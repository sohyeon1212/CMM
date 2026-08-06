from __future__ import annotations

import math

import pytest
from cmm.core.condition import Condition, ReactionBound
from cmm.features.response import flux_response

SUCC = "EX_succ_e"
O2 = "EX_o2_e"
BIOMASS = "Biomass_Ecoli_core"


def test_response_defaults_to_the_objective(ecoli_core):
    result = flux_response(ecoli_core, O2, n_steps=6)
    # With no explicit response the scan is the robustness reading: growth vs the target.
    assert result.response == BIOMASS
    assert result.biomass == BIOMASS
    frame = result.to_frame()
    assert list(frame.columns) == [
        "target_flux",
        "response_flux",
        "biomass_flux",
        "status",
    ]
    assert len(frame) == 6
    optimal = frame[frame["status"] == "optimal"]
    assert (optimal["response_flux"] == optimal["biomass_flux"]).all()


def test_oxygen_scan_reproduces_the_known_growth_optimum(ecoli_core):
    # e_coli_core grows fastest around 21.8 mmol/gDW/h oxygen uptake; the scanned optimum
    # must land near it rather than at either end of the range.
    result = flux_response(ecoli_core, O2, n_steps=40)
    optimum = result.optimum()
    assert optimum is not None
    assert optimum.target_flux == pytest.approx(-21.8, abs=2.0)
    assert optimum.response_flux == pytest.approx(0.8739, abs=0.02)
    assert result.wild_type["biomass_flux"] == pytest.approx(0.8739, abs=1e-3)


def test_auto_range_matches_the_reaction_bounds(ecoli_core):
    result = flux_response(ecoli_core, O2, n_steps=4)
    parameters = result.metadata["parameters"]
    assert parameters["auto_range"] is True
    assert parameters["range_outside_bounds"] is False
    # Oxygen uptake is bounded at -1000 by the model but limited to -60 by the network.
    assert parameters["target_min"] == pytest.approx(-60.0, abs=1e-6)


def test_product_response_records_growth_cost(ecoli_core):
    result = flux_response(ecoli_core, "PGI", response=SUCC, n_steps=8)
    assert result.response == SUCC
    assert result.biomass == BIOMASS
    # Maximizing a product with no growth floor sacrifices growth entirely.
    optimum = result.optimum()
    assert optimum is not None
    assert optimum.response_flux > 0
    assert optimum.biomass_flux == pytest.approx(0.0, abs=1e-6)


def test_biomass_fraction_holds_growth_across_the_scan(ecoli_core):
    result = flux_response(
        ecoli_core, "PGI", response=SUCC, biomass_fraction=0.3, n_steps=10
    )
    feasible = result.feasible_points()
    assert feasible
    floor = 0.3 * result.wild_type["biomass_flux"]
    for point in feasible:
        assert point.biomass_flux >= floor - 1e-6
    assert result.metadata["parameters"]["biomass_fraction"] == 0.3


def test_infeasible_points_are_reported_not_raised(ecoli_core):
    # A 30% growth floor is unreachable at the extremes of the PGI range.
    result = flux_response(
        ecoli_core, "PGI", response=SUCC, biomass_fraction=0.3, n_steps=12
    )
    frame = result.to_frame()
    infeasible = frame[frame["status"] != "optimal"]
    assert not infeasible.empty
    assert infeasible["response_flux"].isna().all()
    assert infeasible["biomass_flux"].isna().all()
    # The feasible window is narrower than the scanned window.
    window = result.feasible_range()
    assert window is not None
    assert window[0] > frame["target_flux"].min()


def test_bottleneck_found_on_a_declining_curve(ecoli_core):
    result = flux_response(ecoli_core, O2, n_steps=30)
    bottleneck = result.bottleneck
    # Growth collapses as oxygen uptake is forced beyond the optimum.
    assert bottleneck.found is True
    assert bottleneck.steepest_decline is not None and bottleneck.steepest_decline < 0
    assert bottleneck.sensitivity == pytest.approx(abs(bottleneck.steepest_decline))
    low, high = bottleneck.decline_interval
    assert low < bottleneck.target_flux < high


def test_bottleneck_absent_when_the_response_never_declines(ecoli_core):
    result = flux_response(ecoli_core, "PGI", response=SUCC, n_steps=10)
    assert result.bottleneck.found is False
    assert "never declines" in result.bottleneck.message


def test_flat_response_is_reported_as_insensitive(branched_model):
    # R3 and R2 are parallel routes to the same product, so forcing flux through R3
    # changes nothing about achievable biomass: a genuinely insensitive reaction.
    result = flux_response(branched_model, "R3", n_steps=6)
    assert result.bottleneck.found is False
    assert "insensitive" in result.bottleneck.message
    responses = [p.response_flux for p in result.feasible_points()]
    assert max(responses) == pytest.approx(min(responses), abs=1e-9)


def test_competing_reaction_declines_the_response(toy_model):
    # PRODUCT and BIOMASS draw on the same metabolite, so every unit forced through
    # PRODUCT costs a unit of biomass — a constant, detectable decline.
    result = flux_response(toy_model, "PRODUCT", response="BIOMASS", n_steps=6)
    assert result.bottleneck.found is True
    assert result.bottleneck.steepest_decline == pytest.approx(-1.0, abs=1e-6)


def test_condition_is_applied_and_restored(ecoli_core):
    anaerobic = Condition(
        name="anaerobic",
        bounds=(ReactionBound(reaction_id=O2, lower_bound=0.0),),
    )
    before = ecoli_core.reactions.get_by_id(O2).lower_bound
    result = flux_response(
        ecoli_core, "PGI", response=SUCC, condition=anaerobic, n_steps=6
    )
    assert result.metadata["parameters"]["condition"] == "anaerobic"
    assert ecoli_core.reactions.get_by_id(O2).lower_bound == before
    # Anaerobic growth is far below the aerobic optimum.
    assert result.wild_type["biomass_flux"] < 0.3


def test_scan_does_not_mutate_the_model(ecoli_core):
    bounds_before = {r.id: r.bounds for r in ecoli_core.reactions}
    objective_before = ecoli_core.slim_optimize()
    flux_response(ecoli_core, "PGI", response=SUCC, biomass_fraction=0.5, n_steps=6)
    assert {r.id: r.bounds for r in ecoli_core.reactions} == bounds_before
    assert ecoli_core.slim_optimize() == pytest.approx(objective_before, abs=1e-9)


def test_explicit_range_outside_bounds_is_disclosed(ecoli_core):
    # Glucose uptake is capped at 10 mmol/gDW/h by the model; scanning to 50 is a
    # deliberate what-if that overrides the model's own capacity, so it must be recorded.
    result = flux_response(
        ecoli_core,
        "EX_glc__D_e",
        response=SUCC,
        target_min=-50.0,
        target_max=0.0,
        n_steps=5,
    )
    assert result.metadata["parameters"]["range_outside_bounds"] is True
    assert result.metadata["parameters"]["auto_range"] is False


def test_range_within_bounds_is_not_flagged(ecoli_core):
    result = flux_response(
        ecoli_core, "PGI", response=SUCC, target_min=-5.0, target_max=5.0, n_steps=5
    )
    assert result.metadata["parameters"]["range_outside_bounds"] is False
    assert result.metadata["parameters"]["auto_range"] is False


def test_rejects_unknown_reactions(ecoli_core):
    with pytest.raises(KeyError, match="target reaction"):
        flux_response(ecoli_core, "NOPE")
    with pytest.raises(KeyError, match="response reaction"):
        flux_response(ecoli_core, "PGI", response="NOPE")


def test_rejects_invalid_parameters(ecoli_core):
    with pytest.raises(ValueError, match="n_steps"):
        flux_response(ecoli_core, "PGI", n_steps=1)
    with pytest.raises(ValueError, match="biomass_fraction"):
        flux_response(ecoli_core, "PGI", biomass_fraction=1.5)
    with pytest.raises(ValueError, match="tol"):
        flux_response(ecoli_core, "PGI", tol=-1.0)


def test_rejects_a_fixed_target(ecoli_core):
    ecoli_core.reactions.get_by_id("PGI").bounds = (1.0, 1.0)
    with pytest.raises(ValueError, match="empty"):
        flux_response(ecoli_core, "PGI")


def test_rejects_a_non_growing_model(ecoli_core):
    ecoli_core.reactions.get_by_id("EX_glc__D_e").bounds = (0.0, 0.0)
    with pytest.raises(ValueError, match="does not grow"):
        flux_response(ecoli_core, "PGI", response=SUCC)


def test_provenance_records_the_scan(ecoli_core):
    result = flux_response(ecoli_core, "PGI", response=SUCC, n_steps=5)
    assert result.metadata["parameters"]["method"] == "flux_response"
    assert result.metadata["parameters"]["target"] == "PGI"
    assert result.metadata["parameters"]["response"] == SUCC
    assert len(result.metadata["model_sha256"]) == 64
    assert math.isfinite(result.metadata["parameters"]["target_max"])
