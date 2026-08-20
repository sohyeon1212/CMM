from __future__ import annotations

import pytest
from cmm.core.condition import Condition, ObjectiveSpec, ReactionBound
from cmm.core.simulation import fba, fva, pfba


def test_pfba_matches_cobra(ecoli_core):
    import cobra.flux_analysis as cfa

    mine = pfba(ecoli_core)
    reference = cfa.pfba(ecoli_core)
    assert mine.status == "optimal"
    assert mine.fluxes["Biomass_Ecoli_core"] == pytest.approx(
        reference.fluxes["Biomass_Ecoli_core"], abs=1e-6
    )
    # pFBA minimizes total flux; the objective value is that minimal sum.
    assert sum(abs(v) for v in mine.fluxes.values()) == pytest.approx(
        reference.fluxes.abs().sum(), abs=1e-3
    )


def test_pfba_total_flux_not_greater_than_fba(ecoli_core):
    fba_solution = fba(ecoli_core)
    pfba_solution = pfba(ecoli_core)
    # Same growth, but pFBA's total flux is <= an arbitrary FBA vertex.
    assert pfba_solution.fluxes["Biomass_Ecoli_core"] == pytest.approx(
        fba_solution.fluxes["Biomass_Ecoli_core"], abs=1e-4
    )
    assert (
        sum(abs(v) for v in pfba_solution.fluxes.values())
        <= sum(abs(v) for v in fba_solution.fluxes.values()) + 1e-6
    )


def test_fba_returns_plain_flux_solution(toy_model):
    solution = fba(toy_model)

    assert solution.status == "optimal"
    assert solution.objective_value == pytest.approx(10)
    assert solution.fluxes["BIOMASS"] == pytest.approx(10)


def test_fba_accepts_condition(toy_model):
    condition = Condition(
        bounds=(ReactionBound("SOURCE_A", upper_bound=2),),
        objective=ObjectiveSpec({"PRODUCT": 1.0}),
    )

    solution = fba(toy_model, condition)

    assert solution.objective_value == pytest.approx(2)
    assert solution.fluxes["PRODUCT"] == pytest.approx(2)
    assert toy_model.reactions.SOURCE_A.upper_bound == 10


def test_fva_returns_flux_ranges(toy_model):
    ranges = fva(
        toy_model,
        reactions=("BIOMASS", "PRODUCT"),
        fraction_of_optimum=0,
    )

    assert ranges["BIOMASS"].minimum == pytest.approx(0)
    assert ranges["BIOMASS"].maximum == pytest.approx(10)
    assert ranges["PRODUCT"].minimum == pytest.approx(0)
    assert ranges["PRODUCT"].maximum == pytest.approx(10)


def test_fva_carries_provenance_like_its_siblings(toy_model):
    result = fva(toy_model, reactions=("BIOMASS",), fraction_of_optimum=0.5)

    assert len(result.metadata["model_sha256"]) == 64
    assert result.metadata["parameters"]["method"] == "fva"
    assert result.metadata["parameters"]["fraction_of_optimum"] == 0.5
    assert result.metadata["parameters"]["loopless"] is False
    assert result.metadata["timestamp_utc"].endswith("Z")
    frame = result.to_frame()
    assert list(frame.columns) == ["reaction_id", "minimum", "maximum"]
    assert list(frame["reaction_id"]) == ["BIOMASS"]


def test_fva_forwards_loopless(ecoli_core):
    """Loopless FVA removes the FRD7/SUCDi thermodynamic loop; plain FVA does not."""

    looped = fva(ecoli_core, reactions=("FRD7",), fraction_of_optimum=0.0)
    loopless = fva(
        ecoli_core, reactions=("FRD7",), fraction_of_optimum=0.0, loopless=True
    )

    assert looped["FRD7"].maximum == pytest.approx(1000.0)
    assert loopless["FRD7"].maximum < looped["FRD7"].maximum
    assert loopless.metadata["parameters"]["loopless_algorithm"] == "cycleFreeFlux"
    assert looped.metadata["parameters"]["loopless_algorithm"] is None


def test_fva_rejects_an_unknown_loopless_algorithm(toy_model):
    with pytest.raises(ValueError, match="fastSNP"):
        fva(toy_model, loopless="nonsense")


def test_fva_does_not_spawn_a_pool_for_a_small_model(toy_model):
    """Small models run in-process: a pool is pure overhead and breaks on macOS spawn."""

    result = fva(toy_model)
    assert result.metadata["parameters"]["processes"] == 1

    forced = fva(toy_model, processes=2)
    assert forced.metadata["parameters"]["processes"] == 2


def test_fva_validates_fraction(toy_model):
    with pytest.raises(ValueError, match="fraction_of_optimum"):
        fva(toy_model, fraction_of_optimum=1.5)
