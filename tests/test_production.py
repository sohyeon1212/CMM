from __future__ import annotations

import pytest
from cmm.core.condition import Condition, ReactionBound
from cmm.features.production import (
    fseof,
    fvseof,
    production_envelope,
    theoretical_yield,
)

GLC = "EX_glc__D_e"
SUCC = "EX_succ_e"
BIOMASS = "Biomass_Ecoli_core"

#: The anaerobic condition, stated the only way 0.4.0 offers: a Condition closing the oxygen
#: exchange, with CO2 uptake closed alongside it exactly as the ``glucose_anaerobic`` preset
#: does. The ``aerobic=False`` flag this replaces was removed in 0.4.0.
ANAEROBIC = Condition(
    name="anaerobic (O2 and CO2 uptake closed)",
    bounds=(
        ReactionBound("EX_o2_e", lower_bound=0.0),
        ReactionBound("EX_co2_e", lower_bound=0.0),
    ),
)


def test_theoretical_yield_succinate_anaerobic(ecoli_core):
    result = theoretical_yield(ecoli_core, SUCC, condition=ANAEROBIC)
    assert result.status == "optimal"
    assert result.substrate == GLC
    assert result.substrate_uptake == pytest.approx(10.0, abs=1e-6)
    # 1.2000 exactly, not the 1.3906 reported before the anaerobic path went through the
    # media layer: the old path closed oxygen by hand and left CO2 uptake open at 1000, so
    # 6.9529 mmol gDW-1 h-1 of CO2 - not available in a closed anaerobic fermentation -
    # inflated the yield by 15.9%.
    assert result.molar_yield == pytest.approx(1.2000, abs=1e-6)
    assert result.co2_exchange >= 0.0  # CO2 is released, never taken up
    assert not result.co2_fixed
    assert result.co2_carbon_fraction == 0.0
    assert not result.exceeds_carbon_ceiling
    assert not result.carbon_imbalance


def test_theoretical_yield_reports_the_condition_it_solved_under(ecoli_core):
    """``condition=`` is the only way to state aeration as of 0.4.0."""

    by_condition = theoretical_yield(ecoli_core, SUCC, condition=ANAEROBIC)
    assert by_condition.molar_yield == pytest.approx(1.2000, abs=1e-6)
    assert by_condition.aerobic is False
    assert by_condition.metadata["parameters"]["condition"] == ANAEROBIC.name


@pytest.mark.parametrize(
    "function, kwargs",
    [
        (theoretical_yield, {}),
        (production_envelope, {"points": 4}),
        (fseof, {"n_steps": 3}),
        (fvseof, {"n_steps": 3}),
    ],
)
def test_aerobic_keyword_was_removed_in_0_4_0(ecoli_core, function, kwargs):
    """Decision 3: no deprecation period - the keyword is gone, not merely discouraged."""

    with pytest.raises(TypeError, match="aerobic"):
        function(ecoli_core, SUCC, aerobic=False, **kwargs)


def test_theoretical_yield_records_the_applied_condition(ecoli_core):
    """A result file must state its own conditions, not just a preset key."""

    result = theoretical_yield(ecoli_core, SUCC, condition=ANAEROBIC)
    applied = result.metadata["parameters"]["applied_condition"]
    assert applied["aerobic"] is False
    assert applied["oxygen_exchange"]["reaction_id"] == "EX_o2_e"
    assert applied["oxygen_exchange"]["lower_bound"] == pytest.approx(0.0)
    assert applied["medium"]["medium"]  # a name, not an empty string
    assert applied["medium"]["applied"][GLC] == pytest.approx(10.0)
    assert {u["reaction_id"] for u in applied["carbon_uptake"]} == {GLC}


def test_theoretical_yield_aerobic_higher(ecoli_core):
    aerobic = theoretical_yield(ecoli_core, SUCC)
    anaerobic = theoretical_yield(ecoli_core, SUCC, condition=ANAEROBIC)
    assert aerobic.molar_yield == pytest.approx(1.638, abs=0.02)
    assert aerobic.molar_yield > anaerobic.molar_yield


def test_theoretical_yield_does_not_mutate_model(ecoli_core):
    before = ecoli_core.reactions.EX_o2_e.lower_bound
    theoretical_yield(ecoli_core, SUCC, condition=ANAEROBIC)
    assert ecoli_core.reactions.EX_o2_e.lower_bound == before
    # Objective restored to biomass.
    assert ecoli_core.slim_optimize() == pytest.approx(0.8739, abs=1e-3)


def test_theoretical_yield_raises_on_closed_substrate(ecoli_core):
    # A closed substrate (lower_bound 0) would divide by zero -> NaN; fail loudly instead.
    ecoli_core.reactions.get_by_id(GLC).bounds = (0.0, 1000.0)
    with pytest.raises(ValueError, match="no uptake capacity"):
        theoretical_yield(ecoli_core, SUCC, substrate=GLC)


def test_production_envelope_is_growth_coupled(ecoli_core):
    envelope = production_envelope(ecoli_core, SUCC, points=15)
    assert envelope.max_growth == pytest.approx(0.8739, abs=1e-2)
    frame = envelope.to_frame()
    # Growth falls as enforced succinate rises (growth-coupled tradeoff).
    assert frame["growth_max"].iloc[0] > frame["growth_max"].iloc[-1]
    assert frame["product_flux"].iloc[-1] > frame["product_flux"].iloc[0]


def test_production_envelope_respects_the_anaerobic_condition(ecoli_core):
    aerobic = production_envelope(ecoli_core, SUCC, points=8)
    before_o2 = ecoli_core.reactions.EX_o2_e.lower_bound
    anaerobic = production_envelope(ecoli_core, SUCC, condition=ANAEROBIC, points=8)

    assert anaerobic.max_growth < aerobic.max_growth
    assert anaerobic.max_product < aerobic.max_product
    assert ecoli_core.reactions.EX_o2_e.lower_bound == before_o2


def test_fseof_finds_known_succinate_targets(ecoli_core):
    result = fseof(ecoli_core, SUCC, BIOMASS, n_steps=8, condition=ANAEROBIC)
    amplification = set(result.amplification_targets())
    # FRD7 (fumarate reductase) and PPC (PEP carboxylase) are textbook succinate
    # amplification targets and must appear.
    assert "FRD7" in amplification
    assert "PPC" in amplification
    # The enforced exchange is a tautological scan signal, not an actionable gene target.
    assert SUCC not in amplification
    assert SUCC in result.amplification_targets(actionable_only=False)


def test_fseof_criterion_string_describes_what_is_computed(ecoli_core):
    """The provenance label claimed monotonicity that ``_classify_trend`` never tested."""

    result = fseof(ecoli_core, SUCC, BIOMASS, n_steps=4, condition=ANAEROBIC)
    criterion = result.metadata["criterion"]
    assert "monotonic" not in criterion
    assert criterion == (
        "endpoint_difference_and_linear_slope_of_flux_magnitude_"
        "without_direction_reversal"
    )
    # The criterion is CMM's, deliberately kept, and must not be presented as Choi's.
    assert result.metadata["criterion_source"] == "CMM (not Choi et al. 2010)"


def test_fseof_trends_shape(ecoli_core):
    result = fseof(ecoli_core, SUCC, BIOMASS, n_steps=6, condition=ANAEROBIC)
    assert result.trends.shape[0] == len(ecoli_core.reactions)
    assert "classification" in result.trends.columns
    assert len(result.enforced_levels) == 6


_FVSEOF_RXNS = ["FRD7", "PPC", "MDH", "FUM", "SUCCt3", "EX_succ_e", "CS", "ACALD"]


def test_fvseof_finds_robust_succinate_targets(ecoli_core):
    result = fvseof(
        ecoli_core,
        SUCC,
        BIOMASS,
        n_steps=4,
        condition=ANAEROBIC,
        reactions=_FVSEOF_RXNS,
    )
    amplification = set(result.amplification_targets())
    # The reductive succinate pathway must be amplification targets...
    assert {"FRD7", "FUM", "MDH"} <= amplification
    # ...and robustly forced (FVA minimum rises, not just the mean).
    robust = set(result.robust_targets())
    assert "FRD7" in robust
    assert (
        robust <= amplification
    )  # robust targets are a subset of amplification targets


def test_fvseof_classifies_on_park_nine_types(ecoli_core):
    """Park et al. classify on the *joint* sign of delta V_avg and delta l_sol.

    ``capacity_slope`` (= delta l_sol) was computed and exported but never selected on.
    """

    result = fvseof(
        ecoli_core,
        SUCC,
        BIOMASS,
        n_steps=4,
        condition=ANAEROBIC,
        reactions=_FVSEOF_RXNS,
    )
    assert set(result.park_type.unique()) <= set(range(1, 10))
    # The band is fixed by delta V_avg and must agree with the three-way classification.
    for rid, park_type in result.park_type.items():
        expected = (
            "amplify" if park_type <= 3 else "knockdown" if park_type <= 6 else "none"
        )
        assert result.classification[rid] == expected
    # The within-band index is fixed by the sign of delta l_sol.
    for rid, park_type in result.park_type.items():
        capacity_slope = result.capacity_slope[rid]
        if park_type in (1, 4, 7):
            assert capacity_slope < 0
        elif park_type in (3, 6, 9):
            assert capacity_slope > 0
    assert result.metadata["criterion_source"] == "Park et al. 2012"
    # The robust flag is CMM's own and must never be labelled Park's.
    assert result.metadata["robust_flag_source"] == "CMM (not in Park et al. 2012)"


def test_fvseof_amplification_targets_use_park_priority_order(ecoli_core):
    """Park: "reactions with smaller values of l_sol received higher priorities"."""

    result = fvseof(
        ecoli_core,
        SUCC,
        BIOMASS,
        n_steps=4,
        condition=ANAEROBIC,
        reactions=_FVSEOF_RXNS,
    )
    targets = result.amplification_targets()
    widths = [result.capacity.loc[rid].mean() for rid in targets]
    assert widths == sorted(widths)
    assert set(targets) == set(result.targets_of_type(1, 2, 3))


def test_fvseof_default_n_steps_meets_parks_minimum(ecoli_core):
    """Park specify n >= 10; the default was 8."""

    import inspect

    assert inspect.signature(fvseof).parameters["n_steps"].default == 10


def test_fvseof_ranges_shape_and_columns(ecoli_core):
    result = fvseof(
        ecoli_core,
        SUCC,
        BIOMASS,
        n_steps=3,
        condition=ANAEROBIC,
        reactions=_FVSEOF_RXNS,
    )
    assert result.mean.shape == (len(_FVSEOF_RXNS), 3)
    assert result.forced.shape == (len(_FVSEOF_RXNS), 3)
    assert result.capacity.shape == (len(_FVSEOF_RXNS), 3)
    # The forced minimum magnitude is never negative (it is a |flux| lower bound).
    assert (result.forced.to_numpy() >= -1e-9).all()
    assert len(result.enforced_levels) == 3


def test_fvseof_zero_yield_product_is_consistent(ecoli_core):
    # A blocked product collapses every enforced level to 0; result stays shape-consistent.
    ecoli_core.reactions.EX_succ_e.bounds = (0.0, 0.0)
    result = fvseof(
        ecoli_core,
        SUCC,
        BIOMASS,
        n_steps=4,
        condition=ANAEROBIC,
        reactions=["FRD7", "PPC"],
    )
    assert result.metadata["max_product"] == pytest.approx(0.0, abs=1e-9)
    assert (
        len(result.enforced_levels) == result.mean.shape[1]
    )  # no level/column mismatch


def test_fseof_zero_yield_product_is_consistent(ecoli_core):
    ecoli_core.reactions.EX_succ_e.bounds = (0.0, 0.0)
    result = fseof(
        ecoli_core,
        SUCC,
        BIOMASS,
        n_steps=4,
        condition=ANAEROBIC,
        reactions=["FRD7", "PPC"],
    )
    scan_columns = [c for c in result.trends.columns if isinstance(c, float)]
    assert result.enforced_levels == (0.0,)
    assert scan_columns == [0.0]
    assert result.amplification_targets() == []


def test_fvseof_accepts_explicit_linear_flux_couplings(ecoli_core):
    # Renamed from "grouping-reaction constraints": Park et al.'s GR constraints are
    # STRING-derived on/off pairs plus a normalised-flux inequality, a different object.
    result = fvseof(
        ecoli_core,
        SUCC,
        BIOMASS,
        n_steps=3,
        condition=ANAEROBIC,
        reactions=["FRD7", "FUM"],
        linear_flux_couplings=[{"FRD7": 1.0, "FUM": -1.0}],
    )
    assert result.metadata["n_linear_flux_couplings"] == 1
    assert set(result.slope.index) == {"FRD7", "FUM"}
    assert result.metadata["coupling_constraint_source"].startswith("CMM")


@pytest.mark.parametrize(
    ("function", "kwargs", "message"),
    [
        (fseof, {"n_steps": 1}, "n_steps"),
        (fseof, {"fraction_min": 0.9, "fraction_max": 0.1}, "fractions"),
        (fvseof, {"biomass_fraction": 0.0}, "biomass_fraction"),
    ],
)
def test_flux_scans_validate_parameters(ecoli_core, function, kwargs, message):
    with pytest.raises(ValueError, match=message):
        function(ecoli_core, SUCC, BIOMASS, condition=ANAEROBIC, **kwargs)


def test_fvseof_does_not_mutate_model(ecoli_core):
    objective = str(ecoli_core.objective.expression)
    growth = ecoli_core.slim_optimize()
    fvseof(
        ecoli_core,
        SUCC,
        BIOMASS,
        n_steps=3,
        condition=ANAEROBIC,
        reactions=_FVSEOF_RXNS,
    )
    assert str(ecoli_core.objective.expression) == objective
    assert ecoli_core.slim_optimize() == pytest.approx(growth, abs=1e-6)


def test_fseof_classifies_reductive_pathway_as_amplify(ecoli_core):
    # MDH/FUM operate in reverse (negative flux) whose magnitude rises with enforced
    # succinate; classifying on magnitude (not signed flux) must label them amplify, and
    # the ethanol pathway (magnitude falls to zero) as knockdown.
    result = fseof(ecoli_core, SUCC, BIOMASS, n_steps=8, condition=ANAEROBIC)
    assert result.trends.loc["MDH", "classification"] == "amplify"
    assert result.trends.loc["FUM", "classification"] == "amplify"
    assert result.trends.loc["ALCD2x", "classification"] == "knockdown"


def test_fseof_knockout_targets_are_disjoint_and_nonempty(ecoli_core):
    result = fseof(ecoli_core, SUCC, BIOMASS, n_steps=8, condition=ANAEROBIC)
    amplify = set(result.amplification_targets())
    knockdown = set(result.knockout_targets())
    assert knockdown  # there are reactions to turn down
    assert amplify.isdisjoint(knockdown)


def test_fseof_does_not_mutate_model(ecoli_core):
    objective = str(ecoli_core.objective.expression)
    direction = ecoli_core.objective_direction
    succ_bounds = ecoli_core.reactions.EX_succ_e.bounds
    growth = ecoli_core.slim_optimize()
    fseof(ecoli_core, SUCC, BIOMASS, n_steps=5, condition=ANAEROBIC)
    assert str(ecoli_core.objective.expression) == objective
    assert ecoli_core.objective_direction == direction
    assert ecoli_core.reactions.EX_succ_e.bounds == succ_bounds
    assert ecoli_core.slim_optimize() == pytest.approx(growth, abs=1e-6)


def test_production_envelope_does_not_mutate_model(ecoli_core):
    growth = ecoli_core.slim_optimize()
    production_envelope(ecoli_core, SUCC, points=8)
    assert ecoli_core.slim_optimize() == pytest.approx(growth, abs=1e-6)


def test_fseof_results_hold_independent_dataframes(ecoli_core):
    # FseofResult.__post_init__ copies the trends frame, so mutating one result's frame
    # cannot corrupt another's (the false-immutability finding).
    r1 = fseof(ecoli_core, SUCC, BIOMASS, n_steps=5, condition=ANAEROBIC)
    r2 = fseof(ecoli_core, SUCC, BIOMASS, n_steps=5, condition=ANAEROBIC)
    r1.trends.loc["FRD7", "classification"] = "TAMPERED"
    assert r2.trends.loc["FRD7", "classification"] != "TAMPERED"


def test_theoretical_yield_discloses_co2_fixation(ecoli_core):
    aerobic = theoretical_yield(ecoli_core, SUCC)
    # Aerobic succinate yield (1.64) exceeds the glucose carbon ceiling (1.5) only via CO2
    # fixation; the result must expose that rather than presenting 1.64 as carbon-from-glucose.
    assert aerobic.carbon_ceiling == pytest.approx(1.5, abs=1e-6)
    assert aerobic.exceeds_carbon_ceiling
    assert aerobic.co2_fixed
    assert aerobic.co2_exchange < 0
    # ...and quantify it, which the bare co2_fixed boolean never did.
    assert aerobic.co2_carbon_fraction > 0.0
    assert aerobic.co2_explains_excess
    assert not aerobic.carbon_imbalance


def test_theoretical_yield_warns_when_co2_uptake_holds_up_the_yield(ecoli_core):
    """The guard that never fired: a CO2-inflated yield *below* the ceiling passed silently.

    Under the condition ``_set_anaerobic`` used to create - oxygen closed by hand, CO2
    uptake left open - the yield is 1.3906 against a ceiling of 1.5, so
    ``exceeds_carbon_ceiling`` is False and nothing complained, while 12.5% of the product
    carbon came from 6.9529 mmol gDW-1 h-1 of CO2 uptake.
    """

    open_co2 = Condition(
        name="oxygen closed, CO2 uptake left open",
        bounds=(ReactionBound("EX_o2_e", lower_bound=0.0),),
    )
    with pytest.warns(UserWarning, match="product carbon from net CO2 uptake"):
        result = theoretical_yield(ecoli_core, SUCC, condition=open_co2)
    assert result.molar_yield == pytest.approx(1.3905777777777781, abs=1e-9)
    assert not result.exceeds_carbon_ceiling  # the old guard still would not fire...
    assert result.co2_carbon_fraction == pytest.approx(0.125, abs=1e-9)  # ...this does


def test_theoretical_yield_carbon_ceiling_counts_every_carbon_uptake(ecoli_core):
    """A co-substrate must enter the denominator, not be attributed to CO2 fixation."""

    ecoli_core.reactions.EX_ac_e.lower_bound = -5.0
    with pytest.warns(UserWarning):
        result = theoretical_yield(ecoli_core, SUCC, substrate=GLC)
    assert result.molar_yield == pytest.approx(1.8571666666666666, abs=1e-9)
    # (10 glucose x 6 C + 5 acetate x 2 C) / (4 C succinate x 10 glucose) = 1.75, not 1.5.
    assert result.carbon_ceiling == pytest.approx(1.75, abs=1e-9)
    assert [u.reaction_id for u in result.co_substrates] == ["EX_ac_e"]
    assert result.co_substrates[0].uptake == pytest.approx(5.0, abs=1e-9)
    assert result.co_substrates[0].carbon_flux == pytest.approx(10.0, abs=1e-9)
    # With the co-substrate counted, the residual excess is exactly the CO2 uptake.
    assert result.excess_carbon == pytest.approx(result.co2_uptake, abs=1e-6)
    assert result.co2_explains_excess
    assert not result.carbon_imbalance


def test_theoretical_yield_rejects_a_non_boundary_reaction(ecoli_core):
    """``theoretical_yield(model, "PFK")`` silently returned 17.66 against a 0.1875 ceiling."""

    with pytest.raises(ValueError, match="not a boundary"):
        theoretical_yield(ecoli_core, "PFK")


def test_carbon_count_uses_stoichiometric_coefficients(ecoli_core):
    """COBRApy's ``_reaction_elements`` weights element counts by coefficient; so must this."""

    from cmm.features.production import _carbon_count

    assert _carbon_count(ecoli_core.reactions.EX_glc__D_e) == pytest.approx(6.0)
    assert _carbon_count(ecoli_core.reactions.EX_succ_e) == pytest.approx(4.0)
    assert _carbon_count(ecoli_core.reactions.EX_ac_e) == pytest.approx(2.0)
    # A mass-balanced internal reaction carries no *net* carbon; the unweighted sum scored
    # PFK as carbon-rich and handed it a meaningless carbon ceiling.
    assert _carbon_count(ecoli_core.reactions.PFK) == pytest.approx(0.0)


def test_detect_substrate_ignores_non_carbon_nutrient(ecoli_core):
    from cmm.features.production import _detect_substrate

    ecoli_core.reactions.EX_pi_e.lower_bound = -15  # tighter than glucose but no carbon
    assert _detect_substrate(ecoli_core) == GLC
