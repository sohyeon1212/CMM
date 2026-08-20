from __future__ import annotations

import pytest
from cmm.core.media import (
    PRESET_MEDIA,
    Medium,
    apply_medium,
    carbon_minimal,
    glucose_minimal,
    preset_medium,
)


def test_glucose_aerobic_growth(ecoli_core):
    apply_medium(ecoli_core, "glucose_aerobic")
    assert ecoli_core.reactions.EX_glc__D_e.lower_bound == pytest.approx(-10.0)
    assert ecoli_core.reactions.EX_o2_e.lower_bound < 0  # oxygen available
    assert ecoli_core.slim_optimize() == pytest.approx(0.8739, abs=1e-3)


def test_glucose_aerobic_reproduces_the_shipped_default_growth(ecoli_core):
    """The applied medium is the model's shipped default minus CO2 uptake.

    Measured: growth is bit-identical to the shipped default (0.8739215070), because CO2
    is a product at the growth optimum and the open uptake bound is never exercised.
    """

    shipped_medium = dict(ecoli_core.medium)
    shipped_growth = ecoli_core.slim_optimize()
    applied = apply_medium(ecoli_core, "glucose_aerobic")
    assert dict(applied) == {
        exchange_id: rate
        for exchange_id, rate in shipped_medium.items()
        if exchange_id != "EX_co2_e"
    }
    assert ecoli_core.slim_optimize() == pytest.approx(shipped_growth, abs=1e-12)


def test_glucose_anaerobic_closes_oxygen(ecoli_core):
    apply_medium(ecoli_core, "glucose_anaerobic")
    assert ecoli_core.reactions.EX_o2_e.lower_bound == pytest.approx(0.0)
    assert ecoli_core.slim_optimize() == pytest.approx(0.2111, abs=1e-3)


def test_anaerobic_oxygen_zero_is_explicit():
    """Anaerobiosis is declared, not left to be inferred from an absent key."""

    medium = glucose_minimal(aerobic=False)
    assert medium.uptake["EX_o2_e"] == 0.0


def test_acetate_medium_switches_carbon_source(ecoli_core):
    apply_medium(ecoli_core, "acetate_aerobic")
    assert ecoli_core.reactions.EX_glc__D_e.lower_bound == pytest.approx(
        0.0
    )  # glucose closed
    assert ecoli_core.reactions.EX_ac_e.lower_bound == pytest.approx(-10.0)
    assert ecoli_core.slim_optimize() == pytest.approx(0.1733, abs=1e-3)


def test_presets_close_co2_uptake_but_not_secretion(ecoli_core):
    apply_medium(ecoli_core, "glucose_aerobic")
    lower, upper = ecoli_core.reactions.EX_co2_e.bounds
    assert lower == pytest.approx(0.0)  # no CO2 fixation
    assert upper > 0  # secretion stays free
    assert ecoli_core.optimize().fluxes["EX_co2_e"] > 0  # respiration unaffected


def test_co2_uptake_is_opt_in(ecoli_core):
    apply_medium(ecoli_core, glucose_minimal(co2_uptake=True))
    assert ecoli_core.reactions.EX_co2_e.lower_bound == pytest.approx(-1000.0)


def test_closing_co2_costs_no_aerobic_growth(ecoli_core):
    """Measured: 0.000% aerobic, -0.247% anaerobic."""

    opened = ecoli_core.copy()
    apply_medium(opened, glucose_minimal(co2_uptake=True))
    apply_medium(ecoli_core, glucose_minimal(co2_uptake=False))
    assert ecoli_core.slim_optimize() == pytest.approx(opened.slim_optimize(), abs=1e-9)


def test_medium_resolve_drops_absent_exchanges(ecoli_core):
    medium = Medium(name="x", uptake={"EX_glc__D_e": 5.0, "EX_nonexistent_e": 99.0})
    resolved = medium.resolve(ecoli_core)
    assert resolved == {"EX_glc__D_e": 5.0}


def test_apply_records_applied_and_dropped_components(ecoli_core):
    with pytest.warns(UserWarning, match="were not applied"):
        applied = apply_medium(ecoli_core, "glucose_aerobic")
    assert applied["EX_glc__D_e"] == pytest.approx(10.0)
    # e_coli_core carries no trace-metal exchanges; they drop, and that is recorded.
    assert "EX_fe2_e" in applied.dropped
    assert "EX_glc__D_e" not in applied.dropped
    record = applied.to_provenance()
    assert record["medium"] == applied.medium
    assert record["applied"] == dict(applied)
    assert set(record["dropped"]) == set(applied.dropped)
    frame = applied.to_frame()
    assert list(frame.columns) == ["exchange_id", "uptake", "status"]
    assert len(frame) == len(applied) + len(applied.dropped)


def test_missing_carbon_source_raises_and_names_it(ecoli_core):
    """A glycerol medium on a model with no glycerol exchange is a different experiment."""

    with pytest.raises(ValueError, match="EX_glyc_e"):
        apply_medium(ecoli_core, "glycerol_aerobic")


def test_required_components_are_inferred_from_growth_limiting_bounds():
    medium = carbon_minimal("EX_glyc_e")
    assert medium.required == frozenset({"EX_glyc_e"})
    inferred = Medium(name="x", uptake={"EX_glc__D_e": 5.0, "EX_nh4_e": 1000.0})
    assert inferred.required == frozenset({"EX_glc__D_e"})


def test_glucose_minimal_uptake_parameter():
    medium = glucose_minimal(glucose=8.0, aerobic=False)
    assert medium.uptake["EX_glc__D_e"] == 8.0
    assert medium.uptake["EX_o2_e"] == 0.0


def test_preset_lookup_unknown_raises():
    with pytest.raises(KeyError, match="unknown preset medium"):
        preset_medium("nonsense")
    assert set(PRESET_MEDIA) >= {"glucose_aerobic", "glucose_anaerobic"}


def test_preset_names_do_not_claim_m9_or_a_standard_medium():
    for medium in PRESET_MEDIA.values():
        assert "M9" not in medium.name
        assert "iJO1366/BiGG mineral set" in medium.name


def test_medium_rejects_negative_uptake(ecoli_core):
    with pytest.raises(ValueError, match="non-negative"):
        Medium("bad", {"EX_glc__D_e": -1.0}).apply_to(ecoli_core)
