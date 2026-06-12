from __future__ import annotations

import pytest
from cmm.core.media import (
    PRESET_MEDIA,
    Medium,
    apply_medium,
    glucose_minimal,
    preset_medium,
)


def test_glucose_aerobic_growth(ecoli_core):
    apply_medium(ecoli_core, "glucose_aerobic")
    assert ecoli_core.reactions.EX_glc__D_e.lower_bound == pytest.approx(-10.0)
    assert ecoli_core.reactions.EX_o2_e.lower_bound < 0  # oxygen available
    assert ecoli_core.slim_optimize() == pytest.approx(0.8739, abs=1e-3)


def test_glucose_anaerobic_closes_oxygen(ecoli_core):
    apply_medium(ecoli_core, "glucose_anaerobic")
    assert ecoli_core.reactions.EX_o2_e.lower_bound == pytest.approx(0.0)
    assert ecoli_core.slim_optimize() == pytest.approx(0.2117, abs=1e-3)


def test_acetate_medium_switches_carbon_source(ecoli_core):
    apply_medium(ecoli_core, "acetate_aerobic")
    assert ecoli_core.reactions.EX_glc__D_e.lower_bound == pytest.approx(0.0)  # glucose closed
    assert ecoli_core.reactions.EX_ac_e.lower_bound == pytest.approx(-10.0)
    assert ecoli_core.slim_optimize() == pytest.approx(0.1733, abs=1e-3)


def test_medium_resolve_drops_absent_exchanges(ecoli_core):
    medium = Medium(name="x", uptake={"EX_glc__D_e": 5.0, "EX_nonexistent_e": 99.0})
    resolved = medium.resolve(ecoli_core)
    assert resolved == {"EX_glc__D_e": 5.0}


def test_glucose_minimal_uptake_parameter():
    medium = glucose_minimal(glucose=8.0, aerobic=False)
    assert medium.uptake["EX_glc__D_e"] == 8.0
    assert "EX_o2_e" not in medium.uptake


def test_preset_lookup_unknown_raises():
    with pytest.raises(KeyError, match="unknown preset medium"):
        preset_medium("nonsense")
    assert set(PRESET_MEDIA) >= {"glucose_aerobic", "glucose_anaerobic"}
