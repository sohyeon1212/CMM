"""Validation: CMM services must equal direct cobrapy computations on e_coli_core.

These guard against the platform silently diverging from the underlying cobra reference —
the GUI runs exactly these services, so matching cobra here means the GUI numbers are right.
"""

from __future__ import annotations

import cobra.flux_analysis as cfa
import pytest
from cmm.core import fba, fva, pfba
from cmm.core.flux_state import FluxState
from cmm.features.comparison import moma
from cmm.features.production import theoretical_yield

GLC = "EX_glc__D_e"
SUCC = "EX_succ_e"
BIOMASS = "Biomass_Ecoli_core"


def test_fba_matches_cobra(ecoli_core):
    cmm_solution = fba(ecoli_core)
    cobra_solution = ecoli_core.optimize()
    assert cmm_solution.objective_value == pytest.approx(cobra_solution.objective_value, abs=1e-6)
    assert cmm_solution.fluxes[BIOMASS] == pytest.approx(cobra_solution.fluxes[BIOMASS], abs=1e-6)


def test_pfba_matches_cobra(ecoli_core):
    cmm_solution = pfba(ecoli_core)
    cobra_solution = cfa.pfba(ecoli_core)
    assert cmm_solution.fluxes[BIOMASS] == pytest.approx(cobra_solution.fluxes[BIOMASS], abs=1e-6)
    assert sum(abs(v) for v in cmm_solution.fluxes.values()) == pytest.approx(
        cobra_solution.fluxes.abs().sum(), abs=1e-3
    )


def test_fva_matches_cobra(ecoli_core):
    reactions = [BIOMASS, GLC, SUCC, "PFK", "PYK"]
    cmm_ranges = fva(ecoli_core, reactions=reactions, fraction_of_optimum=0.9)
    cobra_table = cfa.flux_variability_analysis(
        ecoli_core, reaction_list=reactions, fraction_of_optimum=0.9
    )
    for rid in reactions:
        assert cmm_ranges[rid].minimum == pytest.approx(cobra_table.loc[rid, "minimum"], abs=1e-5)
        assert cmm_ranges[rid].maximum == pytest.approx(cobra_table.loc[rid, "maximum"], abs=1e-5)


def test_moma_matches_cobra(ecoli_core, unrestricted_qp_solver):
    reference = ecoli_core.optimize()
    cmm_reference = FluxState(dict(reference.fluxes), name="wt")
    with ecoli_core:
        ecoli_core.reactions.PFK.bounds = (0.0, 0.0)
        cmm_result = moma(ecoli_core, cmm_reference, linear=False)
    with ecoli_core:
        ecoli_core.reactions.PFK.bounds = (0.0, 0.0)
        cobra_result = cfa.moma(ecoli_core, solution=reference, linear=False)
    assert cmm_result.status == "optimal"
    assert cmm_result.fluxes[BIOMASS] == pytest.approx(cobra_result.fluxes[BIOMASS], abs=1e-4)


def test_theoretical_yield_matches_manual_cobra(ecoli_core):
    cmm_yield = theoretical_yield(ecoli_core, SUCC, aerobic=True)
    # Manual cobra: fix glucose at -10, maximize succinate, divide by 10.
    with ecoli_core:
        ecoli_core.reactions.get_by_id(GLC).bounds = (-10.0, -10.0)
        ecoli_core.objective = ecoli_core.reactions.get_by_id(SUCC)
        manual_max = ecoli_core.slim_optimize()
    assert cmm_yield.molar_yield == pytest.approx(manual_max / 10.0, abs=1e-4)


def test_apply_medium_matches_cobra_medium(ecoli_core):
    from cmm.core import apply_medium

    apply_medium(ecoli_core, "glucose_anaerobic")
    # cobra's medium view should report the same open uptakes.
    assert "EX_o2_e" not in ecoli_core.medium  # anaerobic: oxygen closed
    assert ecoli_core.medium[GLC] == pytest.approx(10.0)
