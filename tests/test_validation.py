"""Validation: CMM services must equal independent, direct COBRApy computations.

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
from cmm.omics.expression import eflux2

GLC = "EX_glc__D_e"
SUCC = "EX_succ_e"
BIOMASS = "Biomass_Ecoli_core"


def test_fba_matches_cobra(ecoli_core):
    cmm_solution = fba(ecoli_core)
    cobra_solution = ecoli_core.optimize()
    assert cmm_solution.objective_value == pytest.approx(
        cobra_solution.objective_value, abs=1e-6
    )
    assert cmm_solution.fluxes[BIOMASS] == pytest.approx(
        cobra_solution.fluxes[BIOMASS], abs=1e-6
    )


def test_pfba_matches_cobra(ecoli_core):
    cmm_solution = pfba(ecoli_core)
    cobra_solution = cfa.pfba(ecoli_core)
    assert cmm_solution.fluxes[BIOMASS] == pytest.approx(
        cobra_solution.fluxes[BIOMASS], abs=1e-6
    )
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
        assert cmm_ranges[rid].minimum == pytest.approx(
            cobra_table.loc[rid, "minimum"], abs=1e-5
        )
        assert cmm_ranges[rid].maximum == pytest.approx(
            cobra_table.loc[rid, "maximum"], abs=1e-5
        )


def test_moma_matches_cobra(branched_model):
    reference = cfa.pfba(branched_model)
    cmm_reference = FluxState(dict(reference.fluxes), name="wt")
    with branched_model:
        branched_model.reactions.R2.bounds = (0.0, 0.0)
        cmm_result = moma(branched_model, cmm_reference, linear=False)
    with branched_model:
        branched_model.reactions.R2.bounds = (0.0, 0.0)
        cobra_result = cfa.moma(branched_model, solution=reference, linear=False)
    assert cmm_result.status == "optimal"
    for rid in ("BIOMASS", "R1", "R2", "R3", "R5"):
        assert cmm_result.fluxes[rid] == pytest.approx(
            cobra_result.fluxes[rid], abs=1e-6
        )


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


def test_eflux2_matches_direct_two_stage_qp(branched_model):
    weights = {
        reaction.id: 1.0 for reaction in branched_model.reactions if reaction.genes
    }
    result = eflux2(branched_model, weights)

    with branched_model:
        optimum = branched_model.slim_optimize()
        problem = branched_model.problem
        objective_expression = branched_model.objective.expression
        optimum_constraint = problem.Constraint(
            objective_expression,
            lb=optimum,
            name="validation_eflux2_optimum",
        )
        branched_model.add_cons_vars([optimum_constraint])
        l2 = sum(
            reaction.flux_expression * reaction.flux_expression
            for reaction in branched_model.reactions
        )
        branched_model.objective = problem.Objective(l2, direction="min")
        direct = branched_model.optimize()

    assert result.status == direct.status == "optimal"
    assert result.objective_value == pytest.approx(optimum, abs=1e-7)
    for rid in ("BIOMASS", "R1", "R2", "R3", "R5"):
        assert result.fluxes[rid] == pytest.approx(direct.fluxes[rid], abs=1e-5)
