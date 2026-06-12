from __future__ import annotations

import pytest
from cmm.core.solvers import SolverCapabilityError
from cmm.features.strain_design import optknock, robustknock

SUCC = "EX_succ_e"


@pytest.fixture
def anaerobic_ecoli(ecoli_core):
    ecoli_core.reactions.EX_o2_e.lower_bound = 0.0
    ecoli_core.reactions.EX_glc__D_e.lower_bound = -10.0
    return ecoli_core


def test_optknock_finds_growth_coupled_succinate_design(anaerobic_ecoli):
    result = optknock(anaerobic_ecoli, SUCC, max_knockouts=3, max_solutions=4)
    assert result.method == "optknock"
    assert len(result.designs) >= 1
    best = result.best()
    assert 1 <= len(best.knockouts) <= 3
    assert best.max_product > 5.0  # succinate is forced to a high flux
    assert best.growth > 0.0


def test_robustknock_returns_only_guaranteed_designs(anaerobic_ecoli):
    result = robustknock(anaerobic_ecoli, SUCC, max_knockouts=3, max_solutions=8)
    assert result.method == "robustknock"
    assert len(result.designs) >= 1
    # Every robust design must guarantee product at maximum growth (worst case > 0).
    assert all(d.guaranteed_product > 1e-6 for d in result.designs)
    assert all(d.growth_coupled for d in result.designs)
    # Ranked by guaranteed product (descending).
    guaranteed = [d.guaranteed_product for d in result.designs]
    assert guaranteed == sorted(guaranteed, reverse=True)


def test_strain_design_does_not_mutate_model(anaerobic_ecoli):
    growth = anaerobic_ecoli.slim_optimize()
    optknock(anaerobic_ecoli, SUCC, max_knockouts=2, max_solutions=2)
    assert anaerobic_ecoli.slim_optimize() == pytest.approx(growth, abs=1e-6)


def test_strain_design_handles_uncouplable_product(anaerobic_ecoli):
    # A product that cannot be produced in the medium yields no designs, not a crash.
    anaerobic_ecoli.reactions.EX_succ_e.bounds = (0.0, 0.0)
    result = optknock(anaerobic_ecoli, SUCC, max_knockouts=2, max_solutions=2)
    assert result.designs == ()
    assert result.best() is None


def test_strain_design_requires_milp(anaerobic_ecoli):
    anaerobic_ecoli.solver = "glpk_exact"  # LP only, no MILP
    with pytest.raises(SolverCapabilityError) as exc:
        optknock(anaerobic_ecoli, SUCC, max_knockouts=2)
    assert exc.value.capability == "MILP"
