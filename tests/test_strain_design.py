from __future__ import annotations

import pytest

# OptKnock/RobustKnock delegate the bilevel MILP to the optional ``straindesign`` package
# (the ``design`` extra). Skip this module gracefully where it is not installed so CI on a
# platform without it stays green instead of erroring at call time.
pytest.importorskip("straindesign")

from cmm.core.solvers import SolverCapabilityError  # noqa: E402
from cmm.features.strain_design import optknock, robustknock  # noqa: E402

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


def test_optknock_and_robustknock_use_distinct_nested_searches(
    anaerobic_ecoli, monkeypatch
):
    """Prevent RobustKnock from regressing to post-filtered OptKnock candidates."""

    import straindesign as sd

    seen: list[str] = []
    real_module = sd.SDModule

    def recording_module(model, module_type, *args, **kwargs):
        seen.append(module_type)
        return real_module(model, module_type, *args, **kwargs)

    monkeypatch.setattr(sd, "SDModule", recording_module)
    monkeypatch.setattr(sd, "compute_strain_designs", lambda *args, **kwargs: None)

    optknock(anaerobic_ecoli, SUCC, max_knockouts=1, max_solutions=1)
    robustknock(anaerobic_ecoli, SUCC, max_knockouts=1, max_solutions=1)

    assert seen == [sd.OPTKNOCK, sd.ROBUSTKNOCK]


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"max_knockouts": 0}, "max_knockouts"),
        ({"max_solutions": 0}, "max_solutions"),
        ({"min_growth": -0.1}, "min_growth"),
    ],
)
def test_strain_design_validates_search_parameters(anaerobic_ecoli, kwargs, message):
    with pytest.raises(ValueError, match=message):
        optknock(anaerobic_ecoli, SUCC, **kwargs)


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
