from __future__ import annotations

import pytest
from cmm.core.flux_state import reference_state_pfba
from cmm.core.solvers import SolverCapabilityError
from cmm.features import revert as revert_module
from cmm.features._perturbation import apply_perturbation, gene_perturbations
from cmm.features.revert import revert_targets
from cmm.omics.differential import differential_expression

SOURCE = {"g1": 50.0, "g2": 100.0, "g3": 1.0, "g5": 1.0, "gb": 50.0}
TARGET = {"g1": 50.0, "g2": 1.0, "g3": 100.0, "g5": 100.0, "gb": 50.0}


def _setup(model):
    reference = reference_state_pfba(model, name="disease")
    direction = differential_expression(model, SOURCE, TARGET, reference=reference)
    return reference, direction


def test_rmta_ranks_disease_branch_gene_first(branched_model):
    reference, direction = _setup(branched_model)
    ranking = revert_targets(
        branched_model, None, reference, direction, method="rmta", perturbation="gene"
    )
    best = ranking.best()
    assert best.target_id == "g2"  # knocking out the disease branch reverts to healthy
    assert best.score > 0
    # Healthy-branch knockouts must not be predicted as normalization targets.
    scores = {t.target_id: t.score for t in ranking}
    assert scores["g3"] == pytest.approx(0, abs=1e-9)
    assert scores["g5"] == pytest.approx(0, abs=1e-9)
    assert scores["g2"] > scores["g1"]


def test_rmta_reaction_perturbation_ranks_r2_first(branched_model):
    reference, direction = _setup(branched_model)
    ranking = revert_targets(
        branched_model, None, reference, direction, method="rmta", perturbation="reaction"
    )
    assert ranking.best().target_id == "R2"


def test_ranking_is_deterministic(branched_model):
    reference, direction = _setup(branched_model)
    first = revert_targets(branched_model, None, reference, direction)
    second = revert_targets(branched_model, None, reference, direction)
    assert [(t.target_id, round(t.score, 6)) for t in first] == [
        (t.target_id, round(t.score, 6)) for t in second
    ]


def test_robustness_gate_zeros_non_robust_knockout(branched_model):
    # g3/g5 help only the best case (best_ts>0) but not the worst case, so rTS must be 0.
    reference, direction = _setup(branched_model)
    target_rxns = revert_module._target_reactions(branched_model, reference)
    perts = {p.target_id: p for p in gene_perturbations(branched_model)}
    with apply_perturbation(branched_model, perts["g3"]):
        scores = revert_module._score_knockout(
            branched_model, reference, direction, target_rxns, "rmta", 0.9
        )
    assert scores.best > 0
    assert scores.robust == pytest.approx(0, abs=1e-9)


def test_mta_single_method_also_ranks_g2_first(branched_model):
    reference, direction = _setup(branched_model)
    ranking = revert_targets(branched_model, None, reference, direction, method="mta")
    assert ranking.best().target_id == "g2"
    assert ranking.metadata["alpha"] == 0.9


def test_mta_miqp_requires_miqp_and_ranks_g2(branched_model):
    reference, direction = _setup(branched_model)
    ranking = revert_targets(
        branched_model, None, reference, direction, method="mta_miqp", perturbation="gene"
    )
    assert ranking.method == "revert_mta_miqp"
    assert ranking.best().target_id == "g2"


def test_revert_gates_on_qp_solver(branched_model):
    branched_model.solver = "glpk"  # LP+MILP only, no QP
    reference, direction = _setup(branched_model)
    with pytest.raises(SolverCapabilityError) as exc:
        revert_targets(branched_model, None, reference, direction, method="rmta")
    assert exc.value.capability == "QP"
