from __future__ import annotations

import pytest
from cmm.core.flux_state import FluxState, reference_state_pfba
from cmm.core.solvers import SolverCapabilityError
from cmm.features import revert as revert_module
from cmm.features._perturbation import apply_perturbation, gene_perturbations
from cmm.features.revert import revert_targets
from cmm.omics.differential import DirectionMap, differential_expression

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
    assert ranking.metadata["formulation"] == "published_valcarcel_2019"
    assert {"bTS", "mTS", "wTS"} <= set(best.detail)


def test_rmta_reaction_perturbation_ranks_r2_first(branched_model):
    reference, direction = _setup(branched_model)
    ranking = revert_targets(
        branched_model,
        None,
        reference,
        direction,
        method="rmta",
        perturbation="reaction",
    )
    assert ranking.best().target_id == "R2"


def test_ranking_is_deterministic(branched_model):
    reference, direction = _setup(branched_model)
    first = revert_targets(branched_model, None, reference, direction)
    second = revert_targets(branched_model, None, reference, direction)
    assert [(t.target_id, round(t.score, 6)) for t in first] == [
        (t.target_id, round(t.score, 6)) for t in second
    ]


def test_published_robustness_score_equation_9():
    assert revert_module._robust_score(2.0, 3.0, -1.0) == pytest.approx(800.0)
    # Outside bTS>0, mTS>0, wTS<0, published rMTA falls back to mTS.
    assert revert_module._robust_score(2.0, 3.0, 1.0) == pytest.approx(2.0)


def test_non_robust_knockout_is_not_ranked_above_disease_target(branched_model):
    reference, direction = _setup(branched_model)
    target_rxns = revert_module._target_reactions(branched_model, reference)
    perts = {p.target_id: p for p in gene_perturbations(branched_model)}
    with apply_perturbation(branched_model, perts["g3"]):
        scores = revert_module._score_knockout(
            branched_model, reference, direction, target_rxns, "rmta", 0.66
        )
    with apply_perturbation(branched_model, perts["g2"]):
        disease_target_scores = revert_module._score_knockout(
            branched_model,
            reference,
            direction,
            target_rxns,
            "rmta",
            0.66,
        )
    assert scores.robust < disease_target_scores.robust


def test_mta_single_method_also_ranks_g2_first(branched_model):
    reference, direction = _setup(branched_model)
    ranking = revert_targets(branched_model, None, reference, direction, method="mta")
    assert ranking.best().target_id == "g2"
    assert ranking.metadata["alpha"] == 0.66


def test_mta_miqp_requires_miqp_and_ranks_g2(branched_model):
    reference, direction = _setup(branched_model)
    ranking = revert_targets(
        branched_model,
        None,
        reference,
        direction,
        method="mta_miqp",
        perturbation="gene",
    )
    assert ranking.method == "revert_mta_miqp"
    assert ranking.best().target_id == "g2"


def test_published_revert_gates_on_miqp_solver(branched_model):
    branched_model.solver = "glpk"  # LP+MILP only, no QP
    reference, direction = _setup(branched_model)
    with pytest.raises(SolverCapabilityError) as exc:
        revert_targets(branched_model, None, reference, direction, method="rmta")
    assert exc.value.capability == "MIQP"


def test_continuous_approximation_gates_on_qp_solver(branched_model):
    branched_model.solver = "glpk"
    reference, direction = _setup(branched_model)
    with pytest.raises(SolverCapabilityError) as exc:
        revert_targets(
            branched_model, None, reference, direction, method="rmta_continuous"
        )
    assert exc.value.capability == "QP"


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"method": "bogus"}, "unknown method"),
        ({"perturbation": "bogus"}, "perturbation"),
        ({"alpha": -0.1}, "alpha"),
        ({"epsilon": -1.0}, "epsilon"),
        ({"parameter_k": 0.0}, "parameter_k"),
    ],
)
def test_revert_validates_public_parameters(branched_model, kwargs, message):
    reference, direction = _setup(branched_model)
    with pytest.raises(ValueError, match=message):
        revert_targets(branched_model, None, reference, direction, **kwargs)


def _published_reference_and_direction():
    reference = FluxState(
        {"r1": 10, "r2": 1, "r3": 8, "r4": 8, "r5": 1, "r6": 9, "r7": 10},
        name="official_cobra_reference",
    )
    direction = DirectionMap(
        {"r1": 0, "r2": 1, "r3": 0, "r4": -1, "r5": 0, "r6": 0, "r7": 0}
    )
    return reference, direction


def test_mta_matches_official_cobra_toolbox_expected_signs(published_mta_model):
    reference, direction = _published_reference_and_direction()
    ranking = revert_targets(
        published_mta_model,
        None,
        reference,
        direction,
        method="mta",
        alpha=0.66,
        epsilon=0.01,
        transcript_separator=".",
    )
    scores = {target.target_id: target.score for target in ranking}
    assert scores["g2"] < 0
    assert scores["g4"] > 0


def test_rmta_matches_official_cobra_toolbox_positive_targets(published_mta_model):
    reference, direction = _published_reference_and_direction()
    gene_ranking = revert_targets(
        published_mta_model,
        None,
        reference,
        direction,
        method="rmta",
        alpha=0.4,
        epsilon=0.01,
        transcript_separator=".",
        targets=["g2", "g4"],
    )
    gene_scores = {target.target_id: target.score for target in gene_ranking}
    assert gene_scores["g4"] > 0

    reaction_ranking = revert_targets(
        published_mta_model,
        None,
        reference,
        direction,
        method="rmta",
        alpha=0.4,
        epsilon=0.01,
        perturbation="reaction",
        targets=["r2", "r3", "r6"],
    )
    reaction_scores = {target.target_id: target.score for target in reaction_ranking}
    assert reaction_scores["r3"] > 0
