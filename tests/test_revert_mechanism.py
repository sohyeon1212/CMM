from __future__ import annotations

import pytest
from cmm.core.flux_state import reference_state_pfba
from cmm.features import knockout_mechanism, revert_targets
from cmm.omics.differential import differential_expression

SOURCE = {"g1": 50.0, "g2": 100.0, "g3": 1.0, "g5": 1.0, "gb": 50.0}
TARGET = {"g1": 50.0, "g2": 1.0, "g3": 100.0, "g5": 100.0, "gb": 50.0}


def _setup(model):
    reference = reference_state_pfba(model, name="disease")
    direction = differential_expression(model, SOURCE, TARGET, reference=reference)
    return reference, direction


@pytest.mark.requires_qp
def test_the_disease_branch_knockout_moves_the_changed_reactions_the_desired_way(
    branched_model,
):
    reference, direction = _setup(branched_model)
    mechanism = knockout_mechanism(branched_model, reference, direction, ["g2"])
    (explanation,) = mechanism.explanations
    assert explanation.status == "optimal"
    assert explanation.blocked_reactions == ("R2",)
    outcomes = {row["reaction_id"]: row["outcome"] for row in explanation.changes}
    # R2 is asked down and R3/R5 up; knocking out g2 reroutes flux through the long branch.
    assert outcomes == {"R2": "desired", "R3": "desired", "R5": "desired"}
    counts = explanation.counts()
    assert counts["n_desired"] == 3 and counts["n_opposite"] == 0
    frame = mechanism.to_frame()
    assert set(frame["target_id"]) == {"g2"}
    assert list(mechanism.summary_frame()["blocked_reactions"]) == ["R2"]


@pytest.mark.requires_miqp
def test_recomputed_mts_matches_the_ranking_it_explains(branched_model):
    reference, direction = _setup(branched_model)
    ranking = revert_targets(branched_model, None, reference, direction, method="rmta")
    ranked = {t.target_id: t.detail["mTS"] for t in ranking}
    mechanism = knockout_mechanism(branched_model, reference, direction, list(ranked))
    for explanation in mechanism.explanations:
        if explanation.status == "optimal":
            assert explanation.mts == pytest.approx(
                ranked[explanation.target_id], abs=1e-6
            )


@pytest.mark.requires_qp
def test_a_lethal_or_unknown_target_stays_in_the_result_with_its_status(branched_model):
    reference, direction = _setup(branched_model)
    # With growth required, losing g1 (the only route from the substrate) is lethal.
    branched_model.reactions.get_by_id("BIOMASS").lower_bound = 1.0
    mechanism = knockout_mechanism(
        branched_model, reference, direction, ["g1", "not_a_gene"]
    )
    statuses = {e.target_id: e.status for e in mechanism.explanations}
    assert statuses["not_a_gene"] == "not_a_knockout"
    assert statuses["g1"] == "infeasible"
    summary = mechanism.summary_frame()
    assert list(summary["target_id"]) == ["g1", "not_a_gene"]
    assert mechanism.to_frame().empty


def test_knockout_mechanism_rejects_an_unknown_perturbation_level(branched_model):
    reference, direction = _setup(branched_model)
    with pytest.raises(ValueError, match="perturbation"):
        knockout_mechanism(
            branched_model,
            reference,
            direction,
            ["g2"],
            perturbation="pathway",  # type: ignore[arg-type]
        )
