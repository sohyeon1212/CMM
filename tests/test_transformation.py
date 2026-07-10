from __future__ import annotations

import pytest
from cmm.core.flux_state import FluxState, reference_state_pfba
from cmm.core.solvers import SolverCapabilityError
from cmm.features.transformation import direction_from_states, transformation_targets


def _states(branched_model):
    source = reference_state_pfba(branched_model, name="A")  # routes through R2
    target = FluxState(
        {"SUP_A": 10, "R1": 10, "R2": 0, "R3": 10, "R5": 10, "BIOMASS": 10}, name="B"
    )
    return source, target


def test_direction_from_states():
    source = FluxState({"R1": 10.0, "R2": 0.0, "R3": 5.0})
    target = FluxState({"R1": 0.0, "R2": 8.0, "R3": 5.0})
    direction = direction_from_states(source, target)
    assert direction["R1"] == -1  # decreases
    assert direction["R2"] == 1  # increases
    assert direction["R3"] == 0  # unchanged


def test_moma_transformation_ranks_disease_branch_first(branched_model):
    source, target = _states(branched_model)
    ranking = transformation_targets(branched_model, source, target, method="moma")
    assert ranking.method == "transform_moma"
    assert ranking.best().target_id == "g2"
    assert ranking.best().score > 0  # genuinely moves flux toward the target
    scores = {t.target_id: t.score for t in ranking}
    assert scores["g3"] == pytest.approx(0, abs=1e-9)


def test_mta_transformation_ranks_disease_branch_first(branched_model):
    source, target = _states(branched_model)
    ranking = transformation_targets(branched_model, source, target, method="mta")
    assert ranking.method == "transform_mta"
    assert ranking.best().target_id == "g2"


def test_transformation_does_not_mutate_model(branched_model):
    source, target = _states(branched_model)
    growth = branched_model.slim_optimize()
    transformation_targets(branched_model, source, target, method="moma")
    assert branched_model.slim_optimize() == pytest.approx(growth, abs=1e-6)


def test_moma_transformation_requires_qp(branched_model):
    branched_model.solver = "glpk"
    source, target = _states(branched_model)
    with pytest.raises(SolverCapabilityError) as exc:
        transformation_targets(branched_model, source, target, method="moma")
    assert exc.value.capability == "QP"


def test_unknown_method_raises(branched_model):
    source, target = _states(branched_model)
    with pytest.raises(ValueError, match="unknown method"):
        transformation_targets(branched_model, source, target, method="bogus")


def test_unknown_perturbation_raises(branched_model):
    source, target = _states(branched_model)
    with pytest.raises(ValueError, match="perturbation"):
        transformation_targets(
            branched_model, source, target, method="moma", perturbation="bogus"
        )
