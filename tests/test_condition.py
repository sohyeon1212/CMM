from __future__ import annotations

import pytest
from cmm.core.condition import Condition, ObjectiveSpec, ReactionBound, with_condition


def test_condition_applies_bounds_and_objective(toy_model):
    condition = Condition(
        name="product-run",
        bounds=(ReactionBound("SOURCE_A", upper_bound=4),),
        objective=ObjectiveSpec({"PRODUCT": 1.0}),
    )

    condition.apply_to(toy_model)

    assert toy_model.reactions.SOURCE_A.upper_bound == 4
    assert toy_model.objective_direction == "max"
    assert toy_model.slim_optimize() == pytest.approx(4)


def test_with_condition_does_not_mutate_original_model(toy_model):
    condition = Condition(bounds=(ReactionBound("SOURCE_A", upper_bound=3),))

    copied = with_condition(toy_model, condition)

    assert copied.reactions.SOURCE_A.upper_bound == 3
    assert toy_model.reactions.SOURCE_A.upper_bound == 10


def test_invalid_bound_is_rejected():
    condition = Condition(bounds=(ReactionBound("R1", lower_bound=5, upper_bound=1),))

    with pytest.raises(ValueError, match="lower_bound"):
        condition.validate()


def test_condition_bounds_set_atomically_when_crossing(toy_model):
    # SOURCE_A starts at (0, 10). Raising both bounds above the current upper yields a valid
    # final pair (20, 30), but assigning lower=20 before upper would transiently exceed the
    # old upper=10 and raise mid-assignment. apply_to must assign bounds atomically (R11).
    condition = Condition(
        name="shift-up",
        bounds=(ReactionBound("SOURCE_A", lower_bound=20.0, upper_bound=30.0),),
    )

    condition.apply_to(toy_model)  # must not raise

    assert toy_model.reactions.SOURCE_A.bounds == (20.0, 30.0)
