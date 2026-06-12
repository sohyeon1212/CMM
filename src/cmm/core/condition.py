"""Condition data model for applying bounds and objectives to cobra models."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Literal

from cobra import Model

ObjectiveDirection = Literal["max", "min"]


@dataclass(frozen=True)
class ReactionBound:
    """A bound override for one reaction."""

    reaction_id: str
    lower_bound: float | None = None
    upper_bound: float | None = None

    def validate(self) -> None:
        if not self.reaction_id:
            raise ValueError("reaction_id must not be empty")
        if (
            self.lower_bound is not None
            and self.upper_bound is not None
            and self.lower_bound > self.upper_bound
        ):
            raise ValueError(
                f"lower_bound must be <= upper_bound for reaction {self.reaction_id!r}"
            )


@dataclass(frozen=True)
class ObjectiveSpec:
    """A linear objective definition."""

    coefficients: Mapping[str, float]
    direction: ObjectiveDirection = "max"

    def validate(self) -> None:
        if self.direction not in {"max", "min"}:
            raise ValueError("objective direction must be 'max' or 'min'")
        if not self.coefficients:
            raise ValueError("objective coefficients must not be empty")
        for reaction_id, coefficient in self.coefficients.items():
            if not reaction_id:
                raise ValueError("objective reaction_id must not be empty")
            if coefficient == 0:
                raise ValueError(
                    f"objective coefficient for reaction {reaction_id!r} must not be 0"
                )


@dataclass(frozen=True)
class Condition:
    """A reusable set of model constraints for a simulation run."""

    name: str = "default"
    bounds: tuple[ReactionBound, ...] = field(default_factory=tuple)
    objective: ObjectiveSpec | None = None
    notes: str = ""

    def validate(self) -> None:
        if not self.name:
            raise ValueError("condition name must not be empty")
        for bound in self.bounds:
            bound.validate()
        if self.objective is not None:
            self.objective.validate()

    def apply_to(self, model: Model) -> None:
        """Apply this condition to a cobra model in place."""

        self.validate()

        for bound in self.bounds:
            reaction = model.reactions.get_by_id(bound.reaction_id)
            new_lower = bound.lower_bound if bound.lower_bound is not None else reaction.lower_bound
            new_upper = bound.upper_bound if bound.upper_bound is not None else reaction.upper_bound
            # Assign atomically: cobra validates lower <= upper on each scalar assignment, so a
            # new value that crosses the current opposite bound would raise mid-assignment (R11).
            reaction.bounds = (new_lower, new_upper)

        if self.objective is not None:
            reactions = {
                model.reactions.get_by_id(reaction_id): coefficient
                for reaction_id, coefficient in self.objective.coefficients.items()
            }
            model.objective = reactions
            model.objective_direction = self.objective.direction


def with_condition(model: Model, condition: Condition | None) -> Model:
    """Return a copy of model with the condition applied."""

    copied = model.copy()
    if condition is not None:
        condition.apply_to(copied)
    return copied
