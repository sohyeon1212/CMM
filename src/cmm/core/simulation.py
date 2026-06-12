"""Flux simulation services."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

from cobra import Model
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis import pfba as _cobra_pfba

from cmm.core.condition import Condition


@dataclass(frozen=True)
class FluxSolution:
    """Serializable result of an optimization run."""

    status: str
    objective_value: float | None
    fluxes: dict[str, float]


@dataclass(frozen=True)
class FluxRange:
    """Minimum and maximum feasible flux for a reaction."""

    minimum: float
    maximum: float


def fba(model: Model, condition: Condition | None = None) -> FluxSolution:
    """Run flux balance analysis and return a plain Python result."""

    with model:
        if condition is not None:
            condition.apply_to(model)
        solution = model.optimize()

        objective_value = solution.objective_value
        if objective_value is not None:
            objective_value = float(objective_value)

        return FluxSolution(
            status=solution.status,
            objective_value=objective_value,
            fluxes={reaction_id: float(value) for reaction_id, value in solution.fluxes.items()},
        )


def pfba(
    model: Model,
    condition: Condition | None = None,
    fraction_of_optimum: float = 1.0,
) -> FluxSolution:
    """Run parsimonious FBA (minimal total flux at the given fraction of the optimum)."""

    with model:
        if condition is not None:
            condition.apply_to(model)
        solution = _cobra_pfba(model, fraction_of_optimum=fraction_of_optimum)
        objective_value = solution.objective_value
        if objective_value is not None:
            objective_value = float(objective_value)
        return FluxSolution(
            status=solution.status,
            objective_value=objective_value,
            fluxes={reaction_id: float(value) for reaction_id, value in solution.fluxes.items()},
        )


def fva(
    model: Model,
    condition: Condition | None = None,
    reactions: Iterable[str] | None = None,
    fraction_of_optimum: float = 1.0,
) -> dict[str, FluxRange]:
    """Run flux variability analysis for selected reactions."""

    if fraction_of_optimum < 0 or fraction_of_optimum > 1:
        raise ValueError("fraction_of_optimum must be between 0 and 1")

    with model:
        if condition is not None:
            condition.apply_to(model)

        reaction_list = None
        if reactions is not None:
            reaction_list = [model.reactions.get_by_id(reaction_id) for reaction_id in reactions]

        table = flux_variability_analysis(
            model,
            reaction_list=reaction_list,
            fraction_of_optimum=fraction_of_optimum,
        )

        return {
            reaction_id: FluxRange(
                minimum=float(row["minimum"]),
                maximum=float(row["maximum"]),
            )
            for reaction_id, row in table.iterrows()
        }
