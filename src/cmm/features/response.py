"""Flux-response analysis: how one reaction's flux constrains the rest of the network.

The target reaction's flux is fixed at each point of a scan and a *response* reaction is
maximized, so the result is a response curve rather than a single operating point. Two
questions share this scan and therefore share one service:

* ``response=None`` maximizes the model objective (biomass), answering "how sensitive is
  growth to this reaction, and where does it break?" — the robustness reading.
* ``response="EX_succ_e"`` maximizes a product while biomass is recorded alongside,
  answering "how much product can this reaction's flux buy me, and at what growth cost?"

Both readings are LP-only, so they run on any solver, and both return the same
serializable result so the visualization layer and the GUI never re-solve.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from cobra import Model
from cobra.flux_analysis import flux_variability_analysis as _cobra_fva

from cmm.core.condition import Condition
from cmm.core.provenance import run_provenance


@dataclass(frozen=True)
class ResponsePoint:
    """One scan point: the target flux held fixed and what the network could do with it."""

    target_flux: float
    response_flux: float  # nan when the fixed target flux is infeasible
    biomass_flux: float  # nan when infeasible
    status: str

    @property
    def feasible(self) -> bool:
        return self.status == "optimal"


@dataclass(frozen=True)
class ResponseBottleneck:
    """Where the response curve falls away fastest.

    A bottleneck is the interval over which forcing more flux through the target costs the
    most response. ``found`` is False when the curve never declines or is flat within
    tolerance, which is itself a finding: the objective does not care about this reaction
    over the scanned range.
    """

    found: bool
    message: str
    target_flux: float | None = None  # midpoint of the steepest-decline interval
    response_flux: float | None = None
    steepest_decline: float | None = None  # d(response)/d(target), most negative
    decline_interval: tuple[float, float] | None = None
    sensitivity: float | None = None  # magnitude of the steepest decline


@dataclass(frozen=True)
class FluxResponseResult:
    """Response of one reaction to an enforced flux in another."""

    target: str
    response: str
    biomass: str
    points: tuple[ResponsePoint, ...]
    bottleneck: ResponseBottleneck
    wild_type: dict[str, float]  # target/response/biomass at the growth optimum
    metadata: dict[str, object] = field(default_factory=dict)

    def to_frame(self) -> pd.DataFrame:
        """Deterministic export table, one row per scan point."""

        return pd.DataFrame(
            [
                (p.target_flux, p.response_flux, p.biomass_flux, p.status)
                for p in self.points
            ],
            columns=["target_flux", "response_flux", "biomass_flux", "status"],
        )

    def feasible_points(self) -> tuple[ResponsePoint, ...]:
        return tuple(p for p in self.points if p.feasible)

    def feasible_range(self) -> tuple[float, float] | None:
        """Lowest and highest scanned target flux that kept the model solvable."""

        feasible = self.feasible_points()
        if not feasible:
            return None
        return (feasible[0].target_flux, feasible[-1].target_flux)

    def optimum(self) -> ResponsePoint | None:
        """The scanned point with the highest response flux."""

        feasible = self.feasible_points()
        if not feasible:
            return None
        return max(feasible, key=lambda p: p.response_flux)


def _objective_reaction(model: Model) -> str:
    for reaction in model.reactions:
        if reaction.objective_coefficient != 0:
            return reaction.id
    raise ValueError("model has no objective reaction; set model.objective first")


def _require_reaction(model: Model, reaction_id: str, label: str) -> None:
    if reaction_id not in model.reactions:
        raise KeyError(f"{label} reaction {reaction_id!r} is not in the model")


def _scan_range(model: Model, target: str) -> tuple[float, float]:
    """Full feasible flux interval of the target, ignoring the objective (FVA at 0%)."""

    table = _cobra_fva(
        model,
        reaction_list=[model.reactions.get_by_id(target)],
        fraction_of_optimum=0.0,
    )
    row = table.loc[target]
    return float(row["minimum"]), float(row["maximum"])


def _bottleneck(points: tuple[ResponsePoint, ...], tol: float) -> ResponseBottleneck:
    feasible = [p for p in points if p.feasible]
    if len(feasible) < 2:
        return ResponseBottleneck(
            found=False,
            message="fewer than two feasible scan points; no curve to differentiate",
        )

    xs = np.asarray([p.target_flux for p in feasible], dtype=float)
    ys = np.asarray([p.response_flux for p in feasible], dtype=float)
    steps = np.diff(xs)
    if not np.all(np.abs(steps) > 1e-12):
        return ResponseBottleneck(
            found=False, message="scan points are not distinct; cannot take a gradient"
        )

    gradients = np.diff(ys) / steps
    if float(np.max(np.abs(gradients))) < tol:
        return ResponseBottleneck(
            found=False,
            message="response is flat within tolerance; it is insensitive to this reaction",
        )

    index = int(np.argmin(gradients))
    steepest = float(gradients[index])
    if steepest >= -tol:
        return ResponseBottleneck(
            found=False,
            message="response never declines over the scanned range; no bottleneck here",
        )

    return ResponseBottleneck(
        found=True,
        message="steepest decline in response per unit target flux",
        target_flux=float((xs[index] + xs[index + 1]) / 2.0),
        response_flux=float((ys[index] + ys[index + 1]) / 2.0),
        steepest_decline=steepest,
        decline_interval=(float(xs[index]), float(xs[index + 1])),
        sensitivity=abs(steepest),
    )


def flux_response(
    model: Model,
    target: str,
    response: str | None = None,
    *,
    biomass: str | None = None,
    condition: Condition | None = None,
    target_min: float | None = None,
    target_max: float | None = None,
    n_steps: int = 20,
    biomass_fraction: float | None = None,
    tol: float = 1e-9,
) -> FluxResponseResult:
    """Scan an enforced flux through ``target`` and maximize ``response`` at each point.

    ``response`` defaults to the biomass/objective reaction, which gives the robustness
    reading (how growth responds to the target). Passing a product exchange instead gives
    the production reading, with biomass recorded at every point so the growth cost of each
    target flux is visible.

    Maximizing a product with no floor on growth returns solutions that sacrifice growth
    entirely, which is a theoretical ceiling rather than a strain. ``biomass_fraction``
    holds biomass at that fraction of the wild-type optimum throughout the scan, making the
    curve the product a *viable* cell could reach — use it whenever the response is a
    product.

    The scan range defaults to the target's full feasible interval (FVA at 0% of the
    optimum). An explicit ``target_min``/``target_max`` may lie outside the reaction's own
    bounds — the scan then overrides them point by point, which is recorded in the result
    metadata as ``range_outside_bounds`` because it is a deliberate what-if, not the model
    as loaded.

    Infeasible points are kept in the result with ``status`` set and NaN fluxes: the flux
    at which a scan stops being solvable is a result, not an error.
    """

    if n_steps < 2:
        raise ValueError("n_steps must be at least 2")
    if tol < 0:
        raise ValueError("tol must be non-negative")
    if biomass_fraction is not None and not 0.0 <= biomass_fraction <= 1.0:
        raise ValueError("biomass_fraction must be in [0, 1]")

    _require_reaction(model, target, "target")
    if response is not None:
        _require_reaction(model, response, "response")
    if biomass is not None:
        _require_reaction(model, biomass, "biomass")

    with model:
        if condition is not None:
            condition.apply_to(model)

        biomass_id = biomass or _objective_reaction(model)
        response_id = response or biomass_id
        target_reaction = model.reactions.get_by_id(target)
        original_bounds = (
            float(target_reaction.lower_bound),
            float(target_reaction.upper_bound),
        )

        # Wild-type reference: the growth optimum the scan is read against.
        model.objective = model.reactions.get_by_id(biomass_id)
        model.objective_direction = "max"
        wild_type_solution = model.optimize()
        if wild_type_solution.status != "optimal":
            raise ValueError(
                f"wild-type solve is {wild_type_solution.status}; the model does not grow "
                "under this condition, so a response scan has no reference"
            )
        wild_type = {
            "target_flux": float(wild_type_solution.fluxes[target]),
            "response_flux": float(wild_type_solution.fluxes[response_id]),
            "biomass_flux": float(wild_type_solution.fluxes[biomass_id]),
        }

        auto_range = target_min is None or target_max is None
        if target_min is None or target_max is None:
            detected_min, detected_max = _scan_range(model, target)
            lower = detected_min if target_min is None else float(target_min)
            upper = detected_max if target_max is None else float(target_max)
        else:
            lower, upper = float(target_min), float(target_max)

        if not np.isfinite(lower) or not np.isfinite(upper):
            raise ValueError(
                f"scan range for {target!r} is not finite ({lower}, {upper}); "
                "give target_min/target_max explicitly"
            )
        if upper - lower <= tol:
            raise ValueError(
                f"scan range for {target!r} is empty ({lower:g}, {upper:g}); the reaction "
                "is fixed under this condition"
            )

        levels = np.linspace(lower, upper, n_steps)

        # Every scan point maximizes the response under the fixed target flux.
        model.objective = model.reactions.get_by_id(response_id)
        model.objective_direction = "max"

        biomass_floor = (
            None
            if biomass_fraction is None
            else biomass_fraction * wild_type["biomass_flux"]
        )

        points: list[ResponsePoint] = []
        for level in levels:
            with model:
                model.reactions.get_by_id(target).bounds = (float(level), float(level))

                if biomass_floor is not None:
                    biomass_reaction = model.reactions.get_by_id(biomass_id)
                    if biomass_floor > float(biomass_reaction.upper_bound):
                        # The growth floor is unreachable for this reaction's capacity;
                        # record it as infeasible rather than letting cobra raise.
                        points.append(
                            ResponsePoint(
                                target_flux=float(level),
                                response_flux=float("nan"),
                                biomass_flux=float("nan"),
                                status="infeasible",
                            )
                        )
                        continue
                    biomass_reaction.bounds = (
                        max(float(biomass_reaction.lower_bound), float(biomass_floor)),
                        float(biomass_reaction.upper_bound),
                    )

                solution = model.optimize()
                if solution.status == "optimal":
                    points.append(
                        ResponsePoint(
                            target_flux=float(level),
                            response_flux=float(solution.fluxes[response_id]),
                            biomass_flux=float(solution.fluxes[biomass_id]),
                            status="optimal",
                        )
                    )
                else:
                    points.append(
                        ResponsePoint(
                            target_flux=float(level),
                            response_flux=float("nan"),
                            biomass_flux=float("nan"),
                            status=str(solution.status),
                        )
                    )

    scanned = tuple(points)
    provenance = run_provenance(
        model,
        method="flux_response",
        target=target,
        response=response_id,
        biomass=biomass_id,
        condition=condition.name if condition else None,
        n_steps=n_steps,
        biomass_fraction=biomass_fraction,
        target_min=lower,
        target_max=upper,
        auto_range=auto_range,
        range_outside_bounds=bool(
            lower < original_bounds[0] - tol or upper > original_bounds[1] + tol
        ),
        target_bounds=list(original_bounds),
    )
    return FluxResponseResult(
        target=target,
        response=response_id,
        biomass=biomass_id,
        points=scanned,
        bottleneck=_bottleneck(scanned, tol=max(tol, 1e-9)),
        wild_type=wild_type,
        metadata=provenance,
    )
