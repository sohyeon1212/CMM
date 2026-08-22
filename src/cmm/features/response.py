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

**Method and attribution.** The scan is *in silico flux response analysis* as described by
Lee KH, Park JH, Kim TY, Kim HU & Lee SY (2007), "Systems metabolic engineering of
*Escherichia coli* for L-threonine production", *Mol Syst Biol* 3:149
(doi:10.1038/msb4100196). When the response is biomass the same scan is the robustness
analysis of Edwards JS & Palsson BO (2000), "Robustness analysis of the *Escherichia coli*
metabolic network", *Biotechnol Prog* 16(6):927-939.

**What is reported from the curve.** ``r(L)``, the maximum response at an enforced target
flux ``L``, is the optimal-value function of an LP parameterised in one right-hand side, so
it is concave and piecewise linear. Its derivative is the *shadow price* of the constraint
``v_target = L`` — a dual variable the LP already returns exactly, not something to estimate
by differencing a grid. Following phenotype phase plane analysis (Edwards JS, Ramakrishna R
& Palsson BO (2002), "Characterizing the metabolic phenotype: a phenotype phase plane
analysis", *Biotechnol Bioeng* 77(1):27-36), where phases are regions of constant shadow
price and phase boundaries are the levels at which the shadow price changes, this module
reports the phases of ``r`` (:class:`ResponsePhase`) and the *response limit*
(:class:`ResponseLimit`): the phase boundary at which the shadow price first falls below
``-limit_threshold``. Because ``r`` is concave the shadow price is non-increasing, so that
crossing is unique and there is no tie to break.

At a phase boundary the LP is dually degenerate and the dual returned depends on the optimal
basis. Shadow prices are therefore always read from the interior of a piece, never at a
kink, and a boundary is located as the intersection of the two neighbouring tangents.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from cobra import Model
from cobra.exceptions import OptimizationError
from cobra.flux_analysis import flux_variability_analysis as _cobra_fva

from cmm.core.condition import Condition
from cmm.core.provenance import run_provenance

_FIXED_TARGET_CONSTRAINT = "cmm_flux_response_fixed_target"


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
class ResponsePhase:
    """One linear piece of the response curve: an interval of constant shadow price.

    ``shadow_price`` is the exact ``d(response)/d(target)`` on this piece, read as the dual
    of the enforced-flux constraint from the interior of the piece. It is a property of the
    network, not of the scan grid.
    """

    target_low: float
    target_high: float
    shadow_price: float
    response_low: float
    response_high: float

    @property
    def width(self) -> float:
        return self.target_high - self.target_low


@dataclass(frozen=True)
class ResponseLimit:
    """The phase boundary at which the response starts to fall faster than ``threshold``.

    ``found`` is False when the shadow price never falls below ``-threshold`` over the
    feasible domain, which is itself a finding: within the stated threshold the objective
    does not pay for flux through this reaction. The location is the intersection of the two
    neighbouring tangents, so it is exact rather than a grid midpoint and does not move with
    ``n_steps``.
    """

    found: bool
    message: str
    target_flux: float | None = None
    response_flux: float | None = None
    shadow_price_before: float | None = None
    shadow_price_after: float | None = None
    threshold: float = 1e-6


@dataclass(frozen=True)
class FluxResponseResult:
    """Response of one reaction to an enforced flux in another."""

    target: str
    response: str
    biomass: str
    points: tuple[ResponsePoint, ...]
    phases: tuple[ResponsePhase, ...]
    limit: ResponseLimit
    wild_type: dict[str, float]  # target/response/biomass at the growth optimum
    feasible_domain: tuple[float, float] | None = None  # exact, from FVA
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

    def phases_frame(self) -> pd.DataFrame:
        """Deterministic export table, one row per response phase."""

        return pd.DataFrame(
            [
                (
                    p.target_low,
                    p.target_high,
                    p.response_low,
                    p.response_high,
                    p.shadow_price,
                )
                for p in self.phases
            ],
            columns=[
                "target_low",
                "target_high",
                "response_low",
                "response_high",
                "shadow_price",
            ],
        )

    def summary_frame(self) -> pd.DataFrame:
        """One-row export of the domain, reference state and response limit.

        The curve and phase tables deliberately stay separate; this frame carries every
        scalar a report may quote so none has to be reconstructed from plotted pixels.
        """

        domain_low, domain_high = self.feasible_domain or (None, None)
        return pd.DataFrame(
            [
                {
                    "target": self.target,
                    "response": self.response,
                    "biomass": self.biomass,
                    "feasible_domain_low": domain_low,
                    "feasible_domain_high": domain_high,
                    "wild_type_target_flux": self.wild_type.get("target_flux"),
                    "wild_type_response_flux": self.wild_type.get("response_flux"),
                    "wild_type_biomass_flux": self.wild_type.get("biomass_flux"),
                    "limit_found": self.limit.found,
                    "limit_message": self.limit.message,
                    "limit_target_flux": self.limit.target_flux,
                    "limit_response_flux": self.limit.response_flux,
                    "limit_shadow_price_before": self.limit.shadow_price_before,
                    "limit_shadow_price_after": self.limit.shadow_price_after,
                    "limit_threshold": self.limit.threshold,
                    "n_points": len(self.points),
                    "n_phases": len(self.phases),
                }
            ]
        )

    def feasible_points(self) -> tuple[ResponsePoint, ...]:
        return tuple(p for p in self.points if p.feasible)

    def feasible_range(self) -> tuple[float, float] | None:
        """Exact feasible interval of the scanned reaction.

        Computed by FVA on the target under the scan's own constraints — condition, growth
        floor, and any explicit ``target_min``/``target_max`` window — so it is a real bound
        rather than "the first and last grid point that happened to solve". It does not move
        with ``n_steps``. ``None`` when nothing was feasible.
        """

        return self.feasible_domain

    def shadow_price_at(self, target_flux: float) -> float:
        """Exact ``d(response)/d(target)`` at ``target_flux``, from the phase structure.

        The phases are built from LP duals, so this is a lookup rather than an estimate and
        needs no further solve. A level on a phase boundary takes the price of the phase to
        its right (the derivative is one-sided at a kink); the last phase includes its own
        upper endpoint. Raises when ``target_flux`` lies outside the feasible domain.
        """

        if not self.phases:
            raise ValueError("no response phases were computed; nothing is feasible")
        low = self.phases[0].target_low
        high = self.phases[-1].target_high
        if not low - 1e-9 <= target_flux <= high + 1e-9:
            raise ValueError(
                f"{target_flux:g} is outside the feasible domain ({low:g}, {high:g})"
            )
        for phase in self.phases:
            if target_flux < phase.target_high:
                return phase.shadow_price
        return self.phases[-1].shadow_price

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


def _scan_range(model: Model, target: str) -> tuple[float, float] | None:
    """Exact feasible flux interval of the target under the model *as currently constrained*.

    FVA at 0% of the optimum, so the objective is free; the caller is responsible for having
    already applied the condition and any growth floor, because both narrow this interval.
    Returns ``None`` when the model is infeasible as constrained.
    """

    try:
        table = _cobra_fva(
            model,
            reaction_list=[model.reactions.get_by_id(target)],
            fraction_of_optimum=0.0,
            processes=1,
        )
    except (OptimizationError, ValueError):
        return None
    row = table.loc[target]
    minimum, maximum = float(row["minimum"]), float(row["maximum"])
    if not np.isfinite(minimum) or not np.isfinite(maximum):
        return None
    return minimum, maximum


@dataclass(frozen=True)
class _Solve:
    """One parametric LP solve: the response value and its exact slope at ``level``."""

    level: float
    value: float
    shadow_price: float
    status: str


def _solve_at(model: Model, target: str, level: float) -> _Solve:
    """Fix ``v_target = level``, maximize the response, and read the dual of that equality.

    The dual of ``v_target == level`` *is* ``d(response)/d(level)``: no finite differencing,
    no neighbouring point, no grid. The fixing is expressed as a named constraint rather
    than a bound so the dual is addressable, and the reaction's own box is widened so the
    named constraint is the binding one — the flux is pinned to ``level`` either way, so the
    feasible set is unchanged.
    """

    with model:
        reaction = model.reactions.get_by_id(target)
        constraint = model.problem.Constraint(
            reaction.flux_expression,
            lb=float(level),
            ub=float(level),
            name=_FIXED_TARGET_CONSTRAINT,
        )
        reaction.bounds = (-1e6, 1e6)
        model.add_cons_vars([constraint])
        model.solver.update()
        solution = model.optimize()
        if solution.status != "optimal":
            return _Solve(
                float(level), float("nan"), float("nan"), str(solution.status)
            )
        try:
            dual = float(model.constraints[_FIXED_TARGET_CONSTRAINT].dual)
        except Exception:  # pragma: no cover - solver without duals for this constraint
            dual = float("nan")
        return _Solve(float(level), float(solution.objective_value), dual, "optimal")


def _phases(
    model: Model,
    target: str,
    lower: float,
    upper: float,
    *,
    slope_tol: float = 1e-6,
    value_tol: float = 1e-6,
    x_tol: float = 1e-9,
    max_depth: int = 40,
) -> tuple[ResponsePhase, ...]:
    """Locate every kink of the concave piecewise-linear ``r(L)`` on ``[lower, upper]``.

    On any sub-interval whose endpoint shadow prices agree, ``r`` is a single line and there
    is no kink inside. Where they differ, the two tangents intersect at a candidate boundary
    ``L*``; if ``r(L*)`` equals the tangent intersection value the interval held exactly one
    kink, otherwise there is more than one and the two halves are recursed. The number of
    solves is set by the number of boundaries, not by any ``n_steps``.
    """

    breaks: list[tuple[float, float, float]] = []  # (level, value, slope_after)

    def recurse(a: _Solve, b: _Solve, depth: int) -> None:
        if depth > max_depth or (b.level - a.level) <= x_tol:
            return
        if a.status != "optimal" or b.status != "optimal":
            return
        if not np.isfinite(a.shadow_price) or not np.isfinite(b.shadow_price):
            return
        if abs(a.shadow_price - b.shadow_price) <= slope_tol:
            return  # one linear piece, no kink inside

        def split() -> None:
            middle = _solve_at(model, target, 0.5 * (a.level + b.level))
            recurse(a, middle, depth + 1)
            recurse(middle, b, depth + 1)

        denominator = a.shadow_price - b.shadow_price
        crossing = (
            b.value - a.value - b.shadow_price * b.level + a.shadow_price * a.level
        ) / denominator
        if not a.level - x_tol <= crossing <= b.level + x_tol:
            split()
            return
        predicted = a.value + a.shadow_price * (crossing - a.level)
        boundary = _solve_at(model, target, float(np.clip(crossing, a.level, b.level)))
        if boundary.status == "optimal" and abs(boundary.value - predicted) <= (
            value_tol * max(1.0, abs(predicted))
        ):
            breaks.append((boundary.level, boundary.value, b.shadow_price))
            return
        split()

    left = _solve_at(model, target, lower)
    right = _solve_at(model, target, upper)
    if left.status != "optimal" or right.status != "optimal":
        return ()
    recurse(left, right, 0)
    breaks.sort(key=lambda item: item[0])

    # Recursion can locate the same kink from two sides; keep the first of each cluster.
    unique: list[tuple[float, float, float]] = []
    for item in breaks:
        if unique and abs(item[0] - unique[-1][0]) < 1e-6:
            continue
        unique.append(item)

    nodes: list[tuple[float, float, float]] = [
        (left.level, left.value, left.shadow_price),
        *unique,
        (right.level, right.value, float("nan")),
    ]
    return tuple(
        ResponsePhase(
            target_low=nodes[i][0],
            target_high=nodes[i + 1][0],
            shadow_price=nodes[i][2],
            response_low=nodes[i][1],
            response_high=nodes[i + 1][1],
        )
        for i in range(len(nodes) - 1)
        if nodes[i + 1][0] - nodes[i][0] > x_tol
    )


def _response_limit(
    phases: tuple[ResponsePhase, ...], threshold: float
) -> ResponseLimit:
    """The first phase boundary at which the shadow price crosses below ``-threshold``."""

    if not phases:
        return ResponseLimit(
            found=False,
            message="no feasible target flux; there is no response curve to read",
            threshold=threshold,
        )
    for index, phase in enumerate(phases):
        if phase.shadow_price >= -threshold:
            continue
        if index == 0:
            return ResponseLimit(
                found=True,
                message=(
                    "the response is already declining at the lower end of the feasible "
                    "domain"
                ),
                target_flux=phase.target_low,
                response_flux=phase.response_low,
                shadow_price_before=phase.shadow_price,
                shadow_price_after=phase.shadow_price,
                threshold=threshold,
            )
        return ResponseLimit(
            found=True,
            message="shadow price crosses below -threshold at this phase boundary",
            target_flux=phase.target_low,
            response_flux=phase.response_low,
            shadow_price_before=phases[index - 1].shadow_price,
            shadow_price_after=phase.shadow_price,
            threshold=threshold,
        )
    if all(abs(phase.shadow_price) <= threshold for phase in phases):
        return ResponseLimit(
            found=False,
            message=(
                "the shadow price is zero throughout; the response is insensitive to this "
                "reaction over its whole feasible domain"
            ),
            threshold=threshold,
        )
    return ResponseLimit(
        found=False,
        message="the shadow price never falls below -threshold over the feasible domain",
        threshold=threshold,
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
    limit_threshold: float = 1e-6,
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
    product. The floor is applied before the scan range is determined, so it narrows the
    feasible domain rather than leaving guaranteed-infeasible points in the grid.

    The scan range defaults to the target's exact feasible interval under those constraints
    (FVA at 0% of the optimum). An explicit ``target_min``/``target_max`` may lie outside the
    reaction's own bounds — the scan then overrides them for the whole scan, which is
    recorded in the result metadata as ``range_outside_bounds`` because it is a deliberate
    what-if, not the model as loaded. Points inside an explicit window that the network
    cannot reach are kept with ``status`` set and NaN fluxes: the flux at which a scan stops
    being solvable is a result, not an error.

    Alongside the scanned points the result carries the *phase structure* of the response
    curve — the intervals of constant LP shadow price — and the *response limit*, the phase
    boundary at which that price first falls below ``-limit_threshold``. Both come from LP
    duals and a tangent intersection, so neither moves with ``n_steps``.
    """

    if n_steps < 2:
        raise ValueError("n_steps must be at least 2")
    if tol < 0:
        raise ValueError("tol must be non-negative")
    if limit_threshold < 0:
        raise ValueError("limit_threshold must be non-negative")
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

        # Wild-type reference: the growth optimum the scan is read against. Taken before the
        # growth floor exists, because the floor is defined as a fraction of it.
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

        # The growth floor is part of the scan's constraint set, so it is applied once, here,
        # rather than per point: the feasible domain, the grid and the phase structure must
        # all be read under the same model.
        biomass_floor = (
            None
            if biomass_fraction is None
            else biomass_fraction * wild_type["biomass_flux"]
        )
        if biomass_floor is not None:
            biomass_reaction = model.reactions.get_by_id(biomass_id)
            if biomass_floor > float(biomass_reaction.upper_bound) + tol:
                raise ValueError(
                    f"the growth floor {biomass_floor:g} exceeds the upper bound of "
                    f"{biomass_id!r}; it can never be met"
                )
            biomass_reaction.bounds = (
                max(float(biomass_reaction.lower_bound), float(biomass_floor)),
                float(biomass_reaction.upper_bound),
            )

        auto_range = target_min is None or target_max is None
        if not auto_range:
            lower, upper = float(target_min), float(target_max)  # type: ignore[arg-type]
        else:
            detected = _scan_range(model, target)
            if detected is None:
                raise ValueError(
                    f"no feasible flux interval for {target!r} under this condition and "
                    "growth floor; give target_min/target_max explicitly"
                )
            lower = detected[0] if target_min is None else float(target_min)
            upper = detected[1] if target_max is None else float(target_max)

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

        # An explicit window is a deliberate what-if that may reach past the reaction's own
        # bounds; apply it once so the feasible domain and the phases are read inside the
        # window the caller actually asked for.
        model.reactions.get_by_id(target).bounds = (lower, upper)

        # Every scan point maximizes the response under the fixed target flux.
        model.objective = model.reactions.get_by_id(response_id)
        model.objective_direction = "max"

        feasible_domain = _scan_range(model, target)

        points: list[ResponsePoint] = []
        for level in np.linspace(lower, upper, n_steps):
            with model:
                model.reactions.get_by_id(target).bounds = (float(level), float(level))
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

        phases = (
            ()
            if feasible_domain is None or feasible_domain[1] - feasible_domain[0] <= tol
            else _phases(model, target, feasible_domain[0], feasible_domain[1])
        )

    limit = _response_limit(phases, threshold=limit_threshold)
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
        limit_threshold=limit_threshold,
        feasible_min=None if feasible_domain is None else feasible_domain[0],
        feasible_max=None if feasible_domain is None else feasible_domain[1],
        n_phases=len(phases),
    )
    return FluxResponseResult(
        target=target,
        response=response_id,
        biomass=biomass_id,
        points=tuple(points),
        phases=phases,
        limit=limit,
        wild_type=wild_type,
        feasible_domain=feasible_domain,
        metadata=provenance,
    )
