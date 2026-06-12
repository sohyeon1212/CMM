"""Reference-distance perturbation responses: MOMA and ROOM.

Given a reference flux distribution (a :class:`~cmm.core.flux_state.FluxState`, typically
the wild-type/source state) and a model that has been perturbed (e.g. a knockout applied),
these services predict the perturbed flux state as the one closest to the reference under a
chosen metric:

- MOMA (Segre 2002): minimal metabolic adjustment. L2 = quadratic (QP); L1 = linear (LP).
- ROOM (Shlomi 2005): minimal number of significantly changed reactions (MILP).

The numerics delegate to cobra's maintained implementations; this layer adapts a
solver-neutral `FluxState` reference (any provenance) into the cobra `Solution` they expect
and returns a plain, serializable result.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import pandas as pd
from cobra import Model
from cobra.core import Solution
from cobra.flux_analysis import moma as _cobra_moma
from cobra.flux_analysis import room as _cobra_room

from cmm.core import solvers
from cmm.core.flux_state import FluxState
from cmm.core.simulation import fba as _fba
from cmm.core.simulation import pfba as _pfba

ComparisonMethod = Literal["moma_l1", "moma_l2", "room"]
ReferenceMethod = Literal["fba", "pfba", "lad", "eflux2"]


def reference_flux(
    model: Model,
    method: ReferenceMethod = "pfba",
    *,
    gene_expression: dict[str, float] | None = None,
    name: str | None = None,
) -> FluxState:
    """Build a reference/template flux state for MOMA/ROOM from a chosen method.

    ``fba``/``pfba`` need only the model; ``lad``/``eflux2`` additionally require single-state
    ``gene_expression`` (the omics-predicted distribution becomes the template).
    """

    if method == "fba":
        fluxes = _fba(model).fluxes
    elif method == "pfba":
        fluxes = _pfba(model).fluxes
    elif method in ("lad", "eflux2"):
        if gene_expression is None:
            raise ValueError(f"method {method!r} requires gene_expression")
        from cmm.omics.expression import integrate_expression

        fluxes = integrate_expression(model, gene_expression, method=method).fluxes
    else:
        raise ValueError(f"unknown reference method {method!r}")
    return FluxState(fluxes, name=name or f"{method}_reference", provenance="imported")


@dataclass(frozen=True)
class ComparisonResult:
    """Outcome of a reference-distance solve."""

    method: ComparisonMethod
    status: str
    distance: float
    fluxes: dict[str, float] = field(default_factory=dict)

    def to_flux_state(self, name: str = "perturbed") -> FluxState:
        return FluxState(self.fluxes, name=name, provenance="imported")


def _reference_solution(model: Model, reference: FluxState) -> Solution:
    """Build a cobra Solution covering every model reaction from the reference state.

    The objective value must be the model objective evaluated at the reference fluxes, not a
    placeholder: ROOM pins the original objective to ``solution.objective_value`` via its
    ``room_old_objective`` constraint, so a zero here would force the whole network to zero
    flux.
    """

    fluxes = pd.Series(
        {r.id: reference.get(r.id) for r in model.reactions}, dtype=float
    )
    objective_value = float(
        sum(r.objective_coefficient * reference.get(r.id) for r in model.reactions)
    )
    return Solution(objective_value=objective_value, status="optimal", fluxes=fluxes)


def moma(model: Model, reference: FluxState, *, linear: bool = False) -> ComparisonResult:
    """Minimal-adjustment flux state nearest the reference.

    ``linear=True`` solves the L1 variant (LP, runs on any solver). ``linear=False`` solves
    the L2 variant (QP, requires a QP-capable solver).
    """

    method: ComparisonMethod = "moma_l1" if linear else "moma_l2"
    if not linear:
        solvers.require("QP", model.solver.interface, feature="L2 MOMA")

    reference_solution = _reference_solution(model, reference)
    try:
        solution = _cobra_moma(model, solution=reference_solution, linear=linear)
    except Exception:  # cobra raises when the perturbed model is infeasible (e.g. lethal KO)
        return ComparisonResult(method=method, status="infeasible", distance=float("nan"))
    return _result(method, solution)


def room(
    model: Model,
    reference: FluxState,
    *,
    linear: bool = False,
    delta: float = 0.03,
    epsilon: float = 1e-3,
) -> ComparisonResult:
    """Minimize the number of reactions that significantly deviate from the reference.

    The default exact formulation is a MILP (binary on/off switches). ``linear=True`` uses
    the LP relaxation.
    """

    if not linear:
        solvers.require("MILP", model.solver.interface, feature="ROOM")

    reference_solution = _reference_solution(model, reference)
    try:
        solution = _cobra_room(
            model, solution=reference_solution, linear=linear, delta=delta, epsilon=epsilon
        )
    except Exception:  # infeasible perturbed model
        return ComparisonResult(method="room", status="infeasible", distance=float("nan"))
    return _result("room", solution)


def _result(method: ComparisonMethod, solution: Solution) -> ComparisonResult:
    distance = (
        float(solution.objective_value) if solution.objective_value is not None else float("nan")
    )
    fluxes = {rid: float(v) for rid, v in solution.fluxes.items()}
    return ComparisonResult(method=method, status=solution.status, distance=distance, fluxes=fluxes)
