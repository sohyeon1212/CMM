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

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import Literal

import pandas as pd
from cobra import Model
from cobra.core import Solution
from cobra.exceptions import OptimizationError
from cobra.flux_analysis import moma as _cobra_moma
from cobra.flux_analysis import room as _cobra_room

from cmm.core import solvers
from cmm.core.flux_state import FluxState, Provenance
from cmm.core.simulation import fba as _fba
from cmm.core.simulation import pfba as _pfba
from cmm.features._perturbation import Perturbation

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
        simulation_result = _fba(model)
        if simulation_result.status != "optimal":
            raise ValueError(f"FBA reference solve is {simulation_result.status}")
        fluxes = simulation_result.fluxes
        provenance: Provenance = "fba"
        metadata = simulation_result.metadata
    elif method == "pfba":
        simulation_result = _pfba(model)
        if simulation_result.status != "optimal":
            raise ValueError(f"pFBA reference solve is {simulation_result.status}")
        fluxes = simulation_result.fluxes
        provenance = "pfba"
        metadata = simulation_result.metadata
    elif method in ("lad", "eflux2"):
        if gene_expression is None:
            raise ValueError(f"method {method!r} requires gene_expression")
        from cmm.omics.expression import integrate_expression

        omics_result = integrate_expression(model, gene_expression, method=method)
        if omics_result.status != "optimal" or not omics_result.fluxes:
            raise ValueError(f"{method} reference solve is {omics_result.status}")
        fluxes = omics_result.fluxes
        provenance = "imported"
        metadata = omics_result.metadata
    else:
        raise ValueError(f"unknown reference method {method!r}")
    return FluxState(
        fluxes,
        name=name or f"{method}_reference",
        provenance=provenance,
        metadata={"source_method": method, **metadata},
    )


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

    missing = [
        reaction.id
        for reaction in model.reactions
        if reaction.id not in reference.fluxes
    ]
    if missing:
        raise ValueError(
            f"reference flux is missing {len(missing)} model reactions: {', '.join(missing[:5])}"
        )
    fluxes = pd.Series(
        {r.id: reference.get(r.id) for r in model.reactions}, dtype=float
    )
    objective_value = float(
        sum(r.objective_coefficient * reference.get(r.id) for r in model.reactions)
    )
    return Solution(objective_value=objective_value, status="optimal", fluxes=fluxes)


def moma(
    model: Model, reference: FluxState, *, linear: bool = False
) -> ComparisonResult:
    """Minimal-adjustment flux state nearest the reference.

    ``linear=True`` solves the L1 variant (LP, runs on any solver). ``linear=False`` solves
    the L2 variant (QP, requires a QP-capable solver).
    """

    method: ComparisonMethod = "moma_l1" if linear else "moma_l2"
    if not linear:
        solvers.require("QP", model.solver.interface, feature="L2 MOMA")
    # COBRApy/Gurobi may leak a backend ``GurobiError`` while reading primal values from an
    # infeasible solve.  Probe feasibility without retrieving a full solution first.  Pass a
    # NaN error_value so a lethal knockout returns non-optimal here instead of raising
    # ``Infeasible`` (which would abort a whole batch run over many knockouts).
    model.slim_optimize(error_value=float("nan"))
    if model.solver.status != "optimal":
        return ComparisonResult(
            method=method, status="infeasible", distance=float("nan")
        )

    reference_solution = _reference_solution(model, reference)
    try:
        solution = _cobra_moma(model, solution=reference_solution, linear=linear)
    except (
        OptimizationError
    ):  # cobra raises for an infeasible perturbed model (e.g. lethal KO)
        return ComparisonResult(
            method=method, status="infeasible", distance=float("nan")
        )
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

    if delta < 0 or epsilon < 0:
        raise ValueError("ROOM delta and epsilon must be non-negative")
    if not linear:
        solvers.require("MILP", model.solver.interface, feature="ROOM")
    # NaN error_value so a lethal knockout returns non-optimal instead of raising Infeasible
    # (which would abort a batch run over many knockouts).
    model.slim_optimize(error_value=float("nan"))
    if model.solver.status != "optimal":
        return ComparisonResult(
            method="room", status="infeasible", distance=float("nan")
        )

    reference_solution = _reference_solution(model, reference)
    try:
        solution = _cobra_room(
            model,
            solution=reference_solution,
            linear=linear,
            delta=delta,
            epsilon=epsilon,
        )
    except OptimizationError:  # infeasible perturbed model
        return ComparisonResult(
            method="room", status="infeasible", distance=float("nan")
        )
    return _result("room", solution)


def _result(method: ComparisonMethod, solution: Solution) -> ComparisonResult:
    distance = (
        float(solution.objective_value)
        if solution.objective_value is not None
        else float("nan")
    )
    fluxes = {rid: float(v) for rid, v in solution.fluxes.items()}
    return ComparisonResult(
        method=method, status=solution.status, distance=distance, fluxes=fluxes
    )


def _objective_reaction_id(model: Model) -> str | None:
    for rxn in model.reactions:
        if rxn.objective_coefficient != 0:
            return rxn.id
    return None


def knockout_comparison(
    model: Model,
    reference: FluxState,
    reaction_ids: Iterable[str],
    *,
    method: ComparisonMethod = "moma_l2",
    delta: float = 0.03,
    epsilon: float = 1e-3,
) -> ComparisonResult:
    """Run MOMA/ROOM against ``reference`` with ``reaction_ids`` forced to zero.

    Accepts one or many reactions, so it covers a single-reaction knockout, a multi-reaction
    knockout, and a gene knockout (pass the reactions the gene disables via its GPR). The model
    is restored on return. Method keys: ``moma_l2`` (QP), ``moma_l1`` (LP), ``room`` (MILP).
    """

    if method not in {"moma_l1", "moma_l2", "room"}:
        raise ValueError("method must be 'moma_l1', 'moma_l2', or 'room'")
    with model:
        for rid in reaction_ids:
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        if method == "room":
            return room(model, reference, delta=delta, epsilon=epsilon)
        return moma(model, reference, linear=method == "moma_l1")


@dataclass(frozen=True)
class BatchComparisonRow:
    """One knockout's MOMA/ROOM outcome in a batch run."""

    target_id: str
    kind: str  # "gene" | "reaction"
    status: str
    distance: float
    objective: float  # objective-reaction flux at the perturbed state (e.g. growth)
    n_reactions: int  # how many reactions the knockout blocked
    product_flux: float = float("nan")  # target-product flux at the perturbed state (if any)


def batch_comparison(
    model: Model,
    reference: FluxState,
    perturbations: Sequence[Perturbation],
    *,
    method: ComparisonMethod = "moma_l2",
    delta: float = 0.03,
    epsilon: float = 1e-3,
    objective_reaction: str | None = None,
    product_reaction: str | None = None,
) -> list[BatchComparisonRow]:
    """Run MOMA/ROOM for each perturbation against one reference, returning a batch table.

    This is the batch perturbation-response job runner: enumerate gene or reaction knockouts
    (via :mod:`cmm.features._perturbation`), score each against the same wild-type reference,
    and collect the distance, the objective (growth) at the perturbed state, and the status.
    When ``product_reaction`` is given, each row also carries that reaction's perturbed flux.
    """

    objective = objective_reaction or _objective_reaction_id(model)
    rows: list[BatchComparisonRow] = []
    for pert in perturbations:
        result = knockout_comparison(
            model,
            reference,
            pert.reaction_ids,
            method=method,
            delta=delta,
            epsilon=epsilon,
        )
        growth = (
            float(result.fluxes.get(objective, float("nan")))
            if (result.fluxes and objective is not None)
            else float("nan")
        )
        product_flux = (
            float(result.fluxes.get(product_reaction, float("nan")))
            if (result.fluxes and product_reaction is not None)
            else float("nan")
        )
        rows.append(
            BatchComparisonRow(
                target_id=pert.target_id,
                kind=pert.kind,
                status=result.status,
                distance=result.distance,
                objective=growth,
                n_reactions=len(pert.reaction_ids),
                product_flux=product_flux,
            )
        )
    return rows
