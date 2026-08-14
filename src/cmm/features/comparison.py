"""Reference-distance perturbation responses: MOMA and ROOM.

Given a reference flux distribution (a :class:`~cmm.core.flux_state.FluxState`, typically
the wild-type/source state) and a model that has been perturbed (e.g. a knockout applied),
these services predict the perturbed flux state as the one closest to the reference under a
chosen metric:

- **MOMA, L2** — minimization of metabolic adjustment, Segrè D, Vitkup D & Church GM (2002),
  *PNAS* 99(23):15112-15117. The QP minimizes ``Sum_i d_i^2``; the paper's Eq. (4) defines
  the distance itself as ``D = sqrt(Sum_i d_i^2)``, which is what
  :attr:`ComparisonResult.distance` reports.
- **MOMA, L1** — the linear (LP) variant as implemented in COBRApy's ``moma(linear=True)``.
  Segrè et al. define no L1 variant; this is the ecosystem's linear relaxation of the same
  idea and is attributed to COBRApy, not to the 2002 paper.
- **ROOM** — regulatory on/off minimization, Shlomi T, Berkman O & Ruppin E (2005), *PNAS*
  102(21):7695-7700. The MILP minimizes the *number* of significantly changed reactions,
  which is a count and not a distance at all.

Because the three objectives are three different quantities, the result keeps the raw solver
objective under :attr:`ComparisonResult.objective_value` and exposes a
:attr:`ComparisonResult.distance` that means one thing, labelled by
:attr:`ComparisonResult.distance_kind`. ROOM's switch count has its own field.

The numerics delegate to cobra's maintained implementations; this layer adapts a
solver-neutral `FluxState` reference (any provenance) into the cobra `Solution` they expect
and returns a plain, serializable result.

**Two properties inherited from COBRApy, disclosed rather than changed.**

1. ROOM as implemented by COBRApy adds a constraint capping the perturbed objective at the
   reference objective (``room_old_objective``). It is not in Shlomi et al.'s Eqs 1-3 and is
   confirmed binding on ``e_coli_core`` (perturbed growth 0.3822 against a reference 0.43696).
2. MOMA and ROOM can return a different flux vector from a different COBRApy version or
   solver run — 12 of 95 and 54 of 95 reactions respectively on ``e_coli_core``. These are
   alternate optima, not disagreements: recomputing each objective from the returned vector
   reproduces the reported optimum (to 4e-14 for MOMA's QP, and exactly 30.0 = 30.0 for
   ROOM's switch count).
"""

from __future__ import annotations

import math
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
DistanceKind = Literal["euclidean_l2", "l1", "none"]
RoomUseCase = Literal["flux_prediction", "lethality"]

#: Shlomi et al. (2005) report two tolerance pairs and use them for two different questions:
#: ``delta=0.03 / epsilon=1e-3`` for **flux prediction** and ``delta=0.1 / epsilon=1e-2`` for
#: **lethality prediction**. The choice is not cosmetic — on ``e_coli_core`` the two pairs
#: give 531 against 401 total switches over the same 35 gene knockouts, a 24% shift in the
#: ranking score — so the use case selects the pair rather than one pair serving both.
ROOM_TOLERANCES: dict[str, tuple[float, float]] = {
    "flux_prediction": (0.03, 1e-3),
    "lethality": (0.1, 1e-2),
}


def _room_tolerances(
    use_case: RoomUseCase, delta: float | None, epsilon: float | None
) -> tuple[float, float]:
    """Resolve ROOM's (delta, epsilon), with explicit values overriding the preset."""

    if use_case not in ROOM_TOLERANCES:
        raise ValueError(
            f"unknown ROOM use case {use_case!r}; use "
            f"{' or '.join(repr(k) for k in ROOM_TOLERANCES)}"
        )
    preset_delta, preset_epsilon = ROOM_TOLERANCES[use_case]
    resolved_delta = preset_delta if delta is None else float(delta)
    resolved_epsilon = preset_epsilon if epsilon is None else float(epsilon)
    if resolved_delta < 0 or resolved_epsilon < 0:
        raise ValueError("ROOM delta and epsilon must be non-negative")
    return resolved_delta, resolved_epsilon


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

    **The default is ``pfba``, which is a deliberate, measured deviation from Segrè et al.
    (2002), who specify an FBA reference.** FBA's optimum is a face, not a point, so the
    reference is not unique in practice: eight identical ``reference_flux(..., "fba")`` calls
    on *i*JO1366 returned three distinct flux vectors. pFBA's minimal-total-flux tie-break
    makes the reference reproducible. The cost of the deviation was measured and is small —
    on *i*JO1366 the FBA-vs-pFBA reference changes MOMA growth by at most 0.00096 h^-1
    (0.098% of wild type) and reclassifies 0 of 24 knockouts; on ``e_coli_core`` it changes
    nothing at all. Pass ``method="fba"`` for Segrè's own protocol.
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
    """Outcome of a reference-distance solve.

    ``objective_value`` is the raw solver objective and therefore means something different
    for each method: ``Sum d^2`` for L2 MOMA, ``Sum |d|`` for L1 MOMA, a reaction count for
    ROOM. ``distance`` is the one quantity that is comparable across runs of the *same*
    method and is always the distance ``distance_kind`` names — the Euclidean distance of
    Segrè et al. Eq. (4) for ``moma_l2`` (i.e. the square root of the QP objective), the L1
    distance for ``moma_l1``, and ``None`` for ROOM, whose objective is a count reported by
    ``n_changed_reactions``. A ``NaN`` distance means the perturbed model had no solution.
    """

    method: ComparisonMethod
    status: str
    objective_value: float
    distance: float | None
    distance_kind: DistanceKind
    n_changed_reactions: float | None = None
    fluxes: dict[str, float] = field(default_factory=dict)
    metadata: dict[str, object] = field(default_factory=dict)

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


def _reference_metadata(reference: FluxState) -> dict[str, object]:
    """Which reference state a distance was measured against — it changes the answer."""

    return {
        "reference": reference.name,
        "reference_provenance": reference.provenance,
        "reference_method": reference.metadata.get("source_method"),
    }


def moma(
    model: Model, reference: FluxState, *, linear: bool = False
) -> ComparisonResult:
    """Minimal-adjustment flux state nearest the reference (Segrè et al. 2002).

    ``linear=False`` solves the published L2 variant (QP, requires a QP-capable solver) and
    reports the Euclidean distance of Eq. (4). ``linear=True`` solves the linear variant as
    implemented in COBRApy (LP, runs on any solver) and reports the L1 distance; the 2002
    paper defines no such variant.
    """

    method: ComparisonMethod = "moma_l1" if linear else "moma_l2"
    metadata = _reference_metadata(reference)
    if not linear:
        solvers.require("QP", model.solver.interface, feature="L2 MOMA")
    # COBRApy/Gurobi may leak a backend ``GurobiError`` while reading primal values from an
    # infeasible solve.  Probe feasibility without retrieving a full solution first.  Pass a
    # NaN error_value so a lethal knockout returns non-optimal here instead of raising
    # ``Infeasible`` (which would abort a whole batch run over many knockouts).
    model.slim_optimize(error_value=float("nan"))
    if model.solver.status != "optimal":
        return _infeasible(method, metadata)

    reference_solution = _reference_solution(model, reference)
    try:
        solution = _cobra_moma(model, solution=reference_solution, linear=linear)
    except (
        OptimizationError
    ):  # cobra raises for an infeasible perturbed model (e.g. lethal KO)
        return _infeasible(method, metadata)
    return _result(method, solution, metadata)


def room(
    model: Model,
    reference: FluxState,
    *,
    linear: bool = False,
    use_case: RoomUseCase = "flux_prediction",
    delta: float | None = None,
    epsilon: float | None = None,
) -> ComparisonResult:
    """Minimize the number of reactions that significantly deviate from the reference.

    The default exact formulation is a MILP (binary on/off switches). ``linear=True`` uses
    the LP relaxation, whose objective is a continuous lower bound on the switch count rather
    than an integer.

    ``use_case`` selects Shlomi et al.'s published tolerance pair: ``"flux_prediction"``
    (delta 0.03, epsilon 1e-3) or ``"lethality"`` (delta 0.1, epsilon 1e-2). Passing ``delta``
    or ``epsilon`` explicitly overrides the preset. See :data:`ROOM_TOLERANCES`.

    The solve is delegated to COBRApy, which additionally caps the perturbed objective at the
    reference objective — a constraint that is not in Shlomi et al. and is confirmed binding
    on ``e_coli_core``. See the module docstring.
    """

    resolved_delta, resolved_epsilon = _room_tolerances(use_case, delta, epsilon)
    metadata: dict[str, object] = {
        **_reference_metadata(reference),
        "room_use_case": use_case,
        "delta": resolved_delta,
        "epsilon": resolved_epsilon,
    }
    if not linear:
        solvers.require("MILP", model.solver.interface, feature="ROOM")
    # NaN error_value so a lethal knockout returns non-optimal instead of raising Infeasible
    # (which would abort a batch run over many knockouts).
    model.slim_optimize(error_value=float("nan"))
    if model.solver.status != "optimal":
        return _infeasible("room", metadata)

    reference_solution = _reference_solution(model, reference)
    try:
        solution = _cobra_room(
            model,
            solution=reference_solution,
            linear=linear,
            delta=resolved_delta,
            epsilon=resolved_epsilon,
        )
    except OptimizationError:  # infeasible perturbed model
        return _infeasible("room", metadata)
    return _result("room", solution, metadata)


def _quantities(
    method: ComparisonMethod, objective_value: float
) -> tuple[float | None, DistanceKind, float | None]:
    """Split a raw solver objective into (distance, distance_kind, switch count)."""

    if method == "moma_l2":
        # Segrè et al. Eq. (4): D = sqrt(Sum d^2). The QP objective is the sum itself, so a
        # value reported without this square root is off by the square.
        distance = (
            float("nan")
            if math.isnan(objective_value)
            else math.sqrt(max(objective_value, 0.0))
        )
        return distance, "euclidean_l2", None
    if method == "moma_l1":
        return objective_value, "l1", None
    # ROOM's objective is a count of significantly changed reactions, not a distance.
    return None, "none", objective_value


def _infeasible(
    method: ComparisonMethod, metadata: dict[str, object]
) -> ComparisonResult:
    distance, kind, changed = _quantities(method, float("nan"))
    return ComparisonResult(
        method=method,
        status="infeasible",
        objective_value=float("nan"),
        distance=distance,
        distance_kind=kind,
        n_changed_reactions=changed,
        metadata=metadata,
    )


def _result(
    method: ComparisonMethod, solution: Solution, metadata: dict[str, object]
) -> ComparisonResult:
    objective_value = (
        float(solution.objective_value)
        if solution.objective_value is not None
        else float("nan")
    )
    distance, kind, changed = _quantities(method, objective_value)
    fluxes = {rid: float(v) for rid, v in solution.fluxes.items()}
    return ComparisonResult(
        method=method,
        status=solution.status,
        objective_value=objective_value,
        distance=distance,
        distance_kind=kind,
        n_changed_reactions=changed,
        fluxes=fluxes,
        metadata=metadata,
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
    room_use_case: RoomUseCase = "flux_prediction",
    delta: float | None = None,
    epsilon: float | None = None,
) -> ComparisonResult:
    """Run MOMA/ROOM against ``reference`` with ``reaction_ids`` forced to zero.

    Accepts one or many reactions, so it covers a single-reaction knockout, a multi-reaction
    knockout, and a gene knockout (pass the reactions the gene disables via its GPR). The model
    is restored on return. Method keys: ``moma_l2`` (QP), ``moma_l1`` (LP), ``room`` (MILP).

    ``room_use_case`` defaults to ``"flux_prediction"`` because this entry point predicts the
    flux state of one named perturbation; the batch screen in :func:`batch_comparison`
    defaults to ``"lethality"`` instead.
    """

    if method not in {"moma_l1", "moma_l2", "room"}:
        raise ValueError("method must be 'moma_l1', 'moma_l2', or 'room'")
    with model:
        for rid in reaction_ids:
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        if method == "room":
            return room(
                model,
                reference,
                use_case=room_use_case,
                delta=delta,
                epsilon=epsilon,
            )
        return moma(model, reference, linear=method == "moma_l1")


@dataclass(frozen=True)
class BatchComparisonRow:
    """One knockout's MOMA/ROOM outcome in a batch run.

    The three quantities the old single ``distance`` column conflated are kept apart, with
    the same meanings as on :class:`ComparisonResult`.
    """

    target_id: str
    kind: str  # "gene" | "reaction"
    status: str
    objective_value: float  # raw solver objective, method-dependent
    distance: float | None  # Euclidean for moma_l2, L1 for moma_l1, None for room
    distance_kind: DistanceKind
    n_changed_reactions: float | None  # ROOM switch count; None for MOMA
    objective: float  # objective-reaction flux at the perturbed state (e.g. growth)
    n_reactions: int  # how many reactions the knockout blocked
    product_flux: float = float(
        "nan"
    )  # target-product flux at the perturbed state (if any)


def batch_comparison(
    model: Model,
    reference: FluxState,
    perturbations: Sequence[Perturbation],
    *,
    method: ComparisonMethod = "moma_l2",
    room_use_case: RoomUseCase = "lethality",
    delta: float | None = None,
    epsilon: float | None = None,
    objective_reaction: str | None = None,
    product_reaction: str | None = None,
) -> list[BatchComparisonRow]:
    """Run MOMA/ROOM for each perturbation against one reference, returning a batch table.

    This is the batch perturbation-response job runner: enumerate gene or reaction knockouts
    (via :mod:`cmm.features._perturbation`), score each against the same wild-type reference,
    and collect the distance, the objective (growth) at the perturbed state, and the status.
    When ``product_reaction`` is given, each row also carries that reaction's perturbed flux.

    A knockout screen asks which deletions are lethal, so ``room_use_case`` defaults to
    Shlomi et al.'s **lethality** tolerance pair (delta 0.1, epsilon 1e-2) rather than their
    flux-prediction pair. Pass ``room_use_case="flux_prediction"`` to score the screen as a
    set of flux predictions instead; the two disagree by roughly a quarter of the switch
    count. The parameter has no effect on the MOMA methods.

    Lethal knockouts stay in the table as ``status="infeasible"`` rows with NaN fluxes — the
    essentiality result is the point of the screen, not an error to filter out.
    """

    objective = objective_reaction or _objective_reaction_id(model)
    rows: list[BatchComparisonRow] = []
    for pert in perturbations:
        result = knockout_comparison(
            model,
            reference,
            pert.reaction_ids,
            method=method,
            room_use_case=room_use_case,
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
                objective_value=result.objective_value,
                distance=result.distance,
                distance_kind=result.distance_kind,
                n_changed_reactions=result.n_changed_reactions,
                objective=growth,
                n_reactions=len(pert.reaction_ids),
                product_flux=product_flux,
            )
        )
    return rows
