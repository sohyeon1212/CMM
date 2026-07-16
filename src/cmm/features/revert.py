"""Revert-metabolism: normalization-target prediction via published MTA/rMTA.

Predicts gene/reaction knockouts that move a perturbed (e.g. disease) metabolic state back
toward a reference (e.g. healthy) state. This is the inverse of the production-design
services: instead of pushing flux toward a product, it ranks interventions by how well they
transform the source flux distribution in the direction implied by two-state expression.

``rmta`` follows Valcárcel et al. (2019) and the COBRA Toolbox ``rMTA.m`` workflow: an MTA
MIQP for the best direction, a MOMA QP, an MTA MIQP with the direction reversed for the worst
case, the published L1 transformation score, and Equation 9 with ``parameterK=100``. ``mta``
is the single published MTA MIQP. The previous all-continuous approximation remains available
only under the explicit name ``rmta_continuous`` so its results cannot be mistaken for rMTA.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

from cobra import Model

from cmm.core import solvers
from cmm.core.condition import Condition
from cmm.core.flux_state import FluxState
from cmm.core.results import TargetRanking, TargetScore
from cmm.core.provenance import run_provenance
from cmm.features._perturbation import (
    Perturbation,
    apply_perturbation,
    gene_perturbations,
    grouped_gene_perturbations,
    reaction_perturbations,
)
from cmm.omics.differential import DirectionMap

RevertMethod = Literal["rmta", "mta", "mta_miqp", "rmta_continuous"]


@dataclass(frozen=True)
class _Scores:
    moma: float
    best: float
    worst: float
    robust: float
    status: str


def _sign(x: float, tol: float = 1e-9) -> int:
    if x > tol:
        return 1
    if x < -tol:
        return -1
    return 0


def _target_reactions(model: Model, reference: FluxState) -> list[str]:
    present = {r.id for r in model.reactions}
    return [rid for rid in reference.reactions() if rid in present]


def _mta_qp(
    model: Model,
    reference: FluxState,
    direction: DirectionMap,
    target_rxns: list[str],
    alpha: float,
) -> tuple[dict[str, float], str]:
    """Solve the MOMA-based transformation QP and return (fluxes, status).

    Objective (minimize): sum_i (v_i - v_ref_i)^2  -  alpha * sum_{nonsteady} d_i (v_i - v_ref_i)
    The quadratic term is full MOMA (keeps the solve bounded/unique); the linear term rewards
    moving target reactions in their desired direction (alpha<0 rewards the wrong way).
    """

    prob = model.problem
    saved_objective = model.objective
    saved_direction = model.objective_direction
    target_set = set(target_rxns)
    try:
        expr = 0.0
        # Full MOMA quadratic over EVERY model reaction (reactions absent from the reference
        # contribute (v_i - 0)^2), matching comparison.py's reference solution and keeping the
        # QP bounded and unique. The linear transformation reward applies only to directional
        # reactions that are part of the reference target set.
        for rxn in model.reactions:
            rid = rxn.id
            vref = reference.get(rid)
            diff = rxn.flux_expression - vref
            expr = expr + diff * diff
            if alpha != 0 and rid in target_set:
                d = direction.get(rid, 0)
                if d != 0:
                    expr = expr - alpha * d * diff
        model.objective = prob.Objective(expr, direction="min")
        # Probe feasibility first: reading primals from an infeasible solve raises a backend
        # GurobiError, so bail out on a non-optimal (e.g. lethal-knockout) solve instead.
        model.slim_optimize(error_value=float("nan"))
        if model.solver.status != "optimal":
            return {}, model.solver.status
        solution = model.optimize()
        fluxes = {rid: float(v) for rid, v in solution.fluxes.items()}
        return fluxes, solution.status
    finally:
        model.objective = saved_objective
        model.objective_direction = saved_direction


def _mta_miqp(
    model: Model,
    reference: FluxState,
    direction: DirectionMap,
    target_rxns: list[str],
    alpha: float,
    eps: float,
) -> tuple[dict[str, float], str]:
    """Published MTA MIQP objective from COBRA Toolbox ``buildMTAproblemFromModel``.

    The objective is ``(1-alpha) * steady_L2 + alpha/2 * failed_changes``.  Constants are
    removed, so the implemented linear binary term is ``-alpha/2 * successful_changes``.
    """

    from optlang.symbolics import Zero

    # Work on a copy: isolates the binary vars from the caller's model/history and avoids
    # cross-knockout name collisions. The copy captures the current (KO-applied) bounds.
    work = model.copy()
    prob = work.problem

    # Phase 1: add the binary variables and FLUSH, so the constraints that reference them in
    # phase 2 do not make optlang re-queue the same variable (the duplicate-add trap).
    switches: list = []
    specs: list[tuple[str, int, float, float]] = []
    for rid in target_rxns:
        d = direction.get(rid, 0)
        if d == 0:
            continue
        rxn = work.reactions.get_by_id(rid)
        vref = reference.get(rid)
        required = vref + eps if d > 0 else vref - eps
        big_m = (
            max(1.0, required - rxn.lower_bound)
            if d > 0
            else max(1.0, rxn.upper_bound - required)
        )
        z = prob.Variable(f"mta_z_{rid}", type="binary")
        switches.append(z)
        specs.append((rid, d, vref, big_m))
    if switches:
        work.add_cons_vars(switches)
        work.solver.update()

    # Phase 2: constraints referencing the now-existing binaries.
    constraints = []
    for (rid, d, vref, big_m), z in zip(specs, switches, strict=True):
        flux = work.reactions.get_by_id(rid).flux_expression
        if d > 0:
            c = prob.Constraint(
                flux - big_m * z, lb=vref + eps - big_m, name=f"mta_up_{rid}"
            )
        else:
            c = prob.Constraint(
                flux + big_m * z, ub=vref - eps + big_m, name=f"mta_dn_{rid}"
            )
        constraints.append(c)
    if constraints:
        work.add_cons_vars(constraints)
        work.solver.update()

    # Objective: minimize weighted steady distance minus the count of transformed reactions.
    quad = Zero
    for rid in target_rxns:
        if direction.get(rid, 0) == 0:
            diff = work.reactions.get_by_id(rid).flux_expression - reference.get(rid)
            quad = quad + (1.0 - alpha) * diff * diff
    binary_reward = (alpha / 2.0) * sum(switches) if switches else Zero
    work.objective = prob.Objective(quad - binary_reward, direction="min")

    # Probe feasibility without reading primals: for an infeasible MIQP (e.g. a lethal
    # knockout) Gurobi raises when reading ``.X``, which would abort the whole revert run.
    work.slim_optimize(error_value=float("nan"))
    if work.solver.status != "optimal":
        return {}, work.solver.status
    solution = work.optimize()
    fluxes = {rid: float(v) for rid, v in solution.fluxes.items()}
    return fluxes, solution.status


def _transformation_score(
    fluxes: dict[str, float],
    reference: FluxState,
    direction: DirectionMap,
    target_rxns: list[str],
) -> float:
    """Published MTA transformation score (COBRA Toolbox ``MTA_TS.m``).

    ``TS = (successful_L1 - unsuccessful_L1) / steady_L1``.  If steady reactions do not
    move, a beneficial/adverse transformation is represented by positive/negative infinity;
    no movement at all is defined as zero rather than the undefined ``0/0``.
    """

    correct = 0.0
    wrong = 0.0
    steady_dev = 0.0
    for rid in target_rxns:
        vref = reference.get(rid)
        delta = fluxes.get(rid, 0.0) - vref
        d = direction.get(rid, 0)
        if d == 0:
            steady_dev += abs(delta)
        else:
            moved = _sign(delta)
            if moved == d:
                correct += abs(delta)
            elif moved == -d:
                wrong += abs(delta)
    numerator = correct - wrong
    if steady_dev <= 1e-12:
        if abs(numerator) <= 1e-12:
            return 0.0
        return float("inf") if numerator > 0 else float("-inf")
    return numerator / steady_dev


def _continuous_transformation_score(
    fluxes: dict[str, float],
    reference: FluxState,
    direction: DirectionMap,
    target_rxns: list[str],
) -> float:
    """Legacy continuous score retained only for ``rmta_continuous`` compatibility."""

    correct = 0.0
    wrong = 0.0
    steady_dev = 0.0
    for rid in target_rxns:
        delta = fluxes.get(rid, 0.0) - reference.get(rid)
        d = direction.get(rid, 0)
        if d == 0:
            steady_dev += delta * delta
        elif _sign(delta) == d:
            correct += abs(delta)
        elif _sign(delta) == -d:
            wrong += abs(delta)
    return (correct - wrong) / (1.0 + steady_dev)


def _prepared_direction(
    model: Model,
    reference: FluxState,
    direction: DirectionMap,
    target_rxns: list[str],
    *,
    reverse: bool = False,
) -> DirectionMap:
    """Apply the published impossible-change preprocessing and optional F/B swap."""

    values: dict[str, int] = {}
    sign = -1 if reverse else 1
    for rid in target_rxns:
        d = sign * direction.get(rid, 0)
        rxn = model.reactions.get_by_id(rid)
        # COBRA rMTA removes a requested backward move for an inactive irreversible reaction.
        if d < 0 and abs(reference.get(rid)) < 1e-6 and rxn.lower_bound >= 0:
            d = 0
        values[rid] = d
    return DirectionMap(values, metadata={"reversed": reverse})


def _score_knockout(
    model: Model,
    reference: FluxState,
    direction: DirectionMap,
    target_rxns: list[str],
    method: RevertMethod,
    alpha: float,
    epsilon: float = 1e-3,
    parameter_k: float = 100.0,
) -> _Scores:
    best_direction = _prepared_direction(model, reference, direction, target_rxns)

    if method in {"mta", "mta_miqp"}:
        flux, status = _mta_miqp(
            model, reference, best_direction, target_rxns, alpha=alpha, eps=epsilon
        )
        if status != "optimal":
            return _Scores(
                moma=float("-inf"),
                best=float("-inf"),
                worst=float("-inf"),
                robust=float("-inf"),
                status=status,
            )
        ts = _transformation_score(flux, reference, best_direction, target_rxns)
        return _Scores(moma=ts, best=ts, worst=ts, robust=ts, status="optimal")

    if method == "rmta":
        best_flux, best_status = _mta_miqp(
            model, reference, best_direction, target_rxns, alpha=alpha, eps=epsilon
        )
        moma_flux, moma_status = _mta_qp(
            model, reference, best_direction, target_rxns, 0.0
        )
        worst_direction = _prepared_direction(
            model, reference, direction, target_rxns, reverse=True
        )
        worst_flux, worst_status = _mta_miqp(
            model, reference, worst_direction, target_rxns, alpha=alpha, eps=epsilon
        )
        statuses = (best_status, moma_status, worst_status)
        if any(status != "optimal" for status in statuses):
            status = next(s for s in statuses if s != "optimal")
            return _Scores(
                moma=float("-inf"),
                best=float("-inf"),
                worst=float("-inf"),
                robust=float("-inf"),
                status=status,
            )
        best_ts = _transformation_score(
            best_flux, reference, best_direction, target_rxns
        )
        moma_ts = _transformation_score(
            moma_flux, reference, best_direction, target_rxns
        )
        # Published wTS is scored against the swapped F/B direction.
        worst_ts = _transformation_score(
            worst_flux, reference, worst_direction, target_rxns
        )
        robust = _robust_score(moma_ts, best_ts, worst_ts, parameter_k=parameter_k)
        return _Scores(
            moma=moma_ts, best=best_ts, worst=worst_ts, robust=robust, status="optimal"
        )

    # Explicit legacy continuous approximation. It is intentionally not called rMTA.
    moma_flux, status = _mta_qp(model, reference, best_direction, target_rxns, 0.0)
    if status != "optimal":
        return _Scores(
            moma=float("-inf"),
            best=float("-inf"),
            worst=float("-inf"),
            robust=float("-inf"),
            status=status,
        )
    moma_ts = _continuous_transformation_score(
        moma_flux, reference, best_direction, target_rxns
    )
    best_flux, b_status = _mta_qp(model, reference, best_direction, target_rxns, alpha)
    worst_flux, w_status = _mta_qp(
        model, reference, best_direction, target_rxns, -alpha
    )
    best_ts = (
        _continuous_transformation_score(
            best_flux, reference, best_direction, target_rxns
        )
        if b_status == "optimal"
        else moma_ts
    )
    worst_ts = (
        _continuous_transformation_score(
            worst_flux, reference, best_direction, target_rxns
        )
        if w_status == "optimal"
        else moma_ts
    )
    robust = _continuous_robust_score(moma_ts, best_ts, worst_ts)
    return _Scores(
        moma=moma_ts, best=best_ts, worst=worst_ts, robust=robust, status="optimal"
    )


def _robust_score(
    moma_ts: float,
    best_ts: float,
    worst_ts: float,
    *,
    parameter_k: float = 100.0,
) -> float:
    """Published rMTA Equation 9, matching COBRA Toolbox ``rMTA.m``."""

    if best_ts != best_ts or moma_ts != moma_ts or worst_ts != worst_ts:
        return float("-inf")
    if best_ts > 0 and moma_ts > 0 and worst_ts < 0:
        return moma_ts * parameter_k * (best_ts - worst_ts)
    return moma_ts


def _continuous_robust_score(moma_ts: float, best_ts: float, worst_ts: float) -> float:
    """Legacy additive score used only by the explicitly named continuous approximation."""

    if best_ts > 0 and worst_ts >= 0:
        return moma_ts + worst_ts
    if best_ts > 0:
        return max(0.0, moma_ts + worst_ts)
    return 0.0


def revert_targets(
    model: Model,
    source_condition: Condition | None,
    reference_state: FluxState,
    direction: DirectionMap,
    *,
    targets: Iterable[str] | None = None,
    method: RevertMethod = "rmta",
    alpha: float = 0.66,
    epsilon: float = 1e-3,
    parameter_k: float = 100.0,
    perturbation: Literal["gene", "reaction"] = "gene",
    transcript_separator: str | None = None,
) -> TargetRanking:
    """Rank knockouts by how well they revert the source state toward the target.

    Parameters mirror the design doc: ``reference_state`` is the source/disease reference
    flux distribution, ``direction`` is the per-reaction desired flux direction from
    two-state differential expression, and ``alpha`` is the transformation weight.
    """

    if method not in {"rmta", "mta", "mta_miqp", "rmta_continuous"}:
        raise ValueError(
            f"unknown method {method!r}; use 'rmta', 'mta', or 'rmta_continuous'"
        )
    if perturbation not in {"gene", "reaction"}:
        raise ValueError("perturbation must be 'gene' or 'reaction'")
    if not 0.0 <= alpha <= 1.0:
        raise ValueError("alpha must be between 0 and 1")
    if epsilon < 0:
        raise ValueError("epsilon must be non-negative")
    if parameter_k <= 0:
        raise ValueError("parameter_k must be positive")
    reference_state.validate()
    missing = [r.id for r in model.reactions if r.id not in reference_state.fluxes]
    if missing:
        shown = ", ".join(missing[:5])
        raise ValueError(
            f"reference_state is missing {len(missing)} model reactions ({shown}); "
            "provide a complete source-state flux vector"
        )

    if method == "rmta_continuous":
        solvers.require(
            "QP", model.solver.interface, feature="continuous rMTA approximation"
        )
    else:
        solvers.require("MIQP", model.solver.interface, feature="published MTA/rMTA")
    provenance = run_provenance(
        model,
        method=method,
        source_condition=source_condition.name if source_condition else None,
        alpha=alpha,
        epsilon=epsilon,
        parameter_k=parameter_k,
        perturbation=perturbation,
        transcript_separator=transcript_separator,
    )

    with model:
        if source_condition is not None:
            source_condition.apply_to(model)

        target_rxns = _target_reactions(model, reference_state)

        if perturbation == "gene":
            perts: list[Perturbation] = (
                grouped_gene_perturbations(model, transcript_separator, targets)
                if transcript_separator
                else gene_perturbations(model, targets)
            )
        else:
            perts = reaction_perturbations(model, targets)

        scored: list[TargetScore] = []
        nonoptimal = 0
        for pert in perts:
            with apply_perturbation(model, pert):
                scores = _score_knockout(
                    model,
                    reference_state,
                    direction,
                    target_rxns,
                    method,
                    alpha,
                    epsilon,
                    parameter_k,
                )
            if scores.status != "optimal":
                nonoptimal += 1
            scored.append(
                TargetScore(
                    pert.target_id,
                    scores.robust,
                    detail={
                        "bTS": scores.best,
                        "mTS": scores.moma,
                        "wTS": scores.worst,
                    },
                )
            )

    return TargetRanking(
        method=f"revert_{method}",
        targets=tuple(scored),
        metadata={
            **provenance,
            "alpha": alpha,
            "epsilon": epsilon,
            "parameter_k": parameter_k,
            "perturbation": perturbation,
            "transcript_separator": transcript_separator,
            "reference": reference_state.name,
            "reference_provenance": reference_state.provenance,
            "n_targets": len(scored),
            "n_nonoptimal": nonoptimal,
            "formulation": (
                "published_valcarcel_2019"
                if method == "rmta"
                else "published_yizhak_2013"
                if method in {"mta", "mta_miqp"}
                else "continuous_heuristic"
            ),
        },
    ).sorted()
