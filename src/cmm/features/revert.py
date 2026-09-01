"""Revert-metabolism: normalization-target prediction via published MTA/rMTA.

Predicts gene/reaction knockouts that move a perturbed (e.g. disease) metabolic state back
toward a reference (e.g. healthy) state. This is the inverse of the production-design
services: instead of pushing flux toward a product, it ranks interventions by how well they
transform the source flux distribution in the direction implied by two-state expression.

``rmta`` follows Valcárcel et al. (2019) and the COBRA Toolbox ``rMTA.m`` workflow: an MTA
MIQP for the best direction, a MOMA QP, an MTA MIQP with the direction reversed for the worst
case, the published L1 transformation score, and Equation 9 with ``parameterK=100``. ``mta``
is the single published MTA MIQP (Yizhak et al. 2013). The previous all-continuous
approximation remains available only under the explicit name ``rmta_continuous`` so its
results cannot be mistaken for rMTA.

**Deviations from the published pipelines, disclosed rather than changed.**

1. ``epsilon`` is a fixed scalar (:data:`DEFAULT_EPSILON`). Yizhak et al. and Valcárcel et al.
   derive the required change per reaction from a *sampled* reference distribution, or use a
   per-reaction epsilon with 1e-3 only as a floor. Deriving it requires the flux sampling CMM
   deliberately avoids for determinism, so the scalar is kept and stated.
2. For ``rmta``, the impossible-change mask is applied to the *already reversed* direction map
   when building the worst-case problem, whereas the COBRA Toolbox masks in preprocessing
   (``diffexprs2rxnFBS.m``) and reverses the masked map. The two orders disagree for an
   inactive irreversible reaction: a raw ``+1`` becomes ``0`` here and ``-1`` there, and a raw
   ``-1`` becomes ``+1`` here and ``0`` there. Masking after the swap is the physically
   consistent reading — a backward move on an inactive irreversible reaction is impossible in
   either frame — but it is not the reference implementation's order, and it changes ``wTS``
   and therefore Equation 9's branch. Unverified against ``rMTA.m`` itself.
3. The source state CMM's callers usually supply is a deterministic E-Flux2 solve at
   ``objective_fraction=1.0``, i.e. contextualization is replaced by a single growth-maximal
   flux vector rather than the contextualization-plus-sampling of the papers, and the
   diseased/source state is forced to be growth-maximal — an assumption Yizhak et al.
   explicitly argue against for non-proliferating cells. ``reference_state`` is a plain
   :class:`~cmm.core.flux_state.FluxState` argument, so passing an externally sampled state
   restores the published pipeline exactly; that escape hatch is the supported route.
4. The published transformation score divides by the steady-set L1 deviation, which is zero
   for a large fraction of candidates on small models (38 of 69 on ``e_coli_core``) and makes
   the score ``±inf``. The denominator is floored at ``epsilon`` — the method's own threshold
   for a change that counts — so the score stays finite and orderable. See
   :func:`_transformation_score`.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Sequence
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
    perturbation_provenance,
    reaction_perturbations,
    target_display_name,
)
from cmm.omics.differential import DirectionMap

RevertMethod = Literal["rmta", "mta", "mta_miqp", "rmta_continuous"]

#: Flux change a changed reaction must clear, measured from the source reference, before the
#: move counts as a success. Yizhak et al. derive this per data set from a sampled reference
#: distribution; CMM keeps a fixed scalar (see the module docstring, deviation 1). Callers who
#: know their flux scale should set it explicitly.
DEFAULT_EPSILON = 1e-3


def tie_structure(scored: Sequence[TargetScore]) -> dict[str, object]:
    """How much of a ranking is actually ordered, as opposed to alphabetical.

    :meth:`TargetRanking.sorted` breaks equal scores on ``target_id``, so a top-k taken
    inside a tie block is a slice of the alphabet. Reporting the block sizes is what lets a
    reader tell the two apart; it is cheap and belongs in every ranking's provenance.
    """

    if not scored:
        return {
            "n_distinct_scores": 0,
            "largest_tie_block": 0,
            "score_resolution": 0.0,
        }
    counts = Counter(round(target.score, 9) for target in scored)
    return {
        "n_distinct_scores": len(counts),
        "largest_tie_block": max(counts.values()),
        "score_resolution": len(counts) / len(scored),
    }


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
    *,
    steady_floor: float = 1e-3,
) -> float:
    """Published MTA transformation score (COBRA Toolbox ``MTA_TS.m``), made finite.

    ``TS = (successful_L1 - unsuccessful_L1) / steady_L1``. The denominator is floored at
    ``steady_floor``, which callers set to the run's own ``epsilon``: a steady-set deviation
    smaller than the change the method itself calls significant cannot be resolved from zero,
    and dividing by it produced ``±inf``.

    That mattered. On ``e_coli_core`` with the SC-02 condition pair the steady deviation is
    *exactly* zero for 38 of 69 solvable gene knockouts, so 38 candidates shared the single
    score ``+inf`` and :meth:`TargetRanking.sorted` broke the tie on ``target_id`` — the
    reported "top 38" was an alphabetical slice. The floor keeps the published ratio wherever
    the denominator is meaningful (the smallest non-zero steady deviation measured on that
    run is 3.93, four thousand times the floor) and orders the degenerate block by the amount
    of correct movement instead of by gene name. Ties that remain are real ties in the
    optimum, and the ranking's provenance reports how many there are.

    No movement at all is still exactly zero rather than the undefined ``0/0``.
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
    if abs(numerator) <= 1e-12 and steady_dev <= 1e-12:
        return 0.0
    return numerator / max(steady_dev, steady_floor, 1e-12)


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
    """Apply the published impossible-change preprocessing and optional F/B swap.

    The incoming ``direction`` carries its own provenance — notably ``gpr_rule``, which says
    how the two-state expression was resolved through the GPR and therefore what the
    direction sets mean. It is carried onto the derived map rather than dropped, so a
    reversed (worst-case) run is still traceable to the rule that produced it.
    """

    values: dict[str, int] = {}
    sign = -1 if reverse else 1
    masked: list[str] = []
    for rid in target_rxns:
        d = sign * direction.get(rid, 0)
        rxn = model.reactions.get_by_id(rid)
        # COBRA rMTA removes a requested backward move for an inactive irreversible reaction.
        if d < 0 and abs(reference.get(rid)) < 1e-6 and rxn.lower_bound >= 0:
            d = 0
            masked.append(rid)
        values[rid] = d
    return DirectionMap(
        values,
        metadata={
            **dict(direction.metadata),
            "reversed": reverse,
            "n_impossible_masked": len(masked),
        },
    )


def _score_knockout(
    model: Model,
    reference: FluxState,
    best_direction: DirectionMap,
    worst_direction: DirectionMap,
    target_rxns: list[str],
    method: RevertMethod,
    alpha: float,
    epsilon: float = DEFAULT_EPSILON,
    parameter_k: float = 100.0,
) -> _Scores:
    """Score one candidate. Both direction maps are prepared once by the caller.

    They are deliberately *not* rebuilt here. The impossible-change mask depends on reaction
    bounds, and this function runs inside the candidate's knockout context, so rebuilding
    would evaluate the mask against knockout-modified bounds and give each candidate its own
    steady set — and with it its own transformation-score denominator. Yizhak et al. mask once
    against the unperturbed model, which is what the caller now does.
    """

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
        ts = _transformation_score(
            flux, reference, best_direction, target_rxns, steady_floor=epsilon
        )
        return _Scores(moma=ts, best=ts, worst=ts, robust=ts, status="optimal")

    if method == "rmta":
        best_flux, best_status = _mta_miqp(
            model, reference, best_direction, target_rxns, alpha=alpha, eps=epsilon
        )
        moma_flux, moma_status = _mta_qp(
            model, reference, best_direction, target_rxns, 0.0
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
            best_flux, reference, best_direction, target_rxns, steady_floor=epsilon
        )
        moma_ts = _transformation_score(
            moma_flux, reference, best_direction, target_rxns, steady_floor=epsilon
        )
        # Published wTS is scored against the swapped F/B direction.
        worst_ts = _transformation_score(
            worst_flux, reference, worst_direction, target_rxns, steady_floor=epsilon
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
    epsilon: float = DEFAULT_EPSILON,
    parameter_k: float = 100.0,
    perturbation: Literal["gene", "reaction"] = "gene",
    transcript_separator: str | None = None,
) -> TargetRanking:
    """Rank knockouts by how well they revert the source state toward the target.

    Parameters mirror the design doc: ``reference_state`` is the source/disease reference
    flux distribution, ``direction`` is the per-reaction desired flux direction from
    two-state differential expression, and ``alpha`` is the transformation weight.

    ``reference_state`` accepts any :class:`~cmm.core.flux_state.FluxState`, so the published
    contextualization-plus-sampling source state can be supplied from outside CMM; the module
    docstring lists that and the other disclosed deviations.

    The returned ranking's metadata carries the direction map's own provenance (including
    ``gpr_rule``), the enumeration counts — an inert gene that blocks no reaction is not
    scored, and how many were skipped is recorded rather than left invisible — and the tie
    structure of the scores, because ties are broken on ``target_id`` and a top-k slice taken
    inside a tie block is alphabetical rather than meaningful.
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
            perts: Sequence[Perturbation] = (
                grouped_gene_perturbations(model, transcript_separator, targets)
                if transcript_separator
                else gene_perturbations(model, targets)
            )
        else:
            perts = reaction_perturbations(model, targets)

        # A ``rmta_continuous`` run is not rMTA. The method string and the ``formulation``
        # tag say so, but a bare ``to_frame().to_csv()`` carries neither — so the marker also
        # rides on every row as its own column.
        marker: dict[str, float] = (
            {"continuous_heuristic_not_rmta": 1.0}
            if method == "rmta_continuous"
            else {}
        )

        # Prepare both direction maps once, against the unperturbed model. Yizhak et al. apply
        # the impossible-change mask in preprocessing, before any knockout is simulated; doing
        # it per candidate would let the steady set — the transformation score's denominator —
        # differ between candidates whose scores are then ranked against each other.
        best_direction = _prepared_direction(
            model, reference_state, direction, target_rxns
        )
        worst_direction = _prepared_direction(
            model, reference_state, direction, target_rxns, reverse=True
        )

        scored: list[TargetScore] = []
        nonoptimal = 0
        for pert in perts:
            with apply_perturbation(model, pert):
                scores = _score_knockout(
                    model,
                    reference_state,
                    best_direction,
                    worst_direction,
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
                        **marker,
                    },
                    target_name=target_display_name(model, pert.target_id, pert.kind),
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
            "direction_provenance": dict(direction.metadata),
            "n_directional_reactions": len(direction.nonsteady()),
            # Now a single number for the whole run rather than a per-candidate quantity,
            # because the mask is applied once against the unperturbed model. The reversed
            # count is reported only for rmta, the one method that solves the worst-case MIQP.
            "n_impossible_masked": best_direction.metadata.get("n_impossible_masked"),
            **(
                {
                    "n_impossible_masked_reversed": worst_direction.metadata.get(
                        "n_impossible_masked"
                    )
                }
                if method == "rmta"
                else {}
            ),
            **perturbation_provenance(perts),
            **tie_structure(scored),
        },
    ).sorted()
