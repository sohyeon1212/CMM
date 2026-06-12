"""Find knockouts that transform one flux state (condition A) toward another (condition B).

Given a source flux distribution (e.g. condition A, possibly omics-predicted) and a target
distribution (condition B), rank candidate gene/reaction knockouts by how close they bring
the flux to the target. Two methods:

- ``moma``: for each knockout, the cell minimally adjusts from A (MOMA); rank by how much that
  adjusted state moves toward B (reduction in distance-to-B).
- ``mta``: derive a per-reaction desired direction from B − A and run robust MTA (rMTA).

This generalizes revert-metabolism from expression-derived directions to two explicit flux
states, so omics-predicted condition fluxes can be used directly as source/target.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Literal

from cobra import Model

from cmm.core import solvers
from cmm.core.flux_state import FluxState
from cmm.core.results import TargetRanking
from cmm.features._perturbation import (
    apply_perturbation,
    gene_perturbations,
    reaction_perturbations,
)
from cmm.features.comparison import moma
from cmm.features.revert import revert_targets
from cmm.omics.differential import DirectionMap

TransformationMethod = Literal["moma", "mta"]


def direction_from_states(
    source: FluxState,
    target: FluxState,
    *,
    reactions: Iterable[str] | None = None,
    tol: float = 1e-6,
) -> DirectionMap:
    """Per-reaction desired flux-value direction to go from source to target (B − A)."""

    keys = set(reactions) if reactions is not None else set(source.fluxes) | set(target.fluxes)
    directions: dict[str, int] = {}
    for rid in keys:
        delta = target.get(rid) - source.get(rid)
        directions[rid] = 1 if delta > tol else (-1 if delta < -tol else 0)
    return DirectionMap(directions=directions, metadata={"from": source.name, "to": target.name})


def transformation_targets(
    model: Model,
    source_state: FluxState,
    target_state: FluxState,
    *,
    method: TransformationMethod = "moma",
    perturbation: Literal["gene", "reaction"] = "gene",
    targets: Iterable[str] | None = None,
    order: Literal[1, 2] = 2,
    alpha: float = 0.9,
) -> TargetRanking:
    """Rank knockouts that move the source flux state toward the target flux state."""

    if method == "mta":
        direction = direction_from_states(source_state, target_state)
        ranking = revert_targets(
            model, None, source_state, direction,
            targets=targets, method="rmta", alpha=alpha, perturbation=perturbation,
        )
        return TargetRanking(
            method="transform_mta", targets=ranking.targets, metadata=ranking.metadata
        )

    if method != "moma":
        raise ValueError(f"unknown method {method!r}; use 'moma' or 'mta'")

    solvers.require("QP", model.solver.interface, feature="MOMA transformation targets")
    compare_rxns = [r.id for r in model.reactions]
    baseline = source_state.distance(target_state, reactions=compare_rxns, order=order)

    if perturbation == "gene":
        perts = gene_perturbations(model, targets)
    else:
        perts = reaction_perturbations(model, targets)

    scored: list[tuple[str, float]] = []
    for pert in perts:
        with apply_perturbation(model, pert):
            result = moma(model, source_state, linear=False)
        if result.status != "optimal":
            scored.append((pert.target_id, 0.0))
            continue
        predicted = FluxState(result.fluxes, name="perturbed")
        distance_to_target = predicted.distance(
            target_state, reactions=compare_rxns, order=order
        )
        scored.append((pert.target_id, baseline - distance_to_target))

    return TargetRanking.from_scores(
        method="transform_moma",
        scores=scored,
        metadata={
            "perturbation": perturbation,
            "baseline_distance": baseline,
            "source": source_state.name,
            "target": target_state.name,
        },
    )
