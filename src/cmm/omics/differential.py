"""Two-state differential expression -> per-reaction desired flux direction.

Revert-metabolism (and the MTA family generally) needs more than single-state expression
constraints: it needs to know, for each reaction, whether flux should *increase*,
*decrease*, or *stay* when moving from a source state (e.g. disease) to a target state
(e.g. healthy). This module derives that per-reaction direction label from gene expression
in the two states, mapped through the GPR and combined with the source reference flux sign.

Direction codes: +1 increase (forward), -1 decrease (backward), 0 steady/ambiguous.
"""

from __future__ import annotations

import ast
import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

from cobra import Model

from cmm.core.flux_state import FluxState


@dataclass(frozen=True)
class DirectionMap:
    """Per-reaction desired flux-change direction (+1 / -1 / 0)."""

    directions: Mapping[str, int]
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "directions", {str(k): int(v) for k, v in self.directions.items()}
        )

    def __getitem__(self, reaction_id: str) -> int:
        return self.directions[reaction_id]

    def __len__(self) -> int:
        return len(self.directions)

    def get(self, reaction_id: str, default: int = 0) -> int:
        return self.directions.get(reaction_id, default)

    def items(self):
        return self.directions.items()

    def forward(self) -> frozenset[str]:
        return frozenset(r for r, d in self.directions.items() if d > 0)

    def backward(self) -> frozenset[str]:
        return frozenset(r for r, d in self.directions.items() if d < 0)

    def steady(self) -> frozenset[str]:
        return frozenset(r for r, d in self.directions.items() if d == 0)

    def nonsteady(self) -> frozenset[str]:
        return self.forward() | self.backward()


def gene_log2_fold_change(
    source: Mapping[str, float],
    target: Mapping[str, float],
    pseudocount: float = 1.0,
) -> dict[str, float]:
    """log2((target + pseudo) / (source + pseudo)) over genes present in both states."""

    genes = set(source) & set(target)
    out: dict[str, float] = {}
    for g in genes:
        out[g] = math.log2((target[g] + pseudocount) / (source[g] + pseudocount))
    return out


def gene_directions(
    source: Mapping[str, float],
    target: Mapping[str, float],
    *,
    up_threshold: float = 1.0,
    down_threshold: float = 1.0,
    pseudocount: float = 1.0,
) -> dict[str, int]:
    """Discretize gene log2 fold change into +1 (up in target) / -1 (down) / 0."""

    lfc = gene_log2_fold_change(source, target, pseudocount=pseudocount)
    out: dict[str, int] = {}
    for g, value in lfc.items():
        if value >= up_threshold:
            out[g] = 1
        elif value <= -down_threshold:
            out[g] = -1
        else:
            out[g] = 0
    return out


def _eval_gpr_direction(node: ast.AST | None, gene_dirs: Mapping[str, int]) -> int:
    """Evaluate a reaction's expression direction over its GPR AST.

    OR (isozymes) takes the max direction (any branch can carry flux); AND (complex) takes
    the min (the limiting subunit). Genes absent from ``gene_dirs`` contribute 0.
    """

    if node is None:
        return 0
    if isinstance(node, ast.Module):
        # cobra's GPR subclasses ast.Module but stores the expression node directly in
        # `.body` rather than a statement list.
        body = getattr(node, "body", None)
        if isinstance(body, list):
            body = body[0] if body else None
        return _eval_gpr_direction(body, gene_dirs)
    if isinstance(node, ast.Expression):
        return _eval_gpr_direction(node.body, gene_dirs)
    if isinstance(node, ast.Name):
        return int(gene_dirs.get(node.id, 0))
    if isinstance(node, ast.BoolOp):
        values = [_eval_gpr_direction(v, gene_dirs) for v in node.values]
        return min(values) if isinstance(node.op, ast.And) else max(values)
    return 0


def reaction_directions(
    model: Model,
    gene_dirs: Mapping[str, int],
    *,
    reference: FluxState | None = None,
    reactions: Iterable[str] | None = None,
    flux_tol: float = 1e-9,
) -> DirectionMap:
    """Map gene directions to per-reaction desired flux-value directions.

    When a ``reference`` (source) flux state is given, the expression direction is combined
    with the sign of the reference flux so the label is in flux-*value* space: a reaction
    operating in reverse with an up-regulated enzyme should carry *more negative* flux, i.e.
    decrease. Reactions inactive in the source (|v_ref| <= tol) are assumed to operate
    forward, so an up-regulated enzyme means "turn on" (increase).
    """

    rxn_ids = list(reactions) if reactions is not None else [r.id for r in model.reactions]
    directions: dict[str, int] = {}
    for rid in rxn_ids:
        rxn = model.reactions.get_by_id(rid)
        if not rxn.genes:
            directions[rid] = 0
            continue
        expr_dir = _eval_gpr_direction(rxn.gpr, gene_dirs)
        if reference is not None:
            v = reference.get(rid)
            if abs(v) > flux_tol:
                expr_dir = expr_dir * (1 if v > 0 else -1)
        directions[rid] = int(expr_dir)
    return DirectionMap(
        directions=directions,
        metadata={"has_reference": reference is not None},
    )


def differential_expression(
    model: Model,
    source: Mapping[str, float],
    target: Mapping[str, float],
    *,
    reference: FluxState | None = None,
    up_threshold: float = 1.0,
    down_threshold: float = 1.0,
    pseudocount: float = 1.0,
    reactions: Iterable[str] | None = None,
) -> DirectionMap:
    """Convenience: two-state gene expression -> reaction `DirectionMap` in one call."""

    g_dirs = gene_directions(
        source,
        target,
        up_threshold=up_threshold,
        down_threshold=down_threshold,
        pseudocount=pseudocount,
    )
    return reaction_directions(model, g_dirs, reference=reference, reactions=reactions)
