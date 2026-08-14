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
        normalized = {str(k): int(v) for k, v in self.directions.items()}
        invalid = {
            key: value for key, value in normalized.items() if value not in {-1, 0, 1}
        }
        if invalid:
            raise ValueError(f"direction values must be -1, 0, or 1; got {invalid}")
        object.__setattr__(self, "directions", normalized)
        object.__setattr__(self, "metadata", dict(self.metadata))

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

    if pseudocount < 0 or not math.isfinite(pseudocount):
        raise ValueError("pseudocount must be finite and non-negative")
    genes = set(source) & set(target)
    out: dict[str, float] = {}
    for g in genes:
        source_value = float(source[g])
        target_value = float(target[g])
        if (
            not math.isfinite(source_value)
            or not math.isfinite(target_value)
            or source_value < 0
            or target_value < 0
        ):
            raise ValueError(
                f"expression for gene {g!r} must be finite and non-negative"
            )
        numerator = target_value + pseudocount
        denominator = source_value + pseudocount
        if denominator == 0:
            out[g] = 0.0 if numerator == 0 else math.inf
        elif numerator == 0:
            out[g] = -math.inf
        else:
            out[g] = math.log2(numerator / denominator)
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

    if up_threshold < 0 or down_threshold < 0:
        raise ValueError("fold-change thresholds must be non-negative")
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


def _eval_gpr_indicator(node: ast.AST | None, indicator: Mapping[str, int]) -> int:
    """Evaluate one **binary** 0/1 gene indicator over a GPR AST: AND -> min, OR -> max.

    On a 0/1 vector, ``min`` at AND is "all of the genes fire" and ``max`` at OR is "at
    least one of them fires" — the two set-membership tests Yizhak et al. (2013) state.
    Genes absent from ``indicator`` contribute 0.
    """

    if node is None:
        return 0
    if isinstance(node, ast.Module):
        # cobra's GPR subclasses ast.Module but stores the expression node directly in
        # `.body` rather than a statement list.
        body = getattr(node, "body", None)
        if isinstance(body, list):
            body = body[0] if body else None
        return _eval_gpr_indicator(body, indicator)
    if isinstance(node, ast.Expression):
        return _eval_gpr_indicator(node.body, indicator)
    if isinstance(node, ast.Name):
        return 1 if indicator.get(node.id, 0) else 0
    if isinstance(node, ast.BoolOp):
        values = [_eval_gpr_indicator(v, indicator) for v in node.values]
        if not values:
            return 0
        return min(values) if isinstance(node.op, ast.And) else max(values)
    return 0


def _eval_gpr_direction(node: ast.AST | None, gene_dirs: Mapping[str, int]) -> int:
    """Evaluate a reaction's expression direction over its GPR AST (Yizhak et al. 2013).

    The paper's rule, verbatim: a reaction is elevated/reduced *"(a) if it is catalysed by a
    complex of enzymes (an 'and' logical relation) and **all** of the genes encoding them
    were categorized as elevated or reduced, respectively … and (b) if it is catalysed by
    isoenzymes (an 'or' logical relation) and **at least one** of them was categorized as
    elevated or reduced, respectively. If a subset is categorized as elevated and another
    subset as reduced, the reaction is considered unchanged. In addition, in any other cases
    not specifically described here, the reaction is considered unchanged."*

    That is implemented, as in the rMTA reference code (COBRA Toolbox
    ``diffexprs2rxnFBS.m``), by evaluating min-at-AND / max-at-OR **twice over binary
    indicators** and subtracting: ``plus`` fires iff the elevated set satisfies the rule,
    ``minus`` iff the reduced set does, and ``plus - minus`` sends the mixed case to 0.

    A signed ``min``/``max`` over ``{-1, 0, +1}`` is *not* equivalent: it diverges in 4 of
    the 7 two-gene label cases with a directional bias (``AND(+1, -1)`` returns ``-1`` where
    the paper says 0; ``OR(-1, 0)`` returns 0 where the paper says ``-1``). CMM used the
    signed form before 0.4.0.
    """

    plus = _eval_gpr_indicator(node, {g: 1 for g, d in gene_dirs.items() if int(d) > 0})
    minus = _eval_gpr_indicator(
        node, {g: 1 for g, d in gene_dirs.items() if int(d) < 0}
    )
    return plus - minus


def reaction_directions(
    model: Model,
    gene_dirs: Mapping[str, int],
    *,
    reference: FluxState | None = None,
    reactions: Iterable[str] | None = None,
    flux_tol: float = 1e-9,
) -> DirectionMap:
    """Map gene directions to per-reaction desired flux-value directions.

    The GPR is resolved by Yizhak et al.'s (2013) rule — AND requires *all* associated genes
    to share a direction, OR requires *at least one*, and conflicting evidence yields
    *unchanged* — via the two-pass binary decomposition in ``_eval_gpr_direction``. Reactions
    with no GPR are 0, matching the reference code's ``rxnFBS(grRules == '') = 0``.

    When a ``reference`` (source) flux state is given, the expression direction is combined
    with the sign of the reference flux so the label is in flux-*value* space: a reaction
    operating in reverse with an up-regulated enzyme should carry *more negative* flux, i.e.
    decrease. That step is Yizhak's, and the reference code's ``rxnFBS(Vref < 0) = -rxnFBS``.
    Reactions inactive in the source (|v_ref| <= tol) are assumed to operate forward, so an
    up-regulated enzyme means "turn on" (increase) — **that convention is CMM's own**; neither
    Yizhak et al. (2013) nor Valcárcel et al. (2019) specify it.
    """

    rxn_ids = (
        list(reactions) if reactions is not None else [r.id for r in model.reactions]
    )
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
        metadata={
            "has_reference": reference is not None,
            "gpr_rule": "yizhak2013_two_pass_binary",
        },
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
