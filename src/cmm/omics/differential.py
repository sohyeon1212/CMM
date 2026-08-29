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
from typing import TYPE_CHECKING

from cobra import Model

from cmm.core.flux_state import FluxState

if TYPE_CHECKING:  # pandas stays out of import time
    import pandas as pd


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


def gene_directions_from_replicates(
    source: "pd.DataFrame",
    target: "pd.DataFrame",
    *,
    p_value_cutoff: float = 0.05,
) -> "pd.DataFrame":
    """Discretize replicate expression by Student's t-test, as Yizhak et al. (2013) specify.

    Their step 2 applies "the Student's t-test to the source and target gene expression
    measurements" at P < 0.05 and splits genes three ways. :func:`gene_directions` cuts on
    fold change instead, which is all a single measurement per gene supports; this is the
    published route and needs replicates.

    Both frames are indexed by gene with one column per replicate. Values are treated as
    already log-scaled — as microarray and RNA-seq summaries usually are — so the difference
    of means is the log2 fold change.

    Returns the evidence alongside the label, because a direction map on its own cannot be
    audited: the report has to be able to say which genes drove it and how strongly.
    """

    import numpy as np
    import pandas as pd
    from scipy import stats

    if not 0.0 < p_value_cutoff <= 1.0:
        raise ValueError("p_value_cutoff must be in (0, 1]")
    shared = source.index.intersection(target.index)
    if not len(shared):
        raise ValueError(
            "source and target expression share no gene ids; check that both use the same "
            "identifier system"
        )
    if source.shape[1] < 2 or target.shape[1] < 2:
        raise ValueError(
            "the t-test needs at least two replicates on each side; this pair has "
            f"{source.shape[1]} and {target.shape[1]}. Use gene_directions() for a "
            "fold-change cut, and state that the published test was not applied."
        )
    src, tgt = source.loc[shared], target.loc[shared]
    statistic, p_value = stats.ttest_ind(
        tgt.to_numpy(), src.to_numpy(), axis=1, equal_var=True
    )
    frame = pd.DataFrame(
        {
            "log2_fold_change": tgt.mean(axis=1) - src.mean(axis=1),
            "t_statistic": statistic,
            "p_value": p_value,
        },
        index=shared,
    )
    frame["significant"] = frame["p_value"] < p_value_cutoff
    # +1 means elevated in the target, so the source reaction's activity must rise.
    frame["direction"] = np.where(
        ~frame["significant"], 0, np.where(frame["log2_fold_change"] > 0, 1, -1)
    ).astype(int)
    return frame


def gene_directions_by_fold_change(
    source: "pd.DataFrame",
    target: "pd.DataFrame",
    *,
    up_threshold: float = 1.0,
    down_threshold: float = 1.0,
) -> "pd.DataFrame":
    """Discretize log-scale expression on fold change, in the replicate-frame evidence shape.

    The fold-change counterpart of :func:`gene_directions_from_replicates`, for data that
    cannot support the published t-test — one measurement per state, or replicates a caller
    has chosen not to test. Both frames are indexed by gene with one column per replicate and
    are treated as **already log-scaled**, so the difference of column means is the log2 fold
    change; that is the same convention the t-test route uses, which is what lets a caller
    swap between the two without also changing what the numbers mean.

    ``p_value`` is present but NaN so the returned frame has the same columns either way, and
    a report can state plainly that no test was run rather than omitting the column.
    """

    import numpy as np
    import pandas as pd

    if up_threshold < 0 or down_threshold < 0:
        raise ValueError("fold-change thresholds must be non-negative")
    shared = source.index.intersection(target.index)
    if not len(shared):
        raise ValueError(
            "source and target expression share no gene ids; check that both use the same "
            "identifier system"
        )
    src, tgt = source.loc[shared], target.loc[shared]
    frame = pd.DataFrame(
        {"log2_fold_change": tgt.mean(axis=1) - src.mean(axis=1)}, index=shared
    )
    frame["t_statistic"] = np.nan
    frame["p_value"] = np.nan
    frame["significant"] = (frame["log2_fold_change"] >= up_threshold) | (
        frame["log2_fold_change"] <= -down_threshold
    )
    frame["direction"] = np.where(
        ~frame["significant"], 0, np.where(frame["log2_fold_change"] > 0, 1, -1)
    ).astype(int)
    return frame


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


def _reaction_change_strength(
    model: Model, reaction_ids: Iterable[str], gene_log2fc: Mapping[str, float]
) -> dict[str, float]:
    """Strongest absolute gene-level log2 fold change among a reaction's GPR genes.

    Used only to order the changed reactions before the top-N cut, never to decide a
    reaction's direction — that stays with the GPR rule in :func:`_eval_gpr_direction`.
    """

    strength: dict[str, float] = {}
    for rid in reaction_ids:
        values = [
            abs(gene_log2fc[g.id])
            for g in model.reactions.get_by_id(rid).genes
            if g.id in gene_log2fc and math.isfinite(gene_log2fc[g.id])
        ]
        strength[rid] = max(values) if values else 0.0
    return strength


def restrict_to_top_changed(
    model: Model,
    direction: DirectionMap,
    gene_log2fc: Mapping[str, float],
    top_n: int,
    *,
    gene_p_values: Mapping[str, float] | None = None,
) -> DirectionMap:
    """Keep only the ``top_n`` most differentially expressed changed reactions.

    Yizhak et al. (2013), Methods: within the reactions that pass the significance cut,
    *"the top 100-200 most differentially expressed reactions are defined as the set of
    'changed' reactions"*. Everything dropped becomes steady (0), which moves it from the
    MIQP's binary success count into its quadratic stay-put term.

    Two reasons the paper gives for the cut: that many changed reactions already suffice to
    recover the correct perturbation, and each one adds a binary variable, so the cut is what
    keeps the MIQP tractable.

    The paper does not say which statistic orders the reactions. Without replicates the only
    evidence available is fold change, which is the default. Pass ``gene_p_values`` when a
    t-test was run: the smallest P value among a reaction's GPR genes is the closer analogue
    of the paper's gene-level test, with |log2FC| breaking ties. Whichever was used is
    recorded in the returned map's metadata.

    Ties in strength are broken on reaction id so the retained set is reproducible.
    """

    if top_n < 0:
        raise ValueError("top_n must be non-negative")
    changed = sorted(direction.nonsteady())
    if len(changed) <= top_n:
        return DirectionMap(
            directions=dict(direction.directions),
            metadata={
                **dict(direction.metadata),
                "top_n_changed": top_n,
                "n_changed_before_cut": len(changed),
                "n_changed_kept": len(changed),
            },
        )
    strength = _reaction_change_strength(model, changed, gene_log2fc)
    if gene_p_values is None:
        ranking = "largest |log2FC| among the reaction's genes"
        ordered = sorted(changed, key=lambda rid: (-strength[rid], rid))
    else:
        ranking = "smallest gene p-value, then largest |log2FC|"

        def evidence(reaction_id: str) -> tuple[float, float, str]:
            gene_ids = [g.id for g in model.reactions.get_by_id(reaction_id).genes]
            best = min(
                (gene_p_values[g] for g in gene_ids if g in gene_p_values), default=1.0
            )
            return (best, -strength[reaction_id], reaction_id)

        ordered = sorted(changed, key=evidence)
    kept = set(ordered[:top_n])
    directions = {
        rid: (value if rid in kept else 0)
        for rid, value in direction.directions.items()
    }
    return DirectionMap(
        directions=directions,
        metadata={
            **dict(direction.metadata),
            "top_n_changed": top_n,
            "n_changed_before_cut": len(changed),
            "n_changed_kept": len(kept),
            "changed_ranking": ranking,
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
    top_n_changed: int | None = None,
) -> DirectionMap:
    """Convenience: two-state gene expression -> reaction `DirectionMap` in one call.

    ``top_n_changed`` applies Yizhak et al.'s cut to the changed set; see
    :func:`restrict_to_top_changed`. Left at ``None`` every reaction that passes the
    fold-change thresholds stays changed, which is CMM's historical behaviour.
    """

    g_dirs = gene_directions(
        source,
        target,
        up_threshold=up_threshold,
        down_threshold=down_threshold,
        pseudocount=pseudocount,
    )
    direction = reaction_directions(
        model, g_dirs, reference=reference, reactions=reactions
    )
    if top_n_changed is None:
        return direction
    lfc = gene_log2_fold_change(source, target, pseudocount=pseudocount)
    return restrict_to_top_changed(model, direction, lfc, top_n_changed)
