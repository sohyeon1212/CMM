"""Multi-condition omics: predict per-condition flux and compare as log-change.

Given an expression table with several condition columns (one column per condition), predict
a flux distribution per condition with E-Flux2 or LAD, then compare conditions as log2
fold-change of flux magnitude.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass, field

import pandas as pd
from cobra import Model

from cmm.core.provenance import run_provenance
from cmm.omics.expression import OmicsFluxResult, integrate_expression


@dataclass(frozen=True)
class ConditionFluxes:
    """Per-condition predicted flux distributions from one expression table.

    ``metadata`` is the run provenance of the multi-condition job as a whole: one
    :func:`~cmm.core.provenance.run_provenance` block for the model every condition was solved
    on, plus the condition names and the integration method. Each condition's own solve keeps
    its own block under ``results[condition].metadata``, so a single condition's numbers stay
    traceable when they are lifted out of the container. Save ``metadata`` alongside any CSV
    of the fluxes.
    """

    method: str
    results: dict[str, OmicsFluxResult]
    metadata: dict[str, object] = field(default_factory=dict)

    def conditions(self) -> tuple[str, ...]:
        return tuple(self.results.keys())

    def fluxes(self, condition: str) -> dict[str, float]:
        result = self.results[condition]
        if result.status != "optimal" or not result.fluxes:
            raise ValueError(
                f"condition {condition!r} has no valid flux state ({result.status})"
            )
        return result.fluxes


def read_expression_table(path: str, gene_column: str | None = None) -> pd.DataFrame:
    """Read a gene-by-condition expression table; index = gene id, columns = conditions."""

    frame = pd.read_csv(path, sep=None, engine="python")
    if frame.empty or frame.shape[1] < 2:
        raise ValueError(
            "expression table must contain a gene column and at least one condition"
        )
    if frame.columns.duplicated().any():
        raise ValueError("expression table contains duplicate condition columns")
    gene_column = gene_column or frame.columns[0]
    frame = frame.set_index(gene_column)
    if frame.index.duplicated().any():
        raise ValueError("expression table contains duplicate gene identifiers")
    numeric = frame.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    finite_or_missing = numeric.apply(
        lambda column: column.map(
            lambda value: pd.isna(value) or math.isfinite(float(value))
        )
    )
    if not finite_or_missing.to_numpy().all() or (numeric < 0).any().any():
        raise ValueError("expression values must be finite and non-negative")
    return numeric


def predict_condition_fluxes(
    model: Model,
    expression: pd.DataFrame,
    *,
    method: str = "eflux2",
    conditions: Iterable[str] | None = None,
    **kwargs,
) -> ConditionFluxes:
    """Predict a flux distribution for each condition column with E-Flux2 or LAD.

    Conditions are solved independently — one ``integrate_expression`` call per column, each
    inside its own ``with model:`` block — so no bound or objective mutation leaks between
    them and adding or removing a column cannot change another column's fluxes. No
    cross-condition normalisation is applied either, by ``read_expression_table`` or here.

    **Consequence to keep in mind when interpreting the result.** E-Flux2's normalisation
    denominator is computed *within* each condition, so a condition in which every gene is
    uniformly k-fold higher produces identical bounds and therefore an identical flux
    distribution: **a uniform global expression shift produces zero predicted flux change**,
    and only within-condition relative expression influences the prediction. (Kim et al. 2016
    are not blind to it — their bound is the absolute expression value.) Normalise the
    expression table for library size before integration, and do not draw conclusions from a
    global shift in expression level.
    """

    if expression.empty:
        raise ValueError("expression table is empty")
    columns = list(conditions) if conditions is not None else list(expression.columns)
    missing = [
        condition for condition in columns if condition not in expression.columns
    ]
    if missing:
        raise KeyError(f"unknown condition columns: {missing}")
    # Fingerprint the model before the per-condition solves: each runs inside its own
    # ``with model:`` block, so this is the one model state every condition was solved on.
    provenance = run_provenance(
        model,
        method="predict_condition_fluxes",
        integration_method=method,
        conditions=tuple(str(condition) for condition in columns),
        n_genes=int(expression.shape[0]),
        **kwargs,
    )
    results: dict[str, OmicsFluxResult] = {}
    for condition in columns:
        gene_expression = {
            str(gene): float(value)
            for gene, value in expression[condition].items()
            if pd.notna(value)
        }
        results[condition] = integrate_expression(
            model, gene_expression, method=method, **kwargs
        )
    return ConditionFluxes(
        method=method,
        results=results,
        metadata={
            **provenance,
            "integration_method": method,
            "conditions": tuple(str(condition) for condition in columns),
            "n_conditions": len(columns),
            "n_nonoptimal": sum(
                1 for result in results.values() if result.status != "optimal"
            ),
        },
    )


def flux_log_change(
    source: dict[str, float],
    target: dict[str, float],
    *,
    reactions: Iterable[str] | None = None,
    pseudocount: float = 1e-3,
) -> dict[str, float]:
    """log2 fold-change of flux *magnitude* between two conditions (target vs source).

    **A CMM utility, not an implementation of a published method.** A pseudocounted log2
    ratio is generic practice; there is no source paper for this function and it must not be
    cited to one, nor described as a differential-flux *method*.

    Uses ``log2((|v_target| + pseudo) / (|v_source| + pseudo))`` so zero/near-zero fluxes are
    handled gracefully; the pseudocount bounds the change for reactions that switch on/off.

    Two artefacts of that definition, both to report alongside the numbers:

    - **A direction reversal reads as "no change".** The ratio is taken on magnitudes, so a
      reaction going ``+5 -> -5`` returns ``0.0``. Use ``sign_flips`` to detect reversals;
      this function alone will not show them.
    - **The on/off magnitude is set by the pseudocount, not by the data.** For ``0 -> v`` the
      result is ``log2(1 + v / pseudo)``, so at the default ``pseudo = 1e-3`` a reaction
      switching on at ``v = 10`` reports about 13.3. Changing the pseudocount changes every
      on/off value while leaving the ranking intact, so the pseudocount must be reported.
    """

    if pseudocount < 0 or not math.isfinite(pseudocount):
        raise ValueError("pseudocount must be finite and non-negative")
    keys = set(reactions) if reactions is not None else set(source) | set(target)
    out: dict[str, float] = {}
    for rid in keys:
        source_value = float(source.get(rid, 0.0))
        target_value = float(target.get(rid, 0.0))
        if not math.isfinite(source_value) or not math.isfinite(target_value):
            raise ValueError(f"flux values for reaction {rid!r} must be finite")
        a = abs(source_value) + pseudocount
        b = abs(target_value) + pseudocount
        # The zero branches are only reachable when pseudocount=0 (the default 1e-3 bounds
        # every change). Guard both off-states symmetrically: log2(0/0)=0, log2(b/0)=+inf
        # (switch on), log2(0/a)=-inf (switch off) — the last would raise a math domain error.
        if a == 0.0:
            out[rid] = 0.0 if b == 0.0 else math.inf
        elif b == 0.0:
            out[rid] = -math.inf
        else:
            out[rid] = math.log2(b / a)
    return out


def sign_flips(
    source: dict[str, float],
    target: dict[str, float],
    *,
    reactions: Iterable[str] | None = None,
    tol: float = 1e-6,
) -> list[str]:
    """Reactions whose flux direction reverses between the two conditions."""

    keys = set(reactions) if reactions is not None else set(source) | set(target)
    flipped = []
    for rid in keys:
        a, b = source.get(rid, 0.0), target.get(rid, 0.0)
        if abs(a) > tol and abs(b) > tol and (a > 0) != (b > 0):
            flipped.append(rid)
    return sorted(flipped)
