# Design: Revert Metabolism (Normalization-Target Prediction)

## Purpose

Predict gene/reaction interventions that move a perturbed metabolic state (for example a
disease, stress, or off-spec production state) back toward a reference "normal" state.

This is the inverse of the production-optimization workflows already in scope (FSEOF,
OptKnock, RobustKnock), which push metabolism *away* from wild type toward a product.
Revert-metabolism instead asks: which knockout makes the diseased flux distribution look
most like the healthy one? No legacy feature covers this, so it is a net-new
capability rather than a reimplementation.

The method is the Metabolic Transformation Algorithm (MTA, Yizhak et al. 2013) and its
robust variant rMTA (Valcárcel et al. 2019). rMTA is the recommended default because it
uses a continuous QP formulation that scales to genome-scale all-gene knockout screens and
does not require a MIQP solver.

## Inputs

- `model`: a cobra model.
- `source_condition` (`core.Condition`): constraints describing the perturbed/source state.
- `reference_state` (`core.FluxState`): the source-state reference flux vector `v_ref`,
  produced by pFBA on the source condition or by the mean of a sampling run. rMTA is
  sensitive to `v_ref`; sampling-mean is preferred for genome-scale models.
- `direction` (`omics.DirectionMap`): per-reaction desired change derived from two-state
  differential expression (target vs source). Each reaction is labeled:
  - `+1` forward / increase  (reaction should carry more flux in the normal state)
  - `-1` backward / decrease
  - `0`  steady / unchanged
- `targets` (optional): the gene or reaction knockouts to score. Defaults to all genes.
- `params`: `alpha` (transformation weight), `epsilon` (numerical floor), and the
  perturbation mode (`gene` or `reaction`).

`reference_state` comes from Phase 3 and `direction` from Phase 5, so this feature is pure
composition over primitives that already exist by the time it lands.

## Algorithm (rMTA, continuous)

For an unperturbed reference flux `v_ref` and a knockout that fixes a reaction set to zero,
solve a MOMA-style quadratic program that trades off two terms:

1. **Steady fidelity** — keep `steady` reactions (`direction == 0`) close to `v_ref`:
   minimize `sum_{i in steady} (v_i - v_ref_i)^2`.
2. **Desired transformation** — reward movement of `forward`/`backward` reactions in the
   labeled direction, weighted by `alpha`:
   for `forward` reward `(v_i - v_ref_i)`, for `backward` reward `(v_ref_i - v_i)`.

rMTA runs this QP three times per knockout with different `alpha` weightings to obtain:

- `bTS` — best-case transformation score (`alpha`: reward desired direction)
- `mTS` — MOMA baseline score (`alpha = 0`: pure minimal adjustment)
- `wTS` — worst-case score (`-alpha`: reward the wrong direction, adversarial)

The per-knockout transformation score `TS(v)` is the magnitude-weighted net agreement
between the achieved flux change `(v - v_ref)` and the `direction` labels, normalized by the
steady disturbance: `TS = (correct_movement - wrong_movement) / (1 + steady_deviation)`.

They combine into a robust score `rTS` (implemented in `_robust_score`):

```
rTS = mTS + wTS               if bTS > 0 and wTS >= 0   (transforms even adversarially)
rTS = max(0, mTS + wTS)       if bTS > 0                (helps best-case but reversible)
rTS = 0                       if bTS <= 0               (cannot transform at all)
```

A knockout ranks high only when it transforms in the best case (`bTS > 0`) **and** its
worst-case behaviour is still beneficial; the score is anchored on the guaranteed
(worst-case) transformation plus the neutral MOMA signal. Note: this differs from the
textbook `mTS·(bTS - wTS)` gate, which requires `wTS < 0` and therefore zeros out knockouts
that *force* a correct transformation (e.g. a knockout whose only feasible reroute is the
healthy branch). The additive worst-case form is the one validated to rank such forced-
reversion targets correctly, so it is the implemented default.

The original MTA (a MIQP that maximizes the count of successfully transformed reactions) is
offered as an optional high-fidelity mode behind the MILP/MIQP solver gate. The two modes
share inputs, the perturbation runner, and the result type; only the per-KO solve differs.

## Outputs

`core.results.TargetRanking`: an ordered list of `(target_id, score, detail)` rows with
`score = rTS` (or MTA `TS` in MIQP mode), sorted descending, plus run metadata (method,
alpha, reference provenance, solver). Exportable as a deterministic table like every other
result service. The `mta` (single-QP) and `mta_miqp` modes are single-solve and have no
robustness gate, so they score by `TS` directly; non-reverting knockouts there approach but
do not always equal exactly zero (a forced reaction may move by the indicator tolerance).
The clear separation between the top target and the rest is preserved in every mode.

## Module Layout

```
cmm/
  core/
    flux_state.py        # FluxState reference primitive (Phase 1/3)
    results.py           # TargetRanking (shared with FSEOF/OptKnock)
    solvers.py           # capability detection + SolverCapabilityError
  omics/
    differential.py      # two-state diff expression -> DirectionMap (Phase 5)
  features/
    revert.py            # revert_targets(...) orchestration + scoring
    _perturbation.py     # shared KO enumeration/runner (also used by batch MOMA/ROOM)
```

Public surface (solver-neutral, GUI-free, mirrors the existing `fba`/`fva` style):

```python
def revert_targets(
    model: Model,
    source_condition: Condition,
    reference_state: FluxState,
    direction: DirectionMap,
    *,
    targets: Iterable[str] | None = None,
    method: Literal["rmta", "mta"] = "rmta",
    alpha: float = 0.9,
    perturbation: Literal["gene", "reaction"] = "gene",
) -> TargetRanking: ...
```

## Performance

Genome-scale x all-gene KO is the expensive case (one-to-three QP solves per gene).

- Reuse the Phase 3 batch perturbation runner for enumeration and parallel solving.
- Solve inside a single persistent `model` context, mutating only the KO bounds per
  iteration and restoring them, to avoid rebuilding the problem.
- Warm-start each KO from `v_ref`.
- Support a candidate-subset `targets` list so the GUI can score a shortlist interactively
  before committing to a full screen.

## Testing

The existing single-metabolite `toy_model` cannot exercise direction logic. Add a small
branched toy model (one shared precursor feeding two competing branches) so a knockout on
one branch demonstrably shifts flux toward the other. Cover:

- `DirectionMap` construction from a two-state expression table with up/down/steady labels.
- A KO that improves agreement scores above one that does not (ranking correctness).
- `rTS == 0` when a KO helps best-case but hurts worst-case (robustness gate).
- `SolverCapabilityError` raised for `method="mta"` when no MIQP solver is present.
- Determinism of the exported `TargetRanking` table.

## Open Questions

- Default reference provenance: pFBA is reproducible but rMTA papers use sampling means.
  Plan: allow both, default to pFBA for small models and recommend sampling for
  genome-scale, recorded in `FluxState.provenance`.
- Differential-expression thresholding (fold-change/percentile) belongs in the omics layer
  so single-state and two-state methods share normalization; revert only consumes the
  resulting labels.
