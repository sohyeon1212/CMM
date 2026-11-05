# ADR-0005 — `v_ref` and ε come from samples of one iMAT activity pattern

**Status:** Accepted · 2026-11-05 · supersedes [ADR-0002](0002-reference-state-and-epsilon.md)
(both its `v_ref` and its ε decision). An earlier draft of this record, never released, sampled
the E-Flux2-constrained space for ε while keeping the E-Flux2 `v_ref`.

## Context

ADR-0002 replaced the papers' reference state with E-Flux2 or LAD because CMM had no iMAT, and
because neither MTA paper says how to turn continuous expression into iMAT's high/moderate/low
classes. Valcárcel et al. (2019, §2.3) supply that rule: in each sample a gene is high above
mean + δ·SD and low below mean − δ·SD of the sample's log2 expression, and a gene keeps a class
only when β of the samples agree (β = 100%; δ = 0.3 for Alzheimer's, 0.7 for prostate cancer).
The MIT-licensed `imatpy` package implements the iMAT MILP of Shlomi et al. (2008).

With both pieces, CMM can follow Yizhak et al. (2013, Methods): solve iMAT, sample 2,000 flux
distributions from the space consistent with the activity pattern of one optimal solution, take
`v_ref` as their mean, and derive each changed reaction's ε from the same samples.

## Decision

SC-02 and the GUI's Revert tab build `v_ref` this way and only this way; the E-Flux2/LAD
reference is removed from both.

- `cmm.omics.discretize_expression` applies the rMTA rule over the model's measured genes. A
  linear file is first transformed with `log2(x + 1)`. `imat_delta` has no default.
- `cmm.omics.imat` builds the MILP with imatpy's public functions, with ε = 1 and imatpy's
  `threshold = 0`, which makes its low-expression constraint the paper's `v = 0`.
- `cmm.features.imat_reference_state` fixes each matched reaction to its solved pattern (active
  reactions keep their direction with |v| ≥ ε; inactive ones are fixed at 0), draws 2,000
  seeded single-chain COBRApy OptGP samples with thinning 500, and averages them.
- `epsilon="per_reaction"` applies `εᵢ = max(t₀.₉₅,ₙ₋₁ · sᵢ / √n, 0.001)` to those samples.

## Consequences

- Three differences from the published runs remain and are recorded in every run as
  `reference_state_deviations`:
  - The GPR maps AND to min and OR to max. The iMAT paper's Methods text states the reverse,
    which would make a complex as available as its best-expressed subunit.
  - COBRApy builds its own warm-up points; optGpSampler was run with 5,000.
  - Sampling uses one chain, for reproducibility.
- Using one of possibly several iMAT optima is the papers' own procedure ("one of iMAT's optimal
  solutions"), not a difference. The pattern the solver returned is written to
  `02_reference/imat_reference_pattern.csv` and the gene classes to
  `02_reference/gene_expression_states.csv`.
- `v_ref` depends on δ, β, the sample count, thinning and seed. All of these are configuration
  and provenance. Runs made with an E-Flux2 reference are not comparable, and a config that
  still names `reference_method` is rejected with the reason.
- iMAT adds one MILP and the sampling cost.
- The E-Flux2 and LAD services themselves remain in `cmm.omics` for other uses.

## Source

[`src/cmm/omics/imat.py`](../../src/cmm/omics/imat.py),
[`src/cmm/features/revert_reference.py`](../../src/cmm/features/revert_reference.py),
[`src/cmm/features/revert_epsilon.py`](../../src/cmm/features/revert_epsilon.py);
Shlomi et al. (2008) *Nat Biotechnol* 26:1003; Yizhak et al. (2013) *Nat Commun* 4:2632,
Methods; Valcárcel et al. (2019) *Bioinformatics* 35:4350, §2.3 and Supplementary Figs. S8–S9.
