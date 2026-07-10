# Revert metabolism: MTA and rMTA

## Scope

CMM ranks gene or reaction knockouts intended to move a source metabolic state toward a
target state. `mta` implements the MTA mixed-integer quadratic formulation of Yizhak et al.
(2013). `rmta` implements the robust best-case/MOMA/worst-case workflow and Equation 9 of
Valcárcel et al. (2019). `mta_miqp` is retained as an API alias for `mta`.

The former continuous approximation is not labeled as rMTA. It is available only as
`rmta_continuous` and its result metadata reports `formulation="continuous_heuristic"`.

## Inputs

- `model`: a COBRApy model with the source-condition bounds.
- `source_condition`: optional bound/objective overrides applied during target scoring.
- `reference_state`: a complete, finite `FluxState` for the source condition. Every model
  reaction must be present.
- `direction`: a `DirectionMap` with values `+1` (desired increase), `-1` (desired decrease),
  or `0` (steady).
- `targets`: optional candidate genes or reactions. Omitting it screens all candidates.
- `alpha`: MTA weighting, in `[0, 1]`; the rMTA paper default used by CMM is `0.66`.
- `epsilon`: minimum directional change represented by each MTA binary switch.
- `parameter_k`: Equation 9 scale factor, default `100`.

For transcript-style gene identifiers such as `g2.1` and `g2.2`,
`transcript_separator="."` groups transcripts into a single gene perturbation before GPR
evaluation.

## Published MTA solve

For steady reactions, MTA penalizes squared deviation from the source flux. Directional
reactions receive a binary success switch constrained by `epsilon` and reaction-specific
big-M bounds. After dropping constants, CMM minimizes

```text
(1 - alpha) · Σsteady (v_i - vref_i)² - (alpha / 2) · Σdirectional z_i
```

subject to steady-state mass balance, reaction bounds, the active knockout, and the binary
direction constraints. This is an MIQP. The solve is built on a model copy so binary variables
cannot leak between candidates.

## Transformation score and robust score

For an optimal candidate flux `v`, the published L1 transformation score is

```text
TS = (successful directional movement - unsuccessful directional movement)
     / steady-reaction L1 disturbance
```

If the denominator is zero, CMM returns `0` for `0/0`, `+∞` for a beneficial nonzero
numerator, and `-∞` for an adverse one. These cases remain sortable and are not converted to
arbitrary finite constants.

Robust rMTA performs three candidate solves:

1. `bTS`: published MTA in the requested direction.
2. `mTS`: full L2 MOMA relative to the source state.
3. `wTS`: published MTA after swapping forward and backward direction sets, scored against
   that swapped direction.

The final score follows Equation 9:

```text
if bTS > 0 and mTS > 0 and wTS < 0:
    rTS = mTS · K · (bTS - wTS)
else:
    rTS = mTS
```

Non-optimal candidates receive `-∞`; the ranking metadata records the count. Every row exposes
`bTS`, `mTS`, and `wTS` for audit.

## Source-state preprocessing

The original MTA/rMTA studies construct the source model with expression contextualization
and flux sampling. CMM's desktop workflow instead creates a deterministic source state by
running E-Flux2 on the source expression at full objective (`objective_fraction=1.0`) and
derives source→target directions from differential expression plus GPR logic. This is a
documented CMM preprocessing variant; the optimization and score stages described above are
the published formulations.

For a manuscript that claims the complete original pipeline, supply an externally generated
source `FluxState` from the study's contextualization/sampling protocol through the Python
API. Record that provenance in the state metadata.

## Validation

`tests/test_revert.py` reconstructs the official COBRA Toolbox MTA test network independently
and checks the published Equation 9 and expected target signs (`g4 > 0`, `g2 < 0` for MTA;
`g4 > 0` for rMTA). `tests/test_scientific_sensitivity.py` verifies that the positive control
and relative ranking persist across `alpha = 0.3, 0.4, 0.66, 0.8`.

Primary references:

- Yizhak et al. (2013), *Nature Communications* 4:2632,
  <https://doi.org/10.1038/ncomms3632>.
- Valcárcel et al. (2019), *Bioinformatics* 35:4350–4355,
  <https://doi.org/10.1093/bioinformatics/btz231>.
