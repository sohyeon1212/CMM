# Changelog

## 0.3.0 — 2026-07-10

### Scientific correctness

- Connected `robustknock` to StrainDesign's actual three-level `ROBUSTKNOCK` module and
  evaluated maximum and guaranteed product at exactly optimal growth.
- Reimplemented `mta` and `rmta` with the published MIQP/MOMA/MIQP workflow, L1
  transformation score, and rMTA Equation 9. Renamed the historical QP approximation to
  `rmta_continuous`.
- Made E-Flux2 a strict two-stage QP with a full-objective default and an explicitly labeled,
  opt-in L1 fallback.
- Aligned FSEOF/FVSEOF scan origins, enforced product/biomass constraints, trend handling,
  and actionable-target filtering with their documented method contracts.

### Reliability and validation

- Added finite/range/completeness checks across conditions, media, expression, flux states,
  perturbations, and scan parameters.
- Added deterministic model fingerprints and solver/runtime/parameter provenance to
  numerical result objects.
- Added direct COBRApy comparisons, the official MTA test topology, non-optional iJO1366
  checks, scientific sensitivity tests, and malformed-input/GUI-state regressions.
- Added an 80% branch-coverage gate, Ruff formatting/linting, mypy checks for scientific
  services, citation validation, and locked-dependency vulnerability auditing.

### Reproducibility and distribution

- Bumped the package to 0.3.0, constrained supported dependency/API ranges, and committed
  `uv.lock` for Python 3.10–3.12.
- Added `CITATION.cff`, an MIT `LICENSE`, scientific validation notes, and
  corrected solver/method documentation.
- Release builds now verify tag/version agreement, rerun all quality gates, validate wheel
  metadata, and install the built wheel in a clean environment.

### Behavior changes

- Code that relied on E-Flux2 silently falling back to pFBA must now request
  `allow_l1_fallback=True` and handle the `eflux2_l1_fallback` method label.
- `rmta` now requires MIQP and denotes the published robust workflow. Use
  `rmta_continuous` only when the historical QP heuristic is intentionally desired.
- Non-optimal FBA results return no objective value or flux vector instead of exposing
  solver-generated invalid numbers.
