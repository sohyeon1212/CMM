# Changelog

## Unreleased

### Added

- Simulation: FBA and pFBA fluxes are shown in separate columns (Reaction / FBA flux /
  pFBA flux / FVA range) so running pFBA no longer overwrites the FBA result; pFBA's minimal
  total flux is shown directly under the objective value.
- Production: a result table beside each plot — FSEOF/FVSEOF amplify/knockdown targets with
  their low/high enforced-level fluxes, and the production-envelope growth range per product
  flux — with a "Show all reactions" toggle that also lists unchanged reactions.
- Comparison (single run): a "Significant change ≥ X % of reference" threshold (default 3%)
  replacing the fixed 1e-6 cutoff, so alternate-optimum drift (notably ROOM) no longer reads
  as a knockout response; the solve is cached and re-filters without re-solving.
- Comparison (batch): the table reports wild-type and post-knockout biomass, an essentiality
  flag, and — when a target product is selected — that product's wild-type and post-knockout
  flux columns.
- Omics: multi-condition expression tables are supported in one tab, computing one predicted-
  flux column per selected condition, with a "Show all reactions" toggle.
- `cmm-guide` project-local skill: an agent-facing operating protocol describing CMM's
  analyses, a goal→method decision guide, the solver requirement matrix, and pitfalls.

### Changed

- Comparison: a two-panel knockout picker (searchable catalogue on the left, chosen knockout
  set on the right) replaces Ctrl/Shift-click selection, making the selection visible and
  clearable.
- Comparison: LAD/E-Flux2 reference templates are disabled (greyed out, not selectable) until
  their integration method has actually been computed on the Omics tab.
- Omics: loading an expression file only stores it and shows its filename; a separate Compute
  button runs the selected method, so the method can change without reloading.
- Revert / Transform: the loaded source/target expression filename is shown next to each input.
- Production: removed the redundant Run FBA button (duplicated the Simulation tab and produced
  no Production-tab output).

### Fixed

- Comparison batch (MOMA/ROOM) no longer aborts the whole run when a lethal knockout makes the
  model infeasible; such a knockout is recorded as infeasible and the run continues.
- Revert / MTA (`_mta_miqp`, `_mta_qp`) no longer crash under Gurobi when a lethal knockout is
  infeasible (backend "Unable to retrieve attribute 'X'"); feasibility is probed before reading
  primals, so infeasible knockouts are skipped and the ranking completes.

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
