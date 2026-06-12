# CMM Autonomous Build — Completion Report

Self-paced PDCA loop over the updated CMM plan, with the revert-metabolism (normalization
target) feature as the headline. Generator and verifier were always separate agents
(Strategy 1); failures were turned into durable rules (Strategy 2); every story had an
explicit done/fail test (Strategy 3); the lead model (Fable 5) designed and wrote the
interdependent core while cheaper Sonnet agents ran independent verification (Strategy 4).

## Final gate

- `pytest`: **36 passed** (was 7 at session start) — runtime `.venv/bin/python`, offscreen Qt.
- `ruff check src tests`: **clean**.
- Layering: **no Qt imports** in `cmm.core` / `cmm.features` / `cmm.omics` (only `cmm.app`).
- New tested code: ~1730 LOC across core, features, omics, app.
- GUI: 6 offscreen screenshots in `temp_figures/` driving real scenarios.

## Stories

| Story | Scope | Result | Verifier |
|-------|-------|--------|----------|
| G1 | Solver tiers, `FluxState`, `TargetRanking` | DONE | Sonnet: PASS, no bugs |
| G2 | Perturbation runner, MOMA (L1/L2), ROOM | DONE | Sonnet: PASS, no bugs |
| G3 | Two-state differential expression -> `DirectionMap` | DONE | Sonnet: PASS, no bugs |
| G4 | **Revert metabolism (rMTA + MTA + MIQP)** | DONE | Sonnet: PASS (g2 #1; +1 doc-sync fix) |
| G5 | Qt shell + offscreen screenshot scenarios | DONE | Sonnet: PASS (6 non-blank PNGs, clean layering, g2 #1) |
| G6 | Final gate + report | DONE | this report |

## What was built

- `cmm.core.solvers` — LP/QP/MILP/MIQP capability detection + `SolverCapabilityError`.
- `cmm.core.flux_state.FluxState` — reference flux primitive (pFBA / sampling-mean / imported),
  L1/L2 distance, serialize.
- `cmm.core.results.TargetRanking` — shared, deterministic ranked-target result + export.
- `cmm.features._perturbation` — gene/reaction KO enumeration + reverting application.
- `cmm.features.comparison` — MOMA (L1 LP, L2 QP) and ROOM (MILP) vs a `FluxState`.
- `cmm.omics.differential` — two-state expression -> per-reaction direction via GPR AST
  (AND=min, OR=max) combined with reference flux sign.
- `cmm.features.revert.revert_targets` — **the headline**: ranks knockouts that revert a
  source state toward a target state. Methods: `rmta` (robust, QP×3), `mta` (single QP),
  `mta_miqp` (original MTA, MIQP-gated). On the disease toy, the disease-branch knockout
  (g2 / R2) ranks #1 across all three methods.
- `cmm.app` — Qt workbench (model browser, flux-range slider, FBA/FVA, revert tab) and an
  offscreen screenshot harness.

## Screenshots (`temp_figures/`)

`01_model_loaded`, `02_fba` (objective 10, optimal), `03_fva`, `04_flux_range_slider`,
`05_revert_rmta` (g2 #1, score 0.8541), `06_revert_mta_miqp` (g2 #1, score 30).

## Failure rules learned (continual-learning memory, `rules.md`)

R1 use .venv not anaconda (Qt clash); R2 QP/MILP/MIQP only in .venv (Gurobi/CPLEX);
R3 offscreen `grab()` for headless screenshots; R4 cobra 0.31 `knock_out_model_genes`;
R5 ROOM pins the original objective to `Solution.objective_value` (never pass 0.0);
R6 add optlang vars then flush before constraints, build one-off MIQP on `model.copy()`.

## Deferred (explicitly out of this slice, recorded not dropped)

Dynamic FBA, random-sampling sampler service, FSEOF/FVSEOF, OptKnock/RobustKnock, the
enzyme-constrained (GECKO) layer, and the full multi-section desktop product. These are
independent of the revert-metabolism vertical slice; the cross-cutting primitives
(`FluxState`, solver tiers, `TargetRanking`, perturbation runner, differential expression)
were built so they drop in without rework. The design and roadmap for each remain in
`docs/feature-roadmap.md` and `docs/design-revert-metabolism.md`.
