# Loop Success Criteria (Strategy 3: hard termination + test definitions)

> Archived pre-0.3.0 criteria. Use `docs/VALIDATION.md` for current scientific acceptance
> criteria and reproducibility commands.

The autonomous build runs as an ordered set of stories. Each story has an explicit
**done test** and **fail condition**. A story is PASS only when an independent verifier
agent (Strategy 1: verifier != generator) confirms the done test against real command
output. The loop terminates when every committed story is PASS and the final gate is clean.

Runtime: all tests/GUI run via `.venv/bin/python` (uv env: Gurobi+CPLEX+glpk, clean PyQt5).
NOT anaconda python (its Qt collides with PyQt5 and it lacks QP/MILP solvers).

Global invariants (apply to every story):
- All new code lives under `CMM/`. No edits outside `CMM/`.
- `.venv/bin/python -m pytest CMM/tests -q` passes (no regressions).
- `.venv/bin/ruff check CMM/src` clean (or only pre-existing, documented warnings).
- Services are GUI-free and import without Qt; only `cmm.app` imports Qt.

## Committed stories

### G1 — Core primitives & solver tiers
- DONE: `cmm.core.solvers` detects LP/QP/MILP/MIQP and raises `SolverCapabilityError`;
  `cmm.core.flux_state.FluxState` (provenance + distance + serialize) and
  `cmm.core.results.TargetRanking` exist with tests.
- FAIL: any new test red, or capability detection wrong for the active solver.

### G2 — Perturbation runner + MOMA/ROOM
- DONE: `cmm.features._perturbation` enumerates+applies gene/reaction KOs;
  `cmm.features.comparison` runs L1-MOMA (LP), L2-MOMA (QP), and ROOM vs a `FluxState`.
- FAIL: MOMA/ROOM disagree with a hand-checked branched-toy expectation.

### G3 — Omics differential expression
- DONE: `cmm.omics.differential` maps a two-state expression table to a per-reaction
  `DirectionMap` (+1/-1/0) through GPR + thresholds, with tests.
- FAIL: direction labels wrong on the toy two-state table.

### G4 — Revert metabolism (rMTA) [HEADLINE]
- DONE: `cmm.features.revert.revert_targets(...)` implements rMTA (QP best/moma/worst ->
  rTS) and optional MIQP MTA; ranks the known normalization target #1 on the toy disease
  model; robustness gate (`rTS==0` when worst-case hurts) and solver gating tested.
- FAIL: known target not ranked #1, or non-deterministic ranking.

### G5 — Minimal desktop GUI + screenshot scenarios
- DONE: `cmm.app` Qt shell over the services; `cmm.app.screenshots` drives >=4 scenarios
  offscreen and writes non-blank PNGs to `CMM/temp_figures/`; a verifier confirms the
  PNGs are non-trivial (size + pixel variance).
- FAIL: app imports Qt into core, or any scenario screenshot is blank/missing.

### G6 — Final gate
- DONE: full pytest + ruff clean; completion report written; ledger closed.
- FAIL: any committed story not PASS.

## Round 2 — Production design, publication figures, usability

Goal: succinate-production analysis on e_coli_core, publication-quality figures vs reference plots,
theoretical yield, easy bound editing, and an adversarial review. Runtime/invariants as above.
Model: cobra `load_model("textbook")` = e_coli_core (the bundled `*.json` is an Escher map,
not a model). Figures use matplotlib Agg at >=200 DPI.

### R2-G1 — Production-design services (`cmm.features.production`)
- DONE: `theoretical_yield`, `production_envelope`, and `fseof` run on e_coli_core for
  succinate (EX_succ_e from EX_glc__D_e), with tests pinning known values (anaerobic
  succinate yield, growth-coupled envelope shape, FSEOF amplification targets include
  TCA/anaplerotic reactions).
- FAIL: yields or envelope endpoints disagree with hand-checked cobra values.

### R2-G2 — Publication-quality visualization (`cmm.visualization`)
- DONE: matplotlib (Agg) figures — production envelope, FSEOF profile, flux-comparison bar,
  theoretical-yield summary — saved at >=200 DPI with consistent styling; a verifier confirms
  the PNGs are non-blank and compares style against reference publication plots.
- FAIL: blank/garbled figures, or pixelated/low-DPI output unfit for a paper.

### R2-G3 — Easy model bound editing
- DONE: GUI reaction table bounds are editable and apply to the model (round-trip verified);
  a `set_bounds`/`Condition` programmatic path stays in sync; re-running FBA reflects edits.
- FAIL: edits do not persist or do not change FBA results.

### R2-G4 — E. coli succinate scenario + screenshots
- DONE: a harness loads e_coli_core, runs theoretical yield + envelope + FSEOF for succinate,
  generates the publication figures, and drives the GUI (incl. a bound edit that increases
  succinate) capturing >=5 offscreen screenshots to `temp_figures/`.
- FAIL: any scenario errors or produces a blank capture.

### R2-G5 — Adversarial review + usability (Workflow, multi-lens)
- DONE: parallel adversarial reviewers (scientific correctness, figure quality vs reference plots,
  usability/ergonomics) run; their confirmed findings are triaged, fixed or recorded.
- FAIL: an unaddressed correctness finding remains.

## Round 4 — Solver status, condition omics, transformation targets, OptKnock

Driven by user follow-up. Same invariants. Data: `test_data/JS66_exp.csv` (5618 genes ×
conditions JS66_2/ALECO70_2/ALECO2_2) + `test_data/JS66_v15.xml` (1112/1122 genes overlap).
`straindesign` is installed for OptKnock/RobustKnock.

### R4-G1 — Solver status check (is it gurobi/QP-capable?)
- DONE: `cmm.core.solvers.solver_status()` reports active solver + capabilities + a
  `recommended` flag (gurobi/cplex); GUI Config menu / header shows it and warns when the
  solver is not QP-capable. (Also fix the E-Flux2 objective-floor direction latent bug.)
- FAIL: status wrong for the active solver, or no warning when QP is unavailable.

### R4-G2 — Flux map metabolite node size
- DONE: `escher_flux_map` metabolite markers are sized appropriately (not oversized) relative
  to the canvas; regenerated `ecoli_05_flux_map.png` / `fig_escher_*` look balanced.
- FAIL: nodes still dominate the edges.

### R4-G3 — Multi-condition omics flux prediction + log-change
- DONE: load a multi-condition expression table, predict per-condition flux with E-Flux2 AND
  LAD on JS66, and compute/visualize log2 fold-change of flux between conditions; both methods
  tested.
- FAIL: a method errors on JS66, or log-change is computed wrong (e.g. sign/zero handling).

### R4-G4 — Condition A→B transformation target finder (MTA + MOMA)
- DONE: `transformation_targets(model, source_state, target_state, method='moma'|'mta')` ranks
  candidate gene knockouts by how close they bring the flux to the target state; works with
  omics-predicted FluxStates; both methods tested on the branched model and reduce distance.
- FAIL: ranking does not move flux toward the target, or methods disagree with hand-check.

### R4-G5 — OptKnock / RobustKnock
- DONE: `cmm.features.strain_design.optknock(...)` / `robustknock(...)` wrap straindesign,
  return ranked knockout designs that couple a product to growth on e_coli_core (e.g.
  succinate); MILP-gated; tested (a known coupling design is found).
- FAIL: no design found for a known-couplable product, or crashes.

## Round 5 — Media, pFBA, MOMA/ROOM templates, branding, README, validation

User follow-up. Same invariants. Screenshots/tutorial use the PUBLIC e_coli_core model only;
NO mention of the private JS66 strain in any published doc.

### R5-G1 — Easy media composition
- DONE: `cmm.core.media` with a `Medium` type + preset media (glucose minimal aerobic/anaerobic,
  etc.), pattern-tolerant exchange resolution, `apply_to(model)` via cobra's medium API; GUI
  medium selector applies it and re-runs FBA.
- FAIL: applying a medium does not set the right exchange bounds / changes the wrong reactions.

### R5-G2 — pFBA
- DONE: `cmm.core.simulation.pfba` returns a `FluxSolution` (parsimonious FBA); GUI "Run pFBA"
  button; result matches `cobra.flux_analysis.pfba`.
- FAIL: pFBA total flux not minimal / disagrees with cobra.

### R5-G3 — MOMA/ROOM reference-template selection
- DONE: a reference-flux builder lets MOMA/ROOM use a template from FBA, pFBA, LAD, or E-Flux2;
  GUI Comparison tab runs MOMA/ROOM against a chosen template + a knockout; both methods work.
- FAIL: the chosen template is not actually used as the MOMA/ROOM reference.

### R5-G4 — Branding: CMM = Cellular Metabolic Modeling Workbench
- DONE: the acronym is expanded consistently in the GUI title/About, README, and package docs.
- FAIL: inconsistent or missing expansion.

### R5-G5 — Publication-quality README + scrub JS66
- DONE: README has overview, install, a worked tutorial, and embedded screenshots (e_coli_core),
  at paper quality; ALL JS66 references removed from README/temp_figures/scenario/docstrings
  (genome-scale scenario genericized to a model-path argument).
- FAIL: README is thin, screenshots missing, or any JS66 mention remains in published material.

### R5-G6 — GUI-vs-cobrapy correctness validation
- DONE: tests assert the GUI/service results (FBA, pFBA, FVA, MOMA, theoretical yield, medium)
  equal direct cobrapy computations on e_coli_core to numerical tolerance.
- FAIL: any GUI result diverges from the cobrapy reference.

## Deferred (documented, not blocking termination)
Heavy/independent roadmap items not required for the headline vertical slice:
dynamic FBA, random sampling sampler service, FSEOF/FVSEOF, OptKnock/RobustKnock,
enzyme-constrained (GECKO) layer, full multi-section desktop product. Each is recorded
in the ledger as DEFERRED with rationale so scope is explicit, not silently dropped.
