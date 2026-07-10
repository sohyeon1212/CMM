# Build Ledger (append-only)

> Archived implementation ledger. Entries describe the state at the time they were written,
> not the validated 0.3.0 method contracts.

Durable record of the autonomous build loop. Event types: START, GENERATE, VERIFY(PASS|FAIL),
FIX, RULE, CHECKPOINT, DEFERRED, CLOSE. Verifier is always a separate agent from the generator.

| # | event | story | model | detail |
|---|-------|-------|-------|--------|
| 1 | START | loop  | Fable 5 (lead) | Recon done; runtime=.venv; solvers=Gurobi/CPLEX/glpk; GUI=offscreen grab OK. Seeded R1–R3. |
| 2 | START | G1    | Fable 5 | Begin core primitives (solvers, flux_state, results). |
| 3 | GENERATE | G1 | Fable 5 (lead) | solvers/flux_state/results + tests; 16 pass, ruff clean. |
| 4 | VERIFY(PASS) | G1 | Sonnet (Explore) | Independent verifier: no bugs; confirmed gurobi has LP/QP/MILP/MIQP. |
| 5 | CHECKPOINT | G1 | Fable 5 | Core primitives accepted. |
| 6 | START | G2 | Fable 5 | Perturbation runner + MOMA/ROOM. |
| 7 | GENERATE | G2 | Fable 5 (lead) | _perturbation + comparison + tests. |
| 8 | VERIFY(FAIL) | G2 | pytest | ImportError: find_gene_knockout_reactions removed in cobra 0.31. |
| 9 | RULE | G2 | Fable 5 | R4: use knock_out_model_genes in a `with model:` context. |
| 10 | FIX | G2 | Fable 5 | Adopt knock_out_model_genes adapter; 5/7 pass. |
| 11 | VERIFY(FAIL) | G2 | pytest | ROOM returned all-zero fluxes (distance=4 shutdown). |
| 12 | RULE | G2 | Fable 5 | R5: set Solution.objective_value to model-objective-at-reference (ROOM pins old objective). |
| 13 | FIX | G2 | Fable 5 | Compute reference objective value; 7/7 pass, full suite 23 pass, ruff clean. |
| 14 | CHECKPOINT | G2 | Sonnet (Explore) | VERIFY(PASS): bounds restore, MOMA/ROOM reroute, solver gating all confirmed; no bugs. |
| 15 | START | G3 | Fable 5 | Omics two-state differential expression -> DirectionMap. |
| 16 | GENERATE | G3 | Fable 5 (lead) | differential.py (GPR AST AND=min/OR=max) + tests; 4 pass first try. |
| 17 | VERIFY(PASS) | G3 | Sonnet (Explore) | Independent: reverse-flux sign branch + GPR combos confirmed; no bugs. |
| 18 | CHECKPOINT | G3 | Fable 5 | DirectionMap accepted. |
| 19 | START | G4 | Fable 5 | Headline: rMTA revert_targets. |
| 20 | GENERATE | G4 | Fable 5 (lead) | revert.py (rMTA QP + robust score); probe: g2 ranks #1 (rTS 0.85). |
| 21 | VERIFY(FAIL) | G4 | pytest | MIQP path: optlang ContainerAlreadyContains on mta_z_R2. |
| 22 | RULE | G4 | Fable 5 | R6: add vars then flush before constraints; build MIQP on model.copy(). |
| 23 | FIX | G4 | Fable 5 | Two-phase add on copy; 7/7 G4 pass, full suite 34 pass, ruff clean. |
| 24 | VERIFY(PASS) | G4 | Sonnet (Explore) | g2 #1 (rTS 0.854, 13x next); robustness gate, MIQP no-leak, solver gating, non-mutation all confirmed; no code bugs. |
| 24b | FIX | G4 | Fable 5 | Doc-sync: design doc rTS formula was stale (textbook gate would zero every gene); updated to the implemented additive worst-case form + MIQP note. |
| 24c | CHECKPOINT | G4 | Fable 5 | Headline revert-metabolism accepted. |
| 25 | START | G5 | Fable 5 | Minimal Qt shell + offscreen screenshot scenarios. |
| 26 | GENERATE | G5 | Fable 5 (lead) | cmm.app (main_window, screenshots) + smoke test; 6 PNGs rendered. |
| 27 | FIX | G5 | Fable 5 | Header label background bled from global QWidget rule; made QLabel transparent. |
| 28 | VERIFY(PASS) | G5 | Sonnet (Explore) | 6 non-blank PNGs (1180x740, ~600-870 colors each), smoke tests pass, layering clean, g2 #1 in both revert scenarios. |
| 29 | CHECKPOINT | G5 | Fable 5 | GUI + screenshots accepted. |
| 30 | CLOSE | loop-R1 | Fable 5 (lead) | Round 1 (revert metabolism): all 6 stories PASS; 36 tests. Loop terminated on success. |
| 31 | START | R2 | Fable 5 | Round 2: production design, publication figures, bound editing, adversarial review (e_coli succinate). |
| 32 | GENERATE | R2-G1 | Fable 5 (lead) | production.py (theoretical_yield, production_envelope, fseof) + tests. |
| 33 | VERIFY(FAIL) | R2-G1 | pytest | aerobic yield 2.96 (actual-flux denominator) + anaerobic infeasible (substrate=CO2). |
| 34 | RULE | R2-G1 | Fable 5 | R7 (fix substrate uptake for yield) + R8 (detect limiting carbon source, not most-negative lb). |
| 35 | FIX | R2-G1 | Fable 5 | Fixed substrate + detection; 6/6 pass (yields 1.39/1.64). |
| 36 | VERIFY(PASS) | R2-G1 | Sonnet (Explore) | Yields hand-checked, FSEOF FRD7/PPC present, knockdowns biologically sane, no mutation, edges graceful. |
| 37 | GENERATE | R2-G2 | Fable 5 (lead) | visualization (envelope/FSEOF/flux/yield figures, 300 DPI, Okabe-Ito) + tests; 5/5 pass. |
| 38 | GENERATE | R2-G3 | Fable 5 (lead) | GUI editable bounds (in-table + programmatic + invalid revert) + tests; pass. |
| 39 | GENERATE | R2-G4 | Fable 5 (lead) | Production tab w/ embedded figures + e_coli succinate scenario harness. |
| 40 | VERIFY(FAIL) | R2-G4 | self-review | Scenario claimed succinate increase but FBA gave 0->0 (growth mode doesn't excrete succinate). |
| 41 | RULE | R2-G4 | Fable 5 | R9: block competing fermentation secretions (upper=0) to force succinate; 0->9.885 via FRD7. |
| 42 | FIX | R2-G4 | Fable 5 | Growth-coupled design + fixed embedded-figure DPI clipping + slider itemChanged guard. |
| 43 | START | R2-G5 | Fable 5 | Adversarial multi-lens review (Workflow): science, figures vs reference plots, usability, code. |
| 44 | VERIFY(concerns) | R2-G5 | Workflow (opus+3×sonnet) | 35 findings: 4 critical, 8 high, 13 medium, 10 low. All 4 lenses "concerns". |
| 45 | RULE | R2-G5 | Fable 5 | R10 (FSEOF magnitude classification), R11 (atomic bounds), R12 (CO2 yield disclosure + currency-metabolite exclusion). |
| 46 | FIX | R2-G5 | Fable 5 | Fixed all 35: FSEOF signed→magnitude (critical); network_flux_map added (critical); infeasible-FBA + crossing-bounds crashes guarded (critical); CO2 disclosure, carbon-aware substrate, distinguishable FSEOF lines, frozen-df copy, context-guarded fallback (high); +13 medium +10 low. |
| 47 | VERIFY(PASS) | R2-G5 | pytest | 58 tests pass (was 48: +10 production/crash regressions), ruff clean, layering clean. Crash repros locked as tests. |
| 48 | CHECKPOINT | R2-G5 | Fable 5 | Adversarial review findings all resolved + re-verified. |
| 49 | CHECKPOINT | loop-R2 | Fable 5 (lead) | Round 2 adversarial fixes complete. |
| 50 | START | R2-G6 | Fable 5 | User follow-up: real curated Escher map + genome-scale (JS66) test + honest figure comparison. |
| 51 | RULE | R2-G6 | Fable 5 | R13 (render Escher map JSON layout, not force-layout) + R14 (genome-scale medium / substrate). |
| 52 | FIX | R2-G6 | Fable 5 | Added escher_flux_map (renders curated Escher layout coloured by flux) + GUI Flux Map tab + substrate selector; tested JS66 (1480 rxn, all <0.2s, FSEOF→FRD2/SUCCt). |
| 53 | VERIFY(PASS) | R2-G6 | pytest | 61 tests (was 58: +3 map renderers), ruff clean, layering clean. e_coli Escher map matches the curated layout + flux coloring. |
| 54 | CHECKPOINT | loop-R2 | Fable 5 (lead) | Round 2 complete incl. curated-Escher-parity flux map + genome-scale validation. |
| 55 | START | R3 | Fable 5 | User: JS66 GUI scenario + screenshots; desktop-style UI menus; omics integration (E-Flux2/LAD) tested. |
| 56 | GENERATE | R3 | Fable 5 (lead) | cmm.omics.expression (gene_to_reaction_weights, eflux2 [QP+pFBA fallback], lad) mirroring published reference algorithms + 9 tests. |
| 57 | GENERATE | R3 | Fable 5 (lead) | GUI: desktop-style menu bar (Analysis/Model/Config) + Omics tab (E-Flux2/LAD, CSV load, demo) + substrate selector. |
| 58 | GENERATE | R3 | Fable 5 (lead) | genome_scale_scenario.py: JS66 (1480 rxn) GUI run — FBA, glucose bound-edit, production design, omics; 8 screenshots. |
| 59 | VERIFY(pending) | R3 | Sonnet (Explore) | Independent omics verifier launched (E-Flux2/LAD science vs published references + Kim 2016). |
| 60 | VERIFY(PASS) | R3 | pytest | 72 tests (was 61: +9 omics, +2 GUI), ruff clean, layering clean. JS66 GUI: yield 1.5 after glucose bound-edit, E-Flux2 maps 1122 genes. |
| 61 | VERIFY(PASS) | R3 omics | Sonnet (Explore) | E-Flux2/LAD match the published reference to 4 decimals; GPR/QP-fallback/non-mutation confirmed. Found 1 low latent bug. |
| 62 | RULE/FIX | R3 omics | Fable 5 | R15: E-Flux2 objective floor must branch on objective_direction (lb for max, ub for min). Fixed. |
| 63 | START | R4 | Fable 5 | User: solver status, flux-map node size, multi-condition omics+logchange, A→B transformation targets, OptKnock/RobustKnock. |
| 64 | GENERATE | R4-G1 | Fable 5 | solver_status() + GUI Config status/warning; eflux2 floor-direction fix. |
| 65 | GENERATE | R4-G2 | Fable 5 | escher_flux_map metabolite node ms 4.5→2.2 (edges now stand out). |
| 66 | GENERATE | R4-G3 | Fable 5 | cmm.omics.conditions (multi-condition predict + flux_log_change + sign_flips) + flux_log_change_figure; verified on JS66_exp.csv (E-Flux2+LAD, 3 conditions optimal). |
| 67 | GENERATE | R4-G4 | Fable 5 | cmm.features.transformation (A→B targets via MOMA and MTA); both rank the correct transformation gene #1. |
| 68 | GENERATE | R4-G5 | Fable 5 | cmm.features.strain_design (OptKnock + RobustKnock via straindesign + guaranteed-product eval); succinate designs found in <1s. |
| 69 | VERIFY(pending) | R4 | Sonnet (Explore) ×2 | Independent verifiers: (conditions+transformation), (OptKnock/RobustKnock). |
| 70 | VERIFY(PASS) | R4 | pytest | 90 tests (was 72: +18), ruff clean, layering clean (Qt/matplotlib/straindesign out of core/omics). |
| 71 | VERIFY(PASS) | R4-G3/G4 | Sonnet (Explore) | JS66 both methods optimal; MOMA improvement +8.376 independently confirmed; MTA g2 #1; non-mutation/gating ok. Minor: pseudocount=0 div-zero. |
| 72 | VERIFY(PASS) | R4-G5 | Sonnet (Explore) | Independent max/guaranteed product match exactly (9.9149/9.9013); KOs biologically sensible; non-mutation/gating/determinism ok. Medium: uncouplable product crash. |
| 73 | FIX | R4 | Fable 5 | Guard flux_log_change pseudocount=0 (R-edge) + uncouplable-product pre-check in strain design (R17) + matplotlib colormap deprecation. |
| 74 | VERIFY(PASS) | R4 | pytest | 92 tests, ruff clean, 1 warning (expected). Both verifier findings resolved + locked with tests. |
| 75 | CLOSE | loop-R4 | Fable 5 (lead) | Round 4 complete: solver status, flux-map nodes, multi-condition omics+logchange, A→B transformation (MOMA/MTA), OptKnock/RobustKnock. |
| 76 | START | R5 | Fable 5 | User: media presets, pFBA, MOMA/ROOM templates, CMM branding, publication README, GUI-vs-cobra validation, scrub JS66. |
| 77 | GENERATE | R5-G1/G2 | Fable 5 | cmm.core.media (presets + apply via cobra medium API) + pfba in core.simulation; GUI medium selector + Run pFBA. |
| 78 | GENERATE | R5-G3 | Fable 5 | reference_flux(fba/pfba/lad/eflux2) template builder + GUI Comparison tab (MOMA/ROOM + template + KO); MOMA/ROOM now return infeasible gracefully on lethal KO. |
| 79 | GENERATE | R5-G4 | Fable 5 | Branding: CMM = "Cellular Metabolic Modeling Workbench" (GUI title/header + README). |
| 80 | GENERATE | R5-G5 | Fable 5 | Publication README (overview/install/tutorial/screenshots via e_coli_core) + python -m cmm.app entry; scrubbed ALL JS66 from src/README/figures; genericized genome_scale_scenario to a model path. |
| 81 | GENERATE | R5-G6 | Fable 5 | tests/test_validation.py: FBA/pFBA/FVA/MOMA/yield/medium == direct cobra on e_coli_core. |
| 82 | VERIFY(PASS) | R5 | pytest | 112 tests (was 92: +20), ruff clean, layering clean, JS66 absent from published files. |
| 83 | START | R6 | Fable 5 | User: add FVSEOF (reference method). |
| 84 | GENERATE | R6 | Fable 5 | fvseof + FvseofResult (FVA per enforced level → mean + forced-min |flux|; robust targets) + fvseof_figure + GUI button. |
| 85 | RULE/FIX | R6 | Fable 5 | R18: cobra FVA processes=1 in the per-step loop (68 s → 0.2 s). |
| 86 | VERIFY(pending) | R6 | Sonnet (Explore) | Independent FVSEOF verifier (science vs FSEOF/published reference). |
| 87 | VERIFY(PASS) | R6 | pytest | 117 tests (was 112: +5), ruff clean, layering clean. FVSEOF finds FRD7/FUM/MDH/PPC; robust = forced-min rises. |
| 88 | VERIFY(PASS) | R6 | Sonnet (Explore) | mean/forced exact to 1e-8 (incl. negative MDH); FVSEOF∩FSEOF amplify 15/17; robust⊆amplify; non-mutation; 0.15 s. Minor: zero-yield level/column shape mismatch. |
| 89 | FIX | R6 | Fable 5 | enforced_levels uses actual scan keys (zero-yield shape consistency) + test. 118 tests pass. |
| 90 | CLOSE | loop-R6 | Fable 5 (lead) | FVSEOF complete: FVA-per-step robust amplification targets + figure + GUI + validated vs FSEOF/cobra. |
