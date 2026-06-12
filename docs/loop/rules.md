# Failure-Derived Rules (Strategy 2: continual-learning memory)

Append-only. Each rule is a *verified fact* generalized from a real failure, so later
stories reference the rule instead of rediscovering the cause. Format:

`R<n>` — <rule> — *(from: <story/event that produced it>)*

## Environment rules (seeded from recon, verified by command output)

- R1 — Use `.venv/bin/python` (uv) for everything, never anaconda python: anaconda's
  `libQt5Core` collides with PyQt5's bundled Qt and the `offscreen` plugin fails to load.
  *(from: recon offscreen grab crash on anaconda; clean grab on .venv)*
- R2 — QP/MILP/MIQP only exist in `.venv` (Gurobi academic + CPLEX). anaconda has only
  glpk(LP)+scipy. Gate QP/MILP features on `cmm.core.solvers` and provide an LP fallback
  where one exists (e.g. linear MOMA). *(from: recon solver enumeration)*
- R3 — Capture GUI screenshots with `QT_QPA_PLATFORM=offscreen` + `QWidget.grab().save()`;
  it produces real rendered pixmaps headlessly. *(from: recon grab produced a 320x80 PNG)*

## Build rules

- R4 — cobra 0.31 removed `find_gene_knockout_reactions`; use
  `cobra.manipulation.delete.knock_out_model_genes(model, gene_list)`, which returns the
  reactions it zeroed and reverts inside a `with model:` context (so call it there to read
  blocked reactions without mutating the caller's model). Keep the old name only as an
  ImportError fallback. *(from: G2 VERIFY FAIL — ImportError on test_comparison collection)*
- R5 — When adapting a `FluxState` into a cobra `Solution` for ROOM, set
  `objective_value` to the model objective evaluated at the reference fluxes
  (`sum(r.objective_coefficient * flux)`), never a placeholder 0.0: cobra's `add_room`
  constrains the original objective to `solution.objective_value` (the `room_old_objective`
  constraint), so 0.0 silently forces the whole network to zero flux and ROOM "succeeds"
  with a meaningless all-zero distribution. *(from: G2 VERIFY FAIL — ROOM returned all-zero
  fluxes; biomass-floor variant went infeasible)*
- R6 — When adding optlang variables AND constraints that reference them, add the variables
  in their own `add_cons_vars([...])` call followed by `model.solver.update()`, THEN add the
  constraints; adding both in one call (or a constraint referencing a not-yet-flushed var)
  makes optlang's gurobi interface re-queue the variable and raise
  `ContainerAlreadyContains` on the next update. For one-off MIQP builds, do this on
  `model.copy()` to also avoid cross-iteration name collisions and history-context conflicts.
  Prefer cobra's built-in moma/room over hand-rolled optlang plumbing where one exists.
  *(from: G4 — repeated ContainerAlreadyContains on mta_z_* during MIQP build)*
- R7 — Theoretical yield must FIX the substrate uptake exactly (`bounds = (-U, -U)`) and use
  U as the denominator. Maximizing a product exchange has alternate optima with variable
  substrate flux and CO2 co-fixation, so dividing by the *actual* substrate flux is
  non-deterministic (gave 2.96 instead of 1.638 for aerobic succinate on e_coli_core).
  *(from: R2-G1 VERIFY FAIL — aerobic succinate yield mismatch)*
- R8 — Carbon-substrate auto-detection must pick the *limiting* uptake exchange (e.g. glucose
  at lb=-10), not the most-negative bound: CO2/H2O/NH4/O2 sit at the default ~-1000 and would
  win a naive "most negative lb" search. Filter to explicitly-restricted uptakes
  (-999 < lb < 0); only then fall back to actual WT uptake flux.
  *(from: R2-G1 VERIFY FAIL — _detect_substrate returned EX_co2_e, forcing infeasibility)*
- R10 — FSEOF (and any flux-trend classifier) must classify on flux MAGNITUDE (`abs`), not
  signed flux: a reaction operating in reverse (negative flux, e.g. the reductive succinate
  arm MDH/FUM) has rising |flux| but falling signed flux, so signed comparison labels it
  "knockdown" — exactly backwards. *(from: R2-G5 adversarial review CRITICAL — amplify/knockdown inverted)*
- R11 — Set a cobra reaction's bounds atomically (`rxn.bounds = (lo, hi)`), never
  `lower_bound=` then `upper_bound=`: cobra validates each assignment independently, so a new
  value that crosses the current opposite bound raises mid-assignment. Same for GUI bound
  edits. *(from: R2-G5 adversarial review CRITICAL — set_reaction_bounds crash on crossing bounds)*
- R12 — Theoretical yield above the substrate carbon ceiling (e.g. succinate 1.64 > 1.5
  mol/mol glucose) is real but only via net CO2 fixation; always expose the ceiling + CO2
  exchange so the number is not read as carbon-from-substrate alone. And exclude currency
  metabolites (ATP, H2O, CO2, NAD(P)H, ...) from network flux maps or they become a hairball.
  *(from: R2-G5 adversarial review HIGH/CRITICAL — undisclosed CO2 yield; network-map gap)*
- R13 — To match curated Escher metabolic-map quality, do NOT invent a force-directed layout (it
  hairballs). Curated Escher maps are JSON whose nodes carry hand-laid (x, y)
  coordinates and whose reactions carry bezier `segments`. Render THAT layout coloured by
  flux (`escher_flux_map`): identical biochemical layout, plus flux encoding. The map file is
  user/BiGG data passed in by path — CMM bundles no maps. *(from: user — "왜 map이
  안떠"; force-layout network map was a hairball, Escher-layout render matches the curated map)*
- R14 — On a genome-scale model, the carbon source is whatever the *medium* allows: JS66's
  default medium has `EX_glc__D_e` lower_bound = 0 (no glucose), so auto-detect correctly
  returns the available carbon exchange (CO2). To analyse succinate-from-glucose you must open
  glucose uptake first (a bound edit) — then everything works (yield 1.5, FSEOF finds FRD2).
  Always surface the detected substrate and let the user override it. *(from: user — test the
  big JS66 model; substrate "surprise" traced to the model's medium, not a code bug)*
- R15 — E-Flux2's objective-floor constraint must branch on objective direction: hold the
  biological objective with a LOWER bound for a maximized objective (`lb = f*opt`) but an
  UPPER bound for a minimized one (`ub = opt/f`). Always using `lb=` silently breaks
  min-direction models. *(from: omics verifier — latent floor-direction bug)*
- R18 — When calling cobra `flux_variability_analysis` in a loop (e.g. FVSEOF runs one FVA
  per enforced level), pass `processes=1`. The default spawns a multiprocessing pool *per
  call*; the pool-spawn overhead dominated FVSEOF at ~68 s for 6 steps — `processes=1` cut it
  to ~0.2 s. *(from: FVSEOF — per-step FVA pool-spawn overhead)*
- R16 — OptKnock/RobustKnock: delegate the bilevel MILP search to `straindesign`
  (SDModule OPTKNOCK), then EVALUATE each design yourself at MAXIMUM growth — fix biomass to
  its optimum, then max product = optimistic (OptKnock), min product = guaranteed worst case
  (RobustKnock). RobustKnock = keep designs with guaranteed product > 0. straindesign is
  chatty: wrap the solve in redirect_stdout/stderr. *(from: R4-G5 — succinate coupling
  designs on e_coli_core)*
- R9 — Under growth-maximizing FBA, wild-type e_coli_core excretes NO succinate (it ferments
  to acetate/ethanol), so "go anaerobic" alone does not raise succinate. To DEMONSTRATE
  succinate production you must block the competing fermentation secretions — and "block
  secretion" means setting the exchange UPPER bound to 0 (e.g. EX_ac_e=(0,0)), not the lower
  bound. Anaerobic + block ac/etoh/for/lac forces succinate 0 -> ~9.9 via FRD7.
  *(from: R2-G4 self-review — scenario claimed a succinate increase that was actually 0->0)*
