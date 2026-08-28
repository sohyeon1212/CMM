# Architecture

CMM separates numerical services from the desktop interface. The Python API and GUI call the
same functions, and result objects carry enough metadata to reproduce a run.

## Runtime layers

- `cmm.core` owns conditions, media, FBA/pFBA/FVA, immutable flux states, solver-capability
  checks, common ranking types, and deterministic provenance.
- `cmm.omics` maps expression through GPR rules, runs E-Flux2 or LAD, predicts multiple
  conditions, and derives source-to-target reaction directions.
- `cmm.features` owns perturbation resolution, MOMA/ROOM, theoretical yield and production
  scans, flux-response scans, random and reference-constrained flux sampling,
  OptKnock/RobustKnock, MTA/rMTA, and A→B transformation ranking.
- `cmm.visualization` converts already-computed results into matplotlib figures. It does not
  solve metabolic models.
- `cmm.workflows.production` composes existing public services into the fixed SC-01 sequence,
  writes typed exports and provenance, and owns the `01_preflight`–`07_validation` run schema.
  It contains no second implementation of the scientific methods and adds no private ranking
  rule.
- `cmm.reporting` currently reads a completed SC-01 run, renders its Nature Genetics-inspired
  R figure/report layer, and validates production-specific artifact and claim coverage. It never
  invokes a metabolic solver or changes numerical results.
- `cmm.cli` is the designated thin adapter for `production-targets`, `report render`, and
  `report validate`: those commands currently target SC-01 and may perform config parsing and
  exit-code mapping only.
- `cmm.app` is the Qt shell. It validates files and UI state, dispatches long analyses to a
  worker, and renders service results.

Only `cmm.app` depends on Qt. The scientific services, workflow, CLI, and validator are
importable and testable in a headless process. The `nature-r` renderer depends on an external
`Rscript` process and fails before rendering when a declared package is unavailable.

SC-01 and SC-04 are the shipped canonical workflows. Other analyses are recipes
over the public services, not installed workflow APIs, CLI commands, or validated schemas. A
contributor adding another canonical workflow must define a separate schema id/version,
renderer, validator, and public boundary rather than treating production-specific helpers as a
generic engine. The complete extension sequence is in
[Adding a canonical workflow to CMM](tutorials/adding-a-canonical-workflow.md).

## State and data flow

```text
confirmed model path + product + explicit condition
              │
              ▼
       production workflow ──────► core / feature services
              │                              │
              │                              ▼
              │                      typed numerical results
              ▼
       CSV + JSON + manifest ─────► R report renderer
                                             │
                                             ▼
                                  HTML + PNG + PDF/SVG
                                             │
                                             ▼
                                       run validator
```

For a narrow analysis, callers bypass the workflow and call the relevant core, omics, or
feature service directly. The workflow is composition and artifact ownership, not a new
numerical method.

`FluxState` is the shared complete reaction-flux vector used by MOMA, ROOM, MTA/rMTA, and
transformation analyses. It rejects empty or non-finite state vectors and records its
origin. A source state must be regenerated after model, medium, bound, or expression changes.
The GUI clears these derived states whenever a new model is loaded.

## Solver contracts

Every method checks the mathematical capability it actually requires:

| Capability | Methods |
|---|---|
| LP | FBA, pFBA, FVA, LAD, yield, envelope, FSEOF, FVSEOF, flux response, flux sampling |
| MILP | ROOM; OptKnock/RobustKnock through StrainDesign |
| QP | L2 MOMA, E-Flux2, explicitly named `rmta_continuous` heuristic |
| MIQP | published MTA and published rMTA |

GLPK supplies LP and MILP. Gurobi and CPLEX supply all four classes; their license limits
still determine feasible model size. A method never silently changes its formulation when a
capability is missing. In particular, E-Flux2 raises `SolverCapabilityError` unless its
explicitly named `allow_l1_fallback=True` approximation is requested.

The canonical production workflow strengthens that local rule into a run-level gate. Its
matched single-knockout stage requires QP for MOMA-L2 and MILP for ROOM; strain design checks
MILP plus importable `straindesign` and surfaces additional backend requirements. Reporting
independently checks `Rscript` and required renderer-package availability. Package metadata
provides compatible minimum versions at runtime; exact versions come from `renv.lock` and the
three-OS CI comparison. A user may request a narrower run, but an orchestrator may not silently
decide to produce one.

## Scientific boundaries

Implemented and tested services are enumerated in `cmm.features.INCLUDED_FEATURES`. **Flux
sampling (`random_flux_sampling`, `reference_constrained_sampling`) and flux-response analysis
(`flux_response`) are shipped**, with method contracts in `VALIDATION.md`, GUI tabs, and
publication figures. The concrete `production_target_workflow` and `publication_reporting`
services are also shipped; their presence does not imply a generic scenario-template engine,
which remains excluded. Dynamic FBA and enzyme-constrained modeling remain roadmap items.
`docs/feature-roadmap.md` holds the current split.

The following distinctions are intentional:

- OptKnock uses the optimistic two-level formulation; RobustKnock uses the distinct
  three-level worst-case formulation.
- `rmta` is the published best/MOMA/worst pipeline and Equation 9; the historical continuous
  approximation is available only as `rmta_continuous`.
- FSEOF ranks a single biomass-optimal flux at each enforced product level. FVSEOF performs
  FVA at each level and reports midpoint, forced-minimum magnitude, and range-width trends.
  Their top candidate sets are retained independently; overlap is descriptive evidence, not a
  workflow selection gate.
- Boundary reactions, biomass, the target exchange, and reactions without a GPR are retained
  in diagnostic tables but excluded from actionable target lists by default.
- `flux_response` reports sensitivity as the exact LP shadow price and its phase boundaries.
  The finite-difference "bottleneck" it previously reported was removed in 0.4.0 as
  grid-dependent and unpublished, and `feasible_range()` is FVA-derived rather than read off
  the scan grid.
- Solver-backed services return serializable typed result containers carrying
  `run_provenance` and stable frame exports. `fva` returns the mapping-compatible `FvaResult`,
  while `BatchComparisonResult` preserves its legacy list behavior and adds metadata plus
  `to_frame()`. Lightweight arithmetic helpers such as `flux_log_change` and `sign_flips`
  intentionally return dictionaries or lists rather than pretending to be solver runs.
  `Medium.apply_to` returns `MediumApplication`, which records components the loaded model could
  not express instead of dropping them silently.
- One convention states the environment of a run: `condition=`. The `aerobic=True|False`
  parameter was removed from `cmm.features.production` in 0.4.0 and `optknock`/`robustknock`
  gained `condition=`, so no service depends silently on the caller's model state.
- Sampling defaults to `processes=1` so a seeded run is bit-for-bit reproducible; parallel
  chains are seeded independently and provenance records which you got.

See [VALIDATION.md](VALIDATION.md) for reference equations, test evidence, reproducibility
commands, and limitations.
