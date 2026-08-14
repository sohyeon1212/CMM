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
- `cmm.app` is the Qt shell. It validates files and UI state, dispatches long analyses to a
  worker, and renders service results.

Only `cmm.app` depends on Qt. The scientific services are importable and testable in a
headless process.

## State and data flow

```text
COBRA model + Condition/expression
              │
              ▼
       core / omics service ──► FluxState or typed result
              │                         │
              ▼                         ▼
       feature service             provenance metadata
              │                 (model SHA-256, solver,
              ▼                  versions, parameters)
       table / ranking / figure
```

`FluxState` is the shared complete reaction-flux vector used by MOMA, ROOM, MTA/rMTA, and
transformation workflows. It rejects empty or non-finite state vectors and records its
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

## Scientific boundaries

Implemented and tested services are enumerated in `cmm.features.INCLUDED_FEATURES`. **Flux
sampling (`random_flux_sampling`, `reference_constrained_sampling`) and flux-response analysis
(`flux_response`) are shipped**, with method contracts in `VALIDATION.md`, GUI tabs, and
publication figures. Dynamic FBA and enzyme-constrained modeling remain roadmap items; they are
not exposed as shipped capabilities. `docs/feature-roadmap.md` holds the current split.

The following distinctions are intentional:

- OptKnock uses the optimistic two-level formulation; RobustKnock uses the distinct
  three-level worst-case formulation.
- `rmta` is the published best/MOMA/worst pipeline and Equation 9; the historical continuous
  approximation is available only as `rmta_continuous`.
- FSEOF ranks a single biomass-optimal flux at each enforced product level. FVSEOF performs
  FVA at each level and reports midpoint, forced-minimum magnitude, and range-width trends.
- Boundary reactions, biomass, the target exchange, and reactions without a GPR are retained
  in diagnostic tables but excluded from actionable target lists by default.
- `flux_response` reports sensitivity as the exact LP shadow price and its phase boundaries.
  The finite-difference "bottleneck" it previously reported was removed in 0.4.0 as
  grid-dependent and unpublished, and `feasible_range()` is FVA-derived rather than read off
  the scan grid.
- Every analysis returns a frozen dataclass carrying `run_provenance`, with no exceptions as of
  0.4.0: `fva` returns `FvaResult` rather than a bare `dict[str, FluxRange]` (it is still a
  `Mapping`, so callers are unaffected), and `Medium.apply_to` returns `MediumApplication`,
  which records the components the loaded model could not express instead of dropping them
  silently.
- One convention states the environment of a run: `condition=`. The `aerobic=True|False`
  parameter was removed from `cmm.features.production` in 0.4.0 and `optknock`/`robustknock`
  gained `condition=`, so no service depends silently on the caller's model state.
- Sampling defaults to `processes=1` so a seeded run is bit-for-bit reproducible; parallel
  chains are seeded independently and provenance records which you got.

See [VALIDATION.md](VALIDATION.md) for reference equations, test evidence, reproducibility
commands, and limitations.
