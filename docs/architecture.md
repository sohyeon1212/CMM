# Architecture

CMM separates numerical services from the desktop interface. The Python API and GUI call the
same functions, and result objects carry enough metadata to reproduce a run.

## Runtime layers

- `cmm.core` owns conditions, media, FBA/pFBA/FVA, immutable flux states, solver-capability
  checks, common ranking types, and deterministic provenance.
- `cmm.omics` maps expression through GPR rules, runs E-Flux2 or LAD, predicts multiple
  conditions, and derives source-to-target reaction directions.
- `cmm.features` owns perturbation resolution, MOMA/ROOM, theoretical yield and production
  scans, OptKnock/RobustKnock, MTA/rMTA, and A→B transformation ranking.
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
| LP | FBA, pFBA, FVA, LAD, yield, envelope, FSEOF, FVSEOF |
| MILP | ROOM; OptKnock/RobustKnock through StrainDesign |
| QP | L2 MOMA, E-Flux2, explicitly named `rmta_continuous` heuristic |
| MIQP | published MTA and published rMTA |

GLPK supplies LP and MILP. Gurobi and CPLEX supply all four classes; their license limits
still determine feasible model size. A method never silently changes its formulation when a
capability is missing. In particular, E-Flux2 raises `SolverCapabilityError` unless its
explicitly named `allow_l1_fallback=True` approximation is requested.

## Scientific boundaries

Implemented and tested services are enumerated in `cmm.features.INCLUDED_FEATURES`. Dynamic
FBA, random flux sampling, flux-response analysis, and enzyme-constrained modeling remain
roadmap items; they are not exposed as shipped capabilities.

The following distinctions are intentional:

- OptKnock uses the optimistic two-level formulation; RobustKnock uses the distinct
  three-level worst-case formulation.
- `rmta` is the published best/MOMA/worst pipeline and Equation 9; the historical continuous
  approximation is available only as `rmta_continuous`.
- FSEOF ranks a single biomass-optimal flux at each enforced product level. FVSEOF performs
  FVA at each level and reports midpoint, forced-minimum magnitude, and range-width trends.
- Boundary reactions, biomass, the target exchange, and reactions without a GPR are retained
  in diagnostic tables but excluded from actionable target lists by default.

See [VALIDATION.md](VALIDATION.md) for reference equations, test evidence, reproducibility
commands, and limitations.
