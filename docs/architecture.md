# Architecture

CMM is organized around simulation workflows instead of inherited desktop menus.

## Layers

- `cmm.core`: conditions, constraints, objectives, the `FluxState` reference primitive,
  solver capability detection, and solver-neutral result objects (including `TargetRanking`)
- `cmm.io`: SBML import/export and the new project archive format
- `cmm.visualization`: map state, flux coloring, and interactive slider state
- `cmm.features`: workflow services for advanced analysis methods
- `cmm.app`: desktop shell and view models

Only `cmm.app` should depend on Qt. Core and feature modules should be testable without
a GUI process.

## Solver Capability

The cobra default solver is LP-only. Each service declares the solver class it needs (LP,
QP, or MILP/MIQP) and checks it at call time via `core.solvers`. LP and QP run on an open
stack; MILP/MIQP features degrade to a typed `SolverCapabilityError` with an actionable
message when no commercial solver is configured. This keeps the open install usable while
documenting which features need CPLEX/Gurobi.

## Reference Flux State

`core.flux_state.FluxState` is a named, serializable flux vector with provenance (pFBA,
sampling mean, or imported). It is the shared reference for MOMA, ROOM, and revert-
metabolism, so distance-to-reference logic lives in one place.

## Project Model

The primary editable state is a condition:

- name
- reaction bound overrides
- objective coefficients and direction
- optional notes

The condition model intentionally excludes scenario templates, separate scenario files,
media-management state, and legacy annotation conventions.

## Included Feature Boundaries

- FBA/FVA: core services backed by cobra model objects
- Flux visualization slider: maps normalized flux ranges to a rendering state
- Dynamic FBA: time-course service with explicit exchange/update callbacks
- MOMA/ROOM: perturbation services that compare reference and mutant states
- Batch MOMA/ROOM: job runner over perturbation tables
- Random sampling: sampler service returning typed flux sample tables
- FSEOF/FVSEOF: objective-enforced scanning services
- Omics integration: expression table normalization plus method-specific constraints
- Enzyme-constrained modeling: separate model extension layer and import/export helpers
- Flux response analysis: parameter sweep service over chosen bounds or objectives
- OptKnock/RobustKnock: new design services with a dedicated product wizard
- Revert metabolism: normalization-target service that ranks KOs by how well they move a
  source `FluxState` toward a target state, using differential expression direction labels

## Removed Feature Boundaries

The desktop app will not include EFM, thermodynamics, legacy MCS workflows, scenario
templates, media management, or separate scenario file formats.
