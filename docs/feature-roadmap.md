# Feature Roadmap

## Cross-Cutting Prerequisites

These are not a phase; they are dependencies that several phases assume and that the
original plan left implicit. They must land before the features that need them.

- **Solver tiers.** Declare solver capability explicitly instead of relying on the cobra
  default (GLPK, LP-only). FBA/FVA/FSEOF need LP. MOMA, ROOM, and revert-metabolism need
  QP. OptKnock/RobustKnock and original MTA (MIQP) need MILP/MIQP. Plan: ship LP + QP on
  an open stack (HiGHS/OSQP via optlang), gate MILP/MIQP behind an optional CPLEX/Gurobi
  install, and have each service raise a typed `SolverCapabilityError` when the active
  solver cannot run it.
- **Reference flux state.** Add a `core.flux_state.FluxState` primitive (a named,
  serializable flux vector plus provenance: pFBA, sampling mean, or imported). MOMA, ROOM,
  and revert-metabolism all compare a candidate state against a reference; without this
  type each feature would reinvent it.
- **Differential expression.** The omics layer is currently scoped to single-state
  constraint methods (E-Flux2 / LAD). Revert-metabolism needs two-state (source vs target)
  differential expression mapped to a per-reaction desired direction. Add this to the omics
  layer so it is shared, not buried in one feature.
- **Target ranking result.** FSEOF, OptKnock, and revert-metabolism all emit a ranked list
  of intervention targets with scores. Add one `core.results.TargetRanking` type they share.

## Phase 1: Core Engine

- Add SBML load/save services
- Add condition apply/serialize/deserialize
- Add FBA and FVA
- Add the `FluxState` reference primitive (pFBA and imported sources first)
- Add deterministic result export tables

## Phase 2: Visualization

- Add map data model
- Add flux-to-color normalization
- Add slider state for flux range inspection
- Add desktop map view after core rendering state is tested

## Phase 3: Advanced Simulation

- Add dynamic FBA
- Add MOMA and ROOM (consume `FluxState` as the reference)
- Add batch perturbation runner (gene/reaction KO enumeration + parallel solve)
- Add random flux sampling (also feeds `FluxState` sampling-mean references)
- Add flux response analysis

## Phase 4: Engineering Workflows (Production Targets)

- Add FSEOF and FVSEOF
- Add OptKnock and RobustKnock as new implementations
- Emit results as `TargetRanking`

## Phase 5: Omics and Enzyme Constraints

- Add omics table import and normalization
- Add single-state omics-to-constraint methods (E-Flux2 / LAD class)
- Add two-state differential expression -> per-reaction desired direction (source vs target)
- Add enzyme-constrained model extension layer
- Add enzyme-constrained import/export helpers and validation

## Phase 6: Normalization Targets (Revert Metabolism)

New capability with no analogue in the legacy app or the original plan. Predicts gene/
reaction interventions that move a perturbed (e.g. disease) metabolic state back toward a
reference (e.g. healthy) state. Depends on Phase 3 (reference state + perturbation runner)
and Phase 5 (differential expression).

- Add `features.revert` service implementing the robust MTA (rMTA) continuous/QP scoring
- Add optional original-MTA (MIQP) high-fidelity mode, gated on a MILP/MIQP solver
- Reuse the batch perturbation runner for KO enumeration
- Emit ranked normalization targets as `TargetRanking`
- See `docs/design-revert-metabolism.md` for the full design

## Phase 7: Desktop Product

- Build the Qt shell with Project, Model, Simulation, Design, and View sections
- Add a compact model browser
- Add the visualization workspace
- Add task-specific dialogs only after the service APIs are stable
