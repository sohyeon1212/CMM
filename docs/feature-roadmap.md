# Feature roadmap

This roadmap starts from the validated CMM 0.3.0 surface. Implemented features are also
enumerated by `cmm.features.INCLUDED_FEATURES`; this document does not advertise planned work
as shipped functionality.

## Current validated surface

- Core: conditions, media, FBA, pFBA, FVA, `FluxState`, solver contracts, provenance.
- Perturbation: gene/reaction/multiple knockout resolution, L1/L2 MOMA, ROOM, batch runs.
- Production: yield, production envelope, FSEOF, FVSEOF, OptKnock, RobustKnock.
- Omics: GPR mapping, LAD, E-Flux2, multi-condition prediction, differential directions.
- Transformation: MOMA A→B ranking, published MTA, published rMTA, explicit continuous
  heuristic.
- Product: Qt desktop workflows, CSV/table export, publication figures, offscreen tests.

The reference tests and scientific limitations are defined in `docs/VALIDATION.md`.

## Publication follow-up

These items improve external scientific evidence without changing the existing method
contracts:

1. Archive a DOI-bearing release and replace the organization-only `CITATION.cff` author
   entry with the manuscript's final authors, ORCIDs, title, and DOI.
2. Publish immutable input model and omics artifacts with checksums rather than relying only
   on upstream model repositories.
3. Add study-specific biological benchmarks with prespecified metrics and held-out
   experimental interventions.
4. Record wall time, peak memory, solver tolerances, and candidate counts for representative
   genome-scale LP, QP, MILP, and MIQP workflows.
5. Add the original MTA/rMTA contextualization plus sampling preprocessing if a manuscript
   claims the complete paper pipeline rather than CMM's documented E-Flux2 variant.

## Planned scientific services

- Dynamic FBA with explicit exchange/update callbacks.
- Reproducible random flux sampling and sampling-mean `FluxState` construction.
- Flux-response analysis over bound and objective sweeps.
- Enzyme-constrained model extensions and validated import/export.
- Automated construction of FVSEOF grouping-reaction constraints from physiological data.

## Engineering follow-up

- Cancellation and progress reporting for long-running GUI jobs.
- Persistent project archives containing model, conditions, expression, result metadata, and
  checksums.
- Solver-tolerance configuration recorded in provenance.
- Optional parallel candidate screens with deterministic ordering and bounded resource use.
