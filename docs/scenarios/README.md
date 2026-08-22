# CMM scenarios

Step-by-step metabolic-engineering pipelines for driving CMM from an AI coding CLI. Each
scenario takes a goal and produces a report, publication figures, and raw data.

Read `AGENTS.md` first for the router, the solver gate, and the rules. Read
`docs/agent-reference.md` for signatures while writing the calls.

For a complete production request, use the auto-discoverable
`.agents/skills/cmm-production-engineering/` skill and the canonical
`cmm production-targets --config CONFIG` workflow. The detailed steps below explain the
scientific roles and support narrow API calls; they are not a reason to invent a second,
one-off orchestrator. To change workflow parameters or define a new research workflow, follow
[Building a reproducible CMM workflow](../building-custom-workflows.md).

## Index

| ID | Goal | Requires | Minimum solver | Key outputs |
|---|---|---|---|---|
| [SC-01](SC-01-production-target-discovery.md) | Increase production of a target metabolite, and design a strain where it is guaranteed | confirmed model, product exchange, condition | **QP + MILP** for the full workflow; strain design also needs importable `straindesign` | MOMA/ROOM single-deletion candidates, OptKnock/RobustKnock designs, FSEOF/FVSEOF targets, forward validation |
| [SC-02](SC-02-omics-context-engineering.md) | Explain and exploit a difference between conditions | model, expression table | LP (LAD) | per-condition fluxes, ranked differences, context targets |
| [SC-03](SC-03-knockout-screening.md) | Screen every single deletion | model | LP (`moma_l1`) | essentiality classes, beneficial deletions |

Every scenario is complete on its own: each answers its own question and ends with its own
report. None is a prerequisite for another.

Shared, used by all of them:

- [`_preflight.md`](_preflight.md) — checks that run before any scenario. Skipping these is how
  a run produces a confident wrong answer instead of an error.
- [`_reporting.md`](_reporting.md) — the artifact contract every scenario ends with: directory
  layout, provenance, raw data, figures, and the `report.html` skeleton.

## How they fit together

Each box below is an entry point. Start at whichever one matches the user's goal — the arrows
are optional continuations, not required steps.

```
              every scenario: _preflight … _reporting

   ┌──────────────┐                                ┌──────────────┐
   │    SC-02     │                                │    SC-03     │
   │ explain a    │                                │ essentiality │
   │ difference   │                                │ + screening  │
   └──────┬───────┘                                └──────┬───────┘
          │  supplies the condition        supplies an exhaustive│
          │  to search in                   single-deletion view │
          └──────────────►┌──────────────┐◄─────────────────┘
                          │    SC-01     │
                          │  production  │
                          │    strain    │
                          └──────────────┘
                       step 3  single KO (forward, MOMA/ROOM)
                       step 4  design    (inverse, MILP)
                       step 6  checks    (loopless FVA + response + sampling)
```

- **SC-01 is the spine of a production goal.** Its single-knockout step evaluates every
  nominated deletion independently with MOMA and ROOM; its separate strain-design step searches
  knockout sets with OptKnock/RobustKnock and proves coupling through `guaranteed_product`.
  These outputs answer different questions and remain separate tables and figure panels.
- **A coupling-only request enters SC-01 at the strain-design step** and skips single-deletion
  and amplification analyses unless the user asks for them.
- **The canonical full workflow does not silently downgrade.** Missing QP, MILP,
  importable `straindesign`, R, or a required R package is a surfaced capability result. Any
  additional backend requirement is surfaced by that backend. A user may
  explicitly approve a narrower run, but its report must state which claim is no longer made.
- **SC-02 and SC-03 are complete studies, not sub-steps.** "Why does condition A outproduce B"
  and "which genes are essential" are finished answers. They *also* compose with SC-01 when the
  goal is production — SC-02 chooses the condition, SC-03 gives the exhaustive deletion picture
  — but neither is subordinate to it, and neither has to be run before SC-01.

Only chain scenarios when the user's goal needs the second one, and say in the report why you
ran it.

## Scenario document structure

Every `SC-*.md` uses the same shape so a pipeline can be followed without reading prose:

- **YAML front matter** — id, goal, when to use (including Korean phrasings), required and
  optional inputs, solver requirements, step list, expected runtime.
- **Objective** — what the run produces, plus measurable success criteria.
- **Pipeline at a glance** — one row per step: question, method, output.
- **Steps** — each with *Goal*, *Preconditions*, *Call* (runnable code), *Outputs*,
  *Artifacts*, *Decision rule*, *Branch*, *Failure → action*, and *Solver*.
- **Cross-checks** — internal consistency to confirm before writing the report.
- **Do not** — the misuses that produce plausible, wrong conclusions.

The *Decision rule* fields carry numeric thresholds. They are documented starting points, not
community standards: state whichever you used in the report.
