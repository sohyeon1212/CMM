# CMM scenarios

Step-by-step metabolic-engineering pipelines for driving CMM from an AI coding CLI. Each
scenario takes a goal and produces a report, publication figures, and raw data.

Read `AGENTS.md` first for the router, the solver gate, and the rules. Read
`docs/agent-reference.md` for signatures while writing the calls.

## Index

| ID | Goal | Requires | Minimum solver | Key outputs |
|---|---|---|---|---|
| [SC-01](SC-01-production-target-discovery.md) | Increase production of a target metabolite, and design a strain where it is guaranteed | model, product exchange | LP (design step needs **MILP** + `straindesign` + Java) | growth-coupled knockout design + amplification targets, verified |
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
                       step 3  design   (inverse, MILP)
                       step 5a check    (forward, MOMA/ROOM)
```

- **SC-01 is the spine of a production goal**, and it now carries the coupling design itself.
  Step 3 searches knockout sets with OptKnock/RobustKnock and proves coupling through
  `guaranteed_product`; step 5a applies the design and predicts the strain's immediate
  phenotype with MOMA/ROOM. Inverse design and forward check live in one pipeline because they
  answer different questions — see SC-01's *Do not* for what must never be merged.
- **A coupling-only request enters SC-01 at step 3** and skips the amplification half.
- **Without MILP**, step 3 falls back to a single-deletion screen and the report must say that
  coupling was not established. SC-01 spells out the three required disclosures.
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
