# CMM scenarios

Step-by-step metabolic-engineering pipelines for driving CMM from an AI coding CLI. Each
scenario takes a goal and produces a report, publication figures, and raw data.

Read `AGENTS.md` first for the router, the solver gate, and the rules. Read
`docs/agent-reference.md` for signatures while writing the calls.

## Index

| ID | Goal | Requires | Minimum solver | Key outputs |
|---|---|---|---|---|
| [SC-01](SC-01-production-target-discovery.md) | Increase production of a target metabolite | model, product exchange | LP | ranked amplification + knockout targets, verified |
| [SC-02](SC-02-growth-coupled-strain-design.md) | Make production obligatory at maximum growth | model, product exchange | **MILP** + `straindesign` + Java | knockout sets with guaranteed product |
| [SC-03](SC-03-omics-context-engineering.md) | Explain and exploit a difference between conditions | model, expression table | LP (LAD) | per-condition fluxes, ranked differences, context targets |
| [SC-04](SC-04-knockout-screening.md) | Screen every single deletion | model | LP (`moma_l1`) | essentiality classes, beneficial deletions |

Every scenario is complete on its own: each answers its own question and ends with its own
report. None is a prerequisite for another.

Shared, used by all of them:

- [`_preflight.md`](_preflight.md) — checks that run before any scenario. Skipping these is how
  a run produces a confident wrong answer instead of an error.
- [`_reporting.md`](_reporting.md) — the artifact contract every scenario ends with: directory
  layout, provenance, raw data, figures, and the `report.md` skeleton.

## How they fit together

Each box below is an entry point. Start at whichever one matches the user's goal — the arrows
are optional continuations, not required steps.

```
              every scenario: _preflight … _reporting

   ┌──────────────┐                                ┌──────────────┐
   │    SC-03     │                                │    SC-04     │
   │ explain a    │                                │ essentiality │
   │ difference   │                                │ + screening  │
   └──────┬───────┘                                └──────┬───────┘
          │  supplies the condition       supplies screened │
          │  to search in                 knockout candidates│
          └──────────────►┌──────────────┐◄─────────────────┘
                          │    SC-01     │
                          │  production  │
                          │   targets    │
                          └──────┬───────┘
                                 │ candidates worth proving
                                 ▼
                          ┌──────────────┐
                          │    SC-02     │
                          │  guaranteed  │
                          │   coupling   │
                          └──────────────┘
```

- **SC-01 is the spine of a production goal.** When the ask is "make more of X" and nothing
  narrower, start and finish here.
- **SC-02 is its own goal — proving coupling.** Its natural place is after SC-01 (screen
  cheaply, then design rigorously), but it needs nothing from SC-01: a model, a product, and a
  MILP solver are enough to start here.
- **SC-03 and SC-04 are complete studies, not sub-steps.** "Why does condition A outproduce B"
  and "which genes are essential" are finished answers. They *also* compose with SC-01 when the
  goal is production — SC-03 chooses the condition, SC-04 supplies candidates — but neither is
  subordinate to it, and neither has to be run before SC-01.

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
