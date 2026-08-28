# Adding a canonical CMM workflow

This tutorial is for contributors who want to add a first-class, goal-level workflow to CMM.
A canonical workflow is more than a script that calls several analyses: it owns a scientific
contract, a typed Python boundary, a deterministic execution sequence, a versioned artifact
schema, a dedicated validator, a publication renderer, a thin CLI, and regression tests.

> **Current shipped boundary**
>
> CMM ships two canonical workflows, and they are worth reading against each other because
> they answer the same structural questions differently:
>
> | | SC-01 production | SC-04 transformation |
> |---|---|---|
> | renderer | R, via `renderer="nature-r"` | pure Python + matplotlib |
> | validator | `validate_production_run` | none yet |
> | reproduce/render/validate scripts | written into `scripts/` | not written |
>
> **Neither is the template.** SC-01 is the fuller reference — it is the one with a validator
> and a scripts directory — while SC-04 shows the smaller shape a workflow can take when its
> panels have no R counterpart. Names such as `MyWorkflowConfig`, `MyWorkflowResult`,
> `run_my_workflow`, `render_my_workflow_report`, and `validate_my_workflow_run` below are
> deliberately generic placeholders. They are not installed, importable, or available as CLI
> commands in the current package.

Do not begin by copying the production module. First decide whether the new question needs a
canonical workflow at all. Numerical work belongs in public solver-neutral services; a workflow
only composes those services and owns their reportable run contract.

## 1. Decide whether the feature is a workflow

Use the smallest boundary that completely answers the scientific question.

| Requested capability | Correct implementation boundary |
|---|---|
| One FBA, FVA, knockout comparison, omics integration, or other single analysis | An existing public service in `cmm.core`, `cmm.features`, or `cmm.omics` |
| A new numerical formulation or scoring method | A new public service in one of those numerical packages |
| The same production study with different thresholds, seeds, or stage switches | New values in `ProductionWorkflowConfig`, not another workflow |
| A repeatable goal that requires several analyses, fixed hand-offs, a stable run directory, figures, a report, and a completion gate | A new canonical workflow |
| A one-off combination used only in one study | A checked-in analysis script plus declared outputs, not necessarily a package workflow |

A new workflow must answer a scientific question that is distinct from SC-01. Different figure
labels, a new organism, a different target product, or a different parameter set do not qualify.
Neither does a generic scenario-template engine: `scenario_templates` remains in
`cmm.features.EXCLUDED_FEATURES`. Add one explicit workflow for one explicit goal rather than a
framework that guesses how arbitrary analyses compose.

Before coding, obtain agreement on all of the following:

- the user question and the claim the workflow may support;
- the intended audience and the decision they will make from the outputs;
- exact inputs and the biological conditions represented by each input;
- stage order and the data passed between stages;
- candidate-universe and ranking rules fixed before results are inspected;
- solver and optional-package capabilities required by every method;
- failure, unavailable, skipped, partial, and lethal-result behavior;
- the raw tables, figures, narrative sections, and replay materials required at completion;
- the validator checks that distinguish a complete run from a merely generated folder;
- claims explicitly outside the workflow, especially biological or wet-lab validation.

If any item remains undecided, the implementation is not ready to start.

## 2. Write the scientific contract first

Add a scenario document under `docs/scenarios/` before implementing orchestration. Use this
copyable outline and replace every bracketed item:

```markdown
# [Workflow title]

## Scientific question
[One question that the workflow answers.]

## Supported claim
[A narrow statement justified by the model and methods.]

## Excluded claims
- [Causality, wet-lab performance, universal essentiality, or another unsupported claim.]

## Inputs
- Exact model: [SBML path; archived byte-for-byte]
- Conditions: [medium, substrate, oxygen/aeration, objective, changed bounds]
- Additional data: [format, identifiers, units, mapping requirements]

## Method roles and order
1. [Preflight]
2. [Reference or baseline]
3. [Inverse discovery, if applicable]
4. [Forward verification, if applicable]
5. [Export, report, and validation]

## Decision rules fixed before execution
- Candidate universe: [complete definition]
- Ranking: [method-specific field, order, and tie behavior]
- Thresholds: [values, units, and boundary behavior]
- Deduplication: [scientific equivalence key and retained aliases]

## Capability contract
- [Method]: [LP/MILP/QP/MIQP plus optional package]
- Missing capability: [fail, or an explicitly requested narrower mode]

## Completion criteria
- [Required artifact roles]
- [Required figures and source tables]
- [Cross-artifact scientific invariants]
- [Dedicated validator returns valid after rendering]
```

Make conditions explicit. A model's current bounds are data, not a universally appropriate
medium. Record the substrate exchange and uptake, oxygen exchange and aeration bounds, objective,
and every changed reaction bound. If the workflow compares conditions, each condition needs its
own complete definition; never change the model in place for one condition and accidentally use
that state as the next condition's baseline.

Separate method roles. Inverse methods propose an intervention from a goal; forward methods
predict the consequence of a supplied intervention. Preserve independent method rankings and
tie blocks unless the protocol defines a scientifically justified combination rule. Do not add a
cross-method score after observing the results.

Define non-optimal behavior before execution. Infeasible and lethal outcomes are scientific
data, not exceptions to delete. A missing solver capability is different: detect it before an
expensive stage, then fail or record the explicitly configured narrower workflow. Never silently
replace a QP, MIQP, or MILP formulation with an LP approximation.

SC-01 demonstrates this separation: numerical methods remain public services, while
`run_production_target_discovery` fixes stage order, capability gates, candidate coverage,
forward validation, and the exported schema. Read its public behavior in
[`cmm.workflows.production`](../../src/cmm/workflows/production.py), but do not treat its private
helpers as a reusable workflow framework.

## 3. Implement missing numerical services before orchestration

The workflow layer must not contain a second implementation of a scientific method. Inventory
the public functions in `cmm.core`, `cmm.features`, and `cmm.omics`, then read the relevant
section of the [agent reference](../agent-reference.md). If a required numerical operation does
not exist, add and validate that service first.

Every new numerical service follows these rules:

1. Accept a `cobra.Model` plus explicit scientific inputs such as `condition=`. Do not depend on
   GUI state or an undocumented mutation previously applied by a caller.
2. Check its own LP, MILP, QP, or MIQP requirement with CMM's solver-capability helpers before
   building the optimization problem.
3. Return a frozen dataclass rather than a bare DataFrame, tuple, or dictionary.
4. Keep row-level statuses, identifiers, values, units where ambiguity is possible, and all
   information needed for interpretation.
5. Provide `to_frame()` for a stable, tidy export. Its columns become part of the workflow's
   artifact contract, so order them explicitly.
6. Store `run_provenance(...)` in a `metadata` field. Include every scientific parameter and
   pass an explicit seed when the method or backend is stochastic. Do not invent a seed for a
   deterministic method merely to make its metadata resemble another service.
7. Export the service and its result types from the owning package's `__init__.py`.
8. Add the method contract and reference evidence to `docs/VALIDATION.md` and unit-test both
   expected results and failure modes.

The following is an illustrative shape only. `MyAnalysisResult` and `run_my_analysis` do not
exist in CMM:

```python
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass

from cobra import Model
import pandas as pd

from cmm.core import Condition, require, run_provenance


@dataclass(frozen=True)
class MyAnalysisRow:
    reaction_id: str
    status: str
    score: float | None


@dataclass(frozen=True)
class MyAnalysisResult:
    rows: tuple[MyAnalysisRow, ...]
    metadata: Mapping[str, object]

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            [asdict(row) for row in self.rows],
            columns=["reaction_id", "status", "score"],
        )


def run_my_analysis(
    model: Model,
    *,
    condition: Condition | None = None,
    threshold: float,
) -> MyAnalysisResult:
    require("LP", model.solver.interface, feature="my analysis")
    with model:
        if condition is not None:
            condition.apply_to(model)
        # Solve through a solver-neutral formulation and retain every status row.
        rows: tuple[MyAnalysisRow, ...] = ()
        metadata = run_provenance(
            model,
            method="my_analysis",
            condition=asdict(condition) if condition is not None else None,
            threshold=threshold,
        )
    return MyAnalysisResult(
        rows=rows,
        metadata=metadata,
    )
```

The empty solve body is intentional: it prevents this tutorial from presenting a planned method
as implemented. A real service must have method-specific tests before a workflow calls it.

## 4. Define the public workflow boundary

Create one module in `src/cmm/workflows/` named after the workflow's stable snake-case slug. Use
one consistent name family:

| Surface | Convention for the tutorial placeholder |
|---|---|
| Module | `cmm.workflows.my_workflow` |
| Schema id | `cmm.my-workflow` |
| Initial schema version | integer `1` |
| Config | `MyWorkflowConfig` |
| Result | `MyWorkflowResult` |
| Domain error | `MyWorkflowError` |
| Runner | `run_my_workflow(config)` |
| Renderer | `render_my_workflow_report(run_dir, renderer="nature-r")` |
| Validator | `validate_my_workflow_run(run_dir)` |
| Analysis CLI | `cmm my-workflow --config CONFIG` |

These are naming conventions, not current APIs. Replace `my_workflow` and `my-workflow` with a
scientifically meaningful name before implementation.

### 4.1 Frozen, serializable configuration

The config is the complete invocation. At minimum, include:

- `model_path`, `output_dir`, and an optional explicit solver;
- all conditions, media, objectives, targets, input tables, and identifier mappings;
- every threshold, search bound, stage switch, sample count, process count, and seed;
- a safe `overwrite: bool = False` switch;
- no loaded model, DataFrame, callback, GUI object, or other non-serializable state.

Use `Path` internally. `from_json(path)` must read a UTF-8 JSON object and resolve every relative
input and output path from the config file's directory, not the process working directory.
`from_mapping(...)` should construct nested CMM types such as `Medium`, `Condition`, and sampling
settings. `__post_init__` should normalize sequences and paths and then call `validate()`.

Validation must reject bad work before a solver runs. Check non-empty identifiers, numeric
ranges, mutually dependent stage switches, candidate capacity guards, supported method names,
seed types and ranges, and ambiguous configurations. A capacity field is a preflight guard; it
must not become an undocumented runtime top-*N* slice.

Defaults are part of the scientific API. Choose them deliberately, document them, and include
their resolved values in `00_config.json`. Never let an external library invent a hidden seed.

### 4.2 Frozen result and status rows

`MyWorkflowResult` should hold the config, typed outputs from every stage, workflow-level
preflight rows, combined provenance, the exported run path, and an immutable artifact index.
Use frozen row dataclasses for workflow annotations that do not belong to a numerical service.
Store the full service result objects rather than reducing them to selected scalar values.

Provide a derived `summary()` only for headline counts and cross-checks. Every numerical value
in the summary must be derivable from a raw CSV or provenance JSON in the same run. The summary
must not become a hidden second ranking implementation.

Keep status explicit at both levels:

- analysis rows retain solver states such as `optimal`, `infeasible`, or method-specific
  unavailable states;
- artifact records use `complete`, `partial`, `skipped`, or `failed` and require a reason for
  every non-complete state;
- a candidate execution index contains every workflow-defined candidate, including candidates
  whose target-specific file could not be produced.

### 4.3 Public runner

The public runner accepts only the config, loads the exact model, applies the selected solver,
and delegates to a loaded-model helper used by focused tests. Match SC-01's high-level shape:

```python
from pathlib import Path

from cobra.io import read_sbml_model


def run_my_workflow(config: MyWorkflowConfig) -> MyWorkflowResult:
    model_path = Path(config.model_path).expanduser().resolve()
    if not model_path.is_file():
        raise FileNotFoundError(f"model_path is not a file: {model_path}")
    if model_path.suffix.lower() not in {".xml", ".sbml"}:
        raise ValueError("the workflow accepts SBML .xml/.sbml models")

    model = read_sbml_model(str(model_path))
    if config.solver is not None:
        model.solver = config.solver
    return _run_my_workflow(model, config)
```

This snippet is also non-installed tutorial code. In the real loaded-model helper:

1. validate the config again;
2. fingerprint the source model and create isolated working copies;
3. resolve and apply every condition without contaminating another condition;
4. run model, exchange, identifier-overlap, baseline, reachability, and capability preflight;
5. stop on an invalid scientific run before expensive analysis;
6. call public numerical services in the declared order, forwarding the same resolved condition,
   reference state, thresholds, and seeds;
7. preserve method-specific results, ties, failures, and unavailable rows;
8. build validation indexes from the complete declared candidate universe;
9. generate workflow provenance, including the source fingerprint and one conditioned-model
   fingerprint per distinct condition;
10. export only after the in-memory scientific result is internally consistent.

Do not catch all exceptions and convert them into empty successful tables. Catch errors only at
a boundary where the scientific contract explicitly permits a `failed` or `skipped` artifact;
store the exception type and message in that artifact's reason or candidate index.

### 4.4 Output-directory safety

Resolve the run root before writing and reject a non-directory path. If the directory exists and
is non-empty, fail unless `overwrite=True` and its existing authoritative manifest has the same
workflow schema id. Derive owned paths from that validated manifest, remove only those exact
paths, and never recursively delete the requested root. Refuse overwrite when unowned entries
would remain, because the result would no longer be a self-contained exhaustive run. Reject
absolute artifact paths, `..` traversal, and any resolved path outside the run root. Treat stale
symlinks as hostile input and do not follow them during cleanup.

Read the exact source model bytes into memory before cleanup so an in-place reproduction whose
model path points into the previous bundle cannot erase its own input; write the archived copy
after the scoped cleanup succeeds. Also export one conditioned SBML for every distinct condition
used for solving. A fingerprint detects model changes; it does not replace the model file.

## 5. Own a workflow-specific artifact schema

Do not pass a new workflow to `cmm.reporting.validate_run`. Despite its generic name, the current
function validates the SC-01 schema-v2 roles in `ARTIFACT_CONTRACTS`, requires production report
fields such as `product_label`, and enforces production-specific numerical invariants. It checks
the schema version but does not establish workflow compatibility by validating `schema_id`.
Likewise, `build_publication_report`, `render_publication_figures`, the checked-in R script, and
`render_production_report` know SC-01 figure ids and source tables. They are not a generic
reporting engine.

Give the new workflow its own schema id, version, artifact contracts, validator, renderer, and
HTML assembly. Reuse neutral behavior only after it has been deliberately extracted into a
shared module with regression tests proving that SC-01 remains unchanged. Never append the new
workflow's required roles to the production `ARTIFACT_CONTRACTS` union: that would make both
validators accept an incoherent mixture of schemas.

### 5.1 Directory contract

Choose numbered stage directories that match the scientific contract. The following shape is a
non-installed template, not a current CMM run schema:

```text
<run>/
  00_config.json
  00_provenance.json
  00_summary.json
  00_manifest.json
  model/
    <source-model>.xml
    <model-id>__<condition-name>.xml
  01_preflight/
    preflight.csv
    preflight.metadata.json
  02_<stage>/
    <method>.csv
    <method>.metadata.json
  03_<stage>/
    candidate_index.csv
    candidate_index.metadata.json
  scripts/
    <workflow>_config.json
    reproduce.py
    render.py
    validate.py
  figures/
    figure_manifest.json
    <figure-id>.png
    <figure-id>.pdf
    <figure-id>.svg
  report.html
  report_standalone.html
  report_validation.json
```

Use one CSV per numerical result or semantically coherent tidy table. Export through result
objects' `to_frame()` methods; do not hand-transcribe values. Every complete CSV gets a UTF-8
JSON metadata sidecar containing its service provenance and workflow annotations. Keep units in
the artifact contract or in unambiguous column names and report all CMM flux units consistently.

Never interpolate an untrusted condition, method, target, or model id directly into a filename.
Derive a bounded filesystem-safe display slug plus a stable digest, reject collisions before
writing, and retain the complete original id in CSV and manifest metadata. Renderers and
validators still resolve artifacts by semantic role; a slug is only a safe path component, not
an artifact-discovery key.

Generated replay scripts must locate the run relative to their own file, use only public APIs,
and reproduce analysis, rendering, and validation separately. The copied replay config sets
`overwrite` to false by default so executing an archived script cannot silently replace evidence.

### 5.2 Manifest contract

`00_manifest.json` is the only artifact-discovery surface. Renderers and validators must resolve
semantic roles from it rather than globbing for filenames. Use a structure like this:

```json
{
  "schema_id": "cmm.my-workflow",
  "schema_version": 1,
  "authoritative": true,
  "status": "complete",
  "report": {
    "title": "CMM <scientific workflow title>",
    "language": "en"
  },
  "directories": [
    "model",
    "01_preflight",
    "02_<stage>",
    "03_<stage>",
    "scripts",
    "figures"
  ],
  "artifacts": {
    "provenance": {
      "path": "00_provenance.json",
      "stage": "workflow",
      "role": "provenance",
      "media_type": "application/json",
      "status": "complete",
      "sha256": "<lowercase sha256>",
      "size_bytes": 1234
    },
    "analysis_result": {
      "path": "02_<stage>/<method>.csv",
      "stage": "02_<stage>",
      "role": "analysis_result",
      "media_type": "text/csv",
      "status": "complete",
      "method": "<method>",
      "metadata_path": "02_<stage>/<method>.metadata.json",
      "sha256": "<lowercase sha256>",
      "size_bytes": 5678
    }
  },
  "supplementary_artifacts": [
    {
      "path": "02_<stage>/<method>.metadata.json",
      "stage": "02_<stage>",
      "role": "analysis_result_metadata",
      "media_type": "application/json",
      "status": "complete",
      "sha256": "<lowercase sha256>",
      "size_bytes": 2345
    }
  ]
}
```

The placeholder schema id above is not registered. The real manifest rules are:

- `schema_id` is a stable, workflow-specific string; `schema_version` is an integer.
- Breaking role, column, status, or invariant changes require a new schema version. Never change
  the meaning of an existing version in place.
- Primary `artifacts` are keyed by unique semantic role. Supplementary records still carry a
  role and are integrity-checked.
- Every record carries stage, role, media type, and status. A record that declares a file also
  carries its relative POSIX path, SHA-256, and byte count; method and metadata path are present
  when relevant.
- `complete` and `partial` records declare files that exist and are non-empty. `skipped` and
  `failed` records include a reason and may omit path and integrity fields when no file exists.
- Hashes and sizes are computed after final bytes are written. Rewriting an artifact requires
  rebuilding its manifest entry.
- `00_manifest.json` is the unhashed trust root and is not listed in its own artifact inventory;
  a file cannot contain its final digest without changing that digest. An archive may carry a
  separate external checksum file, but the manifest validator must never require a self-hash.
- The top-level status is derived from required stage roles; it is not an optimistic constant.

Keep the resolved config, provenance, summary, source model, every conditioned model, preflight,
raw method results, candidate coverage indexes, and replay scripts in the required-role
contract. Figures and reports are post-analysis outputs and should have their own figure/report
manifests and post-render checks rather than being mistaken for numerical source artifacts.

## 6. Add a dedicated validator

Place the workflow validator in a clearly named reporting module, for example
`cmm.reporting.my_workflow_schema`. Expose two levels:

```python
from __future__ import annotations

from pathlib import Path


def validate_my_workflow_source_run(
    run_dir: str | Path,
) -> ValidatedMyWorkflowRun:
    ...


def validate_my_workflow_run(run_dir: str | Path) -> ValidationReport:
    ...
```

The exact result type may reuse a genuinely neutral validation primitive if one has been
extracted, but the artifact contracts and scientific invariants remain workflow-specific. The
source validator raises an aggregated error containing every discovered issue. The public
completion validator returns a non-raising report with `valid`, `issues`, `warnings`, phase, and
the validated run so the CLI can produce a stable exit code and JSON payload.

### 6.1 Structural and integrity checks

Before rendering, validate all of these:

- run directory and authoritative manifest exist;
- manifest is UTF-8 JSON containing an object;
- exact `schema_id`, supported integer `schema_version`, and `authoritative: true` match;
- report metadata required by this workflow is present and English when the workflow promises
  an English report;
- every required semantic role is declared exactly once with the required status;
- no declared relative path is absolute, contains `..`, or resolves outside the run root;
- available artifacts exist as regular, non-empty files;
- declared lowercase SHA-256 and non-negative integer byte count match actual bytes;
- JSON artifacts and metadata sidecars parse to the required shape;
- CSV headers contain the contract columns in a compatible schema;
- complete CSV artifacts declare an existing metadata sidecar that is itself inventoried and
  integrity-checked;
- skipped or failed artifacts state a reason;
- duplicate paths, duplicate roles where uniqueness is required, and unknown primary roles are
  rejected or handled by an explicitly documented compatibility rule.

Validation is read-only. It must never repair files, rerun solvers, fill missing rows, regenerate
hashes, or downgrade a required artifact to a warning.

### 6.2 Scientific and cross-artifact checks

Structural validity is insufficient. Add checks specific to the workflow's claim, including:

- finite values and allowed status vocabulary per numeric column;
- identifier uniqueness, rank ordering, deterministic tie behavior, and expected row counts;
- bounds, fractions, signs, and units consistent with the declared method;
- model, condition, solver, parameter, and seed agreement between resolved config, provenance,
  summary, and metadata sidecars;
- source and every conditioned-model digest matching their declared provenance;
- complete coverage of the predeclared candidate universe in execution indexes;
- method-specific required stages present independently rather than inferred from a consensus;
- lethal, infeasible, partial, and contradictory outcomes retained rather than filtered;
- summary counts exactly recomputed from raw tables;
- every narrative claim and figure category traceable to an artifact status and source column.

Write one focused test for every invariant. Also tamper with a valid fixture one field at a time
and assert that the validator names the violation: a changed byte, escaped path, stale size,
missing column, missing sidecar, duplicate id, omitted candidate, inconsistent seed, or stale
summary must all fail.

### 6.3 Schema dispatch

Before adding dispatch, harden the SC-01 source validator to require the exact
`cmm.production-target-discovery` schema id as well as its supported schema version. Keep the
public `validate_production_run` alias unchanged, but make it reject a role-compatible foreign
manifest. Then update `cmm report render RUN_DIR` and
`cmm report validate RUN_DIR` to read only `schema_id` and `schema_version` from the manifest,
then dispatch through an explicit mapping to that workflow's renderer or validator. Unknown ids
and versions must fail with the supported values. Do not detect a workflow from directory names
or available CSV filenames, and do not send a non-production manifest through
`validate_production_run`.

Keep the existing Python aliases `render_production_report` and `validate_production_run`
unchanged for SC-01 compatibility. Add equally explicit aliases for the new workflow.

## 7. Build a workflow-specific R renderer and report

A canonical publication renderer consumes validated artifacts; it does not run COBRA models,
select candidates, or make scientific decisions. Add one checked-in R script whose required
roles, columns, and figure ids match only the new workflow schema.

Follow the production renderer's reproducibility boundary:

1. Invoke `Rscript --vanilla` with an argument list, never through an interpolated shell string.
2. Check the minimum R and package versions before reading data. Restore the repository's
   `renv.lock` for the exact publication environment.
3. Read only paths found by semantic role in the authoritative manifest.
4. Treat R warnings as render failures so dropped rows, invalid scales, and clipping-related
   warnings cannot pass unnoticed.
5. Use `ggplot2` for composition and save every rendered figure as a 300-DPI PNG plus editable
   PDF and SVG siblings.
6. Use explicit physical dimensions, font sizes, line widths, color-blind-safe colors, axis
   labels with units, deterministic factor order, and fixed facet order.
7. Give long labels enough margin or use controlled wrapping/repulsion. Inspect raster and
   vector output at final manuscript size for cropped labels, legends, and annotations.
8. Write `figures/figure_manifest.json` with renderer versions, script SHA-256, figure status,
   output paths, dimensions, raster DPI, and every source artifact used by each figure.
9. Mark an optional panel `skipped` or `failed` with a reason. A required panel may not disappear
   because its input table is empty.

`ggplot2` itself is pure R, but graphics dependencies such as `ragg`, `systemfonts`,
`textshaping`, and `svglite` may need platform-specific binaries or compilation support. Preserve
the existing Windows, macOS, and Linux R CI matrix, the locked package restore, and the renderer's
cross-platform font fallback. In Python, set deterministic environment values without assuming
that Windows supports POSIX locales: use `TZ=UTC` and `LANG=C.UTF-8`, set `LC_ALL=C.UTF-8` only
on non-Windows platforms, and pass the restored R library through `R_LIBS_USER` when needed.

Build linked and standalone English HTML only after the figure manifest validates. The report
builder should:

- derive all numbers, method names, statuses, and captions from validated CSV/JSON artifacts;
- include the scientific question, declared conditions, methods, result panels, figure legends,
  provenance, and references required by the scenario contract;
- keep method-specific results separate unless a combination rule was declared in advance;
- avoid synthesizing a target recommendation, causal conclusion, or wet-lab claim that no raw
  artifact supports;
- cite every figure's source tables and link supplementary artifacts;
- escape untrusted text and reject output paths outside the run root;
- produce deterministic linked HTML and a standalone copy with all figures embedded.

Post-render validation checks figure ids and statuses, non-empty PNG/PDF/SVG siblings, declared
dimensions and DPI, source-table existence, HTML figure coverage, standalone embedding, expected
method/reference coverage, and an up-to-date `report_validation.json`. Rendering success alone is
not completion.

## 8. Expose the workflow without broadening the numerical layer

After the Python workflow, schema, validator, and renderer pass their own tests:

1. Export the config, result, error, row types that users need, and runner from
   `cmm.workflows.__init__`.
2. Export the dedicated renderer, validator, and their public result types from
   `cmm.reporting.__init__`.
3. Add `cmm <workflow-command> --config CONFIG` to `cmm.cli`. It may parse arguments, invoke
   the public runner, optionally render, validate, print paths or JSON, and map exceptions to a
   nonzero exit code. It must not contain ranking, solver, or reporting logic.
4. Extend `cmm report render` and `cmm report validate` with manifest-schema dispatch as
   described above.
5. Keep the Qt application optional. If a GUI entry is added, it remains a thin view over the
   same public runner and is not required for headless execution.

The normal analysis command should mirror the production user experience:

```text
cmm <workflow-command> --config CONFIG
cmm <workflow-command> --config CONFIG --analysis-only
cmm report render RUN_DIR
cmm report validate RUN_DIR --json
```

The first command performs analysis, workflow-specific rendering, and final validation unless
`--analysis-only` is present. `output_dir` is required for a reportable CLI run. The Python
runner may allow `output_dir=None` for focused in-memory use if the config and result contract
document that behavior.

## 9. Update routing, documentation, skills, and feature status

Do not advertise a workflow before its imports, CLI, validator, and tests exist. In the same
change that ships it:

- add the goal and composition rules to `AGENTS.md` and `docs/scenarios/README.md`;
- add the public config fields, result types, Python example, and CLI commands to
  `docs/agent-reference.md`;
- document runtime ownership and data flow in `docs/architecture.md`;
- add numerical method and workflow/report integrity evidence to `docs/VALIDATION.md`;
- add a concise user entry point to `README.md` and link the full scenario/tutorial;
- update the feature roadmap without presenting a planned capability as available;
- move the exact feature name from `PLANNED_FEATURES` to `INCLUDED_FEATURES` only when the
  complete tested public surface ships. Do not add a generic workflow-engine feature merely
  because one more concrete workflow exists.

Add a repository skill only when goal-level prompts need reliable routing or clarification that
is not already covered by an existing skill. Put it under `.agents/skills/cmm-<workflow>/` with
an English `SKILL.md` that states:

- narrow trigger and explicit non-trigger examples;
- exact required inputs and which facts to inspect before asking the user;
- clarification decisions, each with an evidence-based recommended option;
- the public Python/CLI boundary, stage order, solver gates, and no-silent-fallback rule;
- run-directory and validation completion requirements;
- claim boundaries and stop conditions;
- links to the scenario and only the reference sections needed for execution.

The skill must call the canonical workflow; it must not recreate it with an ad hoc script or
introduce a second candidate-selection policy. Confirm the referenced imports and commands in
the checked-out version before its description says the workflow is shipped.

## 10. Test in layers

Use deterministic toy models for fast service and workflow tests, then add at least one realistic
integration fixture proportionate to the claim. Tests must not depend on the GUI.

### Numerical service tests

- known optimum or invariant on a small model;
- solver-capability failure before formulation/solve;
- condition propagation and model-state isolation;
- infeasible, lethal, empty, tied, and non-optimal cases;
- deterministic results and metadata for explicit seeds;
- stable `to_frame()` columns and provenance keys.

### Config and orchestration tests

- UTF-8 JSON parsing, nested settings, round-trip, and config-relative path resolution;
- invalid identifiers, ranges, seeds, capacity guards, and stage combinations;
- public file-backed runner and loaded-model test helper agree;
- preflight precedes expensive service calls;
- exact call order and forwarding of conditions, references, thresholds, and seeds;
- complete candidate coverage and declared deduplication with aliases retained;
- no input-model mutation leaks between stages or conditions;
- partial/skipped/failed results preserve status and reason;
- output-directory refusal, scoped overwrite, path containment, and replay config safety.

### Artifact and validator tests

- required roles, stable schema id/version, columns, sidecars, hashes, and byte counts;
- source and every distinct conditioned-model archive plus matching fingerprints;
- summary and provenance cross-checks;
- candidate-index coverage when individual analyses fail;
- one tampering test per structural and scientific invariant;
- old supported schema fixtures remain valid after a version bump;
- production schema fixtures remain valid and reject the new schema.

### Renderer and report tests

- missing `Rscript`, missing/incompatible package, R warning, and nonzero process exit fail;
- renderer reads only manifest-declared sources;
- required figures have non-empty 300-DPI PNG and editable PDF/SVG outputs;
- optional figures carry explicit unavailable status and reason;
- figure manifest records source data, package versions, dimensions, DPI, and script digest;
- HTML contains every required figure and source link; standalone HTML embeds all figures;
- report numbers and method coverage agree with the raw artifacts;
- visual QA at publication size catches clipped labels, truncated legends, overlapping
  annotations, excessive whitespace, and unreadable facets;
- R reporting tests run on Windows, macOS, and Linux with the checked-in `renv.lock`.

### CLI and regression tests

- help lists the new command;
- valid analysis-only, full analysis/render/validate, report render, and JSON validate paths;
- invalid config and validation failure return nonzero without a traceback-dependent interface;
- manifest dispatch selects the correct workflow and rejects unknown schema ids/versions;
- SC-01 Python, CLI, artifacts, renderer, validator, and skill behavior remain unchanged.

Run the repository's complete quality gate exactly as configured in `AGENTS.md`:

```bash
QT_QPA_PLATFORM=offscreen uv run --frozen --all-extras pytest -q --cov=cmm --cov-branch --cov-fail-under=80
uv run --frozen --all-extras ruff check src tests
uv run --frozen --all-extras ruff format --check src tests
uv run --frozen --all-extras mypy src/cmm/core src/cmm/features src/cmm/omics src/cmm/workflows src/cmm/reporting
```

Also run `git diff --check`, validate every new local Markdown link, and inspect the rendered
report and all figure formats rather than relying only on automated dimensions.

## 11. Definition of done

A new canonical workflow is complete only when every item below is true:

- [ ] Its distinct scientific question, supported claim, exclusions, inputs, conditions,
      decision rules, and stop conditions are documented.
- [ ] Every numerical method is a tested public solver-neutral service with a frozen result,
      stable `to_frame()`, and `run_provenance` metadata.
- [ ] A frozen config captures every input, parameter, stage switch, solver choice, and seed;
      JSON paths resolve from the config file.
- [ ] A frozen result retains all typed stage outputs, statuses, provenance, run path, and
      artifact index.
- [ ] The orchestrator performs preflight and capability gates before expensive work, forwards
      conditions consistently, and never silently substitutes or truncates a method.
- [ ] The exact source model and one conditioned model per distinct condition, resolved config,
      raw tables, sidecars, summary, provenance, manifest, and replay scripts are exported
      safely.
- [ ] A unique schema id/version and workflow-specific artifact contracts exist.
- [ ] A dedicated source and post-render validator checks integrity plus scientific invariants
      without mutating the run.
- [ ] A workflow-specific R renderer produces traceable 300-DPI PNG and editable PDF/SVG
      figures with no clipped text at manuscript size.
- [ ] Linked and standalone English reports contain only artifact-backed values and claims.
- [ ] Python exports, a thin analysis CLI, schema-dispatched report commands, and documentation
      are present and agree.
- [ ] Routing and an optional skill use the public canonical boundary and do not claim an
      unshipped feature.
- [ ] Unit, composition, artifact-tampering, rendering, CLI, cross-platform R, and SC-01
      regression tests pass.
- [ ] All four offscreen repository quality gates pass with branch coverage at or above 80%.
- [ ] `INCLUDED_FEATURES` and `docs/VALIDATION.md` are updated only after the complete contract
      above has shipped.

Passing this checklist establishes implementation correctness, reproducibility, and artifact
integrity. It does not establish biological validity. Model-derived interventions and state
predictions remain *in silico* hypotheses until independently tested.

## Reference implementation and supporting documents

- [SC-01 production-target discovery](../scenarios/SC-01-production-target-discovery.md)
- [SC-04 transformation-target discovery](../scenarios/SC-04-transformation-target-discovery.md)
- [Production workflow source](../../src/cmm/workflows/production.py)
- [Transformation workflow source](../../src/cmm/workflows/transformation.py)
- [Transformation report renderer](../../src/cmm/reporting/transformation.py)
- [Production schema validator](../../src/cmm/reporting/schema.py)
- [Production publication layer](../../src/cmm/reporting/publication.py)
- [Production R renderer](../../src/cmm/reporting/render_publication_figures.R)
- [Scenario reporting contract](../scenarios/_reporting.md)
- [Repository architecture](../architecture.md)
- [Public API reference](../agent-reference.md)
- [Scientific and implementation validation](../VALIDATION.md)
- [Clean-room policy](../clean-room-policy.md)
