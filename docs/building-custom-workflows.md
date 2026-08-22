# Building a reproducible CMM workflow

CMM workflows are executable research protocols: they bind a scientific question to an exact
model, one explicit condition, solver requirements, typed numerical results, provenance, and a
validated artifact directory. This guide explains how to customize the shipped production
workflow and how to compose a genuinely new workflow from CMM's public services.

The desktop application is useful for exploration, but a workflow intended for a paper should
run through Python or the thin CLI so that its inputs and outputs can be replayed.

## 1. Choose the smallest correct boundary

There are three ways to use CMM. Choose by scientific scope, not by preferred interface.

| Need | Boundary | Reproducible output |
|---|---|---|
| Complete production-target study | `ProductionWorkflowConfig` and `run_production_target_discovery` | Canonical schema-v2 run, R report, validator |
| Same study with different thresholds, candidate counts, search seed, or sampling settings | Change the canonical config | Same canonical schema and validation contract |
| A different scientific question or a single analysis | Compose documented functions in `cmm.core`, `cmm.features`, or `cmm.omics` | Your own declared result and artifact contract |

Do not create a second SC-01 implementation just to change parameters. The existing production
workflow already exposes stage switches, method limits, the number of method-specific FSEOF and
FVSEOF targets, validation settings, a strain-design seed, and a sampling seed. Use direct
service calls when the question itself is different—for example, a standalone FVA study or a
condition comparison.

Start with the [scenario router](scenarios/README.md), then use the relevant signatures in the
[function reference](agent-reference.md). The [architecture](architecture.md) explains which
layer owns numerical work, workflow composition, rendering, and validation.

## 2. Write the scientific contract first

Before writing code or JSON, state the following in a short protocol:

1. **Question and decision rule.** What result changes the conclusion, and which thresholds are
   declared before looking at the ranking?
2. **Exact inputs.** Identify the SBML path, product exchange where applicable, expression
   files, and any id mapping. A model name is not a model version.
3. **One condition.** Record the medium, substrate exchange and uptake, oxygen exchange and
   bounds, objective, and every other changed bound.
4. **Method roles.** Separate inverse target discovery from forward phenotype prediction. Do
   not merge unlike scores into one ranking.
5. **Capabilities.** State the required LP, MILP, QP, or MIQP support and any optional package
   requirement. Define whether a missing capability stops or explicitly narrows the run.
6. **Artifacts.** Name each table, metadata sidecar, figure, and narrative claim that the run
   must produce.
7. **Completion gate.** Define the validator and tests that must pass. A generated HTML file is
   not itself proof that the workflow completed correctly.

The shared [preflight](scenarios/_preflight.md) and
[reporting contract](scenarios/_reporting.md) provide the detailed checks.

## 3. Customize the canonical production workflow

The recommended entry point for a production-engineering study is a UTF-8 JSON config. The
following generic example is intentionally not tied to one organism. Replace every reaction id
with an id present in the exact SBML model.

```json
{
  "model_path": "../model/organism.xml",
  "product": "EX_product_e",
  "output_dir": "../results/product_aerobic",
  "solver": "gurobi",
  "substrate": "EX_substrate_e",
  "biomass": "BIOMASS",
  "medium": {
    "mode": "explicit",
    "name": "defined_substrate_medium",
    "uptake": {
      "EX_substrate_e": 10.0,
      "EX_o2_e": 20.0
    },
    "required": ["EX_substrate_e", "EX_o2_e"]
  },
  "condition": {
    "name": "defined_aerobic",
    "bounds": [
      {
        "reaction_id": "EX_substrate_e",
        "lower_bound": -10.0,
        "upper_bound": 1000.0
      },
      {
        "reaction_id": "EX_o2_e",
        "lower_bound": -20.0,
        "upper_bound": 1000.0
      }
    ],
    "notes": "Defined substrate uptake and aerobic oxygen uptake"
  },
  "top_single_knockouts_per_method": 5,
  "top_amplification_targets_per_method": 10,
  "max_knockouts": 3,
  "strain_design_seed": 0,
  "fseof_steps": 10,
  "fvseof_steps": 10,
  "amplification_loop_diagnostic_top_n": 20,
  "validation": {
    "enabled": true,
    "max_flux_response_targets": 30,
    "flux_response_steps": 20,
    "flux_response_biomass_fraction": 0.3,
    "sampling_growth_fraction": 0.1,
    "sampling": {
      "enabled": true,
      "n": 1000,
      "method": "achr",
      "thinning": 100,
      "processes": 1,
      "seed": 0,
      "store_raw_samples": true
    }
  }
}
```

Positive values in the explicit medium are uptake allowances; the `ReactionBound` entries show
the corresponding COBRA exchange bounds directly. Keeping both in the resolved config makes the
condition auditable. Relative `model_path` and `output_dir` values resolve from the config
file's directory.

Keep `strain_design_seed` explicit even though its default is `0`. It is forwarded unchanged to
both OptKnock and RobustKnock and must be an integer in `0..2_000_000_000`; booleans and floats
are invalid. It is separate from `validation.sampling.seed`. Omitting explicit seed forwarding
would allow `straindesign` to invent a new random seed per call, changing MILP search paths,
runtime, or finite solution pools in otherwise identical runs.

Run analysis, rendering, and validation together:

```bash
uv run cmm production-targets --config workflows/product.json
uv run cmm report validate results/product_aerobic --json
```

Or separate numerical analysis from the R renderer:

```bash
uv run cmm production-targets --config workflows/product.json --analysis-only
uv run cmm report render results/product_aerobic
uv run cmm report validate results/product_aerobic --json
```

The equivalent Python boundary is useful when the workflow is called by another application:

```python
from cmm.reporting import render_production_report, validate_production_run
from cmm.workflows.production import (
    ProductionWorkflowConfig,
    ProductionWorkflowResult,
    run_production_target_discovery,
)

config = ProductionWorkflowConfig.from_json("workflows/product.json")
result: ProductionWorkflowResult = run_production_target_discovery(config)
if result.run_directory is None:
    raise RuntimeError("output_dir is required for a reportable run")

render_production_report(result.run_directory, renderer="nature-r")
validation = validate_production_run(result.run_directory)
validation.raise_for_errors()
```

FSEOF and FVSEOF remain independent methods. The workflow exports and plots the top configured
number from **each** ranking; overlap is descriptive, not an admission rule. Every
method-specific trajectory retains its diagnostic status, and every unique reaction in the two
report-visible lists proceeds to its own flux-response check. Loop-flagged or unresolved
candidates are still scanned but cannot be promoted to supported recommendations.
Consequently, the diagnostic cap must be at least twice the per-method amplification count, and
the response capacity must also accommodate both per-method D1–D5 single-knockout candidate
lists. With 10 amplification targets and five single knockouts per method, those minima are 20
and 30. These limits are capacity guards before deduplication, not top-*N* selectors. Config
validation rejects smaller enabled capacities, and execution never silently drops the tail of
a method.

The canonical knockout validation universe is the unique blocked-reaction signatures
represented by MOMA D1–D5 and ROOM D1–D5. Every candidate receives knockout-background flux
response and matched wild-type/knockout sampling; beneficial-selection and recommendation
filters are applied only after those analyses. Equivalent signatures are simulated once, with
all represented gene ids retained as candidate-id provenance. Every candidate also remains in
the response and sampling indexes if its analysis is infeasible, unavailable, skipped, or
failed.

## 4. Understand what the canonical run owns

The production workflow writes the exact source and conditioned models, resolved config,
provenance, raw CSVs, per-analysis metadata, and an authoritative manifest before reporting.
Its high-level flow is:

```text
config + exact SBML
        │
        ▼
preflight ─► yield/reference ─► MOMA + ROOM single deletions
        │                           │
        ├────────► OptKnock + RobustKnock
        │
        └────────► FSEOF + FVSEOF (independent rankings)
                              │
                              ▼
             all-candidate response / loop / KO sampling checks
                              │
                              ▼
                  manifest ─► R figures + HTML ─► validator
```

`00_manifest.json` is the authority for artifacts; the renderer does not discover files by
filename patterns. `render_production_report(..., renderer="nature-r")` reads only the declared
source tables and writes 300-DPI PNG plus editable PDF/SVG figures and linked/standalone HTML.
`validate_production_run` recomputes declared sizes and SHA-256 values, checks required methods
and columns, and verifies figure/report coverage without changing results.

The generated `scripts/production_config.json`, `reproduce.py`, `render.py`, and `validate.py`
are the run-local replay boundary. Preserve them with the exact model, `uv.lock`, `renv.lock`,
CMM version or commit, and the raw tables used by the manuscript.

## 5. Compose a new workflow from public services

When the question is outside SC-01, compose the documented public services and keep their typed
results intact. This small example asks whether one FSEOF hypothesis has a supportive
target-to-product response; it is not a replacement for the full production workflow.

```python
from dataclasses import dataclass
from typing import Mapping

from cobra.io import read_sbml_model

from cmm.core import (
    Condition,
    Medium,
    ReactionBound,
    apply_medium,
    fba,
    require,
    run_provenance,
)
from cmm.features import FluxResponseResult, flux_response, fseof, theoretical_yield


@dataclass(frozen=True)
class TargetResponseStudy:
    growth: float
    molar_yield: float
    target: str
    response: FluxResponseResult
    metadata: Mapping[str, object]


def run_target_response_study(
    model_path: str,
    *,
    product: str,
    target: str,
    substrate_exchange: str,
    substrate_uptake: float,
    oxygen_exchange: str,
    oxygen_uptake: float,
) -> TargetResponseStudy:
    model = read_sbml_model(model_path)
    medium = Medium(
        name="declared_defined_medium",
        uptake={
            substrate_exchange: substrate_uptake,
            oxygen_exchange: oxygen_uptake,
        },
        required=frozenset({substrate_exchange, oxygen_exchange}),
    )
    applied = apply_medium(model, medium)
    condition = Condition(
        name="declared_aerobic",
        bounds=(
            ReactionBound(
                reaction_id=substrate_exchange,
                lower_bound=-substrate_uptake,
                upper_bound=1000.0,
            ),
            ReactionBound(
                reaction_id=oxygen_exchange,
                lower_bound=-oxygen_uptake,
                upper_bound=1000.0,
            ),
        ),
        notes="Explicit substrate and oxygen uptake",
    )
    require("LP", model.solver.interface, feature="target-response workflow")

    baseline = fba(model, condition=condition)
    if baseline.status != "optimal" or baseline.objective_value <= 1e-6:
        raise RuntimeError("wild type does not grow in the declared condition")

    ceiling = theoretical_yield(model, product, condition=condition)
    if ceiling.molar_yield <= 1e-6:
        raise RuntimeError("product is unreachable in the declared condition")

    ranking = fseof(model, product, condition=condition, n_steps=10)
    if target not in ranking.amplification_targets():
        raise RuntimeError("declared target is not an actionable FSEOF hypothesis")

    response = flux_response(
        model,
        target,
        product,
        condition=condition,
        biomass_fraction=0.3,
        n_steps=20,
    )
    return TargetResponseStudy(
        growth=baseline.objective_value,
        molar_yield=ceiling.molar_yield,
        target=target,
        response=response,
        metadata=run_provenance(
            model,
            method="target_response_study",
            product=product,
            target=target,
            condition=condition.name,
            applied_medium=applied.to_provenance(),
            substrate=substrate_exchange,
            substrate_uptake=substrate_uptake,
            oxygen_exchange=oxygen_exchange,
            oxygen_bounds=(-oxygen_uptake, 1000.0),
        ),
    )
```

For a real custom workflow, define a frozen config dataclass with JSON parsing and validation,
return concrete result types, and export every numerical result through its `to_frame()`
method. Do not import private modules or copy numerical implementations into the orchestration
layer.

### Custom artifact contract

A custom run needs an explicit schema of its own. At minimum, include:

```text
run/
  00_config.json
  00_provenance.json
  00_summary.json
  00_manifest.json
  model/source.xml
  01_preflight/preflight.csv
  02_analysis/<method>.csv
  02_analysis/<method>.metadata.json
  figures/
  scripts/reproduce.py
  report.html
```

Each manifest record should carry a semantic role, relative path, status, reason for a skipped
or failed artifact, media type, byte size, SHA-256, and metadata-sidecar path. Preserve lethal,
infeasible, unavailable, skipped, failed, and contradictory results rather than filtering them
out during export. Candidate-index tables must cover the full workflow-defined candidate
universe even when a target-specific numerical artifact cannot be produced.

The shipped `nature-r` renderer and `validate_production_run` validate the canonical SC-01
schema; they are not generic validators for an arbitrary folder. If a new workflow is being
added to CMM itself, add a schema-specific renderer/validator or deliberately map it to an
existing documented schema, then expose that boundary from `cmm.workflows` and keep any CLI
command a thin adapter.

## 6. Test the workflow as a research instrument

A new workflow is not complete until tests cover both successful results and scientific
failure modes.

- **Config tests:** JSON round-trip, path resolution, invalid ranges, and explicit stage
  disablement.
- **Preflight tests:** zero growth, missing exchange, zero yield, and missing solver capability
  must stop or produce the documented unavailable state.
- **Composition tests:** assert service call order, shared condition, shared reference state,
  complete candidate-universe coverage without runtime slicing, capacity guards, and
  deterministic forwarding/provenance of both strain-design and sampling seeds.
- **Scientific invariants:** retain infeasible/status rows, keep FSEOF and FVSEOF rankings
  independent, run flux response for their complete report-visible top-10 union, run both flux
  response and sampling for every unique MOMA/ROOM D1–D5 candidate, rank RobustKnock by
  guaranteed product, and verify knockouts in the perturbed model.
- **Artifact tests:** required roles, metadata sidecars, hashes, path containment, replay scripts,
  and no hand-transcribed report values.
- **Rendering tests:** source-data traceability, non-empty PNG/PDF/SVG siblings, HTML figure
  coverage, and visual inspection for clipped labels.
- **Regression tests:** run the repository's full offscreen quality gate, not only the new test
  module.

```bash
QT_QPA_PLATFORM=offscreen uv run --frozen --all-extras pytest -q \
  --cov=cmm --cov-branch --cov-fail-under=80
uv run --frozen --all-extras ruff check src tests
uv run --frozen --all-extras ruff format --check src tests
uv run --frozen --all-extras mypy \
  src/cmm/core src/cmm/features src/cmm/omics src/cmm/workflows src/cmm/reporting
```

## 7. Prepare the publication handoff

Archive the following together:

- exact release tag or commit, `uv.lock`, `renv.lock`, and solver/version information;
- exact input and conditioned models, including their redistribution terms;
- resolved config, provenance, manifest, raw tables, validation result, and replay scripts;
- 300-DPI raster and editable vector figures with their source-data tables;
- CMM's citation plus the original paper for every method used;
- declared medium, aeration, substrate uptake, thresholds, search limits, sample count,
  thinning, strain-design seed, and sampling seed;
- limitations covering alternative optima, solver tolerances, GPR resolution, sampler
  convergence, model scope, and experimental validation.

Passing CMM's tests supports implementation and artifact integrity. It does not establish the
biological validity of a predicted intervention. Present targets as *in silico* hypotheses
until they have been tested experimentally.

Further reading:

- [SC-01 production target discovery](scenarios/SC-01-production-target-discovery.md)
- [Reporting and artifact contract](scenarios/_reporting.md)
- [Scientific validation and method references](VALIDATION.md)
- [Public function and result-object reference](agent-reference.md)
- [Repository architecture](architecture.md)
