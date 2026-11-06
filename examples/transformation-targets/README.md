# MTA/rMTA transformation targets

This example runs SC-02 on CMM's existing six-reaction, five-gene synthetic model. The direction
is explicitly **disease source → healthy target**: `SOURCE_EXPRESSION` to `TARGET_EXPRESSION`
from `cmm.app.screenshots`. These labels describe a constructed branch-rerouting fixture, not
patient data or a Recon model.

Run from the repository root with the frozen dependencies and R renderer installed:

```bash
uv run --frozen --all-extras python examples/transformation-targets/prepare.py
uv run --frozen --all-extras cmm transformation-targets --config examples/transformation-targets/config.json
uv run --frozen --all-extras cmm report validate results/example-transformation-mta --json
```

`prepare.py` only exports the existing model and expression constants; it creates no window
and performs no analysis. The public workflow owns preflight, the iMAT reference state,
direction tests, candidate construction, MTA, the MOMA baseline (enabled explicitly in
`config.json`, since MTA has no internal score decomposition), epsilon sensitivity, rendering
and validation. It requires MIQP-capable Gurobi; this small model fits its restricted license.
The condition caps A supply at 10 and declares all six reaction bounds. There is no oxygen
exchange in this synthetic network.

To execute the same defined experiment with published rMTA:

```bash
uv run --frozen --all-extras cmm transformation-targets --config examples/transformation-targets/rmta.json
uv run --frozen --all-extras cmm report validate results/example-transformation-rmta --json
```

[rmta.json](rmta.json) differs in method and output directory, and leaves the MOMA baseline
at its default (off): rMTA's score decomposition (`06_validation/score_decomposition.csv`) answers
the same question with no extra solve. rMTA uses roughly three solves per candidate and adds the
bTS/mTS/wTS component figure. Both methods retain three
eligible gene candidates in this fixture and rank `g2` first. Two other candidates tie, so
their ordering is alphabetical. Infeasible/failed results and tie warnings remain explicit.

There is one linear measurement per gene, so the config explicitly uses fold-change direction
tests instead of the published replicate t-test. `imat_delta` 0.5 classes g1, g2 and gb high and
g3 and g5 low in this five-gene fixture. Epsilon 0.01 and the sensitivity values 0.001/0.01 are
demonstration settings for this
model's flux scale. These choices are recorded in the run and must not be inherited as
scientific defaults for another dataset.

Open the output directory's `report_standalone.html` after validation succeeds. Each run keeps
the model, both original expression files, resolved configuration, provenance, CSVs, and
300-DPI PNG plus editable PDF/SVG figures. Its `scripts/reproduce.py` can replay the archived
inputs after relocation. Choose a new `output_dir` to repeat a run; overwrite is disabled.

## Using Recon or another model

A Recon1/Recon2/Recon3D model alone does not define a transformation study.

Copy this example into a study directory, then resolve these fields before invoking the same
`cmm transformation-targets --config CONFIG` command:

| Field | What must be supplied for the new study |
|---|---|
| `model_path` | The exact Recon release and SBML file, with its source recorded; do not interchange releases. |
| `source_expression_path`, `target_expression_path` | The original expression tables, on the scale declared by `expression_scale` (default linear: finite, non-negative), and the confirmed source → target direction. Gene ids must overlap this model's genes; map ids explicitly if needed. |
| `medium`, `condition` | The accepted full medium, substrate uptake, oxygen/aeration and all changed bounds, using this model's reaction ids. Remove the synthetic six-reaction condition. A medium closes exchanges only: decide open sink reactions (Recon3D ships 95 at -1000) with explicit bounds. Inspect growth before screening. |
| `direction` | Use `significance="ttest"` for replicated measurements; choose and record the ranking and changed-reaction cutoff. Single measurements need the disclosed fold-change alternative. |
| `method`, `perturbation` | Choose `mta` or `rmta` and gene or reaction interventions for the intended question; do not replace published methods with a continuous heuristic. |
| `imat_delta`, `reference_thinning`, `epsilon`, `validation.epsilon_sweep` | Choose δ for the rMTA discretisation (Valcárcel et al.: 0.3 Alzheimer's, 0.7 prostate cancer). Omit `epsilon` to use the default `"per_reaction"` (each changed reaction's threshold derived from the 2,000 samples behind v_ref). A numeric epsilon is in the actual model's flux units, and `TransformationWorkflowConfig.suggest_epsilon` can inform it. Sensitivity values repeat the full ranking. |
| `solver`, `output_dir` | A full Gurobi/CPLEX license for genome-scale MIQP and a fresh output directory. Relative paths resolve from the JSON file. |

Record data/model sources and mappings beside the config, then preserve the workflow's archived
inputs and provenance. Budget genome-scale MTA/rMTA and each epsilon pass from the candidate
count; a reduced candidate list changes the scope and must be declared. Predictions remain
in silico hypotheses, even when implementation and artifact checks pass.

For the agent execution contract use the
[transformation skill](../../.agents/skills/cmm-transformation-engineering/SKILL.md).
For adapting an existing workflow or adding a new one, use the
[customization guide](../../docs/building-custom-workflows.md) and
[contributor tutorial](../../docs/tutorials/adding-a-canonical-workflow.md).
