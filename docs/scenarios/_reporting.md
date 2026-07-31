---
id: _reporting
title: Reporting — the artifact contract every scenario ends with
---

# Reporting contract

A scenario run is finished when someone else could reproduce it from the directory it left
behind. That means raw numbers, figures, provenance, and a narrative that references them —
not a summary in chat.

---

## Directory layout

```
results/<SC-id>_<model-id>_<YYYYMMDD-HHMMSS>/
  00_provenance.json
  01_preflight/
      wild_type_fluxes.csv
  02_yield/
      theoretical_yield.csv
      production_envelope.csv
  03_baseline/
      reference_<method>.csv          one per reference method used
  04_knockout/
      batch_<reference>_<method>.csv
  05_amplification/
      fseof_trends.csv
      fvseof_mean.csv  fvseof_forced.csv  fvseof_capacity.csv
  06_validation/
      flux_response_<target>.csv
      sampling_statistics.csv
      consensus.csv
  figures/
      *.png
  report.md
```

Step directories are numbered by the scenario's own steps; skip numbers a scenario does not
have rather than renumbering. One CSV per analysis, named for what it holds.

---

## Provenance

`00_provenance.json` is written once at preflight and holds what makes the run repeatable:

```python
import json
from pathlib import Path
from cmm.core import run_provenance

run_dir = Path("results/SC-01_e_coli_core_20260729-140355")
run_dir.mkdir(parents=True, exist_ok=True)

provenance = run_provenance(
    model,
    scenario="SC-01",
    medium="glucose_aerobic",
    product="EX_succ_e",
    aerobic=True,
)
(run_dir / "00_provenance.json").write_text(json.dumps(provenance, indent=2))
```

This records the model SHA-256, model id, active solver, Python/CMM/COBRApy/NumPy/pandas/SciPy
versions, and the parameters you pass. Each individual result additionally carries its own
`result.metadata` — keep those too when a step's parameters differ from the run's.

Also record, in `report.md` if not in JSON: any method substitution forced by the solver, and
the exact CMM version or git commit.

---

## Raw data

Every table comes from the result object, never retyped:

```python
result.to_frame().to_csv(run_dir / "05_amplification" / "fseof_trends.csv")
```

- `ProductionEnvelope.to_frame()`, `FluxResponseResult.to_frame()`,
  `SamplingResult.to_frame()` / `.statistics()`, `TargetRanking.to_frame()` all export
  directly.
- `FseofResult.trends` and the three `FvseofResult` frames are DataFrames already.
- Flux dicts (`FluxSolution.fluxes`, `FluxState.fluxes`) → `pd.Series(fluxes).to_csv(...)`.
- `batch_comparison` returns a list of dataclass rows →
  `pd.DataFrame([vars(r) for r in rows]).to_csv(...)`.

Units, stated once in the report and never converted: flux **mmol gDW⁻¹ h⁻¹**, growth
**h⁻¹**, molar yield **mol/mol**.

---

## Figures

```python
from cmm.visualization import save_figure, fseof_figure

save_figure(fseof_figure(result, top_n=6), run_dir / "figures" / "fseof.png")
```

`save_figure` writes 300 DPI, tight-cropped, white-background files and creates parent
directories. For a manuscript also emit PDF or SVG by changing the suffix. Use
`column_width=1` for single-column figures, `2` (default) for double.

Every figure in `report.md` must be reproducible from a CSV in the same run directory.

---

## report.md

```markdown
# <Scenario title> — <model id>

## Summary
Three to five sentences: the goal, the recommended targets, and the single most important
caveat. Someone reading only this should not be misled.

## Setup
The preflight summary table from `_preflight.md`, plus medium, aeration, substrate, solver,
and any method substitution and why.

## Results
One subsection per pipeline step, in order. Each states what was run, the decisive numbers,
the figure, and the CSV they came from.

## Recommended targets
| Rank | Target | Type | Evidence | Predicted effect | Confidence |
Type is amplify / knockdown / knockout. Evidence names the methods that agree. Predicted
effect gives product flux and growth. Confidence reflects how many independent methods
agreed and whether verification passed.

## Limitations
What the analysis does not establish. At minimum: these are model predictions requiring
experimental validation; the medium and aeration assumed; any method whose assumptions
conflict (MOMA's minimal adjustment vs OptKnock's growth maximization); solver capability
that constrained the analysis.

## Provenance
Model fingerprint, solver, CMM version, run directory, and the command or script that
produced it.
```

---

## Rules

1. **No number in `report.md` that is not in a CSV.** If it is worth reporting it is worth
   exporting.
2. **Report infeasible and lethal outcomes.** An infeasible scan point or an essential-gene
   knockout is a result. Silently dropping them makes a target list look cleaner than it is.
3. **Name the assumptions in the report, not just in your reasoning**: reference flux state,
   aerobic/anaerobic, substrate, solver, sampler seed and count.
4. **Say what is a hypothesis.** FSEOF, rMTA, and sampling-based rankings prioritize
   candidates; they do not prove them. The Limitations section is not optional.
5. **Do not fabricate a run.** If a step failed or was skipped, say so and why, in place.
