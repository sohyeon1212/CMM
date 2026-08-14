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
  report.html
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
    medium="glucose_anaerobic",
    product="EX_succ_e",
    oxygen_bounds=(0.0, 0.0),          # the aeration, stated as bounds, not a boolean
    substrate="EX_glc__D_e",
    substrate_uptake=10.0,
)
(run_dir / "00_provenance.json").write_text(json.dumps(provenance, indent=2))
```

This records the model SHA-256, model id, active solver, Python/CMM/COBRApy/NumPy/pandas/SciPy
versions, the run timestamp, the seed, the solver version and the platform, and the parameters
you pass. Each individual result additionally carries its own `result.metadata` — keep those
too when a step's parameters differ from the run's.

**Record the applied condition in full**, as above: the medium preset name, the oxygen exchange
bounds, the substrate and its uptake rate, plus the medium components actually `applied` and any
`dropped`. A result file must state its own conditions without the reader reconstructing them
from a fingerprint. The `aerobic=True|False` shorthand was removed in 0.4.0 and must not appear
in new provenance records.

Also record, in the report if not in JSON: any method substitution forced by the solver, and
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
  `pd.DataFrame([vars(r) for r in rows]).to_csv(...)`. As of 0.4.0 that is **ten** columns, not
  seven: `objective_value` (the raw solver objective — `Σd²` for `moma_l2`, `Σ|d|` for
  `moma_l1`, a *switch count* for `room`), `distance` (a distance only — Segrè et al. Eq. (4)'s
  Euclidean `√(Σd²)`, `None` for ROOM), `distance_kind`, and `n_changed_reactions` replace the
  single overloaded `distance`. **Never label a column "distance" if it holds the objective**;
  runs before 0.4.0 did exactly that, and for `moma_l2` the two differ by a factor of about 36.
  `FvseofResult` likewise exports `park_type` and `capacity_slope` alongside the flux frames,
  and `FluxResponseResult.phases_frame()` exports the phase structure that replaced the
  removed `bottleneck` field.

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

Every figure in the report must be reproducible from a CSV in the same run directory.

---

## report.html

**The report is written as HTML, not Markdown**, so figures can sit next to the numbers they
explain instead of being listed at the end. A reader opens one file and sees the argument in
order; the PNGs stay separate at 300 DPI for reuse in a manuscript.

Reference figures with a **relative path** — `<img src="figures/fseof.png">` — so the run
directory stays self-contained and portable as a folder. Do not inline base64: it bloats the
file and hides the figures from anyone who wants the originals.

### Sections, in order

| Section | Content |
|---|---|
| **Summary** | Three to five sentences: the goal, what is recommended, and the single most important caveat. Someone reading only this must not be misled. |
| **Setup** | The preflight summary table from `_preflight.md`, plus medium, aeration, substrate, solver, and any method substitution with its reason. |
| **Results** | One subsection per pipeline step, in order. Each states what was run, the decisive numbers, the figure, and the CSV they came from. |
| **Recommended targets** | Target / type / evidence / predicted effect / confidence. Type is amplify, knockdown, or knockout. Evidence names the methods that agree. Confidence reflects how many independent methods agreed and whether verification passed. A scenario may use a richer shape — `SC-01` reports a strain proposal — but the columns above are the floor. |
| **Limitations** | What the analysis does *not* establish. At minimum: model predictions requiring experimental validation; the medium and aeration assumed; any conflicting method assumptions (MOMA's minimal adjustment vs OptKnock's growth maximization); solver capability that constrained the run. |
| **Provenance** | Model fingerprint, solver, CMM version, sampler seed, run directory, and the command or script that produced it. |

### Placing figures

Put each figure **inside the subsection that discusses it**, immediately after the sentence
that states its finding. A figure with no sentence pointing at it is decoration; a finding with
no figure beside it is harder to check. Give every one a caption naming the source CSV:

```html
<figure>
  <img src="figures/production_envelope.png" alt="Growth versus succinate flux">
  <figcaption>
    <b>Figure 1.</b> Growth falls from 0.21 h<sup>-1</sup> to zero as succinate is enforced.
    Source: <code>02_yield/production_envelope.csv</code>
  </figcaption>
</figure>
```

Keep the styling minimal and self-contained in a `<style>` block — no external CSS or fonts, so
the file renders the same anywhere. Tables should be plain `<table>` with hairline rules.

---

## Rules

1. **No number in the report that is not in a CSV.** If it is worth reporting it is worth
   exporting.
2. **Report infeasible and lethal outcomes.** An infeasible scan point or an essential-gene
   knockout is a result. Silently dropping them makes a target list look cleaner than it is.
3. **Name the assumptions in the report, not just in your reasoning**: reference flux state,
   aerobic/anaerobic, substrate, solver, sampler seed and count.
4. **Say what is a hypothesis.** FSEOF, rMTA, and sampling-based rankings prioritize
   candidates; they do not prove them. The Limitations section is not optional.
5. **Do not fabricate a run.** If a step failed or was skipped, say so and why, in place.
