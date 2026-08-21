---
id: _reporting
title: Reporting — the artifact contract every scenario ends with
---

# Reporting contract

A scenario run is finished when someone else could reproduce it from the directory it left
behind. That means raw numbers, figures, provenance, and a narrative that references them —
not a summary in chat.

---

## Where the run directory goes

`CMM_RESULTS_DIR` if it is set, otherwise `results/` beside the working directory:

```python
import os
from pathlib import Path

root = Path(os.environ.get("CMM_RESULTS_DIR") or "results").expanduser()
run_dir = root / f"SC-01_{model.id}_{stamp}"
run_dir.mkdir(parents=True, exist_ok=True)
```

`results/` is in CMM's `.gitignore`, so a run started from a clone does not fill
`git status` with the user's own analyses or risk being committed with them. Set
`CMM_RESULTS_DIR` to keep runs somewhere of your choosing — outside the clone, on a larger
disk, or grouped per project:

```bash
CMM_RESULTS_DIR=~/analyses/succinate  python your_scenario_script.py
```

**Say in the report where the directory is** (§7), as an absolute path. A reader given only
the folder cannot otherwise tell where it came from.

## Directory layout

```
<run directory>/
  00_provenance.json
  00_summary.json
  model/
      <model-id>.xml                  the exact file the run used, not just its name
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
  scripts/
      *.py                            the scripts that produced this run, copied verbatim
  report.html
  report_standalone.html
```

Step directories are numbered by the scenario's own steps; skip numbers a scenario does not
have rather than renumbering. One CSV per analysis, named for what it holds.

**The layout above is exhaustive.** Nothing else belongs at the top of the run directory: every
file is either `00_provenance.json`, the pinned model, or inside a numbered step directory,
`figures/` or `scripts/`. A run that leaves scripts and intermediate JSON loose at the root is
still correct science and still unreadable as a deliverable — two runs of the same scenario
should be diffable folder against folder, and they are not if each one invents its own
arrangement.

Where things go when it is not obvious:

- **The scripts that produced the run** → `scripts/`. §7 requires the report to name them, so
  they have to be somewhere findable, and a path into a session's temporary directory is not a
  provenance record.
- **A summary of one step** → that step's directory, not the root.
- **A summary of the whole run** → `00_summary.json`, described below. This is the one summary
  that belongs at the top, because it is about the run rather than about any step in it.
- **A result that spans steps** — a coupling scan across growth floors, a consensus table —
  → the directory of the step that consumed it, which is usually the verification step. If two
  steps consume it, put it with the one that produced it.
- **Anything you would not hand to a reader** — scratch files, caches, half-written attempts —
  → outside the run directory entirely. The run directory is the deliverable, not the workspace.

---

## Provenance and summary

Two files at the root, with two jobs. `00_provenance.json` answers *how would someone repeat
this*; `00_summary.json` answers *what happened*. Keeping them apart is what lets the first be
compared between runs mechanically — a provenance record that also carried results would differ
on every run for reasons that have nothing to do with reproducibility.

### `00_summary.json`

Written at the end, holding the run's own account of itself: the headline number from each
step, the model fingerprint **before and after** the medium was applied (the pair is the
evidence that the condition was actually applied, not merely declared), any warning the medium
raised, and the result of every cross-check the scenario specifies. It is what a later reader
greps when they want a number without opening six CSVs, and what a second run is diffed against
when its conclusion disagrees with the first.

It is a convenience, not a source. Every number in it must also be in a CSV — see Rules — and
the report cites the CSV, never this file.

### `00_provenance.json`

Written once at preflight and holds what makes the run repeatable:

```python
import json
from pathlib import Path
from cmm.core import run_provenance

run_dir = ...  # as chosen above
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
- `batch_comparison` returns a `BatchComparisonResult` — a list of dataclass rows that also
  carries the screen's provenance → `screen.to_frame().to_csv(...)`, and write
  `screen.metadata` beside it as JSON (it is the only copy: the rows carry numbers, the
  container carries the run). `pd.DataFrame([vars(r) for r in screen])` still works and gives
  the same frame. As of 0.4.0 that is **ten** columns, not
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

**The report is HTML, not Markdown**, so a figure sits next to the number it explains instead
of being listed at the end. A reader opens one file and follows the argument in order; the PNGs
stay separate at 300 DPI for reuse in a manuscript.

### Two copies, and why

Ship both. They differ in one thing — how they reference figures — and using the wrong form
fails **silently**: the file saves without error and every figure is simply blank.

| File | `<img src>` form | For |
|---|---|---|
| `report.html` | relative path, `figures/fseof.png` | reading the run in place. Small, and the 300 DPI originals stay reachable beside it |
| `report_standalone.html` | `data:image/png;base64,…` | sending to someone. Survives leaving the folder |

`report.html` alone is not enough: a reader who downloads only that file gets a page with
every figure missing, and nothing says so. `report_standalone.html` alone is not enough
either: it hides the originals and runs to several megabytes.

In the standalone copy, links to data files cannot resolve, so render them as plain
monospace filenames rather than dead `<a href>`.

### Before you call it done

Check each of these on the standalone copy. Any failure means the report is broken for the
person you send it to, and none of them raises an error on its own.

- [ ] No `figures/` or other relative `src`/`href` survives. Search the file for `src="fig`.
- [ ] Number of `<img>` equals the number of figure captions.
- [ ] Figure numbers run 1, 2, 3 … in document order, with none repeated or skipped.
- [ ] Section numbers are unique, and every `§N` cross-reference points at a section that exists.
- [ ] Opening the file from a different directory still shows every figure.

### Sections

Number every section and cross-reference as `§4`. Sections 1–4, 6 and 7 are required; §5 is
the scenario's own answer and its heading changes with the scenario.

| # | Section | Content |
|---|---|---|
| 1 | **Summary** | The findings themselves, not a description of the work. Three to five sentences: what was found, what is recommended, and the single most important caveat. A reader who stops here must not be misled. |
| 2 | **Setup** | The preflight summary table from `_preflight.md`, plus medium, aeration as bounds, substrate and uptake, solver, and any method substitution with its reason. |
| 3 | **Data and methods** | Model and its fingerprint, sizes, every parameter that changes an answer — reference flux state, `n_steps`, sampler seed and count, thresholds. Enough that someone repeats the run without reading the scenario. |
| 4 | **Results** | One subsection per pipeline step, in order. Each states what was run, the decisive numbers, the figure, and the CSV they came from. |
| 5 | *scenario-specific* | SC-01 → **Recommended targets and strain proposal**; SC-02 → **Interpretation — what explains the difference**; SC-03 → **Essentiality classes**. See below. |
| 6 | **Limitations** | What the analysis does *not* establish. At minimum: predictions needing experimental validation; the medium and aeration assumed; conflicting method assumptions (MOMA's minimal adjustment vs OptKnock's growth maximization); solver capability that constrained the run. |
| 7 | **Provenance** | Model fingerprint, solver and version, CMM version or git commit, sampler seed, run directory, and the command or script that produced it. |

**§5 for SC-01** — a table whose floor is: target / type / **wild-type flux → optimum flux
(ratio)** / evidence / predicted effect / confidence. Type is amplify, knockdown or knockout.
Evidence names the methods that agree. Confidence reflects how many independent methods agreed
and whether verification passed. SC-01 additionally reports a strain proposal.

The flux column is not decoration and is not optional. A target whose optimum sits within a few
per cent of its wild-type flux is **not an intervention** — the cell is already where you would
put it — however highly the amplification scan ranked it, and a target whose optimum is *below*
its wild-type flux is a knockdown wearing an amplification label. Both are common: on anaerobic
`e_coli_core`/succinate, most of what FSEOF calls `amplify` is one or the other. Writing the
ratio into the table is what makes that visible, to the reader and to whoever wrote the row.
Reactions that are off in the wild type and switch on, or that move several-fold, are the ones
worth proposing. `SC-01` step 5b decision rule 3 is the source of this test; the table is where
its answer has to appear.

**§5 for SC-02** — what mechanism the per-condition fluxes support, and how strongly. Separate
what the transcriptome decided from what the declared condition decided; when the condition
dominates, say so, because a reader will otherwise credit the data.

**§5 for SC-03** — the essentiality classes with their counts, and the beneficial deletions.
An essential gene is a result, not a failed row.

Descriptive or purely operational runs may leave §5 empty rather than inventing an
interpretation. Say it is empty and why.

### How it is written

- **Headings are plain and descriptive** — `Data and methods`, `Limitations` — never narrative
  (`What the screen told us`, `Reframing the question`).
- **No process narration.** The report states what is true and how it was established, not the
  order the work happened in. Cut "we then", "in the first pass", "next we examined". §4 is
  organised by finding, not by chronology.
- **No emoji.**
- **Every number carries its basis**: the method that produced it, its unit, and the CSV it
  came from. Values are reported at the precision the CSV holds — do not round in the report
  and do not add digits.
- **Distinguish claim strength.** Keep apart what the numbers show, what they are consistent
  with, and what is conjecture. A hedge in one sentence does not cover an unhedged claim in
  the next.
- **Tables carry units in the header**, not in every cell, and headers are written for a
  reader: `Molar yield (mol mol⁻¹)`, not the exported column name `molar_yield_mol_per_mol`.
  The values themselves are copied from the CSV unchanged.
- **All text inside figures is English** — axis labels, titles, legends — and every figure
  carries a caption stating what it shows and what to take from it, naming its source CSV.
- Keep styling minimal and self-contained in a `<style>` block: no external CSS or fonts, so
  the file renders the same anywhere. Tables are plain `<table>` with hairline rules.

### Placing figures

Put each figure **inside the subsection that discusses it**, immediately after the sentence
that states its finding. A figure with no sentence pointing at it is decoration; a finding with
no figure beside it is harder to check.

```html
<figure>
  <img src="figures/production_envelope.png" alt="Growth versus succinate flux">
  <figcaption>
    <b>Figure 1.</b> Growth falls from 0.21 h<sup>-1</sup> to zero as succinate is enforced,
    so the two cannot be maximised together.
    Source: <code>02_yield/production_envelope.csv</code>
  </figcaption>
</figure>
```

---

## Rules

1. **No number in the report that is not in a CSV.** If it is worth reporting it is worth
   exporting.
2. **Report infeasible and lethal outcomes.** An infeasible scan point or an essential-gene
   knockout is a result. Silently dropping them makes a target list look cleaner than it is.
3. **Name the assumptions in the report, not just in your reasoning**: reference flux state,
   the aeration as bounds, substrate, solver, sampler seed and count. These belong in §3, where
   a reader looks for them, not scattered through §4.
4. **Say what is a hypothesis.** FSEOF, rMTA and sampling-based rankings prioritise candidates;
   they do not prove them. §6 is not optional.
5. **Do not fabricate a run.** If a step failed or was skipped, say so and why, in place. A
   report that omits a failed step reads as though the step succeeded.
6. **Ship both copies.** A run that leaves only `report.html` cannot be sent to anyone, and the
   failure is invisible until they open it.
