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

Choose `ProductionWorkflowConfig.output_dir` explicitly. The workflow does not guess where a
publication run belongs, and when `output_dir=None` it returns typed in-memory results without
claiming that a deliverable was written:

```python
config = ProductionWorkflowConfig(
    model_path="model/e_coli_core.xml",
    product="EX_succ_e",
    output_dir="results/SC-01_e_coli_core_succinate_20260821T071140Z",
    # medium= and condition= are also required run-definition decisions
)
result = run_production_target_discovery(config)
run_dir = result.run_directory
```

`results/` is in CMM's `.gitignore`, but an absolute location outside the clone is preferable
for manuscript runs. The CLI config records the resolved output path in `00_config.json`.

**Say in the report where the directory is** (§7), as an absolute path. A reader given only
the folder cannot otherwise tell where it came from.

## Directory layout

```
<run directory>/
  00_config.json
  00_provenance.json
  00_summary.json
  00_manifest.json
  model/
      <model-id>.xml                  byte-for-byte source SBML
      <model-id>__conditioned.xml     bounds/objective used for all analyses
  01_preflight/
      preflight.csv
      preflight.metadata.json
  02_yield/
      theoretical_yield.csv
      production_envelope.csv
      theoretical_yield.metadata.json  production_envelope.metadata.json
  03_reference/
      wild_type_fluxes.csv
      reference_fluxes.csv
      wild_type_fluxes.metadata.json  reference_fluxes.metadata.json
  04_single_knockout/
      gene_knockout_mapping.csv
      single_knockout_moma.csv
      single_knockout_room.csv
      single_knockout_candidates.csv
      single_knockout_consensus.csv
      <each table stem>.metadata.json
  05_strain_design/
      optknock.csv
      robustknock.csv
      optknock.metadata.json  robustknock.metadata.json
  06_amplification/
      fseof.csv  fseof_tidy.csv
      fvseof.csv  fvseof_tidy.csv
      <each table stem>.metadata.json
  07_validation/
      amplification_loop_diagnostic.csv
      amplification_loop_diagnostic.metadata.json
      flux_response_index.csv
      flux_response_tidy.csv
      flux_response_index.metadata.json  flux_response_tidy.metadata.json
      flux_response__<target>.csv
      flux_response_phases__<target>.csv
      flux_response__<target>.metadata.json
      random_sampling_index.csv
      sampling_tidy.csv
      random_sampling_index.metadata.json  sampling_tidy.metadata.json
      random_sampling__<target>.csv.gz
      random_sampling_statistics__<target>.csv
      random_sampling_comparison__<target>.csv
      random_sampling__<target>.metadata.json
      recommendations.csv
      recommendations.metadata.json
  figures/
      *.png                           300-DPI manuscript raster
      *.pdf and *.svg                 editable vector counterparts
      figure_manifest.json
  scripts/
      production_config.json
      reproduce.py
      render.py
      validate.py
  report.html
  report_standalone.html
```

This is the canonical SC-01 schema. Other scenarios keep their own numbered step directories,
but the root control files, model, figure, script/config, and report rules are shared. One CSV
per analysis is named for the quantity it holds; target-specific validation tables are listed
by the two index files rather than discovered through fragile filename globs.

The index files are coverage ledgers, not success-only summaries. `flux_response_index.csv`
contains every unique MOMA/ROOM D1–D5 knockout candidate and every unique reaction in the
independent FSEOF/FVSEOF top-10 lists. `random_sampling_index.csv` contains every unique
MOMA/ROOM D1–D5 knockout candidate plus the shared wild-type reference. Each row records
execution status and reason;
target-specific result files exist for successful attempts, and equivalent knockout genes remain
listed in the candidate-id provenance for their shared blocked-reaction signature. Skipped,
infeasible, unavailable, and failed attempts remain auditable through their index rows.
Loop-flagged or unresolved amplification candidates still receive a response attempt, with
loop/support eligibility kept as separate evidence.
For knockout candidates, one blocked reaction receives a pre-deletion wild-type scan:
reference↔zero when reference flux is nonzero, otherwise the full feasible reaction domain
marked exploratory and not causal support for deletion. A multi-reaction blocked signature
remains explicit unavailable/skipped because no single flux-response x-axis represents it; no
blocked reaction may be selected silently.

For machine audit, knockout index rows use
`candidate_scope=all_display_ranked_candidates` and identify the simulated phenotype with
`blocked_reaction_signature` while preserving aliases in `candidate_target_ids`. Amplification
rows use `candidate_scope=all_report_selected_candidates`. Response rows keep
`loop_diagnostic_status`, `loop_artifact_flag`, `loop_diagnostic_eligible`, and
`loop_diagnostic_reason` separate from execution `status`. `00_summary.json` and
`00_provenance.json` record `validation_candidate_policy` and the expected, attempted,
completed, failed, and skipped coverage counts used by the completion gate.

**The layout above is exhaustive.** Nothing else belongs at the top of the run directory: every
file is one of the four `00_*.json` controls, a report, the pinned model, or inside a numbered
step directory, `figures/` or `scripts/`. The `<stem>.metadata.json` notation above means one
sidecar for each listed manifest-linked primary table, not an open-ended filename convention.
A run that leaves scripts and intermediate JSON loose at the root is still correct science
and still unreadable as a deliverable — two runs of the same scenario should be diffable
folder against folder, and they are not if each one invents its own arrangement.

Where things go when it is not obvious:

- **The scripts that reproduce, render, and validate the run** → `scripts/`. The canonical
  workflow writes the resolved `production_config.json` plus `reproduce.py`, `render.py`, and
  `validate.py`; §7 names them. A path into a session's temporary directory is not a provenance
  record.
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

### `00_config.json`

The fully resolved workflow input after defaults have been applied: model path, product
exchange, solver, condition bounds, method parameters, the explicit strain-design and sampling
seeds, search limits, and output location. This is the human- and machine-readable answer to
"what was requested?" It is not a replacement for provenance, because it does not establish
what code or solver executed it.

### `00_summary.json`

Written at the end, holding the run's own account of itself: model/product/substrate/biomass,
wild-type and theoretical headline values, and result counts for each enabled stage. Source
and conditioned fingerprints, the applied medium, warnings, solver, versions, and detailed
method parameters live in `00_provenance.json`. `00_summary.json` is what a later reader greps
when they want a headline without opening the step CSVs.

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
    strain_design_seed=0,
    sampling_seed=0,
)
(run_dir / "00_provenance.json").write_text(json.dumps(provenance, indent=2))
```

This records the model SHA-256, model id, active solver, Python/CMM/COBRApy/NumPy/pandas/SciPy
versions, the run timestamp, method seed, solver version and platform, and the parameters you
pass. Each individual result additionally carries its own `result.metadata` — keep those too
when a step's parameters differ from the run's. In SC-01, `strain_design_seed` is explicit in
the resolved config and is forwarded as the recorded `seed` for both OptKnock and RobustKnock;
the sampling seed is recorded separately with sampling artifacts. A backend-generated hidden
seed is an incomplete provenance record.

**Record the applied condition in full**, as above: the medium preset name, the oxygen exchange
bounds, the substrate and its uptake rate, plus the medium components actually `applied` and any
`dropped`. A result file must state its own conditions without the reader reconstructing them
from a fingerprint. The `aerobic=True|False` shorthand was removed in 0.4.0 and must not appear
in new provenance records.

Also record, in the report if not in JSON: any method substitution forced by the solver, and
the exact CMM version or git commit.

### `00_manifest.json`

The authoritative schema-v2 artifact map is keyed by semantic role. Primary entries record
relative path, status, reason when applicable, SHA-256, byte size, and the linked
`metadata_path`; supplementary records also carry stage, method, and media type. Each analysis
sidecar records the artifact role/status/reason and the method-specific provenance or evidence
policy. Report rendering and validation read this map instead of guessing filenames.
The validator rejects missing required roles, path escapes, unreadable/empty data, and missing
required CSV columns. It recomputes each declared SHA-256 and byte size, including metadata
sidecars and supplementary artifacts, so an edited or truncated run fails before rendering.

---

## Raw data

Every table comes from a result-object export, never from a hand-built report table:

```python
result.to_frame().to_csv(run_dir / "06_amplification" / "fseof.csv", index=False)
```

- `ProductionEnvelope.to_frame()`, `FluxResponseResult.to_frame()`,
  `SamplingResult.to_frame()` / `.statistics()`, `TargetRanking.to_frame()` all export
  directly.
- `ProductionYield.to_frame()` and `.carbon_uptake_frame()` separate the headline yield from
  its carbon-balance inputs.
- `FseofResult.to_frame()` / `.trajectory_frame()` and
  `FvseofResult.to_frame()` / `.trajectory_frame()` separate rankings from tidy scan data.
- `StrainDesignResult.to_frame()` includes rank, knockout set, maximum product, guaranteed
  product, growth, and coupling verdict.
- `ComparisonResult.summary_frame()` and `.fluxes_frame()` keep one solve summary distinct
  from its reaction-level flux state.
- `GeneKnockoutMappingResult.to_frame()` preserves gene/reaction names, inert status, blocked
  reactions, and the complete GPR used to turn a computational knockout into an experimental
  intervention.
- `single_knockout_candidates.csv` preserves `validation_target_id`,
  `candidate_source_methods`, and `validation_representative`, so equivalent gene ids remain
  traceable to the one blocked-reaction phenotype that was simulated.
- `AmplificationLoopDiagnosticResult.to_frame()` preserves standard and loopless FVA bounds,
  capacities, their ratio and threshold, the artifact flag, status/reason, and the enforced
  product/biomass floors. It supports the recommendation policy without rewriting the raw
  FSEOF/FVSEOF ranking tables.
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

The canonical SC-01 report is rendered from saved CSVs, not in-memory Python objects:

```python
from cmm.reporting import render_production_report, validate_production_run

render_production_report(run_dir, renderer="nature-r")
validation = validate_production_run(run_dir)
```

`nature-r` is CMM's restrained, Nature Genetics-inspired house style; it is not a journal
endorsement or a substitute for checking the current submission instructions. It runs an
auditable `Rscript` renderer with declared package versions and enforces:

- final widths of 89 mm (single column) or at most 180 mm (double column);
- 5–7 pt sans-serif text at final size, sentence-case labels, and panel labels `a`, `b`, …;
- colour-blind-safe colours plus shape/line-type redundancy; no rainbow scale or decorative
  background;
- comparable panels with shared limits, visible wild-type references, physical units, and no
  plot title duplicated from the caption;
- 300-DPI PNG and editable PDF/SVG counterparts produced from the same source-data CSV.

These values follow the current
[Nature Genetics formatting guidance](https://www.nature.com/ng/submission-guidelines/aip-and-formatting),
which still must be checked at submission time.

The minimum SC-01 figure set is: growth/product envelope; matched MOMA and ROOM
growth-versus-product scatter panels with the five highest-ranked candidates labelled in each;
OptKnock/RobustKnock maximum-versus-guaranteed product; independent top-10 FSEOF and top-10
FVSEOF trajectories; candidate-reaction-flux-versus-target-product response curves for every
representable canonical single-knockout and amplification candidate; and paired
wild-type/knockout sampling
distributions for every canonical single-knockout candidate. Here, the knockout universe
contains the unique blocked-reaction signatures represented by MOMA D1–D5 and ROOM D1–D5,
while the amplification universe is the union of the independent FSEOF and FVSEOF top-10
lists. A missing requested method or an otherwise non-runnable target produces an explicit
unavailable/skipped panel or status row,
not a silently shortened story. Loop-flagged or unresolved amplification candidates remain in
the response figure with their eligibility status; the flag blocks support/recommendation, not
response execution.

For Figure 5, enforced candidate-reaction flux is always the x-axis (`target_flux`) and target-
product flux is always the y-axis (`response_flux`). Growth is the configured minimum-growth
constraint and secondary `biomass_flux` output, not an axis. Amplification curves are
wild-type candidate→product scans. A knockout curve is a pre-deletion wild-type titration of
its single blocked reaction over reference↔zero when reference flux is nonzero. A zero-reference
reaction instead uses its full feasible domain and must be labelled exploratory, not causal
deletion support. A multi-reaction signature has no numeric panel and must remain explicit
unavailable/skipped rather than silently choosing one reaction.

The amplification figure retains all available top-10 rows from each method and marks their
loop-diagnostic status. A loop flag or unresolved diagnostic blocks support/recommendation
promotion, not response execution; it does not erase the hypothesis or its trajectory from the
report.

Every plotted row must be traceable to a manifest-declared CSV in the same run directory. The
runtime R gate checks `Rscript`, required package availability, process exit status, non-empty
output, and matching PNG/PDF/SVG siblings; package dependency metadata supplies compatible
minimum-version constraints. Exact renderer versions come from restoring `renv.lock`, are
asserted against the complete lock in the three-OS CI matrix, and are recorded in the figure
manifest. Any runtime failure makes report rendering incomplete.

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

Run `validate_production_run(run_dir)` and check each of these on the standalone copy. Any
validation error means the run is not complete; warnings remain visible in the report and
handoff.

- [ ] No `figures/` or other relative `src`/`href` survives. Search the file for `src="fig`.
- [ ] Number of `<img>` equals the number of figure captions.
- [ ] Figure numbers run 1, 2, 3 … in document order, with none repeated or skipped.
- [ ] Section numbers are unique, and every `§N` cross-reference points at a section that exists.
- [ ] Opening the file from a different directory still shows every figure.
- [ ] Every manifest entry exists and stays under the run directory; validation recomputes its
  recorded SHA-256 and byte size.
- [ ] Every requested method has a result or an explicit unavailable record.
- [ ] Every unique MOMA/ROOM D1–D5 single-knockout candidate has flux-response and matched
  wild-type/knockout sampling index entries, regardless of recommendation status. A
  single-reaction signature has a pre-deletion WT scan: reference↔zero when reference flux is
  nonzero, otherwise the full feasible domain marked exploratory. A multi-reaction signature has
  an explicit unavailable/skipped reason and no silently selected x-axis reaction.
- [ ] Every candidate in the independent FSEOF/FVSEOF top-10 union has a flux-response index
  entry and attempted scan, including loop-flagged or unresolved candidates. Lethal/infeasible,
  unavailable, and failed non-runnable cases carry explicit status and reason rows rather than
  disappearing.
- [ ] Every recommended amplification target is traceable to a method-specific FSEOF or FVSEOF
  ranking and has a completed loopless diagnostic with no artifact flag; unresolved diagnostics
  remain partial.
- [ ] Every figure has source data plus both its 300-DPI PNG and editable PDF/SVG sibling.
- [ ] No missing R package, nonzero renderer exit, clipped label, empty required panel, or
  non-finite plotted value remains.

### Sections

Number every section and cross-reference as `§4`. Sections 1–4, 6 and 7 are required; §5 is
the scenario's own answer and its heading changes with the scenario.

| # | Section | Content |
|---|---|---|
| 1 | **Summary** | The findings themselves, not a description of the work. Three to five sentences: what was found, what is recommended, and the single most important caveat. A reader who stops here must not be misled. |
| 2 | **Setup** | The preflight summary table from `_preflight.md`, plus the confirmed model/product/condition, aeration as bounds, substrate and uptake, solver/R capabilities, and any unavailable or explicitly substituted method with its reason. |
| 3 | **Data and methods** | Model and its fingerprint, sizes, every parameter that changes an answer — reference flux state, `n_steps`, strain-design seed, sampler seed and count, thresholds. Enough that someone repeats the run without reading the scenario. |
| 4 | **Results** | One subsection per pipeline step, in order. Each states what was run, the decisive numbers, the figure, and the CSV they came from. |
| 5 | *scenario-specific* | SC-01 → **Recommended targets and strain proposal**; SC-02 → **Interpretation — what explains the difference**; SC-03 → **Essentiality classes**. See below. |
| 6 | **Limitations** | What the analysis does *not* establish. At minimum: predictions needing experimental validation; the medium and aeration assumed; conflicting method assumptions (MOMA's minimal adjustment vs OptKnock's growth maximization); solver capability that constrained the run. |
| 7 | **Provenance** | Model fingerprint, solver and version, CMM version or git commit, strain-design seed, sampler seed, run directory, and the command or script that produced it. |

**§5 for SC-01** — keep three evidence tables rather than forcing unlike methods into one
score: every canonical single-knockout candidate (MOMA/ROOM display ranks, growth, product,
response, paired-sampling verdict, and GPR-resolved reaction ids/names), strain designs
(OptKnock/RobustKnock, knockout set, growth, maximum and guaranteed product, coupling verdict),
and every canonical amplification candidate (method-specific FSEOF or FVSEOF rank, direction
and robustness, loopless-capacity verdict, wild-type-to-supported flux range, response
verdict). Each row links to its source and states `support`, `contradict`, `inconclusive`,
`skipped`, or `unavailable`; never turn method count alone into a confidence claim. Tables and
figures may paginate or facet these rows for readability, but must not select a smaller
scientific subset.

`recommendations.csv` contains only claims that pass the declared evidence policy in its
metadata sidecar: beneficial single-gene prediction plus paired-sampling support and retained
growth. A nonzero-reference WT titration may be supporting mechanism evidence; a zero-reference
full-domain scan is exploratory and not causal deletion support. A candidate from either FSEOF
or FVSEOF with supporting response and a
non-flagged loopless diagnostic for amplification; or RobustKnock positive guaranteed product
with retained growth for a multi-knockout. Agreement between FSEOF and FVSEOF is reported but
is not required. OptKnock remains visible but is not worst-case evidence. Combined
knockout-amplification proposals are withheld because this workflow does not simulate the
combined intervention.

For amplification, the wild-type flux and supported response range are not decoration. A
target whose response is flat, whose proposed direction is contradicted, or whose useful point
sits at an artificial bound is not presented as a recommendation however highly an inverse
scan ranked it. For knockouts, the sampling column compares the matched wild-type and deletion
ensembles; a wild-type-only sampling plot is not validation.

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
   the aeration as bounds, substrate, solver, strain-design seed, and sampler seed/count. These
   belong in §3, where a reader looks for them, not scattered through §4. Never leave the
   strain-design seed implicit: a hidden backend seed can change MILP paths and runtime.
4. **Say what is a hypothesis.** FSEOF, rMTA and sampling-based rankings prioritise candidates;
   they do not prove them. §6 is not optional.
5. **Do not fabricate a run.** If a step failed or was skipped, say so and why, in place. A
   report that omits a failed step reads as though the step succeeded.
6. **Ship both copies.** A run that leaves only `report.html` cannot be sent to anyone, and the
   failure is invisible until they open it.
7. **Validate the deliverable.** A report that renders but fails `validate_production_run` is
   incomplete. Keep errors, warnings, and unavailable methods visible; do not delete them to
   obtain a cleaner narrative.
