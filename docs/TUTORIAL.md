# CMM Tutorial — Cellular Metabolic Modeling Platform

A hands-on guide to every feature of the CMM desktop platform and Python library, written
against the shipped `e_coli_core` textbook model so you can follow along end to end. Each GUI
step below was exercised in offscreen mode and produces the screenshots referenced under
`temp_figures_new/`.

> **Scope note.** The GUI now covers every advertised method — simulation, production design,
> **strain design (OptKnock/RobustKnock)**, omics, **multi-condition comparison**,
> perturbation response, revert-metabolism, the **A→B transformation finder**, and flux maps —
> and can **export result tables to CSV and save figures** (File menu + the figure toolbar).
> The Python API (§10) exposes the same services for scripting.

---

## 1. What CMM does

CMM is built on [COBRApy](https://opencobra.github.io/cobrapy/). It answers three questions:

1. **How does the cell distribute flux?** — FBA, pFBA, FVA, growth media, omics integration.
2. **How do we make more of a product?** — theoretical yield, production envelope, FSEOF,
   FVSEOF, OptKnock/RobustKnock.
3. **Which knockouts move one metabolic state toward another?** — MOMA/ROOM perturbation
   response and revert-metabolism (rMTA).

Every analysis runs through a solver-neutral service layer (`cmm.core`, `cmm.features`,
`cmm.omics`, `cmm.visualization`); the GUI (`cmm.app`) only renders. That means anything you
can do in the window you can also script, and results are reproducible.

---

## 2. Install & launch

```bash
git clone https://github.com/jyryu3161/CMM.git && cd CMM
./install.sh                       # macOS / Linux  (.\install.ps1 on Windows PowerShell)

# launch on the built-in textbook model
.venv/bin/python -m cmm.app        # Windows: .venv\Scripts\python -m cmm.app

# or launch on your own SBML model
.venv/bin/python -m cmm.app path/to/model.xml
```

Gurobi or CPLEX unlocks the full feature set. GLPK supports LP and MILP, so it can run
FBA/pFBA/FVA/LAD/production scans plus ROOM and StrainDesign workflows. L2 MOMA and E-Flux2
need QP; published MTA/rMTA need MIQP. The restricted Gurobi license is sufficient for the
small validation models but can be too small for genome-scale mixed-integer analyses. Check
the active capabilities under **Config → Solver status…**.

---

## 3. The window at a glance

```
┌───────────────────────────────────────────────────────────────────────┐
│ CMM — Cellular Metabolic Modeling Platform (menu: File/Analysis/Model/Config)│
│ model id · reactions · metabolites · genes · solver                        │
├──────────────────────┬────────────────────────────────────────────────────┤
│ MODEL PANEL          │ TABS                                               │
│ - objective          │ Simulation | Comparison | Production | Strain Des. │
│ - flux-range slider  │ Omics | Multi-condition | [Flux Map] |             │
│ - reaction table     │ Revert Metabolism | Transform (A→B)                │
│   (edit Lower/Upper) │   (active analysis + result table / figure)        │
├──────────────────────┴────────────────────────────────────────────────────┤
│ status bar                                                                 │
└────────────────────────────────────────────────────────────────────────────┘
```

- The **reaction table** on the left lists every reaction with its Lower/Upper bounds and the
  most recent flux. Double-click a **Lower** or **Upper** cell to edit a bound, then re-run
  FBA. Bounds are clamped so `lower ≤ upper` always holds; the status bar tells you when a
  value was clamped.
- The **flux-range slider** highlights reactions whose `|flux|` is at least the chosen
  fraction of the maximum flux (run FBA first). It's a quick way to see the active backbone.
- The **Flux Map** tab only appears when the window was given a curated Escher map (see §9).
- **Saving results.** The **File** menu has **Open Model…** (load a different SBML/JSON model
  without restarting), **Export Table to CSV…** (writes the current tab's result table), and
  **Save Figure…** (saves the current Production / Flux Map figure at 300 DPI). Production and
  flux-map figures also carry a matplotlib toolbar for interactive zoom / pan / save.
- **Responsiveness.** Heavy analyses run on a background thread with a busy indicator, so the
  window keeps repainting instead of freezing; input is blocked until the run finishes.

---

## 4. Simulation tab — FBA, pFBA, FVA, media

This is the starting point for any model.

1. **Pick a medium.** The *Medium* dropdown offers presets (`glucose_aerobic`,
   `glucose_anaerobic`, `acetate_aerobic`, `glycerol_aerobic`). Click **Apply medium** — it
   opens the listed uptakes and closes every other exchange, then refreshes the table. Media
   are matched tolerantly across id conventions (`EX_glc__D_e` / `EX_glc_e` / `EX_glc(e)`).
2. **Run FBA.** The **Objective** label shows the optimum and status; the reaction table's
   *Flux* column and the *Simulation* result table both populate. On `e_coli_core` glucose
   aerobic the growth rate is **0.8739 h⁻¹**.
3. **Run pFBA** for the unique minimal-total-flux distribution at the same growth (the status
   bar reports the total `|flux|`, ≈518.4 on the textbook model).
4. **Run FVA.** Set the **fraction** spin box (default 0.90 = hold 90 % of the optimum), then
   click **Run FVA**. Each reaction gets a `[min, max]` feasible range. FVA auto-runs FBA
   first if fluxes are stale.

**Editing bounds.** To force anaerobic growth by hand, double-click the **Lower** cell of
`EX_o2_e`, set it to `0`, and re-run FBA. The status bar flags stale fluxes after any edit.

---

## 5. Production tab — making more of a target

Select a **Target product** (any exchange reaction, e.g. `EX_succ_e` for succinate), a
**substrate** (`auto` detects the limiting carbon uptake), and **aerobic/anaerobic**.

- **Theoretical yield** — maximum mol product / mol substrate at fixed uptake. The label
  discloses the substrate **carbon ceiling** and flags **net CO₂ fixation** when the yield
  exceeds it (so a number above the carbon ceiling is never misread). Example: aerobic
  succinate ≈ **1.638 mol/mol glucose**, annotated *needs net CO₂ fixation*.
- **Production envelope** — the growth-vs-product phenotypic phase plane, feasible region
  shaded. If min-growth is zero across the range the title says so.
- **FSEOF targets** — scans enforced product flux and classifies each reaction as amplify /
  knockdown by its flux-magnitude trend. The plot shows the top amplification targets; the
  label lists them (e.g. FRD7 / FUM / PPC for succinate). If the product's theoretical yield
  is zero the tool tells you instead of drawing an empty plot.
- **FVSEOF (robust)** — FSEOF + FVA at each step. A target is *robust* when its **forced
  minimum** `|flux|` also rises (the reaction cannot avoid carrying more flux). Solid line =
  mean flux, dashed = forced minimum. This separates genuinely-forced targets from ones that
  merely *can* increase.

**A worked growth-coupled design (by hand).** Switch to anaerobic and, on the Simulation tab,
set `EX_o2_e`, `EX_ac_e`, `EX_etoh_e`, `EX_for_e`, `EX_lac__D_e` lower/upper to `(0,0)` (block
the competing fermentation routes). Re-run FBA: succinate excretion rises from **0 → ~9.9**
while growth stays near zero — growth-coupled, not a free win.

---

## 5b. Strain Design tab — OptKnock / RobustKnock

Searches for a small reaction-knockout set that **couples product to growth** — at maximum
growth the cell is forced to make the product.

1. Pick a **Target product** exchange, a **Method** (`optknock` = maximize product at max
   growth, optimistic; `robustknock` = keep only designs that *guarantee* product at max
   growth, worst case), **max KOs** (design size), and **solutions** (how many to enumerate).
2. Click **Run design**. Each design row shows its knockouts, the growth rate, the optimistic
   **Max product**, and the worst-case **Guaranteed product**; growth-coupled designs (nonzero
   guaranteed product) are highlighted. For succinate on `e_coli_core`, OptKnock finds designs
   such as `{CO2t, PGI}` (growth ≈0.20 h⁻¹, guaranteed succinate ≈9.6).

This search is a nested MILP (delegated to the `straindesign` package) and **needs a MILP
solver**. OptKnock and RobustKnock are distinct module types; RobustKnock is the three-level
worst-case formulation. It can take a while on larger models.

## 6. Omics tab — expression → flux (E-Flux2 / LAD)

1. Choose a **Method**: `eflux2` (scale reaction bounds by normalized expression, then
   minimize L2 flux) or `lad` (fit fluxes to expression-derived targets, an LP).
2. **Load expression CSV…** — a two-column file, header row, **first column = gene id, second
   column = expression value**:

   ```csv
   gene,expression
   b0008,42.7
   b0114,10.1
   b0116,88.0
   ```

   Gene ids are matched case-insensitively against the model's genes. Or click **Run on demo
   expression** for a deterministic synthetic run over the model's genes.
3. The table lists active predicted fluxes (largest first); the summary reports how many
   reactions were mapped through the GPR and the achieved biological objective. E-Flux2 is
   strict: without QP support it raises a capability error. The Python API exposes an
   explicit `allow_l1_fallback=True` approximation labeled `eflux2_l1_fallback`; it must not
   be reported as E-Flux2.

---

## 6b. Multi-condition tab — compare predicted fluxes across conditions

Predicts a flux distribution for each condition in an expression table, then compares two.

1. **Load expression table CSV…** — a **gene × condition** table: first column = gene id, each
   remaining column = one condition:

   ```csv
   gene,glucose,acetate
   b0008,42.7,15.1
   b0114,10.1,80.0
   ```
2. Pick a **Method** (`eflux2` / `lad`) and the two conditions **A** and **B**.
3. **Compare.** The table lists reactions ranked by `log2( |flux_B| / |flux_A| )` — positive =
   higher flux in B. On `e_coli_core` a glucose-vs-acetate table surfaces the expected
   fermentation/TCA shifts (PFL, PDH, ACKr, …). Export the table with **File → Export Table to
   CSV…**.

## 7. Comparison tab — perturbation response (MOMA / ROOM)

Predicts the flux state after a knockout as the one closest to a reference template.

1. **Method**: `MOMA (L2)` (QP), `MOMA (L1)` (LP), or `ROOM` (MILP).
2. **Reference template**: the wild-type flux the perturbed state is compared against —
   `fba`, `pfba`, `lad`, or `eflux2`.
3. **Knockout level**: `reaction` or `gene` (a gene knockout is resolved to the reactions it
   disables through the model's GPR — a multi-gene selection is resolved jointly, so a complex
   that needs two knocked-out subunits is blocked correctly).
4. **Knock out (select one or more):** pick targets in the list — Ctrl/Shift-click to select
   several.
5. **Run (selected as one KO)** knocks out all selected targets *together* (single- or
   multi-knockout). The table shows every changed reaction (reference vs perturbed flux); a
   lethal knockout is reported as *infeasible* rather than crashing. For a gene knockout the
   summary also reports how many reactions it blocked.
6. **Batch (each separately)** runs MOMA/ROOM once per target as its *own* single knockout —
   the CNApy-style batch deletion / essentiality scan — and fills a table of *Target, Kind,
   #reactions, Status, Distance, Objective* sorted most-disrupted first (with no selection it
   scans every gene/reaction of the chosen level). Export it with **File → Export Table to
   CSV…**. On genome-scale models this is many solves; it runs in the background (§3).

> **LAD / E-Flux2 templates.** These use the gene expression you loaded on the **Omics** tab
> (**Load expression CSV…**). If no expression is loaded, the tab falls back to synthetic
> demo data and **says so in the result summary** (a ⚠ note) so it is never mistaken for a
> data-driven template — load a CSV first for a real comparison, or use the `fba`/`pfba`
> templates.

---

## 8. Revert Metabolism tab — normalization targets (rMTA)

Ranks gene/reaction knockouts that move a **source** (e.g. disease) state toward a **target**
(e.g. healthy) state, derived from two-state differential expression.

1. Set **Method** (`rmta`, `mta`, or `rmta_continuous`), **Knockout level** (`gene` /
   `reaction`), and the **transformation weight α** (0–1, default 0.66). `rmta` is the
   published best/MOMA/worst workflow; `mta` is the published single MTA solve; the continuous
   option is an explicitly labeled historical heuristic.
2. **Load source CSV/TSV…** and **Load target CSV/TSV…** — each is a two-column
   gene/expression file (same format as §6). Both must load before **Run Revert** enables.
3. **Run Revert.** The table ranks targets by score (best row highlighted); the summary names
   the top normalization target. The source reference is generated from the **source**
   expression with E-Flux2 at full objective. Per-reaction desired directions come from the
   source→target expression change and GPR logic.

Published `rmta` and `mta` need MIQP; `rmta_continuous` needs QP. On an unsupported solver
the tab reports the capability error cleanly. The original studies use contextualization and
sampling for source-state preprocessing; the deterministic E-Flux2 GUI variant must be
disclosed in a manuscript (see `docs/design-revert-metabolism.md`).

---

## 8b. Transform (A→B) tab — knockouts that move state A toward state B

A generalization of Revert Metabolism: instead of a differential-expression *direction*, it
works from two explicit predicted flux **states**.

1. **Load source (A)** and **target (B)** expression CSV/TSV files (two-column gene/expression,
   as in §6).
2. Choose **Predict states with** (`eflux2` / `lad`) — each expression vector becomes a flux
   state — a **Method** (`moma` = rank by how far the minimal-adjustment state moves toward B;
   `mta` = the published single MTA MIQP on the A→B direction), the **Knockout level**, and α
   (for `mta`).
3. **Run transformation.** Knockouts are ranked (best highlighted) by how well they move A
   toward B. `moma` needs QP; `mta` needs MIQP.

## 9. Flux Map tab — Escher-layout flux maps

The **Flux Map** tab appears only when the window is constructed with a curated Escher map
JSON (`map_path`). It reuses the map's hand-laid coordinates and bezier segments and colours
each reaction by the current FBA flux (diverging: blue = reverse/negative, red =
forward/positive; width ∝ `|flux|`). Click **Render flux map** after running FBA. A textbook
map for `e_coli_core` is bundled under `test_data/`.

To launch with the map wired in, the scenario harnesses pass `map_path=` to
`CmmMainWindow` (see §11). A dependency-free schematic `network_flux_map` is also available in
the Python API when you don't have a curated map.

---

## 10. Python API

Everything in the GUI is a thin call over the same solver-neutral services, so any workflow is
also scriptable and reproducible.

```python
from cobra.io import load_model
from cmm.core import fba, pfba, fva, apply_medium, solver_status
from cmm.features.production import theoretical_yield, production_envelope, fseof, fvseof
from cmm.features.comparison import (
    moma, room, reference_flux, knockout_comparison, batch_comparison,
)
from cmm.features._perturbation import gene_perturbations, blocked_reactions_for_genes
from cmm.features.strain_design import optknock, robustknock
from cmm.features.transformation import transformation_targets
from cmm.omics.expression import integrate_expression
from cmm.omics.conditions import predict_condition_fluxes, flux_log_change

model = load_model("textbook")
print(solver_status(model).summary())

# Simulation
sol = fba(model);  print(sol.objective_value, sol.status)

# Production design
y = theoretical_yield(model, "EX_succ_e", aerobic=True)
print(y.molar_yield, y.carbon_ceiling, y.co2_fixed)

# Growth-coupled strain design (needs a MILP solver)
result = optknock(model, "EX_succ_e", max_knockouts=3, max_solutions=5)
for d in result.designs:
    print(d.knockouts, d.growth, d.guaranteed_product)

# Omics → flux
expr = {g.id: 50.0 for g in model.genes}
flux_state = integrate_expression(model, expr, method="eflux2").to_flux_state()

# Perturbation response against a real reference
ref = reference_flux(model, "pfba")
with model:
    model.reactions.PFK.knock_out()
    print(moma(model, ref, linear=False).distance)

# Gene / multi / batch knockouts
gene_rxns = blocked_reactions_for_genes(model, ["b0726"])          # gene -> reactions (GPR)
print(knockout_comparison(model, ref, gene_rxns, method="moma_l2").distance)   # ~129.9
print(knockout_comparison(model, ref, ["PFK", "TPI"], method="moma_l2").distance)  # multi-KO
batch = batch_comparison(model, ref, gene_perturbations(model), method="moma_l2")
for row in sorted(batch, key=lambda r: -r.distance)[:5]:            # most-disrupted first
    print(row.target_id, row.status, round(row.distance, 3), round(row.objective, 3))

# Multi-condition omics comparison (log2 fold-change of flux magnitude)
# preds = predict_condition_fluxes(model, expression_dataframe, method="eflux2")
# lc = flux_log_change(preds.fluxes("condA"), preds.fluxes("condB"))
```

Export publication figures directly:

```python
from cmm.visualization import production_envelope_figure, save_figure
env = production_envelope(model, "EX_succ_e", aerobic=True, points=20)
save_figure(production_envelope_figure(env, title="Succinate envelope"), "envelope.png")  # 300 DPI
```

---

## 11. Headless / offscreen testing (the "screen-off" mode)

The project ships three scenario harnesses that drive the real GUI with the Qt **offscreen**
platform and save PNG captures — this is how the platform is tested without a display:

```bash
QT_QPA_PLATFORM=offscreen CMM_OUTPUT_DIR=./temp_figures_new PYTHONPATH=src \
  .venv/bin/python -m cmm.app.screenshots            # branched demo: FBA·FVA·slider·rMTA·MIQP
QT_QPA_PLATFORM=offscreen CMM_OUTPUT_DIR=./temp_figures_new PYTHONPATH=src \
  .venv/bin/python -m cmm.app.succinate_scenario     # e_coli_core: yield·envelope·FSEOF·Escher
QT_QPA_PLATFORM=offscreen CMM_OUTPUT_DIR=./temp_figures_new PYTHONPATH=src \
  .venv/bin/python -m cmm.app.genome_scale_scenario [model.xml]   # your genome-scale model
```

Run the unit + scenario test suite:

```bash
uv run pytest -q -ra --strict-markers
uv run ruff check src tests
```

---

## 12. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| "active solver … does not support QP/MILP/MIQP" | Check the method table in §2. GLPK provides LP/MILP; install/configure a QP/MIQP backend for L2 MOMA, E-Flux2, MTA, or rMTA. |
| Theoretical yield raises "no uptake capacity" | The chosen substrate exchange is closed (lower bound 0). Open its uptake or pick another substrate. |
| FSEOF/FVSEOF says "not meaningful: yield is zero" | The product can't carry flux in the current medium — open its exchange / check reachability. |
| Comparison with lad/eflux2 template shows a ⚠ synthetic note | No expression loaded — load a CSV on the Omics tab (§6/§7), or use the fba/pfba templates. |
| Knockout reported "infeasible" | The perturbation is lethal (e.g. can't meet ATP maintenance) — this is a real result, not an error. |
| A long analysis on a large model | Heavy analyses (FVA, FSEOF/FVSEOF, envelope, strain design, revert, transformation, multi-condition, MOMA/ROOM) run on a background thread with a busy indicator — the window stays responsive; input is blocked until it finishes. There is no cancel yet, so give a genome-scale run time, or script it via the API. |

---

*Feature availability reflects CMM 0.3.0. See `docs/VALIDATION.md` for the publication
evidence and limitations, `docs/feature-roadmap.md` for planned additions, and
`docs/architecture.md` for the layering contract.*
