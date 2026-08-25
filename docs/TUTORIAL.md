# CMM Tutorial — Constraint-based Metabolic Modeling

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
┌──────────────────────────────────────────────────────────────────────────────┐
│ CMM — Constraint-based Metabolic Modeling (menu: File/Analysis/Model/Config) │
│ model id · reactions · metabolites · genes · solver                          │
├─────────────────────┬────────────────────────────────────────────────────────┤
│ MODEL PANEL         │ TABS                                                   │
│ - objective         │ Simulation | Comparison | Production | Strain Des.     │
│ - flux-range slider │ Omics | Multi-condition | Flux Map |                   │
│ - reaction table    │ Revert Metabolism | Transform (A→B)                    │
│  (edit Lower/Upper) │   (active analysis + result table / figure)            │
├─────────────────────┴────────────────────────────────────────────────────────┤
│ status bar                                                                   │
└──────────────────────────────────────────────────────────────────────────────┘
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

## 5c. Flux Response tab — how far does one reaction carry the network?

Fixes the **target** reaction's flux at each point of a scan and maximizes the **response**
reaction there, so you see a curve instead of a single operating point.

1. Pick a **Target (scanned)** reaction — the one you would over-express or throttle.
2. Pick a **Response (maximized)** reaction. Leave it on `(objective)` to ask "how sensitive
   is growth to this reaction, and where does it break?"; pick a product exchange to ask
   "how much product does this reaction's flux buy me?"
3. Set **Min growth (% of wild type)** when the response is a product. At 0 the scan
   maximizes product with no growth floor, which returns a non-growing theoretical ceiling;
   at 30 the curve describes a cell that still grows at 30% of wild type.
4. **Range min**/**max** is the scanned interval. It is filled automatically whenever you pick
   a target with the full flux interval that reaction can carry in the current medium — found
   by FVA with the objective unconstrained — so it never starts from a placeholder or the
   previous reaction's numbers. Edit either box to scan a narrower window, and press
   **Detect range** to reset to the detected interval. You may set values beyond the
   reaction's own declared bounds as a what-if; the result records that under
   `range_outside_bounds`.
5. Click **Run scan**. The curve shows the wild-type flux, the optimum, any shaded range that
   has no solution at all, and the **response limit** — the phase boundary past which the
   shadow price `d(response)/d(target)` turns against you — if there is one. The table beside
   it lists every scan point; infeasible points are shown in red rather than dropped.

A response that never declines means the target does not limit the response over that range.
A flat response means the objective is simply insensitive to that reaction — both are
findings, not failures.

## 5d. Sampling tab — is a predicted flux forced, or one of many?

FBA reports one optimal solution out of many that are equally optimal, so a single predicted
flux cannot tell you whether the network *requires* that value. Sampling answers that.

### What each option controls

| Option | What it does | When to change it |
|---|---|---|
| **Mode** | `uniform` samples the whole feasible space defined by the model and medium — "what can this network do at all?". `around a reference` first narrows every reaction to a window around one predicted flux state, then samples inside that — "how much could this *prediction* vary?" | Use `uniform` to test whether a predicted flux is forced. Use `around a reference` to put an uncertainty range on a specific prediction. |
| **Reference** | Which predicted flux state the windows are centred on: `pfba` (the unique minimal-total-flux solution) or `fba`. CMM computes that state for you from the current model and medium — you are choosing the *method*, not supplying numbers. **Greyed out in `uniform` mode**, which has no reference. | `pfba` unless you specifically want the plain FBA solution. |
| **window ±%** | The half-width of each reaction's window, as a percentage of its reference flux. 20 confines each reaction to ±20% of its predicted value, intersected with its existing bounds; reactions at essentially zero flux get a small window around zero instead. **Greyed out in `uniform` mode.** | Narrow (5–10) to ask how tightly the prediction is determined; wide (30–50) to explore around it. |
| **Samples** | How many flux distributions to draw. | ≥1000 for `optgp`. Fewer is faster but under-represents the space. |
| **Sampler** | `optgp` (parallel-capable hit-and-run, needs a large sample count to mix well) or `achr` (artificial centering hit-and-run, single process, better convergence at small counts). | `achr` for quick runs under ~1000 samples; `optgp` for large ones. |
| **Thinning** | Keep every *n*-th iterate of the chain. Higher values give less correlated samples for the same count, at more compute. | Leave at 100. Lower it only when a run is too slow and you accept more correlation. |
| **Seed** | Initializes the random chain. The same seed, sampler, thinning, and sample count reproduce the ensemble exactly. | Change it to check that your conclusion is stable across seeds — a conclusion that moves with the seed is not converged. |
| **Show top** | How many of the most variable reactions to draw in the violin plot. **This affects the figure only, not the table or the exported data.** | Raise it to see more reactions; the plot gets taller. |
| **Knockouts (optional)** | Sample a deletion strain instead of the wild type. Pick a **level** (`reaction` deletes the listed reactions; `gene` resolves each gene to its reactions through the GPR) and move targets into the chosen set. Leave the set empty for wild type. | Use it to ask what a mutant's flux space looks like — which reroutings open up, and which fluxes become forced. |

### Running it

Click **Run sampling**. Then:

- The **violin plot** shows the `Show top` most variable reactions — the informative ones, since a
  reaction with near-constant sampled flux is exactly the one a single FBA already settled.
- The **table lists every reaction in the model**, sorted by standard deviation descending, with
  mean, std, min, median, and max. It is not a subset.
- **Export samples…** writes the raw ensemble to CSV: one row per sample, one column per
  reaction. Use this for your own statistics; the table's summary is available separately
  through *File → Export Table to CSV*.

A reaction with a wide sampled distribution is one whose FBA value was an arbitrary choice among
alternate optima; a narrow one is genuinely pinned by the constraints. Watch for reactions
spanning hundreds of flux units — those are usually thermodynamically infeasible loops
(`FRD7`/`SUCDi` in `e_coli_core`), not biology. Sampled means can be reused as a reference flux
state for MOMA/ROOM/MTA.

Knockouts are applied as a scoped condition, so the loaded model is never edited — the bounds
you see in the left panel are unchanged after a run. If a knockout set leaves no feasible space,
the tab reports it as a probably-lethal set rather than an error, and clears the previous
ensemble so a stale result cannot be exported under the new settings. In `around a reference`
mode the reference is computed **with the knockouts applied**, since a wild-type reference would
place every deleted reaction outside its own sampling window.

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

The **Flux Map** tab draws the current flux distribution on the network. Click **draw: FBA**
or **pFBA** to solve and draw in one step; colour and width encode flux (diverging: blue = reverse/negative, red =
forward/positive; width ∝ `|flux|`).

**Which flux state is drawn.** The map draws the window's current distribution, and the figure
title names the method that produced it — `e_coli_core — pFBA`, `e_coli_core — LAD · condB`.
Rendering with no distribution loaded runs FBA first; editing a bound marks the distribution
stale, so the next render re-solves. Four things can put a distribution on the map:

| From | How |
|---|---|
| FBA, pFBA | **draw: FBA** / **pFBA** on this tab — solves and draws in one click |
| E-Flux2 / LAD | Omics tab → pick the condition → **Show on flux map** |
| MOMA / ROOM | Comparison tab → run a single knockout → **Show on flux map** |

Switching the layout, or changing the reaction count, redraws the loaded flux state by itself.
Those controls never solve — that is the point of keeping them separate from **draw:**, which
always does. Pressing **FBA** merely to see the schematic would replace an omics or knockout
flux state you had drawn with a fresh FBA solve.

Batch comparison is not on this list: it produces one distribution per target, so there is no
single "the" flux state to draw. Run the target you care about as a single knockout.

The reaction table on the left switches to the same distribution, so the map and the numbers
never disagree.

Two layouts, chosen with the **layout** selector:

| Layout | What it draws | Needs |
|---|---|---|
| **Escher map (curated)** | A published Escher map's hand-laid coordinates and bezier segments — every reaction the map covers | an Escher map JSON |
| **Schematic (top reactions)** | The *n* highest-`|flux|` reactions (4–25), currency metabolites (ATP, H₂O, CO₂, NAD(P)H…) dropped so the layout follows the carbon skeleton. Arrow colour *and* width scale with `|flux|`, read against the colorbar; the arrow points the net direction | nothing |

CMM never invents a map layout. A readable metabolic map is hand-drawn, and an automatic
layout of a genome-scale network is a hairball rather than a figure — which is why the
schematic draws a handful of reactions and says so in its title rather than pretending to
show the whole network.

The reaction count is capped at 25 for the same reason. The schematic folds the carbon
backbone into evenly spaced rows, but a metabolic network branches, and past roughly that many
reactions the arrows between rows dominate the picture. Reading more of the network at once is
what a curated map is for.

**The bundled map.** CMM ships Escher's *E. coli* core map and offers it automatically to any
model containing at least half its reactions. That covers `e_coli_core` (94 of 95) and also
genome-scale reconstructions such as iJO1366 — viewing a genome-scale model on a
central-metabolism map is how Escher itself is used. The caption above the figure says which
map is loaded and how much of it this model can fill. Provenance and license are in
`src/cmm/resources/ATTRIBUTION.md`; cite King *et al.* (2015) for the layout.

**Your own map.** **Map JSON…** on the tab, or **File ▸ Open Escher Map (JSON)…**, takes any Escher
map JSON (escher.github.io or BiGG). A file whose reactions do not appear in the loaded model
is refused with a message rather than drawn as a blank grey map.

**A drawing under the flux.** CMM reuses a map's layout but not Escher's drawing conventions —
its arrowheads, node sizes and typography. **Map image…** loads a picture of the same map and
lays it beneath the flux colouring, so the network is drawn the way Escher draws it and CMM
contributes only the colour and width. Export one from escher.github.io with *Map ▸ Export as
SVG*; PNG and JPEG work too. The picture is placed in the map's own canvas coordinates, so an
export of the map you loaded lines up without adjustment. While a background is shown CMM
leaves out its own labels and metabolite dots — the drawing already has them — and **show**
hides the drawing again without unloading it. A background belongs to one map, so loading a
different map or model drops it.

Both layouts are in the Python API as `escher_flux_map` and `network_flux_map`
(`cmm.visualization`), and `cmm.resources.bundled_map_for(model)` returns the bundled map's
path when it suits a model.

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
from cmm.features.response import flux_response
from cmm.features.sampling import random_flux_sampling, reference_constrained_sampling
from cmm.features.strain_design import optknock, robustknock
from cmm.features.transformation import transformation_targets
from cmm.omics.expression import integrate_expression
from cmm.omics.conditions import predict_condition_fluxes, flux_log_change

model = load_model("textbook")
print(solver_status(model).summary())

# Simulation
sol = fba(model);  print(sol.objective_value, sol.status)

# Production design
y = theoretical_yield(model, "EX_succ_e")     # condition comes from the applied medium
# 1.6384 against a ceiling of 1.5; co2_carbon_fraction says how much of the product carbon
# came from CO2 uptake (8.4% here), which is the number `co2_fixed` alone never gave.
print(y.molar_yield, y.carbon_ceiling, y.co2_carbon_fraction, y.carbon_imbalance)

# Growth-coupled strain design (needs a MILP solver)
result = optknock(model, "EX_succ_e", max_knockouts=3, max_solutions=5)
for d in result.designs:
    print(d.knockouts, d.growth, d.guaranteed_product)

# Omics → flux
expr = {g.id: 50.0 for g in model.genes}
flux_state = integrate_expression(model, expr, method="eflux2").to_flux_state()

# Perturbation response against a real reference. `distance` is Segre et al. Eq. (4)'s
# Euclidean distance; `objective_value` is the raw QP objective (its square), and for ROOM
# `distance` is None because a switch count is not a distance -- see `n_changed_reactions`.
ref = reference_flux(model, "pfba")
with model:
    model.reactions.PFK.knock_out()
    r = moma(model, ref, linear=False)
    print(r.distance, r.distance_kind, r.objective_value)   # 51.321 euclidean_l2 2633.846

# Gene / multi / batch knockouts
gene_rxns = blocked_reactions_for_genes(model, ["b0726"])          # gene -> reactions (GPR)
k = knockout_comparison(model, ref, gene_rxns, method="moma_l2")
print(k.distance, k.objective_value)                       # 11.398 distance, 129.925 objective
print(knockout_comparison(model, ref, ["PFK", "TPI"], method="moma_l2").distance)  # multi-KO
batch = batch_comparison(model, ref, gene_perturbations(model), method="moma_l2")
for row in sorted(batch, key=lambda r: -r.distance)[:5]:            # most-disrupted first
    print(row.target_id, row.status, round(row.distance, 3), round(row.objective, 3))
batch.to_frame().to_csv("batch_screen.csv", index=False)   # rows = the numbers
print(batch.metadata["model_sha256"], batch.metadata["n_inert_dropped"])  # container = the run

# Multi-condition omics comparison (log2 fold-change of flux magnitude)
# preds = predict_condition_fluxes(model, expression_dataframe, method="eflux2")
# lc = flux_log_change(preds.fluxes("condA"), preds.fluxes("condB"))

# Verifying a target: does forcing flux through it buy product, and where does it break?
resp = flux_response(model, "PGI", "EX_succ_e", biomass_fraction=0.3, n_steps=20)
print(resp.optimum(), resp.feasible_range(), resp.limit.found)
resp.to_frame().to_csv("flux_response.csv", index=False)

# Is a predicted flux forced, or one of many alternate optima?
ens = random_flux_sampling(model, n=1000, seed=0)        # add achr for small runs
print(ens.statistics().loc[["EX_succ_e", "PGI"]])
near = reference_constrained_sampling(model, ref, n=1000, seed=0)   # around a prediction
sampled_reference = ens.to_flux_state()                   # reusable MOMA/ROOM/MTA reference
```

Export publication figures directly:

```python
from cmm.visualization import production_envelope_figure, save_figure
env = production_envelope(model, "EX_succ_e", points=20)
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

The captures are not tracked in git — they are reproducible from the commands above, and a
committed copy goes stale as soon as the code changes.
[Scenario figures](scenario-figures.md) is the manifest: it names every figure the three
harnesses write and what each one shows.

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

*Feature availability reflects CMM 0.5.0. See `docs/VALIDATION.md` for the publication
evidence and limitations, `docs/feature-roadmap.md` for planned additions, and
`docs/architecture.md` for the layering contract.*
