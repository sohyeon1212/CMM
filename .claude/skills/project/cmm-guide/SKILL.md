---
name: cmm-guide
classification: capability
classification-reason: "Reference protocol an agent reads to understand what CMM can do and which analysis to pick before acting; it maps goals to functions rather than prescribing one workflow."
deprecation-risk: none
effort: low
description: |
  Agent-facing operating manual for CMM (Cellular Metabolic Modeling Platform). Read this
  BEFORE using CMM so you understand its purpose (mainly metabolic-engineering / production-
  target discovery), the analyses it offers, which function answers which question, the
  solver each needs, and the pitfalls to avoid.
  Triggers: /cmm-guide, "how do I use CMM", "what can CMM do", "which CMM analysis"
  Keywords: CMM, metabolic model, metabolic engineering, production target, overexpression,
  knockdown, knockout, strain design, flux, FBA, pFBA, FVA, FSEOF, FVSEOF, OptKnock,
  RobustKnock, MOMA, ROOM, E-Flux2, LAD, MTA, rMTA
argument-hint: "(no args) — read for orientation before working with CMM"
user-invocable: true
allowed-tools:
  - Read
  - Grep
  - Glob
  - Bash
---

# CMM — Agent Operating Protocol

**Read this first.** CMM is a constraint-based metabolic-modeling toolkit (Python API + Qt
GUI) over any `cobra` model. Its analyses are solver-neutral services in `cmm.core`,
`cmm.features`, `cmm.omics`; the GUI is only a thin view. **As an agent you drive CMM through
the Python API**, not by clicking the GUI.

**Primary purpose:** *metabolic engineering* — finding and evaluating genetic intervention
targets (over-expression, knockdown, knockout) that increase production of a target metabolite,
plus the supporting simulation, omics-integration, and perturbation-response analyses.

Before running anything: (1) know the user's **goal** (which metabolite? increase production?
predict a knockout's effect? reach a target state?), (2) check the **active solver's
capabilities** (§Solvers) because most advanced methods fail without QP/MILP/MIQP, (3) pick the
method from §Decision guide, (4) verify **preconditions** (§Pitfalls).

---

## 1. Setup (always do this first)

```python
from cobra.io import load_model, read_sbml_model
from cmm.core import solver_status, supports, apply_medium, PRESET_MEDIA

model = load_model("textbook")          # or read_sbml_model("model.xml"); any cobra model
model.solver = "gurobi"                  # set a capable solver if available (see §Solvers)

status = solver_status(model)            # .name, .capabilities, .warning
# Growth medium (optional): sets exchange uptake bounds.
apply_medium(model, "glucose_aerobic")   # keys: PRESET_MEDIA
```

- CMM never mutates the model destructively during an analysis — functions use `with model:`
  contexts and restore bounds. Apply media / bound edits deliberately.
- Objective reaction = the model's biomass reaction (nonzero `objective_coefficient`).

## 2. Decision guide — goal → analysis

| The user wants to… | Use | Module |
|---|---|---|
| Baseline growth / flux distribution | `fba`, `pfba` | `cmm.core` |
| Flux ranges at (near-)optimal growth | `fva` | `cmm.core` |
| Max theoretical yield of a product | `theoretical_yield` | `cmm.features.production` |
| Production-vs-growth trade-off curve | `production_envelope` | `cmm.features.production` |
| **Rank over/under-expression targets for a product** | `fseof`, `fvseof` | `cmm.features.production` |
| **Find knockout SETS that couple production to growth** | `optknock`, `robustknock` | `cmm.features.strain_design` |
| Predict the flux state after a specific knockout | `knockout_comparison` (MOMA/ROOM) | `cmm.features.comparison` |
| Screen many single knockouts (impact / essentiality) | `batch_comparison` | `cmm.features.comparison` |
| Turn gene expression into a flux state | `integrate_expression` (E-Flux2 / LAD) | `cmm.omics.expression` |
| Per-condition fluxes from an expression table | `predict_condition_fluxes` | `cmm.omics.conditions` |
| Knockouts that revert disease→healthy | `revert_targets` (MTA/rMTA) | `cmm.features.revert` |
| Knockouts that move state A→state B | `transformation_targets` | `cmm.features.transformation` |

Two families to keep straight:
- **Forward (predict a result):** MOMA/ROOM/batch — "you give the knockout, it predicts the
  flux response."
- **Inverse (find the intervention):** FSEOF/FVSEOF, OptKnock/RobustKnock, revert/transform —
  "you give the goal, it finds the targets."

## 3. Capability reference

### Simulation — `cmm.core`
- `fba(model)` / `pfba(model)` — LP. Max biomass; pFBA also minimizes total flux (unique-ish).
- `fva(model, fraction_of_optimum=0.9)` — LP. Flux ranges holding a fraction of the optimum.
- Returns `FluxSolution` / `FluxRange`. Use for a baseline before any design analysis.

### Production design — `cmm.features.production` (main use)
- `theoretical_yield(model, product, substrate=None, aerobic=True)` — max mol product / mol
  substrate. Reports carbon ceiling / CO₂ fixation flags.
- `production_envelope(model, product, substrate=None, aerobic=True, points=20)` — growth
  (min,max) at each enforced product flux → the growth/production trade-off.
- `fseof(model, product, n_steps=10, aerobic=True)` — **heuristic**. Enforces rising product
  flux; classifies each reaction's |flux| trend as `amplify` (over-express target),
  `knockdown` (down-regulate/delete), or `none`. `result.trends`, `.amplification_targets()`.
- `fvseof(...)` — FSEOF + FVA per step; flags **robust** targets whose forced-min |flux| also
  rises (cannot avoid carrying more flux). LP/FVA only.
- `product` is an exchange reaction id (e.g. `EX_succ_e`). The model must have exchanges.

### Perturbation response — `cmm.features.comparison`
- `reference_flux(model, method, gene_expression=None)` — build the wild-type reference:
  `"fba"`/`"pfba"` (model only) or `"lad"`/`"eflux2"` (need `gene_expression`).
- `knockout_comparison(model, reference, reaction_ids, method="moma_l2"|"moma_l1"|"room")` —
  force reactions to 0, predict the **minimally-adjusted / minimally-rerouted** flux state
  nearest the reference. Covers single, multi, and (via GPR-blocked reactions) gene knockouts.
- `batch_comparison(model, reference, perturbations, method=, product_reaction=None)` — run one
  knockout per target; returns per-target growth, distance, essentiality, optional product.
- Assumption: cell makes a *minimal adjustment* from wild-type (MOMA) or changes the *fewest*
  reactions (ROOM). Contrast with OptKnock (§below), which assumes max-growth.

### Omics integration — `cmm.omics`
- `integrate_expression(model, gene_expression, method="eflux2"|"lad")` — gene expression →
  reaction weights (via GPR) → predicted flux state. E-Flux2 scales bounds + L2-min (QP);
  LAD fits fluxes to expression (LP).
- `read_expression_table(path)` → gene × condition DataFrame; `predict_condition_fluxes(...)`
  gives one flux distribution per condition.
- `gene_expression` is `{gene_id: value}`; gene ids MUST match `model.genes` (§Pitfalls).

### Reversion / transformation — `cmm.features.revert`, `.transformation`
- `revert_targets(model, None, reference_state, direction, method="rmta"|"mta", ...)` — rank
  knockouts that move a source (disease) state toward normal, scored by a transformation score
  (how well flux moves in the desired direction). `direction` comes from two-state
  `differential_expression`. rMTA also checks the worst case (robust).
- `transformation_targets(model, source_state, target_state, method="moma"|"mta")` — generalizes
  revert to two explicit flux states A→B. `moma`: score = reduction in distance-to-B; `mta`:
  direction = B−A, scored like rMTA.
- Source/target flux states usually come from `integrate_expression` on two conditions.

### Strain design — `cmm.features.strain_design` (rigorous, growth-coupled)
- `optknock(model, product, max_knockouts=3, max_solutions=5, min_growth=0.05)` — **bilevel
  MILP**: choose ≤K reaction knockouts (outer: max product) such that when the cell maximizes
  growth (inner FBA) it is forced to make the product.
- `robustknock(...)` — max-min: maximizes the **guaranteed** (worst-case) product at max growth.
- Each `StrainDesign` reports `knockouts`, `growth`, `max_product` (optimistic) and
  `guaranteed_product` (worst-case). **`guaranteed_product > 0` ⇒ growth-coupled** (reliable).
- Needs a MILP solver AND the `straindesign` package (which needs Java/OpenJDK).

## 4. Key concepts (to avoid misuse)

- **Reference vs knockout.** Perturbation methods compare a knocked-out model to a wild-type
  *reference flux state*. Choose the reference method deliberately (fba/pfba/lad/eflux2).
- **Growth coupling.** A design is valuable only if the cell *cannot* maximize growth without
  producing — read `guaranteed_product`, not just `max_product`.
- **FSEOF vs OptKnock.** FSEOF = fast heuristic ranking of over/under-expression targets, no
  coupling guarantee. OptKnock/RobustKnock = rigorous knockout-set search that enforces
  coupling. Screen with FSEOF, design with OptKnock.
- **MOMA/ROOM vs OptKnock assumption.** MOMA/ROOM assume *minimal adjustment* from wild-type;
  OptKnock assumes the cell reaches its *growth-maximizing* state. Different cell models.
- **rMTA is optimistic per knockout** (best achievable move toward the goal) — it is a
  prioritization/hypothesis tool, robustified by the worst-case term, not a proof.
- **Direction (revert/transform mta)** is the *goal*, not the intervention: the top knockout
  need not be one of the differentially-expressed genes.

## 5. Solver requirements — CHECK BEFORE RUNNING

`glpk` (the cobra default) is LP+MILP only. Many methods need QP or MIQP.

| Method | Needs | glpk? |
|---|---|:---:|
| FBA, pFBA, FVA, theoretical_yield, production_envelope, FSEOF, FVSEOF, LAD | LP | ✅ |
| MOMA (L1), ROOM | LP / MILP | ✅ |
| MOMA (L2), E-Flux2, transform `moma`, revert `rmta_continuous` | **QP** | ❌ |
| OptKnock, RobustKnock | **MILP** (+ `straindesign` + Java) | ⚠️ package/Java |
| revert `rmta` / `mta`, transform `mta` | **MIQP** | ❌ |

QP → install `osqp` (or gurobi/cplex). MIQP → `gurobi` or `cplex`. Check with
`supports("QP", model.solver.interface)` / `solver_status(model)`. If a required capability is
missing, either switch to an LP-capable method (e.g. MOMA L1 instead of L2, LAD instead of
E-Flux2) or tell the user which solver to install — do not silently produce nothing.

## 6. Preconditions & common pitfalls

- **Expression gene ids must match `model.genes`.** GEO/`GSM*.csv` often use symbols/probes;
  map to the model's ids (e.g. E. coli `b`-numbers) or the integration is meaningless.
  Expression values must be finite and non-negative; no duplicate genes.
- **Product/strain analyses need exchange reactions.** If `model.exchanges` is empty, production
  design is unavailable.
- **Lethal knockouts** (removing an essential reaction) make the solve infeasible — expected,
  reported as `status="infeasible"` / `-inf` / essential=`yes`, not a crash. Don't treat a
  lethal target as an error.
- **Zero-flux reactions.** Knocking out a reaction that carries no flux in the reference changes
  nothing — MOMA returns the reference; don't over-interpret an "empty" result.
- **ROOM vs MOMA display.** ROOM minimizes the *count* of changed reactions, not deviation, so
  its raw flux vector drifts (alternate optima); a significant-change threshold on
  `|perturbed − reference|` (relative to |reference|) filters that noise.
- **Aerobic/anaerobic and substrate** materially change production results — set them to match
  the user's scenario.

## 7. Typical workflow — "increase production of metabolite X"

1. `apply_medium` + `fba`/`pfba` → confirm the model grows; note wild-type X flux.
2. `theoretical_yield(model, "EX_x_e")` and `production_envelope(...)` → is X reachable, and
   what is the growth/production trade-off?
3. `fseof`/`fvseof(model, "EX_x_e")` → shortlist over-expression (amplify) and knockdown/
   knockout (knockdown) targets; prefer FVSEOF-robust ones.
4. `optknock`/`robustknock(model, "EX_x_e", max_knockouts=…)` → find growth-coupled knockout
   sets; rank by **guaranteed_product** (requires MILP + `straindesign`).
5. Optionally `knockout_comparison`/`batch_comparison` to inspect how a chosen design reroutes
   flux and whether it stays feasible.
6. Report: recommended targets, expected product & growth, whether coupling is guaranteed, and
   any solver limitation that constrained the analysis.

---

**When unsure which analysis fits, re-read §2 and §4, and state the assumption you are making
(reference choice, aerobic/anaerobic, solver) so the user can correct it.**
