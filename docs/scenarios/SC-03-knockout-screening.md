---
id: SC-03
title: Knockout screening and essentiality
goal: Screen every single knockout for lethality, growth impact, and effect on a product
when_to_use:
  - "which genes are essential"
  - "screen all knockouts / single deletion study"
  - "넉아웃 스크리닝, 필수 유전자"
role: standalone            # a complete study; not a sub-step of SC-01
optional_next:
  - "SC-01: when the goal is production, hand it the beneficial-deletion candidates"
requires:
  model: cobra model
optional:
  product: exchange reaction id (adds the production column)
  expression: enables LAD / E-Flux2 baselines
solver:
  minimum: LP (moma_l1)
  full: QP (moma_l2), MILP (room)
steps: [preflight, baseline, screen, classify, verify, report]
runtime: minutes on a core model; hours genome-scale
---

# SC-03 — Knockout screening and essentiality

## Objective

Run every single gene (or reaction) deletion against a wild-type reference and classify each
as essential, impairing, neutral, or beneficial — with its effect on a target product when one
is named.

**This is a complete study on its own.** An essentiality map — which deletions are lethal,
which merely impair growth — is a finished answer, and the usual reason to run it (checking a
new reconstruction, planning a deletion library) has nothing to do with production.

It happens to be `SC-01` step 3 run exhaustively, so when the goal *is* production the two
compose: hand the beneficial-deletion candidates to [`SC-01`](SC-01-production-target-discovery.md)
for verification. That is an option, not an obligation — neither scenario requires the other.

**Success criteria.** A complete table covering every perturbation including the lethal ones,
with an explicit essentiality threshold, and a stated reference and method.

## Pipeline at a glance

| Step | Question | Method | Output |
|---|---|---|---|
| 0 | Usable model? | `_preflight.md` | viability, solver |
| 1 | Compared against what? | `reference_flux` | baseline flux state |
| 2 | What does each deletion do? | `batch_comparison` | full screen table |
| 3 | Which class is each? | thresholds on growth and product | classified table |
| 4 | Are the interesting ones real? | `flux_response`, sampling | verified subset |
| 5 | Write it up | `_reporting.md` | report + figures + raw data |

---

## Step 0 — Preflight

**Call.** `_preflight.md` P1–P3; P4 only if a product is named; P5 only if using an omics
baseline.

**Branch.** Genome-scale model → decide the method now, not after a failed overnight run: a
`moma_l1` (LP) screen is the practical default, with `moma_l2` re-checks on the survivors.

---

## Step 1 — Wild-type reference

**Goal.** Every distance and every growth comparison is relative to this state.

**Call.**

```python
from cmm.features import reference_flux

reference = reference_flux(model, "pfba")
wild_type_growth = reference.get(BIOMASS)
```

**Decision rule.** `pfba` is the reproducible default (unique minimal-total-flux solution). Use
`lad`/`eflux2` when the question is condition-specific, and run the screen once per reference
if you need to compare — never mix references within one table.

**Artifacts.** `02_baseline/reference_pfba.csv`.

**Solver.** LP (pfba/lad) or **QP** (eflux2).

---

## Step 2 — The screen

**Goal.** One solve per deletion.

**Call.**

```python
import pandas as pd
from cmm.features import batch_comparison, gene_perturbations, reaction_perturbations

perturbations = gene_perturbations(model)        # or reaction_perturbations(model)
rows = batch_comparison(
    model, reference, perturbations,
    method="moma_l2",                            # moma_l1 (LP) / room (MILP)
    product_reaction=PRODUCT,                    # omit if no product
)
table = pd.DataFrame([vars(r) for r in rows])
```

**Outputs.** `target_id | kind | status | distance | objective | n_reactions | product_flux`.

**Artifacts.** `03_screen/batch_<reference>_<method>.csv` — **the complete table, lethal rows
included**.

**Decision rule.** Keep every row. The infeasible ones are the essentiality result; filtering
them out at export time destroys the primary finding of this scenario.

**Branch.** `n_reactions == 0` means the gene's deletion blocks no reaction under the GPR, so
the result is trivially the reference. Mark these rather than reporting them as "no effect
knockouts" — the model simply has no information about them.

**Failure → action.** A genome-scale `moma_l2` screen may take hours. Reduce by switching to
`moma_l1`, not by sampling a subset of genes; a partial screen must say which genes it covered.

**Solver.** `moma_l1` LP, `moma_l2` **QP**, `room` **MILP**.

---

## Step 3 — Classify

**Goal.** Turn the raw table into the four classes people actually ask for.

**Decision rule** (state your thresholds in the report; these are defaults, not standards):

| Class | Condition |
|---|---|
| essential | `status != "optimal"` or `objective < 0.01 * wild_type_growth` |
| impairing | `0.01 <= objective / wild_type_growth < 0.9` |
| neutral | `objective >= 0.9 * wild_type_growth` and no product change |
| beneficial | `product_flux > wild_type_product_flux` and growth above your viability floor |

Apply them in this order — later rules overwrite earlier ones, and `essential` must win over
everything, including a deletion that looks beneficial for the product but does not grow:

```python
ratio = table["objective"] / wild_type_growth

table["class"] = "neutral"
table.loc[ratio < 0.9, "class"] = "impairing"

# Only meaningful when the screen was run with product_reaction=; skip this line otherwise.
beneficial = (table["product_flux"] > wild_type_product) & (ratio >= 0.1)
table.loc[beneficial, "class"] = "beneficial"

table.loc[(table["status"] != "optimal") | (ratio < 0.01), "class"] = "essential"
```

`0.1` is the viability floor from step 2's decision rule; keep it the same number in both
places and state it in the report.

**Artifacts.** `04_classified/knockout_classes.csv`, `04_classified/essential_genes.csv`,
`figures/knockout_impact.png` (e.g. `flux_comparison_figure` on the top-impact deletions).

**Branch.** No product → the screen is an essentiality study; stop after this step's report.
With a product → the `beneficial` class is the input to step 4 and to `SC-01`.

---

## Step 4 — Verify the interesting deletions

**Goal.** Check that beneficial deletions are real rather than alternate-optima artifacts.

**Call.**

```python
from cmm.features import flux_response, random_flux_sampling

for target in beneficial_targets:
    blocked = blocked_reactions_for_genes(model, [target])
    with model:
        for rid in blocked:
            model.reactions.get_by_id(rid).knock_out()
        after = flux_response(model, PRODUCT, biomass_fraction=0.3, n_steps=15)

ensemble = random_flux_sampling(model, n=1000, seed=0)
```

**Decision rule.**
- The knocked-out model's product response should exceed the wild-type value across a usable
  growth range, not only at a single point.
- Check the deleted reaction's sampled distribution: a reaction whose wild-type flux was
  already near zero across the ensemble cannot have been doing much, so a large predicted
  benefit from deleting it deserves scepticism.

**Artifacts.** `05_verification/flux_response_<target>.csv`,
`05_verification/sampling_statistics.csv`.

**Solver.** LP.

---

## Step 5 — Report

Follow `_reporting.md`. Scenario-specific requirements:

- State the reference method, the comparison method, and every threshold used for
  classification. The essential-gene list is meaningless without them.
- Report class counts and the full table location; list essential genes explicitly.
- Mark the `n_reactions == 0` genes as uninformative rather than neutral.
- **Limitations** must include: essentiality here is *in silico* under one medium and changes
  with it; MOMA/ROOM assume minimal adjustment from wild type; and single-deletion screens
  cannot find synthetic-lethal pairs or the multi-knockout designs of `SC-01` step 3.

## Cross-checks

- Class counts sum to the number of perturbations screened.
- Any gene the literature calls essential that this screen calls neutral is worth flagging —
  usually a medium difference or a gap in the reconstruction.
- Wild-type growth is the same value everywhere it appears.

## Do not

- Do not drop lethal or infeasible rows from the exported table.
- Do not report essentiality without the medium and threshold.
- Do not silently screen a subset; a partial screen states its coverage.
- Do not read a beneficial single deletion as growth-coupled — coupling is proven by
  `SC-01` step 3's `guaranteed_product`, not by a screen.
