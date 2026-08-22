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
  condition: explicit medium, substrate uptake, and oxygen/aeration bounds
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

When the goal *is* production, the canonical
[`SC-01`](SC-01-production-target-discovery.md) workflow already runs a matched MOMA-L2 and
ROOM single-knockout stage and validates its shortlist. Run SC-03 as well only when the user
wants an exhaustive essentiality/classification study; hand its full table to SC-01 as context,
not as a replacement for OptKnock/RobustKnock or forward validation.

**Success criteria.** A complete table covering every perturbation including the lethal ones,
with an explicit essentiality threshold, and a stated reference and method.

## Pipeline at a glance

| Step | Question | Method | Output |
|---|---|---|---|
| 0 | Usable model? | `_preflight.md` | viability, solver |
| 1 | Compared against what? | `reference_flux` | baseline flux state |
| 2 | What does each deletion do? | `batch_comparison` | full screen table |
| 3 | Which class is each? | thresholds on growth and product | classified table |
| 4 | Are the interesting ones real? | `flux_response`, paired wild-type/knockout sampling | verified subset |
| 5 | Write it up | `_reporting.md` | report + figures + raw data |

---

## Step 0 — Preflight

**Call.** First confirm one condition — medium, substrate uptake, oxygen/aeration bounds and
other changed bounds — then run `_preflight.md` P1–P3; P4 only if a product is named; P5 only
if using an omics baseline.

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
import json
from cmm.features import batch_comparison, gene_perturbations, reaction_perturbations

perturbations = gene_perturbations(model)        # or reaction_perturbations(model)
screen = batch_comparison(
    model, reference, perturbations,
    method="moma_l2",                            # moma_l1 (LP) / room (MILP)
    product_reaction=PRODUCT,                    # omit if no product
)
table = screen.to_frame()
provenance = screen.metadata                     # the whole screen's run record
```

**Save the provenance with the table.** `screen.metadata` is a `run_provenance` block — model
fingerprint, UTC timestamp, solver and solver version, platform, package versions, every
parameter, the reference state's identity, and (for `room`) the tolerance pair — plus what the
enumeration covered: `n_perturbations`, `n_inert_dropped` and `n_candidates_considered`.
`gene_perturbations` omits genes whose deletion blocks no reaction (**66 of 137 on
`e_coli_core`**), so without those counts the screen silently understates its own coverage.
Write it next to the CSV: `json.dumps(screen.metadata, indent=2, default=str)`.

**Outputs.** `target_id | kind | status | objective_value | distance | distance_kind |
n_changed_reactions | objective | n_reactions | product_flux` — ten columns as of 0.4.0, which
split the one overloaded `distance` field.

**Name each column for the quantity it holds.** `objective_value` is the raw solver objective
and means something different per method: `Σd²` for `moma_l2`, `Σ|d|` for `moma_l1`, a **count
of switched reactions** for `room`. `distance` is a distance and only a distance — Segrè et al.
Eq. (4)'s Euclidean `√(Σd²)` for `moma_l2`, the L1 sum for `moma_l1`, and `None` for `room`,
whose count lives on `n_changed_reactions`. `distance_kind` records which per row. A screen
exported before 0.4.0 wrote the objective in a column called a distance; the two differ by a
factor of about 36 for `moma_l2`, so never compare one against the other.

**`method="room"` selects a tolerance pair.** `batch_comparison` defaults to
`room_use_case="lethality"` (δ=0.1, ε=0.01), which is Shlomi et al.'s pair for exactly this
question; `room_use_case="flux_prediction"` (δ=0.03, ε=0.001) is the other published pair and
gives about 24% more switches. State which one the screen used.

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
from cmm.features._perturbation import blocked_reactions_for_genes

for target in beneficial_targets:
    blocked = blocked_reactions_for_genes(model, [target])
    with model:
        for rid in blocked:
            model.reactions.get_by_id(rid).knock_out()
        after = flux_response(model, PRODUCT, biomass_fraction=0.3, n_steps=15)

wild_type_ensemble = random_flux_sampling(model, n=1000, seed=0)

knockout_ensembles = {}
for target in beneficial_targets:
    blocked = blocked_reactions_for_genes(model, [target])
    with model:
        for rid in blocked:
            model.reactions.get_by_id(rid).knock_out()
        knockout_ensembles[target] = random_flux_sampling(model, n=1000, seed=0)
```

**Decision rule.**
- The knocked-out model's product response should exceed the wild-type value across a usable
  growth range, not only at a single point.
- Compare matched wild-type and knockout ensembles using the same condition, objective
  conditioning, sampler, seed policy, count, thinning, and reaction set. Report medians and
  intervals for product, biomass, blocked reactions and mechanistically relevant reactions.
  Correlated sampler draws are not biological replicates, so do not headline an unqualified
  p-value.
- A reaction whose wild-type flux was already near zero across the ensemble cannot have been
  doing much; a large predicted deletion benefit then deserves scepticism.

**Artifacts.** `05_verification/flux_response_<target>.csv`,
`05_verification/sampling_index.csv`, and target-specific wild-type/knockout sampling summaries.

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
  cannot find synthetic-lethal pairs or the multi-knockout designs in SC-01's strain-design
  step.

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
  SC-01's OptKnock/RobustKnock `guaranteed_product`, not by a screen.
