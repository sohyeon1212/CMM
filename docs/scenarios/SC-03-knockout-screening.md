---
id: SC-03
title: Single-deletion capacity and adaptation screening
goal: Screen every gene deletion for growth capacity, then optionally predict adaptation and product effects
when_to_use:
  - "which genes are essential"
  - "screen all knockouts / single deletion study"
  - "screen knockouts and identify essential genes"
role: standalone            # a complete study; not a sub-step of SC-01
optional_next:
  - "SC-01: when the goal is production, hand it the beneficial-deletion candidates"
requires:
  model: cobra model
  biomass: explicit biomass reaction id used as the maximum-growth objective
  condition: explicit medium, substrate uptake, and oxygen/aeration bounds
optional:
  product: exchange reaction id (adds the production column)
solver:
  minimum: LP (FBA growth capacity)
  full: QP (MOMA-L2), MILP (ROOM)
steps: [preflight, baseline, screen, classify, verify, report]
runtime: minutes on a core model; hours genome-scale
---

# SC-03 — Single-deletion capacity and adaptation screening

## Objective

Run every single gene deletion, calculate its maximum feasible growth under the declared
condition, and classify it as essential, impairing, dispensable, or model-uninformative. When a
product is named, add separate MOMA/ROOM adaptation and product-phenotype tables.

Essentiality is defined from FBA maximum growth capacity after deletion. MOMA and ROOM minimize
adjustment or regulatory switches relative to a reference; their predicted biomass is an
adaptation phenotype, not the deletion model's maximum capacity, and must not define the
essentiality class.

**This is a complete study on its own.** An essentiality map — which deletions are lethal,
which merely impair growth — is a finished answer, and the usual reason to run it (checking a
new reconstruction, planning a deletion library) has nothing to do with production.

When the goal *is* production, the canonical
[`SC-01`](SC-01-production-target-discovery.md) workflow already runs a matched MOMA-L2 and
ROOM single-knockout stage and validates its shortlist. Run SC-03 as well only when the user
wants an exhaustive essentiality/classification study; hand its full table to SC-01 as context,
not as a replacement for OptKnock/RobustKnock or forward validation.

**Success criteria.** A complete capacity table covering every model gene, including inert and
lethal rows, with explicit thresholds and condition. Optional MOMA/ROOM outputs remain separate
and method-labelled.

## Pipeline at a glance

| Step | Question | Method | Output |
|---|---|---|---|
| 0 | Usable model? | `_preflight.md` | viability, solver |
| 1 | What can wild type achieve? | FBA capacity plus pFBA reference | baseline growth and flux state |
| 2 | What growth capacity remains after each deletion? | GPR-resolved deletion plus FBA | complete capacity table |
| 3 | Which dependency class is each? | thresholds on FBA growth ratio | classified table |
| 4 | How do viable deletions adapt, and what happens to product? | separate MOMA-L2/ROOM and optional flux response | method-specific phenotype tables |
| 5 | Write it up | `_reporting.md` | report + figures + raw data |

---

## Step 0 — Preflight

**Call.** Run `_preflight.md` P0–P3 to resolve and confirm the biomass reaction plus one
condition—medium, substrate uptake, oxygen/aeration bounds and other changed bounds. P4 applies
only if a product is named; this recipe does not define an omics-reference branch. Do not begin
condition-dependent screening until an interview-created run definition is confirmed.

**Branch.** The exhaustive capacity layer remains LP FBA at any model size. For the optional
adaptation layer, declare its gene coverage and MOMA/ROOM methods before the run. A user-approved
MOMA-L1 subset is a narrower adaptation study, not a silent replacement for MOMA-L2.

---

## Step 1 — Conditioned wild-type baseline

**Goal.** Establish wild-type maximum growth for capacity ratios and a deterministic pFBA state
for the optional MOMA/ROOM layer.

**Call.**

```python
from dataclasses import replace

from cmm.core import ObjectiveSpec, apply_medium, fba
from cmm.features import reference_flux

# MEDIUM, CONDITION, and BIOMASS are the complete run definition resolved in preflight.
conditioned_model = model.copy()
medium_application = apply_medium(conditioned_model, MEDIUM)
capacity_condition = replace(
    CONDITION,
    objective=ObjectiveSpec(coefficients={BIOMASS: 1.0}, direction="max"),
)
capacity_condition.apply_to(conditioned_model)

wild_type_capacity = fba(conditioned_model)
if wild_type_capacity.status != "optimal" or wild_type_capacity.objective_value <= 1e-6:
    raise RuntimeError("wild type does not grow in the declared condition")
wild_type_growth = float(wild_type_capacity.objective_value)
reference = reference_flux(conditioned_model, "pfba")
```

**Decision rule.** FBA maximum biomass is the denominator for essentiality. pFBA supplies a
reproducible reference only for the optional adaptation layer; its particular flux vector does
not change the FBA capacity class. An LAD/E-Flux2 reference is a different, separately declared
study and is not substituted automatically.

**Artifacts.** `02_baseline/wild_type_capacity.csv`, `02_baseline/reference_pfba.csv`, plus the
complete applied-medium record and `capacity_condition`, including the explicit max-biomass
objective.

**Solver.** LP.

---

## Step 2 — Exhaustive single-deletion growth capacity

**Goal.** One maximum-biomass FBA solve per model gene, including genes whose GPR deletion
blocks no reaction.

**Call.**

```python
import pandas as pd

from cmm.core import fba
from cmm.features import gene_perturbations

perturbations = gene_perturbations(conditioned_model, include_inert=True)
capacity_results = {}
rows = []
for perturbation in perturbations:
    with conditioned_model:
        for reaction_id in perturbation.reaction_ids:
            conditioned_model.reactions.get_by_id(reaction_id).knock_out()
        result = fba(conditioned_model)

    capacity_results[perturbation.target_id] = result
    growth = (
        float(result.objective_value)
        if result.status == "optimal" and result.objective_value is not None
        else float("nan")
    )
    rows.append(
        {
            "target_id": perturbation.target_id,
            "blocked_reactions": ";".join(perturbation.reaction_ids),
            "n_reactions": len(perturbation.reaction_ids),
            "status": result.status,
            "growth_capacity": growth,
            "growth_ratio": growth / wild_type_growth,
        }
    )

table = pd.DataFrame(rows)
```

**Preserve typed results and provenance.** `capacity_results` retains each `FluxSolution`, so
the complete flux vector, solve summary, and metadata remain exportable. The combined table is
an index derived from those typed results, not a replacement for them. Record that
`include_inert=True` was used and verify that the row count equals `len(conditioned_model.genes)`.

**Artifacts.** `03_screen/single_deletion_capacity.csv`, target-specific FBA result exports and
metadata, and a coverage record containing total model genes, attempted genes, and status
counts. Lethal and non-optimal rows remain in the main table.

**Decision rule.** Keep every row. Infeasible deletions and optimal deletions whose maximum
growth falls below the declared essentiality threshold establish the predicted essential set;
filtering either out destroys the primary finding.

**Branch.** `n_reactions == 0` means the gene's deletion blocks no reaction under the GPR. Mark
it `uninformative` rather than reporting a safe or neutral knockout—the model has no deletion
phenotype for that gene.

**Failure → action.** A partial gene set is allowed only when explicitly declared. Never call a
subset exhaustive, and never omit non-optimal rows from the coverage count.

**Solver.** LP.

---

## Step 3 — Classify

**Goal.** Turn FBA growth capacity into declared dependency classes without mixing in an
adaptation method or product score.

**Decision rule** (state your thresholds in the report; these are defaults, not standards):

| Class | Condition |
|---|---|
| essential | `status != "optimal"` or `growth_ratio < 0.01` |
| impairing | `0.01 <= growth_ratio < 0.9` |
| dispensable | `growth_ratio >= 0.9` |
| uninformative | `n_reactions == 0` |

Apply the capacity classes first, then let `uninformative` override the apparent reference-like
growth of a gene whose deletion changes no model reaction:

```python
ratio = table["growth_ratio"]

table["class"] = "dispensable"
table.loc[ratio < 0.9, "class"] = "impairing"
table.loc[(table["status"] != "optimal") | (ratio < 0.01), "class"] = "essential"
table.loc[table["n_reactions"] == 0, "class"] = "uninformative"
```

The 1% and 90% cutoffs are declared defaults, not universal biological standards. Preserve the
continuous growth ratio beside the class and record any changed thresholds in provenance.

**Artifacts.** `04_classified/knockout_classes.csv`, `04_classified/essential_genes.csv`,
`figures/knockout_impact.png` (an R scatter or ordered bar chart built from `target_id`,
`growth_ratio`, and `class` in the capacity table).

**Branch.** No product → the capacity study may stop after this step's report. With a product →
step 4 adds method-specific product phenotypes; it does not change the capacity class.

---

## Step 4 — Optional adaptation and product phenotype

**Goal.** Predict how viable deletions redistribute flux relative to pFBA, without redefining
their FBA capacity class. When a product is named, compare product across a usable growth range.

The example below assumes `PRODUCT` is named. For an adaptation-only study, omit
`product_reaction` and the product-response block while retaining the separate MOMA/ROOM tables.

**Call.**

```python
from math import isfinite

from cmm.features import (
    batch_comparison,
    blocked_reactions_for_genes,
    flux_response,
)

# Fix these method-independent decision thresholds before inspecting either ranking.
MIN_GROWTH_RETENTION = 0.1
PRODUCT_IMPROVEMENT_TOLERANCE = 1e-6

moma_screen = batch_comparison(
    conditioned_model,
    reference,
    perturbations,
    method="moma_l2",
    product_reaction=PRODUCT,
)
room_screen = batch_comparison(
    conditioned_model,
    reference,
    perturbations,
    method="room",
    room_use_case="flux_prediction",
    product_reaction=PRODUCT,
)

wild_type_product = reference.get(PRODUCT)
viability_floor = MIN_GROWTH_RETENTION * wild_type_growth


def product_improving_targets(screen):
    return {
        row.target_id
        for row in screen
        if row.status == "optimal"
        and isfinite(row.objective)
        and row.objective >= viability_floor
        and isfinite(row.product_flux)
        and row.product_flux > wild_type_product + PRODUCT_IMPROVEMENT_TOLERANCE
    }


moma_candidates = product_improving_targets(moma_screen)
room_candidates = product_improving_targets(room_screen)
PRODUCT_CANDIDATES = tuple(sorted(moma_candidates | room_candidates))

# Retain each method's proposing set beside the union; never merge their unlike scores.
wild_type_response = flux_response(
    conditioned_model,
    BIOMASS,
    response=PRODUCT,
    biomass=BIOMASS,
    target_min=0.3 * wild_type_growth,
    target_max=wild_type_growth,
    n_steps=20,
)
knockout_responses = {}
for target in PRODUCT_CANDIDATES:
    blocked = blocked_reactions_for_genes(conditioned_model, [target])
    with conditioned_model:
        for reaction_id in blocked:
            conditioned_model.reactions.get_by_id(reaction_id).knock_out()
        knockout_responses[target] = flux_response(
            conditioned_model,
            BIOMASS,
            response=PRODUCT,
            biomass=BIOMASS,
            target_min=0.3 * wild_type_growth,
            target_max=wild_type_growth,
            n_steps=20,
        )
```

**Decision rule.**
- Keep MOMA and ROOM rows in separate tables; their objective and distance fields are not
  interchangeable. Neither method defines essentiality.
- Define product-improving rows separately within each method using the same wild-type product
  value, growth floor, tolerance, and tie policy. Preserve lethal and failed rows.
- Plot growth on the x-axis and target-product flux on the y-axis on the same explicit WT-based
  grid. Retain infeasible KO points and statuses when a deletion cannot reach part of that
  range. A single favorable point is not evidence across the usable range.
- Plain random sampling does not reproduce MOMA/ROOM conditioning and is not required for this
  capacity recipe. If added as a sensitivity analysis, declare its sampled constraints and do
  not treat correlated draws as biological replicates.

**Artifacts.** `05_adaptation/single_deletion_moma.csv`,
`05_adaptation/single_deletion_room.csv`, `05_adaptation/wild_type_product_response.csv`, and
target-specific `05_adaptation/product_response_<target>.csv` files.

**Solver.** **QP** for MOMA-L2, **MILP** for ROOM, LP for flux response. A capacity-only study
stops after step 3 and needs LP only.

---

## Step 5 — Report

Follow the reproducibility principles in `_reporting.md`. SC-03 is a public-service recipe, so
its caller must define and validate its own artifact schema; `validate_production_run` does not
validate this study. Scenario-specific requirements:

- State that essentiality uses FBA maximum biomass, name the pFBA reference used only for
  adaptation, list any MOMA/ROOM methods, and report every classification threshold. The
  essential-gene list is meaningless without the condition and capacity rule.
- Report class counts and the full table location; list essential genes explicitly.
- Mark the `n_reactions == 0` genes as uninformative rather than neutral.
- Keep these caveats beside the affected results: essentiality is *in silico* under one medium
  and changes with it; MOMA/ROOM predict adaptation rather than capacity; and a single-deletion
  recipe cannot find synthetic-lethal pairs or SC-01's multi-knockout designs.

## Cross-checks

- Class counts sum to the number of perturbations screened.
- Any gene the literature calls essential that this screen calls dispensable is worth flagging —
  usually a medium difference or a gap in the reconstruction.
- Wild-type growth is the same value everywhere it appears.

## Do not

- Do not drop lethal or infeasible rows from the exported table.
- Do not report essentiality without the medium and threshold.
- Do not define essentiality from MOMA or ROOM biomass; use the FBA maximum-capacity table.
- Do not silently screen a subset; a partial screen states its coverage.
- Do not read a product-improving single deletion as growth-coupled. Model-certified *in silico*
  coupling comes from SC-01's OptKnock/RobustKnock `guaranteed_product`, not this screen.
