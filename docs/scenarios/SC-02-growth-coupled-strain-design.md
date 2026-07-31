---
id: SC-02
title: Growth-coupled strain design
goal: Find knockout sets that force the cell to produce whenever it grows
when_to_use:
  - "design a production strain"
  - "make production obligatory / growth-coupled"
  - "성장 공역 균주 설계"
role: standalone            # its own goal (coupling); runs without any other scenario
commonly_follows:
  - "SC-01: screen candidates cheaply there, then design rigorously here"
requires:
  model: cobra model
  product: exchange reaction id
optional:
  max_knockouts: design size (default 3)
  min_growth: viability floor (default 0.05)
solver:
  minimum: MILP + the straindesign package + Java/OpenJDK
  full: MILP, QP (for the MOMA verification in step 4)
steps: [preflight, feasibility, design, verify, report]
runtime: minutes on a core model; can be hours on genome scale
---

# SC-02 — Growth-coupled strain design

## Objective

Find a small set of reaction knockouts such that a cell maximizing its own growth **cannot
avoid** making the product, and quantify the guaranteed production rate.

**Coupling is its own goal, and this scenario is complete on its own** — a model, a product,
and a MILP solver are all it needs. `SC-01` is not a prerequisite.

Its natural place is nonetheless *after* `SC-01`: SC-01 ranks candidates cheaply, SC-02 proves
coupling with a bilevel MILP, so screening first narrows what this search has to enumerate.
Run SC-01 first when the user has no candidate targets yet and the MILP search would otherwise
be unconstrained; start here directly when the ask is specifically for a guaranteed-production
strain.

**Success criteria.** At least one design with `guaranteed_product > 0` and growth above the
user's viability floor, verified independently, with its knockouts mapped back to genes.

## Pipeline at a glance

| Step | Question | Method | Output |
|---|---|---|---|
| 0 | Usable model, capable solver? | `_preflight.md` | viability, MILP + straindesign |
| 1 | Can coupling exist here at all? | `production_envelope` | trade-off shape |
| 2 | Which knockout sets couple? | `optknock`, `robustknock` | designs |
| 3 | Does the design hold up? | `knockout_comparison`, `flux_response` | verified design |
| 4 | Write it up | `_reporting.md` | report + figures + raw data |

---

## Step 0 — Preflight

**Goal.** Confirm viability, product reachability, **and** that this scenario can run at all.

**Call.** `_preflight.md` P1–P4, plus:

```python
from cmm.core import supports
print(supports("MILP", model.solver.interface))
import importlib.util; print(importlib.util.find_spec("straindesign") is not None)
```

**Branch.** No MILP solver, no `straindesign` package, or no Java → **stop and tell the user**.
There is no LP substitute for OptKnock. Offer `SC-01` instead, stating clearly that it ranks
candidates but cannot guarantee coupling.

---

## Step 1 — Is coupling plausible?

**Goal.** Avoid spending MILP time on a product that cannot be coupled.

**Call.**

```python
from cmm.features import production_envelope
envelope = production_envelope(model, PRODUCT, aerobic=AEROBIC, points=20)
frame = envelope.to_frame()
```

**Decision rule.** Look at `growth_min` across the range. If it is 0 everywhere, the cell can
always grow without producing, so coupling must be *created* by knockouts — possible, but
harder. If `growth_max` falls steeply with product, growth and production compete, which is the
usual precondition for a good design.

**Artifacts.** `02_feasibility/production_envelope.csv`, `figures/production_envelope.png`.

**Solver.** LP.

---

## Step 2 — Search for designs

**Goal.** Enumerate knockout sets that couple production to growth.

**Call.**

```python
from cmm.features import optknock, robustknock

optimistic = optknock(model, PRODUCT, max_knockouts=3, max_solutions=5, min_growth=0.05)
guaranteed = robustknock(model, PRODUCT, max_knockouts=3, max_solutions=8, min_growth=0.05)
```

**Outputs.** `StrainDesignResult` with `.designs` and `.best()`; each `StrainDesign` has
`.knockouts`, `.growth`, `.max_product`, `.guaranteed_product`, `.growth_coupled`.

**Artifacts.** `03_design/optknock.csv`, `03_design/robustknock.csv` — one row per design with
all of the above.

**Decision rule.** **Rank by `guaranteed_product`, never `max_product`.** `max_product` is what
the cell *could* make if it cooperated; `guaranteed_product` is what it makes at worst among
growth-optimal states. Only `guaranteed_product > 0` means growth-coupled. A design with high
`max_product` and zero `guaranteed_product` is not a design — it is a possibility.

`robustknock` optimizes the worst case directly, so prefer its designs when both return
results.

**Branch.**
- No design found → raise `max_knockouts` (cost grows fast), lower `min_growth`, or switch the
  medium/aeration. If still nothing, report that no coupled design exists under these
  constraints — a legitimate negative result.
- Designs found but all with `guaranteed_product == 0` → report them as uncoupled candidates
  and move to `SC-01` verification rather than presenting them as strains.

**Failure → action.** MILP timeouts on genome-scale models are common. Reduce `max_solutions`
first, then `max_knockouts`, and record what you reduced. Do not silently switch models.

**Solver.** **MILP** + `straindesign` + Java.

---

## Step 3 — Verify the design

**Goal.** Confirm the design behaves as claimed and understand how it reroutes flux.

**Call.**

```python
from cmm.features import flux_response, knockout_comparison, reference_flux
from cmm.visualization import flux_comparison_figure, flux_response_figure

best = guaranteed.best()
reference = reference_flux(model, "pfba")

response = knockout_comparison(model, reference, best.knockouts, method="moma_l2")

with model:                                   # the design applied, then restored
    for rid in best.knockouts:
        model.reactions.get_by_id(rid).knock_out()
    scan = flux_response(model, PRODUCT, biomass_fraction=0.0, n_steps=20)
```

**Decision rule.**
- `response.status == "optimal"` — the design is not lethal under a minimal-adjustment cell.
- Compare `response.fluxes[PRODUCT]` with the wild-type product flux: MOMA asks what the cell
  does *immediately* after the deletion, before adaptive evolution. A design that couples at
  growth optimum but produces little under MOMA typically needs adaptation to reach its
  designed phenotype — worth reporting.
- The `flux_response` scan on the knocked-out model shows the growth/product relationship the
  design creates; growth at zero product should be low or infeasible if coupling is real.

**Artifacts.** `04_verification/knockout_comparison.csv`,
`04_verification/design_flux_response.csv`, `figures/design_flux_comparison.png`,
`figures/design_flux_response.png`.

**Branch.** Map reaction knockouts back to genes for anything experimental — a reaction is not
directly deletable. `model.reactions.get_by_id(rid).genes` and the reaction's GPR tell you
which deletions realize it, and whether the GPR is an `or` (needing several deletions).

**Solver.** **QP** for `moma_l2` (use `moma_l1` on LP-only solvers, and say so); LP for the scan.

---

## Step 4 — Report

Follow `_reporting.md`. Scenario-specific requirements:

- The designs table must show `knockouts | growth | max_product | guaranteed_product |
  growth_coupled`, ranked by `guaranteed_product`.
- State explicitly for each recommended design whether it is growth-coupled, and never quote
  `max_product` without `guaranteed_product` beside it.
- Give the gene-level deletions implied by each reaction knockout, including where an `or` GPR
  requires deleting several genes.
- **Limitations** must include: coupling is a property of the model's growth-maximizing
  assumption; MOMA predicts the immediate post-deletion phenotype, which may differ from the
  designed one until the strain adapts; MILP results depend on solver, tolerances, and any
  reduction in `max_knockouts`/`max_solutions` you made for runtime.

## Cross-checks

- `guaranteed_product <= max_product` for every design; a violation means a reporting error.
- No design's product exceeds the theoretical yield from `_preflight.md` P4.
- The knockouts are all real reaction ids present in the model.

## Do not

- Do not rank or recommend by `max_product`.
- Do not call a design growth-coupled without `guaranteed_product > 0`.
- Do not present reaction knockouts as if they were gene deletions.
- Do not mix these designs with SC-01's MOMA/ROOM candidates in one ranked list — they assume
  different cell behavior. Report them as separate evidence.
