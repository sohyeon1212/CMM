---
id: SC-02
title: Omics-driven context engineering
goal: Explain a metabolic difference between conditions or strains, and turn it into targets
when_to_use:
  - "why does condition A produce more than B"
  - "what changes between these two strains / timepoints"
  - "use omics data to compare conditions and identify targets"
role: standalone            # a complete study; not a sub-step of SC-01
optional_next:
  - "SC-01: when the goal is production, hand it the resolved biological condition"
requires:
  model: cobra model
  expression: gene x condition table with ids matching model.genes
optional:
  product: exchange reaction id (needed for the targeting half)
solver:
  minimum: LP (LAD)
  full: QP (E-Flux2)
steps: [preflight, integrate, compare, interpret, target, report]
runtime: minutes per condition
---

# SC-02 — Omics-driven context engineering

## Objective

Turn a gene expression table into condition-specific flux states, identify where those states
differ metabolically, and — when a product is named — use the better-performing condition as
the context in which to search for engineering targets.

**This is a complete study on its own.** "Why does condition A outproduce B" is a finished
answer, and the mechanism it produces is usually the point. When the goal goes further and
asks for targets, it composes with [`SC-01`](SC-01-production-target-discovery.md): SC-02
chooses *which condition* to search in, SC-01 does the systematic search there. That is an
option, not an obligation — neither scenario requires the other.

**Success criteria.** Each condition has a flux state with optimal status; the differences are
quantified per reaction rather than per gene; and any target proposed is tied to the specific
condition it was derived in.

## Pipeline at a glance

| Step | Question | Method | Output |
|---|---|---|---|
| 0 | Do the ids even map? | `_preflight.md` (P5 is critical) | overlap fraction |
| 1 | What flux does each condition imply? | `predict_condition_fluxes` | per-condition fluxes |
| 2 | Where do they differ? | `flux_log_change`, `sign_flips` | ranked differences |
| 3 | What does that mean? | pathway reading of the top changes | mechanism |
| 4 | What should we change? | `fseof` / `flux_response` in context | targets |
| 5 | Write it up | `_reporting.md` | report + figures + raw data |

---

## Step 0 — Preflight

**Goal.** Confirm the expression table actually maps onto the model. This scenario fails
silently more than any other.

**Call.** `_preflight.md` P0–P3 and **P5**, plus P4 if a product is named. P0 performs read-only
discovery and asks only unresolved scientific decisions; do not begin condition-dependent
integration until an interview-created run definition is confirmed.

**Decision rule.** Gene id overlap must cover most model genes. Below roughly half, stop and
ask for the mapping — an unmapped table yields a confident, meaningless flux state.

**Failure → action.** Report the overlap fraction and sample unmatched ids from both sides.
GEO tables typically carry symbols or probe ids; *E. coli* models use b-numbers.

---

## Step 1 — Predict per-condition fluxes

**Goal.** One flux state per expression column, under that state's fully declared external
condition. Expression alone does not encode medium or aeration.

**Call.**

```python
from dataclasses import asdict

from cmm.core import apply_medium, require
from cmm.omics import predict_condition_fluxes, read_expression_table

expression = read_expression_table("expression.csv")     # gene x condition
method = "eflux2"  # resolved in preflight; use "lad" only when explicitly selected
require("QP" if method == "eflux2" else "LP", model.solver.interface,
        feature=f"{method} expression integration")

# RUN_DEFINITIONS maps each expression column to its complete (Medium, Condition) pair from
# preflight. Medium controls available uptake; Condition carries objective and bound overrides.
integrated = {}
resolved_inputs = {}
expression_conditions = {str(name) for name in expression.columns}
defined_conditions = set(RUN_DEFINITIONS)
if defined_conditions != expression_conditions:
    missing = sorted(expression_conditions - defined_conditions)
    extra = sorted(defined_conditions - expression_conditions)
    raise ValueError(f"run-definition mismatch: missing={missing}, extra={extra}")

for name in map(str, expression.columns):
    medium, condition = RUN_DEFINITIONS[name]
    conditioned_model = model.copy()
    applied_medium = apply_medium(conditioned_model, medium)
    condition.apply_to(conditioned_model)
    one = predict_condition_fluxes(
        conditioned_model,
        expression,
        method=method,
        conditions=[name],
    )
    integrated[name] = one.results[name]
    resolved_inputs[name] = {
        "medium": applied_medium.to_provenance(),
        "condition": asdict(condition),
    }

fluxes = {
    name: result.fluxes
    for name, result in integrated.items()
    if result.status == "optimal"
}
```

**Outputs.** One `OmicsFluxResult` per independently conditioned model. If every expression
column genuinely shares the same external bounds, one `predict_condition_fluxes(model,
expression, method=method)` call may return a combined `ConditionFluxes`; its `.fluxes(name)`
raises when that column did not solve optimally. Do not use that convenience call for columns
with different media or oxygen bounds. Archive `resolved_inputs` beside the per-condition
result metadata so model fingerprints are accompanied by the biological definition they hash.

**Artifacts.** `02_integration/fluxes_<condition>.csv` per condition.

**Decision rule.** Every condition must reach `status == "optimal"`. E-Flux2 (QP) scales bounds
by normalized expression then minimizes total squared flux; LAD (LP) fits fluxes to
expression-derived targets. They answer slightly different questions — name which you used.

**Two 0.4.0 changes make every pre-0.4.0 result here stale.** The GPR `OR` rule for continuous
expression is now `"sum"` (Kim et al. 2016 and Lee et al. 2012 both specify it; CMM's previous
`max` matches no source paper) — measured on `e_coli_core` + GSE41189 it moves 30 of 66 mapped
reaction weights, by up to 2.67×, shifts the normalisation denominator, and changes predicted
growth by 26%. And LAD now fits the **absolute** flux, so a reversible reaction is no longer
penalised for running in reverse. `or_rule=` is a parameter on `integrate_expression` /
`gene_to_reaction_weights`; the rule used is recorded as `gpr_or_rule` in provenance, and
`metadata["cmm_deviations"]` lists where each method departs from its source. Quote both.

**E-Flux2 cannot constrain a reaction with no GPR** — `EX_o2_e` among them — so expression
alone does not switch off respiration. On the GSE41189 aerobic/anaerobic pair, expression
alone gets the fermentation half right (formate 4.3×, acetate up 17%, oxidative
phosphorylation down) and the TCA half wrong (up 13%), because the model still takes up 12.25
mmol gDW⁻¹ h⁻¹ of O₂ in the "anaerobic" sample. Applying the oxygen condition as well —
the corresponding state's complete `RUN_DEFINITIONS` entry — gives the expected physiology
(TCA down 33×, fermentation up 3.3×). This is a property of the method, not a tuning knob:
state which external condition you applied.

**Branch.** No QP → pause before changing the claim. Use LAD only when the run definition
explicitly selects that LP method and record the change; never choose it automatically from
solver availability. Never pass `allow_l1_fallback=True` and call the result E-Flux2.

**Failure → action.** A non-optimal condition means its expression-derived bounds are
infeasible. Report which condition and stop for it; do not substitute another method for one
condition and keep E-Flux2 for the rest — the comparison would be meaningless.

**Solver.** LP (LAD) or **QP** (E-Flux2).

---

## Step 2 — Quantify the difference

**Goal.** Rank reactions by how much their flux changed between two conditions.

**Call.**

```python
from cmm.omics import flux_log_change, sign_flips
from cmm.visualization import flux_log_change_figure, save_figure

source, target = "condition_A", "condition_B"
changes = flux_log_change(fluxes[source], fluxes[target])
flipped = sign_flips(fluxes[source], fluxes[target])
```

**Outputs.** `{reaction_id: log2 fold-change of |flux|}` and a list of reactions that reversed
direction.

**Artifacts.** `03_comparison/flux_log_change.csv`, `03_comparison/sign_flips.csv`,
`figures/flux_log_change.png`.

**Decision rule.** Rank by `abs(log2 change)`. Treat sign flips separately and prominently: a
reaction running backwards is a qualitative rerouting, not a magnitude change, and the log2
value understates it.

**Branch.** More than two conditions → run the comparison pairwise against a chosen reference
condition, and say which is the reference. Do not average flux states across conditions.

**Solver.** None (arithmetic on step 1's output).

---

## Step 3 — Interpret

**Goal.** Turn the ranked list into a mechanism a reader can follow.

**Call.**

```python
annotated = sum(1 for r in model.reactions if (r.subsystem or "").strip())
print(f"{annotated} of {len(model.reactions)} reactions carry a subsystem annotation")

for rid, change in sorted(changes.items(), key=lambda kv: -abs(kv[1]))[:20]:
    reaction = model.reactions.get_by_id(rid)
    print(rid, round(change, 2), reaction.subsystem or "—", reaction.gene_reaction_rule)
```

**Decision rule.** Group the top changes by `subsystem` — a coherent story is several reactions
in one pathway moving together, not one isolated reaction. Cross-check against the exchange
reactions: a change in uptake or secretion pattern usually explains a large internal shift.

**Check the annotation first.** Many models carry no subsystem field at all — `e_coli_core` has
it on **none** of its 95 reactions, while iJO1366 has it on about 87%. On an unannotated model
this grouping silently returns nothing, so fall back to grouping by shared metabolites instead,
and say in the report which you used:

```python
# Reactions that move together and share a metabolite are the same story.
top = [rid for rid, _ in sorted(changes.items(), key=lambda kv: -abs(kv[1]))[:20]]
for rid in top:
    shared = {
        other for other in top
        if other != rid
        and set(model.reactions.get_by_id(rid).metabolites)
        & set(model.reactions.get_by_id(other).metabolites)
    }
    print(rid, "shares metabolites with", sorted(shared) or "nothing else in the top 20")
```

**Branch.** With a named product, check the product exchange's own change first — that is the
result the user asked about.

**Failure → action.** If the top changes are scattered across unrelated subsystems with no
exchange changes, treat the interpretation as inconclusive and inspect condition definitions,
gene mapping, integration status, and method sensitivity. Plain `random_flux_sampling` samples
the base model's feasible space; it does not preserve LAD/E-Flux2 constraints and therefore
cannot validate an expression-derived flux difference.

---

## Step 4 — Targets in context

**Goal.** Only when a product is named: find interventions in the condition that matters.

**Call.**

```python
from dataclasses import asdict

from cmm.core import apply_medium
from cmm.features import flux_response, fseof

# Reuse the exact complete run definition selected in step 1. Do not reconstruct only its
# oxygen bound or fall back to the unconditioned source model.
selected_medium, selected_condition = RUN_DEFINITIONS[CHOSEN_CONDITION]
target_model = model.copy()
selected_medium_application = apply_medium(target_model, selected_medium)
selected_condition.apply_to(target_model)
selected_run_definition = {
    "medium": selected_medium_application.to_provenance(),
    "condition": asdict(selected_condition),
}

scan = fseof(target_model, PRODUCT, condition=selected_condition, n_steps=10)
verified = {}
for candidate in scan.amplification_targets()[:10]:
    verified[candidate] = flux_response(
        target_model,
        candidate,
        response=PRODUCT,
        condition=selected_condition,
        biomass_fraction=0.3,
        n_steps=20,
    )
```

Archive `selected_run_definition` with the target tables. Applying the full medium before the
condition preserves every required nutrient and uptake closure from the selected state; passing
the same condition to each service preserves its objective and explicit bound overrides.

**Decision rule.** Cross the FSEOF targets with step 2's ranked changes. A reaction that is
both an amplification target **and** already up in the better-producing condition is a
context-consistent intersection worth inspecting. The two signals are not independent—both
depend on the same model, and one is expression-constrained—so this intersection is neither
standalone validation nor an automatic recommendation. Report each `verified` response
separately and label a candidate forward-supported only when target-product flux improves over
its feasible scan under the declared biomass floor; retain flat, adverse, infeasible, and failed
responses with their status.

**Branch.** For a full target search, hand the resolved biological condition to `SC-01`; the
canonical production workflow then builds its own supported FBA or pFBA reference. Passing an
LAD/E-Flux2 flux state directly into `ProductionWorkflowConfig` is not currently a shipped
handoff. Adding that boundary requires the typed config, provenance, schema, and validation work
described in the [canonical-workflow tutorial](../tutorials/adding-a-canonical-workflow.md).

**Artifacts.** `04_targets/fseof_trends.csv`, `04_targets/flux_response_<target>.csv`,
`04_targets/context_consensus.csv`.

**Solver.** LP.

---

## Step 5 — Report

Follow the reproducibility principles in `_reporting.md`. SC-02 is a public-service recipe, so
its caller must define and validate its own artifact schema; `validate_production_run` does not
validate this study. Scenario-specific requirements:

- State the integration method (E-Flux2 or LAD) and the gene id overlap fraction up front.
  Every number depends on both.
- Report reaction-level changes, not gene-level — the flux state is the model's answer, the
  expression is only its input.
- Give sign flips their own subsection.
- Tie every proposed target to the condition it was derived in.
- **Limitations** must include: expression-derived fluxes are predictions constrained by
  expression, not measurements; E-Flux2's variant here is CMM's documented implementation
  (`docs/VALIDATION.md`); and unmapped genes contribute nothing, so the overlap fraction bounds
  the whole analysis.

## Cross-checks

- Each condition's growth rate is plausible and reported.
- The product's own flux change is consistent with the direction of the user's observation — if
  the model says B produces less while the user measured more, the integration or the medium is
  wrong, and that must be resolved before anything else in the report is trustworthy.

## Do not

- Do not proceed on a near-zero gene id overlap.
- Do not compare a LAD state against an E-Flux2 state — the difference would be method, not
  biology.
- Do not report gene log2 changes as if they were flux changes.
- Do not average flux states across conditions.
