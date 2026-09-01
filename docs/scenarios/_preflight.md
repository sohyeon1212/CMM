---
id: _preflight
title: Preflight — checks every scenario runs first
solver: LP
runtime: seconds
---

# Preflight

Run this before any scenario. Each check exists because skipping it produces a confident
wrong answer rather than an error: a model that does not grow, a product with no exchange, or
an expression table whose gene ids map to nothing will all "succeed" and mean nothing.

Record every result here in `00_provenance.json` (see `_reporting.md`).

---

## Step P0 — Resolve consequential ambiguity

**Goal.** Establish a shared run definition without asking the user for facts CMM can inspect.

Before creating a config or output directory, inspect the candidate model and environment
read-only. Confirm which model files exist, load the intended model when its path is known,
enumerate exchanges and objectives, inspect current substrate and oxygen bounds, and check
solver and R capabilities. Classify each required input as resolved, uniquely discoverable,
unresolved, or blocked. A uniquely discoverable fact is the agent's work; a biological choice is
the user's.

Ask only about unresolved decisions that would change the result. Resolve them in dependency
order: exact model bytes, product exchange when the scenario requires one, one coherent biological
condition, capability-driven changes in method scope, and finally any experimental constraints
the user actually mentioned.
Every question must state the inspected evidence, explain why the choice matters, and provide
compatible choices with the recommended option listed first and clearly labeled, its rationale,
and its effect on the analysis or claim. Ask dependent questions one at a time. If no option is
scientifically defensible from the available evidence, recommend pausing for the missing
information rather than inventing a default.

Do not ask about a model hash, a unique objective or exchange mapping, installed capabilities,
an output timestamp, or documented deterministic protocol defaults. Conversely, current model
bounds describe the file, not necessarily the experiment: `model-as-loaded` requires explicit
acceptance after those bounds are shown, and "use defaults" does not resolve medium, substrate
uptake, or aeration.

If this step required an interview, present the exact model, any scenario-required product
exchange, full condition, method scope or substitutions, and relevant experimental constraints as
a resolved run-definition summary. Obtain explicit confirmation before starting the numerical
workflow or writing run artifacts. Skip the interview and confirmation round when the original
request was already complete and internally consistent. For SC-01 production requests, the
detailed question and exit contract is in the production skill's
[`clarification-interview.md`](../../.agents/skills/cmm-production-engineering/references/clarification-interview.md).

---

## Step P1 — Load the model and set the solver

**Goal.** Have a model and know what the solver can do.

```python
from cobra.io import load_model, read_sbml_model
from cmm.core import model_fingerprint, solver_status

model = read_sbml_model("model.xml")     # or load_model("textbook") / load_model("iJO1366")
status = solver_status(model)
print(status.summary(), status.warning)
print(model_fingerprint(model))
```

**Decision rule.** Compare `status.capabilities` against the solver gate in `AGENTS.md` §2 for
the methods your scenario needs.

**Branch.**
- Missing QP/MIQP and the scenario needs it → use an LP-capable substitute only when that
  scenario permits it and the user explicitly requests or approves the changed scientific claim.
  Otherwise return to step P0 with a recommended capable-solver or installation option. The
  canonical SC-01 workflow never substitutes MOMA-L1 for MOMA-L2.
- MILP or importable `straindesign` missing and the run needs a growth-coupled design → the
  design search cannot run. Any additional backend requirement is surfaced by that backend;
  decide how to narrow the run explicitly rather than silently omitting the design stage.

**Pin the model.** Copy the model file into the run directory as `model/<model-id>.xml` (or
write it there with `cobra.io.write_sbml_model`) before anything else runs. A model *name* does
not identify a model: `load_model("textbook")` and `load_model("e_coli_core")` both return 95
reactions called `e_coli_core`, but they differ — one calls formate transport `FORti` and its
biomass reaction `Biomass_Ecoli_core`, the other `FORt` and `BIOMASS_Ecoli_core_w_GAM` — and
their fingerprints differ accordingly. Two runs of the same scenario that name the same model
are not comparable unless the file is the same file, and a design naming `FORti` cannot be
built by someone holding the other variant.

**Artifact.** `model/<model-id>.xml`, plus `model_id`, `model_sha256`, solver name and
capabilities → `00_provenance.json`. The fingerprint is recorded again after the medium is
applied (step P2), because it changes with the medium and is therefore evidence of the
condition as well as of the model.

---

## Step P2 — Apply the growth medium

**Goal.** Constrain the model to the environment the user actually means.

```python
from cmm.core import Condition, PRESET_MEDIA, ReactionBound, apply_medium

print(sorted(PRESET_MEDIA))

MEDIUM = "glucose_anaerobic"              # a scientific choice; state it in the report
apply_medium(model, MEDIUM)

CONDITION = Condition(                    # carried into every later call
    name=MEDIUM,
    bounds=(ReactionBound(reaction_id="EX_o2_e", lower_bound=0.0, upper_bound=0.0),),
    notes="anaerobic: oxygen uptake closed",
)
```

**Set the medium and the aeration once, here, before the first solve, and pass `CONDITION` to
every later call that accepts it.** The `aerobic=True|False` parameter was removed in 0.4.0:
it duplicated what the medium already says, and because `optknock`/`robustknock` never accepted
it, an `aerobic=False` argument alongside an aerobic medium silently produced an aerobic
design inside an anaerobic run. For an aerobic run use the aerobic preset and drop the
`bounds=` entry rather than leaving a contradicting oxygen bound behind.

**Decision rule.** The medium is a scientific choice, not a default. If the prompt or config
does not resolve medium, substrate uptake, and aeration, return to step P0 rather than choosing
one. Include a recommended coherent condition and explain its scientific consequence. If the
condition is already explicit and consistent, do not ask for redundant confirmation. Re-take the
model fingerprint *after* `apply_medium` — it changes with the medium, so it is evidence of
which condition the run used.

**Missing components are not silent.** A medium component with no matching exchange in the
loaded model is recorded in provenance under `dropped`; a growth-limiting one raises rather
than producing a quietly different experiment. Record both `applied` and `dropped`, not just
the preset key.

**Failure → action.** `apply_medium` raising on an unknown key means the preset does not
exist; print `PRESET_MEDIA` and ask rather than inventing bounds.

---

## Step P3 — Confirm the model grows

**Goal.** Establish that the wild type is viable before asking anything about engineering it.

```python
from cmm.core import fba, pfba

growth = fba(model)
print(growth.status, growth.objective_value)
minimal = pfba(model)
```

**Decision rule.** `growth.status == "optimal"` and `objective_value > 1e-6`.

**Failure → action.** Zero or infeasible growth means the medium is closed, an essential
exchange is blocked, or the objective is unset. Report which and stop — every downstream
number would be meaningless. Do not "fix" it by opening bounds the user did not ask for.

**Artifacts.** `01_preflight/preflight.csv`; the canonical workflow exports the full wild-type
flux state as `03_reference/wild_type_fluxes.csv` and records growth in the summary/provenance.

---

## Step P4 — Confirm the target product is reachable

**Goal.** Establish that production design is even possible. Skip only when the run has no
target product (an essentiality study, or a condition comparison).

```python
from cmm.features import theoretical_yield

print(len(model.exchanges))
result = theoretical_yield(model, "EX_succ_e", condition=CONDITION)
print(result.status, result.molar_yield, result.exceeds_carbon_ceiling, result.co2_fixed)
```

**Decision rule.** `model.exchanges` non-empty, the product id is among them, and
`molar_yield > 1e-6`.

**Branch.**
- No exchanges → production design unavailable. Stop and tell the user (`AGENTS.md` §5).
- Product not an exchange id → inspect exchange ids, names, metabolites, and compartments. If
  the supplied name maps strictly and uniquely to one extracellular exchange, state that
  discoverable resolution and carry it into the run-definition summary. If zero or multiple
  plausible exchanges remain, ask which measured/exported species the user means; show the
  candidates and a recommended choice only when model evidence supports one.
- `molar_yield == 0` → the product is unreachable in this medium. Return to step P0, recommend
  the next condition change only when reachability evidence supports it, and ask whether to
  change medium, substrate, or aeration. Do not continue the target search.
- `exceeds_carbon_ceiling and not co2_fixed` → suspicious; the model may have an unbalanced
  reaction or an open sink. Flag it in the report.

---

## Step P5 — Check gene id overlap

**Goal.** Only for runs using expression data — `SC-02`, or an omics reference in `SC-01`. An
expression table whose ids do not match the model contributes nothing, silently.

```python
from cmm.omics import read_expression_table

expression = read_expression_table("expression.csv")
model_genes = {g.id for g in model.genes}
overlap = model_genes & set(expression.index)
print(len(overlap), "of", len(model_genes), "model genes matched")
```

**Decision rule.** Overlap should cover most of the model's genes. Below roughly half, treat
the table as unmapped.

**Failure → action.** Report the overlap fraction and a few unmatched ids from each side, then
ask the user for the id mapping. GEO tables typically use symbols or probe ids while an
*E. coli* model uses b-numbers. Do not proceed with a near-zero overlap — the integration will
return a plausible flux state built on almost no data.

Also required: expression values finite, non-negative, no duplicate gene rows.

---

## Preflight summary block

Put this at the top of the report:

| Check | Result |
|---|---|
| Model | id, reactions/metabolites/genes, SHA-256 |
| Solver | name, capabilities, substitutions made |
| Medium | preset or condition applied |
| Wild-type growth | value h⁻¹ |
| Product | exchange id, theoretical yield mol/mol, carbon ceiling flags |
| Expression | genes matched / model genes (if applicable) |
