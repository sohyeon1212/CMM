# CMM function reference for agents

Signatures and result objects for every shipped CMM service. `AGENTS.md` decides *which*
analysis to run; this file tells you how to call it. Read the section for the function you
are about to use rather than the whole file.

Everything runs on a plain `cobra.Model`. No CMM function mutates the model destructively:
each opens a `with model:` context and restores bounds, objective, and medium on exit.

---

## 1. Setup

```python
from cobra.io import load_model, read_sbml_model
from cmm.core import PRESET_MEDIA, apply_medium, solver_status, supports

model = load_model("textbook")            # or read_sbml_model("model.xml")
model.solver = "gurobi"                    # if available; see the solver gate in AGENTS.md

status = solver_status(model)              # .name .capabilities .recommended .available
status.warning                             # human-readable gap, or None if fully capable
status.summary()                           # e.g. "gurobi (recommended): LP, MILP, MIQP, QP"

apply_medium(model, "glucose_aerobic")     # keys: PRESET_MEDIA
```

Capability checks: `supports("QP", model.solver.interface)` returns a bool;
`require("QP", model.solver.interface, feature="L2 MOMA")` raises `SolverCapabilityError`
with an actionable message.

### Conditions — reusable constraint sets

Most services accept `condition=` instead of you editing bounds by hand. This is how you
express anaerobic, a substrate limit, or an alternative objective without touching the model.

```python
from cmm.core import Condition, ObjectiveSpec, ReactionBound

anaerobic = Condition(
    name="anaerobic",
    bounds=(ReactionBound(reaction_id="EX_o2_e", lower_bound=0.0),),
)
```

`Condition(name, bounds=(), objective=None, notes="")`; `ReactionBound(reaction_id,
lower_bound=None, upper_bound=None)` leaves a `None` side untouched; `ObjectiveSpec(
coefficients, direction="max")`.

Note the split convention: `cmm.core` simulation and the newer services take `condition=`,
while `cmm.features.production` functions take `aerobic=True|False` and `substrate=` directly.

### Provenance

```python
from cmm.core import model_fingerprint, run_provenance

model_fingerprint(model)      # SHA-256 over bounds, GPRs, stoichiometry, objective
run_provenance(model, method="my_run", note="...")
```

Every numerical result already carries this in `result.metadata`. Never report a number
without it — see the run contract in `AGENTS.md`.

---

## 2. Simulation — `cmm.core`

```python
from cmm.core import fba, fva, pfba, reference_state_pfba

fba(model, condition=None)                                   # -> FluxSolution
pfba(model, condition=None, fraction_of_optimum=1.0)         # -> FluxSolution
fva(model, condition=None, reactions=None, fraction_of_optimum=1.0)   # -> dict[str, FluxRange]
```

- `FluxSolution`: `.status`, `.objective_value` (None when not optimal), `.fluxes`, `.metadata`.
- `FluxRange`: `.minimum`, `.maximum`.
- `reference_state_pfba(model, condition=None, name="reference", fraction_of_optimum=1.0)`
  returns a `FluxState` directly.

`FluxState` is the currency for every comparison method: `.fluxes`, `.name`, `.provenance`,
`.get(rid, default)`, `.reactions()`, `.to_series()`, `.distance(other, order=2)`,
`.serialize()` / `.deserialize()`.

---

## 3. Production design — `cmm.features.production`

```python
from cmm.features import fseof, fvseof, production_envelope, theoretical_yield

theoretical_yield(model, product, substrate=None, *, aerobic=True)
production_envelope(model, product, *, objective=None, substrate=None, aerobic=True, points=20)
fseof(model, product, biomass=None, *, n_steps=10, fraction_min=0.1, fraction_max=0.9,
      aerobic=True, reactions=None, tol=1e-3)
fvseof(model, product, biomass=None, *, n_steps=8, fraction_min=0.1, fraction_max=0.9,
       biomass_fraction=0.95, aerobic=True, reactions=None, group_constraints=None, tol=1e-3)
```

`product` is an **exchange reaction id** (`EX_succ_e`), not a metabolite id. `substrate=None`
auto-detects the carbon source from the medium.

- `ProductionYield`: `.molar_yield`, `.product_flux`, `.substrate_uptake`, `.status`,
  `.carbon_ceiling`, `.co2_exchange`, `.co2_fixed`, `.exceeds_carbon_ceiling`, `.aerobic`.
  Always report `exceeds_carbon_ceiling` and `co2_fixed`: a yield above the substrate's carbon
  ceiling is only legitimate when CO₂ is being fixed.
- `ProductionEnvelope`: `.points`, `.max_growth`, `.max_product`, `.to_frame()` →
  `product_flux | growth_min | growth_max`.
- `FseofResult`: `.trends` (DataFrame indexed by reaction, one column per enforced level plus
  `classification` and `actionable`), `.enforced_levels`, `.amplification_targets()`,
  `.knockout_targets()`. Both target lists take `actionable_only=True` by default, which
  excludes boundary, objective, and no-GPR reactions — keep it on when proposing interventions.
- `FvseofResult`: `.mean`, `.forced`, `.capacity` (all reaction × level), plus
  `.amplification_targets()`, `.knockout_targets()`, `.robust_targets()`. **`robust_targets()`
  is the one that matters**: it lists reactions whose *forced minimum* flux rises with enforced
  product, meaning the reaction cannot avoid carrying more flux. FSEOF alone cannot tell that.

---

## 4. Flux response — `cmm.features.response`

```python
from cmm.features import flux_response

flux_response(model, target, response=None, *, biomass=None, condition=None,
              target_min=None, target_max=None, n_steps=20, biomass_fraction=None, tol=1e-9)
```

Fixes `target` at each point of a linear scan and maximizes `response` there.

- `response=None` maximizes the objective — "how sensitive is growth to this reaction, and
  where does it break?"
- `response="EX_succ_e"` maximizes a product, recording biomass at every point.
- **`biomass_fraction` is required in practice for a product response.** Without a growth
  floor the solver returns non-growing solutions, so the curve is a theoretical ceiling rather
  than a strain. `biomass_fraction=0.3` holds biomass at 30% of the wild-type optimum.
- Range defaults to the target's full feasible interval (FVA at a zero fraction of the
  optimum). An explicit range may exceed the reaction's declared bounds; that is allowed as a
  what-if and recorded as `range_outside_bounds` in provenance.

`FluxResponseResult`: `.target`, `.response`, `.biomass`, `.points`, `.bottleneck`,
`.wild_type` (dict of target/response/biomass flux at the growth optimum), `.metadata`,
`.to_frame()` → `target_flux | response_flux | biomass_flux | status`, `.feasible_points()`,
`.feasible_range()` → `(low, high)` or `None`, `.optimum()` → best `ResponsePoint` or `None`.

`ResponsePoint`: `.target_flux`, `.response_flux`, `.biomass_flux`, `.status`, `.feasible`.
Infeasible points stay in the result with NaN fluxes — **this is data, not an error.**

`ResponseBottleneck`: `.found`, `.message`, `.target_flux`, `.response_flux`,
`.steepest_decline`, `.decline_interval`, `.sensitivity`. `found=False` is a real finding:
either the response never declines (the target does not limit it) or the curve is flat (the
objective is insensitive to that reaction). Read `.message` and report which it was.

---

## 5. Random flux sampling — `cmm.features.sampling`

```python
from cmm.features import random_flux_sampling, reference_constrained_sampling

random_flux_sampling(model, n=1000, *, condition=None, method="optgp",
                     thinning=100, processes=1, seed=0)
reference_constrained_sampling(model, reference, n=1000, *, condition=None,
                               min_fraction=0.8, max_fraction=1.2,
                               zero_tolerance=1e-6, zero_window=0.1,
                               method="optgp", thinning=100, processes=1, seed=0)
```

- `method="optgp"` needs roughly >1000 samples to mix; `method="achr"` converges better for
  small runs and is single-process only.
- `processes=1` is the default deliberately: parallel chains are seeded independently, so a
  multi-process run is not bit-for-bit reproducible. `metadata["parameters"]["reproducible"]`
  records which you got.
- `reference` is a `FluxState` or a plain `{reaction_id: flux}` mapping. Each reaction is
  narrowed to a window around its reference flux, intersected with existing bounds. A
  reference from a *different* condition can violate the sampled model's bounds; that raises
  with the offending reactions named rather than silently sampling a different space.
- Sampling the model as constrained: apply a medium or `condition=` first, or you sample an
  unbounded space.

`SamplingResult`: `.samples` (DataFrame, rows = samples, columns = reactions), `.method`,
`.seed`, `.n_samples`, `.metadata`, `.to_frame()`, `.statistics()` → per-reaction
`mean | std | minimum | q1 | median | q3 | maximum`, `.correlation(min_std=1e-6)` (drops
constant reactions rather than reporting NaN), `.to_flux_state(name="sampled")`.

CMM provides no post-hoc noise addition on purpose: perturbing sampled fluxes independently
breaks `S · v = 0` and the bounds, so the results would no longer be flux distributions.

---

## 6. Perturbation response — `cmm.features.comparison`

```python
from cmm.features import (
    batch_comparison, gene_perturbations, knockout_comparison,
    reaction_perturbations, reference_flux,
)
from cmm.features._perturbation import blocked_reactions_for_genes

reference_flux(model, method="pfba", *, gene_expression=None, name=None)   # -> FluxState
knockout_comparison(model, reference, reaction_ids, *, method="moma_l2",
                    delta=0.03, epsilon=1e-3)
batch_comparison(model, reference, perturbations, *, method="moma_l2", delta=0.03,
                 epsilon=1e-3, objective_reaction=None, product_reaction=None)
```

- `reference_flux` methods: `"fba"`, `"pfba"` (model only) or `"lad"`, `"eflux2"` (need
  `gene_expression={gene_id: value}`). **The reference choice changes every downstream
  result** — state which one you used.
- `knockout_comparison` methods: `"moma_l2"` (QP), `"moma_l1"` (LP), `"room"` (MILP).
- Perturbation lists: `reaction_perturbations(model)`, `gene_perturbations(model)`,
  `grouped_gene_perturbations(...)`. For a specific gene set,
  `blocked_reactions_for_genes(model, ["b0726"])` resolves GPRs to reaction ids you can pass
  straight to `knockout_comparison`.
- `ComparisonResult`: `.method`, `.status`, `.distance`, `.fluxes`, `.to_flux_state()`.
- `BatchComparisonRow`: `.target_id`, `.kind`, `.status`, `.distance`, `.objective` (growth),
  `.n_reactions`, `.product_flux` (NaN unless `product_reaction=` was given). Pass
  `product_reaction` whenever you are screening for production — without it you only learn
  which knockouts hurt growth.

MOMA/ROOM assume the cell makes a *minimal adjustment* from wild type. OptKnock assumes it
reaches its *growth-maximizing* state. These are different cell models; do not mix their
conclusions without saying so.

---

## 7. Strain design — `cmm.features.strain_design`

```python
from cmm.features import optknock, robustknock

optknock(model, product, *, biomass=None, max_knockouts=3, max_solutions=5, min_growth=0.05)
robustknock(model, product, *, biomass=None, max_knockouts=3, max_solutions=8, min_growth=0.05)
```

`StrainDesignResult`: `.designs`, `.best()`. `StrainDesign`: `.knockouts`, `.growth`,
`.max_product` (optimistic), `.guaranteed_product` (worst case), `.growth_coupled`.

**Rank by `guaranteed_product`, never `max_product`.** A design is only valuable if the cell
*cannot* maximize growth without producing; `max_product` is what the cell could do if it
chose to cooperate. `growth_coupled` is `guaranteed_product > 0`.

Needs a MILP solver **and** the `straindesign` package (which needs Java/OpenJDK).

---

## 8. Omics integration — `cmm.omics`

```python
from cmm.omics import (
    differential_expression, flux_log_change, integrate_expression,
    predict_condition_fluxes, read_expression_table, sign_flips,
)

integrate_expression(model, gene_expression, *, method="eflux2")   # or "lad"
read_expression_table(path, gene_column=None)                       # -> gene x condition frame
predict_condition_fluxes(model, expression, *, method="eflux2", conditions=None)
flux_log_change(source_fluxes, target_fluxes, *, reactions=None, pseudocount=1e-3)
sign_flips(source_fluxes, target_fluxes, *, reactions=None, tol=1e-6)
differential_expression(model, source, target, *, reference=None, up_threshold=1.0,
                        down_threshold=1.0, pseudocount=1.0, reactions=None)
```

- `OmicsFluxResult`: `.method`, `.status`, `.objective_value`, `.fluxes`, `.detail`,
  `.metadata`, `.to_flux_state(name="omics")`.
- `ConditionFluxes`: `.conditions()`, `.fluxes(condition)` (raises if that condition's solve
  was not optimal).
- E-Flux2 scales bounds by normalized expression then minimizes total squared flux (**QP**);
  LAD fits fluxes to expression-derived targets (**LP**). Under a QP-less solver use LAD, and
  say so — do not pass `allow_l1_fallback=True` and call the result E-Flux2.
- `differential_expression` returns a `DirectionMap` for the revert/transform methods.

**Gene ids must match `model.genes`.** GEO tables use symbols or probe ids; an *E. coli* model
uses b-numbers. Unmapped ids silently contribute nothing, so an unmapped table produces a
confident, meaningless answer. Check the overlap before integrating.

---

## 9. Reversion and transformation — `cmm.features.revert`, `.transformation`

```python
from cmm.features import direction_from_states, revert_targets, transformation_targets

revert_targets(model, source_condition, reference_state, direction, *, targets=None,
               method="rmta", alpha=0.66, epsilon=1e-3, parameter_k=100.0,
               perturbation="gene", transcript_separator=None)
transformation_targets(model, source_state, target_state, *, method="moma",
                       perturbation="gene", targets=None, order=2, alpha=0.66)
direction_from_states(source, target, *, reactions=None, tol=1e-6)
```

Both return a `TargetRanking`: `.sorted(descending=True)`, `.top(n)`, `.best()`,
`.to_frame()`, `.to_records()`, iterable of `TargetScore`.

- `revert_targets` methods: `"rmta"`, `"mta"` (**MIQP**), `"rmta_continuous"` (**QP**, an
  explicitly labeled heuristic — never report it as published rMTA).
- `transformation_targets` methods: `"moma"` (score = reduction in distance to the target
  state, **QP**) or `"mta"` (**MIQP**).
- `direction` is the *goal*, not the intervention. The top-ranked knockout need not be one of
  the differentially expressed genes.
- rMTA is optimistic per knockout — a prioritization tool robustified by its worst-case term,
  not a proof.

---

## 10. Figures — `cmm.visualization`

```python
from cmm.visualization import (
    escher_flux_map, flux_comparison_figure, flux_log_change_figure, flux_response_figure,
    fseof_figure, fvseof_figure, network_flux_map, production_envelope_figure,
    sampling_figure, save_figure, yield_figure,
)

save_figure(fig, path, dpi=300)      # creates parent dirs, crops tightly, white facecolor
```

All figures are authored at 300 DPI with a colour-blind-safe palette and take
`column_width=1` (single column, ~3.3 in) or `2` (double, ~6.5 in).

| Figure | Takes | Shows |
|---|---|---|
| `production_envelope_figure(envelope)` | `ProductionEnvelope` | growth/product phase plane |
| `yield_figure([yields])` | list of `ProductionYield` | molar yield bars |
| `fseof_figure(result, top_n=6)` | `FseofResult` | target flux vs enforced product |
| `fvseof_figure(result, top_n=5)` | `FvseofResult` | mean (solid) and forced-min (dashed) flux |
| `flux_response_figure(result)` | `FluxResponseResult` | response curve, infeasible span, optimum, bottleneck, growth axis |
| `sampling_figure(result, top_n=8, reference=…)` | `SamplingResult` | per-reaction flux violins vs a reference |
| `flux_comparison_figure(ref, cmp, reactions)` | two flux dicts | grouped bars |
| `flux_log_change_figure(log_changes)` | dict from `flux_log_change` | ranked log2 changes |
| `network_flux_map(model, fluxes)` | model + fluxes | schematic force-directed network |
| `escher_flux_map(map_path, fluxes)` | an Escher map JSON | curated layout coloured by flux |

`network_flux_map` is a quick schematic, not a curated Escher map — do not present it as a
publication network figure.

---

## 11. Feature manifest

`cmm.features.INCLUDED_FEATURES` lists what actually ships; `PLANNED_FEATURES` lists what does
not. As of this writing `PLANNED_FEATURES` is `("dynamic_fba", "enzyme_constrained_modeling")`.
Check the tuple rather than assuming — and never present a planned feature as available.
