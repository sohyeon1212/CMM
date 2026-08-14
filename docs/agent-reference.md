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

applied = apply_medium(model, "glucose_aerobic")   # keys: PRESET_MEDIA
applied.dropped                            # components this model has no exchange for
applied.to_provenance()                    # {"medium", "applied", "dropped"}
```

`apply_medium` / `Medium.apply_to` return a **`MediumApplication`** (0.4.0; previously a plain
dict). It still behaves as the mapping of applied `{exchange_id: uptake}`, so `len(...)`,
indexing and iteration are unchanged. A missing *growth-limiting* component raises rather than
producing a quietly different experiment; anything else warns and is listed in `.dropped`.
Record `.to_provenance()`, not just the preset key: on `e_coli_core` a preset legitimately
drops 18 of its 24 declared components.

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

**As of 0.4.0 there is one convention: `condition=`.** The `aerobic=True|False` parameter was
removed from `cmm.features.production` — it was a redundant second way of saying what the
medium already says, and because `optknock`/`robustknock` never accepted it, a caller could set
`aerobic=False` and still get an aerobic design with no warning. Set the medium and aeration
once, before the first solve, and pass the same `Condition` to every call. `substrate=` is
unchanged.

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
fva(model, condition=None, reactions=None, fraction_of_optimum=1.0, *,
    loopless=False, processes=None)                          # -> FvaResult
```

- `FluxSolution`: `.status`, `.objective_value` (None when not optimal), `.fluxes`, `.metadata`.
- `FluxRange`: `.minimum`, `.maximum`.
- `FvaResult` (0.4.0; `fva` previously returned a bare `dict[str, FluxRange]` with no
  provenance): `.ranges`, `.metadata` (`run_provenance`, exactly as `fba`/`pfba` carry it),
  `.to_frame()` → `reaction_id | minimum | maximum`. It is a `Mapping`, so `result[rid]`,
  `.get(rid)`, `len(result)`, `dict(result)` and iteration behave as before.
- `loopless` is forwarded (`True` → cobra's `"cycleFreeFlux"`; `"fastSNP"` gives optimal
  loopless bounds). Before 0.4.0 loopless FVA was unreachable.
- `processes=None` decides from the problem size: below 500 analysed reactions it resolves to
  1, which avoids ~2.9 s of pure process-pool overhead on `e_coli_core` **and** the macOS
  failure where an unguarded module-level FVA re-spawns interpreters forever. At genome scale
  the pool is still used, so a full iJO1366 FVA on macOS still needs an
  `if __name__ == "__main__":` guard — or pass `processes=1`.
- `reference_state_pfba(model, condition=None, name="reference", fraction_of_optimum=1.0)`
  returns a `FluxState` directly.

`FluxState` is the currency for every comparison method: `.fluxes`, `.name`, `.provenance`,
`.get(rid, default)`, `.reactions()`, `.to_series()`, `.distance(other, order=2)`,
`.serialize()` / `.deserialize()`.

---

## 3. Production design — `cmm.features.production`

```python
from cmm.features import fseof, fvseof, production_envelope, theoretical_yield

theoretical_yield(model, product, substrate=None, *, condition=None)
production_envelope(model, product, *, objective=None, substrate=None, condition=None,
                    points=20)
fseof(model, product, biomass=None, *, n_steps=10, fraction_min=0.1, fraction_max=0.9,
      condition=None, reactions=None, tol=1e-3)
fvseof(model, product, biomass=None, *, n_steps=10, fraction_min=0.1, fraction_max=0.9,
       biomass_fraction=0.95, condition=None, reactions=None,
       linear_flux_couplings=None, tol=1e-3)
```

`product` is an **exchange reaction id** (`EX_succ_e`), not a metabolite id. `substrate=None`
auto-detects the carbon source from the medium.

- `ProductionYield`: `.molar_yield`, `.product_flux`, `.substrate_uptake`, `.status`,
  `.aerobic`, `.carbon_ceiling`, `.co2_exchange`, `.product_carbon`, `.carbon_uptake`
  (every carbon-bearing uptake as `CarbonUptake(reaction_id, uptake, carbon_atoms)`), and the
  derived properties `.co_substrates`, `.co2_uptake`, `.co2_fixed`, `.product_carbon_flux`,
  `.co2_carbon_fraction`, `.excess_carbon`, `.exceeds_carbon_ceiling`, `.co2_explains_excess`,
  `.carbon_imbalance`.
  **Report `co2_carbon_fraction`, not `co2_fixed`.** `co2_fixed` is a bare boolean and the only
  active check was the carbon ceiling, so a CO₂-inflated yield *below* the ceiling passed every
  guard: anaerobic `e_coli_core`/`EX_succ_e` read 1.3906 against a ceiling of 1.5 with 12.5% of
  the product carbon coming from CO₂ uptake. From 0.4.0 the media presets close CO₂ uptake
  (that case now reads **1.2000**), `co2_carbon_fraction` quantifies any remaining contribution
  and a `UserWarning` is raised when it is non-zero, and `carbon_imbalance` flags an excess that
  CO₂ does *not* explain. The carbon ceiling is computed over **every** carbon-bearing uptake,
  so a co-fed model is measured against the carbon it actually receives.
  `theoretical_yield` also raises on a non-boundary reaction rather than returning a
  meaningless number.
- `ProductionEnvelope`: `.points`, `.max_growth`, `.max_product`, `.to_frame()` →
  `product_flux | growth_min | growth_max`.
- `FseofResult`: `.trends` (DataFrame indexed by reaction, one column per enforced level plus
  `classification` and `actionable`), `.enforced_levels`, `.amplification_targets()`,
  `.knockout_targets()`. Both target lists take `actionable_only=True` by default, which
  excludes boundary, objective, and no-GPR reactions — keep it on when proposing interventions.
- `FvseofResult`: `.mean` (Park's V_avg), `.forced` (CMM's forced-minimum |flux|), `.capacity`
  (Park's l_sol) — all reaction × level — plus `.classification`, `.park_type`, `.robust`,
  `.slope`, `.capacity_slope`, `.actionable`, and the accessors `.amplification_targets()`,
  `.knockout_targets()`, `.targets_of_type(*types)`, `.robust_targets()`.
  **`amplification_targets()` is the published selection**: Park et al.'s types 1–3, taken from
  the *joint* sign of ΔV_avg and Δl_sol, returned in Park's own priority order (ascending mean
  l_sol). `park_type` carries the 1–9 index per reaction.
  **`robust_targets()` is CMM's own construct, not Park's variability criterion**, and must not
  be attributed to them: it lists amplification targets whose *forced minimum* flux also rises
  monotonically, meaning the reaction cannot avoid carrying more flux. It is useful — FSEOF
  alone cannot tell that — but it is not what Park's Δl_sol signal measures.
  `linear_flux_couplings` are caller-supplied linear equalities `Σc·v = 0` (renamed from
  `group_constraints` in 0.4.0); Park's grouping-reaction constraints are STRING-derived on/off
  pairs with a normalised-flux inequality and are **not** implemented.

---

## 4. Flux response — `cmm.features.response`

```python
from cmm.features import flux_response

flux_response(model, target, response=None, *, biomass=None, condition=None,
              target_min=None, target_max=None, n_steps=20, biomass_fraction=None,
              limit_threshold=1e-6, tol=1e-9)
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

`FluxResponseResult`: `.target`, `.response`, `.biomass`, `.points`, `.phases`, `.limit`,
`.wild_type` (dict of target/response/biomass flux at the growth optimum), `.metadata`,
`.to_frame()` → `target_flux | response_flux | biomass_flux | status`, `.phases_frame()`,
`.feasible_points()`, `.feasible_range()` → `(low, high)` or `None`, `.optimum()` → best
`ResponsePoint` or `None`, `.feasible_domain`, `.shadow_price_at(target_flux)` → exact
`d(response)/d(target)`. The phases are built from LP duals, so `shadow_price_at` is a lookup
in that structure and needs no further solve; the result stays serializable and model-free.
It raises when `target_flux` is outside the feasible domain.

`ResponsePoint`: `.target_flux`, `.response_flux`, `.biomass_flux`, `.status`, `.feasible`.
Infeasible points stay in the result with NaN fluxes — **this is data, not an error.**

`ResponsePhase`: `.target_low`, `.target_high`, `.shadow_price`, `.response_low`,
`.response_high` — one linear piece of the response curve, i.e. an interval of constant shadow
price. `ResponseLimit`: `.found`, `.message`, `.target_flux`, `.response_flux`,
`.shadow_price_before`, `.shadow_price_after`, `.threshold` — the phase boundary at which the
response starts to fall faster than `threshold` (set it with `limit_threshold=`, default
`1e-6`). `found=False` is a real finding: either the
response never declines (the target does not limit it) or the curve is flat (the objective is
insensitive to that reaction). Read `.message` and report which it was.

> **Removed in 0.4.0: `.bottleneck` / `ResponseBottleneck`.** It located the steepest decline of
> a finite-difference gradient. The response curve is an LP optimal-value function in one bound
> and is therefore concave piecewise linear, so that argmin finds the edge of the scan grid: the
> reported location moved by up to 29.53 flux units and the `found` flag inverted between
> `n_steps` 6 and 160. There is **no published criterion** defining a bottleneck this way. Use
> `.limit` and `.phases`, and do not cite a paper for the removed field.

`feasible_range()` is FVA-derived from 0.4.0. Earlier versions read it off the scan grid and
understated the range whenever a growth floor was applied, so do not quote one from an older run.

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
                    room_use_case="flux_prediction", delta=None, epsilon=None)
batch_comparison(model, reference, perturbations, *, method="moma_l2",
                 room_use_case="lethality", delta=None, epsilon=None,
                 objective_reaction=None, product_reaction=None)
moma(model, reference, *, linear=False)
room(model, reference, *, linear=False, use_case="flux_prediction",
     delta=None, epsilon=None)
```

- `reference_flux` methods: `"fba"`, `"pfba"` (model only) or `"lad"`, `"eflux2"` (need
  `gene_expression={gene_id: value}`). **The reference choice changes every downstream
  result** — state which one you used.
- `knockout_comparison` methods: `"moma_l2"` (QP), `"moma_l1"` (LP), `"room"` (MILP).
- **ROOM tolerances are a use-case choice, not a constant.** Shlomi et al. give δ=0.03/ε=1e-3
  for *flux prediction* and δ=0.1/ε=1e-2 for *lethality*; `cmm.features.ROOM_TOLERANCES` holds
  both. `delta`/`epsilon` default to `None` and override the preset when given. A single
  named perturbation is a flux prediction, so `room()` and `knockout_comparison()` default to
  that pair; `batch_comparison` is a screen and defaults to `"lethality"`. The choice is not
  cosmetic: 531 against 401 switches over the same 35 genes, a 24% shift in the ranking score.
- Perturbation lists: `reaction_perturbations(model)`, `gene_perturbations(model)`,
  `grouped_gene_perturbations(...)`. For a specific gene set,
  `blocked_reactions_for_genes(model, ["b0726"])` resolves GPRs to reaction ids you can pass
  straight to `knockout_comparison`.
- `ComparisonResult`: `.method`, `.status`, `.objective_value`, `.distance`, `.distance_kind`,
  `.n_changed_reactions`, `.fluxes`, `.metadata`, `.to_flux_state()`.
- `BatchComparisonRow`: `.target_id`, `.kind`, `.status`, `.objective_value`, `.distance`,
  `.distance_kind`, `.n_changed_reactions`, `.objective` (growth), `.n_reactions`,
  `.product_flux` (NaN unless `product_reaction=` was given). Pass `product_reaction` whenever
  you are screening for production — without it you only learn which knockouts hurt growth.
- **`objective_value` and `distance` are different quantities and 0.4.0 separates them.**
  `objective_value` is the raw solver objective and means something different per method: `Σd²`
  for `moma_l2`, `Σ|d|` for `moma_l1`, and a *count of switched reactions* for `room`.
  `distance` is a distance and only a distance — Segrè et al. Eq. (4)'s Euclidean `√(Σd²)` for
  `moma_l2`, the L1 sum for `moma_l1`, and **`None` for ROOM**, which defines none; ROOM's count
  is on `n_changed_reactions`. `distance_kind` is `"euclidean_l2" | "l1" | "none"`. Before
  0.4.0 one `distance` field held all three: on the SC-01 design it reported **1303.99** where
  the Euclidean distance is **36.11**. A NaN `distance` still means "no solution"; `None` means
  "not a distance for this method" — code that formats it must handle `None`.
- The perturbation builders return a `PerturbationList` (a `list` subclass): iteration, `len()`
  and indexing are unchanged, and `.inert_dropped` / `.n_inert_dropped` / `.provenance()` record
  the genes that block no reaction and were dropped (66 of 137 on `e_coli_core`).

MOMA/ROOM assume the cell makes a *minimal adjustment* from wild type. OptKnock assumes it
reaches its *growth-maximizing* state. These are different cell models; do not mix their
conclusions without saying so.

---

## 7. Strain design — `cmm.features.strain_design`

```python
from cmm.features import optknock, robustknock

optknock(model, product, *, biomass=None, max_knockouts=3, max_solutions=5,
         min_growth=0.05, condition=None, actionable_only=True)
robustknock(model, product, *, biomass=None, max_knockouts=3, max_solutions=8,
            min_growth=0.05, condition=None, actionable_only=True)
```

`StrainDesignResult`: `.designs`, `.best()`. `StrainDesign`: `.knockouts`, `.growth`,
`.max_product` (optimistic), `.guaranteed_product` (worst case), `.growth_coupled`.

**Rank by `guaranteed_product`, never `max_product`.** A design is only valuable if the cell
*cannot* maximize growth without producing; `max_product` is what the cell could do if it
chose to cooperate. `growth_coupled` is `guaranteed_product > 0`.

`condition=` is new in 0.4.0 and matters: before it, these two functions took no aeration or
medium argument at all and depended silently on the caller's model state, which is how an
aerobic design came to be documented as an anaerobic result (see `SC-01` step 3).

`actionable_only=True` (0.4.0) restricts the candidate set to gene-associated internal
reactions. Without it the search proposes deleting boundary exchanges with no GPR
(`EX_co2_e`, `EX_ac_e`, `EX_for_e`, `EX_etoh_e`, `EX_lac__D_e`), which are not realisable as
gene deletions; on `e_coli_core` it removes 26 of 95 candidates and 10 of 18 returned designs
without changing the top design or its numbers. Pass `actionable_only=False` only if you
genuinely want exchange knockouts.

`max_solutions` caps **MILP solutions, not distinct designs** — coupled-set members can yield
the same intervention more than once, and one run returned 23 designs at `max_solutions=5`.
The solution-pool size is solver-state dependent and is not a reportable quantity; the top
design is.

Needs a MILP solver **and** the `straindesign` package (which needs Java/OpenJDK). Cite
Schneider et al. (2022) for `straindesign` alongside Burgard et al. (2003) and
Tepper & Shlomi (2010) — the package carries no citation of its own.

---

## 8. Omics integration — `cmm.omics`

```python
from cmm.omics import (
    differential_expression, flux_log_change, integrate_expression,
    predict_condition_fluxes, read_expression_table, sign_flips,
)

integrate_expression(model, gene_expression, *, method="eflux2", or_rule="sum", **kwargs)
gene_to_reaction_weights(model, gene_expression, *, or_rule="sum")
eflux2(model, reaction_weights, *, objective_fraction=1.0, min_scale=1e-3,
       weight_threshold=0.0, normalization_percentile=100.0, exclude_exchange=True,
       excluded_reactions=None, allow_l1_fallback=False, gpr_or_rule=None)
lad(model, reaction_weights, *, scaling_factor=1.0, weight_threshold=0.0,
    reaction_sigma=None, gpr_or_rule=None)
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
- **`or_rule` resolves GPR `OR` (isozymes) and defaults to `"sum"` from 0.4.0**, which is what
  both Kim et al. 2016 (E-Flux2) and Lee et al. 2012 (LAD) specify. `AND` is `min` throughout
  and was already correct. `or_rule="max"` is CMM's pre-0.4.0 behaviour, matches no source
  paper, and exists only to reproduce an old result; the rule used is archived in the result's
  provenance as `gpr_or_rule`. On `e_coli_core` + GSE41189 the change moves 30 of 66 mapped
  reaction weights (up to 2.67×) and predicted growth by 26%.
- `lad`'s `weight_threshold` defaults to **0.0** from 0.4.0 (was 0.01), so a low-expression
  reaction is driven toward zero flux as Lee et al. intend rather than dropped from the
  objective. `reaction_sigma` supplies Lee et al.'s per-reaction `1/σ` weights; the unweighted
  default is a deviation from the paper and is recorded as `sigma_weighted` in metadata.
- Every `OmicsFluxResult` carries `metadata["cmm_deviations"]`
  (`cmm.omics.EFLUX2_DEVIATIONS` / `LAD_DEVIATIONS`) listing where the implementation departs
  from its source. Read it before quoting a number.
- `differential_expression` returns a `DirectionMap` for the revert/transform methods. Its GPR
  rule is Yizhak et al.'s ternary rule — all subunits changed (AND), at least one changed (OR),
  **mixed ⇒ unchanged** — recorded as `metadata["gpr_rule"]`. This is a *different* operation
  from the continuous `or_rule` above and is deliberately not shared with it.
- **`flux_log_change` implements no published method.** It is a CMM utility — a log2 ratio of
  two flux vectors with a pseudocount — and must not be cited to any paper. Report the
  pseudocount, which fixes the value returned for switch-on and switch-off reactions.
- Cite E-Flux2 as Kim MK, Lane A, Kelley JJ, Lun DS (2016) and LAD as Lee D, Smallbone K,
  Dunn WB et al. (2012), *BMC Syst Biol* 6:73. Full entries in `docs/VALIDATION.md`.

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
  state, **QP**) or `"mta"` (**MIQP**). `"moma"` remains the default.
- **Check the tie structure before quoting a top-k.** `TargetRanking.metadata` carries
  `n_distinct_scores`, `largest_tie_block` and `score_resolution` (also available from
  `cmm.features.tie_structure`), because `TargetRanking.sorted` breaks ties alphabetically on
  `target_id`. On the project's own SC-02 pair the `mta` path had a single tie block of 18 of
  71 genes even after 0.4.0 floored the score denominator and removed a 38-gene `+∞` block.
- `direction` is the *goal*, not the intervention. The top-ranked knockout need not be one of
  the differentially expressed genes.
- rMTA is optimistic per knockout — a prioritization tool robustified by its worst-case term,
  not a proof.
- **`transformation_targets` is not a CMM invention.** Both paths map to published Yizhak et al.
  (2013) methods: `mta` is that paper's MIQP applied to an arbitrary source→target pair, and
  `moma` is the distance-reduction scoring it uses as the comparison method. Cite Yizhak et al.
  (2013) for both. The paper reports the MOMA-style scoring as *markedly inferior* to MTA, and
  it is the default here, so state which path produced a ranking.

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
| `flux_response_figure(result)` | `FluxResponseResult` | response curve, infeasible span, optimum, response limit, growth axis |
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
