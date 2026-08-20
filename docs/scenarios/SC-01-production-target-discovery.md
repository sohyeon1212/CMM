---
id: SC-01
title: Production target discovery
goal: Find and verify genetic interventions that increase production of a target metabolite
when_to_use:
  - "increase production of X"
  - "which genes should I over-express / knock out to make more X"
  - "design a strain where production is guaranteed, not merely possible"
  - "생산 증대 표적을 찾아줘"
  - "성장 공역 균주 설계"
role: spine                 # the backbone of a production goal; complete on its own
entry_points:
  - "step 0 — the usual start, when no targets are known yet"
  - "step 3 — when only a growth-coupled knockout design is wanted"
optional_inputs_from:
  - "SC-02: a condition-specific flux state, used as a reference in step 2"
  - "SC-03: an exhaustive single-deletion screen, used as context for step 3"
requires:
  model: "cobra model"
  product: "exchange reaction id (e.g. EX_succ_e)"
optional:
  expression: "gene expression table (enables LAD / E-Flux2 reference states)"
  substrate: "exchange reaction id; default: auto-detected from the medium"
  condition: "cmm.core.Condition — medium and aeration, set once in step 0 and passed to every step"
solver:
  minimum: "LP — step 3 falls back to a single-deletion screen and proves no coupling"
  full: "MILP + straindesign + Java (step 3 design); QP (E-Flux2, L2 MOMA)"
steps: [preflight, yield, reference, design, amplification, validation, report]
runtime: "minutes on a core model; the step 3 MILP search dominates on genome scale"
---

# SC-01 — Production target discovery

## Objective

Produce a strain proposal for a target product: a knockout design whose production is
*guaranteed* at maximum growth, the reactions to over-express alongside it, the predicted
product and growth for each, and evidence that every part survived independent verification.

**Success criteria.**
1. Theoretical yield is known and non-zero, with its carbon balance disclosed.
2. Each proposed knockout design either has `guaranteed_product > 0` — coupling proven — or
   the report states plainly that coupling was **not** established and why.
3. Every recommended intervention passed step 5 verification, not just a step 3 or 4 ranking.
4. The run directory satisfies `_reporting.md`.

## Pipeline at a glance

| Step | Question | Method | Output |
|---|---|---|---|
| 0 | Is this model usable at all? | `_preflight.md` | viability, solver, yield |
| 1 | What is the ceiling, and what does it cost? | `theoretical_yield`, `production_envelope` | yield, trade-off curve |
| 2 | What is the cell doing now? | `reference_flux` | reference flux state(s) |
| 3 | Which knockout set forces production? | `optknock`, `robustknock` | designs with guaranteed product |
| 4 | Which reactions should carry more flux? | `fseof`, `fvseof` | amplification candidates |
| 5 | Do design and candidates survive scrutiny? | MOMA/ROOM, `flux_response`, sampling | verified strain proposal |
| 6 | Write it up | `_reporting.md` | report + figures + raw data |

Steps 3 and 4 are independent and may run in either order. Step 5 needs both.

## The condition is set once, in step 0, and every step inherits it

**This is the single most important procedural rule in this scenario, and the one this
document itself previously broke.** A run that applies an aerobic medium and then asks a later
step for an anaerobic answer produces numbers that are individually valid and collectively
meaningless. Only two things fix that: choose the medium and aeration *before* step 0's first
solve, and make every subsequent call state which condition it is running under.

```python
from cmm.core import Condition, ReactionBound, apply_medium

MEDIUM      = "glucose_anaerobic"        # or "glucose_aerobic"; keys in PRESET_MEDIA
PRODUCT     = "EX_succ_e"
OXYGEN      = "EX_o2_e"

apply_medium(model, MEDIUM)              # constrains the model itself

CONDITION = Condition(                   # restates the aeration explicitly, so every
    name=MEDIUM,                         # call below carries it in its own provenance
    bounds=(ReactionBound(reaction_id=OXYGEN, lower_bound=0.0, upper_bound=0.0),),
    notes="anaerobic: oxygen uptake closed",
)
```

For an aerobic run, use the aerobic preset and drop the `bounds=` entry — do not leave an
oxygen bound behind that contradicts the medium.

**From here on, `CONDITION` is passed to every call that accepts it**, and the calls below show
it explicitly rather than leaving it to be inferred. As of 0.4.0 that covers
`fba`/`pfba`/`fva`, `theoretical_yield`, `production_envelope`, `fseof`, `fvseof`, `optknock`,
`robustknock`, `flux_response`, `random_flux_sampling` and `reference_constrained_sampling`.
**The `aerobic=True|False` parameter no longer exists** — it was a second, redundant way of
saying what the medium already says, and because `optknock`/`robustknock` never accepted it,
passing `aerobic=False` alongside an aerobic medium silently produced an aerobic design inside
an anaerobic report. That is exactly how the wrong answer recorded in earlier versions of this
document was produced.

The remaining calls in this scenario — `reference_flux`, `knockout_comparison`,
`batch_comparison` — take no `condition=`. **They inherit the condition from the model state
`apply_medium` established in step 0**, which is why the medium is applied to the model rather
than only expressed as a `Condition`. Do not open a `with model:` block that changes the medium
around any of them.

Record `MEDIUM`, the oxygen bounds, the substrate and its uptake rate in `00_provenance.json`,
and re-record the model fingerprint *after* `apply_medium` — the fingerprint changes with the
medium, so it is itself evidence of which condition the run used.

**Two method families, two jobs.** Step 3 is *inverse* — you give a goal, CMM searches for the
intervention, and the bilevel MILP proves the cell cannot grow without producing. Step 5 is
*forward* — you give the intervention, CMM predicts what follows. They assume different cell
behaviour (growth-maximizing vs minimal-adjustment), which is exactly why one designs and the
other checks. Do not rank their outputs in a single list; see §Do not.

**Entering at step 3.** If the user only wants a growth-coupled knockout design and has no
interest in over-expression targets, run steps 0–1, then 3, then 5a, then report. Steps 2, 4
and 5b exist to serve the amplification half of the answer.

---

## Step 0 — Preflight

**Goal.** Establish that the model grows, the product is reachable, and the solver is capable —
**under the condition fixed above, not under the model's shipped defaults.**

**Call.** Apply the medium and build `CONDITION` first (see above), then follow `_preflight.md`
(P1–P4; P5 only if you have expression data), then check what step 3 will be able to do:

```python
import importlib.util
from cmm.core import supports

can_design = (
    supports("MILP", model.solver.interface)
    and importlib.util.find_spec("straindesign") is not None
)
print("step 3 mode:", "design (OptKnock)" if can_design else "screen (MOMA/ROOM fallback)")
```

**Branch.**
- Theoretical yield zero → stop and report; the product cannot be made in this medium.
  Everything below assumes a non-zero ceiling.
- The preflight growth number is the wild-type growth **in this condition**. On anaerobic
  `e_coli_core` that is 0.211663 h⁻¹, not the 0.873922 h⁻¹ of the aerobic model. Quoting the
  aerobic figure inside an anaerobic run is the first symptom of a mixed-condition report.
- `can_design` false → step 3 runs its LP fallback. Decide this now, not mid-run, and carry the
  fact into the report: the run will produce candidates rather than a coupled design.
  `straindesign` also needs Java, so a missing `java` on `PATH` has the same effect.

---

## Step 1 — Theoretical yield and the growth/production trade-off

**Goal.** Bound the problem. No intervention can exceed the theoretical yield, and the
envelope shows what growth must be given up to approach it.

**Preconditions.** Step 0 passed.

**Call.**

```python
from cmm.features import production_envelope, theoretical_yield
from cmm.visualization import production_envelope_figure, save_figure, yield_figure

result = theoretical_yield(model, PRODUCT, condition=CONDITION)
print(result.molar_yield, result.carbon_ceiling,
      result.exceeds_carbon_ceiling, result.co2_fixed)

envelope = production_envelope(model, PRODUCT, condition=CONDITION, points=20)
frame = envelope.to_frame()
```

**Outputs.** `ProductionYield`, `ProductionEnvelope`.

**Artifacts.**
`02_yield/theoretical_yield.csv`, `02_yield/production_envelope.csv`,
`figures/yield.png` (`yield_figure([result])`),
`figures/production_envelope.png` (`production_envelope_figure(envelope)`).

**Decision rule.**
- `molar_yield > 1e-6` — otherwise stop (step 0 should have caught this).
- **Read the envelope as step 3's feasibility check.** If `growth_max` falls as product rises,
  production and growth compete, which is the usual precondition for a coupled design. If
  `growth_min` is 0 across the whole range, the cell can always grow without producing, so
  coupling has to be *created* by the knockouts — possible, but a harder search. Either way
  step 3 proceeds; this reading sets the expectation.
- Record `exceeds_carbon_ceiling` and `co2_fixed` together. A yield above the substrate's
  carbon ceiling is only legitimate when CO₂ is being fixed; otherwise suspect the model.
- **`co2_fixed=True` below the carbon ceiling is not a clean bill of health.** The ceiling check
  is the only guard, so a CO₂-inflated yield that stays under the ceiling passes it silently.
  Measured on anaerobic `e_coli_core`/`EX_succ_e`: molar yield **1.3906** with CO₂ uptake open
  against **1.2000** with it closed, at a CO₂ exchange of −6.9529 and a ceiling of 1.5 — 15.9%
  high, and `exceeds_carbon_ceiling` is `False` throughout. From 0.4.0 the media presets close
  CO₂ *uptake* (secretion stays free), so this number changes; if you are reading an older run,
  check whether its yield was obtained by taking up CO₂ that a closed anaerobic fermentation
  does not supply.

**Branch.** Comparing substrates or aeration? That is more than one condition, so it is more
than one run: build one `Condition` per case, loop `theoretical_yield` over them, pass all
results to one `yield_figure`, then **pick one condition, rebuild the model state for it, and
run every remaining step under that one.** Do not carry a second condition forward implicitly.

**Failure → action.** `ValueError: no uptake capacity` means the substrate exchange is closed —
fix the medium in step 0. `no carbon uptake exchange found` means auto-detection failed; pass
`substrate=` explicitly.

**Solver.** LP.

---

## Step 2 — Reference flux state

**Goal.** Describe what the cell does *now*, and fix the wild-type numbers every later step is
compared against.

**Scope.** The design search in step 3 does **not** use a reference state — OptKnock optimizes
over the model directly. So this step is narrower than it looks. What it actually feeds:

- the wild-type product flux and growth rate that every "improvement" is measured against;
- the reference MOMA/ROOM needs in step 5a to predict the designed strain's immediate phenotype;
- the comparison point for the sampled ensemble in step 5b;
- the hand-off from `SC-02`, when the question is condition-specific.

One reference (`pfba`) is enough for all four. Build more only when you have expression data and
the condition matters — and note that extra references no longer strengthen the *design*, only
the interpretation.

**Preconditions.** Step 1 gave a non-zero yield.

**Call.**

```python
from cmm.core import supports
from cmm.features import reference_flux

                                    # inherits the condition from the model state applied in
                                    # step 0 — see the note under **Decision rule**
references = {"pfba": reference_flux(model, "pfba")}

if expression is not None:                      # gene -> value, ids matching model.genes
    references["lad"] = reference_flux(model, "lad", gene_expression=expression)
    if supports("QP", model.solver.interface):
        references["eflux2"] = reference_flux(model, "eflux2", gene_expression=expression)

wild_type_growth = references["pfba"].get(BIOMASS)
wild_type_product = references["pfba"].get(PRODUCT)
```

**Outputs.** `dict[str, FluxState]`.

**Artifacts.** `03_reference/reference_<method>.csv` per state; wild-type product flux and
growth in the report.

**Decision rule.** `pfba` is the reproducible default — a unique minimal-total-flux solution.
Use an omics-derived state instead when the run is about a specific condition, and say which
one the wild-type numbers came from, because the improvement figures are relative to it.

`reference_flux` takes no `condition=`; it reads the model exactly as step 0 constrained it.
That is the inheritance, and it only holds if nothing has changed the medium since. Record the
reference state's own fingerprint alongside it so the report can prove which condition produced
it — on anaerobic `e_coli_core` the wild-type reference gives growth 0.211663 and
`EX_succ_e` = 0.0.

**Branch.**
- No expression data → `pfba` alone. That is the normal case and costs the run nothing.
- No QP solver → LAD instead of E-Flux2, and record the substitution (`AGENTS.md` §3.3).
- Several conditions in one expression table → this is `SC-02`; come back here with the
  condition of interest.

**Failure → action.** Near-zero gene id overlap → stop, per preflight P5. A non-optimal
integration status means the expression-derived bounds are infeasible; report it rather than
falling back silently.

**Solver.** LP for pfba/lad; **QP** for eflux2.

---

## Step 3 — Knockout design

**Goal.** Find a small knockout set such that a cell maximizing its own growth **cannot avoid**
making the product, and quantify the guaranteed production rate.

This is the inverse half of the scenario: you supply the goal, the bilevel MILP searches the
knockout sets. It is also where the coupling claim is either earned or explicitly forgone.

**Preconditions.** Step 1 gave a non-zero yield. Step 0 recorded whether MILP, `straindesign`
and Java are available.

**Call.**

```python
import pandas as pd
from cmm.features import optknock, robustknock

optimistic = optknock(model, PRODUCT, condition=CONDITION,
                      max_knockouts=3, max_solutions=5, min_growth=0.05)
guaranteed = robustknock(model, PRODUCT, condition=CONDITION,
                         max_knockouts=3, max_solutions=8, min_growth=0.05)

designs = pd.DataFrame([
    {
        "knockouts": ", ".join(d.knockouts),
        "growth": d.growth,
        "max_product": d.max_product,
        "guaranteed_product": d.guaranteed_product,
        "growth_coupled": d.growth_coupled,
    }
    for d in guaranteed.designs
]).sort_values("guaranteed_product", ascending=False)
```

**Outputs.** `StrainDesignResult` with `.designs` and `.best()`; each `StrainDesign` carries
`knockouts`, `growth`, `max_product`, `guaranteed_product`, `growth_coupled`.

**Artifacts.** `04_design/optknock.csv`, `04_design/robustknock.csv`.

**Decision rule.** **Rank by `guaranteed_product`, never `max_product`.** `max_product` is what
the cell *could* make if it cooperated; `guaranteed_product` is what it makes at worst among
growth-optimal states. Only `guaranteed_product > 0` is growth-coupled, and only then may a
design be called a design rather than a possibility. `robustknock` optimizes the worst case
directly, so prefer its output when both return results.

**The number of designs is not a reproducible quantity — the top design is.** A MILP solution
pool enumerates near-optimal alternatives in an order that depends on the solver's internal
state, so running `fba`/`pfba` or `production_envelope` earlier in the same process can change
how many designs come back. The reference anaerobic run returned **18 OptKnock and 41
RobustKnock designs** (all 41 growth-coupled); another process may return a different count.
Report the ranked top designs and their guaranteed products; **do not headline the count**, and
if you mention it, say it is pool-dependent.

**The design is a property of the condition, and the condition must be the one from step 0.**
On `glucose_anaerobic` `e_coli_core`, both `optknock` and `robustknock` return the same top
design:

| Condition | Knockouts | Growth | Guaranteed product | Coupled |
|---|---|---|---|---|
| `glucose_anaerobic` | `ACALD, D_LACt2, THD2` | 0.090648 | **9.910758** | yes |
| `glucose_aerobic` | `ACALD, D_LACt2, THD2` | 0.873922 | 0.000000 | no |
| `glucose_anaerobic` | `CO2t, FORti, PGI` | 0.000000 | — | **does not grow** |
| `glucose_aerobic` | `CO2t, FORti, PGI` | 0.143322 | 10.406319 | yes |

Read the table as a warning, not as a menu. **Earlier versions of this document presented
`{CO2t, FORti, PGI}` at `guaranteed_product` 10.4063 as the anaerobic answer. It is an aerobic
result** — that deletion set does not grow anaerobically at all, so its guaranteed product is
undefined under the condition the rest of that run used. The mistake was produced exactly as
described above: an aerobic medium was applied while `aerobic=False` was passed to functions
that had the parameter, and the design search, which never had it, ran with oxygen open. The
same design ranked first under either condition and the numbers looked plausible under both,
which is why nothing caught it. **Print the medium and the oxygen bounds next to every design
table.**

The anaerobic design's gene-level deletions, for reference:
`ACALD` → `b0351 or b1241`; `D_LACt2` → `b2975 or b3603`; `THD2` → `b1602 and b1603`. An `or`
rule requires deleting **all** listed genes; an `and` rule requires deleting **any one**.

**Branch.**
- **No MILP, no `straindesign`, or no Java → run the LP fallback below.** This is the only
  supported substitution for this step.
- No design found → raise `max_knockouts` (cost grows fast) or lower `min_growth`. Changing the
  medium or the aeration is **not** a parameter tweak: it is a different experiment, so go back
  to step 0, reset the condition, and re-run every step under it. If still nothing, report that
  no coupled design exists under these constraints. That is a legitimate negative result, not a
  failed run.
- Designs found but every `guaranteed_product == 0` → report them as *uncoupled candidates* and
  carry them into step 5 as hypotheses, not as strains.
- **Map reaction knockouts back to genes** before anything experimental — a reaction is not
  deletable:
  ```python
  for rid in guaranteed.best().knockouts:
      rxn = model.reactions.get_by_id(rid)
      print(rid, sorted(g.id for g in rxn.genes), rxn.gene_reaction_rule)
  ```
  An `or` rule means every listed gene must go. Step 6 requires this gene-level list.
- Want an exhaustive single-deletion study alongside the design — essentiality classes,
  explicit thresholds, a genome-scale strategy → run
  [`SC-03`](SC-03-knockout-screening.md). It answers a different question and is a complete
  study in its own right, not a precursor to this step.

**Failure → action.** MILP timeouts are common at genome scale. Reduce `max_solutions` first,
then `max_knockouts`, and record what you reduced — a silently shrunk search is a silently
weaker claim.

**Solver.** **MILP** + `straindesign` + Java.

### Step 3 fallback — no MILP available

Run a single-deletion screen instead. It finds candidates, not designs.

```python
from cmm.features import batch_comparison, gene_perturbations

result = batch_comparison(
    model, references["pfba"], gene_perturbations(model),
    method="moma_l1",                       # LP; moma_l2 needs QP, room needs MILP
    product_reaction=PRODUCT,
)
screen = result.to_frame()                  # and save result.metadata beside the CSV
```

**Decision rule.** A candidate must satisfy all three: `status == "optimal"` (not lethal),
`objective >= 0.1 * wild_type_growth` (state your threshold — 10% is a starting point, not a
standard), and `product_flux > wild_type_product`. Rank by `product_flux` descending.

**Artifacts.** `04_design/screen_<reference>_<method>.csv`. Keep the lethal rows — they are the
essentiality result, and filtering them makes the table look cleaner than the biology is.

**Id space.** `gene_perturbations` returns *gene* ids (`b0722`) while step 4 returns *reaction*
ids (`PPC`). They never intersect, so a step 5c consensus built naively over both is silently
empty. Either screen with `reaction_perturbations(model)`, or map through the GPR before
crossing:

```python
from cmm.features._perturbation import blocked_reactions_for_genes

screen_candidates = screen.loc[                 # the DataFrame built just above
    (screen["status"] == "optimal")
    & (screen["objective"] >= 0.1 * wild_type_growth)
    & (screen["product_flux"] > wild_type_product),
    "target_id",
]
ko_reactions = {
    r for t in screen_candidates for r in blocked_reactions_for_genes(model, [t])
}
```

(`screen_candidates` is named apart from step 5b's `candidates`, which is a list of *reaction*
ids from step 4 — the whole point of this paragraph is that the two id spaces do not mix.)

**What the report must say.** This substitution changes the scientific claim, so `AGENTS.md`
§3.3 applies in full. State all three of these in **Setup** and again in **Limitations**:

1. the method substituted and why (which of MILP / `straindesign` / Java was missing);
2. that only **single** deletions were examined — no combinations were searched;
3. that **coupling was not established**, and therefore no target in this report may be
   described as growth-coupled.

**Solver.** LP (`moma_l1`). `moma_l2` is QP, `room` is MILP — if you had MILP you would be
running the design search instead.

---

## Step 4 — Amplification targets

**Goal.** Find reactions that must carry more flux for the product to increase — the
over-expression candidates — and the reactions that must carry less.

**Preconditions.** Step 1 gave a non-zero yield.

**Call.**

```python
from cmm.features import fseof, fvseof
from cmm.visualization import fseof_figure, fvseof_figure

scan = fseof(model, PRODUCT, condition=CONDITION, n_steps=10)
amplify = scan.amplification_targets()          # actionable_only=True by default
knockdown = scan.knockout_targets()

robust_scan = fvseof(model, PRODUCT, condition=CONDITION, n_steps=10, biomass_fraction=0.95)
robust = robust_scan.robust_targets()
```

**Outputs.** `FseofResult`, `FvseofResult`.

**Artifacts.** `05_amplification/fseof_trends.csv`, `05_amplification/fvseof_{mean,forced,
capacity}.csv`, `figures/fseof.png`, `figures/fvseof.png`.

**Decision rule.** Priority order:
1. `set(robust) & set(amplify)` — rising mean flux **and** rising forced minimum. The reaction
   cannot avoid carrying more flux, so it is the strongest amplification evidence CMM offers.
2. `amplify` alone — rising mean flux only; the network *may* route around it.
3. `knockdown` — feed to step 3's candidate list if not already there.

`robust_scan.amplification_targets()` is already in Park et al.'s own priority order — types
1–3 by the joint sign of ΔV_avg and Δl_sol, ordered by ascending mean `l_sol` — and
`robust_scan.park_type` gives the 1–9 index per reaction. Report that order as Park's; the
`robust` intersection in rule 1 is CMM's own extra filter on top of it, not a re-ranking by
them.

Keep `actionable_only=True`: it drops boundary, objective, and no-GPR reactions, which are not
things anyone can engineer.

**Two of these are CMM's constructs, not the source papers', and the report must say so.**
CMM's FSEOF selection rule (endpoint difference, positive linear slope, no sign reversal,
baselined at the 10% scan level) is deliberately stricter than the criterion in Choi et al.
(2010), which selects on `|v_j|max > |v_j^initial|` and `v_j^max · v_j^min ≥ 0`. On anaerobic
`e_coli_core`/succinate Choi's rule additionally admits the acetate-secretion branch
(`ACKr`, `ACt2r`, `EX_ac_e`, `PTAr`) — reactions the design search in step 3 *deletes* — so CMM
keeps its own rule. Likewise FVSEOF's `robust_targets()` (forced FVA minimum rising
monotonically) is CMM's addition; it is not the variability criterion of Park et al. (2012).
Cite Choi et al. and Park et al. for the methods, and attribute these two selection rules to
CMM.

**Branch.** Empty `robust_targets()` → fall back to `amplify`, and say the targets lack the
robustness check. Empty `amplify` too → the product does not respond to enforced flux under
this medium. That is a finding about this condition; report it as such, and if you then try a
different aeration or substrate, restart from step 0 rather than swapping it in here.

**Failure → action.** FSEOF on a zero-yield product returns nothing useful; step 1 gates this.
Genome-scale scans are slow — narrow with `reactions=` to a pathway of interest rather than
reducing `n_steps` below 10, which coarsens the trend classification and, for FVSEOF, falls
below the resolution Park et al. specify.

**Solver.** LP (FVSEOF runs an FVA per step, so it costs more wall clock, not more capability).

---

## Step 5 — Verification

**Goal.** Test what steps 3 and 4 proposed. Those steps are *inverse* — they search for
interventions. This step is *forward*: it applies each intervention and predicts what follows.

Step 3 and step 4 produce different kinds of thing, so verification runs in two tracks:

| Track | Verifies | Method | Open question it answers |
|---|---|---|---|
| **5a** | the knockout design | MOMA/ROOM, `flux_response` on the knocked-out model | what does this strain do the day it is built? |
| **5b** | the amplification targets | `flux_response`, sampling | does pushing flux here actually buy product? |
| **5c** | both together | consensus table | which interventions combine into one strain? |

**Preconditions.** A design (or fallback candidates) from step 3, amplification targets from
step 4, and a reference state from step 2.

### 5a — The design: what does the strain do immediately after deletion?

```python
from cmm.features import flux_response, knockout_comparison
from cmm.visualization import flux_comparison_figure, flux_response_figure

best = guaranteed.best()
immediate = knockout_comparison(model, references["pfba"], best.knockouts, method="moma_l2")

with model:                                   # design applied, then restored;
    for rid in best.knockouts:                # the medium from step 0 is untouched
        model.reactions.get_by_id(rid).knock_out()
    designed_scan = flux_response(model, PRODUCT, condition=CONDITION,
                                  biomass_fraction=0.0, n_steps=20)
```

**Decision rule.** MOMA assumes a *minimal adjustment* from wild type, so it describes the
freshly built strain — before any adaptation. Read it against the design's own numbers:

| MOMA result | Reading |
|---|---|
| `status != "optimal"` | The design is lethal under minimal adjustment. **Discard it.** |
| `fluxes[PRODUCT]` ≈ `guaranteed_product` | The strain works as designed from day one. |
| `fluxes[PRODUCT]` ≪ `guaranteed_product` | The design still holds; the strain needs adaptive evolution to reach it. Report the gap — it is an estimate of how much ALE is required. |

**MOMA cannot falsify the coupling claim.** Coupling was proven in step 3 by the bilevel MILP
at the growth optimum; MOMA answers a different question under a different cell model. A low
MOMA product is a schedule, not a refutation. Do not drop a coupled design on this basis.

The `designed_scan` shows the relationship the knockouts created, and its
`feasible_range()` is the crispest evidence of coupling you will get from a forward method: a
lower bound above zero means the strain **cannot** produce less than that, whatever it does. On
anaerobic `e_coli_core` the `{ACALD, D_LACt2, THD2}` design gives `(4.794286, 13.905778)` over
20 of 20 feasible scan points — zero succinate is simply not a solution any more. Compare that
lower bound against the design's own `guaranteed_product` of 9.910758: the scan bounds the
worst case over *all* feasible states, the MILP bounds it over the growth-optimal ones, so the
scan's floor is the looser of the two and both belong in the report.

**Artifacts.** `06_validation/design_moma.csv`, `06_validation/design_flux_response.csv`,
`figures/design_flux_comparison.png`, `figures/design_flux_response.png`.

**Solver.** **QP** for `moma_l2` — use `moma_l1` (LP) on an LP-only solver and say so; LP for
the scan.

### 5b — Amplification targets: does forcing flux buy product?

```python
from cmm.features import flux_response
from cmm.visualization import flux_response_figure

responses = {}
for target in candidates:                    # amplification candidates from step 4
    responses[target] = flux_response(
        model, target, response=PRODUCT,
        condition=CONDITION,                 # the same condition as every step above
        biomass_fraction=0.3,                # keep the cell viable across the scan
        n_steps=20,
    )
```

**Decision rule per target**, in this order — the first two are gates, and skipping them is how
a meaningless number becomes a recommendation:

1. **Does the response vary at all?** Compare `optimum().response_flux` against the response at
   the wild-type target flux. If the curve is essentially flat, this reaction does not drive
   the product and the "optimum" is an arbitrary point on a plateau. Drop the target.
   ```python
   feasible = response.to_frame().query("status == 'optimal'")
   spread = feasible["response_flux"].max() - feasible["response_flux"].min()
   drives_product = spread > 0.05 * feasible["response_flux"].max()
   ```
2. **Is the target physically meaningful?** Check the optimum against the reaction's own
   bounds and its sampled distribution from 5b. A reaction whose optimum sits at hundreds of
   mmol gDW⁻¹ h⁻¹ is usually part of a thermodynamically infeasible loop (`FRD7`/`SUCDi` in
   `e_coli_core` are the classic pair), not an engineering target. Loops show up in 5b as a
   huge sampled standard deviation.
3. **Which direction, and is there one?** `optimum().target_flux` versus
   `wild_type["target_flux"]`: higher means over-express, lower means knock down, and the ratio
   is roughly how much. Equal means **no intervention** — the cell is already at the best value
   for this reaction, so it is not a target however it ranked in step 4.
   **This routinely contradicts step 4.** On anaerobic `e_coli_core` succinate, **7 of the 17
   reactions FSEOF classifies as `amplify` — `FBA`, `FUM`, `GAPD`, `MDH`, `PGI`, `PGM`, `TPI` —
   have a forward-scan optimum *below* their wild-type flux**: the forward method says knock
   them down. `GLCpts`'s optimum sits exactly at its wild-type flux, so it is no intervention at
   all. FSEOF observed a correlation across enforced product levels; the response scan actually
   tested the intervention. Trust the forward method and report the disagreement rather than
   quoting whichever is more flattering.
4. `feasible_range()` bounds how far the intervention can go before the cell stops solving.
   It is computed by FVA on the scanned reaction, so it is the reaction's true range under the
   applied growth floor and not an artefact of where the grid points happened to land. (Before
   0.4.0 it was read off the scan grid and was inward-biased whenever a growth floor was
   applied — on `PGI` at the documented default `n_steps=20` it returned (−37.3684, 6.8421)
   against a true FVA range of (−38.0997, 9.9463), understating the headroom by 31%. Do not
   quote a pre-0.4.0 `feasible_range` as a bound.)
5. **The shadow price of the response with respect to the target** — `d(response)/d(target)`,
   returned exactly by the LP dual — says how much product one more unit of enforced flux buys,
   and the phase boundaries say where that rate changes. A boundary *inside* the useful range
   means pushing past it costs product, so the intervention has a ceiling worth reporting.
   This replaces the `bottleneck` field removed in 0.4.0: that field located the steepest
   finite-difference decline, which for a piecewise-linear LP response curve is an artefact of
   the grid — its reported location moved by up to 29.53 flux units as `n_steps` went 6 → 160,
   and its `found` flag inverted on `PGI`, `TPI` and `EX_o2_e`. **No published criterion defines
   a bottleneck as the argmin of a finite-difference slope**; regions of constant shadow price
   and the boundaries between them are the published objects (Edwards, Ramakrishna & Palsson
   2002). Do not cite any paper for the removed field, and do not carry its numbers forward.

**`biomass_fraction` is not optional here.** With no growth floor the scan reports what a
non-growing cell could do, which is not a strain.

#### Sampling: is the predicted flux forced, or one of many optima?

```python
from cmm.features import random_flux_sampling
from cmm.visualization import sampling_figure

ensemble = random_flux_sampling(model, n=1000, condition=CONDITION, seed=0)
statistics = ensemble.statistics()                          # method="achr" for small n
```

**Decision rule per target.** Compare the reference's predicted flux against the sampled
distribution for that reaction:
- narrow distribution (small `std` relative to `mean`) → the constraints pin this flux, and the
  prediction is solid;
- wide distribution → the reference value was one arbitrary choice among alternate optima, so a
  target ranked on that value is weak evidence. Say so in the report rather than dropping it.

`ensemble.correlation()` shows which reactions co-vary — a target strongly correlated with the
product exchange is mechanistically plausible.

Use `reference_constrained_sampling(model, references["pfba"], ...)` instead when the question
is "how much could this *prediction* vary", rather than "what can the network do at all".

Sampling also supports 5a: sample inside the same `with model:` block that applies the
knockouts, still passing `condition=CONDITION`, to see the flux space the design creates. A
knockout set that leaves no feasible space is lethal.

### 5c — Combining the two tracks into one strain proposal

**First put everything in the same id space.** Step 3's designs are *reaction* ids; step 4's
targets are *reaction* ids; but the step 3 fallback screens *gene* ids. Crossing gene ids with
reaction ids yields an all-`False` table that reads as "no method agreed" when nothing was ever
compared. Map through the GPR before crossing.

Then build one table per proposed strain, not per isolated target:

| Design | Over-express | Knock out (rxn) | Knock out (genes) | Growth | Guaranteed product | MOMA immediate | Coupled | Evidence |
|---|---|---|---|---|---|---|---|---|

- **Over-express** — step 4 targets that passed 5b's three gates.
- **Knock out** — the step 3 design, with its GPR-resolved gene list.
- **MOMA immediate** — 5a's prediction for the freshly built strain.
- **Coupled** — `guaranteed_product > 0`. On the LP fallback this column is `not established`
  for every row; never leave it blank.
- **Evidence** — which independent methods support each part.

**Decision rule.** A design supported by the inverse method (step 3 coupling proof) *and*
surviving the forward check (5a not lethal) is the strongest claim CMM supports. Amplification
targets are ranked separately by how many of {FSEOF, FVSEOF-robust, 5b response} agree. Keep the
two rankings distinct inside the row — they were produced under different cell models and a
single merged score would hide that.

**Artifacts.** `06_validation/flux_response_<target>.csv` per target,
`06_validation/sampling_statistics.csv`, `06_validation/consensus.csv`,
`figures/flux_response_<target>.png`, `figures/sampling.png`.

**Failure → action.** Infeasible response points are data — keep them. A sampler that will not
converge (`RuntimeError` from `_draw`) usually means the model is over-constrained; loosen the
condition or reduce `thinning`. If `optgp` gives unstable statistics, raise `n` above 1000 or
switch to `achr`, and record which you used.

**Solver.** LP for both.

---

## Step 6 — Report

Follow `_reporting.md`. Scenario-specific requirements:

- **The recommended table is per strain proposal**, in the 5c shape — over-expression targets
  and the knockout design in one row, with growth, guaranteed product, MOMA immediate phenotype,
  and the coupling verdict.
- **Give the gene-level deletions** implied by every reaction knockout, including where an `or`
  GPR requires deleting several genes. A reaction id is not something anyone can delete.
- **Never quote `max_product` without `guaranteed_product` beside it.**
- **State the intervention direction quantitatively** where 5b provides it ("increase PGI flux
  from 4.9 to about 6.8 mmol gDW⁻¹ h⁻¹").
- **If step 3 used the LP fallback**, Setup and Limitations must both carry the three
  disclosures listed there: which capability was missing, that only single deletions were
  examined, and that coupling was not established.
- **State the condition once, in Setup, and never restate it differently.** Give the medium
  preset name, the oxygen exchange bounds, the substrate and its uptake rate, and the model
  fingerprint taken after the medium was applied. Every number in the report comes from that
  one condition.
- **Limitations** must include: predictions are hypotheses requiring experimental validation;
  the medium, substrate, and aeration assumed; that coupling is a property of the model's
  growth-maximizing assumption while MOMA predicts the immediate post-deletion phenotype, so a
  gap between them means adaptation rather than error; and any solver capability that forced a
  substitution.

---

## Cross-checks

Before writing the report, confirm the pipeline is internally consistent:

- `guaranteed_product <= max_product` for every design. A violation is a reporting error.
- No design's product and no recommended target's predicted product exceeds the step 1
  theoretical maximum.
- Step 5b's optimum growth is consistent with the step 1 envelope at that product flux.
- A reaction appearing as both an amplification target and part of the knockout design is a
  contradiction — usually a sign the two came from different aeration or substrate settings.
  Resolve it or report it explicitly.
- The wild-type product flux is the same number everywhere it appears.
- Every knockout in the design is a real reaction id present in the model, and every one has a
  gene-level mapping in the report.
- The consensus table is not all-`False` in its cross-method columns. If it is, the tracks were
  compared in different id spaces (genes vs reactions).
- No recommended target is a known futile-cycle reaction. Check the sampled standard deviation:
  a reaction ranging over hundreds of flux units in an ensemble is a loop artifact.

## Do not

- Do not report an FSEOF rank as a recommendation without step 5. FSEOF is a fast heuristic
  with no coupling guarantee.
- Do not rank or recommend a design by `max_product`.
- Do not call a design growth-coupled without `guaranteed_product > 0` — and on the LP fallback,
  do not call anything growth-coupled at all.
- **Do not merge step 3 and step 5 results into a single score.** They assume different cell
  behaviour: step 3 assumes the cell reaches its growth-maximizing state, MOMA in 5a assumes a
  minimal adjustment from wild type. Report them side by side in the same row, never summed.
- Do not discard a coupled design because MOMA predicts low immediate product. That is an
  adaptation estimate, not a refutation.
- Do not present reaction knockouts as if they were gene deletions.
- Do not drop lethal knockouts or infeasible scan points from the exported tables.
- **Do not change the medium, the aeration or the substrate between steps.** Set the condition
  once in step 0 and pass it down. Every number in one report must come from one condition, or
  the condition must be a reported variable with its own column.
- **Do not describe a design without naming the condition it was found under.** A knockout set
  that is growth-coupled aerobically can fail to grow at all anaerobically, and the same set
  can rank first under both — see the table in step 3.
- Do not quote a `bottleneck` location from a pre-0.4.0 run, and do not attribute the removed
  bottleneck criterion to any publication.
