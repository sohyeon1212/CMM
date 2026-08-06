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
  aerobic: "true | false"
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

**Goal.** Establish that the model grows, the product is reachable, and the solver is capable.

**Call.** Follow `_preflight.md` (P1–P4; P5 only if you have expression data), then check what
step 3 will be able to do:

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

PRODUCT, AEROBIC = "EX_succ_e", False

result = theoretical_yield(model, PRODUCT, aerobic=AEROBIC)
print(result.molar_yield, result.carbon_ceiling,
      result.exceeds_carbon_ceiling, result.co2_fixed)

envelope = production_envelope(model, PRODUCT, aerobic=AEROBIC, points=20)
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

**Branch.** Comparing substrates or aeration? Loop `theoretical_yield` over each and pass all
results to one `yield_figure`, then pick the condition for the rest of the run and say so.

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

optimistic = optknock(model, PRODUCT, max_knockouts=3, max_solutions=5, min_growth=0.05)
guaranteed = robustknock(model, PRODUCT, max_knockouts=3, max_solutions=8, min_growth=0.05)

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
how many designs come back. On anaerobic `e_coli_core` succinate the same run returns 32, 39 or
57 designs depending on what preceded it, while the best design is `{CO2t, FORti, PGI}` at
`guaranteed_product` 10.4063 every single time. Report the ranked top designs and their
guaranteed products; **do not headline the count**, and if you mention it, say it is
pool-dependent.

**Branch.**
- **No MILP, no `straindesign`, or no Java → run the LP fallback below.** This is the only
  supported substitution for this step.
- No design found → raise `max_knockouts` (cost grows fast), lower `min_growth`, or revisit the
  medium and aeration. If still nothing, report that no coupled design exists under these
  constraints. That is a legitimate negative result, not a failed run.
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

rows = batch_comparison(
    model, references["pfba"], gene_perturbations(model),
    method="moma_l1",                       # LP; moma_l2 needs QP, room needs MILP
    product_reaction=PRODUCT,
)
screen = pd.DataFrame([vars(r) for r in rows])
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
ko_reactions = {
    r for t in candidates["target_id"] for r in blocked_reactions_for_genes(model, [t])
}
```

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

scan = fseof(model, PRODUCT, n_steps=10, aerobic=AEROBIC)
amplify = scan.amplification_targets()          # actionable_only=True by default
knockdown = scan.knockout_targets()

robust_scan = fvseof(model, PRODUCT, n_steps=8, biomass_fraction=0.95, aerobic=AEROBIC)
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

Keep `actionable_only=True`: it drops boundary, objective, and no-GPR reactions, which are not
things anyone can engineer.

**Branch.** Empty `robust_targets()` → fall back to `amplify`, and say the targets lack the
robustness check. Empty `amplify` too → the product does not respond to enforced flux under
this medium; revisit aeration and substrate in step 1.

**Failure → action.** FSEOF on a zero-yield product returns nothing useful; step 1 gates this.
Genome-scale scans are slow — narrow with `reactions=` to a pathway of interest rather than
reducing `n_steps` below ~8, which coarsens the trend classification.

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

with model:                                   # design applied, then restored
    for rid in best.knockouts:
        model.reactions.get_by_id(rid).knock_out()
    designed_scan = flux_response(model, PRODUCT, biomass_fraction=0.0, n_steps=20)
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
anaerobic `e_coli_core` the `{CO2t, FORti, PGI}` design gives `(3.47, 15.0)` — zero succinate is
simply not a solution any more.

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
   **This routinely contradicts step 4.** On anaerobic `e_coli_core` succinate, FSEOF classifies
   `ADK1` as `amplify`, but its response optimum sits exactly at the wild-type flux: nothing to
   amplify. FSEOF observed a correlation across enforced product levels; the response scan
   actually tested the intervention. Trust the forward method and report the disagreement
   rather than quoting whichever is more flattering.
4. `feasible_range()` bounds how far the intervention can go before the cell stops solving.
5. `bottleneck.found` — a bottleneck *inside* the useful range means pushing past it costs
   product, so the intervention has a ceiling worth reporting.

**`biomass_fraction` is not optional here.** With no growth floor the scan reports what a
non-growing cell could do, which is not a strain.

#### Sampling: is the predicted flux forced, or one of many optima?

```python
from cmm.features import random_flux_sampling
from cmm.visualization import sampling_figure

ensemble = random_flux_sampling(model, n=1000, seed=0)      # method="achr" for small n
statistics = ensemble.statistics()
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

Sampling also supports 5a: `random_flux_sampling(model, condition=<knockouts>)` shows the flux
space the design creates, and a knockout set that leaves no feasible space is lethal.

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
- Do not change `aerobic` or the substrate between steps. Every number in one report must come
  from one condition, or the condition must be a reported variable.
