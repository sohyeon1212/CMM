---
id: SC-01
title: Production target discovery
goal: Find and verify genetic interventions that increase production of a target metabolite
when_to_use:
  - "increase production of X"
  - "which genes should I over-express / knock out to make more X"
  - "생산 증대 표적을 찾아줘"
role: spine                 # the backbone of a production goal; complete on its own
optional_inputs_from:
  - "SC-03: a condition-specific flux state, used as a baseline in step 2"
  - "SC-04: screened knockout candidates, used in step 3"
optional_next:
  - "SC-02: prove coupling for the candidates worth designing around"
requires:
  model: "cobra model"
  product: "exchange reaction id (e.g. EX_succ_e)"
optional:
  expression: "gene expression table (enables LAD / E-Flux2 baselines)"
  substrate: "exchange reaction id; default: auto-detected from the medium"
  aerobic: "true | false"
solver:
  minimum: "LP"
  full: "QP (E-Flux2, L2 MOMA), MILP (ROOM)"
steps: [preflight, yield, baseline, knockout, amplification, validation, report]
runtime: "minutes on a core model; hours on genome scale with many baselines"
---

# SC-01 — Production target discovery

## Objective

Produce a ranked, verified list of genetic interventions — reactions to amplify and reactions
to knock out — that increase flux to a target product, together with the predicted product and
growth for each, and evidence that each survived independent verification.

**Success criteria.**
1. Theoretical yield is known and non-zero, with its carbon balance disclosed.
2. Targets are derived against **more than one baseline flux state**, and the report
   distinguishes targets that all baselines agree on from those only one produced.
3. Every recommended target passed step 5 verification, not just step 3 or 4 ranking.
4. The run directory satisfies `_reporting.md`.

## Pipeline at a glance

| Step | Question | Method | Output |
|---|---|---|---|
| 0 | Is this model usable at all? | `_preflight.md` | viability, solver, yield |
| 1 | What is the ceiling, and what does it cost? | `theoretical_yield`, `production_envelope` | yield, trade-off curve |
| 2 | What is the cell doing now? | `reference_flux` × {pfba, lad, eflux2} | baseline flux states |
| 3 | Which knockouts help? | `batch_comparison` (MOMA/ROOM) per baseline | knockout candidates |
| 4 | Which reactions should carry more flux? | `fseof`, `fvseof` | amplification candidates |
| 5 | Do the candidates survive scrutiny? | `flux_response`, sampling, consensus | verified targets |
| 6 | Write it up | `_reporting.md` | report + figures + raw data |

Steps 3 and 4 are independent and may run in either order. Step 5 needs both.

---

## Step 0 — Preflight

**Goal.** Establish that the model grows, the product is reachable, and the solver is capable.

**Call.** Follow `_preflight.md` (P1–P4; P5 only if you have expression data).

**Branch.** Theoretical yield zero → stop and report; the product cannot be made in this
medium. Everything below assumes a non-zero ceiling.

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
- Read the envelope: if `growth_max` falls as product rises, production and growth compete and
  a growth-coupled design is plausible (→ consider `SC-02`). If `growth_min` is 0 across the
  whole range, nothing forces the cell to produce, so knockouts alone will not guarantee
  production.
- Record `exceeds_carbon_ceiling` and `co2_fixed` together. A yield above the substrate's
  carbon ceiling is only legitimate when CO₂ is being fixed; otherwise suspect the model.

**Branch.** Comparing substrates or aeration? Loop `theoretical_yield` over each and pass all
results to one `yield_figure`, then pick the condition for the rest of the run and say so.

**Failure → action.** `ValueError: no uptake capacity` means the substrate exchange is closed —
fix the medium in step 0. `no carbon uptake exchange found` means auto-detection failed; pass
`substrate=` explicitly.

**Solver.** LP.

---

## Step 2 — Baseline flux states

**Goal.** Describe what the cell does *now*. Every knockout prediction in step 3 is measured
against this, so it is the single largest source of variation in the whole pipeline — which is
exactly why you build more than one.

**Preconditions.** Step 1 gave a non-zero yield.

**Call.**

```python
from cmm.core import supports
from cmm.features import reference_flux

baselines = {"pfba": reference_flux(model, "pfba")}

if expression is not None:                      # gene -> value, ids matching model.genes
    baselines["lad"] = reference_flux(model, "lad", gene_expression=expression)
    if supports("QP", model.solver.interface):
        baselines["eflux2"] = reference_flux(model, "eflux2", gene_expression=expression)
```

**Outputs.** `dict[str, FluxState]`.

**Artifacts.** `03_baseline/reference_<method>.csv` per baseline; wild-type product flux and
growth per baseline in the report.

**Decision rule.** Keep every baseline that solved. Their disagreement is information: a
knockout target that only appears under one baseline is weaker evidence than one all of them
produce, and step 5 uses exactly that.

**Branch.**
- No expression data → `pfba` alone. Say so; the run has no condition-specific evidence.
- No QP solver → LAD only, and record the substitution (`AGENTS.md` §3.3).
- Several conditions in one expression table → this is `SC-03`; come back here with the
  condition of interest.

**Failure → action.** Near-zero gene id overlap → stop, per preflight P5. A non-optimal
integration status means the expression-derived bounds are infeasible; report it rather than
falling back silently.

**Solver.** LP for pfba/lad; **QP** for eflux2.

---

## Step 3 — Knockout targets, per baseline

**Goal.** Find single deletions that raise product flux while keeping the cell viable.

**Preconditions.** Step 2 produced at least one baseline.

**Call.**

```python
import pandas as pd
from cmm.features import batch_comparison, gene_perturbations

perturbations = gene_perturbations(model)     # or reaction_perturbations(model)
METHOD = "moma_l2"                            # moma_l1 (LP) / room (MILP) also valid

knockout_tables = {}
for name, reference in baselines.items():
    rows = batch_comparison(
        model, reference, perturbations,
        method=METHOD, product_reaction=PRODUCT,
    )
    knockout_tables[name] = pd.DataFrame([vars(r) for r in rows])
```

**Outputs.** One `DataFrame` per baseline with `target_id | kind | status | distance |
objective | n_reactions | product_flux`.

**Artifacts.** `04_knockout/batch_<baseline>_<method>.csv` per baseline.

**Decision rule.** A candidate must satisfy all three:
1. `status == "optimal"` — the knockout is not lethal;
2. `objective >= 0.1 * wild_type_growth` — the strain still grows usefully (state your
   threshold; 10% is a starting point, not a standard);
3. `product_flux > wild_type_product_flux` — production actually improves.

Then rank by `product_flux` descending, and mark how many baselines each candidate satisfies.
**Candidates found under every baseline go to step 5 first.**

**Branch.**
- `product_reaction` omitted → you only learn which knockouts hurt growth, not which help
  production. Re-run with it.
- **Id space matters for step 5c.** `gene_perturbations` returns *gene* ids (`b0722`) while
  step 4 returns *reaction* ids (`PPC`). They never intersect, so a consensus table built
  naively over both is silently empty. Either screen with `reaction_perturbations(model)` so
  both steps speak reaction ids, or keep genes here and map them before crossing:
  ```python
  from cmm.features._perturbation import blocked_reactions_for_genes
  gene_to_reactions = {
      row.target_id: blocked_reactions_for_genes(model, [row.target_id])
      for row in rows if row.kind == "gene"
  }
  ```
  Screening genes is the biologically actionable choice; screening reactions is the one that
  composes with step 4. Pick deliberately and say which.
- No candidate passes → single deletions do not help. Move to `SC-02` (multi-knockout
  growth-coupled design) or rely on the amplification targets from step 4.
- Genome-scale model → `gene_perturbations` is thousands of solves. Use `moma_l1` (LP) for the
  screen and re-check the survivors with `moma_l2`, or restrict the list to genes on the
  product's pathway.
- Want the screen read as a study rather than as one step — essentiality classes, explicit
  thresholds, a genome-scale strategy → run [`SC-04`](SC-04-knockout-screening.md) and bring
  its beneficial-deletion candidates back here. SC-04 is this step done exhaustively; it is a
  complete study in its own right, not a mandatory precursor.

**Failure → action.** Lethal knockouts return `status="infeasible"` and `essential` behavior —
this is expected and **belongs in the exported table**. Do not filter them out of the CSV; they
are the essential-gene result.

**Solver.** `moma_l1` LP, `moma_l2` **QP**, `room` **MILP**.

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

**Goal.** Separate targets that survive independent scrutiny from those that were an artifact
of one method's assumptions. Steps 3 and 4 rank; this step tests.

**Preconditions.** Candidate lists from steps 3 and 4.

### 5a — Flux response: does forcing flux through the target actually buy product?

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

### 5b — Sampling: is the predicted flux forced, or one of many optima?

```python
from cmm.features import random_flux_sampling
from cmm.visualization import sampling_figure

ensemble = random_flux_sampling(model, n=1000, seed=0)      # method="achr" for small n
statistics = ensemble.statistics()
```

**Decision rule per target.** Compare the baseline's predicted flux against the sampled
distribution for that reaction:
- narrow distribution (small `std` relative to `mean`) → the constraints pin this flux, and the
  step 2/3 prediction is solid;
- wide distribution → the baseline value was one arbitrary choice among alternate optima, so a
  target ranked on that value is weak evidence. Say so in the report rather than dropping it.

`ensemble.correlation()` shows which reactions co-vary — a target strongly correlated with the
product exchange is mechanistically plausible.

Use `reference_constrained_sampling(model, baselines["pfba"], ...)` instead when the question is
"how much could this *prediction* vary", rather than "what can the network do at all".

### 5c — Cross-method consensus

**First put every candidate in the same id space** — see step 3's branch note. Crossing gene
ids from step 3 with reaction ids from step 4 yields an all-`False` table that looks like
"no method agreed" when in fact nothing was ever compared.

Then build one table over the union:

| Target | FSEOF | FVSEOF robust | KO (pfba) | KO (lad) | KO (eflux2) | Response ↑ | Sampling |
|---|---|---|---|---|---|---|---|

**Decision rule.** Rank by number of independent methods agreeing. A target confirmed by an
inverse method (step 3/4) *and* a forward method (5a) is the strongest claim CMM supports. A
target from one method only is a hypothesis — report it as one.

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

- **Recommended targets table** separates amplification from knockout targets, and gives for
  each: predicted product flux, predicted growth, which methods agreed, and whether step 5
  verification passed.
- **Report the baselines separately** before the consensus. A reader must be able to see that
  a target came only from the E-Flux2 baseline.
- **State the intervention direction quantitatively** where 5a provides it ("increase PGI flux
  from 4.9 to about 6.8 mmol gDW⁻¹ h⁻¹").
- **Limitations** must include: predictions are hypotheses requiring experimental validation;
  the medium, substrate, and aeration assumed; that MOMA/ROOM assume minimal adjustment from
  wild type while FSEOF assumes the cell tolerates enforced product flux; and any solver
  capability that forced a substitution.

---

## Cross-checks

Before writing the report, confirm the pipeline is internally consistent:

- No recommended target's predicted product flux exceeds the step 1 theoretical maximum.
- Step 5a's optimum growth is consistent with the step 1 envelope at that product flux.
- A reaction appearing as both an amplification target and a knockout candidate is a
  contradiction — usually a sign the two came from different aeration or substrate settings.
  Resolve it or report it explicitly.
- The wild-type product flux is the same number everywhere it appears.
- The consensus table is not all-`False` in its cross-method columns. If it is, steps 3 and 4
  were compared in different id spaces (genes vs reactions).
- No recommended target is a known futile-cycle reaction. Check the 5b sampled standard
  deviation: a reaction ranging over hundreds of flux units in an ensemble is a loop artifact.

## Do not

- Do not report FSEOF rank as a recommendation without step 5. FSEOF is a fast heuristic with
  no coupling guarantee.
- Do not mix conclusions from MOMA/ROOM and OptKnock without noting they assume different cell
  behavior (minimal adjustment vs growth maximization).
- Do not drop lethal knockouts or infeasible scan points from the exported tables.
- Do not present a knockout as growth-coupled based on this scenario. Coupling is `SC-02`'s
  `guaranteed_product`.
- Do not change `aerobic` or the substrate between steps. Every number in one report must come
  from one condition, or the condition must be a reported variable.
