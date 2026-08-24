# CMM — agent operating instructions

CMM is a constraint-based metabolic modeling platform (Python library + Qt desktop app) that
runs on any `cobra` model. **Its primary purpose is metabolic engineering**: finding and
verifying genetic interventions — over-expression, knockdown, knockout — that increase
production of a target metabolite, plus the simulation, omics-integration, and
perturbation-response analyses that support that.

**Drive CMM through the Python API or its thin CLI, not the GUI.** Every numerical analysis is
a solver-neutral service in `cmm.core`, `cmm.features`, `cmm.omics`; the desktop app is a thin
view over the same calls. A complete production-engineering request uses the repository skill
`.agents/skills/cmm-production-engineering/` and the composed workflow described below. A
single-analysis request still calls its documented service directly.

Two document layers sit under this one. Read them on demand, not up front:

- **`docs/scenarios/`** — step-by-step pipelines. Start here when the user has a goal.
- **`docs/agent-reference.md`** — signatures and result objects. Read the section for the
  function you are about to call.

---

## 1. Scenario router

Match the user's goal to a scenario and follow that file. Every scenario begins with
`docs/scenarios/_preflight.md` and ends with `docs/scenarios/_reporting.md`.

**Each scenario answers its own question and finishes on its own.** None is a mandatory
prerequisite for another. They also combine, and the combinations below are the useful ones —
but only run a second scenario when the user's goal actually needs it, and say why.

| The user wants | Scenario | Ends here? |
|---|---|---|
| Increase production of metabolite X; find over-expression and knockout targets | [`SC-01`](docs/scenarios/SC-01-production-target-discovery.md) | yes — this is the spine of a production goal |
| A strain where production is *guaranteed*, not merely possible | [`SC-01`](docs/scenarios/SC-01-production-target-discovery.md), entering at **step 4** | yes — step 4 searches and evaluates growth coupling |
| Explain a metabolic difference between conditions or strains | [`SC-02`](docs/scenarios/SC-02-omics-context-engineering.md) | yes — a complete study on its own |
| Which genes are essential; a single-deletion study | [`SC-03`](docs/scenarios/SC-03-knockout-screening.md) | yes — a complete study on its own |
| Screen deletions *for a production goal* | [`SC-03`](docs/scenarios/SC-03-knockout-screening.md) → `SC-01` | no — the screen feeds SC-01's candidates |
| Find targets *in a specific condition* backed by expression data | [`SC-02`](docs/scenarios/SC-02-omics-context-engineering.md) → `SC-01` | no — SC-02 picks the condition, SC-01 searches it |

Korean phrasings map the same way: 생산 증대·과발현/녹아웃 표적 → SC-01; 성장 공역·균주 설계
→ SC-01 (4단계부터); 오믹스·조건 비교·발현 데이터 → SC-02; 넉아웃 스크리닝·필수 유전자 → SC-03.

How they relate:

- **SC-01 is the spine of a production goal.** Start here when the goal is "make more of X"
  and nothing narrower is being asked. Step 3 is an exhaustive forward single-gene screen with
  MOMA/ROOM; step 4 separately searches multi-reaction growth-coupled designs with
  OptKnock/RobustKnock and evaluates `guaranteed_product`. MOMA/ROOM do not validate those
  multi-knockout designs. A request that is only about coupling enters at step 4 and skips the
  single-knockout and amplification stages.
- **SC-02 and SC-03 are complete studies in their own right**, not sub-steps of SC-01. Each
  also composes with SC-01 when the user's goal is production: SC-02 supplies the condition to
  search in, SC-03 supplies an exhaustive single-deletion picture. Neither is subordinate to it.

`docs/scenarios/README.md` holds the same map with a diagram; read it when a request spans
more than one scenario.

If the request is a single analysis rather than a goal ("run FBA", "what is the theoretical
yield"), skip the scenarios and call the function from `docs/agent-reference.md` directly.
Two capabilities have no scenario and are reached that way: `revert_targets` (MTA/rMTA) and
`transformation_targets`, documented in `docs/agent-reference.md` §9.

### Canonical production-workflow boundary

For an end-to-end SC-01 request, do not hand-assemble a one-off script. Resolve the model,
product exchange, and condition first. When that required clarification, obtain confirmation of
the resolved run definition before using one of these equivalent public boundaries:

```bash
cmm production-targets --config CONFIG
cmm report render RUN_DIR
cmm report validate RUN_DIR
```

```python
from cmm.workflows.production import (
    ProductionWorkflowConfig,
    ProductionWorkflowResult,
    run_production_target_discovery,
)
from cmm.reporting import render_production_report, validate_production_run
```

`run_production_target_discovery(config)` owns the scientific sequence and artifact schema;
`render_production_report(run_dir, renderer="nature-r")` owns publication rendering; and
`validate_production_run(run_dir)` is the completion gate. Verify the imports or CLI commands
exist in the checked-out version before describing this composed workflow as shipped. The
individual functions below remain the public building blocks and the correct entry points for
narrow requests.

### Goal → function, at a glance

| Question | Function | Module |
|---|---|---|
| Baseline growth and flux distribution | `fba`, `pfba` | `cmm.core` |
| Flux ranges at (near-)optimal growth | `fva` | `cmm.core` |
| Ceiling on product per substrate | `theoretical_yield` | `cmm.features` |
| Growth-vs-production trade-off | `production_envelope` | `cmm.features` |
| **Rank over/under-expression targets** | `fseof`, `fvseof` | `cmm.features` |
| **Knockout sets that couple product to growth** | `optknock`, `robustknock` | `cmm.features` |
| Flux state after a specific knockout | `knockout_comparison` (MOMA/ROOM) | `cmm.features` |
| Screen many knockouts at once | `batch_comparison` | `cmm.features` |
| **Does forcing flux through a target buy product?** | `flux_response` | `cmm.features` |
| **Is a predicted flux forced, or one of many optima?** | `random_flux_sampling` | `cmm.features` |
| Expression table → flux state | `integrate_expression` (E-Flux2/LAD) | `cmm.omics` |
| Per-condition fluxes from one table | `predict_condition_fluxes` | `cmm.omics` |
| Knockouts that revert disease→healthy | `revert_targets` (MTA/rMTA) | `cmm.features` |
| Knockouts that move state A→B | `transformation_targets` | `cmm.features` |

Keep the two families straight:

- **Forward — predict a result.** MOMA/ROOM/batch, `flux_response`, sampling: *you* supply the
  intervention, CMM predicts the consequence.
- **Inverse — find the intervention.** FSEOF/FVSEOF, OptKnock/RobustKnock, revert/transform:
  *you* supply the goal, CMM proposes targets.

A complete engineering answer uses both: inverse methods to propose, forward methods to verify.

---

## 2. Solver gate — check before running, not after failing

The cobra default (GLPK) is **LP + MILP only**. Several methods need QP or MIQP and will fail
deep inside a solve otherwise.

| Method | Needs | GLPK? |
|---|---|:---:|
| FBA, pFBA, FVA, theoretical yield, envelope, FSEOF, FVSEOF, LAD, **flux response**, **sampling** | LP | ✅ |
| MOMA (L1), ROOM | LP / MILP | ✅ |
| MOMA (L2), E-Flux2, transform `moma`, `rmta_continuous` | **QP** | ❌ |
| OptKnock, RobustKnock | **MILP** + importable `straindesign` | ⚠️ |
| revert `rmta` / `mta`, transform `mta` | **MIQP** | ❌ |

```python
from cmm.core import solver_status, supports
status = solver_status(model)
supports("QP", model.solver.interface)
```

QP → install `osqp`, gurobi, or cplex. MIQP → gurobi or cplex.

If a capability is missing: switch to an LP-capable equivalent (MOMA L1 instead of L2, LAD
instead of E-Flux2) **and say in the report that you substituted and why**, or tell the user
which solver to install. Never silently produce nothing, and never silently downgrade.

The canonical SC-01 workflow is stricter: its single-knockout comparison requires MOMA-L2 and
ROOM so the two requested methods stay comparable. It fails its capability gate when QP or
MILP is missing; it does not silently replace MOMA-L2 with MOMA-L1. OptKnock/RobustKnock also
require importable `straindesign`; require Java only if the selected backend reports that it
needs it. The `nature-r` report renderer separately requires `Rscript` and loadable renderer
packages; compatible minima come from package metadata, while exact versions come from
`renv.lock` and are asserted in CI.

---

## 3. Rules

1. **Resolve the run definition.** A production solve needs an exact model path, product
   exchange reaction, and one explicit condition containing medium, substrate uptake,
   oxygen/aeration bounds, and any other changed bounds. Preflight computes model id and
   fingerprint after loading; do not ask the user to supply a hash. Inspect the model and local
   capabilities read-only before asking anything, and resolve unique facts yourself. Proceed
   when the inputs are already clear. If a missing or ambiguous user decision would change the
   scientific answer, ask only that decision, include a clearly labeled recommended option with
   its evidence and scientific consequence, and resolve dependent questions in order. If any
   clarification was required, summarize the final run definition and obtain explicit approval
   before starting the workflow. Do not add that confirmation round to an initially complete
   request. The production skill's
   [`clarification-interview.md`](.agents/skills/cmm-production-engineering/references/clarification-interview.md)
   defines the detailed protocol.
2. **Preflight first.** Run `docs/scenarios/_preflight.md` before any scenario. A model that
   does not grow, has no exchanges, or whose gene ids do not match the expression table will
   produce confident nonsense.
3. **Check solver capability before the call**, per §2.
4. **Never substitute a method silently.** If you run LAD because E-Flux2 needs QP, that
   belongs in the report, not just in your reasoning.
5. **Every reported number carries its provenance.** `result.metadata` holds the model
   fingerprint, solver, package versions, and parameters. Save it alongside the numbers.
6. **Never present a planned feature as shipped.** Check `cmm.features.PLANNED_FEATURES`.
7. **A lethal knockout is not a bug.** Removing an essential reaction makes the solve
   infeasible; that is reported as `status="infeasible"` / `essential=yes`, and an infeasible
   scan point is data. Do not treat it as an error or drop it from a table.
8. **Make strain design deterministic and rank the right quantity.** Pass an explicit
   `strain_design_seed` (default `0`) through the canonical workflow to both OptKnock and
   RobustKnock, and record it in config/provenance. Never allow the backend to invent a hidden
   random seed. Rank designs by `guaranteed_product`, not `max_product` (`SC-01` strain-design
   step).
9. **State the assumptions that change the answer**: which reference flux state, aerobic or
   anaerobic, which substrate, which solver. Each of these moves results materially.
10. **Do not over-interpret a zero-flux knockout.** Knocking out a reaction carrying no flux in
   the reference changes nothing; MOMA returns the reference. That is not evidence of safety.
11. **Verify before recommending.** A target from FSEOF or rMTA is a hypothesis. Run the
    method-appropriate verification step of the scenario (`flux_response`, sampling, and
    loop diagnostics where applicable) before presenting it as a recommendation. In SC-01,
    the canonical single-knockout candidate universe is the union of the MOMA and ROOM
    display-ranked D1–D5 rows after deduplicating equivalent blocked-reaction signatures.
    Run an equivalent model phenotype once, but retain every represented gene id in the index.
    Every unique candidate receives matched wild-type/knockout random sampling, whether or not
    it later qualifies as beneficial or recommended. For flux response, retain the pre-deletion
    wild-type background. For exactly one blocked reaction with nonzero reference flux, scan
    reference↔zero; when reference flux is already zero, scan the full feasible reaction domain
    and label it exploratory rather than causal support for deletion. A multi-reaction blocked
    signature cannot be represented on one x-axis and must remain explicit unavailable/skipped
    with a reason; never silently choose a reaction. Complete-knockout effects come from
    MOMA/ROOM plus paired sampling.
    FSEOF and FVSEOF keep independent top-10 rankings, and every candidate in their union
    receives a target-to-product flux-response record; membership in both methods is useful
    provenance, not a prerequisite for validation or promotion. Loop-flagged or unresolved
    amplification candidates are still scanned and retain their diagnostic/eligibility status,
    although they cannot become supported recommendations. A lethal, unavailable, failed, or
    otherwise non-runnable case remains as an explicit skipped/failed/status row in the
    appropriate validation index instead of disappearing from coverage. In this workflow,
    `max_flux_response_targets` is a preflight capacity guard for the whole candidate universe,
    never a runtime top-*N* selector; reject an undersized config rather than slicing targets.
    Figure 5 uses the standard flux-response axes for every completed scan: enforced candidate-
    reaction flux on x (`target_flux`) and target-product flux on y (`response_flux`).
    Amplification is a wild-type candidate→product scan. A single-reaction knockout candidate
    uses the pre-deletion wild type: reference↔zero when its reference is nonzero, otherwise the
    full feasible domain as an exploratory response. Growth is the configured minimum-growth
    constraint and secondary `biomass_flux` output, not a Figure 5 axis.
    Preserve `recommendations.csv` as a machine-readable validation artifact, but the canonical
    publication report must not synthesize it into recommended targets, a strain proposal,
    summary promotion, or a figure category. Present each method's results separately and leave
    intervention selection to the user.

---

## 4. Run contract

Every scenario writes one self-contained directory. Full specification in
`docs/scenarios/_reporting.md`; the shape is:

```
<run directory>/            explicit ProductionWorkflowConfig.output_dir — see _reporting.md
  00_config.json          resolved workflow inputs
  00_provenance.json      model fingerprint, solver, versions, every parameter
  00_summary.json         headline results and cross-checks
  00_manifest.json        authoritative artifact inventory
  model/<model-id>.xml    byte-for-byte source model; conditioned SBML sits beside it
  01_preflight/…          model and capability checks
  02_yield/…              yield and production envelope
  03_reference/…          wild-type/reference flux state
  04_single_knockout/…    separate MOMA and ROOM screens
  05_strain_design/…      OptKnock and RobustKnock results
  06_amplification/…      FSEOF and FVSEOF results
  07_validation/…         loop diagnostic, flux response, knockout sampling, recommendations
  scripts/                resolved config plus reproduce/render/validate entry points
  figures/                300 DPI PNG plus editable PDF/SVG
  report.html             the narrative, with figures placed inline
  report_standalone.html  the same, figures embedded, for sending to someone
```

Units are CMM's throughout: fluxes in mmol gDW⁻¹ h⁻¹, growth in h⁻¹, molar yield in mol/mol.
Raw tables come from result-object exports. The canonical manuscript report is rendered from
those CSVs through R with `renderer="nature-r"`; never hand-transcribe a number or draw a
figure from data that are absent from the run directory. Completion requires a clean
`validate_production_run` result, not merely an HTML file that opens.

---

## 5. Stop and ask

Inspect discoverable facts before asking. When a user decision is still required, give compatible
choices and identify the recommended choice, why it is recommended, and how it changes the
analysis or claim. Ask the user rather than guessing when:

- **No target metabolite is named** for a production goal — everything downstream depends on it.
- **The medium, substrate uptake, or oxygen/aeration bounds are not explicit.** Ask for one
  coherent condition rather than inheriting biologically decisive defaults without confirmation.
  `model-as-loaded` is a valid condition only after the user explicitly accepts the inspected
  bounds; a generic request to "use defaults" does not supply a biological condition.
- **The product has no exchange reaction**, or `model.exchanges` is empty. Production design is
  unavailable; do not silently substitute an internal reaction.
- **Expression gene ids do not overlap `model.genes`** (preflight reports the overlap). Ask for
  the id mapping instead of integrating a table that maps to nothing.
- **A required solver capability is missing** and the LP-capable substitute would change the
  scientific claim (e.g. the user explicitly asked for published rMTA, which needs MIQP).
- **Theoretical yield is zero** for the requested product — the product is unreachable in this
  medium. Report that and ask whether to change medium, substrate, or aeration.
- **The user asks for a wet-lab claim.** CMM's tests establish implementation correctness, not
  biological validity. Predictions are hypotheses to test experimentally.

---

## 6. Repository conventions

If you are modifying CMM itself rather than using it:

- `docs/clean-room-policy.md` governs what may be brought in from other codebases: implement
  from behavior and public documentation, never copy source, fixtures, or UI forms.
- Quality gate, all of which must pass:
  ```bash
  QT_QPA_PLATFORM=offscreen uv run --frozen --all-extras pytest -q --cov=cmm --cov-branch --cov-fail-under=80
  uv run --frozen --all-extras ruff check src tests
  uv run --frozen --all-extras ruff format --check src tests
  uv run --frozen --all-extras mypy src/cmm/core src/cmm/features src/cmm/omics src/cmm/workflows src/cmm/reporting
  ```
- New analyses are solver-neutral services in `cmm.features` / `cmm.omics` returning frozen
  dataclasses with `to_frame()` and `run_provenance` metadata. The GUI stays a thin view.
- A new shipped feature moves from `PLANNED_FEATURES` to `INCLUDED_FEATURES` and gains a
  method contract in `docs/VALIDATION.md`.
