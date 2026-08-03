# CMM — agent operating instructions

CMM is a constraint-based metabolic modeling platform (Python library + Qt desktop app) that
runs on any `cobra` model. **Its primary purpose is metabolic engineering**: finding and
verifying genetic interventions — over-expression, knockdown, knockout — that increase
production of a target metabolite, plus the simulation, omics-integration, and
perturbation-response analyses that support that.

**Drive CMM through the Python API, not the GUI.** Every analysis is a solver-neutral service
in `cmm.core`, `cmm.features`, `cmm.omics`; the desktop app is a thin view over the same
calls. Scripts are what make a result reproducible and reportable.

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
| A strain where production is *guaranteed*, not merely possible | [`SC-01`](docs/scenarios/SC-01-production-target-discovery.md), entering at **step 3** | yes — step 3 designs it, step 5a checks it |
| Explain a metabolic difference between conditions or strains | [`SC-02`](docs/scenarios/SC-02-omics-context-engineering.md) | yes — a complete study on its own |
| Which genes are essential; a single-deletion study | [`SC-03`](docs/scenarios/SC-03-knockout-screening.md) | yes — a complete study on its own |
| Screen deletions *for a production goal* | [`SC-03`](docs/scenarios/SC-03-knockout-screening.md) → `SC-01` | no — the screen feeds SC-01's candidates |
| Find targets *in a specific condition* backed by expression data | [`SC-02`](docs/scenarios/SC-02-omics-context-engineering.md) → `SC-01` | no — SC-02 picks the condition, SC-01 searches it |

Korean phrasings map the same way: 생산 증대·과발현/녹아웃 표적 → SC-01; 성장 공역·균주 설계
→ SC-01 (3단계부터); 오믹스·조건 비교·발현 데이터 → SC-02; 넉아웃 스크리닝·필수 유전자 → SC-03.

How they relate:

- **SC-01 is the spine of a production goal.** Start here when the goal is "make more of X"
  and nothing narrower is being asked. Its step 3 is the growth-coupled design search
  (OptKnock/RobustKnock) and its step 5a checks that design with MOMA/ROOM — one runs the
  inverse direction, the other the forward. A request that is only about coupling enters at
  step 3 and skips the amplification half.
- **SC-02 and SC-03 are complete studies in their own right**, not sub-steps of SC-01. Each
  also composes with SC-01 when the user's goal is production: SC-02 supplies the condition to
  search in, SC-03 supplies an exhaustive single-deletion picture. Neither is subordinate to it.

`docs/scenarios/README.md` holds the same map with a diagram; read it when a request spans
more than one scenario.

If the request is a single analysis rather than a goal ("run FBA", "what is the theoretical
yield"), skip the scenarios and call the function from `docs/agent-reference.md` directly.
Two capabilities have no scenario and are reached that way: `revert_targets` (MTA/rMTA) and
`transformation_targets`, documented in `docs/agent-reference.md` §9.

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
| OptKnock, RobustKnock | **MILP** + `straindesign` + Java | ⚠️ |
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

---

## 3. Rules

1. **Preflight first.** Run `docs/scenarios/_preflight.md` before any scenario. A model that
   does not grow, has no exchanges, or whose gene ids do not match the expression table will
   produce confident nonsense.
2. **Check solver capability before the call**, per §2.
3. **Never substitute a method silently.** If you run LAD because E-Flux2 needs QP, that
   belongs in the report, not just in your reasoning.
4. **Every reported number carries its provenance.** `result.metadata` holds the model
   fingerprint, solver, package versions, and parameters. Save it alongside the numbers.
5. **Never present a planned feature as shipped.** Check `cmm.features.PLANNED_FEATURES`.
6. **A lethal knockout is not a bug.** Removing an essential reaction makes the solve
   infeasible; that is reported as `status="infeasible"` / `essential=yes`, and an infeasible
   scan point is data. Do not treat it as an error or drop it from a table.
7. **Rank strain designs by `guaranteed_product`, not `max_product`** (`SC-01` step 3).
8. **State the assumptions that change the answer**: which reference flux state, aerobic or
   anaerobic, which substrate, which solver. Each of these moves results materially.
9. **Do not over-interpret a zero-flux knockout.** Knocking out a reaction carrying no flux in
   the reference changes nothing; MOMA returns the reference. That is not evidence of safety.
10. **Verify before recommending.** A target from FSEOF or rMTA is a hypothesis. Run the
    verification step of the scenario (`flux_response`, sampling, cross-method agreement)
    before presenting it as a recommendation.

---

## 4. Run contract

Every scenario writes one self-contained directory. Full specification in
`docs/scenarios/_reporting.md`; the shape is:

```
results/<SC-id>_<model>_<timestamp>/
  00_provenance.json      model fingerprint, solver, versions, every parameter
  01_<step>/…             raw CSV per step, one table per analysis
  figures/                300 DPI PNG (+ PDF/SVG when asked)
  report.html             the narrative, with figures placed inline
```

Units are CMM's throughout: fluxes in mmol gDW⁻¹ h⁻¹, growth in h⁻¹, molar yield in mol/mol.
Raw tables come from `result.to_frame().to_csv(...)`; figures from
`cmm.visualization` + `save_figure(fig, path)` at 300 DPI. Never hand-transcribe numbers into
the report that are not also in a CSV.

---

## 5. Stop and ask

Ask the user rather than guessing when:

- **No target metabolite is named** for a production goal — everything downstream depends on it.
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
  uv run --frozen --all-extras mypy src/cmm/core src/cmm/features src/cmm/omics
  ```
- New analyses are solver-neutral services in `cmm.features` / `cmm.omics` returning frozen
  dataclasses with `to_frame()` and `run_provenance` metadata. The GUI stays a thin view.
- A new shipped feature moves from `PLANNED_FEATURES` to `INCLUDED_FEATURES` and gains a
  method contract in `docs/VALIDATION.md`.
