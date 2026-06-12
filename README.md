# CMM — Cellular Metabolic Modeling Workbench

CMM is a desktop workbench and Python library for genome-scale metabolic modeling: flux
simulation, omics integration, metabolic-engineering design, and publication-quality
visualization. It is built on [COBRApy](https://opencobra.github.io/cobrapy/) and runs every
analysis through small, solver-neutral services so the same code powers both the GUI and
scripts.

![CMM workbench](docs/images/overview.png)

## Features

- **Simulation** — FBA, parsimonious FBA (pFBA), and FVA.
- **Growth media** — preset media (glucose/acetate/glycerol, aerobic and anaerobic) applied
  with one click; easy reaction-bound editing in the model table.
- **Perturbation response** — MOMA (L1/L2) and ROOM against a reference *template* you choose
  from FBA, pFBA, LAD, or E-Flux2.
- **Omics integration** — single-state expression → flux with **E-Flux2** (Kim 2016) and
  **LAD**; multi-condition prediction with log2 fold-change comparison between conditions.
- **Production design** — theoretical yield (with CO₂-fixation disclosure), production
  envelopes, FSEOF and **FVSEOF** (variability-aware: robustly-forced amplification targets),
  and **OptKnock / RobustKnock** growth-coupled strain designs.
- **Normalization targets** — robust MTA (rMTA) and a MOMA/MTA transformation finder that
  ranks gene knockouts moving flux from one condition toward another.
- **Visualization** — Escher-layout flux maps coloured by flux, plus paper-ready production
  envelopes, FSEOF profiles, and flux comparisons (300 DPI, colour-blind-safe).

## Installation

CMM requires Python ≥ 3.10. A QP/MILP solver (Gurobi or CPLEX) is recommended — the open
GLPK solver runs FBA/pFBA/FVA, but MOMA, ROOM, revert-metabolism (QP), OptKnock (MILP), and
original-MTA (MIQP) need a commercial solver.

```bash
# core library
python -m pip install -e .

# desktop GUI (Qt + matplotlib) and strain design
python -m pip install -e ".[desktop,design]"

# development (tests)
python -m pip install -e ".[dev]"
```

Launch the GUI:

```bash
python -m cmm.app
```

The Config menu reports the active solver and warns when it cannot run the full toolkit.

## Quick start (Python API)

```python
from cobra.io import load_model
from cmm.core import fba, pfba, apply_medium
from cmm.features import theoretical_yield, optknock

model = load_model("textbook")          # e_coli_core

apply_medium(model, "glucose_aerobic")  # preset medium
print(fba(model).objective_value)       # 0.8739  (growth rate, 1/h)
print(pfba(model).fluxes["Biomass_Ecoli_core"])

# theoretical succinate yield from glucose
y = theoretical_yield(model, "EX_succ_e")
print(f"{y.molar_yield:.3f} mol/mol {y.substrate}")   # 1.638 (aerobic)

# growth-coupled succinate knockout design
designs = optknock(model, "EX_succ_e", max_knockouts=3)
print(designs.best().knockouts)
```

Omics integration and condition comparison:

```python
import pandas as pd
from cmm.omics import predict_condition_fluxes, flux_log_change

expression = pd.read_csv("expression.csv").set_index("gene")  # genes × conditions
fluxes = predict_condition_fluxes(model, expression, method="eflux2")
log_fc = flux_log_change(fluxes.fluxes("condition_A"), fluxes.fluxes("condition_B"))
```

## GUI tutorial

1. **Load a model** — start `python -m cmm.app`; the left panel lists every reaction with
   editable Lower/Upper bounds.
2. **Pick a medium** — on the *Simulation* tab choose a preset medium and click **Apply
   medium**, then **Run FBA** / **Run pFBA**. Edit any bound in the table and re-run.
3. **Compare perturbations** — on the *Comparison* tab pick a reference template
   (FBA/pFBA/LAD/E-Flux2), a reaction to knock out, and MOMA or ROOM.
4. **Design for a product** — on the *Production* tab compute theoretical yield, a production
   envelope, or FSEOF amplification targets for a target exchange.
5. **Integrate omics** — on the *Omics* tab load an expression CSV/TSV and run E-Flux2 or LAD.
6. **Rank normalization targets** — on the *Revert Metabolism* tab load source and target
   expression CSV/TSV files, then run rMTA/MTA to rank candidate gene or reaction knockouts.
7. **Visualize** — the *Flux Map* tab renders an Escher-layout map coloured by the current
   flux distribution.

| Production envelope | Escher flux map |
|---|---|
| ![envelope](docs/images/production_envelope.png) | ![flux map](docs/images/flux_map.png) |

| FSEOF targets | Multi-condition log-change |
|---|---|
| ![fseof](docs/images/fseof.png) | ![log change](docs/images/log_change.png) |

## Correctness

Every workbench service is validated against direct COBRApy computations on `e_coli_core`
(see `tests/test_validation.py`): FBA, pFBA, FVA, MOMA, theoretical yield, and media all match
the cobra reference to numerical tolerance, so GUI results equal scripted results.

```bash
python -m pytest          # full test suite
```

## License

Proprietary. See the project owner for terms.
