# CMM — Cellular Metabolic Modeling Platform

[![CI](https://github.com/jyryu3161/CMM/actions/workflows/ci.yml/badge.svg)](https://github.com/jyryu3161/CMM/actions/workflows/ci.yml)

CMM is a Python library and Qt desktop application for constraint-based metabolic modeling.
The same solver-neutral services power scripts and the GUI: FBA/pFBA/FVA, omics integration,
perturbation response, production scans, growth-coupled strain design, MTA/rMTA target
ranking, and publication figures.

![CMM platform](docs/images/overview.png)

## Availability and implementation

Source code, documentation, test data, and reproducibility workflows are freely available
at <https://github.com/jyryu3161/CMM> under the [MIT License](LICENSE). CMM is implemented in
Python and supports Python 3.10–3.12 on Linux, macOS, and Windows. It provides both a Python
API and a Qt desktop interface; installation and a complete test run require no registration.
Tagged releases and their test data will remain available for at least two years after
publication, with issue reporting through the repository's
[GitHub Issues](https://github.com/jyryu3161/CMM/issues).

Before journal submission, the exact `v0.3.0` release must additionally be archived in
Zenodo or an equivalent long-term repository and its DOI added to this section,
`CITATION.cff`, and the manuscript's Availability and Implementation statement.

## Implemented methods

- Simulation: FBA, pFBA, FVA, editable conditions, and growth-media presets.
- Perturbations: reaction/gene/multiple knockouts, L1/L2 MOMA, ROOM, and batch screens.
- Omics: LAD, strict two-stage E-Flux2, multi-condition flux prediction, and log2 changes.
- Production: theoretical yield with carbon/CO₂ disclosure, production envelopes, FSEOF,
  and FVA-based FVSEOF with optional grouping-reaction constraints.
- Strain design: distinct OptKnock and three-level RobustKnock modules through StrainDesign,
  followed by independent maximum/guaranteed-product evaluation.
- Normalization: published MTA MIQP, published rMTA best/MOMA/worst scoring, and an explicitly
  labeled legacy continuous heuristic.
- Auditability: deterministic model fingerprints and solver/package/parameter provenance on
  numerical results.

## Reproducible installation

CMM requires Python 3.10–3.12. The publication environment is locked in `uv.lock`:

```bash
git clone https://github.com/jyryu3161/CMM.git
cd CMM
uv sync --frozen --all-extras
uv run python -m cmm.app
```

For a conventional editable installation:

```bash
python -m pip install -e ".[desktop,design,solver-gurobi]"
python -m cmm.app
```

The cross-platform installers create `.venv` and install the desktop, strain-design, and
Gurobi extras by default:

```bash
./install.sh                 # macOS / Linux / WSL
# .\install.ps1             # Windows PowerShell
./install.sh --dev           # also install test, coverage, lint, and type-check tools
./install.sh --no-gurobi     # GLPK LP/MILP only
```

Tagged wheels and source archives are published on the
[GitHub Releases page](https://github.com/jyryu3161/CMM/releases). The current source version
is 0.3.0.

## Solver requirements

Solver capability is checked before a solve; a method never silently changes formulation.

| Class | Methods |
|---|---|
| LP | FBA, pFBA, FVA, LAD, yield, envelope, FSEOF, FVSEOF |
| MILP | ROOM, OptKnock, RobustKnock |
| QP | L2 MOMA, E-Flux2, `rmta_continuous` |
| MIQP | published MTA and rMTA |

GLPK supports LP/MILP. Gurobi and CPLEX support the full table, subject to their licenses and
model-size limits. The free restricted Gurobi license is suitable for CMM's small QP/MIQP
validation models; genome-scale mixed-integer studies generally require an appropriate
full solver license.

## Python quick start

```python
from cobra.io import load_model
from cmm.core import apply_medium, fba, pfba
from cmm.features import fseof, theoretical_yield

model = load_model("textbook")
apply_medium(model, "glucose_aerobic")

growth = fba(model)
minimal = pfba(model)
yield_result = theoretical_yield(model, "EX_succ_e")
scan = fseof(model, "EX_succ_e", n_steps=8, aerobic=False)

print(growth.objective_value, minimal.status)
print(yield_result.molar_yield, yield_result.metadata["model_sha256"])
print(scan.amplification_targets())
```

Expression integration:

```python
import pandas as pd
from cmm.omics import flux_log_change, predict_condition_fluxes

expression = pd.read_csv("expression.csv").set_index("gene")
predicted = predict_condition_fluxes(model, expression, method="eflux2")
change = flux_log_change(
    predicted.fluxes("condition_A"),
    predicted.fluxes("condition_B"),
)
```

## Validation

The suite includes direct COBRApy cross-checks, a non-optional iJO1366 genome-scale test,
the official COBRA Toolbox MTA test topology, E-Flux2's independent two-stage QP, target
sensitivity checks, offscreen GUI workflows, static analysis, coverage, and distribution
build verification.

```bash
uv sync --frozen --all-extras
QT_QPA_PLATFORM=offscreen uv run pytest -q -ra --strict-markers \
  --durations=10 --cov=cmm --cov-branch --cov-report=term-missing --cov-fail-under=80
uv run ruff check src tests
uv run ruff format --check src tests
uv run mypy src/cmm/core src/cmm/features src/cmm/omics
uvx --from cffconvert==2.0.0 cffconvert --validate
uv build && uvx twine check dist/*
```

The exact evidence, method contracts, references, provenance schema, and limitations are in
[Scientific validation and reproducibility](docs/VALIDATION.md). Passing these checks
supports implementation correctness; it does not constitute wet-lab validation of a new
biological prediction.

## Documentation

- [Desktop and Python tutorial](docs/TUTORIAL.md)
- [Scientific validation and reproducibility](docs/VALIDATION.md)
- [MTA/rMTA design and equations](docs/design-revert-metabolism.md)
- [Architecture and solver contracts](docs/architecture.md)
- [Release changes](CHANGELOG.md)

## Citation and license

CMM is open source under the [MIT License](LICENSE). Citation metadata are machine-readable
in [CITATION.cff](CITATION.cff); replace the contributor placeholder with the final manuscript
authors and add the archived release DOI before submission. A manuscript should also cite the
original papers for every method it uses, listed in `docs/VALIDATION.md`.

## Release process

Every push and pull request installs the frozen lockfile and runs the cross-platform quality
gates. A tag must exactly match `pyproject.toml`; the release workflow reruns all checks,
builds the wheel and sdist, validates them, installs the wheel in a clean environment, and
then attaches the artifacts to a GitHub Release.

```bash
git tag v0.3.0
git push origin v0.3.0
```
