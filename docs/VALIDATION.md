# Scientific validation and reproducibility

This document defines what CMM 0.3.0 has been validated against, how to reproduce the
evidence, and what the tests do not establish. A green test suite supports numerical and
implementation correctness; it is not a substitute for experimental validation of a new
biological target.

CMM source code, documentation, and test data are freely available under the MIT License at
<https://github.com/jyryu3161/CMM>. For a Bioinformatics Application Note submission, archive
the exact tagged release in Zenodo or an equivalent long-term repository and add its DOI to
the README, `CITATION.cff`, and manuscript Availability and Implementation statement.

## Reproduce the publication environment

The repository tracks `uv.lock`, including hashes and resolutions for every supported Python
version. From a clean checkout:

```bash
uv sync --frozen --all-extras
QT_QPA_PLATFORM=offscreen uv run pytest -q -ra --strict-markers \
  --durations=10 --cov=cmm --cov-branch --cov-report=term-missing --cov-fail-under=80
uv run ruff check src tests
uv run ruff format --check src tests
uv run mypy src/cmm/core src/cmm/features src/cmm/omics
uvx --from cffconvert==2.0.0 cffconvert --validate
uv export --frozen --all-extras --no-emit-project -o /tmp/cmm-requirements-audit.txt
uvx --from pip-audit==2.10.1 pip-audit \
  -r /tmp/cmm-requirements-audit.txt --disable-pip --require-hashes
uv build
uvx twine check dist/*
```

CI executes the same locked environment on Linux (Python 3.10 and 3.12), Windows 3.12, and
macOS 3.12. QP/MIQP checks use small models that fit the bundled Gurobi restricted license;
they are required tests, not environment-dependent skips. Genome-scale LP validation uses
GLPK. Pytest's ten slowest durations are retained in every CI log as a lightweight runtime
regression record; solver/model/license details remain part of result provenance.

## Evidence matrix

| Area | Independent/reference evidence | Automated checks |
|---|---|---|
| FBA, pFBA, FVA | Direct COBRApy calls on `e_coli_core` and iJO1366 | `test_validation.py`, `test_genome_scale.py` |
| MOMA | Direct COBRApy L2 MOMA on an independently constructed branched network | `test_validation.py` |
| E-Flux2 | Independent two-stage QP: maximize biological objective, fix it, minimize total squared flux | `test_validation.py`, `test_expression.py` |
| Yield and media | Direct objective/bound calculations in COBRApy | `test_validation.py`, `test_production.py`, `test_media.py` |
| FSEOF/FVSEOF | Enforced product levels, biomass optimization/FVA, regression and boundary filtering | `test_production.py`; scan-resolution sensitivity in `test_scientific_sensitivity.py` |
| MTA/rMTA | Official COBRA Toolbox test topology, expected score signs, published Equation 9 | `test_revert.py`; alpha sensitivity in `test_scientific_sensitivity.py` |
| OptKnock/RobustKnock | Distinct StrainDesign module types (`OPTKNOCK` and three-level `ROBUSTKNOCK`) plus post-solve max/min product evaluation | `test_strain_design.py` |
| GUI/state | Real offscreen Qt workflows, invalid-file rejection, model reload state invalidation | `test_app_smoke.py`, `test_scenarios.py` |
| Provenance | Deterministic SHA-256 fingerprint changes with model bounds and accompanies numerical results | `test_core_primitives.py` and feature tests |

## Genome-scale check

The non-optional genome-scale fixture loads `iJO1366` through COBRApy's model repository and
asserts the imported artifact's exact dimensions (2,583 reactions and 1,367 genes) before any
comparison. CMM FBA, pFBA, and selected-reaction FVA must match direct COBRApy results; a
succinate FSEOF scan must complete without failed levels and must not return biomass or the
product exchange as an actionable intervention.

The reconstruction derives from Orth et al. (2011),
<https://doi.org/10.1038/msb.2011.65>. Exact imported dimensions include repository-level
boundary and annotation representation, so the dimension assertion also detects a changed
upstream artifact.

## Method contracts

### E-Flux2

CMM normalizes non-negative expression-derived reaction weights, constrains reaction bounds,
maximizes the model's biological objective, enforces the requested fraction of that optimum,
and minimizes the L2 norm of all fluxes. The default objective fraction is `1.0`. Without QP
support the method raises; the opt-in `allow_l1_fallback=True` result is labeled
`eflux2_l1_fallback` and must not be reported as E-Flux2.

Reference: Kim et al. (2016), <https://doi.org/10.1371/journal.pone.0157101>.

### FSEOF and FVSEOF

Scans begin at product flux observed under maximum growth and end at a selected fraction of
the theoretical maximum. FSEOF fixes product flux and maximizes biomass at each level. FVSEOF
fixes both product and an exact fraction of maximum biomass, computes reaction FVA, and reports
midpoint, forced-minimum magnitude, range width, and their slopes. Optional grouping-reaction
equalities are caller-supplied. Diagnostic tables contain all requested reactions; actionable
lists exclude boundary/objective/no-GPR reactions.

References: Choi et al. (2010), <https://doi.org/10.1128/AEM.00115-10>; Park et al. (2012),
<https://doi.org/10.1186/1752-0509-6-106>.

### OptKnock and RobustKnock

OptKnock calls StrainDesign's optimistic two-level module. RobustKnock calls its distinct
three-level module, which guards the minimum product flux among growth-optimal states. CMM
then independently fixes each design at its optimal growth and reports both maximum and
minimum product flux. A nonzero `guaranteed_product` is therefore an evaluated property, not
an alias of the optimistic value.

References: Burgard et al. (2003), <https://doi.org/10.1002/bit.10803>; Tepper and Shlomi
(2010), <https://doi.org/10.1093/bioinformatics/btp704>.

### MTA and rMTA

The optimization, transformation score, and robust Equation 9 are documented in
[design-revert-metabolism.md](design-revert-metabolism.md). The GUI uses deterministic
E-Flux2 source-state preprocessing, which differs from the original contextualization plus
sampling protocol and must be disclosed in manuscripts.

References: Yizhak et al. (2013), <https://doi.org/10.1038/ncomms3632>; Valcárcel et al.
(2019), <https://doi.org/10.1093/bioinformatics/btz231>.

## Run provenance

Numerical result metadata includes:

- deterministic model SHA-256 over reaction bounds, objective coefficients, GPR rules, and
  stoichiometry;
- model id and active solver;
- Python, CMM, COBRApy, NumPy, pandas, and SciPy versions;
- method parameters, source-state provenance, and counts of non-optimal scan levels where
  applicable.

Save result tables together with their metadata, the exact input model and omics files, the
Git commit, and the `uv.lock` used for the analysis. Model fingerprints detect scientific
model changes; they do not replace archival storage of the model file.

## Known limits

- No current automated test validates novel target predictions against new wet-lab data.
- MIQP/MILP genome-scale runtime and feasibility depend strongly on the solver license,
  tolerances, candidate set, and network compression.
- rMTA source-state preprocessing and optional FVSEOF grouping constraints must be described
  exactly; CMM does not silently claim the full preprocessing protocol of either paper.
- Floating-point optima can vary within solver tolerance. Manuscripts should report the
  solver, tolerance, model fingerprint, parameters, and acceptance tolerance.
