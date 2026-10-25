"""Post-install smoke check for an installed CMM (installer or wheel), not the source tree.

Run with the interpreter of the environment under test:

    <venv python> -I .github/scripts/install_smoke.py

It imports CMM, runs the console script, and performs real optimisations on cobra's bundled
``textbook`` (e_coli_core) model through CMM's public services, checking the results against
the known aerobic optimum. Every solve uses GLPK so the check is independent of any licence;
when gurobipy is importable, its FBA optimum is cross-checked against GLPK's.
"""

from __future__ import annotations

import csv
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import cobra
from cobra.io import load_model, read_sbml_model, write_sbml_model

import cmm
from cmm.core import fba, fva, pfba, solver_status

BIOMASS = "Biomass_Ecoli_core"
TEXTBOOK_GROWTH = 0.8739215069684  # cobra's documented aerobic textbook optimum
TOL = 1e-6


def check(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"FAIL: {message}")
    print(f"ok: {message}")


def main() -> None:
    print(f"python {sys.version.split()[0]} on {sys.platform}")
    print(f"cmm {cmm.__version__} from {Path(cmm.__file__).parent}")
    print(f"cobra {cobra.__version__}")

    executable = shutil.which("cmm", path=str(Path(sys.executable).parent))
    check(
        executable is not None,
        "the cmm console script is installed beside the interpreter",
    )
    version = subprocess.run(
        [executable, "--version"], capture_output=True, text=True, check=True
    ).stdout.strip()
    check(cmm.__version__ in version, f"`cmm --version` reports {version!r}")

    model = load_model("textbook")
    model.solver = "glpk"
    status = solver_status(model)
    print(f"solver: {status.summary()}")
    check({"LP", "MILP"} <= set(status.capabilities), "GLPK provides LP and MILP")

    solution = fba(model)
    check(solution.status == "optimal", "FBA is optimal")
    check(
        abs(solution.objective_value - TEXTBOOK_GROWTH) < TOL,
        f"FBA growth {solution.objective_value:.10f} matches {TEXTBOOK_GROWTH}",
    )
    check(
        solution.metadata.get("solver") == "glpk"
        and solution.metadata.get("model_sha256"),
        f"FBA provenance records solver and model fingerprint {solution.metadata['model_sha256']}",
    )

    parsimonious = pfba(model)
    check(
        abs(parsimonious.fluxes[BIOMASS] - TEXTBOOK_GROWTH) < TOL,
        "pFBA keeps the FBA growth optimum",
    )
    check(
        sum(abs(v) for v in parsimonious.fluxes.values())
        <= sum(abs(v) for v in solution.fluxes.values()) + TOL,
        "pFBA total flux does not exceed the FBA solution's",
    )

    ranges = fva(model, reactions=[BIOMASS, "EX_glc__D_e"], processes=1)
    check(
        abs(ranges[BIOMASS].minimum - TEXTBOOK_GROWTH) < TOL
        and abs(ranges[BIOMASS].maximum - TEXTBOOK_GROWTH) < TOL,
        "FVA pins biomass at the optimum",
    )

    # File I/O through a path with a space and non-ASCII characters, as users' folders have.
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp) / "CMM smoke ü"
        folder.mkdir()
        sbml = folder / "e_coli_core.xml"
        write_sbml_model(model, str(sbml))
        reloaded = read_sbml_model(str(sbml))
        reloaded.solver = "glpk"
        check(
            len(reloaded.reactions) == len(model.reactions),
            "SBML round trip keeps reactions",
        )
        check(
            abs(fba(reloaded).objective_value - TEXTBOOK_GROWTH) < TOL,
            "the reloaded SBML model reproduces the FBA optimum",
        )

        table = folder / "fva.csv"
        ranges.to_frame().to_csv(table, index=False)
        with table.open(encoding="utf-8", newline="") as handle:
            rows = {row["reaction_id"]: row for row in csv.DictReader(handle)}
        check(
            abs(float(rows[BIOMASS]["maximum"]) - TEXTBOOK_GROWTH) < TOL,
            "FVA CSV export round-trips",
        )

    try:
        import gurobipy  # noqa: F401
    except ImportError:
        print("gurobipy not installed; GLPK-only smoke")
        return
    gurobi_model = load_model("textbook")
    gurobi_model.solver = "gurobi"
    print(f"solver: {solver_status(gurobi_model).summary()}")
    gurobi_growth = fba(gurobi_model).objective_value
    check(
        abs(gurobi_growth - solution.objective_value) < TOL,
        f"Gurobi FBA growth {gurobi_growth:.10f} agrees with GLPK",
    )


if __name__ == "__main__":
    main()
