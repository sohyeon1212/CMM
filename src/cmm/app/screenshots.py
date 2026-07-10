"""Drive the CMM GUI through scenarios and capture offscreen screenshots.

Run with the project venv and the offscreen Qt platform::

    QT_QPA_PLATFORM=offscreen CMM_OUTPUT_DIR=./temp_figures_new \
        PYTHONPATH=src .venv/bin/python -m cmm.app.screenshots

Each scenario performs real GUI actions (FBA, FVA, the flux-range slider, revert-metabolism
prediction) and saves a rendered PNG to ``CMM_OUTPUT_DIR`` or CMM/temp_figures/. This is
the project's GUI test: the captures are inspected to confirm the platform renders and the
services drive the views.
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from cobra import Metabolite, Model, Reaction  # noqa: E402
from qtpy.QtWidgets import QApplication  # noqa: E402

from cmm.app.main_window import CmmMainWindow  # noqa: E402

# Disease (source) vs healthy (target) demo expression for the branched model.
SOURCE_EXPRESSION = {"g1": 50.0, "g2": 100.0, "g3": 1.0, "g5": 1.0, "gb": 50.0}
TARGET_EXPRESSION = {"g1": 50.0, "g2": 1.0, "g3": 100.0, "g5": 100.0, "gb": 50.0}


def _output_dir() -> Path:
    configured = os.environ.get("CMM_OUTPUT_DIR")
    if configured:
        return Path(configured).expanduser()
    return Path(__file__).resolve().parents[3] / "temp_figures"


OUTPUT_DIR = _output_dir()


def build_demo_model() -> Model:
    """Two-branch disease/healthy demo model (mirrors the test fixture)."""

    model = Model("demo_disease")
    a = Metabolite("A_c", compartment="c")
    b = Metabolite("B_c", compartment="c")
    d = Metabolite("D_c", compartment="c")
    p = Metabolite("P_c", compartment="c")

    def make(rid, stoich, gene, lb=0.0, ub=1000.0):
        r = Reaction(rid)
        r.lower_bound = lb
        r.upper_bound = ub
        r.add_metabolites(stoich)
        if gene:
            r.gene_reaction_rule = gene
        return r

    model.add_reactions(
        [
            make("SUP_A", {a: 1}, "", ub=10.0),
            make("R1", {a: -1, b: 1}, "g1"),
            make("R2", {b: -1, p: 1}, "g2"),
            make("R3", {b: -1, d: 1}, "g3"),
            make("R5", {d: -1, p: 1}, "g5"),
            make("BIOMASS", {p: -1}, "gb"),
        ]
    )
    model.objective = model.reactions.BIOMASS
    model.objective_direction = "max"
    return model


def _capture(window: CmmMainWindow, app: QApplication, name: str) -> Path:
    window.show()
    app.processEvents()
    app.processEvents()
    path = OUTPUT_DIR / f"{name}.png"
    pixmap = window.grab()
    pixmap.save(str(path))
    return path


def run_scenarios() -> list[Path]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    app = QApplication.instance() or QApplication([])

    model = build_demo_model()
    window = CmmMainWindow(model)
    saved: list[Path] = []

    # 1. Model loaded.
    saved.append(_capture(window, app, "01_model_loaded"))

    # 2. FBA on the Simulation tab.
    window._goto_tab("Simulation")
    window.run_fba()
    saved.append(_capture(window, app, "02_fba"))

    # 3. FVA.
    window.run_fva()
    saved.append(_capture(window, app, "03_fva"))

    # 4. Flux-range slider highlighting high-flux reactions in the model panel.
    window.flux_slider.setValue(50)
    saved.append(_capture(window, app, "04_flux_range_slider"))

    # 5. Revert-metabolism prediction with rMTA (gene knockouts).
    window._goto_tab("Revert Metabolism")
    window.method_combo.setCurrentText("rmta")
    window.run_revert(SOURCE_EXPRESSION, TARGET_EXPRESSION)
    saved.append(_capture(window, app, "05_revert_rmta"))

    # 6. Revert-metabolism with the published single-MTA MIQP mode.
    window.method_combo.setCurrentText("mta")
    window.run_revert(SOURCE_EXPRESSION, TARGET_EXPRESSION)
    saved.append(_capture(window, app, "06_revert_mta_miqp"))

    return saved


def main() -> int:
    saved = run_scenarios()
    for path in saved:
        size = path.stat().st_size if path.exists() else 0
        print(f"saved {path}  ({size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
