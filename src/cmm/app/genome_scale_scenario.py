"""Genome-scale GUI scenario on an arbitrary SBML model.

Run with the venv and the offscreen platform; pass any SBML model path (defaults to the
public e_coli_core textbook model so the demo is self-contained)::

    QT_QPA_PLATFORM=offscreen CMM_OUTPUT_DIR=./temp_figures_new \
        PYTHONPATH=src .venv/bin/python -m cmm.app.genome_scale_scenario [model.xml]

Drives the GUI through FBA, a glucose bound edit, succinate production design, and omics
integration, capturing offscreen screenshots to ``CMM_OUTPUT_DIR`` or CMM/temp_figures/.
This validates the platform at the scale of a user-supplied genome-scale reconstruction.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import cobra  # noqa: E402
from cobra.io import load_model  # noqa: E402
from qtpy.QtWidgets import QApplication  # noqa: E402

from cmm.app.main_window import CmmMainWindow  # noqa: E402

SUCC = "EX_succ_e"
GLC = "EX_glc__D_e"


def _output_dir() -> Path:
    configured = os.environ.get("CMM_OUTPUT_DIR")
    if configured:
        return Path(configured).expanduser()
    return Path(__file__).resolve().parents[3] / "temp_figures"


OUTPUT_DIR = _output_dir()


def _load(model_path: str | None):
    if model_path:
        return cobra.io.read_sbml_model(model_path)
    return load_model("textbook")


def _capture(window: CmmMainWindow, app: QApplication, name: str) -> Path:
    window.show()
    app.processEvents()
    app.processEvents()
    path = OUTPUT_DIR / f"{name}.png"
    window.grab().save(str(path))
    return path


def run_scenarios(model_path: str | None = None, prefix: str = "genome") -> list[Path]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    app = QApplication.instance() or QApplication([])
    model = _load(model_path)
    default_product = SUCC if SUCC in model.reactions else None
    window = CmmMainWindow(model, default_product=default_product)
    saved: list[Path] = []

    saved.append(_capture(window, app, f"{prefix}_01_model_loaded"))

    window._goto_tab("Simulation")
    window.run_fba()
    saved.append(_capture(window, app, f"{prefix}_02_fba"))

    # Open glucose uptake via bound editing if the substrate exists but is closed.
    if GLC in model.reactions:
        window.set_reaction_bounds(GLC, lower=-10.0)
        if "EX_co2_e" in model.reactions:
            window.set_reaction_bounds("EX_co2_e", lower=0.0)
        window.run_fba()
        window.flux_slider.setValue(15)
        saved.append(_capture(window, app, f"{prefix}_03_glucose_medium"))

    if default_product is not None:
        window._goto_tab("Production")
        window.run_theoretical_yield()
        saved.append(_capture(window, app, f"{prefix}_04_theoretical_yield"))
        window.run_production_envelope_plot()
        saved.append(_capture(window, app, f"{prefix}_05_production_envelope"))

    if model.genes:
        window._goto_tab("Omics")
        window.omics_method_combo.setCurrentText("eflux2")
        window.run_omics_demo()
        saved.append(_capture(window, app, f"{prefix}_06_omics_eflux2"))

    print(
        f"model {model.id}: {len(model.reactions)} reactions, {len(model.genes)} genes"
    )
    return saved


def main() -> int:
    model_path = sys.argv[1] if len(sys.argv) > 1 else None
    for path in run_scenarios(model_path):
        size = path.stat().st_size if path.exists() else 0
        print(f"saved {path.name}  ({size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
