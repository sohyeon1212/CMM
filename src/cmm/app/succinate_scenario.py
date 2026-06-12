"""E. coli succinate production-design scenario: figures + GUI screenshots.

Run with the venv and the offscreen platform::

    QT_QPA_PLATFORM=offscreen CMM_OUTPUT_DIR=./temp_figures_new \
        PYTHONPATH=src .venv/bin/python -m cmm.app.succinate_scenario

Loads e_coli_core, exercises every production-design feature for succinate (theoretical
yield, production envelope, FSEOF, and an anaerobic bound edit that increases succinate), and
writes both standalone publication figures and GUI screenshots to ``CMM_OUTPUT_DIR`` or
CMM/temp_figures/.
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from cobra.io import load_model  # noqa: E402
from qtpy.QtWidgets import QApplication  # noqa: E402

from cmm.app.main_window import CmmMainWindow  # noqa: E402
from cmm.core import fba  # noqa: E402
from cmm.features.production import (  # noqa: E402
    fseof,
    production_envelope,
    theoretical_yield,
)
from cmm.visualization import (  # noqa: E402
    escher_flux_map,
    flux_comparison_figure,
    fseof_figure,
    network_flux_map,
    production_envelope_figure,
    save_figure,
    yield_figure,
)

SUCC = "EX_succ_e"
BIOMASS_ID = "Biomass_Ecoli_core"
# Project root is the CMM package root (src/cmm/app -> CMM): parents[3]. Keep all data
# under CMM/test_data so the scenario is self-contained and never reaches a parent folder.
_REPO_ROOT = Path(__file__).resolve().parents[3]
# Curated Escher map for e_coli_core (user data in test_data; CMM bundles no maps).
ESCHER_MAP = _REPO_ROOT / "test_data" / "e_coli_core.Core metabolism.json"

# Growth-coupled succinate design: go anaerobic and block the competing fermentation
# secretions (upper bound 0), leaving succinate (via FRD7) as the only NADH sink. Under
# growth maximization this forces succinate excretion from 0 up to ~10.
ENGINEERING_BOUNDS = {
    "EX_o2_e": (0.0, 1000.0),    # anaerobic
    "EX_ac_e": (0.0, 0.0),       # block acetate secretion
    "EX_etoh_e": (0.0, 0.0),     # block ethanol
    "EX_for_e": (0.0, 0.0),      # block formate
    "EX_lac__D_e": (0.0, 0.0),   # block lactate
}


def _output_dir() -> Path:
    configured = os.environ.get("CMM_OUTPUT_DIR")
    if configured:
        return Path(configured).expanduser()
    return Path(__file__).resolve().parents[3] / "temp_figures"


OUTPUT_DIR = _output_dir()


def generate_publication_figures(model) -> list[Path]:
    """Standalone paper-ready figures for the succinate scenario."""

    saved: list[Path] = []

    aerobic_yield = theoretical_yield(model, SUCC, aerobic=True)
    anaerobic_yield = theoretical_yield(model, SUCC, aerobic=False)
    saved.append(
        save_figure(
            yield_figure([aerobic_yield, anaerobic_yield], title="Succinate theoretical yield"),
            OUTPUT_DIR / "fig_succinate_yield.png",
        )
    )

    envelope = production_envelope(model, SUCC, aerobic=True, points=20)
    saved.append(
        save_figure(
            production_envelope_figure(envelope, title="Succinate production envelope"),
            OUTPUT_DIR / "fig_succinate_envelope.png",
        )
    )

    result = fseof(model, SUCC, n_steps=10, aerobic=False)
    saved.append(
        save_figure(
            fseof_figure(result, top_n=6, title="Succinate FSEOF amplification targets"),
            OUTPUT_DIR / "fig_succinate_fseof.png",
        )
    )

    # WT (aerobic, growth mode) vs the growth-coupled succinate design.
    wt = fba(model).fluxes
    with model:
        for rid, bounds in ENGINEERING_BOUNDS.items():
            model.reactions.get_by_id(rid).bounds = bounds
        engineered = fba(model).fluxes
    # Compare design *outputs* only; EX_o2_e is an input constraint (its large WT bar would
    # dwarf the product changes), so it is intentionally excluded.
    reactions = ["EX_succ_e", "FRD7", "PPC", "EX_ac_e", "EX_etoh_e"]
    reactions = [r for r in reactions if r in model.reactions]
    saved.append(
        save_figure(
            flux_comparison_figure(
                dict(wt), dict(engineered), reactions,
                reference_label="wild type", comparison_label="succinate design",
                title="Growth-coupled succinate design (anaerobic + byproduct knockouts)",
            ),
            OUTPUT_DIR / "fig_succinate_flux_shift.png",
        )
    )

    # Dependency-free schematic carbon-backbone network (fallback when no map exists).
    saved.append(
        save_figure(
            network_flux_map(model, dict(engineered), top_n=12,
                             title="Succinate design carbon-flux network (schematic)"),
            OUTPUT_DIR / "fig_succinate_network.png",
        )
    )

    # Escher-layout flux maps: wild type and the succinate design.
    if ESCHER_MAP.exists():
        saved.append(
            save_figure(
                escher_flux_map(ESCHER_MAP, dict(wt),
                                title="E. coli core — wild-type flux (Escher layout)"),
                OUTPUT_DIR / "fig_escher_wt.png",
            )
        )
        saved.append(
            save_figure(
                escher_flux_map(ESCHER_MAP, dict(engineered),
                                title="E. coli core — succinate design flux (Escher layout)"),
                OUTPUT_DIR / "fig_escher_succinate.png",
            )
        )
    return saved


def _capture(window: CmmMainWindow, app: QApplication, name: str) -> Path:
    window.show()
    app.processEvents()
    app.processEvents()
    path = OUTPUT_DIR / f"{name}.png"
    window.grab().save(str(path))
    return path


def run_gui_scenarios() -> list[Path]:
    app = QApplication.instance() or QApplication([])
    model = load_model("textbook")
    map_path = str(ESCHER_MAP) if ESCHER_MAP.exists() else None
    window = CmmMainWindow(model, default_product=SUCC, map_path=map_path)
    saved: list[Path] = []

    # Production tab: theoretical yield, envelope, FSEOF for succinate.
    window._goto_tab("Production")
    window.anaerobic_combo.setCurrentText("aerobic")
    window.run_theoretical_yield()
    saved.append(_capture(window, app, "ecoli_01_theoretical_yield"))

    window.run_production_envelope_plot()
    saved.append(_capture(window, app, "ecoli_02_production_envelope"))

    window.anaerobic_combo.setCurrentText("anaerobic")
    window.run_fseof_plot()
    saved.append(_capture(window, app, "ecoli_03_fseof_targets"))

    # Simulation tab + bound edits: WT growth mode, then the growth-coupled design
    # (anaerobic + byproduct knockouts) raises succinate from 0 to ~10.
    window._goto_tab("Simulation")
    window.run_fba()
    wt_succ = window._fluxes.get(SUCC, 0.0)
    for rid, (lower, upper) in ENGINEERING_BOUNDS.items():
        window.set_reaction_bounds(rid, lower=lower, upper=upper)
    window.run_fba()
    design_succ = window._fluxes.get(SUCC, 0.0)
    window.flux_slider.setValue(20)
    saved.append(_capture(window, app, "ecoli_04_succinate_bound_edits"))

    # Flux Map tab: Escher-layout map coloured by the current (engineered) flux.
    if map_path:
        fluxmap_index = next(
            (i for i in range(window.tabs.count()) if window.tabs.tabText(i) == "Flux Map"),
            None,
        )
        if fluxmap_index is not None:
            window.tabs.setCurrentIndex(fluxmap_index)
            window.render_flux_map()
            saved.append(_capture(window, app, "ecoli_05_flux_map"))

    growth = window._fluxes.get(BIOMASS_ID, 0.0)
    print(
        f"succinate flux: wild type={wt_succ:.3f} -> succinate design={design_succ:.3f} "
        f"(growth {growth:.3f} h^-1 — near-zero, growth-coupled, not a free win)"
    )

    return saved


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    model = load_model("textbook")
    figures = generate_publication_figures(model)
    screenshots = run_gui_scenarios()
    for path in figures + screenshots:
        size = path.stat().st_size if path.exists() else 0
        print(f"saved {path.name}  ({size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
