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

from cmm.app.main_window import CmmMainWindow
from cmm.resources import BUNDLED_MAPS  # noqa: E402
from cmm.core import Condition, ReactionBound, fba  # noqa: E402
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
# The curated Escher map CMM bundles for e_coli_core. Held as a Path (rather than the string
# bundled_map_for returns) because this module tests it with .exists() before rendering.
ESCHER_MAP = Path(BUNDLED_MAPS[0])

# Growth-coupled succinate design: go anaerobic and block the competing fermentation
# secretions (upper bound 0), leaving succinate (via FRD7) as the only NADH sink. Under
# growth maximization this forces succinate excretion from 0 up to ~10.
#
# Every bound this design depends on is declared here, including the CO2 supply. That is
# deliberate: reductive-TCA succinate runs through PEP carboxylase (pep + co2 -> oaa), so the
# route consumes CO2 and cannot carry flux without a source of it. Industrial anaerobic
# succinate fermentation supplies exactly that, by CO2 sparging or bicarbonate. Since 0.4.0
# the media layer closes CO2 uptake by default so that a product's carbon is attributed to the
# named substrate unless stated otherwise, and this design states otherwise.
#
# Measured on e_coli_core with glucose_aerobic applied, then these bounds:
#   CO2 supply declared (below)  -> growth 0.020154, succinate 9.8848, CO2 uptake 0.0147
#   CO2 uptake left closed       -> growth 0.000000, succinate 8.3900, CO2 uptake 0.0
# Fluxes are mmol gDW-1 h-1 and growth is h-1. This is a demonstration of bound editing; the
# yields reported by ``theoretical_yield`` elsewhere in this module use the default closed-CO2
# condition and are not affected by the supply declared here.
ENGINEERING_BOUNDS = {
    "EX_o2_e": (0.0, 1000.0),  # anaerobic
    "EX_co2_e": (
        -1000.0,
        1000.0,
    ),  # CO2 supplied (sparging/bicarbonate); PPC requires it
    "EX_ac_e": (0.0, 0.0),  # block acetate secretion
    "EX_etoh_e": (0.0, 0.0),  # block ethanol
    "EX_for_e": (0.0, 0.0),  # block formate
    "EX_lac__D_e": (0.0, 0.0),  # block lactate
}


#: 0.4.0 removed ``aerobic=`` from cmm.features.production: the condition is stated once, the
#: same way every other caller states it. CO2 uptake is closed alongside oxygen, matching what
#: the ``glucose_anaerobic`` preset does through the media layer.
ANAEROBIC = Condition(
    name="anaerobic (O2 and CO2 uptake closed)",
    bounds=(
        ReactionBound(reaction_id="EX_o2_e", lower_bound=0.0),
        ReactionBound(reaction_id="EX_co2_e", lower_bound=0.0),
    ),
)


def _output_dir() -> Path:
    configured = os.environ.get("CMM_OUTPUT_DIR")
    if configured:
        return Path(configured).expanduser()
    return Path(__file__).resolve().parents[3] / "temp_figures"


OUTPUT_DIR = _output_dir()


def generate_publication_figures(model) -> list[Path]:
    """Standalone paper-ready figures for the succinate scenario."""

    saved: list[Path] = []

    aerobic_yield = theoretical_yield(model, SUCC)
    anaerobic_yield = theoretical_yield(model, SUCC, condition=ANAEROBIC)
    saved.append(
        save_figure(
            yield_figure(
                [aerobic_yield, anaerobic_yield], title="Succinate theoretical yield"
            ),
            OUTPUT_DIR / "fig_succinate_yield.png",
        )
    )

    envelope = production_envelope(model, SUCC, points=20)
    saved.append(
        save_figure(
            production_envelope_figure(envelope, title="Succinate production envelope"),
            OUTPUT_DIR / "fig_succinate_envelope.png",
        )
    )

    result = fseof(model, SUCC, n_steps=10, condition=ANAEROBIC)
    saved.append(
        save_figure(
            fseof_figure(
                result, top_n=6, title="Succinate FSEOF amplification targets"
            ),
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
                dict(wt),
                dict(engineered),
                reactions,
                reference_label="wild type",
                comparison_label="succinate design",
                title="Growth-coupled succinate design (anaerobic + byproduct knockouts)",
            ),
            OUTPUT_DIR / "fig_succinate_flux_shift.png",
        )
    )

    # Dependency-free schematic carbon-backbone network (fallback when no map exists).
    saved.append(
        save_figure(
            network_flux_map(
                model,
                dict(engineered),
                top_n=12,
                title="Succinate design carbon-flux network (schematic)",
            ),
            OUTPUT_DIR / "fig_succinate_network.png",
        )
    )

    # Escher-layout flux maps: wild type and the succinate design.
    if ESCHER_MAP.exists():
        saved.append(
            save_figure(
                escher_flux_map(
                    ESCHER_MAP,
                    dict(wt),
                    title="E. coli core — wild-type flux (Escher layout)",
                ),
                OUTPUT_DIR / "fig_escher_wt.png",
            )
        )
        saved.append(
            save_figure(
                escher_flux_map(
                    ESCHER_MAP,
                    dict(engineered),
                    title="E. coli core — succinate design flux (Escher layout)",
                ),
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
            (
                i
                for i in range(window.tabs.count())
                if window.tabs.tabText(i) == "Flux Map"
            ),
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
