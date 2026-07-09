"""Headless smoke test for the Qt shell. Runs the real widgets under the offscreen platform."""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402

pytest.importorskip("qtpy")

from cmm.app.main_window import CmmMainWindow, _read_expression_vector  # noqa: E402
from cmm.app.screenshots import (  # noqa: E402
    SOURCE_EXPRESSION,
    TARGET_EXPRESSION,
    build_demo_model,
)
from qtpy.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def app():
    try:
        instance = QApplication.instance() or QApplication([])
    except Exception as exc:  # pragma: no cover - no usable Qt platform
        pytest.skip(f"no Qt platform available: {exc}")
    return instance


def test_window_builds_and_drives_services(app):
    window = CmmMainWindow(build_demo_model())

    window.run_fba()
    assert "optimal" in window.objective_label.text()
    assert window.reaction_table.rowCount() == 6

    window.run_fva()
    assert window.sim_table.rowCount() == 6

    window.run_revert(SOURCE_EXPRESSION, TARGET_EXPRESSION)
    # The disease-branch gene must be the top-ranked normalization target in the GUI table.
    assert window.revert_table.item(0, 1).text() == "g2"
    assert "g2" in window.revert_summary.text()


def test_bound_editing_updates_model_and_fba(app):
    window = CmmMainWindow(build_demo_model())
    window.run_fba()
    assert "10" in window.objective_label.text()

    # Programmatic bound edit reflects in the model and the table.
    window.set_reaction_bounds("SUP_A", upper=4.0)
    assert window.model.reactions.SUP_A.upper_bound == 4.0
    window.run_fba()
    assert "4" in window.objective_label.text()  # growth now substrate-limited to 4

    # In-table edit path: changing a cell triggers _on_bound_edited.
    sup_row = next(
        r
        for r in range(window.reaction_table.rowCount())
        if window.reaction_table.item(r, 0).text() == "SUP_A"
    )
    window.reaction_table.item(sup_row, 2).setText("7")
    assert window.model.reactions.SUP_A.upper_bound == 7.0
    window.run_fba()
    assert "7" in window.objective_label.text()

    # Invalid input is rejected and the cell reverts.
    window.reaction_table.item(sup_row, 2).setText("not a number")
    assert window.model.reactions.SUP_A.upper_bound == 7.0


def test_infeasible_fba_does_not_crash(app):
    window = CmmMainWindow(build_demo_model())
    # Force infeasibility: demand uptake the network cannot satisfy.
    window.model.reactions.SUP_A.bounds = (5.0, 5.0)
    window.model.reactions.BIOMASS.lower_bound = 1000.0
    window.run_fba()  # must not raise
    assert "infeasible" in window.objective_label.text() or "error" in window.objective_label.text()


def test_crossing_bound_edit_does_not_crash(app):
    window = CmmMainWindow(build_demo_model())
    # Setting lower above the current upper must not raise (atomic assignment + clamp).
    window.set_reaction_bounds("SUP_A", lower=2000.0)
    assert window.model.reactions.SUP_A.lower_bound <= window.model.reactions.SUP_A.upper_bound
    # In-table crossing edit clamps and reflects the stored value in the cell.
    sup_row = next(
        r for r in range(window.reaction_table.rowCount())
        if window.reaction_table.item(r, 0).text() == "SUP_A"
    )
    window.reaction_table.item(sup_row, 1).setText("99999")
    lo = window.model.reactions.SUP_A.lower_bound
    assert lo <= window.model.reactions.SUP_A.upper_bound
    assert window.reaction_table.item(sup_row, 1).text() == f"{lo:g}"


def test_production_buttons_disabled_without_exchanges(app):
    from cobra import Metabolite, Model, Reaction

    model = Model("no_exchange")
    a = Metabolite("a_c", compartment="c")
    b = Metabolite("b_c", compartment="c")
    r = Reaction("R")  # internal (2 metabolites) -> not a boundary/exchange reaction
    r.add_metabolites({a: -1, b: 1})
    r.gene_reaction_rule = "g"
    model.add_reactions([r])
    model.objective = r
    assert len(model.exchanges) == 0
    window = CmmMainWindow(model)
    assert all(not btn.isEnabled() for btn in window._production_buttons)
    # Even if invoked directly, the guarded runner must not raise.
    window.run_theoretical_yield()
    assert "failed" in window.yield_label.text() or "unavailable" in window.yield_label.text()


def test_menu_bar_has_expected_menus(app):
    window = CmmMainWindow(build_demo_model())
    titles = [a.text().replace("&", "") for a in window.menuBar().actions()]
    for expected in ("Analysis", "Model", "Config"):
        assert expected in titles


def test_omics_integration_in_gui(app):
    window = CmmMainWindow(build_demo_model())
    window.omics_method_combo.setCurrentText("eflux2")
    window.run_omics_demo()
    assert window.omics_table.rowCount() > 0
    assert "EFLUX2" in window.omics_summary.text()

    window.omics_method_combo.setCurrentText("lad")
    window.run_omics_demo()
    assert "LAD" in window.omics_summary.text()


def test_expression_vector_reads_tsv(tmp_path):
    path = tmp_path / "expression.tsv"
    path.write_text("gene\texpression\ng1\t2.5\ng2\t7\n")

    assert _read_expression_vector(str(path)) == {"g1": 2.5, "g2": 7.0}


def test_revert_tab_runs_loaded_expressions(app):
    window = CmmMainWindow(build_demo_model())
    assert not window.revert_run_btn.isEnabled()

    window._revert_source_expression = dict(SOURCE_EXPRESSION)
    window._revert_target_expression = dict(TARGET_EXPRESSION)
    window._update_revert_run_state()
    assert window.revert_run_btn.isEnabled()

    window.run_loaded_revert()
    assert window.revert_table.item(0, 1).text() == "g2"
    assert "g2" in window.revert_summary.text()


def test_media_pfba_comparison_in_gui(app, ecoli_core, unrestricted_qp_solver):
    window = CmmMainWindow(ecoli_core)

    # Apply a preset medium and confirm it changed the model's exchange bounds.
    window.medium_combo.setCurrentText("glucose_anaerobic")
    window.apply_selected_medium()
    assert window.model.reactions.EX_o2_e.lower_bound == 0.0

    # pFBA runs and reports growth.
    window.tabs.setCurrentIndex(0)
    window.run_pfba()
    assert "pFBA" in window.objective_label.text()

    # Comparison: MOMA against a pFBA template after a (feasible) knockout.
    window.medium_combo.setCurrentText("glucose_aerobic")
    window.apply_selected_medium()
    window.comparison_method_combo.setCurrentText("MOMA (L2)")
    window.template_combo.setCurrentText("pfba")
    window.ko_combo.setCurrentText("PFK")  # reroutable aerobically
    window.run_comparison()
    assert "MOMA" in window.comparison_summary.text()
    assert window.comparison_table.rowCount() > 0

    # A lethal knockout is reported as infeasible, not a crash.
    window.medium_combo.setCurrentText("glucose_anaerobic")
    window.apply_selected_medium()
    window.ko_combo.setCurrentText("PFK")
    window.run_comparison()
    assert "infeasible" in window.comparison_summary.text()


def test_fvseof_in_gui(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window.tabs.setCurrentIndex(window._tab_index("Production"))
    window.product_combo.setCurrentText("EX_succ_e")
    window.anaerobic_combo.setCurrentText("anaerobic")
    window.run_fvseof_plot()
    assert "FVSEOF" in window.yield_label.text() or "robust" in window.yield_label.text()


def test_menu_bar_has_file_menu(app):
    window = CmmMainWindow(build_demo_model())
    titles = [a.text().replace("&", "") for a in window.menuBar().actions()]
    assert "File" in titles


def test_strain_design_tab_present_and_disabled_without_exchanges(app):
    from cobra import Metabolite, Model, Reaction

    # Strain Design tab exists on any model.
    window = CmmMainWindow(build_demo_model())
    assert window._tab_index("Strain Design") is not None

    # ...but Run design is disabled when the model has no exchange reactions.
    model = Model("no_exchange")
    a = Metabolite("a_c", compartment="c")
    b = Metabolite("b_c", compartment="c")
    r = Reaction("R")
    r.add_metabolites({a: -1, b: 1})
    model.add_reactions([r])
    model.objective = r
    window2 = CmmMainWindow(model)
    assert not window2.sd_run_btn.isEnabled()


def test_transformation_tab_runs(app):
    window = CmmMainWindow(build_demo_model())
    assert window._tab_index("Transform (A→B)") is not None
    assert not window.transform_run_btn.isEnabled()

    window._transform_source_expression = dict(SOURCE_EXPRESSION)
    window._transform_target_expression = dict(TARGET_EXPRESSION)
    window._update_transform_run_state()
    assert window.transform_run_btn.isEnabled()

    window.transform_method_combo.setCurrentText("moma")
    window.transform_omics_combo.setCurrentText("eflux2")
    window.run_transformation()
    # Every gene is a candidate knockout; the ranking table is populated and summarized.
    assert window.transform_table.rowCount() > 0
    assert "A→B target" in window.transform_summary.text()


def test_conditions_tab_runs(app, tmp_path):
    from cmm.omics.conditions import read_expression_table

    path = tmp_path / "conditions.csv"
    path.write_text("gene,condA,condB\ng1,50,50\ng2,100,1\ng3,1,100\ng5,1,100\ngb,50,50\n")

    window = CmmMainWindow(build_demo_model())
    assert window._tab_index("Multi-condition") is not None

    window._condition_table = read_expression_table(str(path))
    window.cond_source_combo.addItems(["condA", "condB"])
    window.cond_target_combo.addItems(["condA", "condB"])
    window.cond_source_combo.setCurrentText("condA")
    window.cond_target_combo.setCurrentText("condB")
    window.cond_run_btn.setEnabled(True)
    window.cond_method_combo.setCurrentText("eflux2")
    window.run_condition_comparison()  # must not raise
    assert "condA" in window.cond_summary.text() and "condB" in window.cond_summary.text()

    # Identical conditions are rejected with a clear message, not a crash.
    window.cond_target_combo.setCurrentText("condA")
    window.run_condition_comparison()
    assert "different conditions" in window.cond_summary.text()


def test_export_table_csv(app, tmp_path):
    window = CmmMainWindow(build_demo_model())
    window.run_fba()
    window.tabs.setCurrentIndex(window._tab_index("Simulation"))

    out = tmp_path / "fluxes.csv"
    from qtpy import QtWidgets

    QtWidgets.QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: (str(out), "CSV"))
    window.export_table_csv()
    assert out.exists()
    header = out.read_text().splitlines()[0]
    assert header.startswith("Reaction")


def test_background_execution_runs_off_thread_and_stays_responsive(app):
    import threading
    import time

    from qtpy.QtCore import QTimer

    window = CmmMainWindow(build_demo_model())
    main_tid = threading.get_ident()
    state = {"pumped": False, "worker_tid": None, "running": False}

    # A main-thread timer that should fire WHILE the worker computes — proving the UI event
    # loop keeps running (the window would repaint) instead of freezing on the solve.
    def on_timer():
        if state["running"]:
            state["pumped"] = True

    QTimer.singleShot(40, on_timer)

    def compute():
        state["worker_tid"] = threading.get_ident()
        state["running"] = True
        time.sleep(0.3)
        state["running"] = False
        return "done"

    result = window._run_in_background(compute, label="test")
    assert result == "done"  # synchronous: result ready on return
    assert state["worker_tid"] != main_tid  # ran off the UI thread
    assert state["pumped"]  # the UI event loop pumped during the compute


def test_background_execution_reraises_errors(app):
    window = CmmMainWindow(build_demo_model())

    def boom():
        raise ValueError("kaboom")

    with pytest.raises(ValueError, match="kaboom"):
        window._run_in_background(boom, label="test")

    # The busy flag is cleared even when the compute raised, so the window stays usable.
    assert window._busy is False


def test_window_renders_non_blank(app):
    window = CmmMainWindow(build_demo_model())
    window.run_fba()
    window.show()
    app.processEvents()

    pixmap = window.grab()
    assert pixmap.width() > 400
    assert pixmap.height() > 300

    image = pixmap.toImage()
    assert not image.isNull()
    # Sample pixels across the canvas; a blank render would be a single uniform color.
    colors = {
        image.pixel(x, y)
        for x in range(0, image.width(), max(1, image.width() // 20))
        for y in range(0, image.height(), max(1, image.height() // 20))
    }
    assert len(colors) > 5
