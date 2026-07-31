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
from qtpy.QtCore import Qt  # noqa: E402
from qtpy.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def app():
    try:
        instance = QApplication.instance() or QApplication([])
    except Exception as exc:  # pragma: no cover - no usable Qt platform
        pytest.skip(f"no Qt platform available: {exc}")
    return instance


def _clear_ko(window):
    """Empty the Comparison tab's chosen-knockout set."""

    window._clear_ko_selected()


def _select_ko(window, ids):
    """Set the Comparison tab's chosen-knockout set to exactly the given ids."""

    window._clear_ko_selected()
    window._add_ko_targets(list(ids))


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
    with pytest.warns(UserWarning, match="infeasible"):
        window.run_fba()  # must not raise
    assert (
        "infeasible" in window.objective_label.text()
        or "error" in window.objective_label.text()
    )


def test_crossing_bound_edit_does_not_crash(app):
    window = CmmMainWindow(build_demo_model())
    # Setting lower above the current upper must not raise (atomic assignment + clamp).
    window.set_reaction_bounds("SUP_A", lower=2000.0)
    assert (
        window.model.reactions.SUP_A.lower_bound
        <= window.model.reactions.SUP_A.upper_bound
    )
    # In-table crossing edit clamps and reflects the stored value in the cell.
    sup_row = next(
        r
        for r in range(window.reaction_table.rowCount())
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
    assert (
        "failed" in window.yield_label.text()
        or "unavailable" in window.yield_label.text()
    )


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


def test_expression_vector_rejects_duplicates_and_negative_values(tmp_path):
    duplicate = tmp_path / "duplicate.csv"
    duplicate.write_text("gene,expression\ng1,1\ng1,2\n")
    with pytest.raises(ValueError, match="duplicate"):
        _read_expression_vector(str(duplicate))

    negative = tmp_path / "negative.csv"
    negative.write_text("gene,expression\ng1,-1\n")
    with pytest.raises(ValueError, match="non-negative"):
        _read_expression_vector(str(negative))


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


def test_media_pfba_comparison_in_gui(app, ecoli_core):
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
    window.comparison_method_combo.setCurrentText("MOMA (L1)")
    window.template_combo.setCurrentText("pfba")
    window.ko_level_combo.setCurrentText("reaction")
    _select_ko(window, ["PFK"])  # reroutable aerobically
    window.run_comparison()
    assert "MOMA" in window.comparison_summary.text()
    assert window.comparison_table.rowCount() > 0

    # A lethal knockout is reported as infeasible, not a crash.
    window.medium_combo.setCurrentText("glucose_anaerobic")
    window.apply_selected_medium()
    _select_ko(window, ["PFK"])
    window.run_comparison()
    assert "infeasible" in window.comparison_summary.text()


def test_fvseof_in_gui(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window.tabs.setCurrentIndex(window._tab_index("Production"))
    window.product_combo.setCurrentText("EX_succ_e")
    window.anaerobic_combo.setCurrentText("anaerobic")
    window.run_fvseof_plot()
    assert (
        "FVSEOF" in window.yield_label.text() or "robust" in window.yield_label.text()
    )


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


def test_omics_multi_condition_in_gui(app, tmp_path):
    from cmm.omics.conditions import read_expression_table

    path = tmp_path / "conditions.csv"
    path.write_text(
        "gene,condA,condB\ng1,50,50\ng2,100,1\ng3,1,100\ng5,1,100\ngb,50,50\n"
    )

    window = CmmMainWindow(build_demo_model())
    # Multi-condition is now folded into the Omics tab.
    assert window._tab_index("Omics") is not None
    assert window._tab_index("Multi-condition") is None

    # Loading detects both conditions (checked by default) but does not compute.
    window._set_omics_source(read_expression_table(str(path)), "conditions.csv")
    assert window._checked_omics_conditions() == ["condA", "condB"]
    assert window.omics_table.rowCount() == 0

    # Compute predicts one flux column per checked condition (column name = condition).
    # LAD is an LP so it runs without a QP solver in CI.
    window.omics_method_combo.setCurrentText("lad")
    window.compute_omics()  # must not raise
    headers = [
        window.omics_table.horizontalHeaderItem(i).text()
        for i in range(window.omics_table.columnCount())
    ]
    assert headers == ["Reaction", "condA", "condB"]
    assert window.omics_table.rowCount() > 0

    # Changing the checked set and recomputing reuses the same data (no reload) and drops
    # the unchecked condition's column.
    for i in range(window.omics_cond_list.count()):
        if window.omics_cond_list.item(i).text() == "condB":
            window.omics_cond_list.item(i).setCheckState(Qt.Unchecked)
    window.compute_omics()
    headers = [
        window.omics_table.horizontalHeaderItem(i).text()
        for i in range(window.omics_table.columnCount())
    ]
    assert headers == ["Reaction", "condA"]


def test_export_table_csv(app, tmp_path):
    window = CmmMainWindow(build_demo_model())
    window.run_fba()
    window.tabs.setCurrentIndex(window._tab_index("Simulation"))

    out = tmp_path / "fluxes.csv"
    from qtpy import QtWidgets

    QtWidgets.QFileDialog.getSaveFileName = staticmethod(
        lambda *a, **k: (str(out), "CSV")
    )
    window.export_table_csv()
    assert out.exists()
    header = out.read_text().splitlines()[0]
    assert header.startswith("Reaction")


def test_comparison_gene_knockout(app):
    window = CmmMainWindow(build_demo_model())
    window._goto_tab("Comparison")
    window.comparison_method_combo.setCurrentText("MOMA (L2)")
    window.template_combo.setCurrentText("pfba")
    window.ko_level_combo.setCurrentText("gene")
    _select_ko(
        window, ["g2"]
    )  # disease-branch gene -> blocks R2 -> reroute through R3/R5
    window.run_comparison()
    assert "gene g2" in window.comparison_summary.text()
    assert "1 reactions blocked" in window.comparison_summary.text()
    assert window.comparison_table.rowCount() > 0


def test_comparison_multi_knockout(app):
    window = CmmMainWindow(build_demo_model())
    window._goto_tab("Comparison")
    window.ko_level_combo.setCurrentText("reaction")
    _select_ko(window, ["R1", "R2"])
    window.run_comparison()
    # Both reactions named in the summary; knocking out R1 is lethal (no substrate route).
    text = window.comparison_summary.text()
    assert "R1" in text and "R2" in text


def test_comparison_requires_a_selection(app):
    window = CmmMainWindow(build_demo_model())
    window._goto_tab("Comparison")
    _clear_ko(window)
    window.run_comparison()
    assert "Add one or more" in window.comparison_summary.text()


def test_comparison_batch_over_genes(app):
    window = CmmMainWindow(build_demo_model())
    window._goto_tab("Comparison")
    window.comparison_method_combo.setCurrentText("MOMA (L2)")
    window.template_combo.setCurrentText("pfba")
    window.ko_level_combo.setCurrentText("gene")
    _clear_ko(window)  # nothing selected -> batch over all genes
    window.run_batch_comparison()
    assert window.comparison_table.columnCount() == 4  # no target product selected
    assert window.comparison_table.rowCount() == len(window.model.genes)
    assert "Batch" in window.comparison_summary.text()
    # The batch table reports wild-type vs post-knockout biomass and essentiality per target.
    headers = [
        window.comparison_table.horizontalHeaderItem(c).text()
        for c in range(window.comparison_table.columnCount())
    ]
    assert headers == ["Target", "WT Biomass", "KO Biomass", "Essential"]


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


def test_flux_response_tab_runs(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Flux Response")
    window.fr_target_combo.setCurrentText("PGI")
    window.fr_response_combo.setCurrentText("EX_succ_e")
    window.fr_growth_spin.setValue(30.0)
    window.fr_steps_spin.setValue(10)
    window.run_flux_response()

    assert window.fr_table.rowCount() == 10
    assert "EX_succ_e" in window.fr_summary.text()
    assert "PGI" in window.fr_summary.text()
    # The scan table is the tab's export target and the curve is its figure.
    assert window._active_table() is window.fr_table
    assert window._active_figure() is not None


def test_flux_response_tab_reports_infeasible_points(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Flux Response")
    window.fr_target_combo.setCurrentText("PGI")
    window.fr_response_combo.setCurrentText("EX_succ_e")
    window.fr_growth_spin.setValue(30.0)
    window.fr_steps_spin.setValue(12)
    window.run_flux_response()

    statuses = [
        window.fr_table.item(r, 3).text() for r in range(window.fr_table.rowCount())
    ]
    # A growth floor makes the extremes of the range unsolvable; they are shown, not hidden.
    assert "infeasible" in statuses
    infeasible_row = statuses.index("infeasible")
    assert window.fr_table.item(infeasible_row, 1).text() == "—"


def test_flux_response_tab_surfaces_errors(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Flux Response")
    window.fr_target_combo.setCurrentText("NOT_A_REACTION")
    window.run_flux_response()
    assert "failed" in window.fr_summary.text().lower()


def test_flux_response_manual_range(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Flux Response")
    window.fr_target_combo.setCurrentText("PGI")
    # The range boxes are always editable; narrowing them narrows the scan.
    assert window.fr_min_spin.isEnabled()
    window.fr_min_spin.setValue(-5.0)
    window.fr_max_spin.setValue(5.0)
    window.fr_steps_spin.setValue(6)
    window.run_flux_response()
    assert window.fr_table.rowCount() == 6
    assert float(window.fr_table.item(0, 0).text()) == pytest.approx(-5.0, abs=1e-6)


def test_sampling_tab_runs(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    window.sample_n_spin.setValue(120)
    window.sample_method_combo.setCurrentText("achr")
    window.sample_thinning_spin.setValue(10)
    window.sample_seed_spin.setValue(4)
    window.run_sampling()

    assert window.sample_table.rowCount() == len(ecoli_core.reactions)
    assert "120" in window.sample_summary.text()
    assert "seed 4" in window.sample_summary.text()
    assert window._active_table() is window.sample_table
    assert window._active_figure() is not None


def test_sampling_tab_reference_mode(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    # The reference controls stay disabled until the mode needs them.
    assert not window.sample_reference_combo.isEnabled()
    window.sample_mode_combo.setCurrentText("around a reference")
    assert window.sample_reference_combo.isEnabled()
    assert window.sample_window_spin.isEnabled()

    window.sample_n_spin.setValue(120)
    window.sample_method_combo.setCurrentText("achr")
    window.sample_thinning_spin.setValue(10)
    window.run_sampling()
    assert window.sample_table.rowCount() == len(ecoli_core.reactions)


def test_sampling_table_exports_to_csv(app, ecoli_core, tmp_path):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    window.sample_n_spin.setValue(120)
    window.sample_method_combo.setCurrentText("achr")
    window.sample_thinning_spin.setValue(10)
    window.run_sampling()

    out = tmp_path / "samples.csv"
    from qtpy import QtWidgets

    QtWidgets.QFileDialog.getSaveFileName = staticmethod(
        lambda *a, **k: (str(out), "CSV")
    )
    window.export_table_csv()
    assert out.exists()
    assert out.read_text().splitlines()[0].startswith("Reaction,Mean,Std")


def test_sampling_exports_raw_ensemble(app, ecoli_core, tmp_path):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    # Nothing to export until a run has happened.
    assert not window.sample_export_btn.isEnabled()

    window.sample_n_spin.setValue(120)
    window.sample_method_combo.setCurrentText("achr")
    window.sample_thinning_spin.setValue(10)
    window.run_sampling()
    assert window.sample_export_btn.isEnabled()

    out = tmp_path / "ensemble.csv"
    from qtpy import QtWidgets

    QtWidgets.QFileDialog.getSaveFileName = staticmethod(
        lambda *a, **k: (str(out), "CSV")
    )
    window.export_sampling_samples()
    assert out.exists()

    import pandas as pd

    frame = pd.read_csv(out, index_col="sample")
    # One row per drawn sample, one column per reaction — the raw ensemble, not the summary.
    assert frame.shape == (120, len(ecoli_core.reactions))


def test_sampling_table_covers_every_reaction(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    window.sample_n_spin.setValue(120)
    window.sample_method_combo.setCurrentText("achr")
    window.sample_thinning_spin.setValue(10)
    window.sample_top_spin.setValue(5)
    window.run_sampling()

    # The statistics table is complete; only the figure is restricted to the top-N.
    assert window.sample_table.rowCount() == len(ecoli_core.reactions)
    listed = {
        window.sample_table.item(r, 0).text()
        for r in range(window.sample_table.rowCount())
    }
    assert listed == {r.id for r in ecoli_core.reactions}
    assert len(window._sampling_canvas.figure.axes[0].get_yticklabels()) == 5


def test_reloading_a_model_clears_the_ensemble(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    window.sample_n_spin.setValue(120)
    window.sample_method_combo.setCurrentText("achr")
    window.sample_thinning_spin.setValue(10)
    window.run_sampling()
    assert window._sampling_result is not None

    window.load_model(build_demo_model())
    assert window._sampling_result is None
    assert not window.sample_export_btn.isEnabled()
    assert window.sample_table.rowCount() == 0


def test_range_is_seeded_with_the_theoretical_interval(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Flux Response")
    window.fr_target_combo.setCurrentText("PGI")

    # Not a placeholder: PGI can carry -50..10 in this medium.
    assert window.fr_min_spin.value() == pytest.approx(-50.0, abs=1e-3)
    assert window.fr_max_spin.value() == pytest.approx(10.0, abs=1e-3)

    # Switching target re-seeds rather than keeping the previous reaction's numbers.
    window.fr_target_combo.setCurrentText("EX_o2_e")
    assert window.fr_min_spin.value() == pytest.approx(-60.0, abs=1e-3)
    assert window.fr_max_spin.value() == pytest.approx(0.0, abs=1e-3)

    # A manual edit survives until the user asks for a reset.
    window.fr_min_spin.setValue(-5.0)
    window._prefill_fr_range()
    assert window.fr_min_spin.value() == pytest.approx(-60.0, abs=1e-3)


def test_reference_controls_are_disabled_in_uniform_mode(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    reference_widgets = (
        window.sample_reference_label,
        window.sample_reference_combo,
        window.sample_window_label,
        window.sample_window_spin,
    )

    assert window.sample_mode_combo.currentText() == "uniform"
    # Uniform sampling has no reference, so neither the combo nor its label invites input.
    assert not any(widget.isEnabled() for widget in reference_widgets)

    window.sample_mode_combo.setCurrentText("around a reference")
    assert all(widget.isEnabled() for widget in reference_widgets)

    window.sample_mode_combo.setCurrentText("uniform")
    assert not any(widget.isEnabled() for widget in reference_widgets)


def test_tab_order_groups_related_analyses(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    names = [window.tabs.tabText(i) for i in range(window.tabs.count())]
    assert names == [
        "Simulation",
        "Sampling",
        "Comparison",
        "Flux Response",
        "Production",
        "Strain Design",
        "Omics",
        "Revert Metabolism",
        "Transform (A→B)",
    ]
    # Each new tab sits beside the analysis asking a similar question.
    assert names.index("Sampling") == names.index("Simulation") + 1
    assert names.index("Flux Response") == names.index("Comparison") + 1


def _run_sampling(window, **overrides):
    window.sample_n_spin.setValue(overrides.get("n", 200))
    window.sample_method_combo.setCurrentText("achr")
    window.sample_thinning_spin.setValue(10)
    window.run_sampling()


def test_sampling_applies_reaction_knockouts(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")

    _run_sampling(window)
    wild_type = window._sampling_result.samples["PGI"].mean()
    assert wild_type > 1e-6
    assert "wild type" in window.sample_summary.text()

    window._ko_add("sample_ko", ["PGI"])
    _run_sampling(window)
    # The deleted reaction carries no flux anywhere in the ensemble.
    assert window._sampling_result.samples["PGI"].abs().max() < 1e-6
    assert "PGI" in window.sample_summary.text()
    assert window._sampling_result.metadata["parameters"]["condition"] == "knockout:PGI"
    # The knockout is a scoped condition, not an edit to the loaded model.
    assert ecoli_core.reactions.PGI.bounds == (-1000.0, 1000.0)


def test_sampling_resolves_gene_knockouts_through_the_gpr(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    window.sample_ko_level_combo.setCurrentText("gene")
    # b4025 is pgi; the GPR must resolve it to the PGI reaction.
    window._ko_add("sample_ko", ["b4025"])
    _run_sampling(window)
    assert window._sampling_result.samples["PGI"].abs().max() < 1e-6
    assert "gene knockout" in window.sample_summary.text()


def test_sampling_level_switch_clears_the_selection(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    window._ko_add("sample_ko", ["PGI"])
    assert window._ko_targets("sample_ko") == ["PGI"]
    # Reaction ids are not gene ids, so switching level must not carry them over.
    window.sample_ko_level_combo.setCurrentText("gene")
    assert window._ko_targets("sample_ko") == []
    assert window.sample_ko_selected_label.text() == "Selected (0):"


def test_reference_mode_uses_a_knocked_out_reference(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    window.sample_mode_combo.setCurrentText("around a reference")
    window._ko_add("sample_ko", ["PGI"])
    # A wild-type reference would put PGI outside its own window and raise; the reference
    # has to be computed under the same knockouts.
    _run_sampling(window)
    assert window._sampling_result is not None
    assert window._sampling_result.samples["PGI"].abs().max() < 1e-6


def test_lethal_knockout_is_reported_not_crashed(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    _run_sampling(window)
    assert window._sampling_result is not None

    window.sample_ko_level_combo.setCurrentText("gene")
    window._ko_add("sample_ko", ["b2415"])  # essential in e_coli_core
    _run_sampling(window)

    assert "lethal" in window.sample_summary.text()
    # The previous run's ensemble must not survive as an exportable result.
    assert window._sampling_result is None
    assert not window.sample_export_btn.isEnabled()
    assert window.sample_table.rowCount() == 0


def test_knockout_picker_is_independent_per_tab(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._add_ko_targets(["PGI"])  # Comparison tab picker
    window._ko_add("sample_ko", ["ENO", "FBA"])  # Sampling tab picker
    assert window._selected_ko_targets() == ["PGI"]
    assert window._ko_targets("sample_ko") == ["ENO", "FBA"]
    window._clear_ko_selected()
    assert window._selected_ko_targets() == []
    assert window._ko_targets("sample_ko") == ["ENO", "FBA"]
