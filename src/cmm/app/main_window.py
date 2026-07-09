"""CMM desktop main window.

A thin Qt view over the solver-neutral services. All numerics live in cmm.core / cmm.features
/ cmm.omics; this module only renders results and wires controls, so the GUI inherits none of
the analysis logic.
"""

from __future__ import annotations

import html
from collections.abc import Mapping

from cobra import Model
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from qtpy.QtCore import QEventLoop, QObject, Qt, QThread, QTimer, Signal
from qtpy.QtGui import QColor, QFont
from qtpy.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QListWidget,
    QMainWindow,
    QMessageBox,
    QProgressDialog,
    QPushButton,
    QSlider,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from cmm.core import (
    PRESET_MEDIA,
    active_solver,
    apply_medium,
    fba,
    fva,
    pfba,
    reference_state_pfba,
    solver_status,
)
from cmm.features._perturbation import (
    blocked_reactions_for_genes,
    gene_perturbations,
    reaction_perturbations,
)
from cmm.features.comparison import (
    batch_comparison,
    knockout_comparison,
    reference_flux,
)
from cmm.features.production import (
    fseof,
    fvseof,
    production_envelope,
    theoretical_yield,
)
from cmm.features.revert import revert_targets
from cmm.features.strain_design import optknock, robustknock
from cmm.features.transformation import transformation_targets
from cmm.omics.conditions import (
    flux_log_change,
    predict_condition_fluxes,
    read_expression_table,
)
from cmm.omics.differential import differential_expression
from cmm.omics.expression import gene_to_reaction_weights, integrate_expression
from cmm.visualization import (
    escher_flux_map,
    fseof_figure,
    fvseof_figure,
    production_envelope_figure,
)

# "Oxford Slate" — an academic / scientific-software palette: a deep-navy header with a
# restrained steel-blue accent over a cool-paper canvas. Every interactive widget (menus,
# combos, spin boxes, the flux slider, scrollbars, the status strip, tooltips) is styled so
# the platform reads as one cohesive professional tool rather than default Qt chrome.
#
#   navy  #1b2b44   accent #2f5e8f (hover #274f79 / pressed #1f3f61)
#   canvas #eef1f5  panels #ffffff  text #1a2433  muted #5c6b7e  border #d2dae4
_STYLE = """
QMainWindow, QWidget {
    background: #eef1f5; color: #1a2433; font-size: 12px;
    font-family: "Helvetica Neue", "Segoe UI", "Arial", sans-serif;
}

/* Title band */
QWidget#header { background: #1b2b44; border-bottom: 2px solid #2f5e8f; }
QLabel { background: transparent; color: #1a2433; }
QLabel#title { color: #ffffff; font-size: 18px; font-weight: 700; }
QLabel#subtitle { color: #9fb3cc; font-size: 12px; }
QLabel#statusbar { background: #e3e9f1; color: #4a5a70;
                   border-top: 1px solid #cdd6e1; padding: 5px 12px; font-size: 11px; }

/* Menu bar */
QMenuBar { background: #e8edf3; color: #27374b; border-bottom: 1px solid #d4dce6; padding: 2px 6px; }
QMenuBar::item { background: transparent; padding: 5px 12px; border-radius: 4px; }
QMenuBar::item:selected { background: #d6e3f2; color: #16324f; }
QMenu { background: #ffffff; color: #1a2433; border: 1px solid #cdd6e1; padding: 4px; }
QMenu::item { padding: 6px 26px 6px 18px; border-radius: 4px; }
QMenu::item:selected { background: #2f5e8f; color: #ffffff; }
QMenu::item:disabled { color: #8a98a8; }
QMenu::separator { height: 1px; background: #dde3ea; margin: 4px 8px; }

/* Cards */
QGroupBox { font-weight: 600; color: #27374b; border: 1px solid #d2dae4; border-radius: 8px;
            margin-top: 12px; padding: 12px 10px 10px 10px; background: #ffffff; }
QGroupBox::title { subcontrol-origin: margin; left: 12px; padding: 0 6px; color: #2f5e8f; }

/* Buttons */
QPushButton { background: #2f5e8f; color: #ffffff; border: 1px solid #2a557f;
              border-radius: 6px; padding: 6px 14px; font-weight: 600; }
QPushButton:hover { background: #274f79; }
QPushButton:pressed { background: #1f3f61; }
QPushButton:disabled { background: #c2cedd; color: #eef1f5; border-color: #c2cedd; }
QPushButton:focus { border: 1px solid #7da3c8; }

/* Combo boxes */
QComboBox { background: #ffffff; color: #1a2433; border: 1px solid #cdd6e1;
            border-radius: 6px; padding: 4px 8px; min-height: 22px; }
QComboBox:hover { border-color: #9bb4d0; }
QComboBox:focus { border-color: #2f5e8f; }
QComboBox::drop-down { border: none; width: 22px; }
QComboBox::down-arrow { image: none; width: 0; height: 0; margin-right: 9px;
    border-left: 4px solid transparent; border-right: 4px solid transparent;
    border-top: 5px solid #5c6b7e; }
QComboBox QAbstractItemView { background: #ffffff; color: #1a2433; border: 1px solid #cdd6e1;
    selection-background-color: #2f5e8f; selection-color: #ffffff; outline: none; }

/* Spin boxes */
QDoubleSpinBox, QSpinBox { background: #ffffff; color: #1a2433; border: 1px solid #cdd6e1;
    border-radius: 6px; padding: 3px 6px; min-height: 22px; }
QDoubleSpinBox:focus, QSpinBox:focus { border-color: #2f5e8f; }

/* Flux-range slider */
QSlider::groove:horizontal { height: 5px; background: #d2dae4; border-radius: 3px; }
QSlider::sub-page:horizontal { background: #2f5e8f; border-radius: 3px; }
QSlider::handle:horizontal { background: #ffffff; border: 2px solid #2f5e8f;
    width: 14px; height: 14px; margin: -6px 0; border-radius: 8px; }
QSlider::handle:horizontal:hover { border-color: #274f79; }

/* Tables */
QTableWidget { background: #ffffff; alternate-background-color: #f5f8fb; color: #1a2433;
    gridline-color: #e4e9f0; border: 1px solid #d2dae4; border-radius: 6px;
    selection-background-color: #d6e3f2; selection-color: #16324f; }
QTableWidget::item { padding: 3px 6px; }
QHeaderView::section { background: #e8edf3; color: #34465c; padding: 6px 8px; border: none;
    border-right: 1px solid #d8dfe8; border-bottom: 1px solid #d2dae4; font-weight: 600; }

/* Tabs */
QTabWidget::pane { border: 1px solid #d2dae4; border-radius: 6px; top: -1px; background: #ffffff; }
QTabBar::tab { background: #e6ebf2; color: #54637a; padding: 8px 18px; margin-right: 2px;
    border: 1px solid #d8dfe8; border-bottom: none;
    border-top-left-radius: 6px; border-top-right-radius: 6px; }
QTabBar::tab:hover { background: #eef2f7; color: #2f5e8f; }
QTabBar::tab:selected { background: #ffffff; color: #1b2b44; font-weight: 600;
    border-bottom: 2px solid #2f5e8f; margin-bottom: -1px; }

/* Scrollbars */
QScrollBar:vertical { background: transparent; width: 11px; margin: 0; }
QScrollBar::handle:vertical { background: #c2cedd; border-radius: 5px; min-height: 28px; }
QScrollBar::handle:vertical:hover { background: #a9b9cc; }
QScrollBar:horizontal { background: transparent; height: 11px; margin: 0; }
QScrollBar::handle:horizontal { background: #c2cedd; border-radius: 5px; min-width: 28px; }
QScrollBar::handle:horizontal:hover { background: #a9b9cc; }
QScrollBar::add-line, QScrollBar::sub-line { height: 0; width: 0; }
QScrollBar::add-page, QScrollBar::sub-page { background: transparent; }

/* Splitter + tooltips */
QSplitter::handle { background: #dde3ea; }
QSplitter::handle:horizontal { width: 6px; }
QSplitter::handle:hover { background: #c2cedd; }
QToolTip { background: #1b2b44; color: #ffffff; border: 1px solid #12203a;
    padding: 5px 8px; border-radius: 4px; }
"""


def _read_expression_vector(path: str) -> dict[str, float]:
    """Read a two-column gene-expression CSV/TSV file into ``gene -> value``."""

    import pandas as pd

    frame = pd.read_csv(path, sep=None, engine="python")
    if frame.shape[1] < 2:
        raise ValueError("expression file must have at least two columns: gene, expression")
    gene_col, value_col = frame.columns[:2]
    expression: dict[str, float] = {}
    for gene, value in zip(frame[gene_col], frame[value_col], strict=False):
        if pd.isna(gene) or pd.isna(value):
            continue
        expression[str(gene)] = float(value)
    if not expression:
        raise ValueError("expression file contains no numeric expression values")
    return expression


class _AnalysisWorker(QObject):
    """Runs one compute callable on a worker thread and reports back to the UI thread.

    The callable must be pure computation over the cobra model / services (no Qt or matplotlib
    access) — those stay on the main thread. The result (or the raised exception) is captured
    so the caller can render it (or re-raise) once ``finished`` fires.
    """

    finished = Signal()

    def __init__(self, fn):
        super().__init__()
        self._fn = fn
        self.result = None
        self.error: Exception | None = None

    def run(self) -> None:
        try:
            self.result = self._fn()
        except Exception as exc:  # captured and re-raised on the UI thread by the caller
            self.error = exc
        finally:
            self.finished.emit()


class CmmMainWindow(QMainWindow):
    """Main platform window over a single cobra model."""

    def __init__(self, model: Model, parent=None, *, default_product: str | None = None,
                 map_path: str | None = None):
        super().__init__(parent)
        self.model = model
        self._fluxes: dict[str, float] = {}
        self._loading = False
        self._fluxes_stale = False
        self._default_product = default_product
        self._map_path = map_path
        self._revert_source_expression: dict[str, float] | None = None
        self._revert_target_expression: dict[str, float] | None = None
        # Real gene expression loaded from a CSV on the Omics tab, reused as the LAD/E-Flux2
        # reference template on the Comparison tab so that path is not left to synthetic data.
        self._omics_expression: dict[str, float] | None = None
        # A→B transformation finder: two expression vectors → two predicted flux states.
        self._transform_source_expression: dict[str, float] | None = None
        self._transform_target_expression: dict[str, float] | None = None
        # Multi-condition comparison: a gene × condition expression table.
        self._condition_table = None
        # Guards re-entrant background runs (see _run_in_background).
        self._busy = False
        self.setWindowTitle("CMM — Cellular Metabolic Modeling Platform")
        self.resize(1180, 740)
        self.setStyleSheet(_STYLE)
        self._build_ui()
        self.load_model(model)

    # -- construction -------------------------------------------------------

    def _build_ui(self) -> None:
        self._build_menus()
        root = QWidget()
        self.setCentralWidget(root)
        outer = QVBoxLayout(root)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        outer.addWidget(self._build_header())

        body = QWidget()
        body_layout = QHBoxLayout(body)
        body_layout.setContentsMargins(12, 12, 12, 12)
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self._build_model_panel())
        splitter.addWidget(self._build_tabs())
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 5)
        body_layout.addWidget(splitter)
        outer.addWidget(body, 1)

        self.status_label = QLabel("Ready.")
        self.status_label.setObjectName("statusbar")
        outer.addWidget(self.status_label)

    def _build_menus(self) -> None:
        """A desktop-style menu bar: Analysis / Model / Config drive the platform."""

        bar = self.menuBar()

        file_menu = bar.addMenu("&File")
        file_menu.addAction("Open Model…", self.open_model_dialog)
        file_menu.addSeparator()
        file_menu.addAction("Export Table to CSV…", self.export_table_csv)
        file_menu.addAction("Save Figure…", self.save_active_figure)
        file_menu.addSeparator()
        file_menu.addAction("Quit", self.close)

        analysis = bar.addMenu("&Analysis")
        analysis.addAction("Run FBA", self.run_fba)
        analysis.addAction("Run FVA", self.run_fva)
        analysis.addSeparator()
        analysis.addAction("Theoretical Yield",
                           lambda: self._in_tab("Production", self.run_theoretical_yield))
        analysis.addAction("Production Envelope",
                           lambda: self._in_tab("Production", self.run_production_envelope_plot))
        analysis.addAction("FSEOF Targets",
                           lambda: self._in_tab("Production", self.run_fseof_plot))
        analysis.addAction("Strain Design (OptKnock / RobustKnock)…",
                           lambda: self._goto_tab("Strain Design"))
        analysis.addSeparator()
        analysis.addAction("Omics Integration (E-Flux2 / LAD)…",
                           lambda: self._goto_tab("Omics"))
        analysis.addAction("Multi-condition Comparison…",
                           lambda: self._goto_tab("Multi-condition"))
        analysis.addAction("Revert Metabolism…", lambda: self._goto_tab("Revert Metabolism"))
        analysis.addAction("Transformation Targets (A→B)…",
                           lambda: self._goto_tab("Transform (A→B)"))
        if self._map_path:
            analysis.addSeparator()
            analysis.addAction("Render Flux Map",
                               lambda: self._in_tab("Flux Map", self.render_flux_map))

        model_menu = bar.addMenu("&Model")
        model_menu.addAction("Model Info…", self._show_model_info)

        config = bar.addMenu("&Config")
        status = solver_status(self.model)
        solver_action = config.addAction(f"Solver: {status.summary()}")
        solver_action.setEnabled(False)
        config.addAction("Solver status…", self._show_solver_status)

    def _tab_index(self, name: str) -> int | None:
        return next(
            (i for i in range(self.tabs.count()) if self.tabs.tabText(i) == name), None
        )

    def _goto_tab(self, name: str) -> None:
        index = self._tab_index(name)
        if index is not None:
            self.tabs.setCurrentIndex(index)

    def _in_tab(self, name: str, action) -> None:
        self._goto_tab(name)
        action()

    def _show_solver_status(self) -> None:
        status = solver_status(self.model)
        lines = [
            f"Active solver: {status.name}",
            f"Capabilities: {', '.join(status.capabilities)}",
            f"Recommended (full QP+MILP+MIQP): {'yes' if status.recommended else 'no'}",
            f"Available solvers: {', '.join(status.available)}",
        ]
        if status.warning:
            lines.append("")
            lines.append("⚠ " + status.warning)
            QMessageBox.warning(self, "Solver status", "\n".join(lines))
        else:
            QMessageBox.information(self, "Solver status", "\n".join(lines))

    def _show_model_info(self) -> None:
        objective = ", ".join(
            r.id for r in self.model.reactions if r.objective_coefficient != 0
        )
        QMessageBox.information(
            self, "Model Info",
            f"{self.model.id}\n\n"
            f"Reactions: {len(self.model.reactions)}\n"
            f"Metabolites: {len(self.model.metabolites)}\n"
            f"Genes: {len(self.model.genes)}\n"
            f"Objective: {objective or '—'} ({self.model.objective_direction})\n"
            f"Solver: {active_solver(self.model)}",
        )

    # -- background execution -----------------------------------------------

    def _run_in_background(self, compute, *, label: str = "Working…"):
        """Run ``compute()`` off the UI thread and return its result synchronously.

        The heavy solve runs on a worker thread while a nested event loop keeps the window
        repainting; a modal busy dialog blocks input so the model is not mutated mid-solve.
        Because the nested loop only exits when the worker finishes, callers still see a plain
        synchronous call (tests and the scenario harnesses read results immediately after).

        Re-entrant calls (a backgrounded handler invoking another) run inline so busy dialogs
        and worker threads never stack.
        """

        if self._busy:  # already inside a background run — don't nest threads/dialogs
            return compute()

        self._busy = True
        worker = _AnalysisWorker(compute)
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)

        loop = QEventLoop()
        worker.finished.connect(loop.quit)

        # A modal, no-cancel busy dialog: blocks input (menus, table edits, other analyses) so
        # only the worker touches the model, and auto-shows only if the work exceeds
        # minimumDuration so fast solves do not flash a dialog.
        dialog = QProgressDialog(label, None, 0, 0, self)
        dialog.setWindowTitle("CMM")
        dialog.setWindowModality(Qt.ApplicationModal)
        dialog.setMinimumDuration(400)
        dialog.setAutoClose(False)
        dialog.setAutoReset(False)

        # Start the thread only once the event loop is running, so the worker cannot finish and
        # queue loop.quit before exec_() begins (which would otherwise hang the loop).
        QTimer.singleShot(0, thread.start)
        try:
            loop.exec_()
        finally:
            dialog.reset()
            dialog.deleteLater()
            thread.quit()
            thread.wait()
            worker.deleteLater()
            self._busy = False

        if worker.error is not None:
            raise worker.error
        return worker.result

    # -- File menu (open / export) ------------------------------------------

    def open_model_dialog(self) -> None:
        """Load a different SBML model into the running window."""

        path, _ = QFileDialog.getOpenFileName(
            self, "Open SBML model", "",
            "SBML models (*.xml *.sbml *.xml.gz *.json);;All files (*)",
        )
        if not path:
            return
        try:
            import cobra

            model = (
                cobra.io.load_json_model(path)
                if path.endswith(".json")
                else cobra.io.read_sbml_model(path)
            )
        except Exception as exc:
            self.status_label.setText(f"Could not load model: {html.escape(str(exc))}")
            return
        # A new model invalidates any loaded expression / revert state.
        self._omics_expression = None
        self._revert_source_expression = None
        self._revert_target_expression = None
        self.revert_source_label.setText("not loaded")
        self.revert_target_label.setText("not loaded")
        self._update_revert_run_state()
        self.load_model(model)
        self.status_label.setText(f"Loaded model '{model.id}' from {path}.")

    def _active_table(self) -> QTableWidget | None:
        """The primary result table of the current tab, or None for figure-only tabs."""

        name = self.tabs.tabText(self.tabs.currentIndex())
        return {
            "Simulation": self.sim_table,
            "Comparison": self.comparison_table,
            "Strain Design": self.sd_table,
            "Omics": self.omics_table,
            "Multi-condition": self.cond_table,
            "Revert Metabolism": self.revert_table,
            "Transform (A→B)": self.transform_table,
        }.get(name)

    def _active_figure(self):
        """The matplotlib figure of the current tab, or None if it has no figure."""

        name = self.tabs.tabText(self.tabs.currentIndex())
        if name == "Production" and self._production_canvas is not None:
            return self._production_canvas.figure
        if name == "Flux Map" and self._map_canvas is not None:
            return self._map_canvas.figure
        return None

    def export_table_csv(self) -> None:
        """Write the current tab's result table to a CSV file."""

        import csv

        table = self._active_table()
        if table is None or table.rowCount() == 0:
            self.status_label.setText(
                "No result table to export on this tab — run an analysis first."
            )
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export table to CSV", "", "CSV files (*.csv);;All files (*)"
        )
        if not path:
            return
        headers = [
            table.horizontalHeaderItem(c).text() if table.horizontalHeaderItem(c) else f"col{c}"
            for c in range(table.columnCount())
        ]
        try:
            with open(path, "w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(headers)
                for row in range(table.rowCount()):
                    writer.writerow(
                        [
                            table.item(row, col).text() if table.item(row, col) else ""
                            for col in range(table.columnCount())
                        ]
                    )
        except OSError as exc:
            self.status_label.setText(f"Could not write CSV: {exc}")
            return
        self.status_label.setText(f"Exported {table.rowCount()} rows to {path}.")

    def save_active_figure(self) -> None:
        """Save the current tab's figure (Production or Flux Map) at 300 DPI."""

        fig = self._active_figure()
        if fig is None:
            self.status_label.setText(
                "No figure to save — open the Production or Flux Map tab and run an analysis."
            )
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save figure", "", "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)"
        )
        if not path:
            return
        try:
            fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
        except Exception as exc:
            self.status_label.setText(f"Could not save figure: {html.escape(str(exc))}")
            return
        self.status_label.setText(f"Saved figure to {path}.")

    def _build_header(self) -> QWidget:
        header = QWidget()
        header.setObjectName("header")
        layout = QVBoxLayout(header)
        layout.setContentsMargins(16, 12, 16, 12)
        title = QLabel("CMM — Cellular Metabolic Modeling Platform")
        title.setObjectName("title")
        self.subtitle = QLabel("")
        self.subtitle.setObjectName("subtitle")
        layout.addWidget(title)
        layout.addWidget(self.subtitle)
        return header

    def _build_model_panel(self) -> QWidget:
        box = QGroupBox("Model")
        layout = QVBoxLayout(box)
        self.model_info = QLabel("")
        self.model_info.setWordWrap(True)
        layout.addWidget(self.model_info)

        slider_row = QHBoxLayout()
        slider_row.addWidget(QLabel("Flux range filter:"))
        self.flux_slider = QSlider(Qt.Horizontal)
        self.flux_slider.setRange(0, 100)
        self.flux_slider.setToolTip(
            "Highlight reactions whose |flux| is at least this fraction of the maximum "
            "flux (run FBA first)."
        )
        self.flux_slider.valueChanged.connect(self._apply_flux_threshold)
        self.flux_threshold_label = QLabel("0%")
        slider_row.addWidget(self.flux_slider, 1)
        slider_row.addWidget(self.flux_threshold_label)
        layout.addLayout(slider_row)

        self.reaction_table = QTableWidget(0, 4)
        self.reaction_table.setHorizontalHeaderLabels(["Reaction", "Lower", "Upper", "Flux"])
        # Reaction IDs vary in length (short on toy models, long like "Biomass_Ecoli_core" on
        # genome-scale ones); let the id column take the slack while numeric columns fit content.
        _rxn_header = self.reaction_table.horizontalHeader()
        _rxn_header.setSectionResizeMode(0, QHeaderView.Stretch)
        for _col in (1, 2, 3):
            _rxn_header.setSectionResizeMode(_col, QHeaderView.ResizeToContents)
        self.reaction_table.verticalHeader().setVisible(False)
        self.reaction_table.setAlternatingRowColors(True)
        self.reaction_table.itemChanged.connect(self._on_bound_edited)
        layout.addWidget(self.reaction_table, 1)
        hint = QLabel("Double-click a Lower/Upper cell to edit a bound, then re-run FBA.")
        hint.setStyleSheet("color: #5c6b7e; font-style: italic;")
        layout.addWidget(hint)
        return box

    def _build_tabs(self) -> QWidget:
        self.tabs = QTabWidget()
        self.tabs.addTab(self._build_simulation_tab(), "Simulation")
        self.tabs.addTab(self._build_comparison_tab(), "Comparison")
        self.tabs.addTab(self._build_production_tab(), "Production")
        self.tabs.addTab(self._build_strain_design_tab(), "Strain Design")
        self.tabs.addTab(self._build_omics_tab(), "Omics")
        self.tabs.addTab(self._build_conditions_tab(), "Multi-condition")
        if self._map_path:
            self.tabs.addTab(self._build_fluxmap_tab(), "Flux Map")
        self.tabs.addTab(self._build_revert_tab(), "Revert Metabolism")
        self.tabs.addTab(self._build_transformation_tab(), "Transform (A→B)")
        return self.tabs

    def _build_omics_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Omics integration (expression → flux)")
        row = QHBoxLayout(controls)
        row.addWidget(QLabel("Method:"))
        self.omics_method_combo = QComboBox()
        self.omics_method_combo.addItems(["eflux2", "lad"])
        self.omics_method_combo.setToolTip(
            "E-Flux2: scale bounds by expression, minimize L2 flux. LAD: fit fluxes to expression."
        )
        row.addWidget(self.omics_method_combo)
        load_btn = QPushButton("Load expression CSV…")
        load_btn.clicked.connect(self.load_expression_csv)
        demo_btn = QPushButton("Run on demo expression")
        demo_btn.clicked.connect(self.run_omics_demo)
        row.addWidget(load_btn)
        row.addWidget(demo_btn)
        row.addStretch(1)
        layout.addWidget(controls)

        self.omics_summary = QLabel(
            "Load a gene-expression CSV (columns: gene, expression) or run on demo data."
        )
        self.omics_summary.setWordWrap(True)
        layout.addWidget(self.omics_summary)

        self.omics_table = QTableWidget(0, 2)
        self.omics_table.setHorizontalHeaderLabels(["Reaction", "Predicted flux"])
        self.omics_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.omics_table.verticalHeader().setVisible(False)
        self.omics_table.setAlternatingRowColors(True)
        layout.addWidget(self.omics_table, 1)
        return tab

    def _build_fluxmap_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)
        controls = QHBoxLayout()
        render_btn = QPushButton("Render flux map")
        render_btn.clicked.connect(self.render_flux_map)
        controls.addWidget(render_btn)
        controls.addWidget(QLabel("Escher layout coloured by the current FBA flux."))
        controls.addStretch(1)
        layout.addLayout(controls)
        self.map_canvas_holder = QVBoxLayout()
        holder = QWidget()
        holder.setLayout(self.map_canvas_holder)
        layout.addWidget(holder, 1)
        self._map_canvas = None
        self._map_toolbar = None
        return tab

    def _build_production_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Production design")
        form = QHBoxLayout(controls)
        form.addWidget(QLabel("Target product:"))
        self.product_combo = QComboBox()
        self.product_combo.setMinimumWidth(160)
        self.product_combo.setToolTip("Exchange reaction of the target product")
        self.product_combo.currentTextChanged.connect(self._on_product_changed)
        form.addWidget(self.product_combo, 1)
        form.addWidget(QLabel("substrate:"))
        self.substrate_combo = QComboBox()
        self.substrate_combo.setMinimumWidth(140)
        self.substrate_combo.setToolTip("Carbon substrate (auto = detect from the medium)")
        form.addWidget(self.substrate_combo, 1)
        self.anaerobic_combo = QComboBox()
        self.anaerobic_combo.addItems(["aerobic", "anaerobic"])
        form.addWidget(self.anaerobic_combo)
        prod_fba_btn = QPushButton("Run FBA")
        prod_fba_btn.clicked.connect(self.run_fba)
        yield_btn = QPushButton("Theoretical yield")
        yield_btn.clicked.connect(self.run_theoretical_yield)
        envelope_btn = QPushButton("Production envelope")
        envelope_btn.clicked.connect(self.run_production_envelope_plot)
        fseof_btn = QPushButton("FSEOF targets")
        fseof_btn.clicked.connect(self.run_fseof_plot)
        fvseof_btn = QPushButton("FVSEOF (robust)")
        fvseof_btn.clicked.connect(self.run_fvseof_plot)
        form.addWidget(prod_fba_btn)
        form.addWidget(yield_btn)
        form.addWidget(envelope_btn)
        form.addWidget(fseof_btn)
        form.addWidget(fvseof_btn)
        self._production_buttons = [yield_btn, envelope_btn, fseof_btn, fvseof_btn]
        layout.addWidget(controls)

        self.yield_label = QLabel("Select a product and compute its theoretical yield.")
        self.yield_label.setFont(QFont("", 13, QFont.Bold))
        self.yield_label.setWordWrap(True)
        self.yield_label.setTextFormat(Qt.RichText)
        layout.addWidget(self.yield_label)

        self.production_canvas_holder = QVBoxLayout()
        holder = QWidget()
        holder.setLayout(self.production_canvas_holder)
        layout.addWidget(holder, 1)
        self._production_canvas = None
        self._production_toolbar = None
        return tab

    def _build_strain_design_tab(self) -> QWidget:
        """Growth-coupled strain design: OptKnock / RobustKnock over the library service."""

        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Growth-coupled strain design (OptKnock / RobustKnock)")
        form = QHBoxLayout(controls)
        form.addWidget(QLabel("Target product:"))
        self.sd_product_combo = QComboBox()
        self.sd_product_combo.setMinimumWidth(150)
        self.sd_product_combo.setToolTip("Exchange reaction of the target product to couple to growth")
        form.addWidget(self.sd_product_combo, 1)
        form.addWidget(QLabel("Method:"))
        self.sd_method_combo = QComboBox()
        self.sd_method_combo.addItems(["optknock", "robustknock"])
        self.sd_method_combo.setToolTip(
            "OptKnock maximizes product at max growth (optimistic); RobustKnock keeps only "
            "designs that guarantee product at max growth (worst case)."
        )
        form.addWidget(self.sd_method_combo)
        form.addWidget(QLabel("max KOs:"))
        self.sd_max_ko_spin = QSpinBox()
        self.sd_max_ko_spin.setRange(1, 6)
        self.sd_max_ko_spin.setValue(3)
        self.sd_max_ko_spin.setToolTip("Maximum number of reaction knockouts per design")
        form.addWidget(self.sd_max_ko_spin)
        form.addWidget(QLabel("solutions:"))
        self.sd_max_sol_spin = QSpinBox()
        self.sd_max_sol_spin.setRange(1, 20)
        self.sd_max_sol_spin.setValue(5)
        self.sd_max_sol_spin.setToolTip("Maximum number of designs to enumerate")
        form.addWidget(self.sd_max_sol_spin)
        self.sd_run_btn = QPushButton("Run design")
        self.sd_run_btn.clicked.connect(self.run_strain_design)
        form.addWidget(self.sd_run_btn)
        form.addStretch(1)
        layout.addWidget(controls)

        self.sd_summary = QLabel(
            "OptKnock/RobustKnock search knockout sets that couple product to growth. "
            "Needs a MILP solver (Gurobi/CPLEX)."
        )
        self.sd_summary.setWordWrap(True)
        layout.addWidget(self.sd_summary)

        self.sd_table = QTableWidget(0, 5)
        self.sd_table.setHorizontalHeaderLabels(
            ["Rank", "Knockouts", "Growth", "Max product", "Guaranteed product"]
        )
        _sd_header = self.sd_table.horizontalHeader()
        _sd_header.setSectionResizeMode(1, QHeaderView.Stretch)
        for _col in (0, 2, 3, 4):
            _sd_header.setSectionResizeMode(_col, QHeaderView.ResizeToContents)
        self.sd_table.verticalHeader().setVisible(False)
        self.sd_table.setAlternatingRowColors(True)
        layout.addWidget(self.sd_table, 1)
        return tab

    def run_strain_design(self) -> None:
        product = self.sd_product_combo.currentText()
        if not product:
            self.sd_summary.setText("Select a target product exchange reaction first.")
            return
        method = self.sd_method_combo.currentText()
        max_ko = self.sd_max_ko_spin.value()
        max_sol = self.sd_max_sol_spin.value()
        self.sd_summary.setText(f"Searching {method} designs for {html.escape(product)}…")
        self.status_label.setText(f"Running {method} (this can take a while)…")
        self.sd_run_btn.setEnabled(False)
        try:
            search = optknock if method == "optknock" else robustknock
            result = self._run_in_background(
                lambda: search(self.model, product, max_knockouts=max_ko, max_solutions=max_sol),
                label=f"Running {method}…",
            )
        except Exception as exc:  # surface solver-capability/model errors instead of freezing
            self.sd_table.setRowCount(0)
            self.sd_summary.setText(f"Strain design failed: {html.escape(str(exc))}")
            self.status_label.setText(f"Strain design failed ({method}).")
            return
        finally:
            self.sd_run_btn.setEnabled(True)

        designs = result.designs
        self.sd_table.setRowCount(len(designs))
        for i, design in enumerate(designs):
            self.sd_table.setItem(i, 0, QTableWidgetItem(str(i + 1)))
            self.sd_table.setItem(i, 1, QTableWidgetItem(", ".join(design.knockouts) or "—"))
            self.sd_table.setItem(i, 2, QTableWidgetItem(f"{design.growth:.4g}"))
            self.sd_table.setItem(i, 3, QTableWidgetItem(f"{design.max_product:.4g}"))
            self.sd_table.setItem(i, 4, QTableWidgetItem(f"{design.guaranteed_product:.4g}"))
            if design.growth_coupled:
                self.sd_table.item(i, 4).setBackground(QColor("#dce9d6"))
        if designs:
            best = result.best()
            self.sd_summary.setText(
                f"<b>{method}</b> for {html.escape(product)}: {len(designs)} design(s). "
                f"Best knocks out <b>{html.escape(', '.join(best.knockouts) or '—')}</b> "
                f"(growth {best.growth:.3g} h⁻¹, guaranteed product {best.guaranteed_product:.3g}, "
                f"max product {best.max_product:.3g})."
            )
        else:
            self.sd_summary.setText(
                f"<b>{method}</b> found no growth-coupled design for "
                f"{html.escape(product)} — try a different product, more knockouts, or check "
                "the product can carry flux in the current medium."
            )
        self.status_label.setText(f"Strain design complete ({method}).")

    def _build_comparison_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)
        controls = QGroupBox("Perturbation response (MOMA / ROOM vs a reference template)")
        controls_layout = QVBoxLayout(controls)

        row = QHBoxLayout()
        row.addWidget(QLabel("Method:"))
        self.comparison_method_combo = QComboBox()
        self.comparison_method_combo.addItems(["MOMA (L2)", "MOMA (L1)", "ROOM"])
        row.addWidget(self.comparison_method_combo)
        row.addWidget(QLabel("Reference template:"))
        self.template_combo = QComboBox()
        self.template_combo.addItems(["fba", "pfba", "lad", "eflux2"])
        self.template_combo.setToolTip(
            "Wild-type flux template the perturbed state is compared to. "
            "LAD/E-Flux2 use expression loaded on the Omics tab (synthetic demo data if none)."
        )
        row.addWidget(self.template_combo)
        row.addWidget(QLabel("Knockout level:"))
        self.ko_level_combo = QComboBox()
        self.ko_level_combo.addItems(["reaction", "gene"])
        self.ko_level_combo.setToolTip("Knock out reactions, or genes (resolved to reactions via GPR)")
        self.ko_level_combo.currentTextChanged.connect(self._populate_ko_list)
        row.addWidget(self.ko_level_combo)
        row.addStretch(1)
        controls_layout.addLayout(row)

        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Knock out (select one or more):"))
        run_btn = QPushButton("Run (selected as one KO)")
        run_btn.setToolTip("MOMA/ROOM with all selected targets knocked out together")
        run_btn.clicked.connect(self.run_comparison)
        batch_btn = QPushButton("Batch (each separately)")
        batch_btn.setToolTip(
            "Run MOMA/ROOM once per target as a separate single knockout "
            "(all targets of this level if none selected)."
        )
        batch_btn.clicked.connect(self.run_batch_comparison)
        row2.addWidget(run_btn)
        row2.addWidget(batch_btn)
        row2.addStretch(1)
        controls_layout.addLayout(row2)

        self.ko_list = QListWidget()
        self.ko_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.ko_list.setMaximumHeight(120)
        self.ko_list.setToolTip(
            "Ctrl/Shift-click to select multiple targets. Selecting several then 'Run' knocks "
            "them out together; 'Batch' runs each on its own."
        )
        controls_layout.addWidget(self.ko_list)
        layout.addWidget(controls)

        self.comparison_summary = QLabel(
            "Pick a reference template (FBA/pFBA/LAD/E-Flux2), a knockout level (reaction/gene), "
            "one or more targets, and a method to predict the perturbed flux state."
        )
        self.comparison_summary.setWordWrap(True)
        layout.addWidget(self.comparison_summary)

        self.comparison_table = QTableWidget(0, 3)
        self.comparison_table.setHorizontalHeaderLabels(
            ["Reaction", "Reference flux", "Perturbed flux"]
        )
        self.comparison_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.comparison_table.verticalHeader().setVisible(False)
        self.comparison_table.setAlternatingRowColors(True)
        layout.addWidget(self.comparison_table, 1)
        return tab

    def _populate_ko_list(self) -> None:
        """Fill the knockout list with the model's genes or reactions for the current level."""

        if not hasattr(self, "ko_list"):
            return
        level = self.ko_level_combo.currentText()
        self.ko_list.clear()
        if level == "gene":
            self.ko_list.addItems([g.id for g in self.model.genes])
        else:
            self.ko_list.addItems([r.id for r in self.model.reactions])

    def _selected_ko_targets(self) -> list[str]:
        return [item.text() for item in self.ko_list.selectedItems()]

    def _comparison_method_key(self) -> str:
        label = self.comparison_method_combo.currentText()
        if label == "ROOM":
            return "room"
        return "moma_l1" if label == "MOMA (L1)" else "moma_l2"

    def _comparison_expression(self):
        """Expression for a LAD/E-Flux2 template: real if loaded, else synthetic (flagged)."""

        template = self.template_combo.currentText()
        if template not in ("lad", "eflux2"):
            return None, False
        if self._omics_expression is not None:
            return self._omics_expression, False
        import numpy as np

        rng = np.random.default_rng(0)
        return {g.id: float(rng.uniform(1, 100)) for g in self.model.genes}, True

    @staticmethod
    def _ko_label(targets: list[str]) -> str:
        shown = ", ".join(targets[:4])
        return shown if len(targets) <= 4 else f"{shown}, +{len(targets) - 4} more"

    def run_comparison(self) -> None:
        method_label = self.comparison_method_combo.currentText()
        template = self.template_combo.currentText()
        level = self.ko_level_combo.currentText()
        targets = self._selected_ko_targets()
        if not targets:
            self.comparison_summary.setText(
                f"Select one or more {level}s to knock out (Ctrl/Shift-click for several)."
            )
            return
        method_key = self._comparison_method_key()
        expression, synthetic_expression = self._comparison_expression()
        ko_label = self._ko_label(targets)

        # Restore the single/multi-KO table shape (a prior batch run reshapes it).
        self.comparison_table.setColumnCount(3)
        self.comparison_table.setHorizontalHeaderLabels(
            ["Reaction", "Reference flux", "Perturbed flux"]
        )

        def _compute():
            reference = reference_flux(self.model, template, gene_expression=expression)
            if level == "gene":
                reaction_ids = blocked_reactions_for_genes(self.model, targets)
            else:
                reaction_ids = tuple(targets)
            result = knockout_comparison(self.model, reference, reaction_ids, method=method_key)
            return reference, result, reaction_ids

        try:
            reference, result, reaction_ids = self._run_in_background(
                _compute, label=f"Running {method_label}…"
            )
        except Exception as exc:
            self.comparison_summary.setText(f"Comparison failed: {html.escape(str(exc))}")
            return

        synthetic_note = (
            " <span style='color:#b45309'>⚠ synthetic demo expression — load a CSV on the "
            "Omics tab for a data-driven template.</span>"
            if synthetic_expression else ""
        )
        blocked_note = (
            f" ({len(reaction_ids)} reactions blocked)" if level == "gene" else ""
        )
        if level == "gene" and not reaction_ids:
            self.comparison_table.setRowCount(0)
            self.comparison_summary.setText(
                f"Knocking out gene(s) {html.escape(ko_label)} blocks no reactions "
                "(no GPR effect) — nothing to compare." + synthetic_note
            )
            self.status_label.setText("Comparison: knockout has no effect.")
            return
        if result.status != "optimal":
            self.comparison_table.setRowCount(0)
            self.comparison_summary.setText(
                f"<b>{method_label}</b> vs <b>{template}</b>: knocking out "
                f"{html.escape(ko_label)}{blocked_note} makes the model <b>infeasible</b> "
                "(e.g. it cannot meet maintenance demand)." + synthetic_note
            )
            self.status_label.setText("Comparison: infeasible perturbation.")
            return
        rows = sorted(result.fluxes.items(), key=lambda kv: -abs(kv[1] - reference.get(kv[0])))
        changed = [(r, v) for r, v in rows if abs(v - reference.get(r)) > 1e-6]
        self.comparison_table.setRowCount(len(changed))
        for i, (rid, flux) in enumerate(changed):
            self.comparison_table.setItem(i, 0, QTableWidgetItem(rid))
            self.comparison_table.setItem(i, 1, QTableWidgetItem(f"{reference.get(rid):.4g}"))
            self.comparison_table.setItem(i, 2, QTableWidgetItem(f"{flux:.4g}"))
        self.comparison_summary.setText(
            f"<b>{method_label}</b> vs <b>{template}</b> template after knocking out "
            f"{level} {html.escape(ko_label)}{blocked_note}: status {result.status}, "
            f"distance {result.distance:.4g}, {len(changed)} reactions changed." + synthetic_note
        )
        self.status_label.setText(f"Comparison complete ({method_label}, {template}).")

    def run_batch_comparison(self) -> None:
        """Run MOMA/ROOM once per target (each a separate single knockout) into a batch table."""

        method_label = self.comparison_method_combo.currentText()
        template = self.template_combo.currentText()
        level = self.ko_level_combo.currentText()
        selected = self._selected_ko_targets()  # empty -> all targets of this level
        method_key = self._comparison_method_key()
        expression, synthetic_expression = self._comparison_expression()

        def _compute():
            reference = reference_flux(self.model, template, gene_expression=expression)
            if level == "gene":
                perts = gene_perturbations(self.model, selected or None)
            else:
                perts = reaction_perturbations(self.model, selected or None)
            return batch_comparison(self.model, reference, perts, method=method_key), len(perts)

        self.status_label.setText(f"Running batch {method_label} ({level})…")
        try:
            rows, n = self._run_in_background(
                _compute, label=f"Batch {method_label} over {level}s…"
            )
        except Exception as exc:
            self.comparison_summary.setText(f"Batch comparison failed: {html.escape(str(exc))}")
            self.status_label.setText("Batch comparison failed.")
            return

        # Most-disrupted first; infeasible (NaN distance) sorted last.
        def _sort_key(r):
            import math

            return (0, -r.distance) if math.isfinite(r.distance) else (1, 0.0)

        rows = sorted(rows, key=_sort_key)
        self.comparison_table.setColumnCount(6)
        self.comparison_table.setHorizontalHeaderLabels(
            ["Target", "Kind", "#reactions", "Status", "Distance", "Objective"]
        )
        self.comparison_table.setRowCount(len(rows))
        for i, r in enumerate(rows):
            dist = "infeasible" if r.status != "optimal" else f"{r.distance:.4g}"
            self.comparison_table.setItem(i, 0, QTableWidgetItem(r.target_id))
            self.comparison_table.setItem(i, 1, QTableWidgetItem(r.kind))
            self.comparison_table.setItem(i, 2, QTableWidgetItem(str(r.n_reactions)))
            self.comparison_table.setItem(i, 3, QTableWidgetItem(r.status))
            self.comparison_table.setItem(i, 4, QTableWidgetItem(dist))
            self.comparison_table.setItem(i, 5, QTableWidgetItem(f"{r.objective:.4g}"))
        synthetic_note = (
            " <span style='color:#b45309'>⚠ synthetic demo expression.</span>"
            if synthetic_expression else ""
        )
        lethal = sum(1 for r in rows if r.status != "optimal")
        self.comparison_summary.setText(
            f"<b>Batch {method_label}</b> vs <b>{template}</b> over {len(rows)} {level} "
            f"knockout(s): {lethal} infeasible/lethal. Sorted by distance (most disrupted "
            "first)." + synthetic_note
        )
        self.status_label.setText(f"Batch comparison complete ({method_label}, {level}).")

    def _build_simulation_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)
        media_row = QHBoxLayout()
        media_row.addWidget(QLabel("Medium:"))
        self.medium_combo = QComboBox()
        self.medium_combo.addItems(sorted(PRESET_MEDIA))
        self.medium_combo.setToolTip("Preset growth medium (sets exchange uptake bounds)")
        media_row.addWidget(self.medium_combo)
        apply_medium_btn = QPushButton("Apply medium")
        apply_medium_btn.clicked.connect(self.apply_selected_medium)
        media_row.addWidget(apply_medium_btn)
        media_row.addStretch(1)
        layout.addLayout(media_row)

        controls = QHBoxLayout()
        fba_btn = QPushButton("Run FBA")
        fba_btn.clicked.connect(self.run_fba)
        pfba_btn = QPushButton("Run pFBA")
        pfba_btn.clicked.connect(self.run_pfba)
        fva_btn = QPushButton("Run FVA")
        fva_btn.clicked.connect(self.run_fva)
        controls.addWidget(fba_btn)
        controls.addWidget(pfba_btn)
        controls.addWidget(fva_btn)
        controls.addWidget(QLabel("fraction:"))
        self.fva_fraction_spin = QDoubleSpinBox()
        self.fva_fraction_spin.setRange(0.0, 1.0)
        self.fva_fraction_spin.setSingleStep(0.1)
        self.fva_fraction_spin.setValue(0.9)
        self.fva_fraction_spin.setToolTip("FVA: fraction of the optimum objective to hold")
        controls.addWidget(self.fva_fraction_spin)
        self.objective_label = QLabel("Objective: —")
        self.objective_label.setFont(QFont("", 13, QFont.Bold))
        controls.addStretch(1)
        controls.addWidget(self.objective_label)
        layout.addLayout(controls)

        self.sim_table = QTableWidget(0, 3)
        self.sim_table.setHorizontalHeaderLabels(["Reaction", "Flux", "FVA range"])
        _sim_header = self.sim_table.horizontalHeader()
        _sim_header.setSectionResizeMode(0, QHeaderView.Stretch)
        _sim_header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        _sim_header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.sim_table.verticalHeader().setVisible(False)
        self.sim_table.setAlternatingRowColors(True)
        layout.addWidget(self.sim_table, 1)
        return tab

    def _build_revert_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Normalization-target prediction (revert metabolism)")
        form = QFormLayout(controls)
        self.method_combo = QComboBox()
        self.method_combo.addItems(["rmta", "mta", "mta_miqp"])
        self.perturbation_combo = QComboBox()
        self.perturbation_combo.addItems(["gene", "reaction"])
        self.alpha_spin = QDoubleSpinBox()
        self.alpha_spin.setRange(0.0, 1.0)
        self.alpha_spin.setSingleStep(0.1)
        self.alpha_spin.setValue(0.9)
        form.addRow("Method:", self.method_combo)
        form.addRow("Knockout level:", self.perturbation_combo)
        form.addRow("Transformation weight α:", self.alpha_spin)

        source_row = QHBoxLayout()
        source_btn = QPushButton("Load source CSV/TSV…")
        source_btn.clicked.connect(lambda: self.load_revert_expression("source"))
        self.revert_source_label = QLabel("not loaded")
        source_row.addWidget(source_btn)
        source_row.addWidget(self.revert_source_label, 1)
        form.addRow("Source expression:", source_row)

        target_row = QHBoxLayout()
        target_btn = QPushButton("Load target CSV/TSV…")
        target_btn.clicked.connect(lambda: self.load_revert_expression("target"))
        self.revert_target_label = QLabel("not loaded")
        target_row.addWidget(target_btn)
        target_row.addWidget(self.revert_target_label, 1)
        form.addRow("Target expression:", target_row)

        self.revert_run_btn = QPushButton("Run Revert")
        self.revert_run_btn.setEnabled(False)
        self.revert_run_btn.clicked.connect(self.run_loaded_revert)
        self.revert_run_btn.setToolTip("Load both source and target expression files to enable.")
        form.addRow("", self.revert_run_btn)
        layout.addWidget(controls)

        self.revert_summary = QLabel(
            "Load source and target expression files to rank normalization targets."
        )
        self.revert_summary.setWordWrap(True)
        layout.addWidget(self.revert_summary)

        self.revert_table = QTableWidget(0, 3)
        self.revert_table.setHorizontalHeaderLabels(["Rank", "Target", "Score"])
        self.revert_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.revert_table.verticalHeader().setVisible(False)
        self.revert_table.setAlternatingRowColors(True)
        layout.addWidget(self.revert_table, 1)
        return tab

    # -- behavior -----------------------------------------------------------

    def load_model(self, model: Model) -> None:
        self.model = model
        solver = active_solver(model)
        self.subtitle.setText(
            f"{model.id}   ·   {len(model.reactions)} reactions   ·   "
            f"{len(model.metabolites)} metabolites   ·   {len(model.genes)} genes   ·   "
            f"solver: {solver}"
        )
        objective = ", ".join(
            r.id for r in model.reactions if r.objective_coefficient != 0
        )
        self.model_info.setText(
            f"<b>{model.id}</b><br>Objective: {objective or '—'} "
            f"({model.objective_direction})"
        )
        self._loading = True
        self.reaction_table.setRowCount(len(model.reactions))
        for row, rxn in enumerate(model.reactions):
            id_item = QTableWidgetItem(rxn.id)
            id_item.setFlags(id_item.flags() & ~Qt.ItemIsEditable)
            flux_item = QTableWidgetItem("—")
            flux_item.setFlags(flux_item.flags() & ~Qt.ItemIsEditable)
            lower_item = QTableWidgetItem(f"{rxn.lower_bound:g}")
            upper_item = QTableWidgetItem(f"{rxn.upper_bound:g}")
            self.reaction_table.setItem(row, 0, id_item)
            self.reaction_table.setItem(row, 1, lower_item)
            self.reaction_table.setItem(row, 2, upper_item)
            self.reaction_table.setItem(row, 3, flux_item)
        self._loading = False

        # Populate the production target selector with exchange reactions.
        self._loading = True
        self.product_combo.clear()
        exchanges = sorted(r.id for r in model.exchanges)
        self.product_combo.addItems(exchanges)
        self.sd_product_combo.clear()
        self.sd_product_combo.addItems(exchanges)
        self.substrate_combo.clear()
        self.substrate_combo.addItems(["auto", *exchanges])
        self._populate_ko_list()
        self._loading = False
        if self._default_product and self._default_product in exchanges:
            self.product_combo.setCurrentText(self._default_product)
            self.sd_product_combo.setCurrentText(self._default_product)

        # Production and strain design need exchange reactions; disable them otherwise.
        has_exchanges = bool(exchanges)
        for btn in getattr(self, "_production_buttons", []):
            btn.setEnabled(has_exchanges)
        if hasattr(self, "sd_run_btn"):
            self.sd_run_btn.setEnabled(has_exchanges)
        if not has_exchanges:
            self.yield_label.setText("No exchange reactions in this model — production design unavailable.")

        self._fluxes = {}
        self._fluxes_stale = False
        self.status_label.setText(f"Loaded model '{model.id}'.")

    def set_reaction_bounds(
        self, reaction_id: str, lower: float | None = None, upper: float | None = None
    ) -> None:
        """Programmatically update a reaction bound and reflect it in the table."""

        rxn = self.model.reactions.get_by_id(reaction_id)
        # Assign both bounds at once so a value that crosses the current opposite bound does
        # not trip cobra's per-assignment validation (lower <= upper).
        new_lower = float(lower) if lower is not None else rxn.lower_bound
        new_upper = float(upper) if upper is not None else rxn.upper_bound
        if new_lower > new_upper:
            new_lower, new_upper = new_upper, new_lower
        rxn.bounds = (new_lower, new_upper)
        self._mark_fluxes_stale()
        for row in range(self.reaction_table.rowCount()):
            if self.reaction_table.item(row, 0).text() == reaction_id:
                self._loading = True
                self.reaction_table.item(row, 1).setText(f"{rxn.lower_bound:g}")
                self.reaction_table.item(row, 2).setText(f"{rxn.upper_bound:g}")
                self._loading = False
                break

    def _on_bound_edited(self, item) -> None:
        """Apply an in-table edit of a Lower/Upper cell to the model."""

        if self._loading or item.column() not in (1, 2):
            return
        row = item.row()
        reaction_id = self.reaction_table.item(row, 0).text()
        rxn = self.model.reactions.get_by_id(reaction_id)
        try:
            value = float(item.text())
        except ValueError:
            self._loading = True
            current = rxn.lower_bound if item.column() == 1 else rxn.upper_bound
            item.setText(f"{current:g}")
            self._loading = False
            self.status_label.setText(f"Invalid bound for {reaction_id}; reverted.")
            return
        if item.column() == 1:
            clamped = min(value, rxn.upper_bound)
            rxn.bounds = (clamped, rxn.upper_bound)
        else:
            clamped = max(value, rxn.lower_bound)
            rxn.bounds = (rxn.lower_bound, clamped)
        # Reflect the actually-stored (possibly clamped) value back into the cell.
        self._loading = True
        item.setText(f"{clamped:g}")
        self._loading = False
        note = "" if clamped == value else "  (clamped to keep lower <= upper)"
        self.status_label.setText(
            f"{reaction_id} bounds set to ({rxn.lower_bound:g}, {rxn.upper_bound:g}).{note}"
        )
        self._mark_fluxes_stale()

    def _mark_fluxes_stale(self) -> None:
        """Flag that bounds changed so the displayed flux column is out of date."""

        self._fluxes_stale = True
        if self._fluxes:
            self.status_label.setText("Bounds changed — fluxes are stale; click Run FBA.")

    def run_fba(self) -> None:
        try:
            solution = fba(self.model)
        except Exception as exc:  # surface solver/model errors instead of freezing
            self.objective_label.setText("Objective: error")
            self.status_label.setText(f"FBA failed: {exc}")
            return
        self._fluxes = dict(solution.fluxes)
        self._fluxes_stale = False
        obj = solution.objective_value
        obj_text = f"{obj:.4g}" if obj is not None else "infeasible"
        self.objective_label.setText(f"Objective: {obj_text} ({solution.status})")
        if obj is None or solution.status != "optimal":
            self.status_label.setText(f"FBA status: {solution.status} — model may be infeasible.")
        else:
            self.status_label.setText("FBA complete.")
        self._populate_flux_column()
        self._fill_sim_table()

    def run_pfba(self) -> None:
        try:
            solution = pfba(self.model)
        except Exception as exc:
            self.objective_label.setText("Objective: error")
            self.status_label.setText(f"pFBA failed: {exc}")
            return
        self._fluxes = dict(solution.fluxes)
        self._fluxes_stale = False
        # pFBA's objective value is the minimal total flux; show the growth (objective rxn).
        growth = next(
            (self._fluxes[r.id] for r in self.model.reactions if r.objective_coefficient != 0),
            None,
        )
        growth_text = f"{growth:.4g}" if growth is not None else "—"
        total = sum(abs(v) for v in self._fluxes.values())
        self.objective_label.setText(f"Objective: {growth_text} (pFBA)")
        self.status_label.setText(f"pFBA complete (minimal total flux = {total:.1f}).")
        self._populate_flux_column()
        self._fill_sim_table()

    def apply_selected_medium(self) -> None:
        name = self.medium_combo.currentText()
        try:
            applied = apply_medium(self.model, name)
        except Exception as exc:
            self.status_label.setText(f"Could not apply medium: {exc}")
            return
        # Reflect the new exchange bounds in the reaction table and mark fluxes stale.
        self.load_model(self.model)
        self.status_label.setText(
            f"Applied medium '{name}' ({len(applied)} open exchanges). Re-run FBA."
        )

    def run_fva(self) -> None:
        if not self._fluxes or self._fluxes_stale:
            self.run_fba()
        fraction = self.fva_fraction_spin.value()
        try:
            ranges = self._run_in_background(
                lambda: fva(self.model, fraction_of_optimum=fraction),
                label="Running FVA…",
            )
        except Exception as exc:
            self.status_label.setText(f"FVA failed: {exc}")
            return
        self.sim_table.setHorizontalHeaderLabels(["Reaction", "Flux", f"FVA range (f={fraction:g})"])
        self.sim_table.setRowCount(len(self._fluxes))
        for row, (rid, flux) in enumerate(sorted(self._fluxes.items())):
            self.sim_table.setItem(row, 0, QTableWidgetItem(rid))
            self.sim_table.setItem(row, 1, QTableWidgetItem(f"{flux:.3g}"))
            rng = ranges.get(rid)
            text = f"[{rng.minimum:.3g}, {rng.maximum:.3g}]" if rng else "—"
            self.sim_table.setItem(row, 2, QTableWidgetItem(text))
        self.status_label.setText(f"FVA complete (fraction {fraction:g}).")

    def _fill_sim_table(self) -> None:
        self.sim_table.setRowCount(len(self._fluxes))
        for row, (rid, flux) in enumerate(sorted(self._fluxes.items())):
            self.sim_table.setItem(row, 0, QTableWidgetItem(rid))
            self.sim_table.setItem(row, 1, QTableWidgetItem(f"{flux:.3g}"))
            self.sim_table.setItem(row, 2, QTableWidgetItem("—"))

    def _populate_flux_column(self) -> None:
        for row in range(self.reaction_table.rowCount()):
            rid = self.reaction_table.item(row, 0).text()
            flux = self._fluxes.get(rid, 0.0)
            self.reaction_table.setItem(row, 3, QTableWidgetItem(f"{flux:.3g}"))
        self._apply_flux_threshold(self.flux_slider.value())

    def _apply_flux_threshold(self, value: int) -> None:
        self.flux_threshold_label.setText(f"{value}%")
        if not self._fluxes:
            return
        max_abs = max((abs(v) for v in self._fluxes.values()), default=0.0) or 1.0
        threshold = (value / 100.0) * max_abs
        # Recolouring emits itemChanged; guard so it is not mistaken for a bound edit.
        self._loading = True
        for row in range(self.reaction_table.rowCount()):
            rid = self.reaction_table.item(row, 0).text()
            active = abs(self._fluxes.get(rid, 0.0)) >= threshold and threshold >= 0
            color = QColor("#dce9d6") if active and value > 0 else QColor("white")
            for col in range(4):
                item = self.reaction_table.item(row, col)
                if item is not None:
                    item.setBackground(color)
        self._loading = False

    # -- production design --------------------------------------------------

    def _current_product(self) -> str:
        return self.product_combo.currentText()

    def _current_substrate(self) -> str | None:
        text = self.substrate_combo.currentText()
        return None if text in ("", "auto") else text

    def _is_aerobic(self) -> bool:
        return self.anaerobic_combo.currentText() == "aerobic"

    def render_flux_map(self) -> None:
        """Render the Escher-layout flux map coloured by the current FBA flux."""

        if self._map_path is None:
            return
        if not self._fluxes or self._fluxes_stale:
            self.run_fba()
        try:
            fig = escher_flux_map(self._map_path, self._fluxes,
                                  title=f"{self.model.id} — flux map")
        except Exception as exc:
            self.status_label.setText(f"Flux map failed: {exc}")
            return
        if self._map_canvas is not None:
            old = self._map_canvas.figure
            self.map_canvas_holder.removeWidget(self._map_canvas)
            self._map_canvas.setParent(None)
            if self._map_toolbar is not None:
                self.map_canvas_holder.removeWidget(self._map_toolbar)
                self._map_toolbar.setParent(None)
            old.clf()
        fig.set_dpi(110)
        self._map_canvas = FigureCanvas(fig)
        self._map_toolbar = NavigationToolbar(self._map_canvas, self)
        self.map_canvas_holder.addWidget(self._map_toolbar)
        self.map_canvas_holder.addWidget(self._map_canvas)
        self._map_canvas.draw()
        self.status_label.setText("Flux map rendered (Escher layout).")

    def _on_product_changed(self, product: str) -> None:
        if not self._loading and product:
            self.yield_label.setText(
                f"Compute the theoretical yield, envelope, or FSEOF targets for "
                f"{html.escape(product)}."
            )

    def _set_production_figure(self, fig) -> None:
        if self._production_canvas is not None:
            old_fig = self._production_canvas.figure
            self.production_canvas_holder.removeWidget(self._production_canvas)
            self._production_canvas.setParent(None)
            if self._production_toolbar is not None:
                self.production_canvas_holder.removeWidget(self._production_toolbar)
                self._production_toolbar.setParent(None)
            old_fig.clf()  # release the previous figure so it is not leaked
        # Figures are authored at 300 DPI for file export; render on screen at a lower DPI
        # so the canvas fits the panel instead of clipping.
        fig.set_dpi(100)
        self._production_canvas = FigureCanvas(fig)
        # A matplotlib toolbar gives interactive zoom/pan and a "Save the figure" button.
        self._production_toolbar = NavigationToolbar(self._production_canvas, self)
        self.production_canvas_holder.addWidget(self._production_toolbar)
        self.production_canvas_holder.addWidget(self._production_canvas)
        self._production_canvas.draw()

    def _run_production(self, action) -> None:
        """Run a production analysis, surfacing any error instead of crashing the UI."""

        try:
            action()
        except Exception as exc:
            self.yield_label.setText(f"Production analysis failed: {html.escape(str(exc))}")
            self.status_label.setText(f"Production analysis failed: {exc}")

    def run_theoretical_yield(self) -> None:
        def _do():
            product = self._current_product()
            result = theoretical_yield(self.model, product,
                                       substrate=self._current_substrate(),
                                       aerobic=self._is_aerobic())
            ceiling = (
                f", carbon ceiling {result.carbon_ceiling:.2f}"
                if result.carbon_ceiling is not None else ""
            )
            co2 = (
                f"; needs net CO₂ fixation ({result.co2_exchange:.1f})"
                if result.co2_fixed and result.exceeds_carbon_ceiling else ""
            )
            self.yield_label.setText(
                f"Theoretical yield of {html.escape(product)}: "
                f"<b>{result.molar_yield:.3f}</b> mol/mol {html.escape(result.substrate)} "
                f"({'aerobic' if result.aerobic else 'anaerobic'}; "
                f"max flux {result.product_flux:.2f}, uptake {result.substrate_uptake:.1f}{ceiling}{co2})."
            )
            self.status_label.setText("Theoretical yield computed.")

        self._run_production(_do)

    def run_production_envelope_plot(self) -> None:
        def _do():
            product = self._current_product()
            aerobic = self._is_aerobic()
            condition = "aerobic" if aerobic else "anaerobic"
            substrate = self._current_substrate()
            envelope = self._run_in_background(
                lambda: production_envelope(self.model, product, substrate=substrate,
                                            aerobic=aerobic, points=20),
                label="Computing production envelope…",
            )
            fig = production_envelope_figure(
                envelope, title=f"Production envelope — {product} ({condition})"
            )
            self._set_production_figure(fig)
            self.status_label.setText(f"Production envelope computed ({condition}).")

        self._run_production(_do)

    def run_fseof_plot(self) -> None:
        def _do():
            product = self._current_product()
            aerobic = self._is_aerobic()
            result = self._run_in_background(
                lambda: fseof(self.model, product, n_steps=10, aerobic=aerobic),
                label="Running FSEOF…",
            )
            if result.metadata.get("max_product", 0.0) <= 1e-9:
                self.yield_label.setText(
                    f"FSEOF not meaningful: theoretical yield of {html.escape(product)} is zero "
                    "— check the exchange is open and the product is reachable."
                )
                self.status_label.setText("FSEOF: zero theoretical yield.")
                return
            fig = fseof_figure(result, top_n=6, title=f"FSEOF targets — {product}")
            self._set_production_figure(fig)
            # Rank the listed targets the same way the figure does (by flux increase).
            levels = list(result.enforced_levels)
            ranked = sorted(
                result.amplification_targets(),
                key=lambda rid: abs(result.trends.loc[rid, levels[-1]])
                - abs(result.trends.loc[rid, levels[0]]),
                reverse=True,
            )[:6]
            self.yield_label.setText(
                "FSEOF amplification targets: " + html.escape(", ".join(ranked))
            )
            self.status_label.setText("FSEOF complete.")

        self._run_production(_do)

    def run_fvseof_plot(self) -> None:
        def _do():
            product = self._current_product()
            aerobic = self._is_aerobic()
            result = self._run_in_background(
                lambda: fvseof(self.model, product, n_steps=8, aerobic=aerobic),
                label="Running FVSEOF…",
            )
            if result.metadata.get("max_product", 0.0) <= 1e-9:
                self.yield_label.setText(
                    f"FVSEOF not meaningful: theoretical yield of {html.escape(product)} is zero."
                )
                self.status_label.setText("FVSEOF: zero theoretical yield.")
                return
            fig = fvseof_figure(result, top_n=5, title=f"FVSEOF robust targets — {product}")
            self._set_production_figure(fig)
            robust = result.robust_targets()[:8]
            self.yield_label.setText(
                "FVSEOF robust amplification targets (forced up): "
                + html.escape(", ".join(robust) if robust else "none")
            )
            self.status_label.setText("FVSEOF complete.")

        self._run_production(_do)

    # -- omics integration --------------------------------------------------

    def run_omics(self, gene_expression: Mapping[str, float]) -> None:
        """Predict a flux distribution from gene expression (E-Flux2 or LAD)."""

        method = self.omics_method_combo.currentText()
        try:
            weights = gene_to_reaction_weights(self.model, dict(gene_expression))
            result = integrate_expression(self.model, dict(gene_expression), method=method)
        except Exception as exc:
            self.omics_summary.setText(f"Omics integration failed: {html.escape(str(exc))}")
            return
        if result.status != "optimal":
            self.omics_summary.setText(
                f"{method.upper()} returned status '{result.status}' "
                f"({len(weights)} reactions mapped from expression)."
            )
            self.omics_table.setRowCount(0)
            return
        ranked = sorted(result.fluxes.items(), key=lambda kv: -abs(kv[1]))
        active = [(rid, v) for rid, v in ranked if abs(v) > 1e-6]
        self.omics_table.setRowCount(len(active))
        for i, (rid, flux) in enumerate(active):
            self.omics_table.setItem(i, 0, QTableWidgetItem(rid))
            self.omics_table.setItem(i, 1, QTableWidgetItem(f"{flux:.4g}"))
        detail = f" [{result.detail}]" if result.detail else ""
        self.omics_summary.setText(
            f"<b>{method.upper()}{detail}</b>: {len(weights)} reactions mapped from "
            f"{len(gene_expression)} genes; {len(active)} active fluxes. "
            f"Objective (FBA optimum) = {result.objective_value:.4g}."
        )
        self.status_label.setText(f"Omics integration complete ({method}).")

    def run_omics_demo(self) -> None:
        """Run omics integration on deterministic demo expression over the model's genes."""

        import numpy as np

        rng = np.random.default_rng(0)
        expression = {g.id: float(rng.uniform(1.0, 100.0)) for g in self.model.genes}
        if not expression:
            self.omics_summary.setText("Model has no genes — omics integration is not applicable.")
            return
        self.run_omics(expression)

    def load_expression_csv(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load gene-expression CSV", "", "CSV files (*.csv *.tsv);;All files (*)"
        )
        if not path:
            return
        try:
            expression = _read_expression_vector(path)
        except Exception as exc:
            self.omics_summary.setText(f"Could not read expression CSV: {html.escape(str(exc))}")
            return
        # Remember real (CSV-loaded) expression so the Comparison tab's LAD/E-Flux2 template
        # uses it instead of synthetic data. Demo expression is intentionally not stored.
        self._omics_expression = dict(expression)
        self.run_omics(expression)

    def _update_revert_run_state(self) -> None:
        ready = (
            self._revert_source_expression is not None
            and self._revert_target_expression is not None
        )
        self.revert_run_btn.setEnabled(ready)

    def load_revert_expression(self, role: str) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, f"Load {role} expression", "", "CSV/TSV files (*.csv *.tsv);;All files (*)"
        )
        if not path:
            return
        try:
            expression = _read_expression_vector(path)
        except Exception as exc:
            self.revert_summary.setText(
                f"Could not read {html.escape(role)} expression: {html.escape(str(exc))}"
            )
            return
        if role == "source":
            self._revert_source_expression = expression
            self.revert_source_label.setText(f"{len(expression)} genes")
        else:
            self._revert_target_expression = expression
            self.revert_target_label.setText(f"{len(expression)} genes")
        self._update_revert_run_state()
        self.status_label.setText(f"Loaded {role} expression ({len(expression)} genes).")

    def run_loaded_revert(self) -> None:
        if self._revert_source_expression is None or self._revert_target_expression is None:
            self.revert_summary.setText("Load both source and target expression files first.")
            self._update_revert_run_state()
            return
        self.run_revert(self._revert_source_expression, self._revert_target_expression)

    def run_revert(
        self,
        source_expression: Mapping[str, float],
        target_expression: Mapping[str, float],
    ) -> None:
        method = self.method_combo.currentText()
        perturbation = self.perturbation_combo.currentText()
        alpha = self.alpha_spin.value()
        def _compute():
            reference = reference_state_pfba(self.model, name="source")
            direction = differential_expression(
                self.model, source_expression, target_expression, reference=reference
            )
            return revert_targets(
                self.model,
                None,
                reference,
                direction,
                method=method,
                alpha=alpha,
                perturbation=perturbation,
            )

        try:
            ranking = self._run_in_background(_compute, label=f"Running revert ({method})…")
        except Exception as exc:  # surface solver-capability/model errors instead of freezing
            self.revert_summary.setText(f"Revert prediction failed: {exc}")
            self.status_label.setText(f"Revert prediction failed ({method}).")
            return
        rows = ranking.to_records()
        self.revert_table.setRowCount(len(rows))
        for i, record in enumerate(rows):
            self.revert_table.setItem(i, 0, QTableWidgetItem(str(record["rank"])))
            self.revert_table.setItem(i, 1, QTableWidgetItem(record["target_id"]))
            self.revert_table.setItem(i, 2, QTableWidgetItem(f"{record['score']:.4g}"))
            if i == 0:
                for col in range(3):
                    self.revert_table.item(i, col).setBackground(QColor("#dce9d6"))
        best = ranking.best()
        if best is not None:
            self.revert_summary.setText(
                f"<b>Top normalization target: {best.target_id}</b> "
                f"(score {best.score:.4g}, method {method}, {perturbation} knockout). "
                f"Knocking it out moves the source state most toward the target state."
            )
        self.status_label.setText(f"Revert prediction complete ({method}).")

    # -- transformation finder (A -> B) -------------------------------------

    def _build_transformation_tab(self) -> QWidget:
        """Rank knockouts that move a source flux state toward a target flux state (A→B)."""

        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Transformation targets (move source → target flux state)")
        form = QFormLayout(controls)
        self.transform_method_combo = QComboBox()
        self.transform_method_combo.addItems(["moma", "mta"])
        self.transform_method_combo.setToolTip(
            "moma: rank by how much the minimal-adjustment state moves toward the target. "
            "mta: robust MTA on the source→target direction."
        )
        self.transform_omics_combo = QComboBox()
        self.transform_omics_combo.addItems(["eflux2", "lad"])
        self.transform_omics_combo.setToolTip("How each expression vector is turned into a flux state")
        self.transform_perturbation_combo = QComboBox()
        self.transform_perturbation_combo.addItems(["gene", "reaction"])
        self.transform_alpha_spin = QDoubleSpinBox()
        self.transform_alpha_spin.setRange(0.0, 1.0)
        self.transform_alpha_spin.setSingleStep(0.1)
        self.transform_alpha_spin.setValue(0.9)
        form.addRow("Method:", self.transform_method_combo)
        form.addRow("Predict states with:", self.transform_omics_combo)
        form.addRow("Knockout level:", self.transform_perturbation_combo)
        form.addRow("Transformation weight α (mta):", self.transform_alpha_spin)

        source_row = QHBoxLayout()
        source_btn = QPushButton("Load source (A) CSV/TSV…")
        source_btn.clicked.connect(lambda: self.load_transform_expression("source"))
        self.transform_source_label = QLabel("not loaded")
        source_row.addWidget(source_btn)
        source_row.addWidget(self.transform_source_label, 1)
        form.addRow("Source expression:", source_row)

        target_row = QHBoxLayout()
        target_btn = QPushButton("Load target (B) CSV/TSV…")
        target_btn.clicked.connect(lambda: self.load_transform_expression("target"))
        self.transform_target_label = QLabel("not loaded")
        target_row.addWidget(target_btn)
        target_row.addWidget(self.transform_target_label, 1)
        form.addRow("Target expression:", target_row)

        self.transform_run_btn = QPushButton("Run transformation")
        self.transform_run_btn.setEnabled(False)
        self.transform_run_btn.clicked.connect(self.run_transformation)
        self.transform_run_btn.setToolTip("Load both source and target expression files to enable.")
        form.addRow("", self.transform_run_btn)
        layout.addWidget(controls)

        self.transform_summary = QLabel(
            "Load source (A) and target (B) expression files; each is predicted to a flux "
            "state, then knockouts are ranked by how well they move A toward B."
        )
        self.transform_summary.setWordWrap(True)
        layout.addWidget(self.transform_summary)

        self.transform_table = QTableWidget(0, 3)
        self.transform_table.setHorizontalHeaderLabels(["Rank", "Target", "Score"])
        self.transform_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.transform_table.verticalHeader().setVisible(False)
        self.transform_table.setAlternatingRowColors(True)
        layout.addWidget(self.transform_table, 1)
        return tab

    def _update_transform_run_state(self) -> None:
        ready = (
            self._transform_source_expression is not None
            and self._transform_target_expression is not None
        )
        self.transform_run_btn.setEnabled(ready)

    def load_transform_expression(self, role: str) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, f"Load {role} expression", "", "CSV/TSV files (*.csv *.tsv);;All files (*)"
        )
        if not path:
            return
        try:
            expression = _read_expression_vector(path)
        except Exception as exc:
            self.transform_summary.setText(
                f"Could not read {html.escape(role)} expression: {html.escape(str(exc))}"
            )
            return
        if role == "source":
            self._transform_source_expression = expression
            self.transform_source_label.setText(f"{len(expression)} genes")
        else:
            self._transform_target_expression = expression
            self.transform_target_label.setText(f"{len(expression)} genes")
        self._update_transform_run_state()
        self.status_label.setText(f"Loaded {role} (A→B) expression ({len(expression)} genes).")

    def run_transformation(self) -> None:
        if (
            self._transform_source_expression is None
            or self._transform_target_expression is None
        ):
            self.transform_summary.setText("Load both source (A) and target (B) expression files first.")
            return
        method = self.transform_method_combo.currentText()
        omics_method = self.transform_omics_combo.currentText()
        perturbation = self.transform_perturbation_combo.currentText()
        alpha = self.transform_alpha_spin.value()
        source_expr = self._transform_source_expression
        target_expr = self._transform_target_expression

        def _compute():
            source_state = integrate_expression(
                self.model, source_expr, method=omics_method
            ).to_flux_state("source")
            target_state = integrate_expression(
                self.model, target_expr, method=omics_method
            ).to_flux_state("target")
            return transformation_targets(
                self.model, source_state, target_state,
                method=method, perturbation=perturbation, alpha=alpha,
            )

        try:
            ranking = self._run_in_background(
                _compute, label=f"Running transformation ({method})…"
            )
        except Exception as exc:  # surface solver-capability/model errors instead of freezing
            self.transform_summary.setText(f"Transformation failed: {html.escape(str(exc))}")
            self.status_label.setText(f"Transformation failed ({method}).")
            return
        rows = ranking.to_records()
        self.transform_table.setRowCount(len(rows))
        for i, record in enumerate(rows):
            self.transform_table.setItem(i, 0, QTableWidgetItem(str(record["rank"])))
            self.transform_table.setItem(i, 1, QTableWidgetItem(record["target_id"]))
            self.transform_table.setItem(i, 2, QTableWidgetItem(f"{record['score']:.4g}"))
            if i == 0:
                for col in range(3):
                    self.transform_table.item(i, col).setBackground(QColor("#dce9d6"))
        best = ranking.best()
        if best is not None:
            self.transform_summary.setText(
                f"<b>Top A→B target: {html.escape(best.target_id)}</b> "
                f"(score {best.score:.4g}, method {method}, {omics_method} states, "
                f"{perturbation} knockout). Higher score = moves the source state more toward "
                "the target."
            )
        else:
            self.transform_summary.setText("No candidate knockouts scored.")
        self.status_label.setText(f"Transformation complete ({method}).")

    # -- multi-condition comparison -----------------------------------------

    def _build_conditions_tab(self) -> QWidget:
        """Predict per-condition fluxes from a gene×condition table and compare two as log2FC."""

        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Multi-condition flux comparison (expression table → per-condition flux)")
        row = QHBoxLayout(controls)
        row.addWidget(QLabel("Method:"))
        self.cond_method_combo = QComboBox()
        self.cond_method_combo.addItems(["eflux2", "lad"])
        row.addWidget(self.cond_method_combo)
        load_btn = QPushButton("Load expression table CSV…")
        load_btn.clicked.connect(self.load_condition_table)
        row.addWidget(load_btn)
        row.addWidget(QLabel("A:"))
        self.cond_source_combo = QComboBox()
        self.cond_source_combo.setMinimumWidth(120)
        row.addWidget(self.cond_source_combo)
        row.addWidget(QLabel("B:"))
        self.cond_target_combo = QComboBox()
        self.cond_target_combo.setMinimumWidth(120)
        row.addWidget(self.cond_target_combo)
        self.cond_run_btn = QPushButton("Compare")
        self.cond_run_btn.setEnabled(False)
        self.cond_run_btn.clicked.connect(self.run_condition_comparison)
        row.addWidget(self.cond_run_btn)
        row.addStretch(1)
        layout.addWidget(controls)

        self.cond_summary = QLabel(
            "Load a gene × condition table (first column = gene id, remaining columns = "
            "conditions), pick two conditions A and B, and compare their predicted fluxes."
        )
        self.cond_summary.setWordWrap(True)
        layout.addWidget(self.cond_summary)

        self.cond_table = QTableWidget(0, 2)
        self.cond_table.setHorizontalHeaderLabels(["Reaction", "log2( |flux B| / |flux A| )"])
        self.cond_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.cond_table.verticalHeader().setVisible(False)
        self.cond_table.setAlternatingRowColors(True)
        layout.addWidget(self.cond_table, 1)
        return tab

    def load_condition_table(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load expression table", "", "CSV/TSV files (*.csv *.tsv);;All files (*)"
        )
        if not path:
            return
        try:
            table = read_expression_table(path)
        except Exception as exc:
            self.cond_summary.setText(f"Could not read expression table: {html.escape(str(exc))}")
            return
        if table.shape[1] < 2:
            self.cond_summary.setText(
                "Expression table needs at least two condition columns to compare."
            )
            return
        self._condition_table = table
        conditions = [str(c) for c in table.columns]
        self._loading = True
        self.cond_source_combo.clear()
        self.cond_source_combo.addItems(conditions)
        self.cond_target_combo.clear()
        self.cond_target_combo.addItems(conditions)
        self.cond_target_combo.setCurrentIndex(min(1, len(conditions) - 1))
        self._loading = False
        self.cond_run_btn.setEnabled(True)
        self.cond_summary.setText(
            f"Loaded {table.shape[0]} genes × {table.shape[1]} conditions "
            f"({', '.join(conditions)}). Pick A and B, then Compare."
        )
        self.status_label.setText(f"Loaded expression table ({table.shape[1]} conditions).")

    def run_condition_comparison(self) -> None:
        if self._condition_table is None:
            self.cond_summary.setText("Load an expression table first.")
            return
        source = self.cond_source_combo.currentText()
        target = self.cond_target_combo.currentText()
        if source == target:
            self.cond_summary.setText("Pick two different conditions for A and B.")
            return
        method = self.cond_method_combo.currentText()
        table = self._condition_table

        def _compute():
            predictions = predict_condition_fluxes(
                self.model, table, method=method, conditions=[source, target]
            )
            return flux_log_change(predictions.fluxes(source), predictions.fluxes(target))

        try:
            log_change = self._run_in_background(
                _compute, label=f"Comparing conditions ({method})…"
            )
        except Exception as exc:
            self.cond_summary.setText(f"Condition comparison failed: {html.escape(str(exc))}")
            self.status_label.setText("Condition comparison failed.")
            return
        ranked = sorted(
            ((rid, v) for rid, v in log_change.items() if v == v),  # drop NaN
            key=lambda kv: -abs(kv[1]),
        )
        changed = [(rid, v) for rid, v in ranked if abs(v) > 1e-6]
        self.cond_table.setRowCount(len(changed))
        for i, (rid, value) in enumerate(changed):
            self.cond_table.setItem(i, 0, QTableWidgetItem(rid))
            text = "+inf" if value == float("inf") else "-inf" if value == float("-inf") else f"{value:.4g}"
            self.cond_table.setItem(i, 1, QTableWidgetItem(text))
        self.cond_summary.setText(
            f"<b>{method}</b>: {source} → {target}. {len(changed)} reactions changed flux "
            f"magnitude (log2 fold-change; positive = higher in {html.escape(target)})."
        )
        self.status_label.setText(f"Condition comparison complete ({method}).")
