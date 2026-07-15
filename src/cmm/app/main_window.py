"""CMM desktop main window.

A thin Qt view over the solver-neutral services. All numerics live in cmm.core / cmm.features
/ cmm.omics; this module only renders results and wires controls, so the GUI inherits none of
the analysis logic.
"""

from __future__ import annotations

import html
import math
from collections.abc import Mapping
from pathlib import Path

_ASSETS_DIR = Path(__file__).resolve().parent / "assets"

from cobra import Model
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from qtpy.QtCore import QEvent, QEventLoop, QObject, Qt, QThread, QTimer, Signal
from qtpy.QtGui import QColor, QFont
from qtpy.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QProgressDialog,
    QPushButton,
    QSlider,
    QSpinBox,
    QSplitter,
    QStyledItemDelegate,
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
    predict_condition_fluxes,
    read_expression_table,
)
from cmm.omics.differential import differential_expression
from cmm.omics.expression import integrate_expression
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
QComboBox::down-arrow { image: url(__ASSETS__/down_arrow.svg);
    width: 10px; height: 7px; margin-right: 8px; }
QComboBox QAbstractItemView { background: #ffffff; color: #1a2433; border: 1px solid #cdd6e1;
    selection-background-color: #2f5e8f; selection-color: #ffffff; outline: none; }

/* Spin boxes */
QDoubleSpinBox, QSpinBox { background: #ffffff; color: #1a2433; border: 1px solid #cdd6e1;
    border-radius: 6px; padding: 3px 6px; min-height: 22px; }
QDoubleSpinBox:focus, QSpinBox:focus { border-color: #2f5e8f; }
/* Styling a spin box in QSS drops the native step arrows, so draw them explicitly. */
QDoubleSpinBox::up-button, QSpinBox::up-button { subcontrol-origin: border;
    subcontrol-position: top right; width: 16px; border-left: 1px solid #cdd6e1; }
QDoubleSpinBox::down-button, QSpinBox::down-button { subcontrol-origin: border;
    subcontrol-position: bottom right; width: 16px; border-left: 1px solid #cdd6e1; }
QDoubleSpinBox::up-button:hover, QSpinBox::up-button:hover,
QDoubleSpinBox::down-button:hover, QSpinBox::down-button:hover { background: #eef2f7; }
QDoubleSpinBox::up-arrow, QSpinBox::up-arrow {
    image: url(__ASSETS__/up_arrow.svg); width: 9px; height: 6px; }
QDoubleSpinBox::down-arrow, QSpinBox::down-arrow {
    image: url(__ASSETS__/down_arrow.svg); width: 9px; height: 6px; }

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
/* No font-weight here: Qt's tabSizeHint ignores a stylesheet font-weight, so a bold label is
   drawn wider than the tab Qt reserved for it and clips. The tab-bar font (incl. its DemiBold
   weight) is set on the widget in code so the size hint measures the real text; the selected
   tab is distinguished by colour + underline. */
QTabBar::tab { background: #e6ebf2; color: #54637a; padding: 7px 13px; margin-right: 2px;
    border: 1px solid #d8dfe8; border-bottom: none;
    border-top-left-radius: 6px; border-top-right-radius: 6px; }
QTabBar::tab:hover { background: #eef2f7; color: #2f5e8f; }
QTabBar::tab:selected { background: #ffffff; color: #1b2b44;
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
        raise ValueError(
            "expression file must have at least two columns: gene, expression"
        )
    gene_col, value_col = frame.columns[:2]
    expression: dict[str, float] = {}
    for gene, value in zip(frame[gene_col], frame[value_col], strict=False):
        if pd.isna(gene) or pd.isna(value):
            continue
        gene_id = str(gene)
        if gene_id in expression:
            raise ValueError(f"expression file contains duplicate gene {gene_id!r}")
        numeric = float(value)
        if not math.isfinite(numeric) or numeric < 0:
            raise ValueError(
                f"expression for gene {gene_id!r} must be finite and non-negative"
            )
        expression[gene_id] = numeric
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
        except (
            Exception
        ) as exc:  # captured and re-raised on the UI thread by the caller
            self.error = exc
        finally:
            self.finished.emit()


class _ComboItemDelegate(QStyledItemDelegate):
    """Left-align combo popup entries and paint disabled ones as clearly inert.

    QSS ``::item:disabled`` styling is unreliable for combo popups (notably on macOS, where it
    does not render), so disabled rows are painted here: a muted neutral fill with faint text
    and no hover/selection highlight. Enabled rows defer to normal styled rendering.
    """

    _DISABLED_TEXT = QColor("#b3bcc7")
    _DISABLED_FILL = QColor("#edeef1")
    _ALIGN = int(Qt.AlignLeft | Qt.AlignVCenter)

    def initStyleOption(self, option, index):  # noqa: N802 (Qt override)
        super().initStyleOption(option, index)
        option.displayAlignment = Qt.AlignLeft | Qt.AlignVCenter

    def paint(self, painter, option, index):
        if not (index.flags() & Qt.ItemIsEnabled):
            painter.save()
            painter.fillRect(option.rect, self._DISABLED_FILL)
            painter.setPen(self._DISABLED_TEXT)
            text = index.data(Qt.DisplayRole)
            # Indent to match the styled left padding of enabled rows.
            rect = option.rect.adjusted(6, 0, -6, 0)
            painter.drawText(rect, self._ALIGN, "" if text is None else str(text))
            painter.restore()
            return
        super().paint(painter, option, index)


class _DisabledHoverCursor(QObject):
    """Show a not-allowed cursor over disabled entries in a combo-box popup."""

    def __init__(self, view):
        super().__init__(view)
        self._view = view

    def eventFilter(self, viewport, event):  # noqa: N802 (Qt override)
        if event.type() == QEvent.MouseMove:
            try:
                pos = event.position().toPoint()  # Qt6
            except AttributeError:
                pos = event.pos()  # Qt5
            index = self._view.indexAt(pos)
            enabled = (
                bool(self._view.model().flags(index) & Qt.ItemIsEnabled)
                if index.isValid() else True
            )
            viewport.setCursor(Qt.ArrowCursor if enabled else Qt.ForbiddenCursor)
        return False


class CmmMainWindow(QMainWindow):
    """Main platform window over a single cobra model."""

    def __init__(
        self,
        model: Model,
        parent=None,
        *,
        default_product: str | None = None,
        map_path: str | None = None,
    ):
        super().__init__(parent)
        self.model = model
        self._fluxes: dict[str, float] = {}  # the last-run distribution (drives left panel/FVA/map)
        self._fba_fluxes: dict[str, float] = {}  # FBA column of the simulation table
        self._pfba_fluxes: dict[str, float] = {}  # pFBA column of the simulation table
        self._loading = False
        self._fluxes_stale = False
        self._default_product = default_product
        self._map_path = map_path
        self._revert_source_expression: dict[str, float] | None = None
        self._revert_target_expression: dict[str, float] | None = None
        # Real gene expression loaded from a CSV on the Omics tab, reused as the LAD/E-Flux2
        # reference template on the Comparison tab so that path is not left to synthetic data.
        self._omics_expression: dict[str, float] | None = None
        # Which integration methods ("lad"/"eflux2") were actually computed on the Omics tab —
        # only those are offered as a Comparison reference template.
        self._omics_computed_methods: set[str] = set()
        # Last single-run MOMA/ROOM result, kept so the change-threshold can re-filter the
        # table without re-solving.
        self._comparison_cache: dict | None = None
        # Last production result (kind, result) so the FSEOF/FVSEOF "show all" toggle can
        # re-render the table without recomputing.
        self._production_result: tuple[str, object] | None = None
        # Keep hover-cursor event filters alive for combos that grey out entries.
        self._cursor_filters: list[_DisabledHoverCursor] = []
        # A→B transformation finder: two expression vectors → two predicted flux states.
        self._transform_source_expression: dict[str, float] | None = None
        self._transform_target_expression: dict[str, float] | None = None
        # Omics integration source: a gene × condition expression table (1+ conditions).
        # Loading only stores it; Compute predicts a flux column per selected condition.
        self._omics_table_df = None
        # Last computed per-condition fluxes (condition -> {reaction: flux}), kept so the
        # active-only / show-all toggle can re-render without recomputing.
        self._omics_fluxes_by_condition: dict[str, dict[str, float]] | None = None
        self._omics_conditions_order: list[str] = []
        # Guards re-entrant background runs (see _run_in_background).
        self._busy = False
        self.setWindowTitle("CMM — Cellular Metabolic Modeling Platform")
        self.resize(1160, 760)
        # Qt QSS ``image: url(...)`` needs a real path; forward slashes work cross-platform.
        self.setStyleSheet(_STYLE.replace("__ASSETS__", _ASSETS_DIR.as_posix()))
        self._build_ui()
        self._style_all_combo_popups()
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
        analysis.addAction(
            "Theoretical Yield",
            lambda: self._in_tab("Production", self.run_theoretical_yield),
        )
        analysis.addAction(
            "Production Envelope",
            lambda: self._in_tab("Production", self.run_production_envelope_plot),
        )
        analysis.addAction(
            "FSEOF Targets", lambda: self._in_tab("Production", self.run_fseof_plot)
        )
        analysis.addAction(
            "Strain Design (OptKnock / RobustKnock)…",
            lambda: self._goto_tab("Strain Design"),
        )
        analysis.addSeparator()
        analysis.addAction(
            "Omics Integration (E-Flux2 / LAD)…", lambda: self._goto_tab("Omics")
        )
        analysis.addAction(
            "Revert Metabolism…", lambda: self._goto_tab("Revert Metabolism")
        )
        analysis.addAction(
            "Transformation Targets (A→B)…", lambda: self._goto_tab("Transform (A→B)")
        )
        if self._map_path:
            analysis.addSeparator()
            analysis.addAction(
                "Render Flux Map",
                lambda: self._in_tab("Flux Map", self.render_flux_map),
            )

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
            self,
            "Model Info",
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
            self,
            "Open SBML model",
            "",
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
        self._omics_computed_methods = set()
        self._refresh_comparison_templates()  # disable LAD/E-Flux2 until recomputed
        self._revert_source_expression = None
        self._revert_target_expression = None
        self._transform_source_expression = None
        self._transform_target_expression = None
        self._omics_table_df = None
        self._omics_fluxes_by_condition = None
        self._omics_conditions_order = []
        self.revert_source_label.setText("not loaded")
        self.revert_target_label.setText("not loaded")
        self.transform_source_label.setText("not loaded")
        self.transform_target_label.setText("not loaded")
        self.omics_file_label.setText("no file loaded")
        self.omics_compute_btn.setEnabled(False)
        self.omics_cond_list.clear()
        self.omics_table.setRowCount(0)
        self.transform_table.setRowCount(0)
        self._update_revert_run_state()
        self._update_transform_run_state()
        self.load_model(model)
        self.status_label.setText(f"Loaded model '{model.id}' from {path}.")

    def _active_table(self) -> QTableWidget | None:
        """The primary result table of the current tab, or None for figure-only tabs."""

        name = self.tabs.tabText(self.tabs.currentIndex())
        return {
            "Simulation": self.sim_table,
            "Comparison": self.comparison_table,
            "Production": self.fseof_table,
            "Strain Design": self.sd_table,
            "Omics": self.omics_table,
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
            table.horizontalHeaderItem(c).text()
            if table.horizontalHeaderItem(c)
            else f"col{c}"
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
        self.reaction_table.setHorizontalHeaderLabels(
            ["Reaction", "Lower", "Upper", "Flux"]
        )
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
        hint = QLabel(
            "Double-click a Lower/Upper cell to edit a bound, then re-run FBA."
        )
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
        if self._map_path:
            self.tabs.addTab(self._build_fluxmap_tab(), "Flux Map")
        self.tabs.addTab(self._build_revert_tab(), "Revert Metabolism")
        self.tabs.addTab(self._build_transformation_tab(), "Transform (A→B)")
        return self.tabs

    def _build_omics_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Omics integration (expression → per-condition flux)")
        cbox = QVBoxLayout(controls)

        # Row 1: method + load (loading only stores the file; Compute runs the integration).
        row = QHBoxLayout()
        row.addWidget(QLabel("Method:"))
        self.omics_method_combo = QComboBox()
        self.omics_method_combo.addItems(["eflux2", "lad"])
        self.omics_method_combo.setToolTip(
            "E-Flux2: scale bounds by expression, minimize L2 flux. LAD: fit fluxes to expression."
        )
        row.addWidget(self.omics_method_combo)
        load_btn = QPushButton("Load expression CSV…")
        load_btn.clicked.connect(self.load_expression_csv)
        row.addWidget(load_btn)
        demo_btn = QPushButton("Use demo expression")
        demo_btn.clicked.connect(self.load_omics_demo)
        row.addWidget(demo_btn)
        row.addStretch(1)
        cbox.addLayout(row)

        # Row 2: the loaded filename and Compute (separate so the method can change and re-run
        # on the same loaded data).
        run_row = QHBoxLayout()
        self.omics_file_label = QLabel("no file loaded")
        self.omics_file_label.setStyleSheet("color: #5c6b7e; font-style: italic;")
        run_row.addWidget(self.omics_file_label, 1)
        self.omics_compute_btn = QPushButton("Compute")
        self.omics_compute_btn.setEnabled(False)
        self.omics_compute_btn.setToolTip(
            "Load an expression source first, then compute the checked conditions."
        )
        self.omics_compute_btn.clicked.connect(self.compute_omics)
        run_row.addWidget(self.omics_compute_btn)
        cbox.addLayout(run_row)

        # Row 2: the conditions detected in the file, each check-selectable.
        cond_row = QHBoxLayout()
        cond_row.addWidget(QLabel("Conditions:"))
        self.omics_cond_list = QListWidget()
        self.omics_cond_list.setSelectionMode(QAbstractItemView.NoSelection)
        self.omics_cond_list.setFlow(QListWidget.LeftToRight)
        self.omics_cond_list.setWrapping(True)
        self.omics_cond_list.setMaximumHeight(64)
        self.omics_cond_list.setToolTip(
            "Columns detected in the file. Check the conditions to predict a flux column for."
        )
        cond_row.addWidget(self.omics_cond_list, 1)
        cbox.addLayout(cond_row)
        layout.addWidget(controls)

        self.omics_summary = QLabel(
            "Load a gene × condition CSV (first column = gene id, remaining columns = "
            "conditions) or use demo data, then check conditions and click Compute."
        )
        self.omics_summary.setWordWrap(True)
        layout.addWidget(self.omics_summary)

        # Result-view controls sit with the table so the toggle reads as adjusting the results.
        table_controls = QHBoxLayout()
        table_controls.addWidget(QLabel("Predicted fluxes:"))
        table_controls.addStretch(1)
        self.omics_show_all_check = QCheckBox("Show all reactions")
        self.omics_show_all_check.setToolTip(
            "Off: only reactions with non-zero flux in some condition. On: every reaction."
        )
        self.omics_show_all_check.toggled.connect(self._render_omics_table)
        table_controls.addWidget(self.omics_show_all_check)
        layout.addLayout(table_controls)

        # Columns are set per run: Reaction + one predicted-flux column per computed condition.
        self.omics_table = QTableWidget(0, 1)
        self.omics_table.setHorizontalHeaderLabels(["Reaction"])
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
        controls_layout = QVBoxLayout(controls)

        # Row 1: the product / substrate / oxygen selectors.
        selectors = QHBoxLayout()
        selectors.addWidget(QLabel("Target product:"))
        self.product_combo = QComboBox()
        self.product_combo.setMinimumWidth(160)
        self.product_combo.setToolTip("Exchange reaction of the target product")
        self.product_combo.currentTextChanged.connect(self._on_product_changed)
        selectors.addWidget(self.product_combo, 1)
        selectors.addWidget(QLabel("substrate:"))
        self.substrate_combo = QComboBox()
        self.substrate_combo.setMinimumWidth(140)
        self.substrate_combo.setToolTip(
            "Carbon substrate (auto = detect from the medium)"
        )
        selectors.addWidget(self.substrate_combo, 1)
        self.anaerobic_combo = QComboBox()
        self.anaerobic_combo.addItems(["aerobic", "anaerobic"])
        selectors.addWidget(self.anaerobic_combo)
        controls_layout.addLayout(selectors)

        # Row 2: the analysis actions, kept on their own line so the controls box does not
        # force the whole window wide.
        actions = QHBoxLayout()
        yield_btn = QPushButton("Theoretical yield")
        yield_btn.clicked.connect(self.run_theoretical_yield)
        envelope_btn = QPushButton("Production envelope")
        envelope_btn.clicked.connect(self.run_production_envelope_plot)
        fseof_btn = QPushButton("FSEOF targets")
        fseof_btn.clicked.connect(self.run_fseof_plot)
        fvseof_btn = QPushButton("FVSEOF (robust)")
        fvseof_btn.clicked.connect(self.run_fvseof_plot)
        actions.addWidget(yield_btn)
        actions.addWidget(envelope_btn)
        actions.addWidget(fseof_btn)
        actions.addWidget(fvseof_btn)
        actions.addStretch(1)
        self._production_buttons = [yield_btn, envelope_btn, fseof_btn, fvseof_btn]
        controls_layout.addLayout(actions)
        layout.addWidget(controls)

        self.yield_label = QLabel("Select a product and compute its theoretical yield.")
        self.yield_label.setFont(QFont("", 13, QFont.Bold))
        self.yield_label.setWordWrap(True)
        self.yield_label.setTextFormat(Qt.RichText)
        layout.addWidget(self.yield_label)

        # Plot on the left, its underlying data table on the right — side by side so the
        # figure and the numbers behind it are visible together (equal 50/50 split).
        results = QSplitter(Qt.Horizontal)
        self.production_canvas_holder = QVBoxLayout()
        holder = QWidget()
        holder.setLayout(self.production_canvas_holder)
        results.addWidget(holder)
        self._production_canvas = None
        self._production_toolbar = None

        # Right side: a header row with the FSEOF/FVSEOF "show all" toggle, above the table.
        right_box = QWidget()
        right_layout = QVBoxLayout(right_box)
        right_layout.setContentsMargins(0, 0, 0, 0)
        table_header = QHBoxLayout()
        table_header.addWidget(QLabel("Result table:"))
        table_header.addStretch(1)
        self.production_show_all_check = QCheckBox("Show all reactions")
        self.production_show_all_check.setToolTip(
            "FSEOF/FVSEOF only. Off: amplify/knockdown targets. On: every reaction "
            "(unchanged ones included)."
        )
        self.production_show_all_check.toggled.connect(self._rerender_production_table)
        table_header.addWidget(self.production_show_all_check)
        right_layout.addLayout(table_header)

        self.fseof_table = QTableWidget(0, 4)
        self._set_production_table_headers(
            ["Reaction", "Flux @ low enforce", "Flux @ high enforce", "Trend"]
        )
        self.fseof_table.verticalHeader().setVisible(False)
        self.fseof_table.setAlternatingRowColors(True)
        right_layout.addWidget(self.fseof_table)
        results.addWidget(right_box)
        results.setStretchFactor(0, 1)
        results.setStretchFactor(1, 1)
        results.setSizes([10000, 10000])
        layout.addWidget(results, 1)
        return tab

    def _build_strain_design_tab(self) -> QWidget:
        """Growth-coupled strain design: OptKnock / RobustKnock over the library service."""

        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Growth-coupled strain design (OptKnock / RobustKnock)")
        form = QVBoxLayout(controls)

        # Row 1: target product + method.
        pick_row = QHBoxLayout()
        pick_row.addWidget(QLabel("Target product:"))
        self.sd_product_combo = QComboBox()
        self.sd_product_combo.setMinimumWidth(150)
        self.sd_product_combo.setToolTip(
            "Exchange reaction of the target product to couple to growth"
        )
        pick_row.addWidget(self.sd_product_combo, 1)
        pick_row.addWidget(QLabel("Method:"))
        self.sd_method_combo = QComboBox()
        self.sd_method_combo.addItems(["optknock", "robustknock"])
        self.sd_method_combo.setToolTip(
            "OptKnock maximizes product at max growth (optimistic); RobustKnock keeps only "
            "designs that guarantee product at max growth (worst case)."
        )
        pick_row.addWidget(self.sd_method_combo)
        form.addLayout(pick_row)

        # Row 2: search limits + run.
        limit_row = QHBoxLayout()
        limit_row.addWidget(QLabel("max KOs:"))
        self.sd_max_ko_spin = QSpinBox()
        self.sd_max_ko_spin.setRange(1, 6)
        self.sd_max_ko_spin.setValue(3)
        self.sd_max_ko_spin.setToolTip(
            "Maximum number of reaction knockouts per design"
        )
        limit_row.addWidget(self.sd_max_ko_spin)
        limit_row.addWidget(QLabel("solutions:"))
        self.sd_max_sol_spin = QSpinBox()
        self.sd_max_sol_spin.setRange(1, 20)
        self.sd_max_sol_spin.setValue(5)
        self.sd_max_sol_spin.setToolTip("Maximum number of designs to enumerate")
        limit_row.addWidget(self.sd_max_sol_spin)
        self.sd_run_btn = QPushButton("Run design")
        self.sd_run_btn.clicked.connect(self.run_strain_design)
        limit_row.addWidget(self.sd_run_btn)
        limit_row.addStretch(1)
        form.addLayout(limit_row)
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
        self.sd_summary.setText(
            f"Searching {method} designs for {html.escape(product)}…"
        )
        self.status_label.setText(f"Running {method} (this can take a while)…")
        self.sd_run_btn.setEnabled(False)
        try:
            search = optknock if method == "optknock" else robustknock
            result = self._run_in_background(
                lambda: search(
                    self.model, product, max_knockouts=max_ko, max_solutions=max_sol
                ),
                label=f"Running {method}…",
            )
        except (
            Exception
        ) as exc:  # surface solver-capability/model errors instead of freezing
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
            self.sd_table.setItem(
                i, 1, QTableWidgetItem(", ".join(design.knockouts) or "—")
            )
            self.sd_table.setItem(i, 2, QTableWidgetItem(f"{design.growth:.4g}"))
            self.sd_table.setItem(i, 3, QTableWidgetItem(f"{design.max_product:.4g}"))
            self.sd_table.setItem(
                i, 4, QTableWidgetItem(f"{design.guaranteed_product:.4g}")
            )
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
        controls = QGroupBox(
            "Perturbation response (MOMA / ROOM vs a reference template)"
        )
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
            "Wild-type flux template the perturbed state is compared to. LAD/E-Flux2 reuse the "
            "expression computed on the Omics tab, and stay disabled until that is available."
        )
        # LAD/E-Flux2 need an omics result; grey them out (delegate + not-allowed cursor) until
        # one exists so they cannot be picked.
        self._style_disabled_items(self.template_combo)
        self._refresh_comparison_templates()
        row.addWidget(self.template_combo)
        row.addWidget(QLabel("Knockout level:"))
        self.ko_level_combo = QComboBox()
        self.ko_level_combo.addItems(["reaction", "gene"])
        self.ko_level_combo.setToolTip(
            "Knock out reactions, or genes (resolved to reactions via GPR)"
        )
        self.ko_level_combo.currentTextChanged.connect(self._populate_ko_list)
        row.addWidget(self.ko_level_combo)
        row.addStretch(1)
        controls_layout.addLayout(row)

        # Row 2: the optional batch product and the single-run significance threshold.
        opts = QHBoxLayout()
        opts.addWidget(QLabel("Batch product (optional):"))
        self.batch_product_combo = QComboBox()
        self.batch_product_combo.setMinimumWidth(140)
        self.batch_product_combo.setToolTip(
            "Optional. When set, the Batch table adds the product's wild-type and post-knockout "
            "flux columns."
        )
        opts.addWidget(self.batch_product_combo, 1)
        # Significant-change threshold for the single-run table: a reaction is shown only when
        # |perturbed - reference| exceeds this fraction of |reference| (plus a small floor).
        # This matches ROOM's own significance criterion and hides alternate-optimum drift.
        opts.addWidget(QLabel("Significant change ≥"))
        self.comparison_threshold_spin = QDoubleSpinBox()
        self.comparison_threshold_spin.setRange(0.0, 100.0)
        self.comparison_threshold_spin.setSingleStep(1.0)
        self.comparison_threshold_spin.setValue(3.0)
        self.comparison_threshold_spin.setSuffix(" % of ref")
        self.comparison_threshold_spin.setToolTip(
            "Show a reaction only when |perturbed − reference| > this % of |reference| + 0.001. "
            "0% shows every reaction that moved at all."
        )
        self.comparison_threshold_spin.valueChanged.connect(self._render_comparison_table)
        opts.addWidget(self.comparison_threshold_spin)
        controls_layout.addLayout(opts)

        # Row 3: the knock-out actions.
        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Knock out:"))
        run_btn = QPushButton("Run (selected as one KO)")
        run_btn.setToolTip("MOMA/ROOM with all selected targets knocked out together")
        run_btn.clicked.connect(self.run_comparison)
        batch_btn = QPushButton("Batch (each separately)")
        batch_btn.setToolTip(
            "Run MOMA/ROOM once per target as a separate single knockout — the selected "
            "targets, or all targets of this level when none are selected."
        )
        batch_btn.clicked.connect(self.run_batch_comparison)
        row2.addWidget(run_btn)
        row2.addWidget(batch_btn)
        row2.addStretch(1)
        controls_layout.addLayout(row2)

        # Two panels: left = searchable catalogue of all targets, right = the chosen knockout
        # set. Add/remove moves ids between them, so the selection is always plainly visible.
        picker = QHBoxLayout()

        left = QVBoxLayout()
        self.ko_filter = QLineEdit()
        self.ko_filter.setPlaceholderText("Search targets…")
        self.ko_filter.setClearButtonEnabled(True)
        self.ko_filter.textChanged.connect(self._filter_ko_list)
        left.addWidget(self.ko_filter)
        self.ko_available = QListWidget()
        self.ko_available.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.ko_available.setMaximumHeight(150)
        self.ko_available.setToolTip("Double-click (or select and Add →) to add a target.")
        self.ko_available.itemDoubleClicked.connect(lambda item: self._add_ko_targets([item.text()]))
        left.addWidget(self.ko_available)
        picker.addLayout(left, 1)

        # Add / remove buttons between the two lists.
        mid = QVBoxLayout()
        mid.addStretch(1)
        add_btn = QPushButton("Add →")
        add_btn.setToolTip("Add the highlighted targets to the knockout set")
        add_btn.clicked.connect(self._add_selected_ko)
        remove_btn = QPushButton("← Remove")
        remove_btn.setToolTip("Remove the highlighted targets from the knockout set")
        remove_btn.clicked.connect(self._remove_selected_ko)
        clear_btn = QPushButton("Clear")
        clear_btn.setToolTip("Empty the knockout set")
        clear_btn.clicked.connect(self._clear_ko_selected)
        mid.addWidget(add_btn)
        mid.addWidget(remove_btn)
        mid.addWidget(clear_btn)
        mid.addStretch(1)
        picker.addLayout(mid)

        right = QVBoxLayout()
        self.ko_selected_label = QLabel("Selected (0):")
        right.addWidget(self.ko_selected_label)
        self.ko_selected = QListWidget()
        self.ko_selected.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.ko_selected.setMaximumHeight(150)
        self.ko_selected.setToolTip("The knockout set. Double-click (or select and ← Remove) to drop one.")
        self.ko_selected.itemDoubleClicked.connect(
            lambda item: self._remove_ko_targets([item.text()])
        )
        right.addWidget(self.ko_selected)
        picker.addLayout(right, 1)

        controls_layout.addLayout(picker)
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
        self.comparison_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch
        )
        self.comparison_table.verticalHeader().setVisible(False)
        self.comparison_table.setAlternatingRowColors(True)
        layout.addWidget(self.comparison_table, 1)
        return tab

    def _populate_ko_list(self) -> None:
        """Fill the available-targets catalogue for the current level; clear the chosen set."""

        if not hasattr(self, "ko_available"):
            return
        level = self.ko_level_combo.currentText()
        ids = (
            [g.id for g in self.model.genes] if level == "gene"
            else [r.id for r in self.model.reactions]
        )
        self.ko_available.clear()
        self.ko_available.addItems(ids)
        self.ko_selected.clear()  # switching level invalidates the previous selection
        if hasattr(self, "ko_filter"):
            self._filter_ko_list(self.ko_filter.text())
        self._update_ko_selected_label()

    def _filter_ko_list(self, text: str) -> None:
        """Hide catalogue targets that do not match the search box (case-insensitive)."""

        query = text.strip().lower()
        for i in range(self.ko_available.count()):
            item = self.ko_available.item(i)
            item.setHidden(bool(query) and query not in item.text().lower())

    def _add_ko_targets(self, ids) -> None:
        """Add target ids to the chosen set, skipping ones already present."""

        existing = {
            self.ko_selected.item(i).text() for i in range(self.ko_selected.count())
        }
        for target_id in ids:
            if target_id not in existing:
                self.ko_selected.addItem(target_id)
                existing.add(target_id)
        self._update_ko_selected_label()

    def _remove_ko_targets(self, ids) -> None:
        drop = set(ids)
        for i in range(self.ko_selected.count() - 1, -1, -1):
            if self.ko_selected.item(i).text() in drop:
                self.ko_selected.takeItem(i)
        self._update_ko_selected_label()

    def _add_selected_ko(self) -> None:
        self._add_ko_targets([item.text() for item in self.ko_available.selectedItems()])

    def _remove_selected_ko(self) -> None:
        self._remove_ko_targets([item.text() for item in self.ko_selected.selectedItems()])

    def _clear_ko_selected(self) -> None:
        self.ko_selected.clear()
        self._update_ko_selected_label()

    def _update_ko_selected_label(self) -> None:
        if hasattr(self, "ko_selected_label"):
            self.ko_selected_label.setText(f"Selected ({self.ko_selected.count()}):")

    def _selected_ko_targets(self) -> list[str]:
        return [self.ko_selected.item(i).text() for i in range(self.ko_selected.count())]

    def _style_all_combo_popups(self) -> None:
        """Centre every combo's popup entries (Qt defaults to left-aligned)."""

        for combo in self.findChildren(QComboBox):
            combo.setItemDelegate(_ComboItemDelegate(combo))

    def _style_disabled_items(self, combo: QComboBox) -> None:
        """Make a combo's disabled entries clearly inert: faint painting + not-allowed cursor."""

        combo.setItemDelegate(_ComboItemDelegate(combo))
        view = combo.view()
        viewport = view.viewport()
        viewport.setMouseTracking(True)
        filt = _DisabledHoverCursor(view)
        viewport.installEventFilter(filt)
        self._cursor_filters.append(filt)

    def _refresh_comparison_templates(self) -> None:
        """Enable LAD/E-Flux2 reference templates only once an omics result exists."""

        combo = getattr(self, "template_combo", None)
        if combo is None:
            return
        for method in ("lad", "eflux2"):
            index = combo.findText(method)
            if index < 0:
                continue
            item = combo.model().item(index)
            if item is None:
                continue
            ready = method in self._omics_computed_methods
            item.setEnabled(ready)
            item.setToolTip(
                "" if ready
                else f"Compute {method.upper()} on the Omics tab to use it as a "
                "reference template."
            )
        # If the current pick just became unavailable, fall back to a template that always runs.
        current = combo.currentText()
        if current in ("lad", "eflux2") and current not in self._omics_computed_methods:
            combo.setCurrentText("fba")

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
        self._comparison_cache = None  # drop any prior result until this run succeeds
        method_label = self.comparison_method_combo.currentText()
        template = self.template_combo.currentText()
        level = self.ko_level_combo.currentText()
        targets = self._selected_ko_targets()
        if not targets:
            self.comparison_summary.setText(
                f"Add one or more {level}s to the knockout set (search on the left, then "
                "Add →) before running."
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
            result = knockout_comparison(
                self.model, reference, reaction_ids, method=method_key
            )
            return reference, result, reaction_ids

        try:
            reference, result, reaction_ids = self._run_in_background(
                _compute, label=f"Running {method_label}…"
            )
        except Exception as exc:
            self.comparison_summary.setText(
                f"Comparison failed: {html.escape(str(exc))}"
            )
            return

        synthetic_note = (
            " <span style='color:#b45309'>⚠ synthetic demo expression — load a CSV on the "
            "Omics tab for a data-driven template.</span>"
            if synthetic_expression
            else ""
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
        # Cache the raw solve so the change-threshold can re-filter without re-solving.
        self._comparison_cache = {
            "reference": reference,
            "fluxes": dict(result.fluxes),
            "method_label": method_label,
            "template": template,
            "level": level,
            "ko_label": ko_label,
            "blocked_note": blocked_note,
            "status": result.status,
            "distance": result.distance,
            "synthetic_note": synthetic_note,
        }
        self._render_comparison_table()

    def _render_comparison_table(self) -> None:
        """Fill the single-run table with reactions whose change exceeds the threshold."""

        cache = self._comparison_cache
        if not cache:
            return
        reference = cache["reference"]
        fluxes = cache["fluxes"]
        delta = self.comparison_threshold_spin.value() / 100.0
        epsilon = 1e-3  # absolute floor so reactions near zero reference aren't over-reported
        rows = sorted(
            fluxes.items(), key=lambda kv: -abs(kv[1] - reference.get(kv[0]))
        )
        changed = [
            (r, v)
            for r, v in rows
            if abs(v - reference.get(r)) > delta * abs(reference.get(r)) + epsilon
        ]
        self.comparison_table.setColumnCount(3)
        self.comparison_table.setHorizontalHeaderLabels(
            ["Reaction", "Reference flux", "Perturbed flux"]
        )
        self.comparison_table.setRowCount(len(changed))
        for i, (rid, flux) in enumerate(changed):
            self.comparison_table.setItem(i, 0, QTableWidgetItem(rid))
            self.comparison_table.setItem(
                i, 1, QTableWidgetItem(f"{reference.get(rid):.4g}")
            )
            self.comparison_table.setItem(i, 2, QTableWidgetItem(f"{flux:.4g}"))
        self.comparison_summary.setText(
            f"<b>{cache['method_label']}</b> vs <b>{cache['template']}</b> template after "
            f"knocking out {cache['level']} {html.escape(cache['ko_label'])}"
            f"{cache['blocked_note']}: status {cache['status']}, distance "
            f"{cache['distance']:.4g}, {len(changed)} reactions changed "
            f"(≥ {self.comparison_threshold_spin.value():g}% of reference)."
            + cache["synthetic_note"]
        )
        self.status_label.setText(
            f"Comparison complete ({cache['method_label']}, {cache['template']})."
        )

    def run_batch_comparison(self) -> None:
        """Run MOMA/ROOM once per target (each a separate single knockout) into a batch table."""

        import math

        # The batch reshapes the table; drop the single-run cache so the change-threshold
        # spinbox can't re-render a stale single-run over the batch table.
        self._comparison_cache = None

        method_label = self.comparison_method_combo.currentText()
        template = self.template_combo.currentText()
        level = self.ko_level_combo.currentText()
        selected = self._selected_ko_targets()  # empty -> all targets of this level
        method_key = self._comparison_method_key()
        expression, synthetic_expression = self._comparison_expression()
        product = self.batch_product_combo.currentText()
        product = None if product in ("", "(none)") else product
        objective_id = next(
            (r.id for r in self.model.reactions if r.objective_coefficient != 0), None
        )

        def _compute():
            reference = reference_flux(self.model, template, gene_expression=expression)
            if level == "gene":
                perts = gene_perturbations(self.model, selected or None)
            else:
                perts = reaction_perturbations(self.model, selected or None)
            rows = batch_comparison(
                self.model, reference, perts, method=method_key, product_reaction=product
            )
            return reference, rows

        self.status_label.setText(f"Running batch {method_label} ({level})…")
        try:
            reference, rows = self._run_in_background(
                _compute, label=f"Batch {method_label} over {level}s…"
            )
        except Exception as exc:
            self.comparison_summary.setText(
                f"Batch comparison failed: {html.escape(str(exc))}"
            )
            self.status_label.setText("Batch comparison failed.")
            return

        # Wild-type reference values (constant across knockouts).
        wt_biomass = float(reference.get(objective_id)) if objective_id else float("nan")
        wt_target = float(reference.get(product)) if product else float("nan")

        def _essential(row) -> bool:
            """A knockout is essential if it leaves the cell unable to grow (or infeasible)."""
            if row.status != "optimal" or not math.isfinite(row.objective):
                return True
            if math.isfinite(wt_biomass) and wt_biomass > 1e-9:
                return row.objective < 0.01 * wt_biomass  # under 1% of wild-type growth
            return row.objective < 1e-6

        # Most impactful first: lowest post-knockout growth (essential/lethal) at the top.
        rows = sorted(
            rows,
            key=lambda r: r.objective if math.isfinite(r.objective) else float("-inf"),
        )

        headers = ["Target", "WT Biomass", "KO Biomass", "Essential"]
        if product:
            headers += ["WT target", "KO target"]
        self.comparison_table.setColumnCount(len(headers))
        self.comparison_table.setHorizontalHeaderLabels(headers)
        self.comparison_table.setRowCount(len(rows))
        for i, r in enumerate(rows):
            ko_biomass = "—" if not math.isfinite(r.objective) else f"{r.objective:.4g}"
            cells = [
                r.target_id,
                f"{wt_biomass:.4g}" if math.isfinite(wt_biomass) else "—",
                ko_biomass,
                "yes" if _essential(r) else "no",
            ]
            if product:
                ko_target = "—" if not math.isfinite(r.product_flux) else f"{r.product_flux:.4g}"
                cells += [
                    f"{wt_target:.4g}" if math.isfinite(wt_target) else "—",
                    ko_target,
                ]
            for col, text in enumerate(cells):
                self.comparison_table.setItem(i, col, QTableWidgetItem(text))

        synthetic_note = (
            " <span style='color:#b45309'>⚠ synthetic demo expression.</span>"
            if synthetic_expression
            else ""
        )
        essential_n = sum(1 for r in rows if _essential(r))
        product_note = f" Target product: {html.escape(product)}." if product else ""
        self.comparison_summary.setText(
            f"<b>Batch {method_label}</b> vs <b>{template}</b> over {len(rows)} {level} "
            f"knockout(s): {essential_n} essential (growth lost). Sorted by post-knockout "
            f"biomass (most impactful first).{product_note}" + synthetic_note
        )
        self.status_label.setText(
            f"Batch comparison complete ({method_label}, {level})."
        )

    def _build_simulation_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)
        media_row = QHBoxLayout()
        media_row.addWidget(QLabel("Medium:"))
        self.medium_combo = QComboBox()
        self.medium_combo.addItems(sorted(PRESET_MEDIA))
        self.medium_combo.setToolTip(
            "Preset growth medium (sets exchange uptake bounds)"
        )
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
        self.fva_fraction_spin.setToolTip(
            "FVA: fraction of the optimum objective to hold"
        )
        controls.addWidget(self.fva_fraction_spin)
        # Objective on top, pFBA's minimal-total-flux directly beneath it (previously only in
        # the easy-to-miss status bar).
        obj_box = QVBoxLayout()
        obj_box.setSpacing(0)
        self.objective_label = QLabel("Objective: —")
        self.objective_label.setFont(QFont("", 13, QFont.Bold))
        self.pfba_total_label = QLabel("")
        self.pfba_total_label.setStyleSheet("color: #5c6b7e; font-size: 11px;")
        obj_box.addWidget(self.objective_label)
        obj_box.addWidget(self.pfba_total_label)
        controls.addStretch(1)
        controls.addLayout(obj_box)
        layout.addLayout(controls)

        # FBA and pFBA fluxes live in separate columns so running pFBA adds to — rather than
        # overwrites — the FBA result and the two can be compared side by side.
        self.sim_table = QTableWidget(0, 4)
        self.sim_table.setHorizontalHeaderLabels(["Reaction", "FBA flux", "pFBA flux", "FVA range"])
        _sim_header = self.sim_table.horizontalHeader()
        _sim_header.setSectionResizeMode(0, QHeaderView.Stretch)
        for _c in (1, 2, 3):
            _sim_header.setSectionResizeMode(_c, QHeaderView.ResizeToContents)
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
        self.method_combo.addItems(["rmta", "mta", "rmta_continuous"])
        self.method_combo.setToolTip(
            "rmta: published best/MOMA/worst workflow (MIQP); "
            "mta: published single MTA (MIQP); "
            "rmta_continuous: historical QP heuristic, explicitly not published rMTA."
        )
        self.perturbation_combo = QComboBox()
        self.perturbation_combo.addItems(["gene", "reaction"])
        self.alpha_spin = QDoubleSpinBox()
        self.alpha_spin.setRange(0.0, 1.0)
        self.alpha_spin.setSingleStep(0.01)
        self.alpha_spin.setValue(0.66)
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
        self.revert_run_btn.setToolTip(
            "Load both source and target expression files to enable."
        )
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
        self.batch_product_combo.clear()
        self.batch_product_combo.addItems(["(none)", *exchanges])  # optional in batch
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
            self.yield_label.setText(
                "No exchange reactions in this model — production design unavailable."
            )

        self._fluxes = {}
        self._fba_fluxes = {}
        self._pfba_fluxes = {}
        self._fluxes_stale = False
        if getattr(self, "pfba_total_label", None) is not None:
            self.pfba_total_label.setText("")
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
            self.status_label.setText(
                "Bounds changed — fluxes are stale; click Run FBA."
            )

    def run_fba(self) -> None:
        try:
            solution = fba(self.model)
        except Exception as exc:  # surface solver/model errors instead of freezing
            self.objective_label.setText("Objective: error")
            self.status_label.setText(f"FBA failed: {exc}")
            return
        self._fluxes = dict(solution.fluxes)
        self._fba_fluxes = dict(solution.fluxes)
        self._fluxes_stale = False
        obj = solution.objective_value
        obj_text = f"{obj:.4g}" if obj is not None else "infeasible"
        self.objective_label.setText(f"Objective: {obj_text} ({solution.status})")
        self.pfba_total_label.setText("")  # FBA does not minimise total flux
        if obj is None or solution.status != "optimal":
            self.status_label.setText(
                f"FBA status: {solution.status} — model may be infeasible."
            )
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
        self._pfba_fluxes = dict(solution.fluxes)
        self._fluxes_stale = False
        # pFBA's objective value is the minimal total flux; show the growth (objective rxn).
        growth = next(
            (
                self._fluxes[r.id]
                for r in self.model.reactions
                if r.objective_coefficient != 0
            ),
            None,
        )
        growth_text = f"{growth:.4g}" if growth is not None else "—"
        total = sum(abs(v) for v in self._fluxes.values())
        self.objective_label.setText(f"Objective: {growth_text} (pFBA)")
        self.pfba_total_label.setText(f"Minimal total flux: {total:.1f}")
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
        self._fill_sim_table(fva_ranges=ranges, fraction=fraction)
        self.status_label.setText(f"FVA complete (fraction {fraction:g}).")

    def _fill_sim_table(self, *, fva_ranges=None, fraction=None) -> None:
        """Fill the simulation table with FBA/pFBA flux columns (and an FVA range if given)."""

        fva_header = "FVA range" if fraction is None else f"FVA range (f={fraction:g})"
        self.sim_table.setHorizontalHeaderLabels(
            ["Reaction", "FBA flux", "pFBA flux", fva_header]
        )
        ids = sorted(set(self._fba_fluxes) | set(self._pfba_fluxes))
        self.sim_table.setRowCount(len(ids))
        for row, rid in enumerate(ids):
            self.sim_table.setItem(row, 0, QTableWidgetItem(rid))
            fba = self._fba_fluxes.get(rid)
            pfba = self._pfba_fluxes.get(rid)
            self.sim_table.setItem(row, 1, QTableWidgetItem(f"{fba:.3g}" if fba is not None else "—"))
            self.sim_table.setItem(row, 2, QTableWidgetItem(f"{pfba:.3g}" if pfba is not None else "—"))
            if fva_ranges is not None:
                rng = fva_ranges.get(rid)
                text = f"[{rng.minimum:.3g}, {rng.maximum:.3g}]" if rng else "—"
            else:
                text = "—"
            self.sim_table.setItem(row, 3, QTableWidgetItem(text))

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
            fig = escher_flux_map(
                self._map_path, self._fluxes, title=f"{self.model.id} — flux map"
            )
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

    def _set_production_table_headers(self, labels: list[str], tooltips=None) -> None:
        """Re-shape the shared production results table for whichever analysis just ran."""

        self.fseof_table.setColumnCount(len(labels))
        self.fseof_table.setHorizontalHeaderLabels(labels)
        self.fseof_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        if tooltips:
            for i, tip in enumerate(tooltips):
                header_item = self.fseof_table.horizontalHeaderItem(i)
                if header_item is not None and tip:
                    header_item.setToolTip(tip)

    def _rerender_production_table(self) -> None:
        """Re-render the cached FSEOF/FVSEOF result when the show-all toggle changes."""

        cached = self._production_result
        if not cached:
            return
        kind, result = cached
        if kind == "fseof":
            self._fill_fseof_table(result)
        elif kind == "fvseof":
            self._fill_fvseof_table(result)

    def _fill_fseof_table(self, result) -> None:
        """Show FSEOF targets ranked by change in |flux|; all reactions when 'show all' is on."""

        self._production_result = ("fseof", result)
        show_all = self.production_show_all_check.isChecked()
        levels = list(result.enforced_levels)
        lo, hi = levels[0], levels[-1]
        rows = []
        for rid in result.trends.index:
            cls = str(result.trends.loc[rid, "classification"])
            if not show_all and cls not in ("amplify", "knockdown"):
                continue
            flo = float(result.trends.loc[rid, lo])
            fhi = float(result.trends.loc[rid, hi])
            rows.append((rid, flo, fhi, cls, abs(abs(fhi) - abs(flo))))
        rows.sort(key=lambda r: -r[4])
        self._set_production_table_headers(
            ["Reaction", "Flux @ low enforce", "Flux @ high enforce", "Trend"],
            [
                "Reaction ID",
                f"Reaction's flux when the product exchange is fixed at the lowest scanned "
                f"level ({lo:.3g}, ~10% of the theoretical maximum) and biomass is maximised",
                f"Reaction's flux when the product exchange is fixed at the highest scanned "
                f"level ({hi:.3g}, ~90% of the theoretical maximum) and biomass is maximised",
                "amplify = |flux| rises as the cell is pushed toward product (over-express "
                "target); knockdown = |flux| falls (down-regulate/delete target)",
            ],
        )
        self.fseof_table.setRowCount(len(rows))
        for i, (rid, flo, fhi, cls, _delta) in enumerate(rows):
            self.fseof_table.setItem(i, 0, QTableWidgetItem(rid))
            self.fseof_table.setItem(i, 1, QTableWidgetItem(f"{flo:.4g}"))
            self.fseof_table.setItem(i, 2, QTableWidgetItem(f"{fhi:.4g}"))
            self.fseof_table.setItem(i, 3, QTableWidgetItem(cls))

    def _fill_fvseof_table(self, result) -> None:
        """Show FVSEOF targets (robust flagged); all reactions when 'show all' is on."""

        self._production_result = ("fvseof", result)
        show_all = self.production_show_all_check.isChecked()
        levels = list(result.enforced_levels)
        lo, hi = levels[0], levels[-1]
        robust_ids = set(result.robust_targets())
        rows = []
        for rid in result.classification.index:
            cls = str(result.classification.loc[rid])
            if not show_all and cls not in ("amplify", "knockdown"):
                continue
            mlo = float(result.mean.loc[rid, lo])
            mhi = float(result.mean.loc[rid, hi])
            is_robust = rid in robust_ids
            tag = f"{cls} · robust" if is_robust else cls
            rows.append((rid, mlo, mhi, tag, is_robust, abs(abs(mhi) - abs(mlo))))
        rows.sort(key=lambda r: (-int(r[4]), -r[5]))
        self._set_production_table_headers(
            ["Reaction", "Mean flux @ low", "Mean flux @ high", "Trend"],
            [
                "Reaction ID",
                f"Midpoint of the reaction's FVA flux range at the lowest enforced product "
                f"level ({lo:.3g})",
                f"Midpoint of the reaction's FVA flux range at the highest enforced product "
                f"level ({hi:.3g})",
                "amplify/knockdown from the |mean-flux| trend; 'robust' = the FVA minimum "
                "|flux| also rises, so the reaction cannot avoid carrying more flux",
            ],
        )
        self.fseof_table.setRowCount(len(rows))
        for i, (rid, mlo, mhi, tag, _robust, _delta) in enumerate(rows):
            self.fseof_table.setItem(i, 0, QTableWidgetItem(rid))
            self.fseof_table.setItem(i, 1, QTableWidgetItem(f"{mlo:.4g}"))
            self.fseof_table.setItem(i, 2, QTableWidgetItem(f"{mhi:.4g}"))
            self.fseof_table.setItem(i, 3, QTableWidgetItem(tag))

    def _fill_envelope_table(self, envelope) -> None:
        """List the production-envelope points: growth range attainable at each product flux."""

        self._production_result = ("envelope", envelope)  # show-all toggle does not apply here
        frame = envelope.to_frame()
        self._set_production_table_headers(
            ["Product flux", "Growth (min)", "Growth (max)"],
            [
                "Enforced product secretion rate along the envelope",
                "Minimum growth rate feasible at this product flux (lower boundary)",
                "Maximum growth rate feasible at this product flux (upper boundary)",
            ],
        )
        self.fseof_table.setRowCount(len(frame))
        for i, (_, record) in enumerate(frame.iterrows()):
            self.fseof_table.setItem(i, 0, QTableWidgetItem(f"{record['product_flux']:.4g}"))
            self.fseof_table.setItem(i, 1, QTableWidgetItem(f"{record['growth_min']:.4g}"))
            self.fseof_table.setItem(i, 2, QTableWidgetItem(f"{record['growth_max']:.4g}"))

    def _run_production(self, action) -> None:
        """Run a production analysis, surfacing any error instead of crashing the UI."""

        try:
            action()
        except Exception as exc:
            self.yield_label.setText(
                f"Production analysis failed: {html.escape(str(exc))}"
            )
            self.status_label.setText(f"Production analysis failed: {exc}")

    def run_theoretical_yield(self) -> None:
        def _do():
            self.fseof_table.setRowCount(0)  # yield is a scalar, not a target scan
            self._production_result = None  # nothing for the show-all toggle to re-render
            product = self._current_product()
            result = theoretical_yield(
                self.model,
                product,
                substrate=self._current_substrate(),
                aerobic=self._is_aerobic(),
            )
            ceiling = (
                f", carbon ceiling {result.carbon_ceiling:.2f}"
                if result.carbon_ceiling is not None
                else ""
            )
            co2 = (
                f"; needs net CO₂ fixation ({result.co2_exchange:.1f})"
                if result.co2_fixed and result.exceeds_carbon_ceiling
                else ""
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
                lambda: production_envelope(
                    self.model, product, substrate=substrate, aerobic=aerobic, points=20
                ),
                label="Computing production envelope…",
            )
            fig = production_envelope_figure(
                envelope, title=f"Production envelope — {product} ({condition})"
            )
            self._set_production_figure(fig)
            self._fill_envelope_table(envelope)
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
            self._fill_fseof_table(result)
            # Rank the listed targets the same way the figure does (by flux increase).
            levels = list(result.enforced_levels)
            ranked = sorted(
                result.amplification_targets(),
                key=lambda rid: (
                    abs(result.trends.loc[rid, levels[-1]])
                    - abs(result.trends.loc[rid, levels[0]])
                ),
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
            fig = fvseof_figure(
                result, top_n=5, title=f"FVSEOF robust targets — {product}"
            )
            self._set_production_figure(fig)
            self._fill_fvseof_table(result)
            robust = result.robust_targets()[:8]
            self.yield_label.setText(
                "FVSEOF robust amplification targets (forced up): "
                + html.escape(", ".join(robust) if robust else "none")
            )
            self.status_label.setText("FVSEOF complete.")

        self._run_production(_do)

    # -- omics integration --------------------------------------------------

    def _set_omics_source(self, table, label: str) -> None:
        """Store a loaded gene × condition table and list its conditions (does not compute)."""

        self._omics_table_df = table
        # Discard any previous result so the show-all toggle can't re-render stale data, and
        # so a new source's methods must be recomputed before they back a Comparison template.
        self._omics_fluxes_by_condition = None
        self._omics_conditions_order = []
        self._omics_expression = None
        self._omics_computed_methods = set()
        self._refresh_comparison_templates()
        self.omics_table.setRowCount(0)
        conditions = [str(c) for c in table.columns]
        self.omics_file_label.setText(
            f"{label} ({table.shape[0]} genes × {len(conditions)} condition(s))"
        )
        self.omics_cond_list.blockSignals(True)
        self.omics_cond_list.clear()
        for cond in conditions:
            item = QListWidgetItem(cond)
            item.setFlags((item.flags() | Qt.ItemIsUserCheckable) & ~Qt.ItemIsSelectable)
            item.setCheckState(Qt.Checked)  # compute all conditions by default
            self.omics_cond_list.addItem(item)
        self.omics_cond_list.blockSignals(False)
        self.omics_compute_btn.setEnabled(True)
        self.omics_summary.setText(
            f"Loaded <b>{html.escape(label)}</b>: {len(conditions)} condition(s) "
            f"({html.escape(', '.join(conditions))}). Check conditions and click Compute."
        )
        self.status_label.setText(f"Loaded expression source: {label}.")

    def _checked_omics_conditions(self) -> list[str]:
        return [
            self.omics_cond_list.item(i).text()
            for i in range(self.omics_cond_list.count())
            if self.omics_cond_list.item(i).checkState() == Qt.Checked
        ]

    def load_omics_demo(self) -> None:
        """Load a deterministic single-condition demo table (does not compute)."""

        import numpy as np
        import pandas as pd

        genes = [g.id for g in self.model.genes]
        if not genes:
            self.omics_summary.setText(
                "Model has no genes — omics integration is not applicable."
            )
            return
        rng = np.random.default_rng(0)
        values = {g: float(rng.uniform(1.0, 100.0)) for g in genes}
        table = pd.DataFrame({"demo": pd.Series(values)})
        self._set_omics_source(table, "demo expression")

    def load_expression_csv(self) -> None:
        """Load a gene × condition expression CSV/TSV as the source (compute is separate)."""

        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load gene-expression CSV",
            "",
            "CSV/TSV files (*.csv *.tsv);;All files (*)",
        )
        if not path:
            return
        try:
            table = read_expression_table(path)
        except Exception as exc:
            self.omics_summary.setText(
                f"Could not read expression CSV: {html.escape(str(exc))}"
            )
            return
        import os

        self._set_omics_source(table, os.path.basename(path))

    def run_omics_demo(self) -> None:
        """Load the demo expression and compute immediately (used by scenarios/tests)."""

        self.load_omics_demo()
        if self._omics_table_df is not None:
            self.compute_omics()

    def compute_omics(self) -> None:
        """Predict a flux column for each checked condition (E-Flux2 or LAD)."""

        import pandas as pd

        if self._omics_table_df is None:
            self.omics_summary.setText("Load an expression source (CSV or demo) first.")
            return
        conditions = self._checked_omics_conditions()
        if not conditions:
            self.omics_summary.setText("Check at least one condition to compute.")
            return
        method = self.omics_method_combo.currentText()
        table = self._omics_table_df

        def _compute():
            return predict_condition_fluxes(
                self.model, table, method=method, conditions=conditions
            )

        try:
            predictions = self._run_in_background(
                _compute, label=f"Computing {method} for {len(conditions)} condition(s)…"
            )
        except Exception as exc:
            self.omics_summary.setText(
                f"Omics integration failed: {html.escape(str(exc))}"
            )
            self.status_label.setText("Omics integration failed.")
            return

        # Keep only conditions that produced an optimal flux state; note the rest.
        per_condition: dict[str, dict[str, float]] = {}
        failed: list[str] = []
        for cond in conditions:
            try:
                per_condition[cond] = predictions.fluxes(cond)
            except Exception:
                failed.append(cond)
        ok = [c for c in conditions if c in per_condition]
        if not ok:
            self._omics_fluxes_by_condition = None
            self._omics_conditions_order = []
            self.omics_table.setColumnCount(1)
            self.omics_table.setHorizontalHeaderLabels(["Reaction"])
            self.omics_table.setRowCount(0)
            self.omics_summary.setText(
                f"{method.upper()} produced no optimal flux state for the selected conditions."
            )
            return

        # Cache the result so the show-all toggle can re-render without recomputing.
        self._omics_fluxes_by_condition = per_condition
        self._omics_conditions_order = ok
        n_active = len(
            {rid for cond in ok for rid, v in per_condition[cond].items() if abs(v) > 1e-6}
        )
        self._render_omics_table()

        # Reuse the first computed condition as the Comparison tab's LAD/E-Flux2 reference,
        # and record this method so only it is offered as a reference template.
        first = ok[0]
        self._omics_expression = {
            str(gene): float(value)
            for gene, value in table[first].items()
            if pd.notna(value)
        }
        self._omics_computed_methods.add(method)
        self._refresh_comparison_templates()  # this method is now available as a reference
        note = (
            f" ({len(failed)} infeasible: {html.escape(', '.join(failed))})" if failed else ""
        )
        self.omics_summary.setText(
            f"<b>{method.upper()}</b>: {len(ok)} condition(s) × {n_active} active "
            f"reactions{note}. Reference for Comparison = '{html.escape(first)}'."
        )
        self.status_label.setText(f"Omics integration complete ({method}).")

    def _render_omics_table(self) -> None:
        """Render the cached per-condition fluxes, honouring the active-only / show-all toggle."""

        per = self._omics_fluxes_by_condition
        ok = self._omics_conditions_order
        if not per or not ok:
            return
        if self.omics_show_all_check.isChecked():
            rxn_set = {rid for cond in ok for rid in per[cond]}
        else:
            rxn_set = {
                rid for cond in ok for rid, v in per[cond].items() if abs(v) > 1e-6
            }
        rxns = sorted(
            rxn_set, key=lambda rid: -max(abs(per[c].get(rid, 0.0)) for c in ok)
        )
        self.omics_table.setColumnCount(1 + len(ok))
        self.omics_table.setHorizontalHeaderLabels(["Reaction", *ok])
        self.omics_table.setRowCount(len(rxns))
        for i, rid in enumerate(rxns):
            self.omics_table.setItem(i, 0, QTableWidgetItem(rid))
            for j, cond in enumerate(ok, start=1):
                self.omics_table.setItem(
                    i, j, QTableWidgetItem(f"{per[cond].get(rid, 0.0):.4g}")
                )

    def _update_revert_run_state(self) -> None:
        ready = (
            self._revert_source_expression is not None
            and self._revert_target_expression is not None
        )
        self.revert_run_btn.setEnabled(ready)

    def load_revert_expression(self, role: str) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            f"Load {role} expression",
            "",
            "CSV/TSV files (*.csv *.tsv);;All files (*)",
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
        import os

        name = os.path.basename(path)
        if role == "source":
            self._revert_source_expression = expression
            self.revert_source_label.setText(f"{name} ({len(expression)} genes)")
        else:
            self._revert_target_expression = expression
            self.revert_target_label.setText(f"{name} ({len(expression)} genes)")
        self._update_revert_run_state()
        self.status_label.setText(
            f"Loaded {role} expression ({len(expression)} genes)."
        )

    def run_loaded_revert(self) -> None:
        if (
            self._revert_source_expression is None
            or self._revert_target_expression is None
        ):
            self.revert_summary.setText(
                "Load both source and target expression files first."
            )
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
            source_result = integrate_expression(
                self.model,
                source_expression,
                method="eflux2",
                objective_fraction=1.0,
            )
            if source_result.status != "optimal" or not source_result.fluxes:
                raise ValueError(
                    "source expression could not produce a valid source-state flux "
                    f"({source_result.status})"
                )
            reference = source_result.to_flux_state("source_expression_state")
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
            ranking = self._run_in_background(
                _compute, label=f"Running revert ({method})…"
            )
        except (
            Exception
        ) as exc:  # surface solver-capability/model errors instead of freezing
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
            "mta: published single MTA MIQP on the source→target direction."
        )
        self.transform_omics_combo = QComboBox()
        self.transform_omics_combo.addItems(["eflux2", "lad"])
        self.transform_omics_combo.setToolTip(
            "How each expression vector is turned into a flux state"
        )
        self.transform_perturbation_combo = QComboBox()
        self.transform_perturbation_combo.addItems(["gene", "reaction"])
        self.transform_alpha_spin = QDoubleSpinBox()
        self.transform_alpha_spin.setRange(0.0, 1.0)
        self.transform_alpha_spin.setSingleStep(0.01)
        self.transform_alpha_spin.setValue(0.66)
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
        self.transform_run_btn.setToolTip(
            "Load both source and target expression files to enable."
        )
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
        self.transform_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch
        )
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
            self,
            f"Load {role} expression",
            "",
            "CSV/TSV files (*.csv *.tsv);;All files (*)",
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
        import os

        name = os.path.basename(path)
        if role == "source":
            self._transform_source_expression = expression
            self.transform_source_label.setText(f"{name} ({len(expression)} genes)")
        else:
            self._transform_target_expression = expression
            self.transform_target_label.setText(f"{name} ({len(expression)} genes)")
        self._update_transform_run_state()
        self.status_label.setText(
            f"Loaded {role} (A→B) expression ({len(expression)} genes)."
        )

    def run_transformation(self) -> None:
        if (
            self._transform_source_expression is None
            or self._transform_target_expression is None
        ):
            self.transform_summary.setText(
                "Load both source (A) and target (B) expression files first."
            )
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
                self.model,
                source_state,
                target_state,
                method=method,
                perturbation=perturbation,
                alpha=alpha,
            )

        try:
            ranking = self._run_in_background(
                _compute, label=f"Running transformation ({method})…"
            )
        except (
            Exception
        ) as exc:  # surface solver-capability/model errors instead of freezing
            self.transform_summary.setText(
                f"Transformation failed: {html.escape(str(exc))}"
            )
            self.status_label.setText(f"Transformation failed ({method}).")
            return
        rows = ranking.to_records()
        self.transform_table.setRowCount(len(rows))
        for i, record in enumerate(rows):
            self.transform_table.setItem(i, 0, QTableWidgetItem(str(record["rank"])))
            self.transform_table.setItem(i, 1, QTableWidgetItem(record["target_id"]))
            self.transform_table.setItem(
                i, 2, QTableWidgetItem(f"{record['score']:.4g}")
            )
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
