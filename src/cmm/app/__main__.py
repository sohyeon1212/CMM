"""Launch the CMM desktop workbench: ``python -m cmm.app [model.xml]``.

Loads an SBML model from the first argument, or the public e_coli_core textbook model when no
path is given, and opens the workbench window.
"""

from __future__ import annotations

import sys


def main() -> int:
    from qtpy.QtWidgets import QApplication

    from cmm.app.main_window import CmmMainWindow

    if len(sys.argv) > 1:
        import cobra

        model = cobra.io.read_sbml_model(sys.argv[1])
    else:
        from cobra.io import load_model

        model = load_model("textbook")

    app = QApplication.instance() or QApplication(sys.argv)
    window = CmmMainWindow(model)
    window.show()
    return app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())
