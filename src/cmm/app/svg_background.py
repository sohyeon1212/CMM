"""Rasterise an SVG so a flux map can be drawn over a picture of the same network.

CNApy's own map format is a background drawing with reaction boxes positioned on top of it,
and the same idea suits an Escher map: the JSON gives the coordinates, and an SVG exported
from that map gives Escher's drawing of it — arrowheads, node sizes, its own typography —
which CMM's renderer does not reproduce. Colouring by flux then goes on top.

This lives in ``cmm.app`` because it needs a Qt renderer and ``cmm.visualization`` stays free
of Qt. The output is a plain RGBA array, so the figure code takes no Qt dependency from it.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from qtpy.QtCore import Qt
from qtpy.QtGui import QImage, QPainter
from qtpy.QtSvg import QSvgRenderer

__all__ = ["svg_background", "MAX_BACKGROUND_PIXELS"]

#: Longest edge of the rasterised background, in pixels. Enough to stay sharp when a reader
#: zooms a 300 dpi export; beyond this the figure costs memory it cannot show.
MAX_BACKGROUND_PIXELS = 4000


def svg_background(
    path: str | Path, max_pixels: int = MAX_BACKGROUND_PIXELS
) -> np.ndarray:
    """Render ``path`` to an RGBA array at the SVG's own aspect ratio.

    Raises ``ValueError`` if the file is not a renderable SVG, rather than returning a blank
    image that would silently look like a map with nothing on it.
    """

    renderer = QSvgRenderer(str(path))
    if not renderer.isValid():
        raise ValueError(f"Not a readable SVG: {path}")

    size = renderer.defaultSize()
    width, height = size.width(), size.height()
    if width <= 0 or height <= 0:
        raise ValueError(f"SVG declares no usable size: {path}")

    scale = min(1.0, max_pixels / max(width, height))
    width, height = max(1, round(width * scale)), max(1, round(height * scale))

    image = QImage(width, height, QImage.Format_RGBA8888)
    image.fill(Qt.transparent)  # keep the figure's white paper behind the drawing
    painter = QPainter(image)
    try:
        renderer.render(painter)
    finally:
        painter.end()  # a QPainter left open on a QImage crashes on some Qt builds

    buffer = image.constBits()
    buffer.setsize(image.sizeInBytes())
    # Copy: the array must outlive the QImage that owns the buffer.
    return np.frombuffer(buffer, dtype=np.uint8).reshape(height, width, 4).copy()
