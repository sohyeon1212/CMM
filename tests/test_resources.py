"""The bundled Escher map: that it is intact, attributed, and offered to the right models."""

from __future__ import annotations

import hashlib
import json
import re

import cobra
from cmm.resources import (
    BUNDLED_MAPS,
    bundled_map_for,
    map_coverage,
    map_reaction_ids,
)

CORE_MAP = BUNDLED_MAPS[0]
ATTRIBUTION = CORE_MAP.parent / "ATTRIBUTION.md"


def test_bundled_map_is_byte_identical_to_its_attributed_digest():
    """The file must stay exactly what ATTRIBUTION.md says CMM redistributes.

    CMM redistributes someone else's work under their license. If the bytes drift from the
    recorded digest the attribution silently becomes false, so the digest is a test, not a
    comment.
    """

    recorded = re.search(r"`([0-9a-f]{64})`", ATTRIBUTION.read_text())
    assert recorded, "ATTRIBUTION.md records no SHA-256 for the bundled map"
    actual = hashlib.sha256(CORE_MAP.read_bytes()).hexdigest()
    assert actual == recorded.group(1)


def test_attribution_names_the_license_and_the_paper_to_cite():
    text = ATTRIBUTION.read_text()
    assert "MIT" in text
    assert "Regents of the University of California" in text
    assert "10.1371/journal.pcbi.1004321" in text  # King et al. 2015, the Escher paper
    assert "escher.github.io" in text


def test_map_is_a_schema_1_escher_document():
    document = json.loads(CORE_MAP.read_text(encoding="utf-8"))
    assert document[0]["schema"].endswith("1-0-0#")
    assert len(map_reaction_ids(str(CORE_MAP))) == 95


def test_core_map_is_offered_for_the_core_model(ecoli_core):
    assert map_coverage(CORE_MAP, ecoli_core) > 0.95
    assert bundled_map_for(ecoli_core) == str(CORE_MAP)


def test_no_map_is_offered_for_an_unrelated_model():
    """A map of a different organism must not be drawn as if it described this model."""

    model = cobra.Model("unrelated")
    a, b = cobra.Metabolite("a_c"), cobra.Metabolite("b_c")
    reaction = cobra.Reaction("MADE_UP_RXN")
    reaction.add_metabolites({a: -1, b: 1})
    model.add_reactions([reaction])

    assert map_coverage(CORE_MAP, model) == 0.0
    assert bundled_map_for(model) is None


def _draw(fig, width_in, height_in, dpi=100):
    """Render at a given panel size the way a resized GUI canvas would."""

    from matplotlib.backends.backend_agg import FigureCanvasAgg

    FigureCanvasAgg(fig)
    fig.set_dpi(dpi)
    fig.set_size_inches(width_in, height_in)
    fig.canvas.draw()
    return fig.canvas.get_renderer()


def _flux_map(ecoli_core):
    from cmm.core import fba
    from cmm.visualization import escher_flux_map

    return escher_flux_map(
        str(CORE_MAP), dict(fba(ecoli_core).fluxes), title="e_coli_core — flux map"
    )


def test_flux_map_survives_a_wide_panel(ecoli_core):
    """The GUI stretches the figure to its panel; the layout has to hold at that shape.

    Three regressions live here, all invisible at the authored figure size and all obvious in
    the application: a title clipped off the top edge, the map drifting right because
    ``colorbar`` re-anchors its parent axes, and a colorbar squeezed to a hairline.
    """

    fig = _flux_map(ecoli_core)
    width, height = 12.42, 6.87  # the panel of a maximised window
    renderer = _draw(fig, width, height)
    ax, cax = fig.axes[0], fig.axes[1]

    title = ax.title.get_window_extent(renderer)
    assert title.y1 <= height * 100, "the title is clipped against the top of the panel"

    box, bar = ax.get_position(), cax.get_position()
    centre = (box.x0 + box.x1) / 2
    # Centred in the room left over once the colorbar and its label have taken theirs.
    assert abs(centre - bar.x0 / 2) < 0.06, f"map is off-centre at x={centre:.3f}"

    long_to_short = (bar.height * height) / (bar.width * width)
    assert long_to_short < 16, f"colorbar is a hairline (1:{long_to_short:.0f})"


def test_flux_map_layout_is_stable_across_panel_sizes(ecoli_core):
    """Re-solved on every draw, so a resize must not move the map or reshape the bar."""

    centres, bars = [], []
    for size in [(12.0, 9.84), (12.42, 6.87), (8.0, 6.0)]:
        fig = _flux_map(ecoli_core)
        _draw(fig, *size)
        box, bar = fig.axes[0].get_position(), fig.axes[1].get_position()
        centres.append((box.x0 + box.x1) / 2)
        bars.append((bar.height * size[1]) / (bar.width * size[0]))

    assert max(centres) - min(centres) < 0.08, f"map centre drifts: {centres}"
    assert max(bars) - min(bars) < 4, f"colorbar aspect drifts: {bars}"
