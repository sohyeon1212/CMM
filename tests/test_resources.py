"""The bundled Escher map: that it is intact, attributed, and offered to the right models."""

from __future__ import annotations

import hashlib
import json
import re

import cobra
import pytest
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


def test_attribution_records_the_chain_the_licence_claim_rests_on():
    """Naming a licence is not enough — the file has to show why that licence applies.

    The map is served from a URL, and a URL carries no licence. What makes redistribution
    sound is that the served bytes are the bytes of a file tracked in a repository whose
    LICENSE is MIT, and each link of that chain is recorded here rather than assumed.
    """

    raw = ATTRIBUTION.read_text()
    # The notice is a wrapped block quote, so match against it with the quote markers and
    # line breaks flattened rather than against however it happens to be laid out today.
    text = " ".join(raw.replace("\n>", " ").split())
    assert "escher/escher.github.io" in text  # the repository the licence comes from
    assert "1-0-0/6/maps" in text  # the path the file is tracked at within it
    assert "byte-identical" in text  # why that repository's file is this file
    assert "MIT" in text
    assert "Copyright © 2019 The Regents of the University of California" in text
    assert "BiGG" in text  # the same map under terms that would not permit this
    assert "10.1371/journal.pcbi.1004321" in text  # King et al. 2015, the Escher paper


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


# --- map background -------------------------------------------------------------------
# The JSON gives the layout; a drawing exported from that same map gives Escher's rendering
# of it, which this renderer does not reproduce. Laying one under the other is optional.


def _probe_svg(tmp_path, width=400, height=200):
    path = tmp_path / "probe.svg"
    path.write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}"><rect width="{width}" height="{height}" '
        'fill="#204080"/></svg>'
    )
    return path


def test_svg_background_rasterises_at_the_svg_aspect_ratio(tmp_path):
    from cmm.app.svg_background import svg_background

    image = svg_background(_probe_svg(tmp_path, 400, 200))
    assert (
        image.shape[2] == 4
    )  # RGBA: the drawing must not paint over the figure's paper
    assert image.shape[1] / image.shape[0] == 2.0
    assert tuple(image[image.shape[0] // 2, image.shape[1] // 2][:3]) == (
        0x20,
        0x40,
        0x80,
    )


def test_svg_background_caps_its_own_size(tmp_path):
    from cmm.app.svg_background import MAX_BACKGROUND_PIXELS, svg_background

    image = svg_background(_probe_svg(tmp_path, 20000, 10000))
    assert max(image.shape[:2]) == MAX_BACKGROUND_PIXELS


def test_svg_background_refuses_a_file_that_is_not_an_svg(tmp_path):
    """A blank image would look like a map with nothing drawn on it — say so instead."""

    from cmm.app.svg_background import svg_background

    junk = tmp_path / "not.svg"
    junk.write_text("this is not markup")
    with pytest.raises(ValueError, match="Not a readable SVG"):
        svg_background(junk)


def test_flux_map_places_a_background_in_the_maps_own_coordinates(ecoli_core):
    """The canvas the layout was drawn on is what a drawing of it lines up with."""

    import json

    import numpy as np

    from cmm.core import fba
    from cmm.visualization import escher_flux_map

    canvas = json.loads(CORE_MAP.read_text(encoding="utf-8"))[1]["canvas"]
    # A drawing with the canvas's own proportions — an export of this map is exactly that.
    height = 400
    width = round(height * canvas["width"] / canvas["height"])
    drawing = np.zeros((height, width, 4), dtype=np.uint8)

    fig = escher_flux_map(
        str(CORE_MAP), dict(fba(ecoli_core).fluxes), background=drawing
    )
    images = fig.axes[0].images
    assert len(images) == 1
    left, right, bottom, top = images[0].get_extent()
    assert left == pytest.approx(canvas["x"], abs=1.0)
    assert right == pytest.approx(canvas["x"] + canvas["width"], abs=1.0)
    # Escher's y grows downward, matching how the axes are inverted.
    assert top == pytest.approx(canvas["y"], abs=1.0)
    assert bottom == pytest.approx(canvas["y"] + canvas["height"], abs=1.0)
    assert images[0].zorder < 1  # under the flux strokes, not over them
    assert images[0].axes.get_aspect() == 1.0  # imshow must not undo the equal aspect


def test_a_background_is_never_stretched_out_of_proportion(ecoli_core):
    """Any picture keeps its own proportions; a mismatch shows as margin, not as distortion."""

    import numpy as np

    from cmm.core import fba
    from cmm.visualization import escher_flux_map

    fluxes = dict(fba(ecoli_core).fluxes)
    for height, width in [(400, 400), (100, 900), (900, 100)]:
        fig = escher_flux_map(
            str(CORE_MAP),
            fluxes,
            background=np.zeros((height, width, 4), dtype=np.uint8),
        )
        left, right, bottom, top = fig.axes[0].images[0].get_extent()
        drawn = abs(right - left) / abs(bottom - top)
        assert drawn == pytest.approx(width / height, rel=1e-6), (
            f"{width}x{height} drawn at aspect {drawn:.4f}"
        )


def test_flux_map_says_it_cannot_rasterise_svg_itself(ecoli_core, tmp_path):
    """cmm.visualization stays free of Qt, so it must not pretend to render SVG."""

    from cmm.core import fba
    from cmm.visualization import escher_flux_map

    with pytest.raises(ValueError, match="cannot rasterise SVG"):
        escher_flux_map(
            str(CORE_MAP), dict(fba(ecoli_core).fluxes), background=_probe_svg(tmp_path)
        )


def test_a_background_silences_the_labels_it_already_carries(ecoli_core):
    """Escher's drawing has the labels; drawing CMM's on top turns the text to mush."""

    import numpy as np

    from cmm.core import fba
    from cmm.visualization import escher_flux_map

    fluxes = dict(fba(ecoli_core).fluxes)
    drawing = np.zeros((10, 10, 4), dtype=np.uint8)

    plain = escher_flux_map(str(CORE_MAP), fluxes)
    assert plain.axes[0].texts, "without a background the map labels itself"

    over = escher_flux_map(str(CORE_MAP), fluxes, background=drawing)
    assert not over.axes[0].texts  # only the title, which lives outside axes.texts
    assert not over.axes[0].lines  # nor the metabolite dots the drawing already shows

    # An explicit request still wins over the default.
    forced = escher_flux_map(
        str(CORE_MAP), fluxes, background=drawing, label_reactions=True
    )
    assert forced.axes[0].texts
