from __future__ import annotations

import json

from cmm.features.production import (
    fseof,
    fvseof,
    production_envelope,
    theoretical_yield,
)
from cmm.features.response import flux_response
from cmm.features.sampling import random_flux_sampling
from cmm.visualization import (
    escher_flux_map,
    flux_comparison_figure,
    flux_response_figure,
    fseof_figure,
    fvseof_figure,
    network_flux_map,
    production_envelope_figure,
    sampling_figure,
    save_figure,
    yield_figure,
)

SUCC = "EX_succ_e"
BIOMASS = "Biomass_Ecoli_core"


def _nonblank(path):
    # A 300-DPI multi-element figure is well over 5 KB; a blank canvas is tiny.
    return path.exists() and path.stat().st_size > 5000


def test_production_envelope_figure(ecoli_core, tmp_path):
    envelope = production_envelope(ecoli_core, SUCC, points=15)
    fig = production_envelope_figure(envelope)
    assert len(fig.axes) == 1
    ax = fig.axes[0]
    assert ax.get_xlabel() and ax.get_ylabel()
    assert len(ax.lines) >= 1
    path = save_figure(fig, tmp_path / "envelope.png")
    assert _nonblank(path)


def test_fseof_figure(ecoli_core, tmp_path):
    result = fseof(ecoli_core, SUCC, BIOMASS, n_steps=8, aerobic=False)
    fig = fseof_figure(result, top_n=5)
    ax = fig.axes[0]
    assert 1 <= len(ax.lines) <= 5
    legend_texts = [t.get_text() for t in ax.get_legend().get_texts()]
    # FRD7 (a top amplification target) should be plotted.
    assert any("FRD7" in t for t in legend_texts)
    path = save_figure(fig, tmp_path / "fseof.png")
    assert _nonblank(path)


def test_fvseof_figure(ecoli_core, tmp_path):
    result = fvseof(
        ecoli_core,
        SUCC,
        BIOMASS,
        n_steps=4,
        aerobic=False,
        reactions=["FRD7", "PPC", "MDH", "FUM", "SUCCt3", "EX_succ_e"],
    )
    fig = fvseof_figure(result, top_n=4)
    ax = fig.axes[0]
    assert len(ax.lines) >= 2  # mean + forced lines
    path = save_figure(fig, tmp_path / "fvseof.png")
    assert _nonblank(path)


def test_flux_comparison_figure(tmp_path):
    ref = {"PPC": 1.0, "FRD7": 0.0, "PDH": 9.0}
    eng = {"PPC": 8.0, "FRD7": 12.0, "PDH": 1.0}
    fig = flux_comparison_figure(ref, eng, ["PPC", "FRD7", "PDH"])
    ax = fig.axes[0]
    assert len(ax.patches) == 6  # two groups x three reactions
    path = save_figure(fig, tmp_path / "compare.png", dpi=200)
    assert _nonblank(path)


def test_yield_figure(ecoli_core, tmp_path):
    yields = [
        theoretical_yield(ecoli_core, SUCC, aerobic=True),
        theoretical_yield(ecoli_core, SUCC, aerobic=False),
    ]
    fig = yield_figure(yields)
    ax = fig.axes[0]
    assert len(ax.patches) == 2
    path = save_figure(fig, tmp_path / "yield.png")
    assert _nonblank(path)


def test_flux_response_figure(ecoli_core, tmp_path):
    result = flux_response(ecoli_core, "EX_o2_e", n_steps=12)
    fig = flux_response_figure(result)
    # Response equals the objective here, so there is no twin biomass axis.
    assert len(fig.axes) == 1
    ax = fig.axes[0]
    assert "EX_o2_e" in ax.get_xlabel()
    legend_texts = [t.get_text() for t in ax.get_legend().get_texts()]
    assert any("wild type" in t for t in legend_texts)
    assert any("optimum" in t for t in legend_texts)
    path = save_figure(fig, tmp_path / "response.png")
    assert _nonblank(path)


def test_flux_response_figure_adds_a_biomass_axis_for_a_product(ecoli_core, tmp_path):
    result = flux_response(
        ecoli_core, "PGI", response=SUCC, biomass_fraction=0.3, n_steps=12
    )
    fig = flux_response_figure(result)
    # A product response gets a secondary growth axis anchored at zero.
    assert len(fig.axes) == 2
    twin = fig.axes[1]
    assert "growth" in twin.get_ylabel()
    assert twin.get_ylim()[0] == 0.0
    path = save_figure(fig, tmp_path / "response_product.png")
    assert _nonblank(path)


def test_flux_response_figure_shades_infeasible_range(ecoli_core):
    result = flux_response(
        ecoli_core, "PGI", response=SUCC, biomass_fraction=0.3, n_steps=12
    )
    frame = result.to_frame()
    assert (frame["status"] != "optimal").any()
    fig = flux_response_figure(result)
    # Infeasible stretches are drawn as shaded spans, not per-point lines.
    patches = fig.axes[0].patches
    assert len(patches) >= 1


def test_sampling_figure(ecoli_core, tmp_path):
    from cmm.core import pfba

    result = random_flux_sampling(ecoli_core, n=80, method="achr", thinning=10, seed=2)
    fig = sampling_figure(result, top_n=5, reference=pfba(ecoli_core).fluxes)
    ax = fig.axes[0]
    assert len(ax.get_yticklabels()) == 5
    # Heading on the figure, run parameters on their own line under it: stacked so neither
    # collides with the other or with the axis labels when the panel is short.
    assert fig.get_suptitle() == "Sampled flux distributions"
    params = ax.get_title()
    assert "n = 80" in params
    assert "seed 2" in params
    path = save_figure(fig, tmp_path / "sampling.png")
    assert _nonblank(path)


def test_sampling_figure_accepts_explicit_reactions(ecoli_core):
    result = random_flux_sampling(ecoli_core, n=80, method="achr", thinning=10, seed=2)
    fig = sampling_figure(result, reactions=[BIOMASS, SUCC])
    labels = [t.get_text() for t in fig.axes[0].get_yticklabels()]
    assert labels == [BIOMASS, SUCC]


def test_sampling_figure_rejects_unknown_reactions(ecoli_core):
    import pytest

    result = random_flux_sampling(ecoli_core, n=80, method="achr", thinning=10, seed=2)
    with pytest.raises(KeyError, match="not present"):
        sampling_figure(result, reactions=["NOPE"])


def test_figures_are_high_dpi(ecoli_core):
    envelope = production_envelope(ecoli_core, SUCC, points=10)
    fig = production_envelope_figure(envelope)
    assert fig.dpi >= 300


def test_network_flux_map(ecoli_core, tmp_path):
    from cmm.core import fba

    fig = network_flux_map(ecoli_core, dict(fba(ecoli_core).fluxes), top_n=10)
    assert len(fig.axes) == 1
    path = save_figure(fig, tmp_path / "network.png")
    assert _nonblank(path)


def _mini_escher_map(path):
    """A minimal valid Escher map JSON: two metabolites linked by reaction R1."""

    body = {
        "reactions": {
            "1": {
                "bigg_id": "R1",
                "label_x": 50,
                "label_y": 50,
                "metabolites": [
                    {"bigg_id": "a_c", "coefficient": -1},
                    {"bigg_id": "b_c", "coefficient": 1},
                ],
                "segments": {
                    "s1": {
                        "from_node_id": "10",
                        "to_node_id": "11",
                        "b1": None,
                        "b2": None,
                    }
                },
            }
        },
        "nodes": {
            "10": {
                "node_type": "metabolite",
                "x": 0,
                "y": 0,
                "bigg_id": "a_c",
                "label_x": 0,
                "label_y": 0,
            },
            "11": {
                "node_type": "metabolite",
                "x": 100,
                "y": 100,
                "bigg_id": "b_c",
                "label_x": 100,
                "label_y": 100,
            },
        },
        "text_labels": {},
        "canvas": {"x": -10, "y": -10, "width": 120, "height": 120},
    }
    path.write_text(json.dumps([{"map_name": "mini"}, body]))
    return path


def test_escher_flux_map_renders_layout(tmp_path):
    map_path = _mini_escher_map(tmp_path / "mini.json")
    fig = escher_flux_map(str(map_path), {"R1": 7.5}, title="mini")
    assert len(fig.axes) >= 1
    # The reaction edge (a PathPatch) is drawn coloured by flux.
    from matplotlib.patches import PathPatch

    assert any(isinstance(p, PathPatch) for p in fig.axes[0].patches)
    path = save_figure(fig, tmp_path / "escher.png")
    assert _nonblank(path)


def test_escher_flux_map_handles_missing_flux(tmp_path):
    map_path = _mini_escher_map(tmp_path / "mini.json")
    # No flux for R1 -> the reaction is drawn in the zero-flux style, no error.
    fig = escher_flux_map(str(map_path), {})
    assert len(fig.axes) >= 1
