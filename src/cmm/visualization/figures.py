"""Publication-quality figures for production-design results.

Figures are built with `matplotlib.figure.Figure` directly (no pyplot global state), so the
module is import-safe and headless. Styling targets print/paper use: configurable single- or
double-column font presets, clean spines, a colour-blind-safe palette, and readable type.

`network_flux_map` is a schematic flux network (force-directed, edge width ∝ |flux|). It is
intentionally NOT a curated Escher map: hand-laid Escher metabolic maps encode biochemical
layout that this lightweight, dependency-free renderer does not attempt. Use it for a quick
flux overview; use a curated Escher map for a publication network figure.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from matplotlib import cm, colormaps
from matplotlib.colors import TwoSlopeNorm
from matplotlib.figure import Figure
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath

from cmm.features.production import (
    FseofResult,
    FvseofResult,
    ProductionEnvelope,
    ProductionYield,
)

# Okabe-Ito colour-blind-safe palette.
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#999999"]
_LINESTYLES = ["-", "--", "-.", ":"]
_FONT_PRESETS = {
    1: {"title": 10, "label": 9, "tick": 8, "legend": 8},   # single column (~3.3 in)
    2: {"title": 13, "label": 12, "tick": 10, "legend": 9},  # double column (~6.5 in)
}


def _new_figure(width: float = 6.0, height: float = 4.0, column_width: int = 2):
    fig = Figure(figsize=(width, height), dpi=300)
    fig.set_facecolor("white")
    ax = fig.subplots()
    return fig, ax, _FONT_PRESETS.get(column_width, _FONT_PRESETS[2])


def _style(ax, font, *, xlabel: str, ylabel: str, title: str) -> None:
    ax.set_xlabel(xlabel, fontsize=font["label"])
    ax.set_ylabel(ylabel, fontsize=font["label"])
    ax.set_title(title, fontsize=font["title"], fontweight="bold", pad=10)
    ax.tick_params(labelsize=font["tick"])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.set_axisbelow(True)


def save_figure(fig: Figure, path: str | Path, dpi: int = 300) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    return path


def production_envelope_figure(
    envelope: ProductionEnvelope,
    *,
    product_label: str | None = None,
    title: str = "Production envelope",
    column_width: int = 2,
) -> Figure:
    """Growth-vs-product phenotypic phase plane with the feasible region shaded."""

    fig, ax, font = _new_figure(column_width=column_width)
    frame = envelope.to_frame().sort_values("product_flux")
    x = frame["product_flux"].to_numpy()
    gmax = frame["growth_max"].to_numpy()
    gmin = frame["growth_min"].to_numpy()
    ax.fill_between(x, gmin, gmax, color=PALETTE[0], alpha=0.18, label="feasible region")
    ax.plot(x, gmax, color=PALETTE[0], linewidth=2.0, label="max growth")
    ax.plot(x, gmin, color=PALETTE[0], linewidth=1.0, linestyle="--", alpha=0.7,
            label="min growth")
    label = product_label or envelope.product
    suffix = "  (min growth = 0 across range)" if float(np.max(gmin)) <= 1e-9 else ""
    _style(ax, font,
           xlabel=f"{label} flux (mmol gDW$^{{-1}}$ h$^{{-1}}$)",
           ylabel="growth rate (h$^{-1}$)",
           title=title + suffix)
    ax.legend(fontsize=font["legend"], frameon=False)
    fig.tight_layout()
    return fig


def fseof_figure(
    result: FseofResult,
    *,
    top_n: int = 6,
    title: str = "FSEOF amplification targets",
    column_width: int = 2,
) -> Figure:
    """Flux of the top amplification-target reactions vs enforced product flux."""

    fig, ax, font = _new_figure(width=7.5, height=4.2, column_width=column_width)
    levels = list(result.enforced_levels)
    targets = result.amplification_targets()
    ranked = sorted(
        targets,
        key=lambda rid: abs(result.trends.loc[rid, levels[-1]]) - abs(result.trends.loc[rid, levels[0]]),
        reverse=True,
    )[:top_n]
    for i, rid in enumerate(ranked):
        ys = [result.trends.loc[rid, lv] for lv in levels]
        # Vary colour AND linestyle/marker so near-coincident lines stay distinguishable.
        ax.plot(levels, ys,
                marker=["o", "s", "^", "D", "v", "P"][i % 6], markersize=6,
                linewidth=2.4 if i < 2 else 1.6,
                linestyle=_LINESTYLES[i % len(_LINESTYLES)],
                color=PALETTE[i % len(PALETTE)], label=rid)
    _style(ax, font,
           xlabel=f"enforced {result.product} flux (mmol gDW$^{{-1}}$ h$^{{-1}}$)",
           ylabel="reaction flux (mmol gDW$^{-1}$ h$^{-1}$)",
           title=title)
    ax.legend(fontsize=font["legend"], frameon=False,
              bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.tight_layout()
    return fig


def fvseof_figure(
    result: FvseofResult,
    *,
    top_n: int = 5,
    title: str = "FVSEOF robust amplification targets",
    column_width: int = 2,
) -> Figure:
    """Mean flux (solid) and forced-minimum |flux| (dashed) of robust targets vs product.

    A target whose forced-minimum rises with enforced product is *robustly* amplified — the
    reaction cannot avoid carrying more flux. That dashed line rising is the FVSEOF signal.
    """

    fig, ax, font = _new_figure(width=7.5, height=4.2, column_width=column_width)
    levels = list(result.enforced_levels)
    level_cols = list(result.mean.columns)
    robust = result.robust_targets() or result.amplification_targets()
    ranked = sorted(
        robust,
        key=lambda rid: abs(result.mean.loc[rid, level_cols[-1]]) - abs(result.mean.loc[rid, level_cols[0]]),
        reverse=True,
    )[:top_n]
    for i, rid in enumerate(ranked):
        colour = PALETTE[i % len(PALETTE)]
        mean = [result.mean.loc[rid, c] for c in level_cols]
        forced = [result.forced.loc[rid, c] for c in level_cols]
        ax.plot(levels, mean, marker="o", markersize=5, linewidth=2.0, color=colour,
                label=f"{rid} (mean)")
        ax.plot(levels, forced, linestyle="--", linewidth=1.4, color=colour, alpha=0.8)
    _style(ax, font,
           xlabel=f"enforced {result.product} flux (mmol gDW$^{{-1}}$ h$^{{-1}}$)",
           ylabel="flux (mmol gDW$^{-1}$ h$^{-1}$)",
           title=title)
    ax.legend(fontsize=font["legend"], frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.text(0.0, -0.16, "solid = mean flux, dashed = forced minimum |flux|",
            transform=ax.transAxes, fontsize=font["legend"], color="#555555")
    fig.tight_layout()
    return fig


def flux_comparison_figure(
    reference: dict[str, float],
    comparison: dict[str, float],
    reactions: list[str],
    *,
    reference_label: str = "wild type",
    comparison_label: str = "engineered",
    title: str = "Flux comparison",
    column_width: int = 2,
) -> Figure:
    """Grouped bar chart comparing two flux distributions over selected reactions."""

    fig, ax, font = _new_figure(width=max(6.0, 0.7 * len(reactions) + 2), height=4.0,
                                column_width=column_width)
    idx = np.arange(len(reactions))
    width = 0.38
    ref = [reference.get(r, 0.0) for r in reactions]
    cmp = [comparison.get(r, 0.0) for r in reactions]
    ax.bar(idx - width / 2, ref, width, color=PALETTE[6], label=reference_label)
    ax.bar(idx + width / 2, cmp, width, color=PALETTE[1], label=comparison_label)
    ax.set_xticks(idx)
    ax.set_xticklabels(reactions, rotation=45, ha="right", fontsize=font["tick"])
    ax.axhline(0, color="black", linewidth=0.8)
    _style(ax, font, xlabel="reaction", ylabel="flux (mmol gDW$^{-1}$ h$^{-1}$)", title=title)
    ax.legend(fontsize=font["legend"], frameon=False)
    fig.tight_layout()
    return fig


def yield_figure(
    yields: list[ProductionYield],
    *,
    substrate_label: str | None = None,
    title: str = "Theoretical molar yield",
    column_width: int = 2,
) -> Figure:
    """Bar chart of theoretical molar yields (e.g. aerobic vs anaerobic).

    A dashed line marks the substrate carbon ceiling; bars above it are annotated with
    ``+CO2`` to flag that the yield is achievable only with net CO2 fixation.
    """

    fig, ax, font = _new_figure(width=5.2, height=4.2, column_width=column_width)
    labels = [f"{'aerobic' if y.aerobic else 'anaerobic'}" for y in yields]
    values = [y.molar_yield for y in yields]
    colors = [PALETTE[0] if y.aerobic else PALETTE[2] for y in yields]
    bars = ax.bar(labels, values, color=colors, width=0.55)
    for bar, y in zip(bars, yields, strict=True):
        text = f"{y.molar_yield:.2f}" + (" +CO$_2$" if y.exceeds_carbon_ceiling else "")
        ax.text(bar.get_x() + bar.get_width() / 2, y.molar_yield, text,
                ha="center", va="bottom", fontsize=font["tick"])
    ceiling = next((y.carbon_ceiling for y in yields if y.carbon_ceiling is not None), None)
    if ceiling is not None:
        ax.axhline(ceiling, color="#999999", linestyle="--", linewidth=1.0)
        ax.text(ax.get_xlim()[1], ceiling, f" carbon ceiling {ceiling:.2f}",
                va="bottom", ha="right", fontsize=font["legend"], color="#555555")
    substrate = substrate_label or (yields[0].substrate if yields else "substrate")
    product = yields[0].product if yields else "product"
    _style(ax, font, xlabel="condition",
           ylabel=f"mol {product} / mol {substrate}", title=title)
    fig.tight_layout()
    return fig


# Currency / cofactor metabolites that connect to nearly everything and turn a metabolic
# network drawing into a hairball. Excluded so the map shows the carbon backbone.
_CURRENCY = frozenset({
    "h2o", "h", "atp", "adp", "amp", "pi", "ppi", "co2", "o2", "nh4",
    "nad", "nadh", "nadp", "nadph", "fad", "fadh2", "coa", "q8", "q8h2",
    "so4", "h2", "pyr_h",
})


def _base_id(mid: str) -> str:
    return mid.rsplit("_", 1)[0]


def flux_log_change_figure(
    log_changes: dict[str, float],
    *,
    top_n: int = 20,
    title: str = "Flux log2 fold-change",
    source_label: str = "A",
    target_label: str = "B",
    column_width: int = 2,
) -> Figure:
    """Horizontal diverging bar chart of the largest flux log2 fold-changes between conditions."""

    ranked = sorted(log_changes.items(), key=lambda kv: -abs(kv[1]))[:top_n]
    ranked = ranked[::-1]  # largest at the top after barh
    labels = [rid for rid, _ in ranked]
    values = [v for _, v in ranked]

    fig, ax, font = _new_figure(width=6.5, height=max(3.5, 0.32 * len(ranked) + 1),
                                column_width=column_width)
    colors = [PALETTE[1] if v > 0 else PALETTE[0] for v in values]
    ax.barh(range(len(values)), values, color=colors)
    ax.set_yticks(range(len(values)))
    ax.set_yticklabels(labels, fontsize=font["tick"])
    ax.axvline(0, color="black", linewidth=0.8)
    _style(ax, font,
           xlabel=f"log$_2$( |flux {target_label}| / |flux {source_label}| )",
           ylabel="reaction", title=title)
    ax.grid(True, axis="y", alpha=0.0)
    fig.tight_layout()
    return fig


def network_flux_map(
    model,
    fluxes: dict[str, float],
    *,
    top_n: int = 12,
    title: str = "Flux network (top reactions)",
    seed: int = 0,
) -> Figure:
    """Schematic carbon-backbone flux network: top-|flux| reactions as edges, width ∝ |flux|.

    Currency metabolites (ATP, H2O, CO2, NAD(P)H, ...) are excluded so the layout follows the
    carbon skeleton. A deterministic force-directed layout; this is a quick flux overview, not
    a curated Escher map (see the module docstring).
    """

    fig, ax, font = _new_figure(width=8.0, height=4.8)
    ax.set_axis_off()
    ax.set_title(title, fontsize=font["title"], fontweight="bold")

    met_ids: list[str] = []
    edges: list[tuple[int, int, float, str]] = []
    index: dict[str, int] = {}

    def _node(mid: str) -> int:
        if mid not in index:
            index[mid] = len(met_ids)
            met_ids.append(mid)
        return index[mid]

    def _carbon_mets(mets):
        return [m for m in mets if _base_id(m.id) not in _CURRENCY]

    ranked = sorted(
        ((rid, flux) for rid, flux in fluxes.items() if abs(flux) > 1e-6),
        key=lambda kv: -abs(kv[1]),
    )
    for rid, flux in ranked:
        rxn = model.reactions.get_by_id(rid)
        if rxn.boundary:
            continue
        reactants = _carbon_mets([m for m, c in rxn.metabolites.items() if c < 0])
        products = _carbon_mets([m for m, c in rxn.metabolites.items() if c > 0])
        if not reactants or not products:
            continue  # pure cofactor reaction — nothing on the carbon backbone to draw
        if flux < 0:
            reactants, products = products, reactants
        edges.append((_node(reactants[0].id), _node(products[0].id), abs(flux), rid))
        if len(edges) >= top_n:
            break

    n = len(met_ids)
    if n == 0:
        ax.text(0.5, 0.5, "no internal carbon flux to display", ha="center", va="center")
        return fig

    _ = seed  # retained for API compatibility; the row layout is deterministic
    pos = _component_row_layout(n, [(s, d) for s, d, _, _ in edges])
    max_w = max((w for _, _, w, _ in edges), default=1.0)
    cmap = colormaps["viridis"]
    reaction_font = max(6, font["legend"] - 2)
    for s, d, w, rid in edges:
        x0, y0 = pos[s]
        x1, y1 = pos[d]
        ax.annotate(
            "", xy=(x1, y1), xytext=(x0, y0),
            arrowprops=dict(arrowstyle="-|>", color=cmap(0.15 + 0.85 * w / max_w),
                            lw=1.2 + 5.0 * w / max_w, alpha=0.9, shrinkA=14, shrinkB=14),
        )
        ax.text((x0 + x1) / 2, (y0 + y1) / 2 - 0.055, rid, fontsize=reaction_font,
                ha="center", va="center", color="#222222",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="#cccccc", alpha=0.85))
    ax.scatter(pos[:, 0], pos[:, 1], s=130, color=PALETTE[0], zorder=3, edgecolors="white")
    metabolite_font = max(7, font["tick"] - 2)
    for i, (mid, (x, y)) in enumerate(zip(met_ids, pos, strict=True)):
        above = i % 2 == 0
        y_offset = 0.055 if above else -0.09
        ax.text(
            x, y + y_offset, _display_metabolite_id(mid),
            fontsize=metabolite_font, fontweight="bold",
            ha="center", va="bottom" if above else "top", color="#11243b",
        )
    # Colour-bar proxy for flux magnitude.
    ax.text(0.0, -0.08, f"edge width / colour ∝ |flux|  (max {max_w:.1f} mmol gDW$^{{-1}}$ h$^{{-1}}$)",
            transform=ax.transAxes, fontsize=font["legend"], color="#555555")
    ax.set_xlim(-0.12, 1.12)
    ax.set_ylim(-0.18, 1.06)
    return fig


def _component_row_layout(n: int, edges: list[tuple[int, int]]) -> np.ndarray:
    """Pack disconnected carbon-backbone components into compact horizontal rows."""

    adjacency = {i: set() for i in range(n)}
    for source, dest in edges:
        adjacency[source].add(dest)
        adjacency[dest].add(source)

    seen: set[int] = set()
    components: list[list[int]] = []
    for node in range(n):
        if node in seen:
            continue
        stack = [node]
        seen.add(node)
        component: list[int] = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbour in adjacency[current]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        components.append(sorted(component))

    components.sort(key=lambda nodes: (-len(nodes), nodes[0]))
    y_values = np.linspace(0.82, 0.22, len(components)) if len(components) > 1 else [0.55]
    pos = np.zeros((n, 2))
    for component, y in zip(components, y_values, strict=True):
        ordered = _ordered_component_nodes(component, edges)
        if len(ordered) == 1:
            xs = [0.5]
        else:
            xs = np.linspace(0.08, 0.92, len(ordered))
        for node, x in zip(ordered, xs, strict=True):
            pos[node] = (float(x), float(y))
    return pos


def _ordered_component_nodes(component: list[int], edges: list[tuple[int, int]]) -> list[int]:
    """Order one component from likely upstream sources to downstream products."""

    nodes = set(component)
    outgoing = {node: [] for node in component}
    indegree = {node: 0 for node in component}
    for source, dest in edges:
        if source in nodes and dest in nodes:
            outgoing[source].append(dest)
            indegree[dest] += 1

    sources = [node for node in component if indegree[node] == 0 and outgoing[node]]
    if not sources:
        sources = [component[0]]

    depths: dict[int, int] = {}
    queue = list(sorted(sources))
    for source in queue:
        depths[source] = 0
    while queue:
        current = queue.pop(0)
        for dest in sorted(outgoing[current]):
            if dest not in depths:
                depths[dest] = depths[current] + 1
                queue.append(dest)

    fallback = max(depths.values(), default=0) + 1
    return sorted(component, key=lambda node: (depths.get(node, fallback), node))


def _display_metabolite_id(mid: str) -> str:
    return _base_id(mid).replace("__", "-")


def escher_flux_map(
    map_path: str | Path,
    fluxes: dict[str, float],
    *,
    title: str | None = None,
    abs_max: float | None = None,
    label_metabolites: bool = True,
    label_reactions: bool = True,
    width: float = 12.0,
) -> Figure:
    """Render an Escher map (curated node/segment layout) coloured by flux.

    This reuses a standard Escher map JSON's hand-laid coordinates and bezier segments — the
    same layout a curated Escher map uses — so the figure matches that biochemical layout
    instead of an invented force layout. Reaction edge width and colour encode flux
    (diverging: blue = negative/reverse, red = positive/forward). ``map_path`` is supplied by
    the caller (CMM bundles no maps).
    """

    with open(map_path) as handle:
        data = json.load(handle)
    body = data[1] if isinstance(data, list) else data
    nodes = body["nodes"]
    reactions = body["reactions"]
    node_xy = {nid: (float(n["x"]), float(n["y"])) for nid, n in nodes.items()}

    fig = Figure(figsize=(width, width * 0.82), dpi=200)
    fig.set_facecolor("white")
    ax = fig.subplots()
    ax.set_axis_off()
    ax.set_aspect("equal")

    signed = [fluxes.get(r["bigg_id"], 0.0) for r in reactions.values()]
    amax = abs_max or max((abs(v) for v in signed), default=1.0) or 1.0
    norm = TwoSlopeNorm(vmin=-amax, vcenter=0.0, vmax=amax)
    cmap = colormaps["coolwarm"]

    for r in reactions.values():
        flux = fluxes.get(r["bigg_id"], 0.0)
        mag = abs(flux)
        if mag <= 1e-9:
            color, lw, alpha = "#d7dbe0", 0.7, 0.7
        else:
            color, lw, alpha = cmap(norm(flux)), 0.8 + 5.5 * mag / amax, 0.95
        for seg in r["segments"].values():
            a = node_xy.get(seg["from_node_id"])
            b = node_xy.get(seg["to_node_id"])
            if a is None or b is None:
                continue
            b1, b2 = seg.get("b1"), seg.get("b2")
            if b1 and b2:
                verts = [a, (b1["x"], b1["y"]), (b2["x"], b2["y"]), b]
                codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
            else:
                verts = [a, b]
                codes = [MplPath.MOVETO, MplPath.LINETO]
            ax.add_patch(PathPatch(MplPath(verts, codes), fill=False, edgecolor=color,
                                   lw=lw, alpha=alpha, capstyle="round", joinstyle="round"))

    xs, ys = [], []
    for n in nodes.values():
        if n.get("node_type") == "metabolite":
            ax.plot(n["x"], n["y"], "o", ms=2.2, color="#2b3a4a", zorder=3)
            xs.append(n["x"])
            ys.append(n["y"])
            if label_metabolites:
                ax.text(n.get("label_x", n["x"]), n.get("label_y", n["y"]),
                        str(n.get("bigg_id", "")).rsplit("_", 1)[0], fontsize=4.5,
                        ha="left", va="center", color="#11243b", zorder=4)
    if label_reactions:
        for r in reactions.values():
            ax.text(r.get("label_x", 0), r.get("label_y", 0), r["bigg_id"], fontsize=4.5,
                    ha="left", va="center", color="#7a3b00", zorder=4, fontstyle="italic")

    if xs:
        mx = 0.04 * (max(xs) - min(xs) + 1)
        my = 0.04 * (max(ys) - min(ys) + 1)
        ax.set_xlim(min(xs) - mx, max(xs) + mx)
        ax.set_ylim(max(ys) + my, min(ys) - my)  # Escher y grows downward
    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(mappable, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label("flux (mmol gDW$^{-1}$ h$^{-1}$)", fontsize=9)
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    fig.tight_layout()
    return fig
