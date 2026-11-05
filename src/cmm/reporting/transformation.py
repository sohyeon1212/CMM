"""Report rendering for transformation-target runs.

Figures come from a checked-in R script and the page is assembled here, which is exactly how
the production report is built. The two workflows draw different panels, so they have separate
R scripts, but a reader gets the same artifact contract from either: a linked page, a
standalone page carrying its figures, a 300-dpi raster plus editable vector per figure, and a
``report_validation.json`` recording whether the run passed its completion gate.

The page is deliberately opinionated about what it must say. A transformation ranking is easy
to over-read: it depends on how the reference state was built (and what still differs from the
papers), on an epsilon that is either derived from the reference samples or chosen outright,
and on a candidate set whose size is the denominator of any percentile claim. Those three facts are
rendered as prose in the report body, not left in a provenance file for a reader to find.
"""

from __future__ import annotations

import base64
from dataclasses import dataclass
import hashlib
import html
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

from cmm.reporting.schema import ValidationReport
from cmm.reporting._rscript import (
    FigureManifest,
    invoke_renderer,
    load_figure_manifest,
    relative_path,
)


@dataclass(frozen=True)
class TransformationReport:
    """Where the rendered report and its figures were written."""

    run_directory: Path
    #: References figures by relative path. Small, and the 300-dpi originals stay beside it.
    report_html: Path
    #: Every figure inlined as a data URI. This is the copy to send someone: the relative-path
    #: copy loses all of its figures the moment it leaves the run directory, and says nothing.
    report_standalone_html: Path
    #: What the R renderer produced, and what it declined to draw with a stated reason.
    figure_manifest: Path
    #: The completion gate's verdict, written beside the report so it travels with the run.
    report_validation: Path
    figures: tuple[Path, ...]


class TransformationReportError(RuntimeError):
    """The transformation report could not be rendered from the run bundle."""


_STYLE = """
body{font:15px/1.6 -apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;
 color:#1a2433;max-width:960px;margin:0 auto;padding:32px 24px 72px}
h1{font-size:26px;margin:0 0 4px} h2{font-size:19px;margin:36px 0 10px;
 padding-bottom:6px;border-bottom:1px solid #e2e8f0}
.sub{color:#5a6b80;margin:0 0 28px}
table{border-collapse:collapse;width:100%;margin:12px 0;font-size:14px}
th,td{border-bottom:1px solid #e2e8f0;padding:7px 10px;text-align:left}
th{background:#f7fafc;font-weight:600}
tr.mark td{background:#f0fff4}
figure{margin:18px 0}img{max-width:100%;height:auto;border:1px solid #e2e8f0;border-radius:4px}
figcaption{color:#5a6b80;font-size:13px;margin-top:6px}
.note{background:#fffaf0;border-left:3px solid #dd6b20;padding:12px 16px;margin:16px 0}
.note strong{color:#9c4221}
code{background:#f7fafc;padding:1px 5px;border-radius:3px;font-size:13px}
.meta{color:#718096;font-size:13px}
"""


def _read_json(path: Path) -> Mapping[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_csv(path: Path) -> list[dict[str, Any]]:
    import csv

    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _number(value: object, digits: int = 4) -> str:
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return html.escape(str(value))
    if number != number:
        return "—"
    return f"{number:.{digits}g}"


def _table(
    rows: Sequence[Mapping[str, object]],
    columns: Sequence[str],
    mark: str | None = None,
) -> str:
    if not rows:
        return "<p class='meta'>No rows.</p>"
    head = "".join(f"<th>{html.escape(c.replace('_', ' '))}</th>" for c in columns)
    body = []
    for row in rows:
        cells = "".join(f"<td>{_number(row.get(c))}</td>" for c in columns)
        highlighted = mark is not None and str(row.get("target_id")) == mark
        body.append(
            f"<tr class='mark'>{cells}</tr>" if highlighted else f"<tr>{cells}</tr>"
        )
    return (
        f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"
    )


_FIGURE_ORDER = (
    "fig01_component_scores",
    "fig02_transformation_ranking",
    "fig03_ranking_vs_moma",
    "fig04_epsilon_sensitivity",
)
#: Only the ranking is unconditional. The component-score panel needs the best/worst contrast
#: that only rMTA computes; the rank comparison needs rMTA's score decomposition or the optional
#: MOMA baseline; the epsilon sweep is a configurable stage. A run that switched one off, or
#: ran MTA, is a valid run whose figure is legitimately absent — the manifest records why,
#: which is what keeps "absent" distinguishable from "broken". Every id here is still expected
#: *in* the manifest: an id the renderer never recorded is a broken renderer, not a choice.
_REQUIRED_FIGURES = frozenset({"fig02_transformation_ranking"})


def renderer_script_path() -> Path:
    """Return the checked-in R renderer used by :func:`render_transformation_figures`."""

    return Path(__file__).with_name("render_transformation_figures.R")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def render_transformation_figures(
    run_dir: str | Path,
    *,
    highlight: str | None = None,
    rscript: str | Path = "Rscript",
    renderer: str | Path | None = None,
    output_dir: str = "figures",
) -> FigureManifest:
    """Render the run's artwork with checked-in R code and return its manifest."""

    root = _validated_root(run_dir)
    script = (
        Path(renderer) if renderer is not None else renderer_script_path()
    ).expanduser()
    script = script.resolve()
    figures_dir = relative_path(
        root,
        output_dir,
        label="figure output directory",
        error_type=TransformationReportError,
    )
    figures_dir.mkdir(parents=True, exist_ok=True)
    figure_manifest_path = figures_dir / "figure_manifest.json"
    invoke_renderer(
        script=script,
        script_sha256=_sha256(script) if script.is_file() else "unknown",
        run_dir=root,
        manifest_path=root / "00_manifest.json",
        figures_dir=figures_dir,
        figure_manifest_path=figure_manifest_path,
        rscript=rscript,
        error_type=TransformationReportError,
        label="R transformation renderer",
        # Which candidate to mark is a reading choice, so it reaches R the same way the font
        # does rather than being written into the run bundle.
        extra_environment={"CMM_TRANSFORMATION_HIGHLIGHT": highlight or ""},
    )
    return load_figure_manifest(
        figure_manifest_path,
        root,
        required=_REQUIRED_FIGURES,
        order=_FIGURE_ORDER,
        error_type=TransformationReportError,
    )


def _validated_root(run_dir: str | Path) -> Path:
    """Resolve a run directory, refusing anything that is not a transformation bundle."""

    root = Path(run_dir).expanduser().resolve()
    manifest_path = root / "00_manifest.json"
    if not manifest_path.is_file():
        raise TransformationReportError(f"no 00_manifest.json in {root}")
    manifest = _read_json(manifest_path)
    if manifest.get("workflow") != "transformation_target_discovery":
        raise TransformationReportError(
            f"{root} is a {manifest.get('workflow')!r} run, not a transformation run"
        )
    return root


def _load_bundle(root: Path) -> dict[str, Any]:
    """Read every table and sidecar the page is built from, in one place."""

    manifest = _read_json(root / "00_manifest.json")
    artifacts = manifest["artifacts"]

    def path_for(role: str) -> Path | None:
        entry = artifacts.get(role)
        return root / entry["path"] if entry else None

    def required_path(role: str) -> Path:
        """A role the page cannot be written without.

        Naming the missing role beats failing later on a ``None``, because the usual cause is
        an interrupted run whose manifest was written before the artifact was.
        """

        path = path_for(role)
        if path is None:
            raise TransformationReportError(
                f"{root}/00_manifest.json has no {role!r} artifact"
            )
        return path

    ranking = _read_csv(required_path("transformation_ranking"))
    if not ranking:
        raise TransformationReportError("the run contains an empty ranking")
    for row in ranking:
        row["score"] = float(row["score"])
        row["rank"] = int(row["rank"])

    baseline_path = path_for("moma_baseline")
    baseline = (
        _read_csv(baseline_path)
        if baseline_path and baseline_path.stat().st_size
        else []
    )
    for row in baseline:
        row["rank"] = int(row["rank"])
        row["moma_score"] = float(row["moma_score"])
    sweep_path = path_for("epsilon_sensitivity")
    sweep = _read_csv(sweep_path) if sweep_path and sweep_path.stat().st_size else []
    for row in sweep:
        row["epsilon"] = float(row["epsilon"])
        row["rank"] = int(row["rank"])
    decomposition_path = path_for("score_decomposition")
    decomposition = (
        _read_csv(decomposition_path)
        if decomposition_path
        and (artifacts.get("score_decomposition") or {}).get("status") == "complete"
        else []
    )
    for row in decomposition:
        for key in ("bTS", "mTS", "wTS", "rTS"):
            row[key] = float(row[key])
        for key in ("rank_rTS", "rank_mTS", "rank_delta"):
            row[key] = int(row[key])
        row["amplified"] = str(row["amplified"]).strip().lower() == "true"
    preflight_path = path_for("preflight")
    preflight = (
        _read_csv(preflight_path) if preflight_path and preflight_path.is_file() else []
    )

    return {
        "manifest": manifest,
        "artifacts": artifacts,
        "summary": _read_json(required_path("summary")),
        "provenance": _read_json(required_path("provenance")),
        "ranking": ranking,
        "baseline": baseline,
        "sweep": sweep,
        "decomposition": decomposition,
        "preflight": preflight,
    }


def render_transformation_report(
    run_dir: str | Path,
    *,
    highlight: str | None = None,
    top_n: int = 15,
    rscript: str | Path = "Rscript",
    renderer: str | Path | None = None,
) -> TransformationReport:
    """Render figures through R, then build both report variants and the validation record.

    ``highlight`` names a candidate to mark throughout — the knockout under test in a
    validation run, or the one a user came to check.
    """

    root = _validated_root(run_dir)
    figure_manifest = render_transformation_figures(
        root, highlight=highlight, rscript=rscript, renderer=renderer
    )
    bundle = _load_bundle(root)

    html_text = _compose(
        root=root,
        summary=bundle["summary"],
        provenance=bundle["provenance"],
        artifacts=bundle["artifacts"],
        ranking=bundle["ranking"],
        baseline=bundle["baseline"],
        sweep=bundle["sweep"],
        decomposition=bundle["decomposition"],
        preflight=bundle["preflight"],
        figure_manifest=figure_manifest,
        highlight=highlight,
        top_n=top_n,
        standalone=False,
    )
    report_path = root / "report.html"
    report_path.write_text(html_text, encoding="utf-8")

    standalone_text = _compose(
        root=root,
        summary=bundle["summary"],
        provenance=bundle["provenance"],
        artifacts=bundle["artifacts"],
        ranking=bundle["ranking"],
        baseline=bundle["baseline"],
        sweep=bundle["sweep"],
        decomposition=bundle["decomposition"],
        preflight=bundle["preflight"],
        figure_manifest=figure_manifest,
        highlight=highlight,
        top_n=top_n,
        standalone=True,
    )
    standalone_path = root / "report_standalone.html"
    standalone_path.write_text(_assert_embedded(standalone_text), encoding="utf-8")

    rendered = tuple(
        root / str(figure["outputs"]["png"])
        for figure in figure_manifest.figures
        if figure.get("status") == "rendered"
    )
    validation_path = _write_report_validation(root, figure_manifest)
    return TransformationReport(
        run_directory=root,
        report_html=report_path,
        report_standalone_html=standalone_path,
        figure_manifest=figure_manifest.path,
        report_validation=validation_path,
        figures=rendered,
    )


def _write_report_validation(root: Path, figure_manifest: FigureManifest) -> Path:
    """Record what was rendered, so a later validator can detect a stale report."""

    rendered = [
        figure["id"]
        for figure in figure_manifest.figures
        if figure.get("status") == "rendered"
    ]
    unavailable = {
        str(figure["id"]): str(figure.get("reason", ""))
        for figure in figure_manifest.figures
        if figure.get("status") != "rendered"
    }
    payload = {
        "schema_version": 2,
        "workflow": "transformation_target_discovery",
        "render": "succeeded",
        "figures_rendered": rendered,
        "figures_unavailable": unavailable,
        "renderer": dict(figure_manifest.renderer),
    }
    path = root / "report_validation.json"
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return path


def _assert_embedded(html_text: str) -> str:
    """Refuse to write a standalone page that would render with a figure missing.

    This is the failure the reporting contract singles out as invisible: a relative ``src``
    resolves against wherever the copy landed, so the browser draws blank space and reports
    nothing. Checking the image sources — not the whole document — is the point, because the
    page deliberately prints artifact paths as text and those must stay readable.
    """

    for source in re.findall(r"<img[^>]*\ssrc='([^']*)'", html_text):
        if not source.startswith("data:image/png;base64,"):
            raise TransformationReportError(
                f"standalone report references an image it does not carry: {source[:60]}"
            )
    return html_text


def _reference_text(provenance: Mapping[str, Any]) -> str:
    """How v_ref was built, in one line, from the parameters the run recorded."""

    if provenance.get("reference_method") != "imat":
        return str(provenance.get("reference_method"))
    return (
        f"iMAT (delta {provenance.get('imat_delta')}, beta {provenance.get('imat_beta')}, "
        f"epsilon {provenance.get('imat_epsilon')}, GPR {provenance.get('imat_gpr_rule')}); "
        f"mean of {provenance.get('reference_samples')} OptGP samples, thinning "
        f"{provenance.get('reference_thinning')}, seed {provenance.get('reference_seed')}"
    )


def _medium_text(provenance: Mapping[str, Any], coverage: Mapping[str, Any]) -> str:
    """The medium as a reader must read it: a preset name, or what "as loaded" leaves open."""

    open_count = coverage.get("n_exchanges_open_for_uptake")
    total = coverage.get("n_exchanges")
    openness = (
        f"{open_count} of {total} exchanges permit uptake"
        if open_count is not None and total is not None
        else None
    )
    name = provenance.get("medium")
    if name:
        return f"{name} ({openness})" if openness else str(name)
    return f"model as loaded — {openness}" if openness else "model as loaded"


def _coverage_text(coverage: Mapping[str, Any]) -> str:
    measured = coverage.get("genes_measured_in_model")
    total = coverage.get("model_genes")
    if measured is None or not total:
        return "not recorded"
    return f"{measured} of {total} ({100.0 * float(measured) / float(total):.1f}%)"


def _masked_text(structure: Mapping[str, Any]) -> str:
    text = str(structure.get("n_impossible_masked"))
    reversed_count = structure.get("n_impossible_masked_reversed")
    if reversed_count is not None:
        text += f" (worst-case map: {reversed_count})"
    return text


def _applied_medium_text(applied: Mapping[str, Any]) -> str:
    uptake = applied.get("applied") or {}
    listed = ", ".join(
        f"{exchange} {float(rate):g}" for exchange, rate in sorted(uptake.items())
    )
    dropped = applied.get("dropped") or []
    text = f"{applied.get('medium')}: {listed or 'no exchanges'}"
    if dropped:
        text += f"; not in the model: {', '.join(map(str, dropped))}"
    return text


def _compose(
    *,
    root: Path,
    summary: Mapping[str, Any],
    provenance: Mapping[str, Any],
    artifacts: Mapping[str, Any],
    ranking: Sequence[Mapping[str, Any]],
    baseline: Sequence[Mapping[str, Any]],
    sweep: Sequence[Mapping[str, Any]],
    figure_manifest: FigureManifest,
    highlight: str | None,
    decomposition: Sequence[Mapping[str, Any]] = (),
    preflight: Sequence[Mapping[str, Any]] = (),
    top_n: int,
    standalone: bool,
) -> str:
    """Build one variant of the page.

    The two variants differ only in how they reference things that live outside the file:
    ``standalone`` embeds each figure as a data URI and renders artifact paths as plain
    filenames, because a dead link is worse than no link once the page has been moved.
    """

    e = html.escape
    method = str(summary.get("method", "?")).upper()
    construction = summary.get("candidate_construction") or {}
    direction = summary.get("direction_construction") or {}
    structure = summary.get("ranking_structure") or {}
    coverage = summary.get("data_coverage") or {}
    amplified_rows = sorted(
        (row for row in decomposition if row["amplified"]),
        key=lambda row: int(row["rank_rTS"]),
    )
    decomposed = {str(row["target_id"]): row for row in decomposition}
    by_id = {str(figure["id"]): figure for figure in figure_manifest.figures}

    def panel(figure_id: str) -> str:
        """Render one figure, or say plainly why it is not there."""

        figure = by_id.get(figure_id)
        if figure is None:
            # The one path that used to fail silently: an id the renderer never recorded —
            # it errored before recording, or the id was renamed — left no figure and no note.
            return (
                f"<div class='note'><strong>{e(figure_id)} was not recorded by the "
                "renderer.</strong> The figure manifest has no entry for it, so it was "
                "neither drawn nor skipped with a reason; re-render the run.</div>"
            )
        label = e(str(figure.get("label", figure_id)))
        caption = e(str(figure.get("caption", "")))
        if figure.get("status") != "rendered":
            reason = e(str(figure.get("reason", "no reason recorded")))
            return f"<div class='note'><strong>{label} was not drawn.</strong> {reason}</div>"
        png = root / str(figure["outputs"]["png"])
        if standalone:
            encoded = base64.b64encode(png.read_bytes()).decode("ascii")
            source = f"data:image/png;base64,{encoded}"
        else:
            source = e(str(figure["outputs"]["png"]))
        alt = e(str(figure.get("alt", caption)))
        vectors = " ".join(
            f"<code>{e(str(figure['outputs'][suffix]))}</code>"
            for suffix in ("svg", "pdf")
            if suffix in figure.get("outputs", {})
        )
        return (
            f"<figure><img src='{source}' alt='{alt}'>"
            f"<figcaption><b>{label}.</b> {caption} "
            f"<span class='meta'>Editable vector: {vectors}</span></figcaption></figure>"
        )

    def sources(*roles: str) -> str:
        """Name the artifacts a section was read from, so no number is unattributable."""

        available = [
            str(artifacts[role]["path"]) for role in roles if role in artifacts
        ]
        if not available:
            return ""
        rendered = " ".join(
            f"<code>{e(path)}</code>"
            if standalone
            else f"<a href='{e(path)}'><code>{e(path)}</code></a>"
            for path in available
        )
        return f"<p class='meta'>Source: {rendered}</p>"

    parts: list[str] = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>Transformation targets — {e(str(root.name))}</title>"
        f"<style>{_STYLE}</style></head><body>",
        f"<h1>Transformation target discovery — {method}</h1>",
        "<p class='sub'>Knockouts ranked by how far they move the <b>source</b> metabolic "
        "state toward the <b>target</b> state.</p>",
    ]

    # --- 1. Summary --------------------------------------------------------------------
    parts.append("<h2>1. Summary</h2>")
    top = ranking[0]
    top_branch = ""
    top_decomposed = decomposed.get(str(top.get("target_id")))
    if top_decomposed is not None:
        top_branch = (
            " Its score comes from Equation 9&rsquo;s <b>amplification branch</b> "
            "(bTS&nbsp;&gt;&nbsp;0, mTS&nbsp;&gt;&nbsp;0, wTS&nbsp;&lt;&nbsp;0), so it is "
            "mTS multiplied by k&nbsp;&times;&nbsp;(bTS&nbsp;&minus;&nbsp;wTS)."
            if top_decomposed["amplified"]
            else " The amplification branch did not fire for it, so its score is its MOMA "
            "component (rTS&nbsp;=&nbsp;mTS)."
        )
    parts.append(
        f"<p>{method} scored <b>{len(ranking)}</b> knockout candidates. The highest-scoring "
        f"candidate is <b>{e(str(top.get('target_id')))}</b> at "
        f"{_number(top.get('score'))}.{top_branch} <b>This candidate count is the denominator "
        "of any &ldquo;ranked in the top <i>N</i>%&rdquo; statement about this run</b>, so it "
        "belongs beside any percentile quoted from it.</p>"
    )
    if decomposition:
        parts.append(
            f"<p><b>{len(amplified_rows)}</b> of {len(decomposition)} scored candidates "
            f"{'satisfies' if len(amplified_rows) == 1 else 'satisfy'} rMTA&rsquo;s "
            "amplification branch (Equation 9: bTS&nbsp;&gt;&nbsp;0, "
            "mTS&nbsp;&gt;&nbsp;0, wTS&nbsp;&lt;&nbsp;0). Every other candidate is ranked by "
            "its MOMA component alone.</p>"
        )
    parts.append(
        "<div class='note'><strong>The direction is an input, not a finding.</strong> "
        "The same pair of files asks a different question when exchanged: "
        "source&nbsp;→&nbsp;target asks which perturbation <em>produces</em> the second state, "
        "and the reverse asks which <em>reverts</em> it. Nothing in the model detects a swap."
        "</div>"
    )
    deviations = provenance.get("reference_state_deviations") or []
    parts.append(
        "<div class='note'><strong>How the reference state was built.</strong> As in "
        "Yizhak et al. (2013) and Valcárcel et al. (2019): source genes were classed high, "
        "moderate or low by the rMTA rule, iMAT found one activity pattern consistent with "
        f"those classes, and v_ref is the mean of {e(str(provenance.get('reference_samples')))} "
        "flux distributions sampled from that pattern&rsquo;s space. iMAT places no objective "
        "on growth. What still differs from the papers: "
        + "; ".join(e(str(item)) for item in deviations)
        + ".</div>"
    )
    sweep_note = (
        " Section 4 reports how far the ranking moves across other values."
        if sweep
        else " No sensitivity sweep was run, so this run does not show whether the value "
        "decided the answer."
    )
    if provenance.get("epsilon") == "per_reaction":
        eps_summary = direction.get("epsilon_summary") or {}
        parts.append(
            "<div class='note'><strong>Epsilon was derived per reaction.</strong> As in "
            "Yizhak et al. (2013), each changed reaction's threshold is the change that a "
            "one-sided t-test over the "
            f"{e(str(provenance.get('reference_samples')))} samples behind v_ref calls "
            "significant (P&nbsp;&lt;&nbsp;0.05), with a minimum of 0.001. Across "
            f"{e(str(eps_summary.get('n_reactions')))} reactions the median epsilon is "
            f"{_number(eps_summary.get('median'))}; "
            f"{e(str(eps_summary.get('n_floored')))} sit at the minimum."
            + sweep_note
            + "</div>"
        )
    else:
        parts.append(
            "<div class='note'><strong>Epsilon was chosen, not derived.</strong> It is how "
            "far a reaction's flux must move from v_ref before the move counts, measured in "
            "this model's flux units. The source papers derive it per reaction from a sampled "
            "reference distribution; this run used one value for every reaction instead."
            + sweep_note
            + "</div>"
        )
    if direction.get("changed_set_in_published_range") is False:
        parts.append(
            "<div class='note'><strong>The changed set is outside the published range."
            f"</strong> {e(str(direction.get('n_reactions_changed')))} reactions were "
            "labelled changed; Yizhak et al. and Valcárcel et al. tune their differential-"
            "expression filter until 100&ndash;200 are.</div>"
        )

    # --- 2. Setup ----------------------------------------------------------------------
    parts.append("<h2>2. Setup</h2>")
    parts.append(
        _table(
            [
                {
                    "field": "source (moving away from)",
                    "value": provenance.get("source_expression"),
                },
                {
                    "field": "target (moving toward)",
                    "value": provenance.get("target_expression"),
                },
                {"field": "model", "value": provenance.get("model_path")},
                {"field": "medium", "value": _medium_text(provenance, coverage)},
                {"field": "condition", "value": provenance.get("condition") or "none"},
                {"field": "method", "value": provenance.get("method")},
                {"field": "knockout level", "value": provenance.get("perturbation")},
                {
                    "field": "reference state (v_ref)",
                    "value": _reference_text(provenance),
                },
                {"field": "alpha", "value": provenance.get("alpha")},
                {"field": "epsilon", "value": provenance.get("epsilon")},
            ],
            ["field", "value"],
        )
    )
    warnings = [row for row in preflight if row.get("status") == "warning"]
    if warnings:
        # A warning is a fact about this run's inputs; left in the CSV, no reader sees it.
        parts.append(
            "<div class='note'><strong>Preflight warnings.</strong> The run continued, but "
            "these checks did not pass.</div>"
        )
        parts.append(_table(warnings, ["check", "message"]))
    parts.append(sources("workflow_configuration", "preflight"))

    # --- 3. Data and methods -----------------------------------------------------------
    parts.append("<h2>3. Data and methods</h2>")
    parts.append("<h3>Which reactions were asked to change</h3>")
    parts.append(
        _table(
            [
                {
                    "setting": "expression scale",
                    "value": provenance.get("expression_scale"),
                },
                {
                    "setting": "significance test",
                    "value": direction.get("significance"),
                },
                {
                    "setting": "P-value adjustment",
                    "value": direction.get("p_adjust") or "none",
                },
                {
                    "setting": "minimum |log2 fold change|",
                    "value": direction.get("min_abs_log2fc", 0.0),
                },
                {
                    "setting": "model genes with an expression value",
                    "value": _coverage_text(coverage),
                },
                {
                    "setting": "genes compared",
                    "value": direction.get("n_genes_compared"),
                },
                {
                    "setting": "genes significant",
                    "value": direction.get("n_genes_significant"),
                },
                {
                    "setting": "reactions labelled changed",
                    "value": direction.get("n_reactions_labelled"),
                },
                {
                    "setting": "changed set kept",
                    "value": direction.get("n_reactions_changed"),
                },
                {
                    "setting": "changed-set ranking",
                    "value": direction.get("changed_ranking")
                    or direction.get("ranking"),
                },
                *(
                    [
                        {
                            "setting": "requested moves impossible under the reaction "
                            "bounds (masked)",
                            "value": _masked_text(structure),
                        }
                    ]
                    if structure.get("n_impossible_masked") is not None
                    else []
                ),
            ],
            ["setting", "value"],
        )
    )
    parts.append(
        "<p class='meta'>A model gene with no expression value is treated as "
        "neither high nor low: iMAT classes an unmeasured gene as moderate, as the published "
        "pipeline does, and the gene contributes no direction. A requested move is masked, as "
        "in the published rMTA preprocessing, when it asks an inactive irreversible reaction "
        "to run backwards.</p>"
    )
    parts.append(
        sources(
            "source_reference_fluxes",
            "gene_expression_states",
            "imat_reference_pattern",
            "gene_differential_expression",
            "reaction_direction_map",
            "reaction_epsilons",
        )
    )
    parts.append("<h3>Candidate universe</h3>")
    parts.append(
        _table(
            [
                {"step": k.replace("_", " "), "value": v}
                for k, v in construction.items()
                if not isinstance(v, dict)
            ],
            ["step", "value"],
        )
    )
    coupling = construction.get("coupling")
    if isinstance(coupling, Mapping):
        parts.append(
            f"<p class='meta'>Coupling: {e(str(coupling.get('coupling')))}, "
            f"{e(str(coupling.get('n_sets')))} sets from "
            f"{e(str(coupling.get('n_reactions')))} reactions. CMM computes full coupling, "
            "which is stronger than partial coupling: it can split a partially coupled set "
            "but never merge two, so this count is an upper bound on a partial-coupling "
            "one.</p>"
        )
    parts.append(sources("transformation_candidates"))

    # --- 4. Results --------------------------------------------------------------------
    parts.append("<h2>4. Results</h2>")
    parts.append("<h3>Ranking</h3>")
    parts.append(panel("fig01_component_scores"))
    parts.append(panel("fig02_transformation_ranking"))
    columns = ["rank", "target_id"]
    if "target_name" in ranking[0]:
        columns.append("target_name")
    columns.append("score")
    for optional in ("bTS", "mTS", "wTS"):
        if optional in ranking[0]:
            columns.append(optional)
    parts.append(_table(list(ranking)[:top_n], columns, mark=highlight))
    if "bTS" in ranking[0]:
        parts.append(
            "<p class='meta'>rMTA's Equation 9 branches on the signs of bTS, mTS and wTS, so "
            "the components are reported alongside the combined score: which branch fired "
            "cannot be recovered from the score alone.</p>"
        )
    if structure.get("n_distinct_scores") is not None:
        parts.append(
            _table(
                [
                    {"quantity": "scored candidates", "value": len(ranking)},
                    {
                        "quantity": "distinct scores",
                        "value": structure.get("n_distinct_scores"),
                    },
                    {
                        "quantity": "largest tie block",
                        "value": structure.get("largest_tie_block"),
                    },
                ],
                ["quantity", "value"],
            )
        )
        parts.append(
            "<p class='meta'>Ties break on target id, so a top-<i>k</i> that ends inside a "
            "tie block is alphabetical past the block&rsquo;s first member.</p>"
        )
    parts.append(sources("transformation_ranking"))
    if decomposition:
        parts.append("<h3>Amplified candidates</h3>")
        if amplified_rows:
            amplified_columns = ["target_id"]
            if "target_name" in amplified_rows[0]:
                amplified_columns.append("target_name")
            amplified_columns += ["rTS", "bTS", "mTS", "wTS", "rank_rTS", "rank_mTS"]
            parts.append(
                _table(
                    amplified_rows,
                    amplified_columns,
                    mark=highlight,
                )
            )
        else:
            parts.append(
                "<p>No candidate satisfies the amplification branch, so every score is the "
                "candidate&rsquo;s MOMA component and the ranking equals the mTS ranking.</p>"
            )
        parts.append(sources("score_decomposition"))

    parts.append("<h3>Score decomposition</h3>")
    parts.append(panel("fig03_ranking_vs_moma"))
    if decomposition:
        parts.append(
            "<p>When the amplification branch does not fire, rTS&nbsp;=&nbsp;mTS, so this "
            "panel shows how much of the ranking is the MOMA component alone. Spearman&rsquo;s "
            f"&rho; between rTS and mTS is <b>{_number(structure.get('spearman_rTS_vs_mTS'), 3)}"
            "</b>. It is computed from the scores already in the ranking, with no extra "
            "solve.</p>"
        )
        parts.append(sources("score_decomposition"))
    if baseline:
        parts.append(
            "<p>An external MOMA baseline was also run: target expression put through the "
            "same discretisation, iMAT and sampling as v_ref, and each knockout scored by how close MOMA brings it to "
            "that target flux state. Yizhak et al. report such a baseline as markedly inferior "
            "for this task; here it is a cross-check and negative control rather than a "
            "superiority claim, and a high rank still says a knockout <i>could</i> produce the "
            "difference, not that MOMA <i>would</i> fail to find it. Spearman&rsquo;s &rho; "
            "between the transformation rank and the MOMA-baseline rank is "
            f"<b>{_number(structure.get('spearman_rank_vs_moma_baseline'), 3)}</b>.</p>"
        )
        failed_baseline = [
            row
            for row in baseline
            if row.get("status") not in (None, "", "optimal")
            or not math.isfinite(float(row["moma_score"]))
        ]
        if failed_baseline:
            parts.append(
                "<p>Unsuccessful MOMA solves remain in the CSV with their status and "
                "are excluded from the rank comparison.</p>"
            )
            parts.append(_table(failed_baseline, ["target_id", "status"]))
        parts.append(sources("moma_baseline"))
    elif decomposition:
        parts.append(
            "<p class='meta'>The external MOMA baseline was not run. The decomposition above "
            "answers the same question without extra solves; the baseline is an optional "
            "cross-check (<code>validation.run_moma_baseline</code>).</p>"
        )
    else:
        parts.append(
            "<p class='meta'>MTA computes one score per candidate, so there is no internal "
            "decomposition, and the optional external MOMA baseline "
            "(<code>validation.run_moma_baseline</code>) was not run.</p>"
        )

    parts.append("<h3>Epsilon sensitivity</h3>")
    if sweep:
        parts.append(panel("fig04_epsilon_sensitivity"))
        parts.append(sources("epsilon_sensitivity"))
    else:
        parts.append(
            "<div class='note'><strong>Not run.</strong> "
            + (
                "Epsilon was derived per reaction, so no single value was chosen; a sweep "
                "would scale that derivation rather than test a choice.</div>"
                if provenance.get("epsilon") == "per_reaction"
                else "Epsilon was chosen rather than derived, so nothing here shows whether "
                "that choice decided the ranking.</div>"
            )
        )

    parts.append("<h3>Interpretation limits</h3>")
    parts.append(
        "<ul>"
        "<li>These are <i>in silico</i> hypotheses that prioritise experiments. They do not "
        "establish a wet-lab phenotype.</li>"
        "<li>A high rank says a knockout <i>could</i> produce the observed difference, not that "
        "it did.</li>"
        "<li>The ranking is conditioned on one medium and one reference state. A different "
        "condition is a different run.</li>"
        "<li>Candidates ranked near the true answer are typically those sharing its "
        "consequences, so the top of the list is a neighbourhood rather than a single call.</li>"
        "<li>Ties break on candidate id, so a slice taken inside a tie block is alphabetical "
        "rather than meaningful.</li>"
        + (
            "<li>Scores are not comparable across runs or implementations. The k multiplier "
            "makes rTS discontinuous at the amplification branch, so the gap between the "
            "amplified tier and the rest is set by k, not by biology. Compare ranks and branch "
            "membership.</li>"
            if decomposition
            else ""
        )
        + "<li>Genes without an expression value are unconstrained rather than low, so the "
        "scores of the reactions they encode rest on the network alone.</li>"
        "</ul>"
    )

    # --- 5. References -----------------------------------------------------------------
    parts.append("<h2>5. References</h2>")
    parts.append(
        "<p>MTA — Yizhak K, Gabay O, Cohen H, Ruppin E (2013) "
        "<i>Nature Communications</i> <b>4</b>:2632. rMTA — Valc&aacute;rcel LV <i>et al.</i> "
        "(2019) <i>Bioinformatics</i> <b>35</b>(21):4350&ndash;4355. The MOMA baseline is "
        "defined in the same 2013 paper. Neither method originates in CMM.</p>"
    )

    # --- 6. Provenance -----------------------------------------------------------------
    parts.append("<h2>6. Provenance</h2>")
    renderer = figure_manifest.renderer
    parts.append(
        _table(
            [
                {"field": "run directory", "value": root.name},
                {"field": "model fingerprint", "value": provenance.get("model_id")},
                {"field": "solver", "value": provenance.get("solver")},
                *(
                    [
                        {
                            "field": "medium as applied",
                            "value": _applied_medium_text(provenance["medium_applied"]),
                        }
                    ]
                    if isinstance(provenance.get("medium_applied"), Mapping)
                    else []
                ),
                {"field": "figure renderer", "value": renderer.get("engine")},
                {"field": "R", "value": renderer.get("r")},
                {"field": "renderer sha256", "value": renderer.get("script_sha256")},
            ],
            ["field", "value"],
        )
    )
    parts.append("<h3>Artifact inventory</h3>")
    parts.append(
        _table(
            [
                {
                    "role": role,
                    "path": str(entry.get("path", "—")),
                    "status": str(entry.get("status", "—")),
                }
                for role, entry in sorted(artifacts.items())
            ],
            ["role", "path", "status"],
        )
    )
    parts.append(
        "<p class='meta'>Every number on this page is read from that inventory&rsquo;s CSVs; "
        "<code>00_manifest.json</code> is the only path-discovery surface.</p>"
    )
    parts.append("</body></html>")
    return "\n".join(parts)


# --------------------------------------------------------------------------------------
# Completion gate
# --------------------------------------------------------------------------------------

#: Roles the workflow always writes. A run missing one of these did not finish, whatever else
#: is present. The optional roles below are configurable stages, so their absence is a choice.
_REQUIRED_ROLES: Mapping[str, tuple[str, ...]] = {
    "provenance": (),
    "summary": (),
    "workflow_configuration": (),
    "model": (),
    "source_expression": (),
    "target_expression": (),
    "reproduction_config": (),
    "reproduce_script": (),
    "preflight": ("check", "status"),
    "source_reference_fluxes": ("reaction_id", "flux"),
    "reaction_direction_map": ("reaction_id", "direction"),
    "transformation_candidates": ("target_id",),
    "transformation_ranking": ("target_id", "score", "rank"),
}
_OPTIONAL_ROLES: Mapping[str, tuple[str, ...]] = {
    "gene_differential_expression": ("gene_id",),
    "moma_baseline": ("target_id", "moma_score", "rank"),
    "epsilon_sensitivity": ("epsilon", "target_id", "rank"),
    "score_decomposition": (
        "target_id",
        "bTS",
        "mTS",
        "wTS",
        "rTS",
        "amplified",
        "rank_rTS",
        "rank_mTS",
    ),
    "gene_expression_states": ("gene_id", "state"),
    "imat_reference_pattern": ("reaction_id", "expression_state", "activity", "v_ref"),
    "reaction_epsilons": ("reaction_id", "epsilon"),
}


def validate_transformation_run(
    run_dir: str | Path,
    *,
    manifest_name: str = "00_manifest.json",
) -> ValidationReport:
    """Check that a transformation run is complete, self-consistent and readable.

    This is the completion gate. An HTML file that opens is not evidence the run finished: the
    figures it points at can be absent, the ranking it quotes can disagree with the CSV it
    claims to come from, and a manifest can name an artifact that was never written. Each of
    those reads as success in a browser.

    Never raises for a failed run — the issues are returned so a CLI or an agent can report
    them. Only a malformed argument raises.
    """

    root = Path(run_dir).expanduser().resolve()
    issues: list[str] = []
    warnings: list[str] = []

    manifest_path = root / manifest_name
    if not manifest_path.is_file():
        return ValidationReport(
            valid=False,
            issues=(f"required manifest is missing: {manifest_path}",),
            warnings=(),
            phase="pre-render",
        )
    try:
        manifest = _read_json(manifest_path)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        return ValidationReport(
            valid=False,
            issues=(f"{manifest_name} is not valid UTF-8 JSON: {exc}",),
            warnings=(),
            phase="pre-render",
        )
    if manifest.get("workflow") != "transformation_target_discovery":
        return ValidationReport(
            valid=False,
            issues=(
                f"{root} is a {manifest.get('workflow')!r} run, "
                "not a transformation_target_discovery run",
            ),
            warnings=(),
            phase="pre-render",
        )
    if manifest.get("schema_version") != 2:
        issues.append(
            f"schema_version must be the integer 2, got {manifest.get('schema_version')!r}"
        )
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, Mapping):
        return ValidationReport(
            valid=False,
            issues=("manifest artifacts must be an object keyed by semantic role",),
            warnings=(),
            phase="pre-render",
        )

    csv_rows: dict[str, list[dict[str, Any]]] = {}
    artifact_paths: dict[str, Path] = {}
    for role, columns in {**_REQUIRED_ROLES, **_OPTIONAL_ROLES}.items():
        entry = artifacts.get(role)
        if entry is None:
            if role in _REQUIRED_ROLES:
                issues.append(f"manifest is missing required artifact role {role!r}")
            continue
        if not isinstance(entry, Mapping) or not isinstance(entry.get("path"), str):
            issues.append(f"artifact {role!r} has no usable path")
            continue
        try:
            path = relative_path(
                root, entry["path"], label=role, error_type=TransformationReportError
            )
        except TransformationReportError as exc:
            issues.append(str(exc))
            continue
        if not path.is_file():
            issues.append(f"artifact {role!r} is declared but missing: {entry['path']}")
            continue
        _check_integrity(role, entry, path, issues=issues, warnings=warnings)
        artifact_paths[role] = path
        status = entry.get("status")
        if status not in {"complete", "partial"}:
            # A stage switched off writes an empty table and says so. Demanding its columns
            # would turn a deliberate configuration into a validation failure; the reason is
            # what has to be present, so the absence stays legible rather than silent.
            if role in _REQUIRED_ROLES:
                issues.append(f"required artifact {role!r} is {status!r}")
            elif not str(entry.get("reason", "")).strip():
                issues.append(
                    f"artifact {role!r} is {status!r} but records no reason, so a reader "
                    "cannot tell a disabled stage from a failed one"
                )
            continue
        if columns:
            rows = _read_csv(path)
            csv_rows[role] = rows
            missing = [name for name in columns if not rows or name not in rows[0]]
            if missing:
                issues.append(
                    f"artifact {role!r} is missing required column(s): {sorted(missing)}"
                )

    _check_reproduction_inputs(artifact_paths, issues=issues)
    _check_ranking_invariants(csv_rows, issues=issues, warnings=warnings)

    outputs = _check_report_outputs(root, issues=issues, warnings=warnings)
    checks: dict[str, Any] = {
        "artifact_contract": {
            "status": "fail" if issues else "pass",
            "roles_checked": len(csv_rows) + len(_REQUIRED_ROLES),
        },
        **outputs,
    }
    return ValidationReport(
        valid=not issues,
        issues=tuple(issues),
        warnings=tuple(warnings),
        phase="post-render" if (root / "report.html").is_file() else "pre-render",
        checks=checks,
        report_validation=(root / "report_validation.json")
        if (root / "report_validation.json").is_file()
        else None,
    )


def _check_reproduction_inputs(
    artifact_paths: Mapping[str, Path], *, issues: list[str]
) -> None:
    """Both configs must refer to the archived inputs, relative to their own location."""

    for role in ("workflow_configuration", "reproduction_config"):
        path = artifact_paths.get(role)
        if path is None:
            continue
        try:
            config = _read_json(path)
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
            issues.append(f"{role} is not valid UTF-8 JSON: {exc}")
            continue
        if not isinstance(config, Mapping):
            issues.append(f"{role} must contain a JSON object")
            continue
        for field, input_role in (
            ("model_path", "model"),
            ("source_expression_path", "source_expression"),
            ("target_expression_path", "target_expression"),
        ):
            value = config.get(field)
            if (
                not isinstance(value, str)
                or not value
                or Path(value).is_absolute()
                or (path.parent / value).resolve() != artifact_paths.get(input_role)
            ):
                issues.append(
                    f"{role}.{field} must reference the manifest's {input_role} "
                    "using a relative path; external inputs make the bundle unreproducible"
                )


def _check_integrity(
    role: str,
    entry: Mapping[str, Any],
    path: Path,
    *,
    issues: list[str],
    warnings: list[str],
) -> None:
    """Confirm the file on disk is the file the manifest describes.

    A stale hash is how a partially re-run bundle looks: the manifest was written by one run
    and the CSV beside it by another, and every number downstream is then attributed wrongly.
    """

    declared_size = entry.get("size_bytes")
    actual_size = path.stat().st_size
    if isinstance(declared_size, int) and declared_size != actual_size:
        issues.append(
            f"artifact {role!r} is {actual_size} bytes; the manifest declares {declared_size}"
        )
    declared_hash = entry.get("sha256")
    if isinstance(declared_hash, str) and declared_hash:
        if _sha256(path) != declared_hash:
            issues.append(
                f"artifact {role!r} does not match its manifest sha256; "
                "the bundle mixes files from more than one run"
            )
    else:
        warnings.append(
            f"artifact {role!r} declares no sha256, so it cannot be verified"
        )


def _check_ranking_invariants(
    csv_rows: Mapping[str, Sequence[Mapping[str, Any]]],
    *,
    issues: list[str],
    warnings: list[str],
) -> None:
    """Check the ranking against itself and against the artifacts that cite it."""

    ranking = csv_rows.get("transformation_ranking")
    if not ranking:
        return
    try:
        ranks = [int(row["rank"]) for row in ranking]
    except (KeyError, TypeError, ValueError):
        issues.append("transformation_ranking has a non-integer rank")
        return
    if ranks != list(range(1, len(ranks) + 1)):
        issues.append(
            "transformation_ranking is not ranked 1..N in file order; "
            "a reader taking the top rows would not get the top candidates"
        )
    scores = [float(row["score"]) for row in ranking]
    if scores != sorted(scores, reverse=True):
        issues.append("transformation_ranking rows are not in descending score order")

    distinct = len(set(scores))
    if distinct < len(scores):
        largest_tie = max(scores.count(score) for score in set(scores))
        warnings.append(
            f"the ranking holds {distinct} distinct scores across {len(scores)} candidates; "
            f"its largest tie block is {largest_tie}, and ties break on target id, so a "
            "top-k slice taken inside one is alphabetical rather than meaningful"
        )

    baseline = csv_rows.get("moma_baseline")
    if baseline:
        ranked = {str(row["target_id"]) for row in ranking}
        extra = {str(row["target_id"]) for row in baseline} - ranked
        if extra:
            issues.append(
                f"moma_baseline scores {len(extra)} candidate(s) absent from the ranking; "
                "the two are meant to cover the same candidate universe"
            )

    decomposition = csv_rows.get("score_decomposition")
    if decomposition:
        ranked = {str(row["target_id"]) for row in ranking}
        extra = {str(row["target_id"]) for row in decomposition} - ranked
        if extra:
            issues.append(
                f"score_decomposition holds {len(extra)} candidate(s) absent from the ranking"
            )


def _check_report_outputs(
    root: Path, *, issues: list[str], warnings: list[str]
) -> dict[str, Any]:
    """Audit the rendered outputs, including the failures a browser will not show."""

    figure_manifest_path = root / "figures" / "figure_manifest.json"
    linked = root / "report.html"
    standalone = root / "report_standalone.html"
    present = [
        path for path in (figure_manifest_path, linked, standalone) if path.is_file()
    ]
    if not present:
        warnings.append(
            "the run has not been rendered; run cmm.reporting.render_transformation_report"
        )
        return {"publication_outputs": {"status": "not-rendered"}}
    for path in (figure_manifest_path, linked, standalone):
        if not path.is_file():
            issues.append(f"required report output is missing: {path.name}")
    if issues:
        return {"publication_outputs": {"status": "failed"}}

    try:
        figure_manifest = load_figure_manifest(
            figure_manifest_path,
            root,
            required=_REQUIRED_FIGURES,
            order=_FIGURE_ORDER,
            error_type=TransformationReportError,
        )
    except TransformationReportError as exc:
        issues.append(str(exc))
        return {"publication_outputs": {"status": "failed"}}

    # Iterate the expected ids, not the manifest: an id the renderer never recorded is
    # otherwise invisible to every check below, because each of them walks the manifest.
    recorded_ids = {str(figure["id"]) for figure in figure_manifest.figures}
    for figure_id in _FIGURE_ORDER:
        if figure_id not in recorded_ids:
            issues.append(
                f"figure manifest has no entry for {figure_id!r}; the renderer neither drew "
                "it nor recorded why it was skipped"
            )
    rendered = [
        figure
        for figure in figure_manifest.figures
        if figure.get("status") == "rendered"
    ]
    linked_text = linked.read_text(encoding="utf-8")
    standalone_text = standalone.read_text(encoding="utf-8")
    for figure in rendered:
        png = str(figure["outputs"]["png"])
        if png not in linked_text:
            issues.append(
                f"report.html does not place rendered figure {figure['id']!r}; "
                "the page is stale with respect to the figure manifest"
            )
    # The whole reason two variants exist. A relative src in the standalone copy renders as
    # blank space in silence once the file has been moved off this machine.
    for source in re.findall(r"<img[^>]*\ssrc='([^']*)'", standalone_text):
        if not source.startswith("data:image/png;base64,"):
            issues.append(
                f"report_standalone.html references an image it does not carry: {source[:60]}"
            )
    for source in re.findall(r"<img[^>]*\ssrc=\"?'?([^'\">]*)", linked_text):
        if source.startswith("data:"):
            warnings.append(
                "report.html embeds an image; the linked variant is meant to stay small and "
                "reference figures/ so the 300-dpi originals are the ones read"
            )
            break

    validation_path = root / "report_validation.json"
    if validation_path.is_file():
        record = _read_json(validation_path)
        recorded = set(record.get("figures_rendered", ()))
        if recorded != {str(figure["id"]) for figure in rendered}:
            issues.append(
                "report_validation.json disagrees with the figure manifest about what was "
                "rendered; the run was re-rendered without re-validating"
            )
    else:
        issues.append("report_validation.json is missing")

    return {
        "publication_outputs": {
            "status": "fail" if issues else "pass",
            "figures_rendered": len(rendered),
            "figures_unavailable": len(figure_manifest.figures) - len(rendered),
        }
    }
