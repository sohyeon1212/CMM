"""Report rendering for transformation-target runs.

Figures are matplotlib and the page is built here in Python, so this renderer has no R
dependency — the production report's ``nature-r`` backend draws production-specific panels
(yield, envelope, strain design) that have no counterpart in a transformation run.

The page is deliberately opinionated about what it must say. A transformation ranking is easy
to over-read: it is produced from a reference state that is not the published one, with an
epsilon that was chosen rather than derived, over a candidate set whose size is the denominator
of any percentile claim. Those three facts are rendered as prose in the report body, not left
in a provenance file for a reader to find.
"""

from __future__ import annotations

import base64
from dataclasses import dataclass
import html
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


@dataclass(frozen=True)
class TransformationReport:
    """Where the rendered report and its figures were written."""

    run_directory: Path
    #: References figures by relative path. Small, and the 300-dpi originals stay beside it.
    report_html: Path
    #: Every figure inlined as a data URI. This is the copy to send someone: the relative-path
    #: copy loses all of its figures the moment it leaves the run directory, and says nothing.
    report_standalone_html: Path
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


def render_transformation_report(
    run_dir: str | Path,
    *,
    highlight: str | None = None,
    top_n: int = 15,
) -> TransformationReport:
    """Render figures and a self-contained HTML page from a transformation run bundle.

    ``highlight`` names a candidate to mark throughout — the knockout under test in a
    validation run, or the one a user came to check.
    """

    from cmm.visualization.figures import (
        epsilon_sensitivity_figure,
        save_figure,
        transformation_ranking_figure,
        transformation_vs_moma_figure,
    )

    root = Path(run_dir).resolve()
    manifest_path = root / "00_manifest.json"
    if not manifest_path.is_file():
        raise TransformationReportError(f"no 00_manifest.json in {root}")
    manifest = _read_json(manifest_path)
    if manifest.get("workflow") != "transformation_target_discovery":
        raise TransformationReportError(
            f"{root} is a {manifest.get('workflow')!r} run, not a transformation run"
        )
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

    summary = _read_json(required_path("summary"))
    provenance = _read_json(required_path("provenance"))
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

    figures_dir = root / "figures"
    figures_dir.mkdir(exist_ok=True)
    written: list[Path] = []

    def emit(figure, stem: str) -> None:
        """Write one panel as raster and as vector.

        The PNG is what the page shows; the SVG and PDF are what a manuscript needs, and a
        figure that only exists as a raster has to be redrawn to be used anywhere else.
        """

        for suffix in ("png", "svg", "pdf"):
            path = save_figure(figure, figures_dir / f"{stem}.{suffix}")
            if suffix == "png":
                written.append(path)

    emit(
        transformation_ranking_figure(ranking, highlight=highlight, top_n=10), "ranking"
    )
    if baseline:
        emit(
            transformation_vs_moma_figure(ranking, baseline, highlight=highlight),
            "ranking_vs_moma",
        )
    if sweep:
        emit(
            epsilon_sensitivity_figure(sweep, highlight=highlight),
            "epsilon_sensitivity",
        )

    html_text = _compose(
        root=root,
        summary=summary,
        provenance=provenance,
        ranking=ranking,
        baseline=baseline,
        sweep=sweep,
        figures=written,
        highlight=highlight,
        top_n=top_n,
    )
    report_path = root / "report.html"
    report_path.write_text(html_text, encoding="utf-8")
    standalone_path = root / "report_standalone.html"
    standalone_path.write_text(_standalone(html_text, root), encoding="utf-8")
    return TransformationReport(
        run_directory=root,
        report_html=report_path,
        report_standalone_html=standalone_path,
        figures=tuple(written),
    )


def _standalone(html_text: str, root: Path) -> str:
    """Rewrite the page so every figure travels inside it.

    ``report.html`` points at ``figures/*.png`` by relative path, which is right for reading the
    run in place and useless the moment the file is sent to someone: the images resolve against
    wherever the copy landed, and a browser renders the missing ones as blank space without
    saying anything. This produces the copy that survives being detached from its directory.

    Raises if any relative reference survives, because that is precisely the failure that is
    invisible in the rendered output.
    """

    out = html_text
    for png in sorted((root / "figures").glob("*.png")):
        encoded = base64.b64encode(png.read_bytes()).decode("ascii")
        out = out.replace(
            f"src='figures/{png.name}'", f"src='data:image/png;base64,{encoded}'"
        )
    if "figures/" in out:
        raise ValueError(
            "standalone report still references a figure by relative path; "
            "it would render with that figure silently missing"
        )
    return out


def _compose(
    *,
    root: Path,
    summary: Mapping[str, Any],
    provenance: Mapping[str, Any],
    ranking: Sequence[Mapping[str, Any]],
    baseline: Sequence[Mapping[str, Any]],
    sweep: Sequence[Mapping[str, Any]],
    figures: Sequence[Path],
    highlight: str | None,
    top_n: int,
) -> str:
    e = html.escape
    method = str(summary.get("method", "?")).upper()
    construction = summary.get("candidate_construction") or {}
    direction = summary.get("direction_construction") or {}
    figure_names = {path.name for path in figures}

    parts: list[str] = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>Transformation targets — {e(str(root.name))}</title>"
        f"<style>{_STYLE}</style></head><body>",
        f"<h1>Transformation target discovery — {method}</h1>",
        "<p class='sub'>Knockouts ranked by how far they move the <b>source</b> metabolic "
        "state toward the <b>target</b> state.</p>",
    ]

    # --- what was asked ----------------------------------------------------------------
    parts.append("<h2>1. What was asked</h2>")
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
                {
                    "field": "medium",
                    "value": provenance.get("medium") or "model as loaded",
                },
                {"field": "condition", "value": provenance.get("condition") or "none"},
                {"field": "method", "value": provenance.get("method")},
                {"field": "knockout level", "value": provenance.get("perturbation")},
            ],
            ["field", "value"],
        )
    )
    parts.append(
        "<div class='note'><strong>The direction is an input, not a finding.</strong> "
        "The same pair of files asks a different question when exchanged: source&nbsp;→&nbsp;target "
        "asks which perturbation <em>produces</em> the second state, and the reverse asks which "
        "<em>reverts</em> it. Nothing in the model detects a swap.</div>"
    )

    # --- how it was set up -------------------------------------------------------------
    parts.append("<h2>2. How the run was set up</h2>")
    parts.append(
        _table(
            [
                {
                    "setting": "reference state (v_ref)",
                    "value": provenance.get("reference_method"),
                },
                {
                    "setting": "significance test",
                    "value": direction.get("significance"),
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
                {"setting": "alpha", "value": provenance.get("alpha")},
                {"setting": "epsilon", "value": provenance.get("epsilon")},
            ],
            ["setting", "value"],
        )
    )
    parts.append(
        "<div class='note'><strong>The reference state is not the published one.</strong> "
        f"{e(str(provenance.get('reference_state_deviation', '')))}. iMAT places no objective "
        "on growth, so the ranking below is conditioned on the estimator that was used.</div>"
    )
    parts.append(
        "<div class='note'><strong>Epsilon was chosen, not derived.</strong> It is how far a "
        "reaction's flux must move from v_ref before the move counts, measured in this model's "
        "flux units. The source papers derive it per data set from a sampled reference "
        "distribution, which CMM does not reproduce."
        + (
            " Section 6 reports how far the ranking moves across other values."
            if sweep
            else " No sensitivity sweep was run, so this run does not show whether the value "
            "decided the answer."
        )
        + "</div>"
    )

    # --- candidates ---------------------------------------------------------------------
    parts.append("<h2>3. Candidate universe</h2>")
    rows = [
        {"step": k.replace("_", " "), "value": v}
        for k, v in construction.items()
        if not isinstance(v, dict)
    ]
    parts.append(_table(rows, ["step", "value"]))
    parts.append(
        f"<p><b>{len(ranking)}</b> candidates were scored. "
        "<b>This count is the denominator of any &ldquo;ranked in the top <i>N</i>%&rdquo; "
        "statement about this run</b>, so it belongs beside any percentile quoted from it.</p>"
    )
    coupling = construction.get("coupling")
    if isinstance(coupling, Mapping):
        parts.append(
            f"<p class='meta'>Coupling: {e(str(coupling.get('coupling')))}, "
            f"{e(str(coupling.get('n_sets')))} sets from "
            f"{e(str(coupling.get('n_reactions')))} reactions. CMM computes full coupling, "
            "which is stronger than the partial coupling the source paper uses: it can split "
            "one of their sets but never merge two.</p>"
        )

    # --- ranking -------------------------------------------------------------------------
    parts.append("<h2>4. Ranking</h2>")
    if "ranking.png" in figure_names:
        parts.append(
            "<figure><img src='figures/ranking.png' alt='transformation score by rank'>"
            "<figcaption>Transformation score against rank. The band marks the leading "
            "candidates; a score that decays smoothly rather than falling away means the cut "
            "between &ldquo;top&rdquo; and &ldquo;rest&rdquo; is a choice, not a "
            "boundary.</figcaption></figure>"
        )
    columns = ["rank", "target_id", "score"]
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

    # --- MOMA baseline -------------------------------------------------------------------
    parts.append("<h2>5. MOMA baseline</h2>")
    if baseline:
        if "ranking_vs_moma.png" in figure_names:
            parts.append(
                "<figure><img src='figures/ranking_vs_moma.png' alt='transformation rank "
                "against MOMA rank'><figcaption>Each candidate's rank under the two methods. "
                "Points on the diagonal are candidates the two agree on.</figcaption></figure>"
            )
        parts.append(
            "<p>Yizhak et al. compare their method against a MOMA baseline and report it as "
            "<i>markedly inferior</i> for this task. The comparison is here because a ranking "
            "that merely reproduces MOMA's ordering has not shown that its signal comes from "
            "the method rather than from the inputs.</p>"
        )
    else:
        parts.append(
            "<div class='note'><strong>Not run.</strong> Without it this run does not "
            "demonstrate that its ordering differs from the baseline the source paper reports "
            "as inferior.</div>"
        )

    # --- epsilon sweep --------------------------------------------------------------------
    if sweep:
        parts.append("<h2>6. Epsilon sensitivity</h2>")
        if "epsilon_sensitivity.png" in figure_names:
            parts.append(
                "<figure><img src='figures/epsilon_sensitivity.png' alt='rank against "
                "epsilon'><figcaption>How each leading candidate's rank moves with epsilon. "
                "A candidate whose rank is flat across the sweep was not produced by the value "
                "chosen.</figcaption></figure>"
            )

    # --- limits ----------------------------------------------------------------------------
    parts.append("<h2>Interpretation limits</h2>")
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
        "</ul>"
    )
    parts.append(
        "<h2>Methods to cite</h2><p>MTA — Yizhak K, Gabay O, Cohen H, Ruppin E (2013) "
        "<i>Nature Communications</i> <b>4</b>:2632. rMTA — Valc&aacute;rcel LV <i>et al.</i> "
        "(2019) <i>Bioinformatics</i> <b>35</b>(21):4350&ndash;4355. The MOMA baseline is "
        "defined in the same 2013 paper. Neither method originates in CMM.</p>"
    )
    parts.append(
        f"<p class='meta'>Run directory <code>{e(str(root))}</code>. "
        "Every number here is read from that bundle's CSVs; "
        "<code>00_manifest.json</code> is the only path-discovery surface.</p>"
    )
    parts.append("</body></html>")
    return "\n".join(parts)
