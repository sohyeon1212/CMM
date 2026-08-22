"""R-rendered figures, journal-style HTML, and strict post-render validation."""

from __future__ import annotations

import base64
import csv
from dataclasses import dataclass
import hashlib
import html
from html.parser import HTMLParser
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import struct
import subprocess
from typing import Any, Mapping, Sequence
from urllib.parse import urlparse
import xml.etree.ElementTree as ET

from cmm.reporting.schema import (
    RUN_SCHEMA_VERSION,
    RunValidationError,
    ValidationReport,
    ValidatedRun,
    validate_run,
)

_FIGURE_ORDER = (
    "fig01_yield_envelope",
    "fig02_single_knockout",
    "fig03_strain_design",
    "fig04_amplification",
    "fig05_flux_response",
    "fig06_sampling_shift",
)
_REQUIRED_FIGURES = frozenset(_FIGURE_ORDER[:4])
_REFERENCES = (
    (
        "MOMA",
        "Analysis of optimality in natural and perturbed metabolic networks",
        "Proceedings of the National Academy of Sciences (2002)",
        "10.1073/pnas.232349399",
    ),
    (
        "ROOM",
        "Regulatory on/off minimization of metabolic flux changes after genetic perturbations",
        "Proceedings of the National Academy of Sciences (2005)",
        "10.1073/pnas.0406346102",
    ),
    (
        "OptKnock",
        "OptKnock: a bilevel programming framework for identifying gene knockout strategies for microbial strain optimization",
        "Biotechnology and Bioengineering (2003)",
        "10.1002/bit.10803",
    ),
    (
        "RobustKnock",
        "Predicting metabolic engineering knockout strategies for chemical production: accounting for competing pathways",
        "Bioinformatics (2010)",
        "10.1093/bioinformatics/btp704",
    ),
    (
        "FSEOF",
        "In silico identification of gene amplification targets for improvement of lycopene production",
        "Applied and Environmental Microbiology (2010)",
        "10.1128/AEM.00115-10",
    ),
    (
        "FVSEOF",
        "Flux variability scanning based on enforced objective flux for identifying gene amplification targets",
        "BMC Systems Biology (2012)",
        "10.1186/1752-0509-6-106",
    ),
    (
        "StrainDesign",
        "StrainDesign: a comprehensive Python package for computational design of metabolic networks",
        "Bioinformatics (2022)",
        "10.1093/bioinformatics/btac632",
    ),
    (
        "Fast-SNP",
        "Fast-SNP: a fast matrix pre-processing algorithm for efficient loopless flux optimization of metabolic models",
        "Bioinformatics (2016)",
        "10.1093/bioinformatics/btw555",
    ),
    (
        "Randomized flux-space sampling",
        "Use of randomized sampling for analysis of metabolic networks",
        "Journal of Biological Chemistry (2009)",
        "10.1074/jbc.R800048200",
    ),
)
_EXPECTED_DOIS = tuple(reference[3] for reference in _REFERENCES)
_REPORT_VALIDATION_NAME = "report_validation.json"


@dataclass(frozen=True)
class FigureManifest:
    path: Path
    renderer: Mapping[str, Any]
    figures: tuple[Mapping[str, Any], ...]


@dataclass(frozen=True)
class ReportBuildResult:
    report_html: Path
    report_standalone_html: Path
    figure_manifest: Path
    report_validation: Path
    validation_warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class PublicationBundle:
    validated_run: ValidatedRun
    figures: FigureManifest
    report: ReportBuildResult


class FigureRenderError(RuntimeError):
    """The external renderer or generated publication bundle failed validation."""


def renderer_script_path() -> Path:
    """Return the checked-in R renderer used by :func:`render_publication_figures`."""

    return Path(__file__).with_name("render_publication_figures.R")


def _relative_path(root: Path, value: object, *, label: str) -> Path:
    if not isinstance(value, str) or not value:
        raise FigureRenderError(f"{label} must be a non-empty relative path")
    pure = PurePosixPath(value)
    if pure.is_absolute() or ".." in pure.parts:
        raise FigureRenderError(f"{label} escapes the run directory: {value!r}")
    path = (root / Path(*pure.parts)).resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise FigureRenderError(
            f"{label} escapes the run directory: {value!r}"
        ) from exc
    return path


def _load_figure_manifest(path: Path, root: Path) -> FigureManifest:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FigureRenderError(
            f"figure manifest is not valid JSON: {path}: {exc}"
        ) from exc
    if not isinstance(data, Mapping):
        raise FigureRenderError("figure manifest must contain a JSON object")
    schema_version = data.get("schema_version")
    if type(schema_version) is not int or schema_version != RUN_SCHEMA_VERSION:
        raise FigureRenderError(
            f"figure manifest schema_version must be the integer {RUN_SCHEMA_VERSION}"
        )
    renderer = data.get("renderer", {})
    figures_value = data.get("figures")
    if not isinstance(renderer, Mapping) or not isinstance(figures_value, list):
        raise FigureRenderError(
            "figure manifest needs renderer object and figures list"
        )
    script_sha256 = renderer.get("script_sha256")
    if (
        not isinstance(script_sha256, str)
        or re.fullmatch(r"[0-9a-f]{64}", script_sha256) is None
    ):
        raise FigureRenderError(
            "figure manifest renderer.script_sha256 must be a lowercase SHA-256 digest"
        )

    figures: list[Mapping[str, Any]] = []
    seen: set[str] = set()
    rendered_required: set[str] = set()
    for index, raw in enumerate(figures_value):
        if not isinstance(raw, Mapping):
            raise FigureRenderError(f"figure manifest entry {index} is not an object")
        figure_id = raw.get("id")
        status = raw.get("status")
        if not isinstance(figure_id, str) or status not in {
            "rendered",
            "skipped",
            "failed",
        }:
            raise FigureRenderError(
                f"figure manifest entry {index} has invalid id/status"
            )
        if figure_id in seen:
            raise FigureRenderError(f"figure manifest repeats id {figure_id!r}")
        seen.add(figure_id)
        if status == "rendered":
            outputs = raw.get("outputs")
            sources = raw.get("sources")
            if not isinstance(outputs, Mapping) or not isinstance(sources, list):
                raise FigureRenderError(
                    f"rendered figure {figure_id!r} needs outputs and sources"
                )
            for suffix in ("png", "pdf", "svg"):
                output = _relative_path(
                    root, outputs.get(suffix), label=f"{figure_id}.{suffix}"
                )
                if not output.is_file() or output.stat().st_size == 0:
                    raise FigureRenderError(
                        f"rendered figure {figure_id!r} is missing non-empty {suffix}: {output}"
                    )
            for source in sources:
                source_path = _relative_path(root, source, label=f"{figure_id} source")
                if not source_path.is_file():
                    raise FigureRenderError(
                        f"figure {figure_id!r} cites a missing source artifact: {source_path}"
                    )
            if figure_id in _REQUIRED_FIGURES:
                rendered_required.add(figure_id)
        else:
            reason = raw.get("reason")
            if not isinstance(reason, str) or not reason.strip():
                raise FigureRenderError(
                    f"unavailable figure {figure_id!r} must state a reason"
                )
            if figure_id in _REQUIRED_FIGURES:
                raise FigureRenderError(
                    f"required figure {figure_id!r} was {status}: {reason}"
                )
        figures.append(raw)

    missing = sorted(_REQUIRED_FIGURES - rendered_required)
    if missing:
        raise FigureRenderError(
            f"figure manifest omitted required figure(s): {missing}"
        )
    if seen - set(_FIGURE_ORDER):
        raise FigureRenderError(
            f"figure manifest contains unknown figure id(s): {sorted(seen - set(_FIGURE_ORDER))}"
        )
    return FigureManifest(path=path, renderer=renderer, figures=tuple(figures))


def _renderer_environment(
    script_sha256: str, *, platform_name: str | None = None
) -> dict[str, str]:
    """Build a deterministic renderer environment without assuming POSIX locales on Windows."""

    platform = os.name if platform_name is None else platform_name
    environment = dict(os.environ)
    environment.update(
        {
            "TZ": "UTC",
            "LANG": "C.UTF-8",
            "CMM_RENDERER_SHA256": script_sha256,
        }
    )
    if platform != "nt":
        environment["LC_ALL"] = "C.UTF-8"
    else:
        environment.pop("LC_ALL", None)
    return environment


def _decode_renderer_stream(value: bytes | str | None) -> str:
    """Decode external renderer output explicitly and loss-tolerantly as UTF-8."""

    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


def render_publication_figures(
    run_dir: str | Path,
    *,
    manifest_name: str = "00_manifest.json",
    rscript: str | Path = "Rscript",
    renderer: str | Path | None = None,
    output_dir: str = "figures",
) -> FigureManifest:
    """Validate source artifacts and render publication artwork with checked-in R code."""

    validated = validate_run(run_dir, manifest_name=manifest_name)
    script = Path(renderer) if renderer is not None else renderer_script_path()
    script = script.expanduser().resolve()
    if not script.is_file():
        raise FigureRenderError(f"R publication renderer is missing: {script}")
    executable = shutil.which(os.fspath(rscript))
    if executable is None:
        raise FigureRenderError(
            f"Rscript executable {os.fspath(rscript)!r} was not found; publication figures require R"
        )

    figures_dir = _relative_path(
        validated.root, output_dir, label="figure output directory"
    )
    figures_dir.mkdir(parents=True, exist_ok=True)
    figure_manifest_path = figures_dir / "figure_manifest.json"
    command = [
        executable,
        "--vanilla",
        os.fspath(script),
        os.fspath(validated.root),
        os.fspath(validated.manifest_path),
        os.fspath(figures_dir),
        os.fspath(figure_manifest_path),
    ]
    environment = _renderer_environment(_sha256(script))
    completed = subprocess.run(
        command,
        cwd=validated.root,
        env=environment,
        check=False,
        capture_output=True,
    )
    renderer_output = "\n".join(
        output
        for value in (completed.stdout, completed.stderr)
        if (output := _decode_renderer_stream(value).strip())
    )
    if completed.returncode != 0:
        raise FigureRenderError(
            f"R publication renderer failed ({completed.returncode}): "
            f"{renderer_output or 'no renderer output'}"
        )
    if re.search(r"(^|\n)Warning(?: message)?", renderer_output, flags=re.IGNORECASE):
        raise FigureRenderError(
            f"R publication renderer emitted a warning: {renderer_output}"
        )
    if not figure_manifest_path.is_file():
        raise FigureRenderError(
            "R renderer completed without writing figures/figure_manifest.json"
        )
    return _load_figure_manifest(figure_manifest_path, validated.root)


def _read_json_object(path: Path) -> Mapping[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    return value if isinstance(value, Mapping) else {}


def _read_csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or ())
        return fields, [dict(row) for row in reader]


def _rows(validated: ValidatedRun, role: str) -> list[dict[str, str]]:
    path = validated.artifact(role, required=False)
    return _read_csv_rows(path)[1] if path is not None else []


def _json(validated: ValidatedRun, role: str) -> Mapping[str, Any]:
    path = validated.artifact(role, required=False)
    return _read_json_object(path) if path is not None else {}


def _escape(value: object) -> str:
    return html.escape(str(value), quote=True)


def _human_label(value: str) -> str:
    return value.replace("_", " ").strip().capitalize()


def _display_value(value: object) -> object:
    """Return a deterministic, readable representation for nested setup values."""

    if isinstance(value, (Mapping, list)):
        return json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(", ", ": "),
        )
    return value


def _format_number(value: object, *, digits: int = 4) -> str:
    """Format a finite numeric value concisely without inventing precision."""

    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return "—" if value is None or value == "" else str(value)
    if not (float("-inf") < number < float("inf")):
        return "—"
    return f"{number:.{digits}g}"


def _format_growth_retained(value: object) -> str:
    try:
        fraction = float(str(value).strip())
    except (TypeError, ValueError):
        return "Not quantified"
    if not (0 <= fraction <= 1):
        return "Not quantified"
    return f"{100 * fraction:.3g}%"


def _format_signed_number(value: object, *, digits: int = 4) -> str:
    """Format a finite effect size with an explicit sign when it is positive."""

    rendered = _format_number(value, digits=digits)
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return rendered
    return f"+{rendered}" if number > 0 and rendered != "—" else rendered


def _plural_noun(count: int, singular: str, plural: str | None = None) -> str:
    """Return publication-ready singular/plural wording without placeholder syntax."""

    return singular if count == 1 else (plural or f"{singular}s")


def _publication_prose(value: object) -> str:
    """Resolve machine-facing ``noun(s)`` placeholders before publication display."""

    text = str(value)
    pattern = re.compile(r"([A-Za-z][A-Za-z-]*)\(s\)")

    def replace(match: re.Match[str]) -> str:
        clause = re.split(r"[.;:]", text[: match.start()])[-1]
        counts = re.findall(r"\b\d+\b", clause)
        count = int(counts[-1]) if counts else 2
        return _plural_noun(count, match.group(1))

    return pattern.sub(replace, text)


def _semantic_warning_key(value: str) -> str:
    """Collapse equivalent condition warnings without hiding distinct risks."""

    normalized = re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()
    words = set(normalized.split())
    if "medium" in words and (
        {"loaded", "retained"} & words and {"model", "bounds"} & words
    ):
        return "condition:medium-as-loaded"
    return normalized


def _deduplicate_warnings(values: Sequence[str]) -> list[str]:
    unique: list[str] = []
    seen: set[str] = set()
    for value in values:
        key = _semantic_warning_key(value)
        if not key or key in seen:
            continue
        seen.add(key)
        unique.append(value)
    return unique


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def _table(
    headers: Sequence[str],
    rows: Sequence[Sequence[object]],
    *,
    css_class: str | None = None,
) -> str:
    head = "".join(f"<th>{_escape(header)}</th>" for header in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{_escape(value)}</td>" for value in row) + "</tr>"
        for row in rows
    )
    if not body:
        body = (
            f'<tr><td colspan="{max(1, len(headers))}" class="muted">'
            "No rows returned.</td></tr>"
        )
    class_attr = f' class="{_escape(css_class)}"' if css_class else ""
    return (
        f"<table{class_attr}><thead><tr>{head}</tr></thead>"
        f"<tbody>{body}</tbody></table>"
    )


def _artifact_ref(validated: ValidatedRun, role: str, *, standalone: bool) -> str:
    artifact = validated.artifacts.get(role)
    relative = artifact.relative_path if artifact is not None else None
    if not relative:
        return f"<code>{_escape(role)} unavailable</code>"
    label = f"<code>{_escape(relative)}</code>"
    return label if standalone else f'<a href="{_escape(relative)}">{label}</a>'


def _source_refs(
    validated: ValidatedRun, roles: Sequence[str], *, standalone: bool
) -> str:
    return ", ".join(
        _artifact_ref(validated, role, standalone=standalone) for role in roles
    )


def _flatten_scalars(
    value: object, *, prefix: str = "", depth: int = 0, limit: int = 120
) -> list[tuple[str, object]]:
    rows: list[tuple[str, object]] = []
    if depth > 4:
        return rows
    if isinstance(value, Mapping):
        for key in sorted(value, key=str):
            child = value[key]
            label = f"{prefix}.{key}" if prefix else str(key)
            rows.extend(
                _flatten_scalars(child, prefix=label, depth=depth + 1, limit=limit)
            )
            if len(rows) >= limit:
                break
    elif isinstance(value, list):
        if all(not isinstance(item, (Mapping, list)) for item in value):
            rows.append((prefix, "; ".join(str(item) for item in value)))
    elif isinstance(value, (str, int, float, bool)) or value is None:
        rows.append((prefix, value))
    return rows[:limit]


def _provenance_table(validated: ValidatedRun) -> str:
    provenance = _json(validated, "provenance")
    preferred = (
        "model_id",
        "model_sha256",
        "source_model_sha256",
        "conditioned_model_sha256",
        "solver",
        "solver_version",
        "python",
        "cmm",
        "cobra",
        "numpy",
        "pandas",
        "scipy",
        "timestamp_utc",
        "seed",
    )
    rows = [
        (_human_label(key), provenance[key]) for key in preferred if key in provenance
    ]
    return _table(("Field", "Value"), rows)


def _artifact_table(validated: ValidatedRun, *, standalone: bool) -> str:
    rows: list[tuple[object, ...]] = []
    for role in sorted(validated.artifacts):
        artifact = validated.artifacts[role]
        path = artifact.relative_path or "—"
        path_html = (
            f"<code>{_escape(path)}</code>"
            if standalone or path == "—"
            else f'<a href="{_escape(path)}"><code>{_escape(path)}</code></a>'
        )
        rows.append(
            (
                _human_label(role),
                artifact.status,
                path_html,
                _publication_prose(artifact.reason) if artifact.reason else "—",
            )
        )
    head = "".join(
        f"<th>{_escape(header)}</th>"
        for header in ("Artifact", "Status", "Path", "Reason")
    )
    body = "".join(
        "<tr>"
        f"<td>{_escape(role)}</td><td>{_escape(status)}</td><td>{path_html}</td>"
        f"<td>{_escape(reason)}</td></tr>"
        for role, status, path_html, reason in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _ko_display_ranks(
    validated: ValidatedRun, rows: Sequence[Mapping[str, str]]
) -> dict[int, int]:
    """Return presentation-only top-five ranks without changing workflow verdicts."""

    declared: dict[int, int] = {}
    for index, row in enumerate(rows):
        try:
            rank = int(float(row.get("display_rank", "")))
        except ValueError:
            continue
        if 1 <= rank <= 5:
            declared[index] = rank
    if declared:
        return declared

    configuration = _json(validated, "reproduction_config") or _json(
        validated, "workflow_configuration"
    )
    try:
        viability_fraction = float(configuration.get("viability_fraction", 0))
    except (TypeError, ValueError):
        viability_fraction = 0
    summary = _json(validated, "summary")
    try:
        wild_type_growth = float(summary.get("wild_type_growth", ""))
    except (TypeError, ValueError):
        wild_type_growth = math.nan

    def number(row: Mapping[str, str], column: str) -> float:
        try:
            value = float(row.get(column, ""))
        except (TypeError, ValueError):
            return math.nan
        return value if math.isfinite(value) else math.nan

    eligible: list[tuple[int, Mapping[str, str], float, float]] = []
    for index, row in enumerate(rows):
        growth = number(row, "objective")
        product = number(row, "product_flux")
        growth_fraction = number(row, "growth_fraction")
        if not math.isfinite(growth_fraction) and wild_type_growth > 0:
            growth_fraction = growth / wild_type_growth
        if (
            row.get("status", "").lower() == "optimal"
            and math.isfinite(growth)
            and math.isfinite(product)
            and math.isfinite(growth_fraction)
            and growth_fraction >= viability_fraction
        ):
            eligible.append((index, row, product, growth))
    eligible.sort(
        key=lambda item: (
            -item[2],
            -item[3],
            item[1].get("target_id", ""),
            item[0],
        )
    )
    return {index: rank for rank, (index, _, _, _) in enumerate(eligible[:5], 1)}


def _ko_mapping_index(
    validated: ValidatedRun,
) -> dict[str, tuple[str, tuple[tuple[str, str, str, str], ...]]]:
    names: dict[str, str] = {}
    reactions: dict[str, list[tuple[str, str, str, str]]] = {}
    for row in _rows(validated, "gene_knockout_mapping"):
        gene = row.get("gene_id", "")
        if not gene:
            continue
        if row.get("gene_name", ""):
            names.setdefault(gene, row["gene_name"])
        reaction = row.get("blocked_reaction", "")
        if reaction:
            item = (
                reaction,
                row.get("reaction_name", ""),
                row.get("reaction_equation", ""),
                row.get("gpr", ""),
            )
            values = reactions.setdefault(gene, [])
            if item not in values:
                values.append(item)
    return {
        gene: (names.get(gene, ""), tuple(reactions.get(gene, ())))
        for gene in names.keys() | reactions.keys()
    }


def _screen_table(validated: ValidatedRun, role: str, label: str) -> str:
    rows = _rows(validated, role)
    display_ranks = _ko_display_ranks(validated, rows)
    mapping = _ko_mapping_index(validated)
    supported = {
        row.get("target", "")
        for row in _rows(validated, "recommendations")
        if row.get("type") == "single_gene_knockout"
        and row.get("verdict", "").lower() == "support"
    }

    def finite_number(row: Mapping[str, str], column: str) -> float:
        try:
            number = float(row.get(column, ""))
        except (TypeError, ValueError):
            return float("-inf")
        return number if math.isfinite(number) else float("-inf")

    ordered = sorted(
        ((index, row) for index, row in enumerate(rows) if index in display_ranks),
        key=lambda indexed: (
            display_ranks.get(indexed[0], 99),
            -finite_number(indexed[1], "product_flux"),
            -finite_number(indexed[1], "objective"),
            indexed[1].get("target_id", ""),
            indexed[0],
        ),
    )
    table_rows: list[tuple[object, ...]] = []
    for index, row in ordered:
        target = row.get("target_id", "")
        mapped_name, mapped_reactions = mapping.get(target, ("", ()))
        gene_name = row.get("target_name", row.get("gene_name", mapped_name))
        gene_label = f"{target} ({gene_name})" if gene_name else target
        blocked_ids = [
            value.strip()
            for value in row.get("blocked_reactions", "").split(";")
            if value.strip()
        ]
        reaction_names = {reaction: name for reaction, name, _, _ in mapped_reactions}
        reaction_equations = {
            reaction: equation for reaction, _, equation, _ in mapped_reactions
        }
        direct_names = [
            value.strip() for value in row.get("blocked_reaction_names", "").split(";")
        ]
        direct_equations = [
            value.strip()
            for value in row.get("blocked_reaction_equations", "").split(";")
        ]
        if not blocked_ids:
            blocked_ids = [reaction for reaction, _, _, _ in mapped_reactions]
        for position, reaction in enumerate(blocked_ids):
            if position < len(direct_names) and direct_names[position]:
                reaction_names[reaction] = direct_names[position]
            if position < len(direct_equations) and direct_equations[position]:
                reaction_equations[reaction] = direct_equations[position]
        reaction_context = (
            "; ".join(
                (
                    f"{reaction} — {reaction_names[reaction]}"
                    if reaction_names.get(reaction)
                    else reaction
                )
                for reaction in blocked_ids
            )
            or "No model reaction blocked"
        )
        equation_context = (
            "; ".join(
                f"{reaction}: {reaction_equations[reaction]}"
                for reaction in blocked_ids
                if reaction_equations.get(reaction)
            )
            or "—"
        )
        evidence_status = (
            "Final support"
            if target in supported
            else "Beneficial screen candidate (forward-tested)"
            if _truthy(row.get("selected", ""))
            else "Display candidate (forward-tested)"
        )
        table_rows.append(
            (
                f"D{display_ranks[index]}" if index in display_ranks else "—",
                gene_label,
                reaction_context,
                equation_context,
                row.get("status", ""),
                _format_number(row.get("objective", "")),
                _format_number(row.get("product_flux", "")),
                evidence_status,
            )
        )
    note = (
        '<p class="muted"><strong>Display-rank definition.</strong> D1–D5 use the '
        "method-specific workflow display rank when exported; current exports choose one "
        "representative per blocked-reaction signature. For legacy runs without that field, "
        "the reporter selects the five highest-product feasible rows retaining the configured "
        "viability fraction, with ties resolved by growth and gene ID. Display rank is not a "
        "recommendation or confidence rank. These are the canonical D1–D5 forward-validation "
        "candidates; final support remains a separate verdict. The complete unmodified screen "
        "remains in the source CSV."
    )
    note += (
        f" Showing all {len(ordered)} {_plural_noun(len(ordered), 'canonical candidate')} "
        "and no screen-only rows; "
        f"the linked CSV contains all {len(rows)} screen outcomes."
    )
    return (
        f"<h4>{_escape(label)}</h4>"
        + '<div class="table-scroll">'
        + _table(
            (
                "Display rank",
                "Gene candidate",
                "Blocked reactions",
                "Reaction equations",
                "Status",
                "Growth rate (h⁻¹)",
                "Product flux (mmol gDW⁻¹ h⁻¹)",
                "Workflow evidence",
            ),
            table_rows,
            css_class="ko-screen-table",
        )
        + "</div>"
        + note
        + "</p>"
    )


def _design_table(validated: ValidatedRun, role: str, label: str) -> str:
    rows = _rows(validated, role)
    return f"<h4>{_escape(label)}</h4>" + _table(
        (
            "Knockout set",
            "Growth rate (h⁻¹)",
            "Maximum product (mmol gDW⁻¹ h⁻¹)",
            "Guaranteed product (mmol gDW⁻¹ h⁻¹)",
            "Growth coupled",
        ),
        [
            (
                row.get("knockouts", ""),
                _format_number(row.get("growth", "")),
                _format_number(row.get("max_product", "")),
                _format_number(row.get("guaranteed_product", "")),
                row.get("growth_coupled", ""),
            )
            for row in rows[:25]
        ],
    )


def _figure_html(
    root: Path,
    figure: Mapping[str, Any],
    *,
    standalone: bool,
    display_number: int | None,
) -> str:
    if figure.get("status") != "rendered":
        return (
            '<div class="panel-unavailable"><strong>Validation panel unavailable.</strong> '
            f"{_escape(figure.get('reason', 'No reason was recorded.'))}</div>"
        )
    if display_number is None:
        raise FigureRenderError(
            f"rendered figure {figure.get('id')!r} has no display number"
        )
    outputs = figure.get("outputs", {})
    if not isinstance(outputs, Mapping):
        raise FigureRenderError(f"figure {figure.get('id')!r} has no output mapping")
    png = _relative_path(root, outputs.get("png"), label=f"{figure.get('id')} PNG")
    if standalone:
        encoded = base64.b64encode(png.read_bytes()).decode("ascii")
        source = f"data:image/png;base64,{encoded}"
    else:
        source = png.relative_to(root).as_posix()
    sources = figure.get("sources", [])
    source_text = (
        ", ".join(str(item) for item in sources) if isinstance(sources, list) else ""
    )
    output_text = ", ".join(
        str(outputs.get(kind)) for kind in ("pdf", "svg", "png") if outputs.get(kind)
    )
    width_mm = figure.get("width_mm")
    if width_mm not in {89, 180}:
        raise FigureRenderError(
            f"rendered figure {figure.get('id')!r} has invalid width_mm {width_mm!r}"
        )
    image_style = f"width:{width_mm}mm;max-width:100%"
    return (
        f'<figure id="{_escape(figure.get("id", "figure"))}">'
        f'<img src="{_escape(source)}" alt="{_escape(figure.get("alt", "CMM result figure"))}" '
        f'style="{image_style}">'
        f"<figcaption><strong>Figure {display_number}.</strong> "
        f"{_escape(figure.get('caption', ''))} "
        f"Sources: <code>{_escape(source_text)}</code>. "
        f"Files: <code>{_escape(output_text)}</code>.</figcaption></figure>"
    )


def _figure_lookup(
    figure_manifest: FigureManifest,
) -> tuple[dict[str, Mapping[str, Any]], dict[str, int]]:
    figures = {str(item.get("id")): item for item in figure_manifest.figures}
    rendered_order = [
        figure_id
        for figure_id in _FIGURE_ORDER
        if figures.get(figure_id, {}).get("status") == "rendered"
    ]
    numbers = {
        figure_id: number for number, figure_id in enumerate(rendered_order, start=1)
    }
    return figures, numbers


def _preflight_html(validated: ValidatedRun) -> tuple[str, list[str]]:
    rows = _rows(validated, "preflight_checks")
    warnings = [
        row.get("message", "")
        for row in rows
        if row.get("status", "").lower() == "warning" and row.get("message", "")
    ]
    table = _table(
        ("Check", "Status", "Value", "Interpretation"),
        [
            (
                row.get("check", ""),
                row.get("status", ""),
                row.get("value", ""),
                row.get("message", ""),
            )
            for row in rows
        ],
    )
    return f'<div class="preflight-table">{table}</div>', warnings


def _setup_table(validated: ValidatedRun) -> str:
    provenance = _json(validated, "provenance")
    summary = _json(validated, "summary")
    parameters = provenance.get("parameters", {})
    medium = provenance.get("medium_application", {})
    rows: list[tuple[object, object]] = []
    for label, value in (
        ("Model", provenance.get("model_id", summary.get("model_id", ""))),
        ("Product exchange", summary.get("product", "")),
        ("Substrate exchange", summary.get("substrate", "")),
        ("Biomass reaction", summary.get("biomass", "")),
        ("Solver", provenance.get("solver", "")),
    ):
        rows.append((label, _display_value(value)))
    if isinstance(parameters, Mapping):
        for key in ("condition", "reference_method", "requested_solver"):
            if key in parameters:
                rows.append((_human_label(key), _display_value(parameters[key])))
    if isinstance(medium, Mapping):
        for key in ("name", "mode", "applied", "dropped"):
            if key in medium:
                rows.append(
                    (f"Medium: {_human_label(key)}", _display_value(medium[key]))
                )
    return _table(("Run definition", "Resolved value"), rows)


def _amplification_support_statement(
    validated: ValidatedRun, supported_targets: Sequence[str]
) -> str:
    if supported_targets:
        return (
            "Computational support was assigned to "
            + _plural_noun(len(supported_targets), "amplification target")
            + " "
            + ", ".join(supported_targets)
            + "."
        )
    return (
        "No independently proposed FSEOF or FVSEOF target passed the workflow's "
        "method-specific forward-validation support rule."
    )


def _summary_html(validated: ValidatedRun, *, standalone: bool) -> str:
    report = validated.report
    findings_value = report.get("findings", ())
    findings = (
        [str(item) for item in findings_value if isinstance(item, str) and item.strip()]
        if isinstance(findings_value, list)
        else []
    )
    summary = _json(validated, "summary")
    yield_rows = _rows(validated, "theoretical_yield")
    yield_row = yield_rows[0] if yield_rows else {}
    recommendations = _rows(validated, "recommendations")
    supported_knockout_rows = sorted(
        (
            row
            for row in recommendations
            if row.get("type") == "single_gene_knockout"
            and row.get("verdict", "").lower() == "support"
            and row.get("target", "")
        ),
        key=lambda row: row.get("target", ""),
    )
    supported_amplifications = sorted(
        row.get("target", "")
        for row in recommendations
        if row.get("type") == "amplification"
        and row.get("verdict", "").lower() == "support"
        and row.get("target", "")
    )
    coupled_designs = [
        row
        for row in recommendations
        if row.get("type") == "multi_knockout"
        and row.get("verdict", "").lower() == "coupled"
    ]

    def product_effect(row: Mapping[str, str]) -> float:
        try:
            return float(row.get("product_effect", ""))
        except ValueError:
            return float("-inf")

    top_design = max(
        coupled_designs,
        key=product_effect,
        default=None,
    )
    findings_html = (
        "<ul>" + "".join(f"<li>{_escape(item)}</li>" for item in findings) + "</ul>"
        if findings
        else ""
    )
    knockout_details: list[str] = []
    for row in supported_knockout_rows:
        details: list[str] = []
        growth = _format_growth_retained(row.get("growth_retained", ""))
        if growth != "Not quantified":
            details.append(f"{growth} WT growth retained")
        effect = _format_signed_number(row.get("product_effect", ""))
        if effect != "—":
            details.append(
                f"conservative paired-sampling product-flux shift {effect} mmol gDW⁻¹ h⁻¹"
            )
        suffix = f" ({'; '.join(details)})" if details else ""
        knockout_details.append(f"{_escape(row.get('target', ''))}{suffix}")
    knockout_statement = (
        "Final computational support was assigned to "
        + _plural_noun(
            len(knockout_details),
            "single-gene knockout",
            "single-gene knockouts",
        )
        + " "
        + "; ".join(knockout_details)
        + "."
        if knockout_details
        else "No single-gene knockout passed the exported support rule."
    )
    amplification_statement = _amplification_support_statement(
        validated, supported_amplifications
    )
    sentences = [
        (
            f"The requested product exchange {_escape(summary.get('product', validated.report.get('product_label', '')))} "
            f"had a theoretical flux ceiling of {_escape(_format_number(yield_row.get('product_flux', '')))} "
            "mmol gDW⁻¹ h⁻¹ and a molar yield of "
            f"{_escape(_format_number(yield_row.get('molar_yield', '')))} mol mol⁻¹ "
            f"({_source_refs(validated, ('theoretical_yield',), standalone=standalone)})."
        ),
        (
            f"MOMA and ROOM evaluated {len(_rows(validated, 'single_knockout_moma'))} and "
            f"{len(_rows(validated, 'single_knockout_room'))} exported deletion outcomes, respectively; "
            "their predictions are reported separately because they encode different adaptation assumptions."
        ),
        (
            f"{knockout_statement} "
            f"{_escape(amplification_statement)} These labels mean that each promoted row "
            "passed the workflow's declared proposal and forward-validation rules, not that efficacy is causal "
            "or experimentally demonstrated."
        ),
    ]
    coverage = summary.get("validation_coverage")
    if isinstance(coverage, Mapping):
        sentences.append(
            "Forward validation covered "
            f"{_escape(coverage.get('flux_response_completed', 0))}/"
            f"{_escape(coverage.get('flux_response_expected', 0))} flux-response "
            "candidate analyses and "
            f"{_escape(coverage.get('sampling_completed', 0))}/"
            f"{_escape(coverage.get('sampling_expected', 0))} sampling executions "
            "including the shared wild-type ensemble; failed or skipped executions are "
            "listed individually in §4.5."
        )
    if top_design is not None:
        top_effect = product_effect(top_design)
        tied_top = [
            row
            for row in coupled_designs
            if math.isclose(
                product_effect(row), top_effect, rel_tol=1e-9, abs_tol=1e-12
            )
        ]
        top_growth = _format_growth_retained(top_design.get("growth_retained", ""))
        growth_clause = (
            f" ({top_growth} WT growth retained)"
            if top_growth != "Not quantified"
            else ""
        )
        tied_targets = [
            row.get("target", "")
            for row in tied_top
            if row is not top_design and row.get("target", "")
        ]
        tie_clause = (
            f" This top value is tied with {_escape(', '.join(tied_targets))}."
            if tied_targets
            else ""
        )
        sentences.append(
            f"The export contains {len(coupled_designs)} growth-coupled reaction-level "
            f"{_plural_noun(len(coupled_designs), 'multi-knockout design')}; "
            "the highest guaranteed product flux is "
            f"{_escape(_format_number(top_design.get('product_effect', '')))} mmol gDW⁻¹ h⁻¹ "
            f"for {_escape(top_design.get('target', ''))}{growth_clause}.{tie_clause} "
            "Reaction interventions still require "
            "GPR-aware resolution to experimental gene edits."
        )
    return findings_html + "".join(f"<p>{sentence}</p>" for sentence in sentences)


def _amplification_table(validated: ValidatedRun) -> str:
    diagnostic = {
        row.get("target", ""): row
        for row in _rows(validated, "amplification_loop_diagnostic")
        if row.get("target", "")
    }

    def ranked_rows(raw_role: str, tidy_role: str) -> list[dict[str, str]]:
        raw = _rows(validated, raw_role)
        tidy = _rows(validated, tidy_role)
        tidy_order = list(dict.fromkeys(row.get("target", "") for row in tidy))
        if not raw:
            raw = []
            for target in tidy_order:
                first = next(
                    (row for row in tidy if row.get("target", "") == target),
                    {},
                )
                fallback = dict(first)
                fallback["reaction_id"] = target
                fallback["report_selected"] = "true"
                raw.append(fallback)
        tidy_rank: dict[str, float] = {}
        for order, target in enumerate(tidy_order, start=1):
            matching = [row for row in tidy if row.get("target", "") == target]
            declared = []
            for row in matching:
                try:
                    declared.append(float(row.get("method_rank", "")))
                except (TypeError, ValueError):
                    continue
            tidy_rank[target] = min(declared) if declared else float(order)

        has_selection = any("report_selected" in row for row in raw)
        ranked: list[tuple[float, int, dict[str, str]]] = []
        for order, original in enumerate(raw):
            row = dict(original)
            target = row.get("reaction_id", row.get("target", ""))
            if not target:
                continue
            if has_selection:
                if not _truthy(row.get("report_selected", "")):
                    continue
            elif row.get("classification", "").lower() != "amplify" or not _truthy(
                row.get("actionable", "")
            ):
                continue
            rank = math.nan
            for column in ("amplification_rank", "method_rank"):
                try:
                    rank = float(row.get(column, ""))
                except (TypeError, ValueError):
                    continue
                if math.isfinite(rank):
                    break
            if not math.isfinite(rank):
                rank = tidy_rank.get(target, float(order + 1))
            row["target"] = target
            row["display_rank"] = _format_number(rank)
            ranked.append((rank, order, row))
        ranked.sort(key=lambda item: (item[0], item[1], item[2]["target"]))
        return [row for _, _, row in ranked[:10]]

    def loop_status(row: Mapping[str, str]) -> str:
        target = row.get("target", "")
        detail = diagnostic.get(target, {})
        flagged = _truthy(
            row.get("loop_artifact_flag", detail.get("loop_artifact_flag", ""))
        )
        if flagged:
            return "Diagnostic-only loop artifact"
        status = row.get("diagnostic_status", detail.get("diagnostic_status", ""))
        return "Cleared" if status == "complete" else status or "Not assessed"

    fseof = ranked_rows("amplification_target_ranking", "fseof_tidy")
    fvseof = ranked_rows("variability_supported_amplification_targets", "fvseof_tidy")
    fseof_table = _table(
        ("Method rank", "Reaction target", "Class", "Slope", "Loop diagnostic"),
        [
            (
                row.get("display_rank", ""),
                row.get("target", ""),
                row.get("classification", ""),
                _format_number(row.get("slope", "")),
                loop_status(row),
            )
            for row in fseof
        ],
    )
    fvseof_table = _table(
        (
            "Method rank",
            "Reaction target",
            "Class",
            "Robust",
            "Slope",
            "Mean capacity",
            "Loop diagnostic",
        ),
        [
            (
                row.get("display_rank", ""),
                row.get("target", ""),
                row.get("classification", ""),
                row.get("robust", ""),
                _format_number(row.get("slope", "")),
                _format_number(row.get("mean_capacity", "")),
                loop_status(row),
            )
            for row in fvseof
        ],
    )
    return (
        "<h4>FSEOF independent top ten</h4>"
        + fseof_table
        + "<h4>FVSEOF independent top ten</h4>"
        + fvseof_table
        + '<p class="muted">Each method contributes its own ranked hypotheses; target '
        "intersection is not required. Loop-flagged rows are displayed for diagnosis and "
        "retained in flux-response validation, but are excluded from support and recommendation "
        "eligibility.</p>"
    )


def _mapping_table(validated: ValidatedRun) -> str:
    rows = _rows(validated, "gene_knockout_mapping")
    recommendations = _rows(validated, "recommendations")
    consensus = _rows(validated, "single_knockout_consensus")
    prioritized: list[str] = []

    def add_priority(target: str) -> None:
        if target and target not in prioritized:
            prioritized.append(target)

    for row in recommendations:
        if (
            row.get("type") == "single_gene_knockout"
            and row.get("verdict", "").lower() == "support"
        ):
            add_priority(row.get("target", ""))
    for role in ("single_knockout_moma", "single_knockout_room"):
        screen_rows = _rows(validated, role)
        display_ranks = _ko_display_ranks(validated, screen_rows)
        for index in sorted(display_ranks, key=lambda value: display_ranks[value]):
            add_priority(screen_rows[index].get("target_id", ""))
    for row in consensus:
        if _truthy(row.get("recommended", "")):
            add_priority(row.get("target_id", ""))
    priority = {target: index for index, target in enumerate(prioritized)}

    def normalized_signature(reactions: Sequence[str]) -> str:
        return ";".join(sorted(value for value in reactions if value))

    mapped_reactions: dict[str, set[str]] = {}
    for row in rows:
        gene_id = row.get("gene_id", "").strip()
        blocked_reaction = row.get("blocked_reaction", "").strip()
        if gene_id and blocked_reaction:
            mapped_reactions.setdefault(gene_id, set()).add(blocked_reaction)
    gene_signatures = {
        gene: normalized_signature(tuple(reactions))
        for gene, reactions in mapped_reactions.items()
    }
    screen_signatures: dict[str, str] = {}
    for role in ("single_knockout_moma", "single_knockout_room"):
        for row in _rows(validated, role):
            target = row.get("target_id", "").strip()
            signature = normalized_signature(
                tuple(
                    value.strip()
                    for value in row.get("blocked_reaction_signature", "").split(";")
                )
            )
            if target and signature:
                screen_signatures[target] = signature

    representative_signatures: dict[str, int] = {}
    for gene, rank in priority.items():
        signature = screen_signatures.get(gene, gene_signatures.get(gene, ""))
        if signature:
            representative_signatures[signature] = min(
                representative_signatures.get(signature, len(priority)), rank
            )
    indexed_displayed = [
        (index, row)
        for index, row in enumerate(rows)
        if row.get("gene_id", "") in priority
        or gene_signatures.get(row.get("gene_id", ""), "") in representative_signatures
    ]
    indexed_displayed.sort(
        key=lambda item: (
            representative_signatures.get(
                gene_signatures.get(item[1].get("gene_id", ""), ""),
                priority.get(item[1].get("gene_id", ""), len(priority)),
            ),
            item[1].get("gene_id", "") not in priority,
            item[0],
        )
    )
    displayed = [row for _, row in indexed_displayed]

    def mapping_relation(row: Mapping[str, str]) -> str:
        gene_id = row.get("gene_id", "")
        signature = gene_signatures.get(gene_id, "")
        if signature and priority.get(
            gene_id, len(priority)
        ) > representative_signatures.get(signature, len(priority)):
            suffix = " (also D1–D5 candidate)" if gene_id in priority else ""
            return f"Signature-equivalent model deletion{suffix}"
        return "Screened validation representative"

    table = _table(
        (
            "Gene",
            "Relation to screened target",
            "Name",
            "Inert",
            "Blocked reaction",
            "Reaction name",
            "Reaction equation",
            "GPR",
        ),
        [
            (
                row.get("gene_id", ""),
                mapping_relation(row),
                row.get("gene_name", ""),
                row.get("inert", ""),
                row.get("blocked_reaction", ""),
                row.get("reaction_name", ""),
                row.get("reaction_equation", ""),
                row.get("gpr", ""),
            )
            for row in displayed
        ],
    )
    artifact = validated.artifacts.get("gene_knockout_mapping")
    path = artifact.relative_path if artifact is not None else None
    note = (
        "Shown rows include supported or canonical D1–D5 validation representatives and every gene "
        "with the same blocked-reaction signature. Signature-equivalent rows are model-equivalent "
        "deletions, not a preferred wet-lab gene; GPR-aware experimental resolution is still required. "
        f"The full gene-to-reaction mapping remains in <code>{_escape(path or 'gene_knockout_mapping unavailable')}</code>."
    )
    return table + f'<p class="muted">{note}</p>'


def _loop_table(validated: ValidatedRun) -> str:
    rows = _rows(validated, "amplification_loop_diagnostic")
    return _table(
        (
            "Target",
            "Standard capacity",
            "Loopless capacity",
            "Loopless/standard",
            "Loop artifact",
            "Status",
            "Reason",
        ),
        [
            (
                row.get("target", ""),
                _format_number(row.get("standard_capacity", "")),
                _format_number(row.get("loopless_capacity", "")),
                _format_number(row.get("loopless_to_standard_capacity_ratio", "")),
                row.get("loop_artifact_flag", ""),
                row.get("diagnostic_status", ""),
                row.get("reason", ""),
            )
            for row in rows
        ],
    )


def _validation_table(validated: ValidatedRun) -> str:
    responses = _rows(validated, "flux_response_tidy")
    index_rows = _rows(validated, "flux_response_validation_index")

    def base_key(row: Mapping[str, str]) -> tuple[str, str, str, str]:
        return (
            row.get("target", ""),
            row.get("scan_reaction", ""),
            row.get("response_reaction", ""),
            row.get("background", ""),
        )

    indexed_by_base: dict[tuple[str, str, str, str], list[Mapping[str, str]]] = {}
    for row in index_rows:
        indexed_by_base.setdefault(base_key(row), []).append(row)

    def response_scope(row: Mapping[str, str]) -> str:
        scope = row.get("candidate_scope", "").strip()
        if scope:
            return scope
        matches = indexed_by_base.get(base_key(row), ())
        scopes = {match.get("candidate_scope", "").strip() for match in matches}
        scopes.discard("")
        return next(iter(scopes)) if len(scopes) == 1 else ""

    grouped: dict[tuple[str, str, str, str, str], dict[str, int]] = {}
    for row in responses:
        key = (*base_key(row), response_scope(row))
        values = grouped.setdefault(key, {"optimal": 0, "infeasible": 0})
        status = row.get("status", "")
        bucket = "optimal" if status == "optimal" else "infeasible"
        values[bucket] += 1
    indexed = {
        (*base_key(row), row.get("candidate_scope", "").strip()): row
        for row in index_rows
    }
    recommendation_verdicts = {
        row.get("target", ""): row.get("verdict", "")
        for row in _rows(validated, "recommendations")
        if row.get("target", "")
    }
    keys = list(grouped)
    for key in indexed:
        if key not in grouped:
            keys.append(key)

    def aliases(row: Mapping[str, str] | None, target: str) -> str:
        if row is None:
            return target
        value = row.get("candidate_target_ids", "").strip()
        return value or target

    table = _table(
        (
            "Simulation target",
            "Candidate target IDs",
            "Scan reaction",
            "Response reaction",
            "Background",
            "Candidate scope",
            "Reference candidate flux (mmol gDW⁻¹ h⁻¹)",
            "Loop diagnostic eligible",
            "Execution status",
            "Recommendation verdict",
            "Optimal points",
            "Other/infeasible",
            "Reason",
        ),
        [
            (
                key[0],
                aliases(indexed.get(key), key[0]),
                key[1],
                key[2],
                key[3],
                key[4] or "legacy/not declared",
                _format_number(indexed.get(key, {}).get("scan_reference_flux", ""))
                if indexed.get(key, {}).get("scan_reference_flux", "") != ""
                else "—",
                indexed.get(key, {}).get("loop_diagnostic_eligible", "—") or "—",
                indexed.get(key, {}).get("status", "complete"),
                recommendation_verdicts.get(key[0], "Not promoted"),
                grouped.get(key, {}).get("optimal", 0),
                grouped.get(key, {}).get("infeasible", 0),
                indexed.get(key, {}).get("reason", "")
                or indexed.get(key, {}).get("error", "")
                or indexed.get(key, {}).get("loop_diagnostic_reason", ""),
            )
            for key in sorted(keys)
        ],
        css_class="validation-table",
    )
    return f'<div class="table-scroll">{table}</div>'


def _display_ranked_signature_candidates(
    validated: ValidatedRun,
) -> dict[str, tuple[str, ...]]:
    """Return the D1-D5 union grouped by model-equivalent reaction signature."""

    def normalized_signature(value: str) -> str:
        return ";".join(
            sorted(part.strip() for part in value.split(";") if part.strip())
        )

    screen_rows: list[Mapping[str, str]] = []
    ranked_signatures: set[str] = set()
    for role in ("single_knockout_moma", "single_knockout_room"):
        for row in _rows(validated, role):
            screen_rows.append(row)
            try:
                rank = int(float(row.get("display_rank", "")))
            except (TypeError, ValueError):
                continue
            if not 1 <= rank <= 5:
                continue
            signature = normalized_signature(row.get("blocked_reaction_signature", ""))
            if signature:
                ranked_signatures.add(signature)

    grouped: dict[str, set[str]] = {signature: set() for signature in ranked_signatures}
    for screen_row in screen_rows:
        signature = normalized_signature(
            screen_row.get("blocked_reaction_signature", "")
        )
        target = screen_row.get("target_id", "").strip()
        if signature in grouped and target:
            grouped[signature].add(target)

    mapped_reactions: dict[str, set[str]] = {}
    for row in _rows(validated, "gene_knockout_mapping"):
        gene = row.get("gene_id", "").strip()
        reaction = row.get("blocked_reaction", "").strip()
        if gene and reaction:
            mapped_reactions.setdefault(gene, set()).add(reaction)
    for gene, reactions in mapped_reactions.items():
        signature = normalized_signature(";".join(reactions))
        if signature in grouped:
            grouped[signature].add(gene)
    return {key: tuple(sorted(values)) for key, values in sorted(grouped.items())}


def _report_selected_amplification_targets(validated: ValidatedRun) -> set[str]:
    targets: set[str] = set()
    for role in (
        "amplification_target_ranking",
        "variability_supported_amplification_targets",
    ):
        for row in _rows(validated, role):
            if (
                _truthy(row.get("report_selected", ""))
                and row.get("reaction_id", "").strip()
            ):
                targets.add(row["reaction_id"].strip())
    return targets


def _status_counts(rows: Sequence[Mapping[str, str]]) -> dict[str, int]:
    counts = {"complete": 0, "failed": 0, "skipped": 0, "other": 0}
    for row in rows:
        status = row.get("status", "").strip().lower()
        counts[status if status in counts else "other"] += 1
    return counts


def _coverage_cell(*, expected: int, indexed: int, exhaustive: bool) -> str:
    if exhaustive:
        return f"{indexed}/{expected}"
    return f"{indexed} indexed; legacy scope not declared exhaustive"


def _validation_execution_summary(validated: ValidatedRun) -> str:
    """Describe complete and unavailable candidate-level forward checks without promotion."""

    response_index = _rows(validated, "flux_response_validation_index")
    sampling_index = _rows(validated, "single_knockout_sampling_validation_index")
    response_rows = _rows(validated, "flux_response_tidy")
    sampling_rows = _rows(validated, "sampling_tidy")
    display_signatures = _display_ranked_signature_candidates(validated)
    expected_amplification = _report_selected_amplification_targets(validated)

    response_amplification: list[Mapping[str, str]] = []
    response_knockout: list[Mapping[str, str]] = []
    for row in response_index:
        scope = row.get("candidate_scope", "")
        if scope == "all_report_selected_candidates":
            response_amplification.append(row)
        elif scope == "all_display_ranked_candidates":
            response_knockout.append(row)
        elif row.get("background") == "gene_knockout":
            response_knockout.append(row)
        else:
            response_amplification.append(row)
    sampling_knockout = [
        row for row in sampling_index if row.get("target_id") != "wild_type"
    ]
    wild_type_sampling = [
        row for row in sampling_index if row.get("target_id") == "wild_type"
    ]

    # Legacy v2 bundles may not declare index roles.  Their completed tidy targets remain
    # reportable, but the reporter must not imply that unindexed candidates were executed.
    if not response_index:
        response_amplification = [
            {"target": target, "status": "complete", "background": "wild_type"}
            for target in sorted(
                {
                    row.get("target", "")
                    for row in response_rows
                    if row.get("background") == "wild_type" and row.get("target", "")
                }
            )
        ]
        response_knockout = [
            {
                "target": target,
                "candidate_target_ids": target,
                "status": "complete",
                "background": "gene_knockout",
            }
            for target in sorted(
                {
                    row.get("target", "")
                    for row in response_rows
                    if row.get("background") == "gene_knockout"
                    and row.get("target", "")
                }
            )
        ]
    if not sampling_index:
        sampling_knockout = [
            {
                "target_id": target,
                "candidate_target_ids": target,
                "status": "complete",
            }
            for target in sorted(
                {
                    row.get("target", "")
                    for row in sampling_rows
                    if row.get("target", "")
                }
            )
        ]

    amplification_exhaustive = bool(response_amplification) and all(
        row.get("candidate_scope") == "all_report_selected_candidates"
        for row in response_amplification
    )
    response_ko_exhaustive = bool(response_knockout) and all(
        row.get("candidate_scope") == "all_display_ranked_candidates"
        for row in response_knockout
    )
    sampling_exhaustive = bool(sampling_knockout) and all(
        row.get("candidate_scope") == "all_display_ranked_candidates"
        for row in sampling_knockout
    )

    coverage_rows: list[tuple[object, ...]] = []
    for label, rows, expected, exhaustive in (
        (
            "Flux response — amplification",
            response_amplification,
            len(expected_amplification),
            amplification_exhaustive,
        ),
        (
            "Flux response — single-gene knockout signatures",
            response_knockout,
            len(display_signatures),
            response_ko_exhaustive,
        ),
        (
            "Random sampling — single-gene knockout signatures",
            sampling_knockout,
            len(display_signatures),
            sampling_exhaustive,
        ),
    ):
        counts = _status_counts(rows)
        coverage_rows.append(
            (
                label,
                _coverage_cell(
                    expected=expected, indexed=len(rows), exhaustive=exhaustive
                ),
                counts["complete"],
                counts["failed"],
                counts["skipped"],
            )
        )

    exceptions: list[tuple[object, ...]] = []
    for assay, target_column, rows in (
        ("Flux response", "target", response_index),
        ("Random sampling", "target_id", sampling_index),
    ):
        for row in rows:
            if row.get("status", "").lower() == "complete":
                continue
            exceptions.append(
                (
                    assay,
                    row.get(target_column, ""),
                    row.get("candidate_target_ids", "") or row.get(target_column, ""),
                    row.get("status", ""),
                    row.get("reason", "") or row.get("error", ""),
                )
            )
    for role, label in (
        ("flux_response_tidy", "Flux response"),
        ("sampling_tidy", "Random sampling"),
    ):
        artifact = validated.artifacts.get(role)
        if artifact is not None and artifact.status in {"skipped", "failed"}:
            exceptions.append(
                (
                    label,
                    "All",
                    "—",
                    artifact.status,
                    artifact.reason or "No reason recorded",
                )
            )

    aliases = sorted(
        {
            (
                row.get("target", row.get("target_id", "")),
                row.get("candidate_target_ids", ""),
                row.get("blocked_reaction_signature", ""),
            )
            for row in (*response_knockout, *sampling_knockout)
            if ";" in row.get("candidate_target_ids", "")
        }
    )
    alias_note = (
        '<p class="muted"><strong>Signature-equivalent candidates.</strong> '
        + "; ".join(
            f"{_escape(representative)} represents {_escape(candidate_ids)} "
            f"({_escape(signature or 'shared blocked-reaction signature')})"
            for representative, candidate_ids, signature in aliases
        )
        + ". These aliases share one model intervention and therefore one forward-validation execution; none is treated as unvalidated.</p>"
        if aliases
        else ""
    )
    legacy_note = (
        '<div class="notice"><strong>Legacy schema-v2 coverage.</strong> One or more validation '
        "indexes do not declare the exhaustive candidate-scope marker. Counts below describe "
        "only indexed or tidy completed analyses and do not imply that every D1-D5 or "
        "report-selected candidate was executed.</div>"
        if not (
            amplification_exhaustive and response_ko_exhaustive and sampling_exhaustive
        )
        else ""
    )
    exception_html = (
        "<h4>Failed or skipped candidate analyses</h4>"
        + _table(
            (
                "Assay",
                "Simulation target",
                "Candidate target IDs",
                "Status",
                "Reason",
            ),
            exceptions,
        )
        if exceptions
        else '<p class="muted">No candidate-level flux-response or knockout-sampling execution was failed or skipped.</p>'
    )
    wild_type_note = (
        f" A shared wild-type sampling execution is also indexed ({len(wild_type_sampling)} row)."
        if sampling_index
        else ""
    )
    response_totals = _status_counts([*response_amplification, *response_knockout])
    sampling_totals = _status_counts(
        [*wild_type_sampling, *sampling_knockout]
        if sampling_index
        else sampling_knockout
    )
    exploratory_knockouts: list[str] = []
    for response_row in response_knockout:
        if response_row.get("status", "").strip().lower() != "complete":
            continue
        try:
            reference_flux = float(response_row.get("scan_reference_flux", ""))
        except (TypeError, ValueError):
            continue
        if math.isfinite(reference_flux) and reference_flux == 0:
            exploratory_knockouts.append(response_row.get("target", ""))
    exploratory_note = (
        '<p class="muted"><strong>Zero-reference knockout candidates.</strong> '
        + _escape(", ".join(sorted(set(exploratory_knockouts))))
        + " were scanned across the full feasible candidate-reaction domain. These scans are "
        "exploratory and cannot causally support deletion.</p>"
        if exploratory_knockouts
        else ""
    )
    totals_html = (
        "<p><strong>Execution totals.</strong> Flux response indexed "
        f"{len(response_amplification) + len(response_knockout)} "
        f"{_plural_noun(len(response_amplification) + len(response_knockout), 'candidate execution')} "
        f"({response_totals['complete']} complete, {response_totals['failed']} failed, "
        f"{response_totals['skipped']} skipped). Random sampling indexed "
        f"{len(wild_type_sampling) + len(sampling_knockout) if sampling_index else len(sampling_knockout)} "
        f"{_plural_noun(len(wild_type_sampling) + len(sampling_knockout) if sampling_index else len(sampling_knockout), 'execution')} "
        "including the shared wild type when declared "
        f"({sampling_totals['complete']} complete, {sampling_totals['failed']} failed, "
        f"{sampling_totals['skipped']} skipped).</p>"
    )
    return (
        "<h4>Forward-validation execution coverage</h4>"
        + totals_html
        + _table(
            (
                "Assay and candidate scope",
                "Indexed/expected",
                "Complete",
                "Failed",
                "Skipped",
            ),
            coverage_rows,
        )
        + f'<p class="muted">Candidate coverage is separate from recommendation status.{wild_type_note} '
        "Every completed tidy target is shown in the corresponding figure; raw scan and sampling rows remain in the source CSVs.</p>"
        + exploratory_note
        + alias_note
        + legacy_note
        + exception_html
    )


def _sampling_validation_table(validated: ValidatedRun) -> str:
    index_rows = [
        row
        for row in _rows(validated, "single_knockout_sampling_validation_index")
        if row.get("target_id") != "wild_type"
    ]
    samples = _rows(validated, "sampling_tidy")
    sample_counts: dict[str, dict[str, int]] = {}
    for row in samples:
        target = row.get("target", "")
        condition = row.get("condition", "")
        if target and condition in {"wild_type", "knockout"}:
            sample_counts.setdefault(target, {"wild_type": 0, "knockout": 0})[
                condition
            ] += 1
    if not index_rows:
        index_rows = [
            {"target_id": target, "candidate_target_ids": target, "status": "complete"}
            for target in sorted(sample_counts)
        ]
    table = _table(
        (
            "Simulation target",
            "Candidate target IDs",
            "Blocked-reaction signature",
            "Status",
            "WT tidy rows",
            "Knockout tidy rows",
            "Reason",
        ),
        [
            (
                row.get("target_id", ""),
                row.get("candidate_target_ids", "") or row.get("target_id", ""),
                row.get("blocked_reaction_signature", ""),
                row.get("status", ""),
                sample_counts.get(row.get("target_id", ""), {}).get("wild_type", 0),
                sample_counts.get(row.get("target_id", ""), {}).get("knockout", 0),
                row.get("reason", "") or row.get("error", ""),
            )
            for row in index_rows
        ],
        css_class="sampling-validation-table",
    )
    return f'<div class="table-scroll">{table}</div>'


def _recommendations_html(validated: ValidatedRun) -> str:
    rows = _rows(validated, "recommendations")
    if not rows:
        return (
            '<div class="notice">No recommendation row was exported. The inverse-method outputs '
            "remain hypotheses and are not promoted to interventions.</div>"
        )
    sections: list[str] = []
    labels = (
        (
            "single_gene_knockout",
            "Single-gene knockouts",
            "Conservative sampling Δ product flux (mmol gDW⁻¹ h⁻¹)",
            "Product effect is the minimum of the positive paired-sampling mean and median knockout-minus-WT product-flux shifts.",
        ),
        (
            "multi_knockout",
            "Reaction-level multi-knockout strain designs",
            "Guaranteed product flux (mmol gDW⁻¹ h⁻¹)",
            "Product effect is the RobustKnock guaranteed product flux at the reported growth optimum; reaction sets require GPR-aware gene resolution.",
        ),
        (
            "amplification",
            "Amplification targets",
            "Response Δ product flux (mmol gDW⁻¹ h⁻¹)",
            "Product effect is the optimized product-flux change across the exported target-flux response scan.",
        ),
    )
    for intervention_type, label, product_header, definition in labels:
        selected = [row for row in rows if row.get("type") == intervention_type]
        sections.append(f"<h3>{_escape(label)}</h3>")
        if intervention_type == "amplification" and not selected:
            sections.append(
                f'<div class="notice">{_escape(_amplification_support_statement(validated, ()))}</div>'
            )
        sections.append(
            _table(
                (
                    "Target/intervention",
                    "Verdict",
                    "Proposal methods",
                    "Validation methods",
                    "WT growth retained (%)",
                    product_header,
                    "Loop/artifact flag",
                    "Reason",
                ),
                [
                    (
                        row.get("target", ""),
                        row.get("verdict", "unavailable") or "unavailable",
                        row.get("proposal_methods", row.get("evidence", "")),
                        row.get("validation_methods", ""),
                        _format_growth_retained(row.get("growth_retained", "")),
                        _format_number(row.get("product_effect", "")),
                        row.get("artifact_flag", ""),
                        row.get("reason", ""),
                    )
                    for row in selected
                ],
            )
        )
        sections.append(
            f'<p class="muted"><strong>Definition.</strong> {_escape(definition)} '
            "WT growth retained is the predicted growth-rate fraction relative to the wild type, shown as a percentage.</p>"
        )
    return "".join(sections)


_CSS = """
:root{color-scheme:light;--ink:#171717;--muted:#626262;--rule:#d7d7d7;--accent:#1769aa}
*{box-sizing:border-box}body{margin:0 auto;padding:36px 28px 64px;max-width:1080px;color:var(--ink);
font:15px/1.55 Arial,Helvetica,sans-serif;background:#fff}header{border-bottom:2px solid var(--ink);
padding-bottom:16px;margin-bottom:28px}h1{font-size:28px;line-height:1.15;margin:0 0 8px}header p{margin:0;
color:var(--muted)}h2{font-size:20px;margin:38px 0 12px;border-bottom:1px solid var(--rule);
padding-bottom:5px}h3{font-size:17px;margin:26px 0 10px}h4{font-size:14px;margin:20px 0 8px}
p{max-width:88ch}table{width:100%;border-collapse:collapse;margin:12px 0 18px;font-size:12.5px;
line-height:1.35}th,td{padding:7px 8px;border-bottom:1px solid var(--rule);text-align:left;
vertical-align:top;overflow-wrap:anywhere}th{font-weight:700;border-bottom:1.5px solid var(--ink)}code{
font:11.5px/1.4 ui-monospace,SFMono-Regular,Menlo,monospace;background:#f5f5f5;padding:1px 3px}
a{color:var(--accent)}.preflight-table th:nth-child(2),.preflight-table td:nth-child(2){min-width:68px;
white-space:nowrap;overflow-wrap:normal}.table-scroll{max-width:100%;overflow-x:auto}.ko-screen-table{
table-layout:fixed;min-width:980px}.ko-screen-table th,.ko-screen-table td{overflow-wrap:normal;
word-break:normal;hyphens:none}.ko-screen-table th:nth-child(1){width:6%}.ko-screen-table th:nth-child(2){width:10%}
.ko-screen-table th:nth-child(3){width:20%}.ko-screen-table th:nth-child(4){width:28%}.ko-screen-table th:nth-child(5){width:7%}
.ko-screen-table th:nth-child(6){width:9%}.ko-screen-table th:nth-child(7){width:11%}.ko-screen-table th:nth-child(8){width:9%}
.ko-screen-table td:nth-child(1),.ko-screen-table td:nth-child(5),.ko-screen-table td:nth-child(6),
.ko-screen-table td:nth-child(7){white-space:nowrap}.validation-table{min-width:1180px}.sampling-validation-table{
min-width:900px}.validation-table th,.validation-table td,.sampling-validation-table th,
.sampling-validation-table td{overflow-wrap:normal;word-break:normal;hyphens:none}.validation-table td:nth-child(1),
.validation-table td:nth-child(3),.validation-table td:nth-child(4),.validation-table td:nth-child(7),
.validation-table td:nth-child(8),.sampling-validation-table td:nth-child(1),
.sampling-validation-table td:nth-child(4){white-space:nowrap}figure{margin:22px 0 30px}figure img{display:block;max-width:100%;height:auto;
margin:0 auto;border:0}figcaption{margin-top:8px;color:#333;font-size:12px;line-height:1.45}.muted{
color:var(--muted)}.notice,.panel-unavailable{border-left:3px solid #888;background:#f6f6f3;padding:10px 13px;
margin:14px 0}.panel-unavailable{color:var(--muted)}ul{padding-left:22px}@media(max-width:700px){body{
padding:22px 15px}table{display:block;overflow-x:auto;white-space:normal}}@media print{body{max-width:none;
padding:0;font-size:10pt}h2,figure{break-inside:avoid}a{color:inherit;text-decoration:none}}
""".strip()


def _document(
    validated: ValidatedRun,
    figure_manifest: FigureManifest,
    *,
    standalone: bool,
) -> str:
    report = validated.report
    title = str(report.get("title", "CMM production-engineering report"))
    subtitle = str(report.get("subtitle", "Schema-v2 reproducible analysis bundle"))
    figures, numbers = _figure_lookup(figure_manifest)
    summary_data = _json(validated, "summary")
    recommendation_rows = _rows(validated, "recommendations")
    supported_single = sorted(
        row.get("target", "")
        for row in recommendation_rows
        if row.get("type") == "single_gene_knockout"
        and row.get("verdict", "").lower() == "support"
        and row.get("target", "")
    )
    supported_amplification = sorted(
        row.get("target", "")
        for row in recommendation_rows
        if row.get("type") == "amplification"
        and row.get("verdict", "").lower() == "support"
        and row.get("target", "")
    )
    coupled_multi = [
        row
        for row in recommendation_rows
        if row.get("type") == "multi_knockout"
        and row.get("verdict", "").lower() == "coupled"
    ]

    def numeric_effect(row: Mapping[str, str]) -> float:
        try:
            return float(row.get("product_effect", ""))
        except ValueError:
            return float("-inf")

    top_multi = max(coupled_multi, key=numeric_effect, default=None)

    def available_roles(*roles: str) -> tuple[str, ...]:
        return tuple(
            role
            for role in roles
            if (artifact := validated.artifacts.get(role)) is not None
            and artifact.relative_path is not None
            and artifact.status in {"complete", "partial"}
        )

    single_sources = available_roles(
        "single_knockout_moma",
        "single_knockout_room",
        "single_knockout_consensus",
        "recommendations",
        "summary",
    )

    def figure(figure_id: str) -> str:
        value = figures.get(figure_id)
        if value is None:
            return '<div class="panel-unavailable">Figure was not declared.</div>'
        return _figure_html(
            validated.root,
            value,
            standalone=standalone,
            display_number=numbers.get(figure_id),
        )

    preflight, condition_warnings = _preflight_html(validated)
    summary_warnings_value = summary_data.get("warnings", ())
    summary_warnings = (
        [
            str(item)
            for item in summary_warnings_value
            if isinstance(item, str) and item.strip()
        ]
        if isinstance(summary_warnings_value, list)
        else []
    )
    warning_items = _deduplicate_warnings(
        [*validated.warnings, *summary_warnings, *condition_warnings]
    )
    warning_html = (
        '<div class="notice"><strong>Condition and validation warnings.</strong><ul>'
        + "".join(f"<li>{_escape(item)}</li>" for item in warning_items)
        + "</ul></div>"
        if warning_items
        else ""
    )
    # The reproduction configuration is deliberately relocatable; the original invocation
    # configuration can contain workstation-specific model and output paths.
    configuration = _json(validated, "reproduction_config") or _json(
        validated, "workflow_configuration"
    )
    parameter_rows = _flatten_scalars(configuration)
    loop_rows = _rows(validated, "amplification_loop_diagnostic")
    flagged = [
        row.get("target", "")
        for row in loop_rows
        if str(row.get("loop_artifact_flag", "")).lower() == "true"
    ]
    unavailable = [
        f"{_human_label(role)}: {artifact.reason}"
        for role, artifact in validated.artifacts.items()
        if artifact.status in {"skipped", "failed"}
    ]
    dynamic_limitations = "".join(f"<li>{_escape(item)}</li>" for item in unavailable)
    if flagged:
        dynamic_limitations += (
            "<li>Loopless diagnostics flagged the following amplification candidates as potential "
            f"loop artifacts; they are not presented as supported targets: {_escape(', '.join(flagged))}.</li>"
        )

    references_html = (
        "<ol>"
        + "".join(
            f"<li><strong>{_escape(method)}.</strong> {_escape(title)}. "
            f"<em>{_escape(venue)}</em>. "
            f'<a href="https://doi.org/{_escape(doi)}">{_escape(doi)}</a>.</li>'
            for method, title, venue, doi in _REFERENCES
        )
        + "</ol>"
    )
    supported_single_text = ", ".join(supported_single) if supported_single else "none"
    amplification_result = _amplification_support_statement(
        validated, supported_amplification
    )
    multi_result = (
        f"{len(coupled_multi)} coupled reaction-level "
        f"{_plural_noun(len(coupled_multi), 'design')} "
        f"{'was' if len(coupled_multi) == 1 else 'were'} exported; the highest "
        f"guaranteed-product design was {_escape(top_multi.get('target', ''))} at "
        f"{_escape(_format_number(top_multi.get('product_effect', '')))} mmol gDW⁻¹ h⁻¹. "
        if top_multi is not None
        else "No coupled multi-knockout recommendation was exported. "
    )

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{_escape(title)}</title><style>{_CSS}</style></head><body>
<header><h1>{_escape(title)}</h1><p>{_escape(subtitle)}</p></header><main>
<section id="section-1"><h2>1. Summary</h2>{_summary_html(validated, standalone=standalone)}{warning_html}
<div class="notice"><strong>Interpretation boundary.</strong> All interventions in this report are <em>in silico</em>
hypotheses. A computational support verdict is not evidence of wet-lab efficacy.</div></section>
<section id="section-2"><h2>2. Setup</h2>{_setup_table(validated)}
<h3>Preflight checks</h3>{preflight}<p>Source: {_source_refs(validated, ("preflight_checks",), standalone=standalone)}.</p></section>
<section id="section-3"><h2>3. Data and methods</h2>
<p>Fluxes are reported in mmol gDW⁻¹ h⁻¹, growth in h⁻¹, and molar yield in mol mol⁻¹. MOMA and ROOM
predict deletion phenotypes from the same reference state; OptKnock and RobustKnock search multi-reaction
growth-coupled designs; FSEOF and FVSEOF nominate flux-direction hypotheses. Flux response and paired
sampling are forward checks and do not turn a model prediction into an experimental claim.</p>
<h3>Resolved workflow parameters</h3>{_table(("Parameter", "Value"), parameter_rows)}</section>
<section id="section-4"><h2>4. Results</h2>
<h3>4.1 Production ceiling and growth trade-off</h3>
<p>The yield and feasible growth/product envelope establish the model-specific production boundary. Sources:
{_source_refs(validated, ("theoretical_yield", "production_envelope"), standalone=standalone)}.</p>{figure("fig01_yield_envelope")}
<h3>4.2 Single-knockout phenotypes</h3>
<p>MOMA and ROOM results remain method-specific; infeasible rows are retained in the source CSVs. The final
supported single-gene set is {_escape(supported_single_text)}. Every D1–D5 display-ranked signature is a
forward-validation candidate, while final support remains a distinct category. Sources:
{_source_refs(validated, single_sources, standalone=standalone)}.</p>
{figure("fig02_single_knockout")}{_screen_table(validated, "single_knockout_moma", "MOMA")}
{_screen_table(validated, "single_knockout_room", "ROOM")}<h4>Gene-to-reaction interpretation</h4>
{_mapping_table(validated)}
<h3>4.3 Multi-knockout strain designs</h3>
<p>{multi_result}Guaranteed product, rather than the cooperative maximum, determines the coupling interpretation.
These are reaction-level interventions and require GPR-aware resolution before any gene-edit proposal. Sources:
{_source_refs(validated, ("optknock", "robustknock"), standalone=standalone)}.</p>{figure("fig03_strain_design")}
{_design_table(validated, "optknock", "OptKnock")}{_design_table(validated, "robustknock", "RobustKnock")}
<h3>4.4 Amplification hypotheses</h3>
<p>FSEOF and FVSEOF independently contribute their method-specific top ten; intersection is not required, and
the trajectories are not gene-expression fold-change prescriptions. Loop-flagged ranks are isolated on a
diagnostic scale and retained in flux-response validation, but excluded from support and recommendation
eligibility. {_escape(amplification_result)} Sources:
{_source_refs(validated, ("amplification_target_ranking", "variability_supported_amplification_targets", "fseof_tidy", "fvseof_tidy"), standalone=standalone)}.</p>
{figure("fig04_amplification")}{_amplification_table(validated)}<h4>Loopless capacity diagnostic</h4>{_loop_table(validated)}
<h3>4.5 Forward validation</h3>
<p>Flux response is a conditional optimization under the stated background, not independent causal evidence.
Figure 5 uses the standard candidate-reaction-to-product flux-response definition in every facet: enforced candidate-reaction
flux (<code>target_flux</code>) is on the x-axis and optimized target-product flux
(<code>response_flux</code>) is on the y-axis. Biomass flux is a secondary value recorded under the configured
minimum-growth constraint, not a plot axis. The knockout-derived block contains wild-type pre-deletion
single-reaction titrations identified by the index candidate scope; a multi-reaction knockout signature remains
an explicit skipped or unavailable index row because it has no single candidate-reaction scan. Legacy schema-v2
product-to-growth rows remain auditable in their source CSV and are not relabelled as product responses. A
single-reaction candidate with zero reference flux is scanned over its full feasible domain and labelled
exploratory; that response cannot causally support deletion. Paired sampling compares feasible-state ensembles
and is not biological replication.
Every completed candidate analysis is displayed independently of whether it became a recommendation. The
sampling figure is restricted to the product exchange and biomass reaction; all other sampled reactions remain
in the source CSV. Sources:
{_source_refs(validated, available_roles("flux_response_validation_index", "flux_response_tidy", "single_knockout_sampling_validation_index", "sampling_tidy"), standalone=standalone)}.</p>
{_validation_execution_summary(validated)}
{figure("fig05_flux_response")}<h4>Flux-response candidate accounting</h4>{_validation_table(validated)}
{figure("fig06_sampling_shift")}<h4>Single-knockout sampling accounting</h4>{_sampling_validation_table(validated)}</section>
<section id="section-5"><h2>5. Recommended targets and strain proposal</h2>
<p>Rows are kept by intervention class. “Support”, “contradict”, “inconclusive”, and “unavailable” describe
agreement among the exported computations only; method count is not a confidence score.</p>{_recommendations_html(validated)}
<p>Source: {_source_refs(validated, ("recommendations",), standalone=standalone)}.</p></section>
<section id="section-6"><h2>6. Limitations</h2><ul>
<li>Predictions require experimental validation under the reported medium, aeration, substrate and strain background.</li>
<li>MOMA minimal adjustment, ROOM regulatory switching and bilevel growth optimisation are distinct assumptions;
agreement is model robustness, not replication.</li><li>Flux-space samples are correlated feasible states, not biological
variance or confidence intervals.</li><li>Reaction interventions require GPR-aware gene mapping; isozymes and
multi-subunit enzymes can change the experimental deletion set.</li><li>Artificial bounds, thermodynamically infeasible
cycles and alternative optima can inflate a response unless the exported diagnostics exclude them.</li>{dynamic_limitations}</ul></section>
<section id="section-7"><h2>7. References</h2>{references_html}</section>
<section id="section-8"><h2>8. Provenance</h2>{_provenance_table(validated)}
<p>Run directory: <code>.</code>. Authoritative manifest:
<code>{_escape(validated.manifest_path.name)}</code>.</p><h3>Reproduction assets</h3>
<p>{_source_refs(validated, ("reproduction_config", "reproduce_script", "render_script", "validate_script"), standalone=standalone)}</p>
<h3>Artifact inventory</h3>{_artifact_table(validated, standalone=standalone)}</section>
</main></body></html>
"""


class _HTMLAudit(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.image_sources: list[str] = []
        self.image_styles: list[str] = []
        self.hrefs: list[str] = []
        self.figure_count = 0
        self.figcaption_count = 0
        self.ids: list[str] = []
        self._in_caption = False
        self._caption_text: list[str] = []
        self.caption_texts: list[str] = []
        self.absolute_local_paths: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        if values.get("id") is not None:
            self.ids.append(str(values["id"]))
        if tag == "img" and values.get("src") is not None:
            self.image_sources.append(str(values["src"]))
            self.image_styles.append(str(values.get("style", "")))
        elif tag == "a" and values.get("href") is not None:
            self.hrefs.append(str(values["href"]))
        elif tag == "figure":
            self.figure_count += 1
        elif tag == "figcaption":
            self.figcaption_count += 1
            self._in_caption = True
            self._caption_text = []

    def handle_endtag(self, tag: str) -> None:
        if tag == "figcaption" and self._in_caption:
            self.caption_texts.append(" ".join(self._caption_text))
            self._in_caption = False

    def handle_data(self, data: str) -> None:
        for expression in (
            r"(?:^|[\s(])(/(?!/)[^\s<]+)",
            r"(?:^|[\s(])([A-Za-z]:[\\/][^\s<]+)",
        ):
            self.absolute_local_paths.extend(
                match.group(1) for match in re.finditer(expression, data)
            )
        if self._in_caption and data.strip():
            self._caption_text.append(data.strip())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _png_properties(path: Path) -> tuple[int, int, float | None, float | None]:
    payload = path.read_bytes()
    if payload[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("invalid PNG signature")
    offset = 8
    width = height = 0
    dpi_x: float | None = None
    dpi_y: float | None = None
    while offset + 12 <= len(payload):
        length = struct.unpack(">I", payload[offset : offset + 4])[0]
        chunk_type = payload[offset + 4 : offset + 8]
        data = payload[offset + 8 : offset + 8 + length]
        if chunk_type == b"IHDR" and len(data) >= 8:
            width, height = struct.unpack(">II", data[:8])
        elif chunk_type == b"pHYs" and len(data) == 9 and data[8] == 1:
            pixels_x, pixels_y = struct.unpack(">II", data[:8])
            dpi_x = pixels_x * 0.0254
            dpi_y = pixels_y * 0.0254
        if chunk_type == b"IEND":
            break
        offset += 12 + length
    if width <= 0 or height <= 0:
        raise ValueError("PNG has no valid IHDR dimensions")
    return width, height, dpi_x, dpi_y


def _svg_number(value: str | None) -> float | None:
    if value is None:
        return None
    matched = re.fullmatch(
        r"\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)\s*(?:px|pt)?\s*",
        value,
    )
    return float(matched.group(1)) if matched else None


def _svg_horizontal_text_overflow(svg: str) -> tuple[str, ...]:
    """Return labels whose declared svglite bounds exceed the SVG viewBox."""

    try:
        root = ET.fromstring(svg)
    except ET.ParseError as exc:
        raise ValueError(f"invalid XML: {exc}") from exc
    view_box = root.get("viewBox", "").split()
    if len(view_box) != 4:
        raise ValueError("SVG root has no four-value viewBox")
    try:
        view_x, _, view_width, _ = (float(value) for value in view_box)
    except ValueError as exc:
        raise ValueError("SVG viewBox is not numeric") from exc
    if view_width <= 0:
        raise ValueError("SVG viewBox width is not positive")
    lower_limit = view_x
    upper_limit = view_x + view_width
    tolerance = 0.75
    clipped: list[str] = []
    for element in root.iter():
        if element.tag.rsplit("}", 1)[-1] != "text":
            continue
        transform = element.get("transform", "")
        if "rotate" in transform.lower():
            continue
        x = _svg_number(element.get("x"))
        text_length = _svg_number(element.get("textLength"))
        if x is None or text_length is None or text_length < 0:
            continue
        anchor = element.get("text-anchor", "start").lower()
        if anchor == "middle":
            left, right = x - text_length / 2, x + text_length / 2
        elif anchor == "end":
            left, right = x - text_length, x
        else:
            left, right = x, x + text_length
        if left < lower_limit - tolerance or right > upper_limit + tolerance:
            label = "".join(element.itertext()).strip() or "<empty>"
            clipped.append(label[:80])
    return tuple(clipped)


def _artwork_checks(
    validated: ValidatedRun,
    figure_manifest: FigureManifest,
    *,
    issues: list[str],
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    specification = figure_manifest.renderer.get("specification", {})
    declared_script_sha256 = figure_manifest.renderer.get("script_sha256")
    expected_script_sha256 = _sha256(renderer_script_path())
    if declared_script_sha256 != expected_script_sha256:
        issues.append(
            "figure renderer script SHA-256 does not match the checked-in renderer"
        )
    declared_dpi = (
        specification.get("raster_dpi") if isinstance(specification, Mapping) else None
    )
    if declared_dpi != 300:
        issues.append(
            f"figure renderer must declare raster_dpi 300, got {declared_dpi!r}"
        )
    for figure in figure_manifest.figures:
        if figure.get("status") != "rendered":
            continue
        figure_id = str(figure.get("id"))
        width_mm = figure.get("width_mm")
        height_mm = figure.get("height_mm")
        dpi = figure.get("dpi")
        if width_mm not in {89, 180}:
            issues.append(
                f"figure {figure_id!r} width_mm must be 89 or 180, got {width_mm!r}"
            )
        if not isinstance(height_mm, (int, float)) or height_mm <= 0:
            issues.append(f"figure {figure_id!r} must declare a positive height_mm")
        if dpi != 300:
            issues.append(f"figure {figure_id!r} dpi must be 300, got {dpi!r}")
        outputs = figure.get("outputs", {})
        if not isinstance(outputs, Mapping):
            issues.append(f"figure {figure_id!r} has no output map")
            continue
        for output_type in ("png", "pdf", "svg"):
            try:
                output = _relative_path(
                    validated.root,
                    outputs.get(output_type),
                    label=f"{figure_id}.{output_type}",
                )
            except FigureRenderError as exc:
                issues.append(str(exc))
                continue
            if not output.is_file() or output.stat().st_size == 0:
                issues.append(
                    f"figure {figure_id!r} has no non-empty {output_type} output"
                )
                continue
            record: dict[str, object] = {
                "figure_id": figure_id,
                "format": output_type,
                "path": output.relative_to(validated.root).as_posix(),
                "sha256": _sha256(output),
                "size_bytes": output.stat().st_size,
            }
            if output_type == "png":
                try:
                    width_px, height_px, dpi_x, dpi_y = _png_properties(output)
                except ValueError as exc:
                    issues.append(f"figure {figure_id!r} PNG failed validation: {exc}")
                else:
                    record.update(
                        {
                            "width_px": width_px,
                            "height_px": height_px,
                            "dpi_x": dpi_x,
                            "dpi_y": dpi_y,
                        }
                    )
                    if isinstance(width_mm, (int, float)):
                        expected_width = round(float(width_mm) / 25.4 * 300)
                        if abs(width_px - expected_width) > 2:
                            issues.append(
                                f"figure {figure_id!r} PNG width is {width_px}px, expected "
                                f"{expected_width}px for {width_mm} mm at 300 DPI"
                            )
                    if dpi_x is None or dpi_y is None:
                        issues.append(
                            f"figure {figure_id!r} PNG has no physical-DPI metadata"
                        )
                    elif abs(dpi_x - 300) > 1 or abs(dpi_y - 300) > 1:
                        issues.append(
                            f"figure {figure_id!r} PNG DPI is {dpi_x:.2f}x{dpi_y:.2f}, expected 300"
                        )
            elif output_type == "pdf" and not output.read_bytes().startswith(b"%PDF-"):
                issues.append(f"figure {figure_id!r} PDF has an invalid signature")
            elif output_type == "svg":
                try:
                    svg = output.read_text(encoding="utf-8")
                except UnicodeDecodeError:
                    issues.append(f"figure {figure_id!r} SVG is not UTF-8 text")
                else:
                    if "<svg" not in svg:
                        issues.append(f"figure {figure_id!r} SVG has no svg root")
                    try:
                        clipped_text = _svg_horizontal_text_overflow(svg)
                    except ValueError as exc:
                        issues.append(
                            f"figure {figure_id!r} SVG failed structural validation: {exc}"
                        )
                    else:
                        if clipped_text:
                            issues.append(
                                f"figure {figure_id!r} SVG has horizontally clipped text: "
                                f"{', '.join(clipped_text)}"
                            )
                    if re.search(r"<(?:[A-Za-z0-9_-]+:)?image\b", svg, re.IGNORECASE):
                        issues.append(f"figure {figure_id!r} SVG embeds a raster image")
                    if re.search(
                        r"(?:data:image|href=[\"'][^\"']+\.(?:png|jpe?g))",
                        svg,
                        re.IGNORECASE,
                    ):
                        issues.append(
                            f"figure {figure_id!r} SVG references raster artwork"
                        )
            records.append(record)
    return records


def _audit_html(
    path: Path,
    *,
    validated: ValidatedRun,
    rendered: Sequence[Mapping[str, Any]],
    standalone: bool,
    issues: list[str],
) -> Mapping[str, object]:
    try:
        document = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        issues.append(f"report {path.name!r} is not readable UTF-8 HTML: {exc}")
        return {"path": path.name, "status": "failed"}
    parser = _HTMLAudit()
    parser.feed(document)
    local_root = str(validated.root)
    if local_root and local_root in document:
        issues.append(f"report {path.name!r} exposes the absolute local run directory")
    if "file://" in document.lower():
        issues.append(f"report {path.name!r} contains a local file URI")
    if parser.absolute_local_paths:
        issues.append(
            f"report {path.name!r} contains absolute local path text: "
            f"{sorted(set(parser.absolute_local_paths))}"
        )
    expected_images = len(rendered)
    if len(parser.image_sources) != expected_images:
        issues.append(
            f"report {path.name!r} has {len(parser.image_sources)} images; expected {expected_images}"
        )
    expected_styles = [
        f"width:{figure.get('width_mm')}mm;max-width:100%" for figure in rendered
    ]
    if parser.image_styles != expected_styles:
        issues.append(
            f"report {path.name!r} image widths do not match the figure manifest"
        )
    if (
        parser.figure_count != expected_images
        or parser.figcaption_count != expected_images
    ):
        issues.append(
            f"report {path.name!r} figure/caption counts do not match rendered artwork"
        )
    numbers: list[int] = []
    for caption in parser.caption_texts:
        matched = re.search(r"Figure\s+(\d+)\.", caption)
        if matched:
            numbers.append(int(matched.group(1)))
    if numbers != list(range(1, expected_images + 1)):
        issues.append(
            f"report {path.name!r} figure numbers are not continuous: {numbers}"
        )
    if len(parser.ids) != len(set(parser.ids)):
        issues.append(f"report {path.name!r} contains duplicate HTML ids")
    expected_sections = {f"section-{number}" for number in range(1, 9)}
    if not expected_sections.issubset(parser.ids):
        issues.append(
            f"report {path.name!r} is missing section id(s): {sorted(expected_sections - set(parser.ids))}"
        )

    if standalone:
        for source in parser.image_sources:
            if not source.startswith("data:image/png;base64,"):
                issues.append(
                    f"standalone report contains a non-embedded image source: {source}"
                )
                continue
            try:
                decoded = base64.b64decode(source.split(",", 1)[1], validate=True)
            except (ValueError, TypeError):
                issues.append("standalone report contains invalid base64 image data")
            else:
                if not decoded.startswith(b"\x89PNG\r\n\x1a\n"):
                    issues.append("standalone report embeds a non-PNG image")
        for href in parser.hrefs:
            parsed = urlparse(href)
            if not parsed.scheme and not href.startswith("#"):
                issues.append(f"standalone report contains a relative href: {href}")
    else:
        if any(source.startswith("data:") for source in parser.image_sources):
            issues.append("linked report contains embedded image data")
        expected_sources = [
            str(figure.get("outputs", {}).get("png"))
            for figure in rendered
            if isinstance(figure.get("outputs"), Mapping)
        ]
        if parser.image_sources != expected_sources:
            issues.append(
                f"linked report image sources do not match figure manifest: {parser.image_sources}"
            )
        for href in parser.hrefs:
            parsed = urlparse(href)
            if parsed.scheme or href.startswith("#"):
                continue
            try:
                linked = _relative_path(validated.root, href, label="report href")
            except FigureRenderError as exc:
                issues.append(str(exc))
            else:
                if not linked.is_file():
                    issues.append(f"linked report href does not exist: {href}")
    for doi in _EXPECTED_DOIS:
        if f"https://doi.org/{doi}" not in parser.hrefs:
            issues.append(f"report {path.name!r} is missing DOI link {doi}")
    return {
        "path": path.name,
        "sha256": _sha256(path),
        "size_bytes": path.stat().st_size,
        "images": len(parser.image_sources),
        "figures": parser.figure_count,
        "captions": parser.figcaption_count,
    }


def _post_render_audit(
    validated: ValidatedRun,
) -> tuple[list[str], list[str], Mapping[str, Any]]:
    issues: list[str] = []
    warnings = list(validated.warnings)
    manifest_path = validated.root / "figures" / "figure_manifest.json"
    linked_path = validated.root / "report.html"
    standalone_path = validated.root / "report_standalone.html"
    for path in (manifest_path, linked_path, standalone_path):
        if not path.is_file():
            issues.append(f"required publication output is missing: {path.name}")
    if issues:
        return issues, warnings, {"publication_outputs": {"status": "failed"}}
    try:
        figure_manifest = _load_figure_manifest(manifest_path, validated.root)
    except FigureRenderError as exc:
        issues.append(str(exc))
        return issues, warnings, {"publication_outputs": {"status": "failed"}}

    rendered = [
        figure
        for figure_id in _FIGURE_ORDER
        for figure in figure_manifest.figures
        if figure.get("id") == figure_id and figure.get("status") == "rendered"
    ]
    artwork = _artwork_checks(validated, figure_manifest, issues=issues)
    linked = _audit_html(
        linked_path,
        validated=validated,
        rendered=rendered,
        standalone=False,
        issues=issues,
    )
    standalone = _audit_html(
        standalone_path,
        validated=validated,
        rendered=rendered,
        standalone=True,
        issues=issues,
    )
    checks: Mapping[str, Any] = {
        "artifact_contract": {
            "status": "pass",
            "roles_checked": len(validated.artifacts),
        },
        "artwork": {
            "status": "fail" if issues else "pass",
            "rendered_figures": len(rendered),
            "files": artwork,
        },
        "linked_report": linked,
        "standalone_report": standalone,
    }
    return issues, warnings, checks


def _write_report_validation(validated: ValidatedRun) -> Path:
    issues, warnings, checks = _post_render_audit(validated)
    path = validated.root / _REPORT_VALIDATION_NAME
    payload = {
        "schema_version": RUN_SCHEMA_VERSION,
        "valid": not issues,
        "issues": issues,
        "warnings": warnings,
        "checks": checks,
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    if issues:
        raise FigureRenderError(
            "post-render publication validation failed:\n"
            + "\n".join(f"- {issue}" for issue in issues)
        )
    return path


def build_publication_report(
    run_dir: str | Path,
    *,
    manifest_name: str = "00_manifest.json",
    figure_manifest_path: str | Path = "figures/figure_manifest.json",
    report_name: str = "report.html",
    standalone_name: str = "report_standalone.html",
) -> ReportBuildResult:
    """Build deterministic linked/standalone English HTML and validate the bundle."""

    validated = validate_run(run_dir, manifest_name=manifest_name)
    raw_manifest = Path(figure_manifest_path)
    manifest_path = (
        raw_manifest.expanduser().resolve()
        if raw_manifest.is_absolute()
        else _relative_path(
            validated.root, raw_manifest.as_posix(), label="figure manifest"
        )
    )
    figure_manifest = _load_figure_manifest(manifest_path, validated.root)
    report_path = _relative_path(validated.root, report_name, label="report output")
    standalone_path = _relative_path(
        validated.root, standalone_name, label="standalone report output"
    )
    report_path.write_text(
        _document(validated, figure_manifest, standalone=False), encoding="utf-8"
    )
    standalone_path.write_text(
        _document(validated, figure_manifest, standalone=True), encoding="utf-8"
    )
    validation_path = _write_report_validation(validated)
    return ReportBuildResult(
        report_html=report_path,
        report_standalone_html=standalone_path,
        figure_manifest=manifest_path,
        report_validation=validation_path,
        validation_warnings=validated.warnings,
    )


def build_publication_bundle(
    run_dir: str | Path,
    *,
    manifest_name: str = "00_manifest.json",
    rscript: str | Path = "Rscript",
    renderer: str | Path | None = None,
) -> PublicationBundle:
    """Validate source data, render figures, and build both report variants."""

    validated = validate_run(run_dir, manifest_name=manifest_name)
    figures = render_publication_figures(
        validated.root,
        manifest_name=manifest_name,
        rscript=rscript,
        renderer=renderer,
    )
    report = build_publication_report(
        validated.root,
        manifest_name=manifest_name,
        figure_manifest_path=figures.path,
    )
    return PublicationBundle(validated_run=validated, figures=figures, report=report)


def validate_production_run(
    run_dir: str | Path,
    *,
    manifest_name: str = "00_manifest.json",
) -> ValidationReport:
    """Return pre-render warnings or a strict post-render validation report."""

    root = Path(run_dir).expanduser().resolve()
    publication_paths = (
        root / "figures" / "figure_manifest.json",
        root / "report.html",
        root / "report_standalone.html",
        root / _REPORT_VALIDATION_NAME,
    )
    phase = (
        "post-render"
        if any(path.exists() for path in publication_paths)
        else "pre-render"
    )
    try:
        validated = validate_run(root, manifest_name=manifest_name)
    except RunValidationError as exc:
        return ValidationReport(
            valid=False,
            issues=exc.issues,
            warnings=(),
            run=None,
            phase=phase,  # type: ignore[arg-type]
        )
    if phase == "pre-render":
        return ValidationReport(
            valid=True,
            issues=(),
            warnings=(
                *validated.warnings,
                "publication figures and reports have not been rendered; source-artifact validation only",
            ),
            run=validated,
            phase="pre-render",
            checks={"artifact_contract": {"status": "pass"}},
        )

    issues, warnings, checks = _post_render_audit(validated)
    record_path = root / _REPORT_VALIDATION_NAME
    if not record_path.is_file():
        issues.append(
            f"required publication output is missing: {_REPORT_VALIDATION_NAME}"
        )
    else:
        try:
            record = json.loads(record_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
            issues.append(f"report_validation.json is not valid UTF-8 JSON: {exc}")
        else:
            if not isinstance(record, Mapping):
                issues.append("report_validation.json must contain an object")
            else:
                if record.get("schema_version") != RUN_SCHEMA_VERSION:
                    issues.append("report_validation.json has the wrong schema_version")
                if record.get("valid") is not True or record.get("issues") != []:
                    issues.append(
                        "report_validation.json does not record a successful render"
                    )
                if record.get("checks") != checks:
                    issues.append("report_validation.json checks are stale")
    return ValidationReport(
        valid=not issues,
        issues=tuple(issues),
        warnings=tuple(warnings),
        run=validated,
        phase="post-render",
        checks=checks,
        report_validation=record_path if record_path.is_file() else None,
    )


def render_production_report(
    run_dir: str | Path,
    renderer: str = "nature-r",
) -> PublicationBundle:
    """Render the canonical production report through the named publication backend."""

    if renderer != "nature-r":
        raise ValueError(
            f"unknown production report renderer {renderer!r}; use 'nature-r'"
        )
    return build_publication_bundle(run_dir)
