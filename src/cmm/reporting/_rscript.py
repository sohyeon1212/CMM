"""Invoking a checked-in R renderer, shared by the workflows that have one.

Extracted from :mod:`cmm.reporting.publication` when a second workflow needed the same
invocation. Only the mechanics live here — locating ``Rscript``, pinning the environment so a
render is reproducible, running the script, and refusing to accept a render that emitted a
warning. Which figures exist, what they are drawn from, and what counts as required are the
calling workflow's business and stay in its own module.

The leading underscore is the whole promise this module carries: it is package-internal, and a
workflow outside CMM must not import it.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
from typing import Any, Mapping

from cmm.reporting.schema import RUN_SCHEMA_VERSION


def renderer_environment(
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


def decode_renderer_stream(value: bytes | str | None) -> str:
    """Decode external renderer output explicitly and loss-tolerantly as UTF-8."""

    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


def invoke_renderer(
    *,
    script: Path,
    script_sha256: str,
    run_dir: Path,
    manifest_path: Path,
    figures_dir: Path,
    figure_manifest_path: Path,
    rscript: str | Path,
    error_type: type[Exception],
    label: str = "R renderer",
    extra_environment: dict[str, str] | None = None,
) -> None:
    """Run ``script`` over one run directory, raising ``error_type`` on any failure.

    A warning is treated as a failure. R reports a dropped row, an invalid scale and a font
    substitution as warnings rather than errors, and each of those changes what the figure
    shows while still producing a file that opens.
    """

    if not script.is_file():
        raise error_type(f"{label} is missing: {script}")
    executable = shutil.which(os.fspath(rscript))
    if executable is None:
        raise error_type(
            f"Rscript executable {os.fspath(rscript)!r} was not found; "
            "publication figures require R"
        )
    environment = renderer_environment(script_sha256)
    if extra_environment:
        environment.update(extra_environment)
    completed = subprocess.run(
        [
            executable,
            "--vanilla",
            os.fspath(script),
            os.fspath(run_dir),
            os.fspath(manifest_path),
            os.fspath(figures_dir),
            os.fspath(figure_manifest_path),
        ],
        cwd=run_dir,
        env=environment,
        check=False,
        capture_output=True,
    )
    renderer_output = "\n".join(
        output
        for value in (completed.stdout, completed.stderr)
        if (output := decode_renderer_stream(value).strip())
    )
    if completed.returncode != 0:
        raise error_type(
            f"{label} failed ({completed.returncode}): "
            f"{renderer_output or 'no renderer output'}"
        )
    if re.search(r"(^|\n)Warning(?: message)?", renderer_output, flags=re.IGNORECASE):
        raise error_type(f"{label} emitted a warning: {renderer_output}")
    if not figure_manifest_path.is_file():
        raise error_type(
            f"{label} completed without writing "
            f"{figure_manifest_path.parent.name}/{figure_manifest_path.name}"
        )


@dataclass(frozen=True)
class FigureManifest:
    path: Path
    renderer: Mapping[str, Any]
    figures: tuple[Mapping[str, Any], ...]


def relative_path(
    root: Path, value: object, *, label: str, error_type: type[Exception]
) -> Path:
    if not isinstance(value, str) or not value:
        raise error_type(f"{label} must be a non-empty relative path")
    pure = PurePosixPath(value)
    if pure.is_absolute() or ".." in pure.parts:
        raise error_type(f"{label} escapes the run directory: {value!r}")
    path = (root / Path(*pure.parts)).resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise error_type(f"{label} escapes the run directory: {value!r}") from exc
    return path


def load_figure_manifest(
    path: Path,
    root: Path,
    *,
    required: frozenset[str],
    order: tuple[str, ...],
    error_type: type[Exception],
) -> FigureManifest:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise error_type(f"figure manifest is not valid JSON: {path}: {exc}") from exc
    if not isinstance(data, Mapping):
        raise error_type("figure manifest must contain a JSON object")
    schema_version = data.get("schema_version")
    if type(schema_version) is not int or schema_version != RUN_SCHEMA_VERSION:
        raise error_type(
            f"figure manifest schema_version must be the integer {RUN_SCHEMA_VERSION}"
        )
    renderer = data.get("renderer", {})
    figures_value = data.get("figures")
    if not isinstance(renderer, Mapping) or not isinstance(figures_value, list):
        raise error_type("figure manifest needs renderer object and figures list")
    script_sha256 = renderer.get("script_sha256")
    if (
        not isinstance(script_sha256, str)
        or re.fullmatch(r"[0-9a-f]{64}", script_sha256) is None
    ):
        raise error_type(
            "figure manifest renderer.script_sha256 must be a lowercase SHA-256 digest"
        )

    figures: list[Mapping[str, Any]] = []
    seen: set[str] = set()
    rendered_required: set[str] = set()
    for index, raw in enumerate(figures_value):
        if not isinstance(raw, Mapping):
            raise error_type(f"figure manifest entry {index} is not an object")
        figure_id = raw.get("id")
        status = raw.get("status")
        if not isinstance(figure_id, str) or status not in {
            "rendered",
            "skipped",
            "failed",
        }:
            raise error_type(f"figure manifest entry {index} has invalid id/status")
        if figure_id in seen:
            raise error_type(f"figure manifest repeats id {figure_id!r}")
        seen.add(figure_id)
        if status == "rendered":
            outputs = raw.get("outputs")
            sources = raw.get("sources")
            if not isinstance(outputs, Mapping) or not isinstance(sources, list):
                raise error_type(
                    f"rendered figure {figure_id!r} needs outputs and sources"
                )
            for suffix in ("png", "pdf", "svg"):
                output = relative_path(
                    root,
                    outputs.get(suffix),
                    label=f"{figure_id}.{suffix}",
                    error_type=error_type,
                )
                if not output.is_file() or output.stat().st_size == 0:
                    raise error_type(
                        f"rendered figure {figure_id!r} is missing non-empty {suffix}: {output}"
                    )
            for source in sources:
                source_path = relative_path(
                    root, source, label=f"{figure_id} source", error_type=error_type
                )
                if not source_path.is_file():
                    raise error_type(
                        f"figure {figure_id!r} cites a missing source artifact: {source_path}"
                    )
            if figure_id in required:
                rendered_required.add(figure_id)
        else:
            reason = raw.get("reason")
            if not isinstance(reason, str) or not reason.strip():
                raise error_type(
                    f"unavailable figure {figure_id!r} must state a reason"
                )
            if figure_id in required:
                raise error_type(
                    f"required figure {figure_id!r} was {status}: {reason}"
                )
        figures.append(raw)

    missing = sorted(required - rendered_required)
    if missing:
        raise error_type(f"figure manifest omitted required figure(s): {missing}")
    if seen - set(order):
        raise error_type(
            f"figure manifest contains unknown figure id(s): {sorted(seen - set(order))}"
        )
    return FigureManifest(path=path, renderer=renderer, figures=tuple(figures))
