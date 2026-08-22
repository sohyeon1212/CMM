"""Schema-v2 validation for reproducible CMM publication runs.

``00_manifest.json`` is the only path-discovery surface. This module validates its stable
role mapping, every declared digest and byte count, CSV metadata sidecars, and a compact set
of scientific invariants before any report text or artwork is generated. Presentation code
therefore never guesses filenames or silently drops a requested method.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
from typing import Any, Literal, Mapping, Sequence

RUN_SCHEMA_VERSION = 2
_SHA256 = re.compile(r"[0-9a-f]{64}")
_FLOAT_TOLERANCE = 1e-9

ArtifactStatus = Literal["complete", "partial", "skipped", "failed"]


@dataclass(frozen=True)
class ArtifactContract:
    kind: Literal["csv", "json", "file"]
    required: bool
    columns: tuple[str, ...] = ()


# Stable roles shared by the workflow, report builder and R renderer. Required method roles
# must be present and complete for a publication-ready SC-01 run. A deliberately disabled
# stage remains a useful partial run, but it is not silently promoted to publication-ready.
ARTIFACT_CONTRACTS: Mapping[str, ArtifactContract] = {
    "provenance": ArtifactContract("json", True),
    "summary": ArtifactContract("json", True),
    "model": ArtifactContract("file", True),
    "wild_type_fluxes": ArtifactContract("csv", True, ("reaction_id", "flux")),
    "theoretical_yield": ArtifactContract(
        "csv",
        True,
        ("status", "molar_yield", "product_flux", "substrate_uptake"),
    ),
    "production_envelope": ArtifactContract(
        "csv", True, ("product_flux", "growth_min", "growth_max")
    ),
    "reference_fluxes": ArtifactContract("csv", True, ("reaction_id", "flux")),
    "single_knockout_moma": ArtifactContract(
        "csv", True, ("target_id", "kind", "status", "objective", "product_flux")
    ),
    "single_knockout_room": ArtifactContract(
        "csv", True, ("target_id", "kind", "status", "objective", "product_flux")
    ),
    "optknock": ArtifactContract(
        "csv",
        True,
        ("knockouts", "growth", "max_product", "guaranteed_product", "growth_coupled"),
    ),
    "robustknock": ArtifactContract(
        "csv",
        True,
        ("knockouts", "growth", "max_product", "guaranteed_product", "growth_coupled"),
    ),
    "fseof_tidy": ArtifactContract(
        "csv", True, ("target", "enforced_product_flux", "reaction_flux")
    ),
    "fvseof_tidy": ArtifactContract(
        "csv",
        True,
        ("target", "enforced_product_flux", "mean_flux", "forced_min_flux"),
    ),
    "amplification_target_ranking": ArtifactContract(
        "csv", False, ("reaction_id", "slope", "classification", "actionable")
    ),
    "variability_supported_amplification_targets": ArtifactContract(
        "csv",
        False,
        ("reaction_id", "classification", "robust", "slope", "actionable"),
    ),
    # Supporting roles are declared when attempted, with a reason for a skipped/failed state.
    # Some, including the required preflight table, are promoted by semantic role from the
    # supplementary inventory rather than duplicated in the primary mapping.
    "preflight_checks": ArtifactContract(
        "csv", True, ("check", "status", "value", "message")
    ),
    "workflow_configuration": ArtifactContract("json", False),
    "single_knockout_consensus": ArtifactContract(
        "csv", False, ("target_id", "recommended")
    ),
    "gene_knockout_mapping": ArtifactContract(
        "csv",
        False,
        ("gene_id", "gene_name", "inert", "blocked_reaction", "reaction_name", "gpr"),
    ),
    "flux_response_tidy": ArtifactContract(
        "csv",
        False,
        (
            "target",
            "target_flux",
            "response_flux",
            "biomass_flux",
            "status",
            "scan_reaction",
            "response_reaction",
            "background",
        ),
    ),
    "sampling_tidy": ArtifactContract(
        "csv", False, ("target", "condition", "reaction_id", "flux")
    ),
    # Candidate-level execution indexes make failed analyses visible even though they have
    # no tidy trajectory/distribution rows.  They are optional additions to schema v2 so
    # preserved runs written before exhaustive validation remain readable and valid.
    "flux_response_validation_index": ArtifactContract(
        "csv",
        False,
        (
            "target",
            "scan_reaction",
            "response_reaction",
            "background",
            "status",
            "error",
            "data_file",
            "phases_file",
            "metadata_file",
        ),
    ),
    "single_knockout_sampling_validation_index": ArtifactContract(
        "csv",
        False,
        (
            "target_id",
            "status",
            "error",
            "samples_file",
            "statistics_file",
            "comparison_file",
            "metadata_file",
        ),
    ),
    "amplification_loop_diagnostic": ArtifactContract(
        "csv",
        False,
        (
            "rank",
            "target",
            "source_methods",
            "standard_minimum",
            "standard_maximum",
            "standard_capacity",
            "loopless_minimum",
            "loopless_maximum",
            "loopless_capacity",
            "loopless_to_standard_capacity_ratio",
            "capacity_ratio_threshold",
            "loop_artifact_flag",
            "diagnostic_status",
            "reason",
            "enforced_product_floor",
            "biomass_floor",
        ),
    ),
    "recommendations": ArtifactContract(
        "csv",
        False,
        (
            "target",
            "type",
            "evidence",
            "verdict",
            "proposal_methods",
            "validation_methods",
            "growth_retained",
            "product_effect",
            "artifact_flag",
            "reason",
        ),
    ),
    "reproduction_config": ArtifactContract("json", False),
    "reproduce_script": ArtifactContract("file", False),
    "render_script": ArtifactContract("file", False),
    "validate_script": ArtifactContract("file", False),
}


@dataclass(frozen=True)
class RunArtifact:
    role: str
    relative_path: str | None
    path: Path | None
    kind: str
    required: bool
    status: ArtifactStatus
    reason: str | None = None
    sha256: str | None = None
    size_bytes: int | None = None
    metadata_relative_path: str | None = None
    metadata_path: Path | None = None
    source: Literal["primary", "supplementary"] = "primary"

    @property
    def available(self) -> bool:
        return self.status in {"complete", "partial"} and self.path is not None


@dataclass(frozen=True)
class ValidatedRun:
    root: Path
    manifest_path: Path
    manifest: Mapping[str, Any]
    artifacts: Mapping[str, RunArtifact]
    warnings: tuple[str, ...] = ()

    def artifact(self, role: str, *, required: bool = True) -> Path | None:
        """Return an available artifact path, optionally accepting an absent optional role."""

        artifact = self.artifacts.get(role)
        if artifact is None or not artifact.available:
            if required:
                raise KeyError(f"artifact role {role!r} is not available")
            return None
        return artifact.path

    @property
    def report(self) -> Mapping[str, Any]:
        value = self.manifest.get("report", {})
        return value if isinstance(value, Mapping) else {}


class RunValidationError(ValueError):
    """Raised with every schema or scientific-contract violation found in one pass."""

    def __init__(self, issues: list[str] | tuple[str, ...]):
        self.issues = tuple(issues)
        detail = "\n".join(f"- {issue}" for issue in self.issues)
        super().__init__(f"CMM run schema v2 validation failed:\n{detail}")


@dataclass(frozen=True)
class ValidationReport:
    """Non-raising validation surface for CLIs and agent workflows."""

    valid: bool
    issues: tuple[str, ...]
    warnings: tuple[str, ...]
    run: ValidatedRun | None = None
    phase: Literal["pre-render", "post-render"] = "pre-render"
    checks: Mapping[str, Any] = field(default_factory=dict)
    report_validation: Path | None = None

    def raise_for_errors(self) -> ValidatedRun:
        if not self.valid or self.run is None:
            raise RunValidationError(
                self.issues or ("validation did not return a run",)
            )
        return self.run


def _read_json(path: Path, *, label: str, issues: list[str]) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        issues.append(f"{label} is not valid UTF-8 JSON: {exc}")
        return None


def _resolve_artifact_path(
    root: Path, value: object, *, role: str, issues: list[str]
) -> Path | None:
    if not isinstance(value, str) or not value.strip():
        issues.append(f"artifact {role!r} must declare a non-empty relative path")
        return None
    pure = PurePosixPath(value)
    if pure.is_absolute() or ".." in pure.parts:
        issues.append(
            f"artifact {role!r} path must stay inside the run directory: {value!r}"
        )
        return None
    candidate = (root / Path(*pure.parts)).resolve()
    try:
        candidate.relative_to(root)
    except ValueError:
        issues.append(
            f"artifact {role!r} resolves outside the run directory: {value!r}"
        )
        return None
    return candidate


def _normalise_entry(
    role: str,
    raw: object,
    contract: ArtifactContract,
    *,
    root: Path,
    issues: list[str],
    source: Literal["primary", "supplementary"] = "primary",
) -> RunArtifact:
    if isinstance(raw, str):
        entry: Mapping[str, object] = {"path": raw, "status": "complete"}
    elif isinstance(raw, Mapping):
        entry = raw
    else:
        issues.append(f"artifact {role!r} must be a path string or object")
        entry = {}

    status_value = entry.get("status", "complete")
    if status_value == "disabled" and source == "supplementary":
        status_value = "skipped"
    allowed = {"complete", "partial", "skipped", "failed"}
    if not isinstance(status_value, str) or status_value not in allowed:
        issues.append(
            f"artifact {role!r} has invalid status {status_value!r}; expected one of {sorted(allowed)}"
        )
        status: ArtifactStatus = "failed"
    else:
        status = status_value  # type: ignore[assignment]

    reason_value = entry.get("reason")
    reason = str(reason_value).strip() if reason_value is not None else None
    if source == "supplementary" and status == "skipped" and not reason:
        reason = "analysis was disabled"
    path_value = entry.get("path")
    relative_path = str(path_value) if isinstance(path_value, str) else None
    path = None
    if status in {"complete", "partial"} or path_value is not None:
        path = _resolve_artifact_path(root, path_value, role=role, issues=issues)
    if status in {"skipped", "failed"} and not reason:
        issues.append(f"artifact {role!r} with status {status!r} must state a reason")
    if contract.required and status != "complete":
        issues.append(
            f"required artifact {role!r} must have status 'complete', not {status!r}"
        )

    metadata_value = entry.get("metadata_path")
    metadata_relative = str(metadata_value) if isinstance(metadata_value, str) else None
    metadata_path = None
    if metadata_value is not None:
        metadata_path = _resolve_artifact_path(
            root, metadata_value, role=f"{role} metadata", issues=issues
        )
    sha_value = entry.get("sha256")
    size_value = entry.get("size_bytes")
    return RunArtifact(
        role=role,
        relative_path=relative_path,
        path=path,
        kind=contract.kind,
        required=contract.required,
        status=status,
        reason=reason,
        sha256=sha_value if isinstance(sha_value, str) else None,
        size_bytes=(
            size_value
            if type(size_value) is int
            and isinstance(size_value, int)
            and size_value >= 0
            else None
        ),
        metadata_relative_path=metadata_relative,
        metadata_path=metadata_path,
        source=source,
    )


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validate_integrity_fields(
    *,
    path: Path,
    role: str,
    raw: Mapping[str, object],
    issues: list[str],
    allow_missing: bool = False,
) -> None:
    expected_hash = raw.get("sha256")
    expected_size = raw.get("size_bytes")
    if not isinstance(expected_hash, str) or _SHA256.fullmatch(expected_hash) is None:
        if not allow_missing:
            issues.append(f"artifact {role!r} must declare a lowercase SHA-256 digest")
    elif path.is_file() and _file_sha256(path) != expected_hash:
        issues.append(f"artifact {role!r} SHA-256 does not match its manifest entry")
    if (
        type(expected_size) is not int
        or not isinstance(expected_size, int)
        or expected_size < 0
    ):
        if not allow_missing:
            issues.append(
                f"artifact {role!r} must declare a non-negative integer size_bytes"
            )
    elif path.is_file() and path.stat().st_size != expected_size:
        issues.append(
            f"artifact {role!r} byte size does not match its manifest entry: "
            f"expected {expected_size}, found {path.stat().st_size}"
        )


def _read_csv_rows(
    path: Path, *, role: str, issues: list[str]
) -> tuple[list[str], list[dict[str, str]]] | None:
    try:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            header = list(reader.fieldnames or ())
            rows = [dict(row) for row in reader]
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        issues.append(f"artifact {role!r} is not readable CSV: {exc}")
        return None
    if not header:
        issues.append(f"artifact {role!r} has no CSV header")
        return None
    return header, rows


def _validate_file(
    artifact: RunArtifact,
    contract: ArtifactContract,
    raw_entry: object,
    *,
    supplementary_by_path: Mapping[str, Mapping[str, object]],
    issues: list[str],
    warnings: list[str],
) -> tuple[list[str], list[dict[str, str]]] | Mapping[str, Any] | None:
    if not artifact.available or artifact.path is None:
        if not artifact.required and artifact.status in {"skipped", "failed"}:
            warnings.append(
                f"optional artifact {artifact.role!r} is {artifact.status}: {artifact.reason}"
            )
        return None
    path = artifact.path
    if not path.is_file():
        issues.append(f"artifact {artifact.role!r} does not exist as a file: {path}")
        return None
    if path.stat().st_size == 0:
        issues.append(f"artifact {artifact.role!r} is empty: {path}")
        return None

    if not isinstance(raw_entry, Mapping):
        issues.append(
            f"artifact {artifact.role!r} must use an object entry with sha256 and size_bytes"
        )
        raw_mapping: Mapping[str, object] = {}
    else:
        raw_mapping = raw_entry
    _validate_integrity_fields(
        path=path, role=artifact.role, raw=raw_mapping, issues=issues
    )

    if contract.kind == "json":
        parsed = _read_json(path, label=f"artifact {artifact.role!r}", issues=issues)
        if parsed is not None and not isinstance(parsed, Mapping):
            issues.append(f"artifact {artifact.role!r} JSON must contain an object")
            return None
        return parsed if isinstance(parsed, Mapping) else None
    if contract.kind != "csv":
        return None

    parsed_csv = _read_csv_rows(path, role=artifact.role, issues=issues)
    if parsed_csv is None:
        return None
    header, rows = parsed_csv
    expected = list(contract.columns)
    declared = raw_mapping.get("columns", ())
    if isinstance(declared, list) and all(isinstance(item, str) for item in declared):
        expected.extend(str(item) for item in declared)
    elif declared not in (None, ()):
        issues.append(f"artifact {artifact.role!r} columns must be a list of strings")
    missing = sorted(set(expected) - set(header))
    if missing:
        issues.append(
            f"artifact {artifact.role!r} is missing CSV column(s) {missing}; found {header}"
        )
    if not rows:
        warnings.append(f"artifact {artifact.role!r} has a header but no data rows")

    if artifact.status == "complete":
        if artifact.metadata_path is None or artifact.metadata_relative_path is None:
            issues.append(
                f"complete CSV artifact {artifact.role!r} must declare metadata_path"
            )
        elif not artifact.metadata_path.is_file():
            issues.append(
                f"metadata sidecar for artifact {artifact.role!r} does not exist: "
                f"{artifact.metadata_path}"
            )
        else:
            metadata = _read_json(
                artifact.metadata_path,
                label=f"metadata sidecar for artifact {artifact.role!r}",
                issues=issues,
            )
            if metadata is not None and not isinstance(metadata, Mapping):
                issues.append(
                    f"metadata sidecar for artifact {artifact.role!r} must contain an object"
                )
            sidecar_record = supplementary_by_path.get(artifact.metadata_relative_path)
            if sidecar_record is None:
                issues.append(
                    f"metadata sidecar for artifact {artifact.role!r} is not listed in "
                    "supplementary_artifacts"
                )
            else:
                _validate_integrity_fields(
                    path=artifact.metadata_path,
                    role=f"{artifact.role} metadata",
                    raw=sidecar_record,
                    issues=issues,
                )
    return parsed_csv


def _finite_number(
    value: object,
    *,
    role: str,
    column: str,
    row_number: int,
    issues: list[str],
    optional: bool = False,
) -> float | None:
    if value is None or str(value).strip() in {"", "NA", "NaN", "nan"}:
        if not optional:
            issues.append(
                f"artifact {role!r} row {row_number} has no finite {column!r} value"
            )
        return None
    try:
        number = float(str(value))
    except ValueError:
        number = math.nan
    if not math.isfinite(number):
        if not optional:
            issues.append(
                f"artifact {role!r} row {row_number} has non-finite {column!r}: {value!r}"
            )
        return None
    return number


def _flag(value: object) -> bool | None:
    normalised = str(value).strip().lower()
    if normalised in {"true", "1", "yes", "y"}:
        return True
    if normalised in {"false", "0", "no", "n"}:
        return False
    return None


def _validate_numeric_rows(
    role: str,
    rows: Sequence[Mapping[str, str]],
    *,
    issues: list[str],
    warnings: list[str],
) -> None:
    if role in {"wild_type_fluxes", "reference_fluxes"}:
        identifiers: set[str] = set()
        for index, row in enumerate(rows, start=2):
            reaction = row.get("reaction_id", "").strip()
            if not reaction:
                issues.append(f"artifact {role!r} row {index} has an empty reaction_id")
            elif reaction in identifiers:
                issues.append(f"artifact {role!r} repeats reaction_id {reaction!r}")
            identifiers.add(reaction)
            _finite_number(
                row.get("flux"),
                role=role,
                column="flux",
                row_number=index,
                issues=issues,
            )
        return

    if role == "theoretical_yield":
        if len(rows) != 1:
            issues.append(
                "artifact 'theoretical_yield' must contain exactly one data row"
            )
        for index, row in enumerate(rows, start=2):
            if row.get("status") != "optimal":
                issues.append("theoretical_yield status must be 'optimal'")
            for column in ("molar_yield", "product_flux", "substrate_uptake"):
                number = _finite_number(
                    row.get(column),
                    role=role,
                    column=column,
                    row_number=index,
                    issues=issues,
                )
                if number is not None and number <= 0:
                    issues.append(f"theoretical_yield {column} must be positive")
        return

    if role == "production_envelope":
        previous_product: float | None = None
        for index, row in enumerate(rows, start=2):
            product = _finite_number(
                row.get("product_flux"),
                role=role,
                column="product_flux",
                row_number=index,
                issues=issues,
            )
            lower = _finite_number(
                row.get("growth_min"),
                role=role,
                column="growth_min",
                row_number=index,
                issues=issues,
            )
            upper = _finite_number(
                row.get("growth_max"),
                role=role,
                column="growth_max",
                row_number=index,
                issues=issues,
            )
            if (
                lower is not None
                and upper is not None
                and lower > upper + _FLOAT_TOLERANCE
            ):
                issues.append(
                    f"production_envelope row {index} has growth_min greater than growth_max"
                )
            if lower is not None and lower < -_FLOAT_TOLERANCE:
                issues.append(
                    f"production_envelope row {index} has negative growth_min"
                )
            if upper is not None and upper < -_FLOAT_TOLERANCE:
                issues.append(
                    f"production_envelope row {index} has negative growth_max"
                )
            if (
                product is not None
                and previous_product is not None
                and product < previous_product
            ):
                issues.append("production_envelope product_flux must be non-decreasing")
            if product is not None:
                previous_product = product
        if len(rows) < 2:
            issues.append("production_envelope needs at least two scan points")
        return

    if role in {"single_knockout_moma", "single_knockout_room"}:
        for index, row in enumerate(rows, start=2):
            if not row.get("target_id", "").strip():
                issues.append(f"artifact {role!r} row {index} has an empty target_id")
            if row.get("status") == "optimal":
                growth = _finite_number(
                    row.get("objective"),
                    role=role,
                    column="objective",
                    row_number=index,
                    issues=issues,
                )
                _finite_number(
                    row.get("product_flux"),
                    role=role,
                    column="product_flux",
                    row_number=index,
                    issues=issues,
                )
                if growth is not None and growth < -_FLOAT_TOLERANCE:
                    issues.append(f"artifact {role!r} row {index} has negative growth")
        return

    if role in {"optknock", "robustknock"}:
        previous_guaranteed: float | None = None
        for index, row in enumerate(rows, start=2):
            growth = _finite_number(
                row.get("growth"),
                role=role,
                column="growth",
                row_number=index,
                issues=issues,
            )
            maximum = _finite_number(
                row.get("max_product"),
                role=role,
                column="max_product",
                row_number=index,
                issues=issues,
            )
            guaranteed = _finite_number(
                row.get("guaranteed_product"),
                role=role,
                column="guaranteed_product",
                row_number=index,
                issues=issues,
            )
            if growth is not None and growth < -_FLOAT_TOLERANCE:
                issues.append(f"artifact {role!r} row {index} has negative growth")
            if maximum is not None and maximum < -_FLOAT_TOLERANCE:
                issues.append(f"artifact {role!r} row {index} has negative max_product")
            if guaranteed is not None and guaranteed < -_FLOAT_TOLERANCE:
                issues.append(
                    f"artifact {role!r} row {index} has negative guaranteed_product"
                )
            if (
                maximum is not None
                and guaranteed is not None
                and guaranteed > maximum + _FLOAT_TOLERANCE
            ):
                issues.append(
                    f"artifact {role!r} row {index} has guaranteed_product above max_product"
                )
            coupled = _flag(row.get("growth_coupled"))
            if coupled is None:
                issues.append(
                    f"artifact {role!r} row {index} has invalid growth_coupled flag"
                )
            elif guaranteed is not None and coupled != (guaranteed > 1e-6):
                issues.append(
                    f"artifact {role!r} row {index} has a growth_coupled verdict inconsistent "
                    "with guaranteed_product"
                )
            if (
                guaranteed is not None
                and previous_guaranteed is not None
                and guaranteed > previous_guaranteed + _FLOAT_TOLERANCE
            ):
                issues.append(f"artifact {role!r} is not ranked by guaranteed_product")
            if guaranteed is not None:
                previous_guaranteed = guaranteed
        return

    trajectory_columns = {
        "fseof_tidy": ("enforced_product_flux", "reaction_flux"),
        "fvseof_tidy": ("enforced_product_flux", "mean_flux", "forced_min_flux"),
    }
    if role in trajectory_columns:
        points_by_target: dict[str, int] = {}
        for index, row in enumerate(rows, start=2):
            target = row.get("target", "").strip()
            if not target:
                issues.append(f"artifact {role!r} row {index} has an empty target")
            else:
                points_by_target[target] = points_by_target.get(target, 0) + 1
            for column in trajectory_columns[role]:
                _finite_number(
                    row.get(column),
                    role=role,
                    column=column,
                    row_number=index,
                    issues=issues,
                )
        undersampled = sorted(
            target for target, count in points_by_target.items() if count < 2
        )
        if undersampled:
            issues.append(
                f"artifact {role!r} has targets with fewer than two scan points: {undersampled}"
            )
        return

    if role == "flux_response_tidy":
        allowed_backgrounds = {"wild_type", "gene_knockout"}
        allowed_scopes = {
            "all_report_selected_candidates",
            "all_display_ranked_candidates",
        }
        for index, row in enumerate(rows, start=2):
            for column in ("target", "scan_reaction", "response_reaction"):
                if not row.get(column, "").strip():
                    issues.append(
                        f"artifact {role!r} row {index} has an empty {column}"
                    )
            if row.get("background") not in allowed_backgrounds:
                issues.append(
                    f"artifact {role!r} row {index} has invalid background {row.get('background')!r}"
                )
            scope = row.get("candidate_scope", "").strip()
            if scope and scope not in allowed_scopes:
                issues.append(
                    f"artifact {role!r} row {index} has invalid candidate_scope {scope!r}"
                )
            _finite_number(
                row.get("target_flux"),
                role=role,
                column="target_flux",
                row_number=index,
                issues=issues,
            )
            if row.get("status") == "optimal":
                _finite_number(
                    row.get("response_flux"),
                    role=role,
                    column="response_flux",
                    row_number=index,
                    issues=issues,
                )
                _finite_number(
                    row.get("biomass_flux"),
                    role=role,
                    column="biomass_flux",
                    row_number=index,
                    issues=issues,
                    optional=True,
                )
        return

    if role == "sampling_tidy":
        conditions_by_target: dict[str, set[str]] = {}
        for index, row in enumerate(rows, start=2):
            target = row.get("target", "").strip()
            condition = row.get("condition", "").strip()
            if condition not in {"wild_type", "knockout"}:
                issues.append(
                    f"artifact {role!r} row {index} has invalid condition {condition!r}"
                )
            if target:
                conditions_by_target.setdefault(target, set()).add(condition)
            _finite_number(
                row.get("flux"),
                role=role,
                column="flux",
                row_number=index,
                issues=issues,
            )
        unmatched = sorted(
            target
            for target, conditions in conditions_by_target.items()
            if conditions != {"wild_type", "knockout"}
        )
        if unmatched:
            issues.append(
                "sampling_tidy must contain matched wild_type and knockout ensembles for "
                f"each target; unmatched: {unmatched}"
            )
        return

    if role == "flux_response_validation_index":
        seen: set[tuple[str, str, str]] = set()
        allowed_scopes = {
            "all_report_selected_candidates",
            "all_display_ranked_candidates",
        }
        for index, row in enumerate(rows, start=2):
            target = row.get("target", "").strip()
            background = row.get("background", "").strip()
            scope = row.get("candidate_scope", "").strip()
            status = row.get("status", "").strip().lower()
            key = (target, background, scope)
            if not target:
                issues.append(f"artifact {role!r} row {index} has an empty target")
            if background not in {"wild_type", "gene_knockout"}:
                issues.append(
                    f"artifact {role!r} row {index} has invalid background {background!r}"
                )
            if scope and scope not in allowed_scopes:
                issues.append(
                    f"artifact {role!r} row {index} has invalid candidate_scope {scope!r}"
                )
            if status not in {"complete", "failed", "skipped"}:
                issues.append(
                    f"artifact {role!r} row {index} has invalid status {status!r}"
                )
            if target and background:
                if key in seen:
                    issues.append(
                        f"artifact {role!r} repeats target/background/candidate_scope {key!r}"
                    )
                seen.add(key)
            if status == "complete":
                for column in ("data_file", "phases_file", "metadata_file"):
                    if not row.get(column, "").strip():
                        issues.append(
                            f"artifact {role!r} row {index} is complete but has no {column}"
                        )
            elif status == "failed" and not row.get("error", "").strip():
                issues.append(
                    f"artifact {role!r} row {index} failed without an error reason"
                )
            elif status == "skipped" and not row.get("reason", "").strip():
                issues.append(
                    f"artifact {role!r} row {index} was skipped without a reason"
                )
        return

    if role == "single_knockout_sampling_validation_index":
        seen_targets: set[str] = set()
        for index, row in enumerate(rows, start=2):
            target = row.get("target_id", "").strip()
            status = row.get("status", "").strip().lower()
            if not target:
                issues.append(f"artifact {role!r} row {index} has an empty target_id")
            elif target in seen_targets:
                issues.append(f"artifact {role!r} repeats target_id {target!r}")
            else:
                seen_targets.add(target)
            if status not in {"complete", "failed", "skipped"}:
                issues.append(
                    f"artifact {role!r} row {index} has invalid status {status!r}"
                )
            if status == "complete":
                for column in ("statistics_file", "metadata_file"):
                    if not row.get(column, "").strip():
                        issues.append(
                            f"artifact {role!r} row {index} is complete but has no {column}"
                        )
                if target != "wild_type" and not row.get("comparison_file", "").strip():
                    issues.append(
                        f"artifact {role!r} row {index} is a complete knockout but has no comparison_file"
                    )
            elif status == "failed" and not row.get("error", "").strip():
                issues.append(
                    f"artifact {role!r} row {index} failed without an error reason"
                )
            elif status == "skipped" and not row.get("reason", "").strip():
                issues.append(
                    f"artifact {role!r} row {index} was skipped without a reason"
                )
        return

    if role == "amplification_loop_diagnostic":
        for index, row in enumerate(rows, start=2):
            if not row.get("target", "").strip():
                issues.append(f"artifact {role!r} row {index} has an empty target")
            for column in (
                "standard_minimum",
                "standard_maximum",
                "standard_capacity",
                "loopless_minimum",
                "loopless_maximum",
                "loopless_capacity",
                "loopless_to_standard_capacity_ratio",
                "capacity_ratio_threshold",
                "enforced_product_floor",
                "biomass_floor",
            ):
                _finite_number(
                    row.get(column),
                    role=role,
                    column=column,
                    row_number=index,
                    issues=issues,
                    optional=row.get("diagnostic_status") != "complete",
                )
            if (
                row.get("diagnostic_status") == "complete"
                and _flag(row.get("loop_artifact_flag")) is None
            ):
                issues.append(
                    f"artifact {role!r} row {index} has invalid loop_artifact_flag"
                )


def _validate_cross_artifact_invariants(
    csv_rows: Mapping[str, Sequence[Mapping[str, str]]],
    json_objects: Mapping[str, Mapping[str, Any]],
    *,
    issues: list[str],
    warnings: list[str],
) -> None:
    moma = {
        row.get("target_id", "") for row in csv_rows.get("single_knockout_moma", ())
    }
    room = {
        row.get("target_id", "") for row in csv_rows.get("single_knockout_room", ())
    }
    if moma != room:
        issues.append(
            "MOMA and ROOM must cover the same knockout universe; "
            f"MOMA-only={sorted(moma - room)}, ROOM-only={sorted(room - moma)}"
        )

    provenance = json_objects.get("provenance", {})
    summary = json_objects.get("summary", {})
    provenance_hash = provenance.get("model_sha256")
    summary_hash = summary.get("model_sha256")
    if provenance_hash and summary_hash and provenance_hash != summary_hash:
        issues.append("summary model_sha256 does not match provenance model_sha256")
    if summary.get("status") not in {None, "complete"}:
        issues.append(
            f"summary status must be 'complete', got {summary.get('status')!r}"
        )

    recommendations = csv_rows.get("recommendations", ())
    responses = {
        row.get("target", "") for row in csv_rows.get("flux_response_tidy", ())
    }
    sampled = {row.get("target", "") for row in csv_rows.get("sampling_tidy", ())}
    response_index = csv_rows.get("flux_response_validation_index", ())
    if response_index:
        completed_responses = {
            row.get("target", "")
            for row in response_index
            if row.get("status", "").strip().lower() == "complete"
        }
        if completed_responses != responses:
            issues.append(
                "complete flux-response index targets must match flux_response_tidy; "
                f"index-only={sorted(completed_responses - responses)}, "
                f"tidy-only={sorted(responses - completed_responses)}"
            )
        indexed_scopes = {
            (row.get("target", ""), row.get("candidate_scope", ""))
            for row in response_index
            if row.get("status", "").strip().lower() == "complete"
            and row.get("candidate_scope", "").strip()
        }
        tidy_scopes = {
            (row.get("target", ""), row.get("candidate_scope", ""))
            for row in csv_rows.get("flux_response_tidy", ())
            if row.get("candidate_scope", "").strip()
        }
        if tidy_scopes and tidy_scopes != indexed_scopes:
            issues.append(
                "flux-response tidy/index candidate scopes must match; "
                f"index-only={sorted(indexed_scopes - tidy_scopes)}, "
                f"tidy-only={sorted(tidy_scopes - indexed_scopes)}"
            )
    sampling_index = csv_rows.get("single_knockout_sampling_validation_index", ())
    if sampling_index:
        completed_sampling = {
            row.get("target_id", "")
            for row in sampling_index
            if row.get("target_id", "") != "wild_type"
            and row.get("status", "").strip().lower() == "complete"
        }
        if completed_sampling != sampled:
            issues.append(
                "complete knockout-sampling index targets must match sampling_tidy; "
                f"index-only={sorted(completed_sampling - sampled)}, "
                f"tidy-only={sorted(sampled - completed_sampling)}"
            )

    def normalized_signature(value: str) -> str:
        return ";".join(
            sorted(part.strip() for part in value.split(";") if part.strip())
        )

    def display_ranked_signature_aliases() -> dict[str, set[str]]:
        ranked_signatures: set[str] = set()
        screen_rows: list[Mapping[str, str]] = []
        for role in ("single_knockout_moma", "single_knockout_room"):
            for row in csv_rows.get(role, ()):
                screen_rows.append(row)
                try:
                    display_rank = int(float(row.get("display_rank", "")))
                except (TypeError, ValueError):
                    continue
                if 1 <= display_rank <= 5:
                    signature = normalized_signature(
                        row.get("blocked_reaction_signature", "")
                    )
                    if signature:
                        ranked_signatures.add(signature)

        aliases: dict[str, set[str]] = {
            signature: set() for signature in ranked_signatures
        }
        for row in screen_rows:
            signature = normalized_signature(row.get("blocked_reaction_signature", ""))
            target = row.get("target_id", "").strip()
            if signature in aliases and target:
                aliases[signature].add(target)

        # The mapping also accounts for signature-equivalent genes that an older screen
        # may have collapsed before export.  Compare the complete blocked-reaction set per
        # gene so a multi-reaction deletion is never mistaken for a single-reaction alias.
        mapped_reactions: dict[str, set[str]] = {}
        for row in csv_rows.get("gene_knockout_mapping", ()):
            gene = row.get("gene_id", "").strip()
            reaction = row.get("blocked_reaction", "").strip()
            if gene and reaction:
                mapped_reactions.setdefault(gene, set()).add(reaction)
        for gene, reactions in mapped_reactions.items():
            signature = normalized_signature(";".join(reactions))
            if signature in aliases:
                aliases[signature].add(gene)
        return aliases

    expected_aliases = display_ranked_signature_aliases()
    expected_signatures = set(expected_aliases)

    def validate_candidate_aliases(
        rows: Sequence[Mapping[str, str]],
        *,
        assay: str,
    ) -> None:
        for row in rows:
            signature = normalized_signature(row.get("blocked_reaction_signature", ""))
            expected = expected_aliases.get(signature)
            if not expected:
                continue
            declared = {
                target.strip()
                for target in row.get("candidate_target_ids", "").split(";")
                if target.strip()
            }
            if declared != expected:
                issues.append(
                    f"{assay} candidate_target_ids for signature {signature!r} must "
                    "list every signature-equivalent gene; "
                    f"declared={sorted(declared)}, expected={sorted(expected)}"
                )

    exhaustive_response_ko = [
        row
        for row in response_index
        if row.get("candidate_scope", "") == "all_display_ranked_candidates"
    ]
    if exhaustive_response_ko:
        indexed_signatures = {
            normalized_signature(row.get("blocked_reaction_signature", ""))
            for row in exhaustive_response_ko
            if normalized_signature(row.get("blocked_reaction_signature", ""))
        }
        if not expected_signatures:
            issues.append(
                "exhaustive knockout flux-response policy was declared but display-ranked "
                "screen signatures were unavailable"
            )
        elif indexed_signatures != expected_signatures:
            issues.append(
                "exhaustive knockout flux-response signatures must match the D1-D5 union; "
                f"index-only={sorted(indexed_signatures - expected_signatures)}, "
                f"screen-only={sorted(expected_signatures - indexed_signatures)}"
            )
        validate_candidate_aliases(
            exhaustive_response_ko,
            assay="exhaustive knockout flux-response",
        )

    exhaustive_sampling = [
        row
        for row in sampling_index
        if row.get("target_id", "") != "wild_type"
        and row.get("candidate_scope", "") == "all_display_ranked_candidates"
    ]
    if exhaustive_sampling:
        indexed_signatures = {
            normalized_signature(row.get("blocked_reaction_signature", ""))
            for row in exhaustive_sampling
            if normalized_signature(row.get("blocked_reaction_signature", ""))
        }
        if not expected_signatures:
            issues.append(
                "exhaustive knockout sampling policy was declared but display-ranked "
                "screen signatures were unavailable"
            )
        elif indexed_signatures != expected_signatures:
            issues.append(
                "exhaustive knockout sampling signatures must match the D1-D5 union; "
                f"index-only={sorted(indexed_signatures - expected_signatures)}, "
                f"screen-only={sorted(expected_signatures - indexed_signatures)}"
            )
        validate_candidate_aliases(
            exhaustive_sampling,
            assay="exhaustive knockout sampling",
        )

    expected_amplification = {
        row.get("reaction_id", "").strip()
        for role in (
            "amplification_target_ranking",
            "variability_supported_amplification_targets",
        )
        for row in csv_rows.get(role, ())
        if _flag(row.get("report_selected")) is True
        and row.get("reaction_id", "").strip()
    }
    exhaustive_amplification = [
        row
        for row in response_index
        if row.get("candidate_scope", "") == "all_report_selected_candidates"
    ]
    if exhaustive_amplification:
        indexed_targets = {
            row.get("target", "").strip()
            for row in exhaustive_amplification
            if row.get("target", "").strip()
        }
        if not expected_amplification:
            issues.append(
                "exhaustive amplification flux-response policy was declared but no "
                "report-selected FSEOF/FVSEOF targets were available"
            )
        elif indexed_targets != expected_amplification:
            issues.append(
                "exhaustive amplification flux-response targets must match the independent "
                "report-selected union; "
                f"index-only={sorted(indexed_targets - expected_amplification)}, "
                f"ranking-only={sorted(expected_amplification - indexed_targets)}"
            )

    coverage = summary.get("validation_coverage")
    if isinstance(coverage, Mapping):
        response_counts = {
            "flux_response_expected": len(response_index),
            "flux_response_attempted": sum(
                row.get("status", "").strip().lower() in {"complete", "failed"}
                for row in response_index
            ),
            "flux_response_completed": sum(
                row.get("status", "").strip().lower() == "complete"
                for row in response_index
            ),
            "flux_response_failed": sum(
                row.get("status", "").strip().lower() == "failed"
                for row in response_index
            ),
        }
        sampling_counts = {
            "sampling_expected": len(sampling_index),
            "sampling_attempted": sum(
                row.get("status", "").strip().lower() in {"complete", "failed"}
                for row in sampling_index
            ),
            "sampling_completed": sum(
                row.get("status", "").strip().lower() == "complete"
                for row in sampling_index
            ),
            "sampling_failed": sum(
                row.get("status", "").strip().lower() == "failed"
                for row in sampling_index
            ),
            "sampling_skipped": sum(
                row.get("status", "").strip().lower() == "skipped"
                for row in sampling_index
            ),
        }
        indexed_counts = {
            **(response_counts if response_index else {}),
            **(sampling_counts if sampling_index else {}),
        }
        for column, actual in indexed_counts.items():
            declared = coverage.get(column)
            if declared is not None and declared != actual:
                issues.append(
                    f"summary validation_coverage.{column}={declared!r} does not match "
                    f"the candidate index count {actual}"
                )
        declared_knockouts = coverage.get("single_knockout_candidates_expected")
        if (
            exhaustive_response_ko
            and declared_knockouts is not None
            and declared_knockouts != len(expected_signatures)
        ):
            issues.append(
                "summary validation_coverage.single_knockout_candidates_expected does not "
                "match the D1-D5 signature union"
            )
        declared_amplification = coverage.get("amplification_candidates_expected")
        if (
            exhaustive_amplification
            and declared_amplification is not None
            and declared_amplification != len(expected_amplification)
        ):
            issues.append(
                "summary validation_coverage.amplification_candidates_expected does not "
                "match the report-selected amplification union"
            )
    loop_checked = {
        row.get("target", "")
        for row in csv_rows.get("amplification_loop_diagnostic", ())
    }
    for row in recommendations:
        target = row.get("target", "").strip()
        intervention_type = row.get("type", "").strip()
        verdict = row.get("verdict", "").strip().lower()
        if not target or intervention_type not in {
            "single_gene_knockout",
            "multi_knockout",
            "amplification",
        }:
            issues.append(
                f"recommendations contains an invalid target/type row: {target!r}/{intervention_type!r}"
            )
            continue
        if verdict and verdict not in {
            "support",
            "contradict",
            "inconclusive",
            "unavailable",
            "coupled",
            "uncoupled",
        }:
            issues.append(f"recommendation {target!r} has invalid verdict {verdict!r}")
        if (
            verdict in {"support", "contradict", "inconclusive"}
            and target not in responses
        ):
            issues.append(
                f"recommendation {target!r} has verdict {verdict!r} but no flux-response rows"
            )
        if (
            intervention_type == "single_gene_knockout"
            and verdict == "support"
            and target not in sampled
        ):
            issues.append(
                f"supported single-knockout recommendation {target!r} has no paired sampling rows"
            )
        if (
            intervention_type == "amplification"
            and verdict == "support"
            and target not in loop_checked
        ):
            warnings.append(
                f"supported amplification recommendation {target!r} has no loop diagnostic row"
            )


def _supplementary_records(
    manifest: Mapping[str, Any],
    *,
    root: Path,
    issues: list[str],
) -> tuple[list[Mapping[str, object]], dict[str, Mapping[str, object]]]:
    raw_records = manifest.get("supplementary_artifacts", [])
    if raw_records is None:
        raw_records = []
    if not isinstance(raw_records, list):
        issues.append("manifest supplementary_artifacts must be a list")
        return [], {}
    records: list[Mapping[str, object]] = []
    by_path: dict[str, Mapping[str, object]] = {}
    for index, raw in enumerate(raw_records):
        if not isinstance(raw, Mapping):
            issues.append(f"supplementary_artifacts entry {index} must be an object")
            continue
        path_value = raw.get("path")
        role_value = raw.get("role")
        role = (
            str(role_value)
            if isinstance(role_value, str)
            else f"supplementary[{index}]"
        )
        path = _resolve_artifact_path(root, path_value, role=role, issues=issues)
        if path is None or not isinstance(path_value, str):
            continue
        if path_value in by_path:
            issues.append(f"supplementary_artifacts repeats path {path_value!r}")
            continue
        by_path[path_value] = raw
        records.append(raw)
        if (
            path_value == "00_manifest.json"
            and role == "authoritative_artifact_manifest"
        ):
            continue
        if not path.is_file():
            issues.append(f"supplementary artifact {role!r} does not exist: {path}")
            continue
        if path.stat().st_size == 0:
            issues.append(f"supplementary artifact {role!r} is empty: {path}")
            continue
        _validate_integrity_fields(path=path, role=role, raw=raw, issues=issues)
    return records, by_path


def validate_run(
    run_dir: str | Path,
    *,
    manifest_name: str = "00_manifest.json",
) -> ValidatedRun:
    """Strictly validate and return a schema-v2 CMM publication input run.

    Required method roles must be complete. Optional roles may be omitted or explicitly marked
    partial/skipped/failed. Report and figure outputs are checked separately by
    :func:`validate_production_run`, allowing this function to gate the R renderer.
    """

    root = Path(run_dir).expanduser().resolve()
    issues: list[str] = []
    warnings: list[str] = []
    if not root.is_dir():
        raise RunValidationError([f"run directory does not exist: {root}"])
    manifest_path = root / manifest_name
    if not manifest_path.is_file():
        raise RunValidationError([f"required manifest is missing: {manifest_path}"])
    manifest_value = _read_json(manifest_path, label=manifest_name, issues=issues)
    if not isinstance(manifest_value, Mapping):
        raise RunValidationError(
            issues or [f"{manifest_name} must contain a JSON object"]
        )
    manifest: Mapping[str, Any] = manifest_value

    schema_version = manifest.get("schema_version")
    if type(schema_version) is not int or schema_version != RUN_SCHEMA_VERSION:
        issues.append(
            f"schema_version must be the integer {RUN_SCHEMA_VERSION}, got {schema_version!r}"
        )
    report = manifest.get("report")
    if not isinstance(report, Mapping):
        issues.append("manifest report must be an object")
    else:
        for report_field in ("title", "product_label"):
            if (
                not isinstance(report.get(report_field), str)
                or not str(report.get(report_field)).strip()
            ):
                issues.append(
                    f"manifest report.{report_field} must be a non-empty string"
                )
        if report.get("language", "en") != "en":
            issues.append("publication report language must be 'en'")

    raw_artifacts = manifest.get("artifacts")
    if not isinstance(raw_artifacts, Mapping):
        issues.append("manifest artifacts must be an object keyed by semantic role")
        raw_artifacts = {}

    supplementary, supplementary_by_path = _supplementary_records(
        manifest, root=root, issues=issues
    )
    supplementary_by_role: dict[str, list[Mapping[str, object]]] = {}
    for record in supplementary:
        role_value = record.get("role")
        if isinstance(role_value, str):
            supplementary_by_role.setdefault(role_value, []).append(record)

    artifacts: dict[str, RunArtifact] = {}
    csv_rows: dict[str, Sequence[Mapping[str, str]]] = {}
    json_objects: dict[str, Mapping[str, Any]] = {}
    for role, contract in ARTIFACT_CONTRACTS.items():
        raw = raw_artifacts.get(role)
        source: Literal["primary", "supplementary"] = "primary"
        if raw is None and len(supplementary_by_role.get(role, ())) == 1:
            raw = supplementary_by_role[role][0]
            source = "supplementary"
        elif raw is None and len(supplementary_by_role.get(role, ())) > 1:
            issues.append(
                f"supplementary_artifacts declares role {role!r} more than once"
            )
            continue
        if raw is None:
            if contract.required:
                issues.append(f"manifest is missing required artifact role {role!r}")
            continue
        artifact = _normalise_entry(
            role, raw, contract, root=root, issues=issues, source=source
        )
        artifacts[role] = artifact
        parsed = _validate_file(
            artifact,
            contract,
            raw,
            supplementary_by_path=supplementary_by_path,
            issues=issues,
            warnings=warnings,
        )
        if contract.kind == "csv" and isinstance(parsed, tuple):
            _, rows = parsed
            csv_rows[role] = rows
            _validate_numeric_rows(role, rows, issues=issues, warnings=warnings)
        elif contract.kind == "json" and isinstance(parsed, Mapping):
            json_objects[role] = parsed

    unknown = sorted(set(raw_artifacts) - set(ARTIFACT_CONTRACTS))
    for role in unknown:
        warnings.append(
            f"manifest declares unrecognised artifact role {role!r}; its integrity is not interpreted"
        )
    _validate_cross_artifact_invariants(
        csv_rows, json_objects, issues=issues, warnings=warnings
    )

    if issues:
        raise RunValidationError(issues)
    return ValidatedRun(
        root=root,
        manifest_path=manifest_path,
        manifest=manifest,
        artifacts=artifacts,
        warnings=tuple(warnings),
    )


def validate_production_run(
    run_dir: str | Path,
    *,
    manifest_name: str = "00_manifest.json",
) -> ValidationReport:
    """Validate core artifacts and publication outputs through the public stable alias."""

    # Delayed import avoids a module cycle: publication rendering itself depends on validate_run.
    from cmm.reporting.publication import validate_production_run as validate_outputs

    return validate_outputs(run_dir, manifest_name=manifest_name)
