"""Schema-v2 run-bundle writing, shared by CMM's canonical workflows.

Extracted verbatim from :mod:`cmm.workflows.production` so a second workflow can produce the
same bundle without duplicating the format. ``00_manifest.json`` being the only path-discovery
surface is what lets one validator serve every workflow; two writers would be two contracts.

The only change made while moving is that the path-escape guard raises an injectable error
type. :mod:`cmm.workflows.production` passes its own ``ProductionWorkflowError``, so the
exception a production run raises is unchanged.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, is_dataclass, replace
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from cmm.reporting.schema import ArtifactStatus


def _jsonable(value: object) -> object:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value) and not isinstance(value, type):
        return _jsonable(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_jsonable(item) for item in value]
    return str(value)


class BundleError(RuntimeError):
    """A run bundle could not be written as specified."""


@dataclass(frozen=True)
class ArtifactRecord:
    path: str
    stage: str
    role: str
    media_type: str
    status: ArtifactStatus = "complete"
    method: str | None = None
    reason: str | None = None
    metadata_path: str | None = None
    sha256: str | None = None
    size_bytes: int | None = None


class _ArtifactWriter:
    def __init__(
        self, root: Path, *, error_type: type[Exception] = BundleError
    ) -> None:
        self.root = root.resolve()
        self._error_type = error_type
        self.records: list[ArtifactRecord] = []
        self.metadata_links: dict[str, str] = {}

    def _path(self, relative: str) -> Path:
        path = self.root / relative
        resolved = path.resolve(strict=False)
        if not resolved.is_relative_to(self.root):
            raise self._error_type(
                f"artifact path escapes the run directory: {relative!r}"
            )
        return path

    def csv(
        self,
        relative: str,
        frame: pd.DataFrame,
        *,
        stage: str,
        role: str,
        method: str | None = None,
        status: ArtifactStatus = "complete",
        reason: str | None = None,
        metadata_path: str | None = None,
    ) -> None:
        path = self._path(relative)
        frame.to_csv(path, index=False)
        self._record(
            path,
            stage,
            role,
            "text/csv",
            method,
            status,
            reason,
            metadata_path,
        )

    def json(
        self,
        relative: str,
        payload: object,
        *,
        stage: str,
        role: str,
        status: ArtifactStatus = "complete",
        reason: str | None = None,
    ) -> None:
        path = self._path(relative)
        path.write_text(
            json.dumps(_jsonable(payload), indent=2, sort_keys=True, ensure_ascii=False)
            + "\n",
            encoding="utf-8",
        )
        self._record(path, stage, role, "application/json", None, status, reason, None)

    def metadata(
        self,
        role: str,
        relative: str,
        payload: object,
        *,
        stage: str,
    ) -> None:
        self.json(
            relative,
            payload,
            stage=stage,
            role=f"{role}_metadata",
        )
        self.metadata_links[role] = relative
        self.records = [
            replace(record, metadata_path=relative)
            if record.role == role and record.media_type == "text/csv"
            else record
            for record in self.records
        ]

    def text(
        self,
        relative: str,
        content: str,
        *,
        stage: str,
        role: str,
        media_type: str = "text/x-python",
        executable: bool = False,
    ) -> None:
        path = self._path(relative)
        path.write_text(content, encoding="utf-8")
        if executable:
            path.chmod(0o755)
        self._record(path, stage, role, media_type, None, "complete", None, None)

    def existing(
        self,
        relative: str,
        *,
        stage: str,
        role: str,
        media_type: str,
    ) -> None:
        self._record(
            self._path(relative),
            stage,
            role,
            media_type,
            None,
            "complete",
            None,
            None,
        )

    def _record(
        self,
        path: Path,
        stage: str,
        role: str,
        media_type: str,
        method: str | None,
        status: ArtifactStatus,
        reason: str | None,
        metadata_path: str | None,
    ) -> None:
        payload = path.read_bytes()
        self.records.append(
            ArtifactRecord(
                path=path.relative_to(self.root).as_posix(),
                stage=stage,
                role=role,
                media_type=media_type,
                status=status,
                method=method,
                reason=reason,
                metadata_path=metadata_path,
                sha256=hashlib.sha256(payload).hexdigest(),
                size_bytes=len(payload),
            )
        )
