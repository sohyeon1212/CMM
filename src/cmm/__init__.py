"""CMM public package interface."""

from importlib.metadata import PackageNotFoundError, version as _version

from cmm.core.condition import Condition, ObjectiveSpec, ReactionBound
from cmm.core.simulation import FluxRange, FluxSolution, fba, fva

__all__ = [
    "Condition",
    "FluxRange",
    "FluxSolution",
    "ObjectiveSpec",
    "ReactionBound",
    "fba",
    "fva",
]

# Single source of truth is pyproject.toml; read it from the installed package metadata so the
# version never drifts from the release. The literal fallback is only used when running from a
# source tree that was never installed (no dist metadata to read).
try:
    __version__ = _version("cmm")
except PackageNotFoundError:  # pragma: no cover - source tree without an install
    __version__ = "0.3.0"
