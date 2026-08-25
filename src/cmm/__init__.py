"""CMM public package interface."""

from importlib.metadata import PackageNotFoundError, version as _version

from cmm.core.condition import Condition, ObjectiveSpec, ReactionBound
from cmm.core.media import Medium, MediumApplication
from cmm.core.simulation import FluxRange, FluxSolution, FvaResult, fba, fva, pfba
from cmm.features.response import ResponseLimit, ResponsePhase

#: The top level carries the result *types* a caller has to name to annotate or unpack a
#: result — the services themselves live in ``cmm.core``, ``cmm.features`` and ``cmm.omics``,
#: which is the boundary ``docs/agent-reference.md`` documents. 0.4.0 adds ``FvaResult``
#: (``fva`` no longer returns a bare dict), ``MediumApplication`` (``Medium.apply_to`` now
#: reports what it dropped) and ``ResponsePhase``/``ResponseLimit`` (which replace the deleted
#: ``ResponseBottleneck``).
__all__ = [
    "Condition",
    "FluxRange",
    "FluxSolution",
    "FvaResult",
    "Medium",
    "MediumApplication",
    "ObjectiveSpec",
    "ReactionBound",
    "ResponseLimit",
    "ResponsePhase",
    "fba",
    "fva",
    "pfba",
]

# Single source of truth is pyproject.toml; read it from the installed package metadata so the
# version never drifts from the release. The literal fallback is only used when running from a
# source tree that was never installed (no dist metadata to read).
try:
    __version__ = _version("cmm")
except PackageNotFoundError:  # pragma: no cover - source tree without an install
    __version__ = "0.5.0"
