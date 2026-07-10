"""Reference flux state primitive.

A :class:`FluxState` is a named, serializable flux vector with provenance. It is the
shared reference object for distance-based methods (MOMA, ROOM, revert-metabolism), so the
"distance from a reference distribution" logic lives in exactly one place.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
import math
from typing import Literal, cast

import pandas as pd
from cobra import Model
from cobra.flux_analysis import pfba

Provenance = Literal["fba", "pfba", "sampling_mean", "imported"]


@dataclass(frozen=True)
class FluxState:
    """A named reference flux distribution."""

    fluxes: Mapping[str, float]
    name: str = "reference"
    provenance: Provenance = "imported"
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Freeze the mapping into a plain dict of floats for hashing-free equality and
        # to reject accidental non-numeric payloads early.
        object.__setattr__(
            self,
            "fluxes",
            {str(rid): float(value) for rid, value in self.fluxes.items()},
        )
        object.__setattr__(self, "metadata", dict(self.metadata))
        self.validate()

    def validate(self) -> None:
        if not self.name:
            raise ValueError("FluxState name must not be empty")
        if not self.fluxes:
            raise ValueError("FluxState must contain at least one flux")
        if self.provenance not in {"fba", "pfba", "sampling_mean", "imported"}:
            raise ValueError(f"unknown FluxState provenance {self.provenance!r}")
        invalid = [
            rid for rid, value in self.fluxes.items() if not math.isfinite(value)
        ]
        if invalid:
            raise ValueError(
                f"FluxState contains non-finite fluxes: {', '.join(invalid[:5])}"
            )

    def get(self, reaction_id: str, default: float = 0.0) -> float:
        return self.fluxes.get(reaction_id, default)

    def reactions(self) -> tuple[str, ...]:
        return tuple(self.fluxes.keys())

    def to_series(self) -> pd.Series:
        return pd.Series(self.fluxes, name=self.name, dtype=float).sort_index()

    def distance(
        self,
        other: FluxState | Mapping[str, float],
        reactions: Iterable[str] | None = None,
        order: Literal[1, 2] = 2,
    ) -> float:
        """Distance to another flux vector over the union (or a subset) of reactions."""

        other_map = other.fluxes if isinstance(other, FluxState) else dict(other)
        if reactions is None:
            keys = set(self.fluxes) | set(other_map)
        else:
            keys = set(reactions)
        if order == 1:
            return float(sum(abs(self.get(k) - other_map.get(k, 0.0)) for k in keys))
        if order == 2:
            return (
                float(sum((self.get(k) - other_map.get(k, 0.0)) ** 2 for k in keys))
                ** 0.5
            )
        raise ValueError("order must be 1 or 2")

    def serialize(self) -> dict:
        return {
            "name": self.name,
            "provenance": self.provenance,
            "fluxes": dict(self.fluxes),
            "metadata": dict(self.metadata),
        }

    @classmethod
    def deserialize(cls, payload: Mapping[str, object]) -> FluxState:
        raw_fluxes = payload.get("fluxes", {})
        raw_metadata = payload.get("metadata", {})
        if not isinstance(raw_fluxes, Mapping):
            raise ValueError("serialized FluxState 'fluxes' must be a mapping")
        if not isinstance(raw_metadata, Mapping):
            raise ValueError("serialized FluxState 'metadata' must be a mapping")
        raw_provenance = str(payload.get("provenance", "imported"))
        if raw_provenance not in {"fba", "pfba", "sampling_mean", "imported"}:
            raise ValueError(f"unknown FluxState provenance {raw_provenance!r}")
        return cls(
            fluxes={str(key): float(value) for key, value in raw_fluxes.items()},
            name=str(payload.get("name", "reference")),
            provenance=cast(Provenance, raw_provenance),
            metadata={str(key): value for key, value in raw_metadata.items()},
        )

    @classmethod
    def from_solution(
        cls, solution, name: str = "reference", provenance: Provenance = "fba"
    ) -> FluxState:
        if getattr(solution, "status", None) != "optimal":
            raise ValueError(
                f"cannot create FluxState from solution status {getattr(solution, 'status', None)!r}"
            )
        return cls(fluxes=dict(solution.fluxes), name=name, provenance=provenance)


def reference_state_pfba(
    model: Model,
    condition=None,
    name: str = "reference",
    fraction_of_optimum: float = 1.0,
) -> FluxState:
    """Build a reference state from a parsimonious FBA (pFBA) solution.

    pFBA yields a unique, minimal-total-flux distribution at the given fraction of the
    optimum, which makes it a reproducible default reference for small models.
    """

    if not 0.0 < fraction_of_optimum <= 1.0:
        raise ValueError("fraction_of_optimum must be in (0, 1]")
    with model:
        if condition is not None:
            condition.apply_to(model)
        solution = pfba(model, fraction_of_optimum=fraction_of_optimum)
    return FluxState.from_solution(solution, name=name, provenance="pfba")


def reference_state_from_samples(
    samples: pd.DataFrame,
    name: str = "reference",
) -> FluxState:
    """Build a reference state from the column means of a flux-sampling table."""

    if samples.empty:
        raise ValueError("sample table is empty")
    numeric = samples.apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        raise ValueError("sample table contains non-numeric or non-finite values")
    means = numeric.mean(axis=0)
    return FluxState(
        fluxes={str(rid): float(v) for rid, v in means.items()},
        name=name,
        provenance="sampling_mean",
        metadata={"n_samples": int(len(samples))},
    )
