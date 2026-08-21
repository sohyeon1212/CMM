"""Data files CMM ships with, and the rules for when one applies to a model.

Only one kind of resource lives here today: a curated Escher pathway map. CMM does not
generate map layouts — a readable metabolic map is hand-drawn, and an automatic layout of a
genome-scale network is a hairball rather than a figure. Escher's published maps are the
layouts the field already reads, so CMM renders those and bundles the *E. coli* core one.

Provenance and license for every bundled file are recorded in ``ATTRIBUTION.md`` beside this
module; the map is redistributed byte-for-byte under Escher's MIT license.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - typing only
    from cobra import Model

__all__ = ["BUNDLED_MAPS", "bundled_map_for", "map_reaction_ids", "MAP_COVERAGE_MIN"]

_HERE = Path(__file__).parent

#: Curated Escher maps shipped with CMM, most specific first.
BUNDLED_MAPS: tuple[Path, ...] = (_HERE / "e_coli_core.Core metabolism.json",)

#: A map is offered for a model when it can draw at least this fraction of its own reactions.
#: A central-metabolism map is a legitimate view of a genome-scale model — that is how Escher
#: itself is used — so the test is how much of the *map* the model can fill, never how much of
#: the model the map covers.
MAP_COVERAGE_MIN = 0.5


@lru_cache(maxsize=8)
def map_reaction_ids(path: str) -> frozenset[str]:
    """Return the BiGG reaction ids an Escher map JSON draws."""

    with open(path, encoding="utf-8") as handle:
        document = json.load(handle)
    # An Escher map is [metadata, {"reactions": {...}, "nodes": {...}, ...}].
    reactions = document[1].get("reactions", {})
    return frozenset(
        entry["bigg_id"] for entry in reactions.values() if entry.get("bigg_id")
    )


def map_coverage(path: Path | str, model: Model) -> float:
    """Fraction of the map's reactions this model actually contains."""

    drawn = map_reaction_ids(str(path))
    if not drawn:
        return 0.0
    present = {rxn.id for rxn in model.reactions}
    return len(drawn & present) / len(drawn)


def bundled_map_for(model: Model) -> str | None:
    """Best bundled Escher map for ``model``, or ``None`` when none is a reasonable view.

    Returns a path string rather than a ``Path`` so callers can pass it straight to
    ``escher_flux_map`` and to Qt, both of which take either.
    """

    best: tuple[float, Path] | None = None
    for candidate in BUNDLED_MAPS:
        if not candidate.exists():  # pragma: no cover - defensive; wheels include it
            continue
        coverage = map_coverage(candidate, model)
        if coverage >= MAP_COVERAGE_MIN and (best is None or coverage > best[0]):
            best = (coverage, candidate)
    return str(best[1]) if best else None
