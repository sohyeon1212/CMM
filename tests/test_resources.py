"""The bundled Escher map: that it is intact, attributed, and offered to the right models."""

from __future__ import annotations

import hashlib
import json
import re

import cobra
from cmm.resources import (
    BUNDLED_MAPS,
    bundled_map_for,
    map_coverage,
    map_reaction_ids,
)

CORE_MAP = BUNDLED_MAPS[0]
ATTRIBUTION = CORE_MAP.parent / "ATTRIBUTION.md"


def test_bundled_map_is_byte_identical_to_its_attributed_digest():
    """The file must stay exactly what ATTRIBUTION.md says CMM redistributes.

    CMM redistributes someone else's work under their license. If the bytes drift from the
    recorded digest the attribution silently becomes false, so the digest is a test, not a
    comment.
    """

    recorded = re.search(r"`([0-9a-f]{64})`", ATTRIBUTION.read_text())
    assert recorded, "ATTRIBUTION.md records no SHA-256 for the bundled map"
    actual = hashlib.sha256(CORE_MAP.read_bytes()).hexdigest()
    assert actual == recorded.group(1)


def test_attribution_names_the_license_and_the_paper_to_cite():
    text = ATTRIBUTION.read_text()
    assert "MIT" in text
    assert "Regents of the University of California" in text
    assert "10.1371/journal.pcbi.1004321" in text  # King et al. 2015, the Escher paper
    assert "escher.github.io" in text


def test_map_is_a_schema_1_escher_document():
    document = json.loads(CORE_MAP.read_text(encoding="utf-8"))
    assert document[0]["schema"].endswith("1-0-0#")
    assert len(map_reaction_ids(str(CORE_MAP))) == 95


def test_core_map_is_offered_for_the_core_model(ecoli_core):
    assert map_coverage(CORE_MAP, ecoli_core) > 0.95
    assert bundled_map_for(ecoli_core) == str(CORE_MAP)


def test_no_map_is_offered_for_an_unrelated_model():
    """A map of a different organism must not be drawn as if it described this model."""

    model = cobra.Model("unrelated")
    a, b = cobra.Metabolite("a_c"), cobra.Metabolite("b_c")
    reaction = cobra.Reaction("MADE_UP_RXN")
    reaction.add_metabolites({a: -1, b: 1})
    model.add_reactions([reaction])

    assert map_coverage(CORE_MAP, model) == 0.0
    assert bundled_map_for(model) is None
