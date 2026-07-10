"""Sensitivity checks for manuscript-facing scientific conclusions."""

from __future__ import annotations

import pytest

from cmm.core.flux_state import FluxState
from cmm.features.production import fseof
from cmm.features.revert import revert_targets
from cmm.omics.differential import DirectionMap
from cmm.omics.expression import eflux2


def _official_mta_inputs():
    reference = FluxState(
        {"r1": 10, "r2": 1, "r3": 8, "r4": 8, "r5": 1, "r6": 9, "r7": 10}
    )
    direction = DirectionMap(
        {"r1": 0, "r2": 1, "r3": 0, "r4": -1, "r5": 0, "r6": 0, "r7": 0}
    )
    return reference, direction


@pytest.mark.parametrize("alpha", [0.3, 0.4, 0.66, 0.8])
def test_published_rmta_positive_control_is_stable_across_alpha(
    published_mta_model, alpha
):
    reference, direction = _official_mta_inputs()
    ranking = revert_targets(
        published_mta_model,
        None,
        reference,
        direction,
        method="rmta",
        alpha=alpha,
        epsilon=0.01,
        transcript_separator=".",
        targets=["g2", "g4"],
    )
    scores = {target.target_id: target.score for target in ranking}
    assert scores["g4"] > 0
    assert scores["g4"] > scores["g2"]


@pytest.mark.parametrize("n_steps", [6, 8, 10, 12])
def test_fseof_headline_targets_are_stable_across_scan_resolution(ecoli_core, n_steps):
    result = fseof(
        ecoli_core,
        "EX_succ_e",
        "Biomass_Ecoli_core",
        n_steps=n_steps,
        aerobic=False,
    )
    targets = set(result.amplification_targets())
    assert {"FRD7", "FUM", "PPC"} <= targets


def test_eflux2_is_invariant_to_global_expression_units(ecoli_core):
    unit = {reaction.id: 1.0 for reaction in ecoli_core.reactions if reaction.genes}
    scaled = {reaction_id: 1000.0 for reaction_id in unit}
    first = eflux2(ecoli_core, unit)
    second = eflux2(ecoli_core, scaled)
    assert first.status == second.status == "optimal"
    for rid in ("Biomass_Ecoli_core", "PFK", "PYK", "ATPM"):
        assert first.fluxes[rid] == pytest.approx(second.fluxes[rid], abs=1e-6)
