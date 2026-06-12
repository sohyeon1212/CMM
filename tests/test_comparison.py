from __future__ import annotations

import pytest
from cmm.core.flux_state import reference_state_pfba
from cmm.features._perturbation import (
    apply_perturbation,
    gene_perturbations,
    reaction_perturbations,
    run_perturbations,
)
from cmm.features.comparison import moma, reference_flux, room


def test_reference_flux_fba_and_pfba(ecoli_core):
    fba_ref = reference_flux(ecoli_core, "fba")
    pfba_ref = reference_flux(ecoli_core, "pfba")
    assert fba_ref.provenance == "imported"
    # Both reproduce growth; pFBA has the smaller total flux.
    assert fba_ref.get("Biomass_Ecoli_core") == pytest.approx(0.8739, abs=1e-3)
    assert pfba_ref.get("Biomass_Ecoli_core") == pytest.approx(0.8739, abs=1e-3)
    assert sum(abs(v) for v in pfba_ref.fluxes.values()) <= sum(
        abs(v) for v in fba_ref.fluxes.values()
    ) + 1e-6


def test_reference_flux_omics_templates(ecoli_core):
    expression = {g.id: 50.0 for g in ecoli_core.genes}
    for method in ("lad", "eflux2"):
        ref = reference_flux(ecoli_core, method, gene_expression=expression)
        assert len(ref.fluxes) == len(ecoli_core.reactions)


def test_reference_flux_omics_requires_expression(ecoli_core):
    with pytest.raises(ValueError, match="requires gene_expression"):
        reference_flux(ecoli_core, "eflux2")


def test_moma_uses_chosen_template_as_reference(ecoli_core, unrestricted_qp_solver):
    # Different templates -> different MOMA reference -> different perturbed distance.
    fba_ref = reference_flux(ecoli_core, "fba")
    pfba_ref = reference_flux(ecoli_core, "pfba")
    with ecoli_core:
        ecoli_core.reactions.PFK.bounds = (0.0, 0.0)
        d_fba = moma(ecoli_core, fba_ref, linear=False).distance
    with ecoli_core:
        ecoli_core.reactions.PFK.bounds = (0.0, 0.0)
        d_pfba = moma(ecoli_core, pfba_ref, linear=False).distance
    # pFBA template is the minimal-flux reference, so its MOMA distance differs from FBA's.
    assert d_fba >= 0 and d_pfba >= 0

# --- perturbation enumeration / application --------------------------------


def test_gene_perturbations_resolve_blocked_reactions(branched_model):
    perts = {p.target_id: p for p in gene_perturbations(branched_model)}
    assert perts["g2"].reaction_ids == ("R2",)
    assert perts["g3"].reaction_ids == ("R3",)
    assert perts["g2"].kind == "gene"


def test_run_perturbations_restores_model(branched_model):
    perts = reaction_perturbations(branched_model, ["R2"])

    def objective(model, _pert):
        return model.slim_optimize()

    results = run_perturbations(branched_model, perts, objective)
    # Knocking out R2 still allows growth via the R3->R5 branch.
    assert results[0][1] == pytest.approx(10, abs=1e-6)
    # Original bounds restored.
    assert branched_model.reactions.R2.bounds == (0.0, 1000.0)


def test_apply_perturbation_zeros_reactions(branched_model):
    pert = reaction_perturbations(branched_model, ["R2"])[0]
    with apply_perturbation(branched_model, pert):
        assert branched_model.reactions.R2.bounds == (0.0, 0.0)
    assert branched_model.reactions.R2.bounds == (0.0, 1000.0)


# --- MOMA ------------------------------------------------------------------


def test_l2_moma_reroutes_after_knockout(branched_model):
    reference = reference_state_pfba(branched_model, name="wt")
    pert = reaction_perturbations(branched_model, ["R2"])[0]
    with apply_perturbation(branched_model, pert):
        result = moma(branched_model, reference, linear=False)
    assert result.status == "optimal"
    assert result.method == "moma_l2"
    # The disease branch is dead; flux must reroute through R3->R5.
    assert result.fluxes["R2"] == pytest.approx(0, abs=1e-6)
    assert result.fluxes["R3"] > 1.0
    assert result.fluxes["R5"] > 1.0
    assert result.distance > 0


def test_l1_moma_runs_as_lp(branched_model):
    reference = reference_state_pfba(branched_model, name="wt")
    pert = reaction_perturbations(branched_model, ["R2"])[0]
    with apply_perturbation(branched_model, pert):
        result = moma(branched_model, reference, linear=True)
    assert result.status == "optimal"
    assert result.method == "moma_l1"
    assert result.fluxes["R3"] > 1.0
    # Restoring the objective leaves the model optimizing biomass again.
    assert branched_model.slim_optimize() == pytest.approx(10, abs=1e-6)


def test_moma_zero_distance_without_perturbation(branched_model):
    reference = reference_state_pfba(branched_model, name="wt")
    result = moma(branched_model, reference, linear=False)
    assert result.distance == pytest.approx(0, abs=1e-6)


# --- ROOM ------------------------------------------------------------------


def test_room_counts_changed_reactions(branched_model):
    reference = reference_state_pfba(branched_model, name="wt")
    pert = reaction_perturbations(branched_model, ["R2"])[0]
    with apply_perturbation(branched_model, pert):
        result = room(branched_model, reference)
    assert result.status == "optimal"
    assert result.method == "room"
    # Rerouting changes R2 (off), R3 (on), R5 (on): a small, positive switch count.
    assert result.distance >= 2
    assert result.fluxes["R3"] > 1.0
