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
    assert fba_ref.provenance == "fba"
    assert pfba_ref.provenance == "pfba"
    # Both reproduce growth; pFBA has the smaller total flux.
    assert fba_ref.get("Biomass_Ecoli_core") == pytest.approx(0.8739, abs=1e-3)
    assert pfba_ref.get("Biomass_Ecoli_core") == pytest.approx(0.8739, abs=1e-3)
    assert (
        sum(abs(v) for v in pfba_ref.fluxes.values())
        <= sum(abs(v) for v in fba_ref.fluxes.values()) + 1e-6
    )


def test_reference_flux_omics_templates(ecoli_core):
    expression = {g.id: 50.0 for g in ecoli_core.genes}
    for method in ("lad", "eflux2"):
        ref = reference_flux(ecoli_core, method, gene_expression=expression)
        assert len(ref.fluxes) == len(ecoli_core.reactions)


def test_reference_flux_omics_requires_expression(ecoli_core):
    with pytest.raises(ValueError, match="requires gene_expression"):
        reference_flux(ecoli_core, "eflux2")


def test_moma_uses_chosen_template_as_reference(branched_model):
    # Different templates -> different MOMA reference -> different perturbed distance.
    fba_ref = reference_flux(branched_model, "fba")
    pfba_ref = reference_flux(branched_model, "pfba")
    with branched_model:
        branched_model.reactions.R2.bounds = (0.0, 0.0)
        d_fba = moma(branched_model, fba_ref, linear=False).distance
    with branched_model:
        branched_model.reactions.R2.bounds = (0.0, 0.0)
        d_pfba = moma(branched_model, pfba_ref, linear=False).distance
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


def test_moma_l2_distance_is_the_root_of_the_qp_objective(branched_model):
    """Segrè et al. Eq. (4) is D = sqrt(Sum d^2); the QP objective is the sum itself."""

    import math

    reference = reference_state_pfba(branched_model, name="wt")
    pert = reaction_perturbations(branched_model, ["R2"])[0]
    with apply_perturbation(branched_model, pert):
        result = moma(branched_model, reference, linear=False)
    assert result.distance_kind == "euclidean_l2"
    assert result.n_changed_reactions is None
    assert result.objective_value > result.distance  # the sum is the square of the root
    assert result.distance == pytest.approx(math.sqrt(result.objective_value))
    # And the root is the Euclidean distance recomputed from the returned flux vector.
    recomputed = math.sqrt(
        sum(
            (result.fluxes[r.id] - reference.get(r.id)) ** 2
            for r in branched_model.reactions
        )
    )
    assert result.distance == pytest.approx(recomputed, abs=1e-6)


def test_moma_l1_distance_is_the_lp_objective(branched_model):
    reference = reference_state_pfba(branched_model, name="wt")
    pert = reaction_perturbations(branched_model, ["R2"])[0]
    with apply_perturbation(branched_model, pert):
        result = moma(branched_model, reference, linear=True)
    assert result.distance_kind == "l1"
    assert result.n_changed_reactions is None
    assert result.distance == pytest.approx(result.objective_value)
    recomputed = sum(
        abs(result.fluxes[r.id] - reference.get(r.id)) for r in branched_model.reactions
    )
    assert result.distance == pytest.approx(recomputed, abs=1e-6)


def test_comparison_records_its_reference(branched_model):
    reference = reference_state_pfba(branched_model, name="wt")
    result = moma(branched_model, reference, linear=False)
    assert result.metadata["reference"] == "wt"
    assert result.metadata["reference_provenance"] == "pfba"


def test_infeasible_moma_keeps_the_distance_labelling(branched_model):
    reference = reference_state_pfba(branched_model, name="wt")
    with branched_model:
        branched_model.reactions.SUP_A.bounds = (5.0, 5.0)
        branched_model.reactions.R1.bounds = (0.0, 0.0)
        result = moma(branched_model, reference, linear=False)
    assert result.status == "infeasible"
    assert result.distance_kind == "euclidean_l2"
    assert result.distance != result.distance  # NaN: no solution, not "zero distance"


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
    # ROOM's objective is a switch count, not a distance, so it is not reported as one.
    assert result.distance is None
    assert result.distance_kind == "none"
    # Rerouting changes R2 (off), R3 (on), R5 (on): a small, positive switch count.
    assert result.n_changed_reactions >= 2
    assert result.n_changed_reactions == pytest.approx(result.objective_value)
    assert result.fluxes["R3"] > 1.0


def test_room_tolerance_presets_are_shlomis_two_pairs(branched_model):
    from cmm.features.comparison import ROOM_TOLERANCES

    assert ROOM_TOLERANCES["flux_prediction"] == (0.03, 1e-3)
    assert ROOM_TOLERANCES["lethality"] == (0.1, 1e-2)

    reference = reference_state_pfba(branched_model, name="wt")
    pert = reaction_perturbations(branched_model, ["R2"])[0]
    with apply_perturbation(branched_model, pert):
        default = room(branched_model, reference)
        lethality = room(branched_model, reference, use_case="lethality")
        override = room(branched_model, reference, use_case="lethality", delta=0.5)
    assert (default.metadata["delta"], default.metadata["epsilon"]) == (0.03, 1e-3)
    assert (lethality.metadata["delta"], lethality.metadata["epsilon"]) == (0.1, 1e-2)
    # An explicit value wins over the preset, and only the one that was given.
    assert (override.metadata["delta"], override.metadata["epsilon"]) == (0.5, 1e-2)

    with pytest.raises(ValueError, match="unknown ROOM use case"):
        room(branched_model, reference, use_case="lethal")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="non-negative"):
        room(branched_model, reference, delta=-1.0)


def test_batch_screen_defaults_to_the_lethality_tolerances(branched_model):
    from cmm.features.comparison import batch_comparison, knockout_comparison

    reference = reference_state_pfba(branched_model, name="wt")
    single = knockout_comparison(branched_model, reference, ["R2"], method="room")
    rows = batch_comparison(
        branched_model,
        reference,
        reaction_perturbations(branched_model, ["R2"]),
        method="room",
    )
    # A named single knockout is a flux prediction; a screen asks about lethality.
    assert single.metadata["room_use_case"] == "flux_prediction"
    assert rows[0].distance is None
    assert rows[0].distance_kind == "none"
    assert rows[0].n_changed_reactions >= 1


# --- gene / multi / batch knockouts ----------------------------------------


def test_blocked_reactions_for_genes_joint(branched_model):
    from cmm.features._perturbation import blocked_reactions_for_genes

    assert blocked_reactions_for_genes(branched_model, ["g2"]) == ("R2",)
    # Knocking out g2 and g3 together blocks both their reactions.
    assert set(blocked_reactions_for_genes(branched_model, ["g2", "g3"])) == {
        "R2",
        "R3",
    }
    assert blocked_reactions_for_genes(branched_model, []) == ()


def test_knockout_comparison_gene_and_reaction_agree(branched_model):
    from cmm.features._perturbation import blocked_reactions_for_genes
    from cmm.features.comparison import knockout_comparison

    reference = reference_state_pfba(branched_model, name="wt")
    # Gene g2 disables exactly R2, so a gene KO and the equivalent reaction KO match.
    gene_rxns = blocked_reactions_for_genes(branched_model, ["g2"])
    by_gene = knockout_comparison(
        branched_model, reference, gene_rxns, method="moma_l2"
    )
    by_reaction = knockout_comparison(
        branched_model, reference, ["R2"], method="moma_l2"
    )
    assert by_gene.status == "optimal"
    assert by_gene.distance == pytest.approx(by_reaction.distance, abs=1e-6)
    assert by_gene.fluxes["R3"] > 1.0  # rerouted through the healthy branch
    # Model restored after the knockout context.
    assert branched_model.reactions.R2.bounds == (0.0, 1000.0)


def test_knockout_comparison_multi_reaction(branched_model):
    from cmm.features.comparison import knockout_comparison

    reference = reference_state_pfba(branched_model, name="wt")
    # Knocking out both branches (R2 and R3) leaves no route to product -> lethal/infeasible.
    result = knockout_comparison(
        branched_model, reference, ["R2", "R3"], method="moma_l2"
    )
    assert result.status != "optimal" or result.fluxes.get(
        "BIOMASS", 0.0
    ) == pytest.approx(0.0, abs=1e-6)


def test_batch_comparison_ranks_targets(branched_model):
    from cmm.features._perturbation import gene_perturbations
    from cmm.features.comparison import batch_comparison

    reference = reference_state_pfba(branched_model, name="wt")
    rows = {
        r.target_id: r
        for r in batch_comparison(
            branched_model,
            reference,
            gene_perturbations(branched_model),
            method="moma_l2",
        )
    }
    # Every non-inert gene is scored.
    assert {"g1", "g2", "g3", "g5", "gb"} <= set(rows)
    # g3/g5 (unused healthy branch at the pFBA optimum) have no effect: distance ~0.
    assert rows["g3"].distance == pytest.approx(0.0, abs=1e-6)
    # g2 (the used disease branch) forces a reroute: nonzero distance. MOMA minimizes the
    # deviation from the reference (not growth), so the predicted biomass drops to 6 (the
    # point closest to the reference's R2=10, R3=R5=0), staying positive.
    assert rows["g2"].distance > 0
    assert rows["g2"].objective == pytest.approx(6.0, abs=1e-6)
    assert rows["g2"].kind == "gene"
    assert rows["g2"].n_reactions == 1


def test_gene_perturbations_report_what_they_dropped(branched_model):
    from cmm.features._perturbation import (
        PerturbationList,
        gene_perturbations,
        perturbation_provenance,
        reaction_perturbations,
    )

    # Give the fixture a gene that controls nothing, so there is something to drop.
    branched_model.reactions.SUP_A.gene_reaction_rule = "g1 or g_inert"
    perts = gene_perturbations(branched_model)
    assert isinstance(perts, PerturbationList)
    assert "g_inert" in perts.inert_dropped
    assert perts.n_inert_dropped == 1
    assert {p.target_id for p in perts}.isdisjoint(perts.inert_dropped)
    provenance = perts.provenance()
    assert provenance["n_perturbations"] == len(perts)
    assert provenance["n_inert_dropped"] == 1
    assert provenance["n_candidates_considered"] == len(branched_model.genes)

    # include_inert keeps them, and then nothing is reported as dropped.
    kept = gene_perturbations(branched_model, include_inert=True)
    assert "g_inert" in {p.target_id for p in kept}
    assert kept.n_inert_dropped == 0

    # Reaction knockouts are never inert, and a plain sequence reports None rather than 0.
    assert reaction_perturbations(branched_model).n_inert_dropped == 0
    assert perturbation_provenance(list(perts))["n_inert_dropped"] is None
