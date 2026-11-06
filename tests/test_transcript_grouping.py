from __future__ import annotations

import pytest
from cmm.core.flux_state import reference_state_pfba
from cmm.features import knockout_mechanism, revert_candidates, transformation_targets
from cmm.omics.differential import DirectionMap


def test_one_isoform_blocks_nothing_but_the_grouped_gene_is_a_candidate(
    transcript_model,
):
    per_transcript, record = revert_candidates(transcript_model)
    grouped, grouped_record = revert_candidates(
        transcript_model, transcript_separator="."
    )
    # Each isoform alone is inert, so per-transcript screening never tests gene t1 at all.
    assert per_transcript == ("t3.1",)
    assert grouped == ("t1", "t3")
    assert grouped_record["transcript_separator"] == "."
    # t2 blocks the only route to biomass, under either enumeration.
    assert record["n_genes_essential_removed"] == 1
    assert grouped_record["n_genes_essential_removed"] == 1


def test_transcript_grouping_is_refused_for_reaction_knockouts(transcript_model):
    with pytest.raises(ValueError, match="gene knockouts only"):
        revert_candidates(
            transcript_model, perturbation="reaction", transcript_separator="."
        )


@pytest.mark.requires_qp
def test_mechanism_knocks_out_every_transcript_of_a_grouped_gene(transcript_model):
    reference = reference_state_pfba(transcript_model, name="source")
    direction = DirectionMap({"R1": -1, "R3": 1})
    mechanism = knockout_mechanism(
        transcript_model,
        reference,
        direction,
        ["t1", "t1.1"],
        transcript_separator=".",
    )
    by_id = {e.target_id: e for e in mechanism.explanations}
    assert by_id["t1"].blocked_reactions == ("R1",)
    # A transcript id is not a knockout unit once transcripts are grouped.
    assert by_id["t1.1"].status == "not_a_knockout"


@pytest.mark.requires_qp
def test_moma_baseline_scores_grouped_genes(transcript_model):
    source = reference_state_pfba(transcript_model, name="source")
    ranking = transformation_targets(
        transcript_model, source, source, transcript_separator="."
    )
    assert {t.target_id for t in ranking} == {"t1", "t2", "t3"}
