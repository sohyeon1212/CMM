from __future__ import annotations

import pytest

from cmm.core.flux_state import reference_state_pfba
from cmm.omics.differential import (
    differential_expression,
    gene_directions,
    reaction_directions,
    restrict_to_top_changed,
)
from cobra import Metabolite, Model, Reaction

# Disease (source) vs healthy (target) gene expression for the branched model.
SOURCE = {"g1": 50.0, "g2": 100.0, "g3": 1.0, "g5": 1.0, "gb": 50.0}
TARGET = {"g1": 50.0, "g2": 1.0, "g3": 100.0, "g5": 100.0, "gb": 50.0}


def test_gene_directions_discretize_fold_change():
    dirs = gene_directions(SOURCE, TARGET)
    assert dirs["g3"] == 1  # strongly up in healthy
    assert dirs["g5"] == 1
    assert dirs["g2"] == -1  # strongly down in healthy
    assert dirs["g1"] == 0  # unchanged
    assert dirs["gb"] == 0


def test_reaction_directions_with_reference(branched_model):
    reference = reference_state_pfba(branched_model, name="disease")
    g_dirs = gene_directions(SOURCE, TARGET)
    dmap = reaction_directions(branched_model, g_dirs, reference=reference)

    # R2 active in source (v>0), enzyme down -> flux should decrease.
    assert dmap["R2"] == -1
    # R3/R5 inactive in source, enzyme up -> turn on (increase).
    assert dmap["R3"] == 1
    assert dmap["R5"] == 1
    # Unchanged enzymes and the gene-less supply stay steady.
    assert dmap["R1"] == 0
    assert dmap["BIOMASS"] == 0
    assert dmap["SUP_A"] == 0

    assert dmap.forward() == frozenset({"R3", "R5"})
    assert dmap.backward() == frozenset({"R2"})
    assert "R1" in dmap.steady()


def test_differential_expression_convenience(branched_model):
    reference = reference_state_pfba(branched_model, name="disease")
    dmap = differential_expression(branched_model, SOURCE, TARGET, reference=reference)
    assert dmap.nonsteady() == frozenset({"R2", "R3", "R5"})


def _two_gene_gpr_model():
    model = Model("gpr")
    x = Metabolite("x_c", compartment="c")
    r_and = Reaction("R_AND")
    r_and.add_metabolites({x: -1})
    r_and.gene_reaction_rule = "gA and gB"
    r_or = Reaction("R_OR")
    r_or.add_metabolites({x: -1})
    r_or.gene_reaction_rule = "gA or gB"
    model.add_reactions([r_and, r_or])
    return model


def _paper_rule(op: str, labels: tuple[int, ...]) -> int:
    """Yizhak et al. 2013 prose, transcribed directly: 'all' / 'at least one' / mixed -> 0."""

    has_up = any(label > 0 for label in labels)
    has_down = any(label < 0 for label in labels)
    if has_up and has_down:
        return 0  # "a subset ... elevated and another subset ... reduced" -> unchanged
    if op == "and":
        if all(label > 0 for label in labels):
            return 1
        if all(label < 0 for label in labels):
            return -1
        return 0
    if has_up:
        return 1
    if has_down:
        return -1
    return 0


@pytest.mark.parametrize("n_genes", [2, 3])
def test_gpr_direction_matches_yizhak_rule_on_every_label_combination(n_genes):
    from itertools import product

    from cmm.omics.differential import _eval_gpr_direction

    genes = [f"g{i}" for i in range(n_genes)]
    model = Model("gpr")
    x = Metabolite("x_c", compartment="c")
    for op in ("and", "or"):
        rxn = Reaction(f"R_{op.upper()}")
        rxn.add_metabolites({x: -1})
        rxn.gene_reaction_rule = f" {op} ".join(genes)
        model.add_reactions([rxn])

    for labels in product((-1, 0, 1), repeat=n_genes):
        g_dirs = dict(zip(genes, labels, strict=True))
        for op in ("and", "or"):
            rxn = model.reactions.get_by_id(f"R_{op.upper()}")
            assert _eval_gpr_direction(rxn.gpr, g_dirs) == _paper_rule(op, labels), (
                f"{op} {labels}"
            )


def test_gpr_and_or_combination():
    # gA up (+1), gB down (-1): conflicting evidence is 'unchanged' on both node types
    # (Yizhak et al. 2013). CMM's pre-0.4.0 signed min/max gave -1 and +1 here.
    model = _two_gene_gpr_model()
    dmap = reaction_directions(model, {"gA": 1, "gB": -1})
    assert dmap["R_AND"] == 0
    assert dmap["R_OR"] == 0


def test_gpr_or_propagates_a_single_reduced_isozyme():
    # OR(-1, 0) is 'at least one reduced' -> -1; signed max would have returned 0.
    model = _two_gene_gpr_model()
    dmap = reaction_directions(model, {"gA": 0, "gB": -1})
    assert dmap["R_OR"] == -1
    assert dmap["R_AND"] == 0  # not all subunits reduced
    assert dmap.metadata["gpr_rule"] == "yizhak2013_two_pass_binary"


def test_differential_expression_rejects_invalid_values():
    with pytest.raises(ValueError, match="non-negative"):
        gene_directions({"g": -1.0}, {"g": 1.0})
    with pytest.raises(ValueError, match="pseudocount"):
        gene_directions({"g": 1.0}, {"g": 2.0}, pseudocount=-1.0)


# --- Yizhak et al.'s top-N cut on the changed set -----------------------------


def _four_reaction_model():
    """Four single-gene reactions, so each gene's label lands on exactly one reaction."""

    model = Model("top_n")
    mets = {name: Metabolite(f"{name}_c", compartment="c") for name in "abcde"}
    for index, (gene, left, right) in enumerate(
        [("g1", "a", "b"), ("g2", "b", "c"), ("g3", "c", "d"), ("g4", "d", "e")]
    ):
        rxn = Reaction(f"R{index + 1}")
        rxn.bounds = (0.0, 100.0)
        rxn.add_metabolites({mets[left]: -1, mets[right]: 1})
        rxn.gene_reaction_rule = gene
        model.add_reactions([rxn])
    return model


def test_top_changed_cut_keeps_the_strongest_and_zeroes_the_rest():
    model = _four_reaction_model()
    # Strength order: R4 (8) > R3 (4) > R2 (2) > R1 (1); all four are labelled -1.
    lfc = {"g1": -1.0, "g2": -2.0, "g3": -4.0, "g4": -8.0}
    direction = reaction_directions(model, {g: -1 for g in lfc})
    assert len(direction.nonsteady()) == 4

    cut = restrict_to_top_changed(model, direction, lfc, 2)
    assert cut.nonsteady() == {"R3", "R4"}
    assert cut["R1"] == 0 and cut["R2"] == 0
    # The retained labels keep their sign; the cut only moves reactions into the steady set.
    assert cut["R3"] == -1 and cut["R4"] == -1
    assert cut.metadata["n_changed_before_cut"] == 4
    assert cut.metadata["n_changed_kept"] == 2
    assert cut.metadata["gpr_rule"] == "yizhak2013_two_pass_binary"


def test_top_changed_cut_is_a_no_op_when_the_set_is_already_small():
    model = _four_reaction_model()
    lfc = {"g1": -1.0, "g2": -2.0, "g3": -4.0, "g4": -8.0}
    direction = reaction_directions(model, {g: -1 for g in lfc})
    cut = restrict_to_top_changed(model, direction, lfc, 10)
    assert cut.nonsteady() == direction.nonsteady()
    assert cut.metadata["n_changed_kept"] == 4


def test_top_changed_ties_break_on_reaction_id():
    model = _four_reaction_model()
    lfc = {"g1": -3.0, "g2": -3.0, "g3": -3.0, "g4": -3.0}
    direction = reaction_directions(model, {g: -1 for g in lfc})
    first = restrict_to_top_changed(model, direction, lfc, 2)
    second = restrict_to_top_changed(model, direction, lfc, 2)
    assert first.nonsteady() == second.nonsteady() == {"R1", "R2"}


def test_differential_expression_applies_the_cut():
    model = _four_reaction_model()
    source = {"g1": 8.0, "g2": 8.0, "g3": 8.0, "g4": 8.0}
    target = {"g1": 3.0, "g2": 1.0, "g3": 0.0, "g4": 0.0}
    uncut = differential_expression(model, source, target)
    cut = differential_expression(model, source, target, top_n_changed=1)
    assert len(cut.nonsteady()) == 1
    assert cut.nonsteady() < uncut.nonsteady()
    assert cut.metadata["top_n_changed"] == 1


def test_top_changed_cut_rejects_a_negative_count():
    model = _four_reaction_model()
    direction = reaction_directions(model, {"g1": -1})
    with pytest.raises(ValueError, match="non-negative"):
        restrict_to_top_changed(model, direction, {"g1": -1.0}, -1)


# --- Student's t-test over replicates ----------------------------------------


def _replicates(values, columns=("r1", "r2", "r3")):
    import pandas as pd

    return pd.DataFrame(values, columns=list(columns))


def test_t_test_labels_genes_by_direction_of_the_shift():
    import numpy as np

    from cmm.omics.differential import gene_directions_from_replicates

    source = _replicates(np.array([[8.0, 8.1, 7.9], [8.0, 8.1, 7.9], [8.0, 8.1, 7.9]]))
    target = _replicates(
        np.array([[11.0, 11.1, 10.9], [5.0, 5.1, 4.9], [8.0, 8.1, 7.9]])
    )
    frame = gene_directions_from_replicates(source, target)
    assert list(frame["direction"]) == [1, -1, 0]
    # The evidence travels with the label so a report can say why a gene was called.
    assert {"log2_fold_change", "t_statistic", "p_value", "significant"} <= set(frame)
    assert frame["log2_fold_change"].iloc[0] == pytest.approx(3.0, abs=0.01)


def test_t_test_needs_replicates_and_names_the_alternative():
    import numpy as np

    from cmm.omics.differential import gene_directions_from_replicates

    single = _replicates(np.array([[8.0], [8.0]]), columns=("only",))
    with pytest.raises(ValueError, match="gene_directions\\(\\)"):
        gene_directions_from_replicates(single, single)


def test_t_test_rejects_disjoint_identifiers():
    import numpy as np
    import pandas as pd

    from cmm.omics.differential import gene_directions_from_replicates

    a = pd.DataFrame(np.ones((2, 3)), index=["a", "b"])
    b = pd.DataFrame(np.ones((2, 3)), index=["x", "y"])
    with pytest.raises(ValueError, match="share no gene ids"):
        gene_directions_from_replicates(a, b)


def test_t_test_cutoff_is_validated():
    import numpy as np

    from cmm.omics.differential import gene_directions_from_replicates

    frame = _replicates(np.ones((2, 3)))
    with pytest.raises(ValueError, match="p_value_cutoff"):
        gene_directions_from_replicates(frame, frame, p_value_cutoff=0.0)


def test_top_changed_cut_can_rank_on_p_value_when_replicates_exist():
    # Fold change alone would keep R1: it has the largest |log2FC|. The p-value ordering
    # keeps R2, whose shift is smaller but far better supported — which is the paper's
    # gene-level test carried to the reaction level.
    model = _four_reaction_model()
    lfc = {"g1": -5.0, "g2": -2.0, "g3": -0.5, "g4": -0.4}
    p_values = {"g1": 0.40, "g2": 0.001, "g3": 0.9, "g4": 0.95}
    direction = reaction_directions(model, {g: -1 for g in lfc})

    by_fold = restrict_to_top_changed(model, direction, lfc, 1)
    by_p = restrict_to_top_changed(model, direction, lfc, 1, gene_p_values=p_values)
    assert by_fold.nonsteady() == {"R1"}
    assert by_p.nonsteady() == {"R2"}
    assert "log2FC" in by_fold.metadata["changed_ranking"]
    assert "p-value" in by_p.metadata["changed_ranking"]
