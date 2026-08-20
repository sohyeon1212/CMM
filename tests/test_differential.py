from __future__ import annotations

import pytest

from cmm.core.flux_state import reference_state_pfba
from cmm.omics.differential import (
    differential_expression,
    gene_directions,
    reaction_directions,
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
