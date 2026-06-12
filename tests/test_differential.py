from __future__ import annotations

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


def test_gpr_and_or_combination():
    model = Model("gpr")
    x = Metabolite("x_c", compartment="c")
    r_and = Reaction("R_AND")
    r_and.add_metabolites({x: -1})
    r_and.gene_reaction_rule = "gA and gB"
    r_or = Reaction("R_OR")
    r_or.add_metabolites({x: -1})
    r_or.gene_reaction_rule = "gA or gB"
    model.add_reactions([r_and, r_or])

    # gA up (+1), gB down (-1): AND takes the limiting subunit (min = -1),
    # OR takes the best isozyme (max = +1).
    g_dirs = {"gA": 1, "gB": -1}
    dmap = reaction_directions(model, g_dirs)
    assert dmap["R_AND"] == -1
    assert dmap["R_OR"] == 1
