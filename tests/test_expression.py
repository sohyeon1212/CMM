from __future__ import annotations

import pytest
from cmm.core.solvers import SolverCapabilityError
from cmm.omics.expression import (
    eflux2,
    gene_to_reaction_weights,
    integrate_expression,
    lad,
)
from cobra import Metabolite, Model, Reaction

BIOMASS = "Biomass_Ecoli_core"


def _gpr_model():
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


def test_gene_to_reaction_weights_gpr_combination():
    model = _gpr_model()
    weights = gene_to_reaction_weights(model, {"gA": 5.0, "gB": 2.0})
    assert weights["R_AND"] == 2.0  # complex -> min subunit
    assert weights["R_OR"] == 5.0  # isozymes -> max


def test_gene_to_reaction_weights_case_insensitive():
    model = _gpr_model()
    weights = gene_to_reaction_weights(model, {"GA": 5.0, "GB": 2.0})
    assert weights["R_OR"] == 5.0


def test_gene_to_reaction_weights_missing_complex_subunit_drops_reaction():
    model = _gpr_model()
    weights = gene_to_reaction_weights(model, {"gA": 5.0})  # gB missing
    assert "R_AND" not in weights  # AND with a missing subunit -> None
    assert weights["R_OR"] == 5.0  # OR keeps the available isozyme


def test_eflux2_runs_and_holds_objective(ecoli_core):
    weights = {r.id: 50.0 for r in ecoli_core.reactions if r.genes}
    result = eflux2(ecoli_core, weights)
    assert result.status == "optimal"
    assert result.detail == "qp"  # gurobi supports the L2 step
    assert result.fluxes[BIOMASS] == pytest.approx(result.objective_value, abs=1e-6)
    assert result.objective_value == pytest.approx(
        result.metadata["unconstrained_optimum"], abs=1e-6
    )


def test_eflux2_does_not_mutate_model(ecoli_core):
    growth = ecoli_core.slim_optimize()
    weights = {r.id: 50.0 for r in ecoli_core.reactions if r.genes}
    eflux2(ecoli_core, weights)
    assert ecoli_core.slim_optimize() == pytest.approx(growth, abs=1e-6)


def test_eflux2_requires_qp_unless_approximation_is_explicit(ecoli_core):
    ecoli_core.solver = "glpk"  # LP+MILP only, no QP
    weights = {r.id: 50.0 for r in ecoli_core.reactions if r.genes}
    with pytest.raises(SolverCapabilityError):
        eflux2(ecoli_core, weights)
    result = eflux2(ecoli_core, weights, allow_l1_fallback=True)
    assert result.status == "optimal"
    assert result.method == "eflux2_l1_fallback"
    assert result.detail == "pfba_l1_approximation"


def test_lad_runs_on_any_solver(ecoli_core):
    ecoli_core.solver = "glpk"
    weights = {r.id: 5.0 for r in ecoli_core.reactions if r.genes}
    result = lad(ecoli_core, weights)
    assert result.status == "optimal"
    assert result.objective_value is not None
    # Model restored afterwards.
    assert ecoli_core.slim_optimize() == pytest.approx(0.8739, abs=1e-3)


def test_eflux2_min_direction_floor_branch(ecoli_core):
    # Exercise the min-objective floor branch (ub instead of lb): it must not raise.
    ecoli_core.objective_direction = "min"
    weights = {r.id: 50.0 for r in ecoli_core.reactions if r.genes}
    result = eflux2(ecoli_core, weights)
    assert result.status in ("optimal", "infeasible")


def test_eflux2_empty_expression_returns_no_targets(ecoli_core):
    assert eflux2(ecoli_core, {}).status == "no_targets"
    assert lad(ecoli_core, {}).status == "no_targets"


def test_integrate_expression_dispatch(ecoli_core):
    expr = {g.id: 20.0 for g in ecoli_core.genes}
    assert integrate_expression(ecoli_core, expr, method="eflux2").status == "optimal"
    assert integrate_expression(ecoli_core, expr, method="lad").status == "optimal"
    with pytest.raises(ValueError, match="unknown omics method"):
        integrate_expression(ecoli_core, expr, method="bogus")


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"objective_fraction": 0.0}, "objective_fraction"),
        ({"objective_fraction": 1.1}, "objective_fraction"),
        ({"min_scale": -0.1}, "min_scale"),
        ({"normalization_percentile": 101}, "normalization_percentile"),
    ],
)
def test_eflux2_validates_parameters(ecoli_core, kwargs, message):
    weights = {r.id: 1.0 for r in ecoli_core.reactions if r.genes}
    with pytest.raises(ValueError, match=message):
        eflux2(ecoli_core, weights, **kwargs)


def test_omics_rejects_negative_or_nonfinite_expression(ecoli_core):
    rid = next(r.id for r in ecoli_core.reactions if r.genes)
    with pytest.raises(ValueError, match="non-negative"):
        eflux2(ecoli_core, {rid: -1.0})
    gene = next(iter(ecoli_core.genes)).id
    with pytest.raises(ValueError, match="finite"):
        gene_to_reaction_weights(ecoli_core, {gene: float("nan")})
