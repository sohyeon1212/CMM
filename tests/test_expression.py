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
    assert weights["R_AND"] == 2.0  # complex -> min subunit (Kim 2016, Lee 2012)
    assert weights["R_OR"] == 7.0  # isozymes -> sum (Kim 2016, Lee 2012)


def test_gene_to_reaction_weights_or_rule_max_reproduces_pre_040():
    model = _gpr_model()
    weights = gene_to_reaction_weights(model, {"gA": 5.0, "gB": 2.0}, or_rule="max")
    assert weights["R_AND"] == 2.0
    assert weights["R_OR"] == 5.0
    with pytest.raises(ValueError, match="or_rule"):
        gene_to_reaction_weights(model, {"gA": 5.0}, or_rule="mean")


def test_gene_to_reaction_weights_case_insensitive():
    model = _gpr_model()
    weights = gene_to_reaction_weights(model, {"GA": 5.0, "GB": 2.0})
    assert weights["R_OR"] == 7.0


def test_gene_to_reaction_weights_missing_complex_subunit_drops_reaction():
    model = _gpr_model()
    weights = gene_to_reaction_weights(model, {"gA": 5.0})  # gB missing
    assert "R_AND" not in weights  # AND with a missing subunit -> None
    assert weights["R_OR"] == 5.0  # OR sums the available isozymes


@pytest.mark.requires_qp
def test_eflux2_runs_and_holds_objective(ecoli_core):
    weights = {r.id: 50.0 for r in ecoli_core.reactions if r.genes}
    result = eflux2(ecoli_core, weights)
    assert result.status == "optimal"
    assert result.detail == "qp"  # gurobi supports the L2 step
    assert result.fluxes[BIOMASS] == pytest.approx(result.objective_value, abs=1e-6)
    assert result.objective_value == pytest.approx(
        result.metadata["unconstrained_optimum"], abs=1e-6
    )


@pytest.mark.requires_qp
def test_eflux2_bounds_follow_each_reactions_own_expression(ecoli_core):
    # Non-uniform weights, so every scaled bound differs from its template by its own ratio.
    import numpy as np

    rng = np.random.default_rng(1)
    gene_reactions = [r.id for r in ecoli_core.reactions if r.genes]
    weights = {rid: float(rng.uniform(1.0, 100.0)) for rid in gene_reactions}
    weights["PFK"] = 0.01  # floored at min_scale: ub 1000 -> 1
    denom = max(weights.values())

    result = eflux2(ecoli_core, weights)
    assert result.status == "optimal"
    n_changed = 0
    for rid, (lb, ub) in result.bounds.items():
        template = ecoli_core.reactions.get_by_id(rid)
        ratio = min(1.0, max(1e-3, weights[rid] / denom))
        expected_lb = (
            template.lower_bound * ratio
            if template.lower_bound < 0
            else template.lower_bound
        )
        expected_ub = (
            template.upper_bound * ratio
            if template.upper_bound > 0
            else template.upper_bound
        )
        assert (lb, ub) == pytest.approx((expected_lb, expected_ub)), rid
        n_changed += (lb, ub) != template.bounds
        assert lb - 1e-6 <= result.fluxes[rid] <= ub + 1e-6
    assert n_changed > len(result.bounds) // 2
    # Exchanges and the biomass reaction are never scaled.
    assert not any(
        ecoli_core.reactions.get_by_id(rid).boundary for rid in result.bounds
    )
    assert BIOMASS not in result.bounds

    # A barely expressed PFK caps glycolytic flux that uniform expression lets through.
    assert result.bounds["PFK"][1] == pytest.approx(1.0)
    uniform = eflux2(ecoli_core, {rid: 50.0 for rid in gene_reactions})
    assert uniform.fluxes["PFK"] > 1.0 + 1e-3
    assert result.fluxes["PFK"] <= 1.0 + 1e-6


@pytest.mark.requires_qp
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
    # The archived provenance must agree with the reported method label.
    assert result.metadata["parameters"]["method"] == "eflux2_l1_fallback"


def _reverse_only_model():
    """Model whose reaction ``R`` can only carry negative flux (bounds -10 .. -5)."""

    model = Model("reverse")
    x = Metabolite("x_c", compartment="c")
    r = Reaction("R", lower_bound=-10, upper_bound=-5)
    r.add_metabolites({x: -1})
    sink = Reaction("SINK", lower_bound=0, upper_bound=1000)
    sink.add_metabolites({x: -1})
    model.add_reactions([r, sink])
    return model


def test_lad_fits_absolute_flux_for_a_reverse_reaction():
    # Lee et al. 2012 fit |v|: a reaction pinned at v = -5 with target 5 deviates by 0.
    # A signed residual (CMM before 0.4.0) would score |-5 - 5| = 10.
    model = _reverse_only_model()
    result = lad(model, {"R": 5.0})
    assert result.status == "optimal"
    assert result.objective_value == pytest.approx(0.0, abs=1e-6)
    assert result.fluxes["R"] == pytest.approx(-5.0, abs=1e-6)


def test_lad_sigma_weights_scale_the_objective(ecoli_core):
    weights = {r.id: 5.0 for r in ecoli_core.reactions if r.genes}
    unweighted = lad(ecoli_core, weights)
    weighted = lad(ecoli_core, weights, reaction_sigma={r: 2.0 for r in weights})
    assert weighted.status == "optimal"
    assert weighted.objective_value == pytest.approx(
        unweighted.objective_value / 2.0, rel=1e-6
    )
    assert weighted.metadata["sigma_weighted"] is True
    assert unweighted.metadata["sigma_weighted"] is False
    with pytest.raises(ValueError, match="reaction_sigma"):
        lad(ecoli_core, weights, reaction_sigma={next(iter(weights)): 0.0})


def test_lad_default_threshold_keeps_zero_expression_targets(ecoli_core):
    # Lee et al. drive zero-expression reactions toward zero flux rather than dropping them.
    rid = next(r.id for r in ecoli_core.reactions if r.genes)
    result = lad(ecoli_core, {rid: 0.0})
    assert result.status == "optimal"
    assert result.metadata["parameters"]["weight_threshold"] == 0.0


def test_integrate_expression_records_the_gpr_or_rule(ecoli_core):
    expr = {g.id: 20.0 for g in ecoli_core.genes}
    result = integrate_expression(ecoli_core, expr, method="lad")
    assert result.metadata["parameters"]["gpr_or_rule"] == "sum"
    legacy = integrate_expression(ecoli_core, expr, method="lad", or_rule="max")
    assert legacy.metadata["parameters"]["gpr_or_rule"] == "max"


def test_lad_runs_on_any_solver(ecoli_core):
    ecoli_core.solver = "glpk"
    weights = {r.id: 5.0 for r in ecoli_core.reactions if r.genes}
    result = lad(ecoli_core, weights)
    assert result.status == "optimal"
    assert result.objective_value is not None
    # Model restored afterwards.
    assert ecoli_core.slim_optimize() == pytest.approx(0.8739, abs=1e-3)


@pytest.mark.requires_qp
def test_eflux2_min_direction_floor_branch(ecoli_core):
    # Exercise the min-objective floor branch (ub instead of lb): it must not raise.
    ecoli_core.objective_direction = "min"
    weights = {r.id: 50.0 for r in ecoli_core.reactions if r.genes}
    result = eflux2(ecoli_core, weights)
    assert result.status in ("optimal", "infeasible")


def test_eflux2_empty_expression_returns_no_targets(ecoli_core):
    assert eflux2(ecoli_core, {}).status == "no_targets"
    assert lad(ecoli_core, {}).status == "no_targets"


@pytest.mark.requires_qp
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


@pytest.mark.requires_qp
def test_eflux2_reports_the_bounds_it_imposed(ecoli_core):
    weights = {r.id: 50.0 for r in ecoli_core.reactions if r.genes}
    weights[next(iter(weights))] = 0.0  # one reaction squeezed to the noise floor
    result = eflux2(ecoli_core, weights)
    assert result.bounds
    for rid, (lb, ub) in result.bounds.items():
        template = ecoli_core.reactions.get_by_id(rid)
        assert template.lower_bound <= lb <= max(0.0, template.lower_bound) or lb == 0
        assert min(0.0, template.upper_bound) <= ub <= template.upper_bound
        assert result.fluxes[rid] >= lb - 1e-6 and result.fluxes[rid] <= ub + 1e-6
    frame = result.bounds_frame()
    assert list(frame.columns) == ["reaction_id", "lower_bound", "upper_bound"]
    assert len(frame) == len(result.bounds)
    # LAD changes no bounds and says so by leaving the field empty.
    assert lad(ecoli_core, weights).bounds == {}


def test_expression_scale_helpers_never_log_twice():
    import numpy as np
    import pandas as pd

    from cmm.omics import linear_expression, log2_expression

    linear = pd.Series([0.0, 1.0, 3.0])
    log2 = np.log2(linear + 1.0)
    assert log2_expression(linear, "linear").tolist() == log2.tolist()
    assert log2_expression(log2, "log2") is log2
    assert linear_expression(pd.Series([0.0, 3.0]), "log2").tolist() == [1.0, 8.0]
    assert linear_expression(linear, "linear") is linear
    for helper in (linear_expression, log2_expression):
        with pytest.raises(ValueError, match="expression scale"):
            helper(linear, "ln")


def test_expression_columns_report_dropped_and_identifier_like_columns():
    import pandas as pd

    from cmm.omics import describe_expression_columns

    frame = pd.DataFrame(
        {
            "symbol": ["pgi", "pfkA", "eno"],
            "entrez": [948535, 948412, 945032],
            "s1": [10.0, 3.0, 3.0],
            "s2": [12, 4, 4],  # integer counts with ties are measurements
        },
        index=["b4025", "b3916", "b2779"],
    )
    columns = describe_expression_columns(frame)
    assert columns.samples == ("entrez", "s1", "s2")
    assert columns.dropped == ("symbol",)
    assert columns.identifier_like == ("entrez",)
    notes = columns.warnings()
    assert any("symbol" in note for note in notes)
    assert any("entrez" in note and "identifiers" in note for note in notes)

    clean = describe_expression_columns(frame[["s1", "s2"]])
    assert clean.warnings() == []
