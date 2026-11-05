"""iMAT (Shlomi et al. 2008) with the rMTA discretisation, and MTA's sampled v_ref."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from cobra.util.array import create_stoichiometric_matrix

from cmm.features import IMATReferenceError, imat_reference_state
from cmm.features.revert_reference import (
    PUBLISHED_REFERENCE_SAMPLES,
    PUBLISHED_REFERENCE_THINNING,
)
from cmm.omics import discretize_expression, imat
from cmm.omics.imat import IMAT_GPR_RULE, PUBLISHED_IMAT_BETA, PUBLISHED_IMAT_EPSILON

# g1, g2 and gb high; g3 and g5 low: the "disease" branch R2 is on, R3 -> R5 is off.
STATES = {"g1": 1, "g2": 1, "g3": -1, "g5": -1, "gb": 1}


def test_published_constants():
    assert PUBLISHED_IMAT_EPSILON == 1.0
    assert PUBLISHED_IMAT_BETA == 1.0
    assert PUBLISHED_REFERENCE_SAMPLES == 2000
    assert PUBLISHED_REFERENCE_THINNING == 500


# --- discretisation -----------------------------------------------------------


def test_discretisation_cuts_each_sample_at_mean_plus_minus_delta_sd():
    values = pd.DataFrame({"s1": [0.0, 2.0, 4.0, 6.0, 8.0]}, index=list("abcde"))
    mean, sd = 4.0, float(np.std([0, 2, 4, 6, 8], ddof=1))  # sd = 3.162...
    states = discretize_expression(values, delta=0.5)
    # high > 4 + 1.58 -> d, e; low < 4 - 1.58 -> a, b; c is moderate.
    assert states.states == {"a": -1, "b": -1, "c": 0, "d": 1, "e": 1}
    assert mean + 0.5 * sd < 6.0
    # A wider band leaves only the extremes classed.
    wide = discretize_expression(values, delta=1.0)
    assert wide.states == {"a": -1, "b": 0, "c": 0, "d": 0, "e": 1}


def test_beta_is_the_share_of_samples_that_must_agree():
    values = pd.DataFrame(
        {
            "s1": [0.0, 5.0, 10.0, 10.0],
            "s2": [0.0, 5.0, 10.0, 10.0],
            "s3": [0.0, 5.0, 0.0, 10.0],  # c drops to low in one sample
        },
        index=list("abcd"),
    )
    strict = discretize_expression(values, delta=0.1)
    assert strict.states["c"] == 0  # beta = 100%: two of three is not enough
    assert strict.states["d"] == 1
    assert strict.high_fraction["c"] == pytest.approx(2 / 3)
    loose = discretize_expression(values, delta=0.1, beta=0.6)
    assert loose.states["c"] == 1
    summary = strict.summary()
    assert summary["n_samples"] == 3 and summary["beta"] == 1.0


def test_thresholds_are_set_on_the_genes_kept():
    # A non-model gene with a huge value moves the mean; restricting to model genes undoes it.
    values = pd.DataFrame({"s": [1.0, 2.0, 3.0, 100.0]}, index=["a", "b", "c", "x"])
    everything = discretize_expression(values, delta=0.3)
    model_only = discretize_expression(values, delta=0.3, genes=["a", "b", "c", "z"])
    assert "x" not in model_only.states and "z" not in model_only.states
    assert everything.states["c"] == -1
    assert model_only.states == {"a": -1, "b": 0, "c": 1}
    frame = model_only.to_frame()
    assert list(frame.columns) == ["gene_id", "state", "high_fraction", "low_fraction"]


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"delta": -0.1}, "delta"),
        ({"delta": 0.3, "beta": 0.0}, "beta"),
        ({"delta": 0.3, "beta": 1.1}, "beta"),
        ({"delta": 0.3, "genes": ["a"]}, "two genes"),
    ],
)
def test_discretisation_rejects_bad_settings(kwargs, message):
    values = pd.DataFrame({"s": [1.0, 2.0]}, index=["a", "b"])
    with pytest.raises(ValueError, match=message):
        discretize_expression(values, **kwargs)


# --- the iMAT MILP --------------------------------------------------------------


def test_imat_matches_every_label_it_can_and_leaves_the_model_alone(branched_model):
    bounds = {r.id: r.bounds for r in branched_model.reactions}
    objective = str(branched_model.objective.expression)

    result = imat(branched_model, STATES)

    assert result.status == "optimal"
    # R_H = R1, R2, BIOMASS and R_L = R3, R5 are all satisfiable at once.
    assert result.objective_value == pytest.approx(5.0)
    assert result.n_matched == 5
    assert result.activity == {
        "R1": "active_forward",
        "R2": "active_forward",
        "BIOMASS": "active_forward",
        "R3": "inactive",
        "R5": "inactive",
    }
    assert result.reaction_states["SUP_A"] == 0  # no GPR: moderate
    assert all(abs(result.fluxes[r]) <= 1e-9 for r in ("R3", "R5"))
    assert result.fluxes["R2"] >= PUBLISHED_IMAT_EPSILON - 1e-9
    assert {r.id: r.bounds for r in branched_model.reactions} == bounds
    assert str(branched_model.objective.expression) == objective
    parameters = result.metadata["parameters"]
    assert parameters["threshold"] == 0.0  # the paper's v = 0, not imatpy's relaxation
    assert parameters["gpr_rule"] == IMAT_GPR_RULE
    assert parameters["implementation"] == "imatpy"


def test_imat_reports_what_it_cannot_match(branched_model):
    # R2 low but BIOMASS high: with R3/R5 also low, P has no producer, so one label must go.
    states = {"g1": 1, "g2": -1, "g3": -1, "g5": -1, "gb": 1}
    result = imat(branched_model, states)
    assert result.status == "optimal"
    assert result.objective_value == pytest.approx(4.0)
    assert list(result.activity.values()).count("not_matched") == 1
    assert result.summary()["n_not_matched"] == 1


def test_gpr_uses_min_for_and_max_for_or(branched_model):
    branched_model.reactions.R2.gene_reaction_rule = "g2 and g3"  # complex
    branched_model.reactions.R5.gene_reaction_rule = "g5 or g2"  # isozymes
    result = imat(branched_model, {"g2": 1, "g3": -1, "g5": -1})
    assert (
        result.reaction_states["R2"] == -1
    )  # a complex is as scarce as its scarcest part
    assert result.reaction_states["R5"] == 1  # isozymes as their most available member


def test_unlisted_genes_are_moderate_and_bad_states_are_refused(branched_model):
    result = imat(branched_model, {"g2": 1})
    assert result.reaction_states == {
        "SUP_A": 0,
        "R1": 0,
        "R2": 1,
        "R3": 0,
        "R5": 0,
        "BIOMASS": 0,
    }
    with pytest.raises(ValueError, match="-1, 0 or 1"):
        imat(branched_model, {"g2": 2})
    with pytest.raises(ValueError, match="positive"):
        imat(branched_model, STATES, epsilon=0.0)


def test_pattern_bounds_keep_the_solved_direction(branched_model):
    branched_model.reactions.R2.lower_bound = -1000.0  # reversible
    result = imat(branched_model, STATES)
    bounds = result.pattern_bounds(branched_model)
    lb, ub = bounds["R2"]
    if result.activity["R2"] == "active_forward":
        assert (lb, ub) == (1.0, 1000.0)
    else:
        assert (lb, ub) == (-1000.0, -1.0)
    assert bounds["R3"] == (0.0, 0.0)
    assert "SUP_A" not in bounds  # moderate reactions keep the model's bounds
    frame = result.to_frame()
    assert set(frame["activity"]) >= {"active_forward", "inactive", "moderate"}


# --- v_ref from the pattern's space -------------------------------------------------


def test_reference_state_is_the_mean_of_samples_inside_the_pattern(branched_model):
    bounds = {r.id: r.bounds for r in branched_model.reactions}
    first = imat_reference_state(branched_model, STATES, n_samples=200, thinning=10)
    second = imat_reference_state(branched_model, STATES, n_samples=200, thinning=10)

    assert first.fluxes == second.fluxes  # seeded single chain
    assert {r.id: r.bounds for r in branched_model.reactions} == bounds
    samples = first.samples
    assert samples is not None and len(samples) == 200 == first.n_samples
    # Every sample honours the pattern: R3/R5 off, the high reactions at least epsilon.
    assert (samples[["R3", "R5"]].abs() <= 1e-9).all().all()
    assert (samples[["R1", "R2", "BIOMASS"]] >= 1.0 - 1e-6).all().all()
    means = samples.mean(axis=0)
    assert first.fluxes == pytest.approx({k: float(v) for k, v in means.items()})
    # A mean of steady states is a steady state.
    order = [r.id for r in branched_model.reactions]
    stoichiometry = create_stoichiometric_matrix(branched_model)
    residual = stoichiometry @ np.array([first.fluxes[r] for r in order])
    assert np.abs(residual).max() < 1e-6

    state = first.to_flux_state()
    assert state.provenance == "sampling_mean"
    parameters = first.metadata["parameters"]
    assert parameters["sampled_space"] == "imat_activity_pattern"
    assert parameters["imat"]["objective_value"] == pytest.approx(5.0)
    assert first.metadata["seed"] == 0
    frame = first.to_frame()
    assert list(frame.columns) == [
        "reaction_id",
        "expression_state",
        "activity",
        "imat_flux",
        "pattern_lower_bound",
        "pattern_upper_bound",
        "v_ref",
    ]
    r3 = frame.set_index("reaction_id").loc["R3"]
    assert (r3["pattern_lower_bound"], r3["pattern_upper_bound"]) == (0.0, 0.0)


def test_no_imat_solution_means_no_reference_state(branched_model):
    branched_model.reactions.SUP_A.lower_bound = 5.0
    branched_model.reactions.BIOMASS.upper_bound = 0.0  # nothing can drain the supply
    result = imat(branched_model, STATES)
    assert result.status != "optimal" and result.activity == {}
    with pytest.raises(IMATReferenceError, match="no optimal solution"):
        imat_reference_state(branched_model, STATES, n_samples=10, thinning=1)
    with pytest.raises(ValueError, match="n_samples"):
        imat_reference_state(branched_model, STATES, n_samples=1)
