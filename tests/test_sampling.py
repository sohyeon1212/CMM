from __future__ import annotations

import pytest
from cmm.core.condition import Condition, ReactionBound
from cmm.core.flux_state import FluxState
from cmm.features.comparison import reference_flux
from cmm.features.sampling import (
    random_flux_sampling,
    reference_constrained_sampling,
)

O2 = "EX_o2_e"
BIOMASS = "Biomass_Ecoli_core"

# 'achr' converges better than 'optgp' at the small sample counts a test can afford.
SMALL = {"n": 120, "method": "achr", "thinning": 10}


def test_samples_are_a_reaction_table(ecoli_core):
    result = random_flux_sampling(ecoli_core, seed=3, **SMALL)
    assert result.n_samples == 120
    assert result.method == "achr"
    assert result.seed == 3
    frame = result.to_frame()
    assert frame.shape[1] == len(ecoli_core.reactions)
    assert BIOMASS in frame.columns

    tidy = result.long_frame()
    assert list(tidy.columns) == ["sample_id", "reaction_id", "flux"]
    assert len(tidy) == result.n_samples * len(ecoli_core.reactions)
    first_sample = tidy[tidy["sample_id"] == result.samples.index[0]].set_index(
        "reaction_id"
    )["flux"]
    assert first_sample.to_dict() == result.samples.iloc[0].to_dict()


def test_same_seed_reproduces_the_same_samples(ecoli_core):
    first = random_flux_sampling(ecoli_core, seed=11, **SMALL)
    second = random_flux_sampling(ecoli_core, seed=11, **SMALL)
    assert first.samples.equals(second.samples)


def test_different_seeds_give_different_samples(ecoli_core):
    first = random_flux_sampling(ecoli_core, seed=1, **SMALL)
    second = random_flux_sampling(ecoli_core, seed=2, **SMALL)
    assert not first.samples.equals(second.samples)


def test_samples_respect_the_model_bounds(ecoli_core):
    result = random_flux_sampling(ecoli_core, seed=5, **SMALL)
    for reaction in ecoli_core.reactions:
        column = result.samples[reaction.id]
        assert column.min() >= reaction.lower_bound - 1e-6
        assert column.max() <= reaction.upper_bound + 1e-6


def test_samples_satisfy_the_steady_state(ecoli_core):
    # S . v = 0 for every drawn sample; this is what makes them flux distributions and
    # what post-hoc noise would destroy.
    import numpy as np
    from cobra.util.array import create_stoichiometric_matrix

    result = random_flux_sampling(ecoli_core, seed=5, **SMALL)
    matrix = np.ascontiguousarray(
        create_stoichiometric_matrix(ecoli_core), dtype=np.float64
    )
    ordered = np.ascontiguousarray(
        result.samples[[r.id for r in ecoli_core.reactions]].to_numpy(),
        dtype=np.float64,
    )
    residual = np.abs(matrix.dot(ordered.T))
    assert np.isfinite(residual).all()
    assert residual.max() < 1e-6


def test_statistics_table_shape(ecoli_core):
    result = random_flux_sampling(ecoli_core, seed=5, **SMALL)
    table = result.statistics()
    assert list(table.columns) == [
        "mean",
        "std",
        "minimum",
        "q1",
        "median",
        "q3",
        "maximum",
    ]
    assert table.index.name == "reaction_id"
    assert len(table) == len(ecoli_core.reactions)
    row = table.loc[BIOMASS]
    assert row["minimum"] <= row["median"] <= row["maximum"]


def test_correlation_drops_constant_reactions(ecoli_core):
    result = random_flux_sampling(ecoli_core, seed=5, **SMALL)
    correlation = result.correlation()
    constant = result.samples.std()
    dropped = set(constant[constant <= 1e-6].index)
    assert dropped.isdisjoint(set(correlation.columns))
    assert correlation.notna().all().all()


def test_flux_state_bridge(ecoli_core):
    result = random_flux_sampling(ecoli_core, seed=5, **SMALL)
    state = result.to_flux_state(name="sampled_wt")
    assert isinstance(state, FluxState)
    assert state.name == "sampled_wt"
    assert state.provenance == "sampling_mean"
    assert state.metadata["n_samples"] == 120
    assert state.get(BIOMASS) == pytest.approx(result.samples[BIOMASS].mean())


def test_result_frame_is_a_copy(ecoli_core):
    result = random_flux_sampling(ecoli_core, seed=5, **SMALL)
    frame = result.to_frame()
    frame.iloc[0, 0] = 12345.0
    assert result.samples.iloc[0, 0] != 12345.0


def test_condition_narrows_the_sampled_space(ecoli_core):
    anaerobic = Condition(
        name="anaerobic",
        bounds=(ReactionBound(reaction_id=O2, lower_bound=0.0),),
    )
    before = ecoli_core.reactions.get_by_id(O2).lower_bound
    result = random_flux_sampling(ecoli_core, condition=anaerobic, seed=5, **SMALL)
    assert result.samples[O2].min() >= -1e-6
    assert result.metadata["parameters"]["condition"] == "anaerobic"
    assert ecoli_core.reactions.get_by_id(O2).lower_bound == before


def test_sampling_does_not_mutate_the_model(ecoli_core):
    bounds_before = {r.id: r.bounds for r in ecoli_core.reactions}
    random_flux_sampling(ecoli_core, seed=5, **SMALL)
    assert {r.id: r.bounds for r in ecoli_core.reactions} == bounds_before


def test_provenance_records_reproducibility(ecoli_core):
    result = random_flux_sampling(ecoli_core, seed=5, processes=1, **SMALL)
    parameters = result.metadata["parameters"]
    assert parameters["method"] == "random_flux_sampling"
    assert parameters["sampler"] == "achr"
    assert parameters["seed"] == 5
    assert parameters["reproducible"] is True


def test_rejects_invalid_sampling_arguments(ecoli_core):
    with pytest.raises(ValueError, match="n must be"):
        random_flux_sampling(ecoli_core, n=0)
    with pytest.raises(ValueError, match="thinning"):
        random_flux_sampling(ecoli_core, thinning=0)
    with pytest.raises(ValueError, match="processes"):
        random_flux_sampling(ecoli_core, processes=0)
    with pytest.raises(ValueError, match="unknown sampling method"):
        random_flux_sampling(ecoli_core, method="gibbs")
    with pytest.raises(ValueError, match="single-process"):
        random_flux_sampling(ecoli_core, method="achr", processes=4)


def test_reference_constrained_sampling_stays_near_the_reference(ecoli_core):
    reference = reference_flux(ecoli_core, "pfba")
    result = reference_constrained_sampling(
        ecoli_core,
        reference,
        min_fraction=0.9,
        max_fraction=1.1,
        seed=5,
        **SMALL,
    )
    # Every reaction carrying real reference flux stays inside its +/-10% window.
    for reaction_id, flux in reference.fluxes.items():
        if abs(flux) < 1e-6:
            continue
        column = result.samples[reaction_id]
        low, high = sorted((flux * 0.9, flux * 1.1))
        assert column.min() >= low - 1e-6
        assert column.max() <= high + 1e-6


def test_reference_constrained_sampling_accepts_a_mapping(ecoli_core):
    reference = reference_flux(ecoli_core, "pfba")
    result = reference_constrained_sampling(
        ecoli_core, dict(reference.fluxes), seed=5, **SMALL
    )
    assert result.n_samples == 120
    assert result.metadata["parameters"]["reference"] == "mapping"
    assert result.metadata["parameters"]["constrained_reactions"] > 0


def test_reference_constrained_sampling_is_narrower_than_uniform(ecoli_core):
    reference = reference_flux(ecoli_core, "pfba")
    uniform = random_flux_sampling(ecoli_core, seed=5, **SMALL)
    constrained = reference_constrained_sampling(
        ecoli_core, reference, min_fraction=0.9, max_fraction=1.1, seed=5, **SMALL
    )
    assert constrained.samples[BIOMASS].std() < uniform.samples[BIOMASS].std()


def test_reference_outside_the_model_bounds_raises(ecoli_core):
    reference = reference_flux(ecoli_core, "pfba")
    # A reference produced under a different condition can exceed the sampled bounds.
    inconsistent = FluxState(
        fluxes={**reference.fluxes, "PGI": 5000.0},
        name="other_condition",
        provenance="pfba",
    )
    with pytest.raises(ValueError, match="outside the model bounds"):
        reference_constrained_sampling(ecoli_core, inconsistent, seed=5, **SMALL)


def test_reference_constrained_rejects_bad_fractions(ecoli_core):
    reference = reference_flux(ecoli_core, "pfba")
    with pytest.raises(ValueError, match="min_fraction"):
        reference_constrained_sampling(ecoli_core, reference, min_fraction=1.5)
    with pytest.raises(ValueError, match="max_fraction"):
        reference_constrained_sampling(ecoli_core, reference, max_fraction=0.5)
    with pytest.raises(ValueError, match="zero_tolerance"):
        reference_constrained_sampling(ecoli_core, reference, zero_window=0.0)


def test_reference_constrained_rejects_an_empty_reference(ecoli_core):
    with pytest.raises(ValueError, match="empty"):
        reference_constrained_sampling(ecoli_core, {}, seed=5, **SMALL)
