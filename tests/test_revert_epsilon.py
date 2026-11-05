from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from cmm.features import ReactionEpsilons, epsilons_from_samples
from cmm.features.revert_epsilon import (
    PUBLISHED_EPSILON_MINIMUM,
    PUBLISHED_EPSILON_SAMPLES,
)


def test_published_constants():
    assert PUBLISHED_EPSILON_MINIMUM == 1e-3
    assert PUBLISHED_EPSILON_SAMPLES == 2000


def test_epsilon_is_the_one_sided_t_test_threshold_with_a_floor():
    rng = np.random.default_rng(0)
    n = 2000
    samples = pd.DataFrame(
        {
            "wide": rng.normal(5.0, 2.0, n),
            "narrow": rng.normal(1.0, 1e-4, n),
            "fixed": np.full(n, 3.0),
        }
    )
    result = epsilons_from_samples(samples, ["wide", "narrow", "fixed", "absent"])

    s = samples["wide"].std(ddof=1)
    expected = stats.t.ppf(0.95, n - 1) * s / math.sqrt(n)
    assert result.epsilons["wide"] == pytest.approx(expected)
    # At this n the threshold is about 0.037 standard deviations.
    assert result.epsilons["wide"] / s == pytest.approx(0.0368, abs=5e-4)
    # The minimum applies to a reaction that barely varies, never varies, or was not sampled.
    for rid in ("narrow", "fixed", "absent"):
        assert result.epsilons[rid] == PUBLISHED_EPSILON_MINIMUM
    assert result.n_floored == 3
    assert result.sample_std["absent"] == 0.0


def test_the_threshold_closes_the_one_sample_t_test_at_the_cutoff():
    # Shifting the reference by epsilon must put the sampled values exactly at P = 0.05.
    rng = np.random.default_rng(1)
    values = rng.normal(0.0, 3.0, 2000)
    eps = epsilons_from_samples(pd.DataFrame({"r": values}), ["r"]).epsilons["r"]
    p = stats.ttest_1samp(values, values.mean() + eps, alternative="less").pvalue
    assert p == pytest.approx(0.05, rel=1e-6)


def test_frame_summary_and_order():
    samples = pd.DataFrame({"b": [0.0, 10.0, 20.0], "a": [1.0, 1.0, 1.0]})
    result = epsilons_from_samples(samples, ["b", "a", "b"])
    frame = result.to_frame()
    assert list(frame["reaction_id"]) == ["a", "b"]
    assert list(frame.columns) == [
        "reaction_id",
        "sample_std",
        "t_test_epsilon",
        "epsilon",
        "floored",
    ]
    assert frame.set_index("reaction_id").loc["a", "floored"]
    summary = result.summary()
    assert summary["n_reactions"] == 2
    assert summary["n_floored"] == 1
    assert summary["minimum"] == PUBLISHED_EPSILON_MINIMUM
    assert summary["maximum"] == result.epsilons["b"]
    assert summary["median"] == pytest.approx(
        (result.epsilons["a"] + result.epsilons["b"]) / 2
    )
    empty = ReactionEpsilons({}, {}, 2, 0.05, 1e-3)
    assert empty.summary() == {"n_reactions": 0, "n_floored": 0}


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [({"p_value": 0.0}, "p_value"), ({"minimum": -1.0}, "minimum")],
)
def test_epsilons_from_samples_validates(kwargs, message):
    samples = pd.DataFrame({"r": [0.0, 1.0]})
    with pytest.raises(ValueError, match=message):
        epsilons_from_samples(samples, ["r"], **kwargs)
    with pytest.raises(ValueError, match="two samples"):
        epsilons_from_samples(pd.DataFrame({"r": [0.0]}), ["r"])


def test_only_the_requested_reactions_samples_are_kept():
    samples = pd.DataFrame(
        {"a": [0.0, 1.0, 2.0], "b": [1.0, 1.0, 1.0], "c": [5.0, 6.0, 7.0]}
    )
    result = epsilons_from_samples(samples, ["c", "a", "missing"])
    kept = result.samples_frame()
    assert list(kept.columns) == ["c", "a"]  # requested order; unsampled ids absent
    assert kept.index.name == "sample_id" and len(kept) == 3
    assert ReactionEpsilons({}, {}, 2, 0.05, 1e-3).samples_frame().empty
