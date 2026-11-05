"""Per-reaction MTA/rMTA significance thresholds from a sampled source state.

MTA counts a changed reaction as transformed when its flux moves past ``v_ref ± epsilon``.
Yizhak et al. (2013, Methods, "Selection of thresholds in the MTA procedure") choose that
epsilon from the sampled source state: for each changed reaction, epsilon is the change at
which ``v_ref,i ± epsilon`` is significantly different (one-sided t-test, P < 0.05) from the
2,000 flux values sampled for that reaction. Valcárcel et al. (2019) use this per-reaction rule
for human and mouse data with a minimum of 0.001 "to avoid numerical issues".

Read as a one-sample one-sided t-test of the ``n`` sampled values against ``v_ref,i ± epsilon``
— where ``v_ref,i`` is the sample mean — the threshold is closed-form::

    epsilon_i = max(t_{1-p, n-1} * s_i / sqrt(n), minimum)

with ``s_i`` the sample standard deviation of reaction ``i``. The papers do not name the
one-sample form explicitly; it is the reading under which their description defines a single
threshold per reaction, and it is stated here so the choice can be checked.

The samples are those :func:`cmm.features.revert_reference.imat_reference_state` draws from
the space consistent with one optimal iMAT activity pattern, the same space the papers
sample, so v_ref and the thresholds come from one ensemble.

Because ``epsilon_i`` scales with ``1 / sqrt(n)``, the sample count is part of the method, not a
precision knob; the default is the papers' 2,000.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
import math

import pandas as pd

#: Valcárcel et al. (2019): "we set a minimum epsilon_i of 0.001 to avoid numerical issues".
PUBLISHED_EPSILON_MINIMUM = 1e-3
#: Yizhak et al. (2013) and Valcárcel et al. (2019) both sample 2,000 flux distributions.
PUBLISHED_EPSILON_SAMPLES = 2000


@dataclass(frozen=True)
class ReactionEpsilons:
    """Per-reaction significance thresholds and the sampled spread they came from."""

    epsilons: dict[str, float]
    sample_std: dict[str, float]
    n_samples: int
    p_value: float
    minimum: float
    metadata: dict[str, object] = field(default_factory=dict)
    #: The sampled fluxes of the requested reactions only (rows = samples), kept so the
    #: thresholds can be checked against their distribution. The full-model sample table is
    #: not kept: at genome scale it runs to hundreds of megabytes.
    samples: pd.DataFrame | None = field(default=None, compare=False, repr=False)

    def samples_frame(self) -> pd.DataFrame:
        """The kept samples, one row per sample (``sample_id``) and one column per reaction."""

        if self.samples is None:
            return pd.DataFrame()
        frame = self.samples.copy()
        frame.index.name = "sample_id"
        return frame

    @property
    def n_floored(self) -> int:
        """Reactions whose threshold is the ``minimum`` rather than the t-test value."""

        return sum(1 for rid in self.epsilons if self._raw(rid) < self.minimum)

    def _raw(self, rid: str) -> float:
        return _raw_epsilon(self.sample_std[rid], self.n_samples, self.p_value)

    def summary(self) -> dict[str, float | int]:
        values = sorted(self.epsilons.values())
        if not values:
            return {"n_reactions": 0, "n_floored": 0}
        middle = len(values) // 2
        median = (
            values[middle]
            if len(values) % 2
            else (values[middle - 1] + values[middle]) / 2.0
        )
        return {
            "n_reactions": len(values),
            "n_floored": self.n_floored,
            "minimum": values[0],
            "median": median,
            "maximum": values[-1],
        }

    def to_frame(self) -> pd.DataFrame:
        rows = [
            {
                "reaction_id": rid,
                "sample_std": self.sample_std[rid],
                "t_test_epsilon": self._raw(rid),
                "epsilon": value,
                "floored": self._raw(rid) < self.minimum,
            }
            for rid, value in sorted(self.epsilons.items())
        ]
        return pd.DataFrame(
            rows,
            columns=[
                "reaction_id",
                "sample_std",
                "t_test_epsilon",
                "epsilon",
                "floored",
            ],
        )


def _raw_epsilon(std: float, n: int, p_value: float) -> float:
    from scipy import stats

    if not math.isfinite(std) or std <= 0.0:
        return 0.0
    return float(stats.t.ppf(1.0 - p_value, n - 1) * std / math.sqrt(n))


def epsilons_from_samples(
    samples: pd.DataFrame,
    reactions: Iterable[str],
    *,
    p_value: float = 0.05,
    minimum: float = PUBLISHED_EPSILON_MINIMUM,
    metadata: Mapping[str, object] | None = None,
) -> ReactionEpsilons:
    """Apply the per-reaction threshold rule to an existing sample table.

    ``samples`` has one row per sampled flux distribution and one column per reaction. A
    requested reaction missing from the table, or one whose sampled flux does not vary, gets
    ``minimum``.
    """

    if not 0.0 < p_value < 1.0:
        raise ValueError("p_value must be in (0, 1)")
    if minimum < 0.0:
        raise ValueError("minimum must be non-negative")
    n = int(len(samples))
    if n < 2:
        raise ValueError("at least two samples are needed to estimate a spread")
    std = samples.std(ddof=1)
    wanted = list(dict.fromkeys(reactions))
    sample_std: dict[str, float] = {}
    epsilons: dict[str, float] = {}
    for rid in wanted:
        s = float(std[rid]) if rid in std.index else 0.0
        if not math.isfinite(s):
            s = 0.0
        sample_std[rid] = s
        epsilons[rid] = max(_raw_epsilon(s, n, p_value), minimum)
    return ReactionEpsilons(
        epsilons=epsilons,
        sample_std=sample_std,
        n_samples=n,
        p_value=p_value,
        minimum=minimum,
        metadata=dict(metadata or {}),
        samples=samples[[rid for rid in wanted if rid in samples.columns]].reset_index(
            drop=True
        ),
    )
