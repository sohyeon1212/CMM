"""MTA/rMTA's source reference state v_ref, built as the papers build it.

Yizhak et al. (2013, Methods) "used iMAT to obtain a solution space consistent with the
source's expression, and then sampled 2,000 flux distributions from the space consistent with
the activity/inactivity pattern of one of the optimal iMAT solutions. v_ref is then defined as
the mean of these samples." Valcárcel et al. (2019) follow the same protocol with
optGpSampler (2,000 points, 500 steps between points).

:func:`imat_reference_state` does exactly that: it solves iMAT
(:func:`cmm.omics.imat.imat`), confines every reaction iMAT matched to its solved pattern
(:meth:`cmm.omics.imat.IMATResult.pattern_bounds`), samples that space with COBRApy's OptGP
sampler — an implementation of the same algorithm as optGpSampler (Megchelenbrink et al.
2014) — and averages. The same samples then give the per-reaction significance thresholds
(:func:`cmm.features.revert_epsilon.epsilons_from_samples`), as both papers derive them.

**Differences from the published runs, disclosed rather than hidden.**

* COBRApy builds its warm-up points itself (two per reaction, from flux-range extremes) and
  does not expose a count; optGpSampler was run with 5,000.
* Sampling is one seeded chain (``processes=1``) so the result is reproducible; optGpSampler
  runs several chains in parallel.

Using one of possibly several iMAT optima is the papers' own procedure, not a difference; the
pattern the solver returned is kept on the result so it can be inspected.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import pandas as pd
from cobra import Model

from cmm.core.flux_state import FluxState
from cmm.core.provenance import run_provenance
from cmm.features.sampling import random_flux_sampling
from cmm.omics.imat import PUBLISHED_IMAT_EPSILON, IMATResult, imat

#: Yizhak et al. (2013) and Valcárcel et al. (2019) both average 2,000 samples.
PUBLISHED_REFERENCE_SAMPLES = 2000
#: Valcárcel et al. used 500 steps between points.
PUBLISHED_REFERENCE_THINNING = 500


@dataclass(frozen=True)
class IMATReferenceState:
    """v_ref as the mean of samples from one iMAT solution's activity pattern."""

    imat: IMATResult
    fluxes: dict[str, float]
    bounds: dict[str, tuple[float, float]]
    n_samples: int
    metadata: dict[str, object] = field(default_factory=dict)
    #: The full sample table (rows = samples). Kept for the epsilon derivation; at genome
    #: scale it is large, so callers export only the columns they need.
    samples: pd.DataFrame | None = field(default=None, compare=False, repr=False)

    def to_flux_state(self, name: str = "source_state_imat") -> FluxState:
        return FluxState(
            fluxes=self.fluxes,
            name=name,
            provenance="sampling_mean",
            metadata=self.metadata,
        )

    def to_frame(self) -> pd.DataFrame:
        """One row per reaction: iMAT class and activity, imposed bounds, and v_ref."""

        frame = self.imat.to_frame()
        frame["pattern_lower_bound"] = frame["reaction_id"].map(
            lambda rid: self.bounds.get(rid, (None, None))[0]
        )
        frame["pattern_upper_bound"] = frame["reaction_id"].map(
            lambda rid: self.bounds.get(rid, (None, None))[1]
        )
        frame["v_ref"] = frame["reaction_id"].map(self.fluxes)
        return frame


class IMATReferenceError(RuntimeError):
    """iMAT found no optimal solution, so there is no pattern to sample."""


def imat_reference_state(
    model: Model,
    gene_states: Mapping[str, int],
    *,
    epsilon: float = PUBLISHED_IMAT_EPSILON,
    n_samples: int = PUBLISHED_REFERENCE_SAMPLES,
    thinning: int = PUBLISHED_REFERENCE_THINNING,
    seed: int = 0,
) -> IMATReferenceState:
    """Solve iMAT, sample its activity pattern, and return the sample mean as v_ref.

    The model is expected to carry the run's medium and condition already; it is left
    unchanged. ``gene_states`` comes from :func:`cmm.omics.imat.discretize_expression`.
    ``n_samples`` and ``thinning`` are part of the method, not precision knobs: v_ref and every
    sampled epsilon depend on them.
    """

    if n_samples < 2:
        raise ValueError("n_samples must be at least 2")
    solution = imat(model, gene_states, epsilon=epsilon)
    if solution.status != "optimal":
        raise IMATReferenceError(
            f"iMAT found no optimal solution ({solution.status}); there is no activity "
            "pattern to sample"
        )
    bounds = solution.pattern_bounds(model)
    provenance = run_provenance(
        model,
        seed=seed,
        method="imat_reference_state",
        imat_epsilon=epsilon,
        n_samples=n_samples,
        thinning=thinning,
        sampler="optgp",
        processes=1,
        sampled_space="imat_activity_pattern",
        n_pattern_bounds=len(bounds),
        imat=solution.summary(),
        imat_provenance=solution.metadata,
    )
    with model:
        for rid, (lb, ub) in bounds.items():
            model.reactions.get_by_id(rid).bounds = (lb, ub)
        sampled = random_flux_sampling(
            model, n_samples, method="optgp", thinning=thinning, processes=1, seed=seed
        )
    samples = sampled.samples
    means = samples.mean(axis=0)
    return IMATReferenceState(
        imat=solution,
        fluxes={str(rid): float(value) for rid, value in means.items()},
        bounds=bounds,
        n_samples=int(len(samples)),
        metadata=provenance,
        samples=samples,
    )
