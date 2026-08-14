"""Reproducible random flux sampling.

FBA returns one optimal flux distribution out of many that are equally optimal, so a single
solution cannot say whether a predicted flux is *forced* by the network or merely one
arbitrary choice among alternate optima. Sampling answers that: it draws flux distributions
from the feasible space and reports the distribution of each reaction's flux, which is the
evidence a predicted engineering target needs before it is believed.

Two services are provided:

* :func:`random_flux_sampling` samples the feasible space as constrained by the model and
  the given condition — the uniform reading.
* :func:`reference_constrained_sampling` first narrows each reaction to a window around a
  reference flux state (from FBA, pFBA, LAD, or E-Flux2) and samples inside it — the
  uncertainty-around-a-prediction reading.

Both are seeded and single-process by default so a result can be reproduced exactly, and
both return samples that satisfy the steady-state and bound constraints. Note that adding
noise to sampled fluxes afterwards would break both, so this module never does it.

**Two deviations from optGpSampler as published (Megchelenbrink W, Huynen M & Marchiori E
(2014), *PLoS ONE* 9(2):e86587), disclosed rather than changed.**

1. ``processes=1`` runs a single chain. optGpSampler's defining feature is many parallel
   chains whose warm-up points come from one shared optimization; a single chain is a
   correct but slower and less well-mixed sampler, not the published algorithm. It is the
   default because parallel chains are seeded independently, so ``processes>1`` cannot be
   reproduced bit-for-bit — and reproducibility is the property this module exists to
   provide. Raise ``processes`` when throughput matters more than an exactly repeatable run,
   and say in the report that you did.
2. ``thinning=100`` is COBRApy's default and is below the thinning the paper needed at
   genome scale, where mixing required on the order of k ≈ 500–2500 steps between retained
   samples. On a genome-scale reconstruction, raise ``thinning`` (and ``n``) and check
   convergence rather than trusting the default; on small models such as ``e_coli_core`` the
   default is ample.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Literal

import pandas as pd
from cobra import Model
from cobra.sampling import sample as _cobra_sample

from cmm.core.condition import Condition
from cmm.core.flux_state import FluxState, reference_state_from_samples
from cmm.core.provenance import run_provenance

SamplingMethod = Literal["optgp", "achr"]

_STATISTIC_COLUMNS = ("mean", "std", "minimum", "q1", "median", "q3", "maximum")


@dataclass(frozen=True)
class SamplingResult:
    """A sampled flux ensemble and the statistics derived from it."""

    samples: pd.DataFrame  # rows = samples, columns = reaction ids
    method: str
    seed: int
    metadata: dict[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        # frozen=True only blocks reassignment; copy so a caller mutating the frame cannot
        # corrupt this result (same contract as FseofResult.trends).
        object.__setattr__(self, "samples", self.samples.copy())

    @property
    def n_samples(self) -> int:
        return int(len(self.samples))

    def to_frame(self) -> pd.DataFrame:
        """The raw sample table — one row per sample, one column per reaction."""

        return self.samples.copy()

    def statistics(self) -> pd.DataFrame:
        """Per-reaction summary: mean, std, min, quartiles, max."""

        described = self.samples.describe()
        table = pd.DataFrame(
            {
                "mean": described.loc["mean"],
                "std": described.loc["std"],
                "minimum": described.loc["min"],
                "q1": described.loc["25%"],
                "median": described.loc["50%"],
                "q3": described.loc["75%"],
                "maximum": described.loc["max"],
            }
        )
        table.index.name = "reaction_id"
        return table[list(_STATISTIC_COLUMNS)]

    def correlation(self, *, min_std: float = 1e-6) -> pd.DataFrame:
        """Flux correlation between reactions that actually vary.

        Reactions whose sampled flux is constant carry no information and would produce
        undefined correlations, so they are dropped rather than reported as NaN.
        """

        if min_std < 0:
            raise ValueError("min_std must be non-negative")
        varying = self.samples.std()
        columns = varying[varying > min_std].index
        if len(columns) == 0:
            return pd.DataFrame()
        return self.samples[columns].corr()

    def to_flux_state(self, name: str = "sampled") -> FluxState:
        """Sample means as a :class:`FluxState` usable as a reference for MOMA/ROOM/MTA."""

        return reference_state_from_samples(self.samples, name=name)


def _validate_sampling_arguments(
    *, n: int, thinning: int, processes: int, method: str
) -> None:
    if n < 1:
        raise ValueError("n must be at least 1")
    if thinning < 1:
        raise ValueError("thinning must be at least 1")
    if processes < 1:
        raise ValueError("processes must be at least 1")
    if method not in ("optgp", "achr"):
        raise ValueError(f"unknown sampling method {method!r}; use 'optgp' or 'achr'")
    if method == "achr" and processes != 1:
        raise ValueError("the 'achr' sampler is single-process; use processes=1")


def _draw(
    model: Model, *, n: int, method: str, thinning: int, processes: int, seed: int
) -> pd.DataFrame:
    try:
        return _cobra_sample(
            model,
            n,
            method=method,
            thinning=thinning,
            processes=processes,
            seed=seed,
        )
    except Exception as error:  # pragma: no cover - solver/geometry dependent
        raise RuntimeError(
            f"flux sampling failed ({error}); the feasible space may be empty or "
            "degenerate under the current bounds"
        ) from error


def random_flux_sampling(
    model: Model,
    n: int = 1000,
    *,
    condition: Condition | None = None,
    method: SamplingMethod = "optgp",
    thinning: int = 100,
    processes: int = 1,
    seed: int = 0,
) -> SamplingResult:
    """Sample flux distributions uniformly from the model's feasible space.

    ``optgp`` needs a large ``n`` (roughly >1000) to mix well; ``achr`` converges better for
    small sample counts. ``processes`` defaults to 1 because parallel chains are seeded
    independently, so raising it makes a run faster but no longer bit-for-bit reproducible —
    which also means the published multi-chain optGpSampler design is switched off by
    default. ``thinning`` defaults to COBRApy's 100, below the k ≈ 500–2500 optGpSampler
    needed at genome scale; raise both on a genome-scale reconstruction. Both deviations are
    set out in the module docstring.

    The sampled space is the model *as constrained*, so it is the caller's job to apply a
    medium or condition first — sampling an unconstrained model samples an unbounded space.
    """

    _validate_sampling_arguments(
        n=n, thinning=thinning, processes=processes, method=method
    )

    provenance = run_provenance(
        model,
        method="random_flux_sampling",
        sampler=method,
        n=n,
        thinning=thinning,
        processes=processes,
        seed=seed,
        condition=condition.name if condition else None,
        reproducible=processes == 1,
    )

    with model:
        if condition is not None:
            condition.apply_to(model)
        samples = _draw(
            model, n=n, method=method, thinning=thinning, processes=processes, seed=seed
        )

    return SamplingResult(
        samples=samples, method=method, seed=seed, metadata=provenance
    )


def reference_constrained_sampling(
    model: Model,
    reference: FluxState | Mapping[str, float],
    n: int = 1000,
    *,
    condition: Condition | None = None,
    min_fraction: float = 0.8,
    max_fraction: float = 1.2,
    zero_tolerance: float = 1e-6,
    zero_window: float = 0.1,
    method: SamplingMethod = "optgp",
    thinning: int = 100,
    processes: int = 1,
    seed: int = 0,
) -> SamplingResult:
    """Sample the flux space narrowed to a window around a reference flux state.

    Each reaction carrying reference flux ``v`` is restricted to ``[v * min_fraction,
    v * max_fraction]`` (mirrored for negative ``v``) intersected with its existing bounds;
    reactions at essentially zero flux are restricted to ``[-zero_window, zero_window]``.
    The result quantifies how much a *predicted* state could vary while staying close to the
    prediction — uncertainty around a solution, not the whole feasible space.

    Because ``min_fraction <= 1 <= max_fraction``, every window contains its reference flux,
    so an empty window can only mean the reference violates the model's own bounds — which
    raises rather than silently sampling a different space.

    ``processes`` and ``thinning`` carry the same defaults, and the same two disclosed
    deviations from published optGpSampler, as :func:`random_flux_sampling`.
    """

    _validate_sampling_arguments(
        n=n, thinning=thinning, processes=processes, method=method
    )
    if not 0.0 <= min_fraction <= 1.0:
        raise ValueError("min_fraction must be in [0, 1]")
    if max_fraction < 1.0:
        raise ValueError("max_fraction must be at least 1")
    if zero_tolerance < 0 or zero_window <= 0:
        raise ValueError("zero_tolerance must be >= 0 and zero_window must be > 0")

    fluxes = reference.fluxes if isinstance(reference, FluxState) else dict(reference)
    if not fluxes:
        raise ValueError("reference flux state is empty")

    provenance = run_provenance(
        model,
        method="reference_constrained_sampling",
        sampler=method,
        n=n,
        thinning=thinning,
        processes=processes,
        seed=seed,
        condition=condition.name if condition else None,
        min_fraction=min_fraction,
        max_fraction=max_fraction,
        zero_tolerance=zero_tolerance,
        zero_window=zero_window,
        reference=getattr(reference, "name", "mapping"),
        reproducible=processes == 1,
    )

    applied: dict[str, tuple[float, float]] = {}
    conflicts: list[str] = []

    with model:
        if condition is not None:
            condition.apply_to(model)

        for reaction_id, flux in fluxes.items():
            if reaction_id not in model.reactions:
                continue
            reaction = model.reactions.get_by_id(reaction_id)
            value = float(flux)

            if abs(value) <= zero_tolerance:
                window = (-zero_window, zero_window)
            elif value > 0:
                window = (value * min_fraction, value * max_fraction)
            else:
                window = (value * max_fraction, value * min_fraction)

            lower = max(float(reaction.lower_bound), window[0])
            upper = min(float(reaction.upper_bound), window[1])
            if lower > upper:
                conflicts.append(reaction_id)
                continue

            reaction.bounds = (lower, upper)
            applied[reaction_id] = (lower, upper)

        if conflicts:
            shown = ", ".join(sorted(conflicts)[:5])
            more = "" if len(conflicts) <= 5 else f" (and {len(conflicts) - 5} more)"
            raise ValueError(
                f"reference flux lies outside the model bounds for {len(conflicts)} "
                f"reaction(s): {shown}{more}. The reference was probably produced under a "
                "different condition than the one being sampled"
            )

        samples = _draw(
            model, n=n, method=method, thinning=thinning, processes=processes, seed=seed
        )

    provenance["parameters"]["constrained_reactions"] = len(applied)  # type: ignore[index]
    return SamplingResult(
        samples=samples, method=method, seed=seed, metadata=provenance
    )
