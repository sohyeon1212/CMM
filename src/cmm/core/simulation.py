"""Flux simulation services."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from typing import Literal

from cobra import Model
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis import pfba as _cobra_pfba
import pandas as pd

from cmm.core.condition import Condition
from cmm.core.provenance import run_provenance

#: Loop-removal algorithms cobra accepts. ``"fastSNP"`` returns optimal loopless bounds;
#: ``"cycleFreeFlux"`` is cheaper but does not guarantee them.
LooplessAlgorithm = Literal["fastSNP", "cycleFreeFlux"]

# Below this many analysed reactions, FVA runs in-process: a pool costs seconds of pure
# overhead on a small model and, on macOS spawn semantics, fails outright unless the caller
# guards the call with ``if __name__ == "__main__":``.
_PROCESS_POOL_MIN_REACTIONS = 500


def _loopless_algorithm(loopless: bool | LooplessAlgorithm) -> LooplessAlgorithm | None:
    """Resolve the ``loopless`` argument to the algorithm name cobra expects.

    ``True`` maps to ``"cycleFreeFlux"``, which is what cobra's deprecated boolean does, so
    the two spellings stay equivalent while CMM avoids the deprecation path.
    """

    if isinstance(loopless, bool):
        return "cycleFreeFlux" if loopless else None
    if loopless not in ("fastSNP", "cycleFreeFlux"):
        raise ValueError(
            f"loopless must be a bool, 'fastSNP' or 'cycleFreeFlux'; got {loopless!r}"
        )
    return loopless


@dataclass(frozen=True)
class FluxSolution:
    """Serializable result of an optimization run."""

    status: str
    objective_value: float | None
    fluxes: dict[str, float]
    metadata: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class FluxRange:
    """Minimum and maximum feasible flux for a reaction."""

    minimum: float
    maximum: float


@dataclass(frozen=True)
class FvaResult(Mapping[str, FluxRange]):
    """Flux ranges with the run provenance that produced them.

    Behaves as the ``{reaction_id: FluxRange}`` mapping :func:`fva` used to return, so
    ``result[reaction_id]`` and ``result.get(reaction_id)`` are unchanged, while
    ``metadata`` carries the same :func:`~cmm.core.provenance.run_provenance` record that
    :class:`FluxSolution` carries for FBA and pFBA.
    """

    ranges: Mapping[str, FluxRange]
    metadata: dict[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "ranges", dict(self.ranges))

    def __getitem__(self, reaction_id: str) -> FluxRange:
        return self.ranges[reaction_id]

    def __iter__(self) -> Iterator[str]:
        return iter(self.ranges)

    def __len__(self) -> int:
        return len(self.ranges)

    def to_frame(self) -> pd.DataFrame:
        """Deterministic export table (reaction_id, minimum, maximum)."""

        records = [
            {
                "reaction_id": reaction_id,
                "minimum": flux_range.minimum,
                "maximum": flux_range.maximum,
            }
            for reaction_id, flux_range in self.ranges.items()
        ]
        return pd.DataFrame(records, columns=["reaction_id", "minimum", "maximum"])


def fba(model: Model, condition: Condition | None = None) -> FluxSolution:
    """Run flux balance analysis and return a plain Python result."""

    provenance = run_provenance(
        model, method="fba", condition=condition.name if condition else None
    )
    with model:
        if condition is not None:
            condition.apply_to(model)
        solution = model.optimize()

        if solution.status != "optimal":
            return FluxSolution(
                status=solution.status,
                objective_value=None,
                fluxes={},
                metadata=provenance,
            )
        objective_value = solution.objective_value
        if objective_value is not None:
            objective_value = float(objective_value)

        return FluxSolution(
            status=solution.status,
            objective_value=objective_value,
            fluxes={
                reaction_id: float(value)
                for reaction_id, value in solution.fluxes.items()
            },
            metadata=provenance,
        )


def pfba(
    model: Model,
    condition: Condition | None = None,
    fraction_of_optimum: float = 1.0,
) -> FluxSolution:
    """Run parsimonious FBA (minimal total flux at the given fraction of the optimum)."""

    if not 0.0 < fraction_of_optimum <= 1.0:
        raise ValueError("fraction_of_optimum must be in (0, 1]")
    provenance = run_provenance(
        model,
        method="pfba",
        condition=condition.name if condition else None,
        fraction_of_optimum=fraction_of_optimum,
    )
    with model:
        if condition is not None:
            condition.apply_to(model)
        solution = _cobra_pfba(model, fraction_of_optimum=fraction_of_optimum)
        objective_value = solution.objective_value
        if objective_value is not None:
            objective_value = float(objective_value)
        return FluxSolution(
            status=solution.status,
            objective_value=objective_value,
            fluxes={
                reaction_id: float(value)
                for reaction_id, value in solution.fluxes.items()
            },
            metadata=provenance,
        )


def fva(
    model: Model,
    condition: Condition | None = None,
    reactions: Iterable[str] | None = None,
    fraction_of_optimum: float = 1.0,
    *,
    loopless: bool | LooplessAlgorithm = False,
    processes: int | None = None,
) -> FvaResult:
    """Run flux variability analysis for selected reactions.

    ``loopless`` is forwarded to cobra, which removes thermodynamically infeasible loops
    from the reported ranges at a substantially higher cost per reaction. ``True`` selects
    ``"cycleFreeFlux"``, the algorithm cobra's own (now deprecated) boolean maps to; pass
    ``"fastSNP"`` for optimal loopless bounds instead. The algorithm actually used is
    recorded in provenance as ``loopless_algorithm``.

    ``processes`` is forwarded too. It defaults to ``None``, meaning *decide from the
    problem size*: fewer than 500 analysed reactions run in-process, because a pool costs
    seconds of pure overhead on a small model and fails on macOS unless the caller wraps the
    call in ``if __name__ == "__main__":``. Pass an integer to force the choice.
    """

    if fraction_of_optimum < 0 or fraction_of_optimum > 1:
        raise ValueError("fraction_of_optimum must be between 0 and 1")

    algorithm = _loopless_algorithm(loopless)
    requested = list(reactions) if reactions is not None else None
    n_reactions = len(requested) if requested is not None else len(model.reactions)
    if processes is None and n_reactions < _PROCESS_POOL_MIN_REACTIONS:
        processes = 1
    provenance = run_provenance(
        model,
        method="fva",
        condition=condition.name if condition else None,
        fraction_of_optimum=fraction_of_optimum,
        loopless=loopless,
        loopless_algorithm=algorithm,
        processes=processes,
        n_reactions=n_reactions,
    )

    with model:
        if condition is not None:
            condition.apply_to(model)

        reaction_list = None
        if requested is not None:
            reaction_list = [
                model.reactions.get_by_id(reaction_id) for reaction_id in requested
            ]

        table = flux_variability_analysis(
            model,
            reaction_list=reaction_list,
            fraction_of_optimum=fraction_of_optimum,
            loopless=algorithm,
            processes=processes,
        )

        return FvaResult(
            ranges={
                str(reaction_id): FluxRange(
                    minimum=float(row["minimum"]),
                    maximum=float(row["maximum"]),
                )
                for reaction_id, row in table.iterrows()
            },
            metadata=provenance,
        )
