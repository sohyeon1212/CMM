"""Computational strain design: OptKnock and RobustKnock.

OptKnock (Burgard 2003) searches for a small reaction-knockout set that couples a target
product to growth — at maximum growth the cell is *forced* to make the product. RobustKnock
(Tepper & Shlomi 2010) additionally requires the coupling to hold under the worst case (the
guaranteed minimum product at maximum growth, not the optimistic maximum).

The bilevel MILP search is delegated to the ``straindesign`` package; each returned design is
then evaluated for its optimistic and guaranteed product at maximum growth so OptKnock and
RobustKnock can be ranked on the right quantity.
"""

from __future__ import annotations

import contextlib
import io
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

from cobra import Model

from cmm.core import solvers


@dataclass(frozen=True)
class StrainDesign:
    """One knockout design and the product it forces at maximum growth."""

    knockouts: tuple[str, ...]
    growth: float
    max_product: float        # optimistic: max product at max growth
    guaranteed_product: float  # robust: min product at max growth

    @property
    def growth_coupled(self) -> bool:
        return self.guaranteed_product > 1e-6


@dataclass(frozen=True)
class StrainDesignResult:
    method: str
    product: str
    designs: tuple[StrainDesign, ...]

    def best(self) -> StrainDesign | None:
        return self.designs[0] if self.designs else None


def _objective_reaction(model: Model) -> str:
    for rxn in model.reactions:
        if rxn.objective_coefficient != 0:
            return rxn.id
    raise ValueError("model has no objective reaction")


def _evaluate_design(
    model: Model, knockouts: Iterable[str], product: str, biomass: str
) -> StrainDesign | None:
    """Growth and optimistic/guaranteed product for a knockout set, at maximum growth."""

    kos = tuple(knockouts)
    with model:
        for rid in kos:
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        growth = model.slim_optimize()
        if growth is None or growth != growth or growth < 1e-6:  # NaN/infeasible/no growth
            return None
        biomass_rxn = model.reactions.get_by_id(biomass)
        biomass_rxn.bounds = (0.999 * growth, growth)
        product_rxn = model.reactions.get_by_id(product)
        model.objective = product_rxn
        model.objective_direction = "max"
        max_product = model.slim_optimize()
        model.objective_direction = "min"
        guaranteed = model.slim_optimize()
    return StrainDesign(
        knockouts=kos,
        growth=float(growth),
        max_product=float(max_product) if max_product is not None else 0.0,
        guaranteed_product=float(guaranteed) if guaranteed is not None else 0.0,
    )


def _search_designs(
    model: Model,
    product: str,
    biomass: str,
    *,
    max_knockouts: int,
    max_solutions: int,
    min_growth: float,
) -> list[tuple[str, ...]]:
    try:
        import straindesign as sd
    except ImportError as exc:  # pragma: no cover - straindesign is a declared dependency
        raise RuntimeError("strain design requires the 'straindesign' package") from exc

    # A product that cannot carry flux in the current medium has no coupling design — and
    # straindesign prunes such reactions and then crashes on the missing outer objective.
    with model:
        model.objective = model.reactions.get_by_id(product)
        model.objective_direction = "max"
        if (model.slim_optimize() or 0.0) < 1e-6:
            return []

    module = sd.SDModule(
        model, sd.OPTKNOCK,
        inner_objective=biomass,
        outer_objective=product,
        constraints=[f"{biomass} >= {min_growth}"],
    )
    # straindesign is chatty; silence its solver logs.
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        solutions = sd.compute_strain_designs(
            model, sd_modules=[module], max_cost=max_knockouts, max_solutions=max_solutions
        )
    reaction_sd = getattr(solutions, "reaction_sd", None) or []
    designs: list[tuple[str, ...]] = []
    for design in reaction_sd:
        kos = tuple(sorted(rid for rid, coeff in design.items() if coeff <= 0))
        if kos and kos not in designs:
            designs.append(kos)
    return designs


def _strain_design(
    model: Model,
    product: str,
    *,
    method: Literal["optknock", "robustknock"],
    biomass: str | None,
    max_knockouts: int,
    max_solutions: int,
    min_growth: float,
) -> StrainDesignResult:
    solvers.require("MILP", model.solver.interface, feature=method)
    biomass = biomass or _objective_reaction(model)
    ko_sets = _search_designs(
        model, product, biomass,
        max_knockouts=max_knockouts, max_solutions=max_solutions, min_growth=min_growth,
    )
    evaluated = [
        d for d in (_evaluate_design(model, kos, product, biomass) for kos in ko_sets)
        if d is not None
    ]
    if method == "robustknock":
        evaluated = [d for d in evaluated if d.growth_coupled]
        evaluated.sort(key=lambda d: (-d.guaranteed_product, -d.growth))
    else:
        evaluated.sort(key=lambda d: (-d.max_product, -d.growth))
    return StrainDesignResult(method=method, product=product, designs=tuple(evaluated))


def optknock(
    model: Model,
    product: str,
    *,
    biomass: str | None = None,
    max_knockouts: int = 3,
    max_solutions: int = 5,
    min_growth: float = 0.05,
) -> StrainDesignResult:
    """OptKnock: knockout sets that maximize the product at maximum growth (optimistic)."""

    return _strain_design(
        model, product, method="optknock", biomass=biomass,
        max_knockouts=max_knockouts, max_solutions=max_solutions, min_growth=min_growth,
    )


def robustknock(
    model: Model,
    product: str,
    *,
    biomass: str | None = None,
    max_knockouts: int = 3,
    max_solutions: int = 8,
    min_growth: float = 0.05,
) -> StrainDesignResult:
    """RobustKnock: keep only designs that guarantee product at maximum growth (worst case)."""

    return _strain_design(
        model, product, method="robustknock", biomass=biomass,
        max_knockouts=max_knockouts, max_solutions=max_solutions, min_growth=min_growth,
    )
