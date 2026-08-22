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
from dataclasses import field
from importlib.metadata import version as package_version
from numbers import Integral
from typing import Literal, cast

from cobra import Model
import pandas as pd

from cmm.core import solvers
from cmm.core.condition import Condition
from cmm.core.provenance import run_provenance
from cmm.features.production import (
    _actionable_reaction,
    _apply_condition,
    _condition_provenance,
    _record_condition,
)


# Gurobi's Seed parameter accepts integers through 2,000,000,000.  Using the same
# conservative range for the solver-neutral API also remains valid for the other
# straindesign backends while preventing a backend error after expensive preprocessing.
STRAIN_DESIGN_SEED_MAX = 2_000_000_000


@dataclass(frozen=True)
class StrainDesign:
    """One reaction-knockout design and the product it forces at maximum growth.

    The reaction ids are a computational design, not a gene recipe.  Experimental use
    requires resolving every reaction's GPR, including isoenzymes and multisubunit enzymes.
    """

    knockouts: tuple[str, ...]
    growth: float
    max_product: float  # optimistic: max product at max growth
    guaranteed_product: float  # robust: min product at max growth

    @property
    def growth_coupled(self) -> bool:
        return self.guaranteed_product > 1e-6


@dataclass(frozen=True)
class StrainDesignResult:
    method: str
    product: str
    designs: tuple[StrainDesign, ...]
    metadata: dict[str, object] = field(default_factory=dict)

    def best(self) -> StrainDesign | None:
        return self.designs[0] if self.designs else None

    def to_frame(self) -> pd.DataFrame:
        """One row per design, in the method's documented rank order.

        ``knockouts`` is stored as a semicolon-delimited stable string for CSV portability;
        ``n_knockouts`` makes filtering independent of parsing that display field.
        """

        return pd.DataFrame(
            [
                {
                    "rank": rank,
                    "method": self.method,
                    "product": self.product,
                    "knockouts": ";".join(design.knockouts),
                    "n_knockouts": len(design.knockouts),
                    "growth": design.growth,
                    "max_product": design.max_product,
                    "guaranteed_product": design.guaranteed_product,
                    "growth_coupled": design.growth_coupled,
                }
                for rank, design in enumerate(self.designs, start=1)
            ],
            columns=[
                "rank",
                "method",
                "product",
                "knockouts",
                "n_knockouts",
                "growth",
                "max_product",
                "guaranteed_product",
                "growth_coupled",
            ],
        )


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
        if (
            growth is None or growth != growth or growth < 1e-6
        ):  # NaN/infeasible/no growth
            return None
        biomass_rxn = model.reactions.get_by_id(biomass)
        # Evaluate the product range on the *growth-optimal face*: exactly Tepper & Shlomi's
        # Omega', the set of growth-maximal flux states, over which RobustKnock minimises.
        # Relaxing the pin to a 99.9% face does not overstate the guarantee - measured, it
        # *lowers* the reported minimum every time, because a slightly sub-optimal face
        # contains extra flux states with less product. Pinning by equality is therefore the
        # faithful reading and the less conservative-looking one; both are reasons to keep it.
        biomass_rxn.bounds = (growth, growth)
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
    method: Literal["optknock", "robustknock"],
    max_knockouts: int,
    max_solutions: int,
    min_growth: float,
    actionable_only: bool,
    seed: int,
) -> tuple[list[tuple[str, ...]], int, int, str]:
    try:
        import straindesign as sd
    except (
        ImportError
    ) as exc:  # pragma: no cover - straindesign is a declared dependency
        raise RuntimeError("strain design requires the 'straindesign' package") from exc

    # A product that cannot carry flux in the current medium has no coupling design — and
    # straindesign prunes such reactions and then crashes on the missing outer objective.
    with model:
        model.objective = model.reactions.get_by_id(product)
        model.objective_direction = "max"
        if (model.slim_optimize() or 0.0) < 1e-6:
            return [], 0, 0, "not_run_product_unreachable"

    module_type = sd.OPTKNOCK if method == "optknock" else sd.ROBUSTKNOCK
    module = sd.SDModule(
        model,
        module_type,
        inner_objective=biomass,
        outer_objective=product,
        constraints=[f"{biomass} >= {min_growth}"],
    )
    search_kwargs: dict[str, object] = {}
    n_candidates = len(model.reactions)
    if actionable_only:
        # Burgard et al. restrict the candidate set to central metabolism. Without any
        # restriction the search returns designs that delete boundary reactions with no GPR
        # (EX_co2_e, EX_ac_e, EX_for_e, EX_etoh_e, EX_lac__D_e) — not realisable as gene
        # deletions. straindesign's ko_cost restricts candidates to the listed reactions.
        candidates = {
            rxn.id: 1.0
            for rxn in model.reactions
            if _actionable_reaction(model, rxn.id, product, biomass)
        }
        n_candidates = len(candidates)
        search_kwargs["ko_cost"] = candidates
    # straindesign is chatty; silence its solver logs.
    with (
        contextlib.redirect_stdout(io.StringIO()),
        contextlib.redirect_stderr(io.StringIO()),
    ):
        solutions = sd.compute_strain_designs(
            model,
            sd_modules=[module],
            max_cost=max_knockouts,
            max_solutions=max_solutions,
            compress=True,
            seed=seed,
            **search_kwargs,
        )
    reaction_sd = getattr(solutions, "reaction_sd", None) or []
    search_status = str(getattr(solutions, "status", "unknown"))
    n_decompressed_designs = len(reaction_sd)
    designs: list[tuple[str, ...]] = []
    for design in reaction_sd:
        kos = tuple(sorted(rid for rid, coeff in design.items() if coeff <= 0))
        if kos and kos not in designs:
            designs.append(kos)
    return designs, n_candidates, n_decompressed_designs, search_status


def _strain_design(
    model: Model,
    product: str,
    *,
    method: Literal["optknock", "robustknock"],
    biomass: str | None,
    max_knockouts: int,
    max_solutions: int,
    min_growth: float,
    condition: Condition | None,
    actionable_only: bool,
    seed: int,
) -> StrainDesignResult:
    if max_knockouts < 1:
        raise ValueError("max_knockouts must be at least 1")
    if max_solutions < 1:
        raise ValueError("max_solutions must be at least 1")
    if min_growth < 0:
        raise ValueError("min_growth must be non-negative")
    if (
        isinstance(seed, bool)
        or not isinstance(seed, Integral)
        or not 0 <= seed <= STRAIN_DESIGN_SEED_MAX
    ):
        raise ValueError(
            "seed must be an integer in the Gurobi-compatible range "
            f"[0, {STRAIN_DESIGN_SEED_MAX}]"
        )
    seed = int(seed)
    solvers.require("MILP", model.solver.interface, feature=method)
    biomass = biomass or _objective_reaction(model)
    provenance = run_provenance(
        model,
        seed=seed,
        method=method,
        product=product,
        biomass=biomass,
        max_knockouts=max_knockouts,
        max_solutions=max_solutions,
        min_growth=min_growth,
        condition=condition.name if condition is not None else None,
        actionable_only=actionable_only,
        strain_design_seed=seed,
    )
    with model:
        _apply_condition(model, condition=condition)
        _record_condition(provenance, _condition_provenance(model, condition))
        ko_sets, n_candidates, n_decompressed_designs, search_status = _search_designs(
            model,
            product,
            biomass,
            method=method,
            max_knockouts=max_knockouts,
            max_solutions=max_solutions,
            min_growth=min_growth,
            actionable_only=actionable_only,
            seed=seed,
        )
        evaluated = [
            d
            for d in (_evaluate_design(model, kos, product, biomass) for kos in ko_sets)
            if d is not None
        ]
    n_evaluated_designs = len(evaluated)
    if method == "robustknock":
        evaluated = [d for d in evaluated if d.growth_coupled]
    # A high optimistic maximum is not a coupling guarantee.  Both methods use the same
    # publication rank: worst-case product first, then optimistic product and growth, with
    # the reaction set as a deterministic final tie-break.
    evaluated.sort(
        key=lambda d: (
            -d.guaranteed_product,
            -d.max_product,
            -d.growth,
            d.knockouts,
        )
    )
    deduplicated = _deduplicate_designs(evaluated)
    provenance["parameters"] = {
        **cast(dict[str, object], provenance["parameters"]),
        "n_knockout_candidates": n_candidates,
        "strain_design_backend": "straindesign",
        "straindesign_version": package_version("straindesign"),
        "straindesign_search_status": search_status,
        "straindesign_search_complete": search_status == "optimal",
        "straindesign_compress": True,
        "max_compressed_milp_solutions_requested": max_solutions,
        "max_solutions_semantics": (
            "caps compressed MILP solutions before straindesign network decompression; "
            "the returned reaction-design count can be larger"
        ),
        "n_compressed_milp_solutions": None,
        "n_compressed_milp_solutions_reason": (
            "straindesign does not expose the generated compressed-solution count"
        ),
        "n_decompressed_designs": n_decompressed_designs,
        "n_unique_decompressed_designs": len(ko_sets),
        "n_evaluated_designs": n_evaluated_designs,
        "n_method_eligible_designs": len(evaluated),
        "n_returned_designs_after_deduplication": len(deduplicated),
        "intervention_level": "reaction",
        "requires_gpr_resolution": True,
        # Backward-compatible historical keys.  ``n_milp_designs`` was a misnomer: it
        # counted evaluated decompressed designs, not compressed MILP solutions.
        "n_milp_designs": len(evaluated),
        "n_milp_designs_semantics": (
            "legacy alias for n_method_eligible_designs; not a compressed MILP count"
        ),
        "n_designs_after_deduplication": len(deduplicated),
    }
    return StrainDesignResult(
        method=method,
        product=product,
        designs=tuple(deduplicated),
        metadata=provenance,
    )


def _deduplicate_designs(designs: Iterable[StrainDesign]) -> list[StrainDesign]:
    """Keep the first design for each distinct knockout *set*, preserving rank order.

    ``max_solutions`` caps MILP solutions, not designs: straindesign decompresses one MILP
    solution into every member of the compressed reaction's coupled set, so the same
    intervention comes back repeatedly under different member ids. Collapsing on the
    knockout set as a set (order-insensitive) leaves one row per genuinely distinct design.
    """

    seen: set[frozenset[str]] = set()
    unique: list[StrainDesign] = []
    for design in designs:
        key = frozenset(design.knockouts)
        if key in seen:
            continue
        seen.add(key)
        unique.append(design)
    return unique


def optknock(
    model: Model,
    product: str,
    *,
    biomass: str | None = None,
    max_knockouts: int = 3,
    max_solutions: int = 5,
    min_growth: float = 0.05,
    condition: Condition | None = None,
    actionable_only: bool = True,
    seed: int = 0,
) -> StrainDesignResult:
    """OptKnock: knockout sets that maximize the product at maximum growth (optimistic).

    ``condition`` states the medium and aeration the search runs under, and is recorded in
    full in the provenance. Passing it is strongly preferred to conditioning the model
    beforehand: these two functions previously took no condition at all and depended
    silently on the caller's model state, which is how an aerobic design came to be
    published as the anaerobic answer.

    ``max_solutions`` caps **MILP solutions, not designs.** straindesign decompresses each
    MILP solution into every member of the compressed reaction's coupled set, so the number
    of designs returned is normally larger — ``max_solutions=5`` has returned 23. Designs
    are de-duplicated by knockout set before being returned; the raw and de-duplicated
    counts are both in the provenance.

    ``actionable_only`` restricts knockout candidates to reactions that have a GPR and are
    not boundary reactions, following Burgard et al.'s restriction to central metabolism.
    Set it to False to allow exchange knockouts, which are not realisable as gene deletions.

    ``seed`` is forwarded explicitly to straindesign's MILP backend.  Keeping it fixed pins
    the solver's search perturbation for reproducible reruns on the same solver/platform;
    omitting it no longer lets straindesign invent a different seed on every call.
    """

    return _strain_design(
        model,
        product,
        method="optknock",
        biomass=biomass,
        max_knockouts=max_knockouts,
        max_solutions=max_solutions,
        min_growth=min_growth,
        condition=condition,
        actionable_only=actionable_only,
        seed=seed,
    )


def robustknock(
    model: Model,
    product: str,
    *,
    biomass: str | None = None,
    max_knockouts: int = 3,
    max_solutions: int = 8,
    min_growth: float = 0.05,
    condition: Condition | None = None,
    actionable_only: bool = True,
    seed: int = 0,
) -> StrainDesignResult:
    """RobustKnock: maximize guaranteed product over all growth-maximal flux states.

    ``condition``, ``max_solutions``, ``actionable_only`` and ``seed`` behave exactly as in
    :func:`optknock`; in particular ``max_solutions`` caps MILP solutions, not designs.
    """

    return _strain_design(
        model,
        product,
        method="robustknock",
        biomass=biomass,
        max_knockouts=max_knockouts,
        max_solutions=max_solutions,
        min_growth=min_growth,
        condition=condition,
        actionable_only=actionable_only,
        seed=seed,
    )
