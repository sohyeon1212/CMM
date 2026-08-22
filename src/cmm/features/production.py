"""Production-design services: theoretical yield, production envelope, FSEOF.

These answer "how do we make more of a target product?" — the metabolic-engineering side of
the platform. They run on plain cobra models and return solver-neutral, serializable
results so the visualization layer and the GUI can render them without re-solving.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import cast
import warnings

import numpy as np
import pandas as pd
from cobra import Model, Reaction
from cobra.flux_analysis import flux_variability_analysis as _cobra_fva
from cobra.flux_analysis import production_envelope as _cobra_envelope

from cmm.core.condition import Condition
from cmm.core.media import MediumApplication
from cmm.core.provenance import run_provenance

_CO2_IDS = ("EX_co2_e", "EX_co2(e)", "EX_co2_e_")
_O2_IDS = ("EX_o2_e", "EX_o2(e)", "EX_o2_e_")

#: Tolerance on a carbon flux, in mmol C gDW-1 h-1, below which it is treated as zero.
_CARBON_TOL = 1e-6


@dataclass(frozen=True)
class CarbonUptake:
    """One carbon-containing uptake flux and the carbon it carries."""

    reaction_id: str
    uptake: float  # positive: mmol gDW-1 h-1 entering the model
    carbon_atoms: float  # carbon atoms per unit flux of the exchange

    @property
    def carbon_flux(self) -> float:
        """Carbon entering through this exchange, in mmol C gDW-1 h-1."""

        return self.uptake * self.carbon_atoms


@dataclass(frozen=True)
class ProductionYield:
    """Maximum product flux and yield per substrate at a fixed substrate uptake.

    The carbon ceiling is a plain carbon balance — the carbon actually taken up divided by
    the carbon in one unit of product — computed over **every** organic uptake, not the one
    nominated substrate, following COBRApy's ``carbon_yield`` convention of summing element
    counts over the exchanges that carry them. It is textbook carbon balance and is
    deliberately cited to no primary paper. CO2 is excluded from the ceiling and accounted
    separately: including it would make a CO2-inflated yield unfalsifiable by construction,
    which is exactly the defect ``co2_carbon_fraction`` exists to expose.
    """

    product: str
    substrate: str
    product_flux: float
    substrate_uptake: float
    molar_yield: float
    status: str
    aerobic: bool
    carbon_ceiling: float | None = (
        None  # max mol product / mol substrate from all organic C uptake
    )
    co2_exchange: float = 0.0  # net CO2 exchange (negative = fixation)
    product_carbon: float = 0.0  # carbon atoms per unit product flux
    carbon_uptake: tuple[
        CarbonUptake, ...
    ] = ()  # every organic carbon uptake, incl. substrate
    metadata: dict[str, object] = field(default_factory=dict)

    @property
    def co_substrates(self) -> tuple[CarbonUptake, ...]:
        """Carbon sources other than the nominated substrate that feed this yield."""

        return tuple(u for u in self.carbon_uptake if u.reaction_id != self.substrate)

    @property
    def co2_uptake(self) -> float:
        """Net CO2 *taken up*, as a positive flux (0.0 when CO2 is being released)."""

        return max(0.0, -self.co2_exchange)

    @property
    def co2_fixed(self) -> bool:
        """Whether net CO2 is being fixed (consumed), which can lift yield above the ceiling."""

        return self.co2_exchange < -1e-6

    @property
    def product_carbon_flux(self) -> float:
        """Carbon leaving in the product, in mmol C gDW-1 h-1."""

        return self.product_flux * self.product_carbon

    @property
    def co2_carbon_fraction(self) -> float:
        """Fraction of the product's carbon supplied by net CO2 uptake.

        This is the number ``co2_fixed`` alone never gave: a bare boolean cannot say whether
        CO2 contributed 1% or 30% of the product carbon.
        """

        if self.product_carbon_flux <= _CARBON_TOL:
            return 0.0
        return self.co2_uptake / self.product_carbon_flux

    @property
    def excess_carbon(self) -> float:
        """Product carbon above what the organic uptake supplies, in mmol C gDW-1 h-1."""

        organic = sum(u.carbon_flux for u in self.carbon_uptake)
        return max(0.0, self.product_carbon_flux - organic)

    @property
    def exceeds_carbon_ceiling(self) -> bool:
        if self.carbon_ceiling is None:
            return False
        return self.molar_yield > self.carbon_ceiling + 1e-6

    @property
    def co2_explains_excess(self) -> bool:
        """Whether CO2 uptake quantitatively accounts for the carbon above the ceiling."""

        return (
            self.exceeds_carbon_ceiling
            and self.co2_uptake >= self.excess_carbon - _CARBON_TOL
        )

    @property
    def carbon_imbalance(self) -> bool:
        """Yield above the ceiling that CO2 uptake does *not* account for.

        True means carbon is entering from somewhere neither the organic uptakes nor the CO2
        exchange explain — a modelling error, not a metabolic result.
        """

        return self.exceeds_carbon_ceiling and not self.co2_explains_excess

    def to_frame(self) -> pd.DataFrame:
        """One-row export of the yield and its carbon-balance diagnostics."""

        return pd.DataFrame(
            [
                {
                    "product": self.product,
                    "substrate": self.substrate,
                    "status": self.status,
                    "product_flux": self.product_flux,
                    "substrate_uptake": self.substrate_uptake,
                    "molar_yield": self.molar_yield,
                    "aerobic": self.aerobic,
                    "carbon_ceiling": self.carbon_ceiling,
                    "co2_exchange": self.co2_exchange,
                    "co2_uptake": self.co2_uptake,
                    "product_carbon": self.product_carbon,
                    "product_carbon_flux": self.product_carbon_flux,
                    "co2_carbon_fraction": self.co2_carbon_fraction,
                    "exceeds_carbon_ceiling": self.exceeds_carbon_ceiling,
                    "co2_explains_excess": self.co2_explains_excess,
                    "carbon_imbalance": self.carbon_imbalance,
                    "n_carbon_uptakes": len(self.carbon_uptake),
                }
            ]
        )

    def carbon_uptake_frame(self) -> pd.DataFrame:
        """One row per organic carbon uptake used for the carbon ceiling."""

        return pd.DataFrame(
            [
                {
                    "reaction_id": uptake.reaction_id,
                    "uptake": uptake.uptake,
                    "carbon_atoms": uptake.carbon_atoms,
                    "carbon_flux": uptake.carbon_flux,
                    "is_nominated_substrate": uptake.reaction_id == self.substrate,
                }
                for uptake in self.carbon_uptake
            ],
            columns=[
                "reaction_id",
                "uptake",
                "carbon_atoms",
                "carbon_flux",
                "is_nominated_substrate",
            ],
        )


@dataclass(frozen=True)
class EnvelopePoint:
    product_flux: float
    growth_min: float
    growth_max: float


@dataclass(frozen=True)
class ProductionEnvelope:
    product: str
    objective: str
    points: tuple[EnvelopePoint, ...]
    metadata: dict[str, object] = field(default_factory=dict)

    @property
    def max_growth(self) -> float:
        return max((p.growth_max for p in self.points), default=0.0)

    @property
    def max_product(self) -> float:
        return max((p.product_flux for p in self.points), default=0.0)

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            [(p.product_flux, p.growth_min, p.growth_max) for p in self.points],
            columns=["product_flux", "growth_min", "growth_max"],
        )


@dataclass(frozen=True)
class FseofResult:
    """Flux Scanning based on Enforced Objective Flux (Choi 2010)."""

    product: str
    biomass: str
    enforced_levels: tuple[float, ...]
    trends: (
        pd.DataFrame
    )  # index reaction_id; columns = enforced levels; + 'classification'
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        # frozen=True only blocks reassignment; copy the DataFrame so a caller mutating it
        # cannot corrupt this (or another caller's) cached result.
        object.__setattr__(self, "trends", self.trends.copy())

    def amplification_targets(self, *, actionable_only: bool = True) -> list[str]:
        """Return amplification candidates ranked by the FSEOF trend signal.

        FSEOF's CMM selection rule already requires a positive slope of flux magnitude.
        Rank the selected reactions by that slope (largest first), with reaction id as a
        deterministic tie-break.  Returning the model's reaction order made workflow
        ``top-N`` selection depend on SBML serialization rather than scientific evidence.
        """

        selected = self.trends["classification"] == "amplify"
        if actionable_only and "actionable" in self.trends:
            selected &= self.trends["actionable"].astype(bool)
        ranked = self.trends.loc[selected, ["slope"]].copy()
        ranked["_score"] = pd.to_numeric(ranked["slope"], errors="coerce").fillna(
            float("-inf")
        )
        ranked["_reaction_id"] = ranked.index.astype(str)
        ranked = ranked.sort_values(
            ["_score", "_reaction_id"],
            ascending=[False, True],
            kind="stable",
        )
        return ranked.index.tolist()

    def knockout_targets(self, *, actionable_only: bool = True) -> list[str]:
        selected = self.trends["classification"] == "knockdown"
        if actionable_only and "actionable" in self.trends:
            selected &= self.trends["actionable"].astype(bool)
        return self.trends.index[selected].tolist()

    def to_frame(self) -> pd.DataFrame:
        """Wide reaction-level export in explicit amplification rank order."""

        frame = self.trends.copy()
        rank_by_reaction = {
            reaction_id: rank
            for rank, reaction_id in enumerate(self.amplification_targets(), start=1)
        }
        frame.insert(
            0,
            "amplification_rank",
            pd.Series(rank_by_reaction, dtype="Int64"),
        )
        frame["_reaction_id"] = frame.index.astype(str)
        frame = frame.sort_values(
            ["amplification_rank", "_reaction_id"],
            ascending=[True, True],
            na_position="last",
            kind="stable",
        ).drop(columns="_reaction_id")
        frame.index.name = "reaction_id"
        return frame.reset_index()

    def trajectory_frame(self) -> pd.DataFrame:
        """Tidy export with one row per reaction and enforced-product level."""

        records: list[dict[str, object]] = []
        for reaction_id, row in self.trends.iterrows():
            for level in self.enforced_levels:
                records.append(
                    {
                        "reaction_id": reaction_id,
                        "enforced_product_flux": level,
                        "flux": float(row[level]),
                        "slope": float(row["slope"]),
                        "classification": str(row["classification"]),
                        "actionable": bool(row.get("actionable", False)),
                    }
                )
        return pd.DataFrame(
            records,
            columns=[
                "reaction_id",
                "enforced_product_flux",
                "flux",
                "slope",
                "classification",
                "actionable",
            ],
        )


@dataclass(frozen=True)
class FvseofResult:
    """Flux Variability Scanning based on Enforced Objective Flux (Park et al. 2012).

    Like FSEOF but runs FVA at each enforced product level, so each reaction has a flux
    *range* per step, summarised by Park's two statistics: ``V_avg = (V'max + V'min)/2``
    (here ``mean``) and ``l_sol = V'max - V'min`` (here ``capacity``).

    **Classification follows Park's nine types**, taken from the *joint* sign of the change
    in V_avg and the change in l_sol across the scan — see :func:`_park_flux_type`. Types
    1-3 (V_avg rising) are the amplification candidates, types 4-6 (V_avg falling) the
    knockdown candidates, types 7-9 (no clear V_avg trend) neither. Within the amplification
    band Park prioritise "reactions with smaller values of l_sol", which is the order
    :meth:`amplification_targets` returns.

    **The ``robust`` flag is CMM's own construct and is not in Park et al.** It marks
    reactions whose *forced minimum* |flux| rises monotonically — the reaction cannot avoid
    carrying more flux, which single-FBA FSEOF cannot distinguish. It is useful and is kept,
    but it is not Park's variability criterion: Park's variability signal is the change in
    l_sol, which is what the nine types encode. The two are related but not equivalent —
    a forced minimum can rise while l_sol also rises, and l_sol can shrink while the forced
    minimum stays at 0 because the interval still straddles zero.
    """

    product: str
    biomass: str
    enforced_levels: tuple[float, ...]
    mean: pd.DataFrame  # reaction x level: V_avg, Park's (V'max + V'min)/2
    forced: (
        pd.DataFrame
    )  # reaction x level: forced minimum |flux| (0 if the range spans 0); CMM's own
    capacity: pd.DataFrame  # reaction x level: l_sol, Park's V'max - V'min
    classification: (
        pd.Series
    )  # reaction -> amplify (types 1-3) / knockdown (4-6) / none (7-9)
    park_type: pd.Series  # reaction -> Park nine-type index, 1-9
    robust: (
        pd.Series
    )  # reaction -> bool (forced-min |flux| monotonically rises); CMM's own
    # Park's q_slope regresses the *signed* flux; this regresses |V_avg|, matching CMM's
    # magnitude convention for reactions that run in reverse. Different quantity, same role.
    slope: pd.Series  # reaction -> linear slope of |V_avg| versus enforced product
    capacity_slope: pd.Series  # reaction -> slope of l_sol (Park's delta l_sol)
    actionable: pd.Series  # reaction -> has GPR and is internal/non-objective
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "mean", self.mean.copy())
        object.__setattr__(self, "forced", self.forced.copy())
        object.__setattr__(self, "capacity", self.capacity.copy())

    def _mean_capacity(self) -> pd.Series:
        """Mean l_sol per reaction across the scan, Park's amplification priority key."""

        return self.capacity.mean(axis=1)

    def amplification_targets(self) -> list[str]:
        """Rank actionable type 1-3 reactions by Park's mean l_sol priority.

        Smaller mean ``l_sol`` ranks first; reaction id is the deterministic tie-break so
        a workflow shortlist does not depend on the model's serialization order.
        """

        selected = (self.classification == "amplify") & self.actionable.astype(bool)
        ranked = self._mean_capacity()[selected].rename("mean_capacity").to_frame()
        ranked["reaction_id"] = ranked.index.astype(str)
        ranked = ranked.sort_values(
            ["mean_capacity", "reaction_id"],
            ascending=[True, True],
            kind="stable",
        )
        return ranked.index.tolist()

    def knockout_targets(self) -> list[str]:
        selected = (self.classification == "knockdown") & self.actionable.astype(bool)
        return self.classification.index[selected].tolist()

    def targets_of_type(self, *types: int) -> list[str]:
        """Actionable reactions in the given Park types, in scan order."""

        selected = self.park_type.isin(types) & self.actionable.astype(bool)
        return self.park_type.index[selected].tolist()

    def robust_targets(self) -> list[str]:
        """Amplification targets also flagged ``robust`` — CMM's flag, not Park's."""

        amp = set(self.amplification_targets())
        return [rid for rid in self.robust.index[self.robust] if rid in amp]

    def to_frame(self) -> pd.DataFrame:
        """Reaction-level classification and ranking quantities."""

        reaction_ids = list(self.mean.index)
        rank_by_reaction = {
            reaction_id: rank
            for rank, reaction_id in enumerate(self.amplification_targets(), start=1)
        }
        frame = pd.DataFrame(
            {
                "reaction_id": reaction_ids,
                "amplification_rank": pd.Series(
                    [rank_by_reaction.get(reaction_id) for reaction_id in reaction_ids],
                    dtype="Int64",
                ),
                "classification": [
                    str(self.classification.loc[r]) for r in reaction_ids
                ],
                "park_type": [int(self.park_type.loc[r]) for r in reaction_ids],
                "robust": [bool(self.robust.loc[r]) for r in reaction_ids],
                "slope": [float(self.slope.loc[r]) for r in reaction_ids],
                "capacity_slope": [
                    float(self.capacity_slope.loc[r]) for r in reaction_ids
                ],
                "mean_capacity": [
                    float(self._mean_capacity().loc[r]) for r in reaction_ids
                ],
                "actionable": [bool(self.actionable.loc[r]) for r in reaction_ids],
            }
        )
        return frame.sort_values(
            ["amplification_rank", "reaction_id"],
            ascending=[True, True],
            na_position="last",
            kind="stable",
            ignore_index=True,
        )

    def trajectory_frame(self) -> pd.DataFrame:
        """Tidy export of V_avg, forced minimum and l_sol across the scan."""

        records: list[dict[str, object]] = []
        for reaction_id in self.mean.index:
            for level in self.enforced_levels:
                records.append(
                    {
                        "reaction_id": reaction_id,
                        "enforced_product_flux": level,
                        "mean_flux": float(self.mean.loc[reaction_id, level]),
                        "forced_min_abs_flux": float(
                            self.forced.loc[reaction_id, level]
                        ),
                        "capacity": float(self.capacity.loc[reaction_id, level]),
                        "classification": str(self.classification.loc[reaction_id]),
                        "park_type": int(self.park_type.loc[reaction_id]),
                        "robust": bool(self.robust.loc[reaction_id]),
                        "actionable": bool(self.actionable.loc[reaction_id]),
                    }
                )
        return pd.DataFrame(
            records,
            columns=[
                "reaction_id",
                "enforced_product_flux",
                "mean_flux",
                "forced_min_abs_flux",
                "capacity",
                "classification",
                "park_type",
                "robust",
                "actionable",
            ],
        )


def _carbon_count(rxn) -> float:
    """Carbon atoms carried per unit flux of ``rxn``, weighted by stoichiometry.

    Element counts are multiplied by their stoichiometric coefficients and summed, as
    COBRApy's ``_reaction_elements`` does, then taken in magnitude so the result is a
    per-unit-flux atom count rather than a signed balance. Ignoring the coefficients — as
    this did before — counts a mass-balanced internal reaction as carbon-rich (PFK scored
    32 carbons, giving it a meaningless "carbon ceiling") and mis-counts any exchange whose
    coefficient is not -1.
    """

    total = 0.0
    for met, coefficient in rxn.metabolites.items():
        try:
            total += float(coefficient) * float(met.elements.get("C", 0))
        except (AttributeError, TypeError, ValueError):
            continue
    return abs(total)


def _detect_substrate(model: Model, unlimited: float = -999.0) -> str:
    """Pick the main carbon source: the limiting *carbon-containing* uptake exchange.

    The substrate is the carbon exchange with an explicitly restricted uptake bound (e.g.
    glucose at -10), not an unconstrained one and not a non-carbon nutrient (phosphate,
    sulfate, ...). Exchanges left at the default ~-1000 (CO2, H2O, NH4, O2) are availability
    constraints, not the carbon source.
    """

    carbon_limited = [
        rxn
        for rxn in model.exchanges
        if unlimited < rxn.lower_bound < 0.0 and _carbon_count(rxn) > 0
    ]
    if carbon_limited:
        return min(carbon_limited, key=lambda r: r.lower_bound).id
    # Fall back to the carbon exchange actually taken up at the wild-type optimum,
    # inside a context guard so the probe does not disturb the caller's model state.
    with model:
        solution = model.optimize()
    if solution.status != "optimal":
        raise ValueError(
            f"cannot detect substrate from a model with solver status {solution.status!r}"
        )
    uptaken = [
        (rxn.id, solution.fluxes[rxn.id])
        for rxn in model.exchanges
        if solution.fluxes[rxn.id] < 0 and _carbon_count(rxn) > 0
    ]
    if not uptaken:
        raise ValueError(
            "no carbon uptake exchange found; specify substrate explicitly"
        )
    return min(uptaken, key=lambda x: x[1])[0]


def _exchange(model: Model, candidates: Iterable[str]):
    """The first of ``candidates`` present in the model, or None."""

    for reaction_id in candidates:
        if reaction_id in model.reactions:
            return model.reactions.get_by_id(reaction_id)
    return None


def _apply_condition(model: Model, *, condition: Condition | None) -> None:
    """Apply the run condition to ``model`` in place (caller supplies a ``with model:``).

    ``condition`` is the *only* way to state the environment of a production run as of
    0.4.0. The ``aerobic=True|False`` flag it replaces was a second, redundant way of saying
    what the medium already says, and because ``optknock``/``robustknock`` never accepted it
    a caller could pass ``aerobic=False`` and still receive an aerobic design in silence.
    Anaerobiosis is now stated the same way as everything else — a
    :class:`~cmm.core.condition.Condition` closing the oxygen exchange, on top of a medium
    (``glucose_anaerobic``) whose CO2 uptake is already closed by ``cmm.core.media``.
    """

    if condition is not None:
        condition.apply_to(model)


def _condition_provenance(
    model: Model, condition: Condition | None
) -> dict[str, object]:
    """Full record of the condition a run actually solved under.

    A result file must state its own conditions rather than leave a reader to reconstruct
    them from a model fingerprint: the medium as applied, the oxygen exchange bounds, and
    every carbon uptake with its rate. The medium half is written through
    :meth:`~cmm.core.media.MediumApplication.to_provenance` so its shape matches what
    ``cmm.core.media`` already records. Call this *after* the condition has been applied.
    """

    application = MediumApplication(
        medium=condition.name if condition is not None else "model as loaded",
        applied={
            reaction_id: float(rate) for reaction_id, rate in model.medium.items()
        },
    )
    oxygen = _exchange(model, _O2_IDS)
    return {
        "condition": condition.name if condition is not None else None,
        "medium": application.to_provenance(),
        "oxygen_exchange": (
            None
            if oxygen is None
            else {
                "reaction_id": oxygen.id,
                "lower_bound": float(oxygen.lower_bound),
                "upper_bound": float(oxygen.upper_bound),
            }
        ),
        "aerobic": _is_aerobic(model),
        "carbon_uptake": [
            {"reaction_id": u.reaction_id, "uptake_rate": u.uptake}
            for u in _carbon_uptake_bounds(model)
        ],
    }


def _record_condition(
    provenance: dict[str, object], applied: dict[str, object]
) -> None:
    """Attach an applied-condition record to a ``run_provenance`` parameter block."""

    cast(dict[str, object], provenance["parameters"])["applied_condition"] = applied


def _is_aerobic(model: Model) -> bool:
    """Whether the oxygen exchange is open for uptake (True when the model has none)."""

    oxygen = _exchange(model, _O2_IDS)
    return oxygen is None or oxygen.lower_bound < -1e-9


def _carbon_uptake_bounds(model: Model) -> tuple[CarbonUptake, ...]:
    """Carbon-containing exchanges open for uptake, at their *bound* (not their flux)."""

    co2 = _exchange(model, _CO2_IDS)
    co2_id = None if co2 is None else co2.id
    return tuple(
        CarbonUptake(rxn.id, abs(float(rxn.lower_bound)), carbon)
        for rxn in model.exchanges
        if rxn.lower_bound < 0.0
        and rxn.id != co2_id
        and (carbon := _carbon_count(rxn)) > 0
    )


def _carbon_uptake_fluxes(model: Model, solution) -> tuple[CarbonUptake, ...]:
    """Carbon-containing exchanges actually *consumed* in a solution, excluding CO2.

    CO2 is excluded deliberately and reported on its own: folding it into the carbon
    ceiling would raise the ceiling by exactly the amount that inflated the yield, so no
    guard could ever fire.
    """

    co2 = _exchange(model, _CO2_IDS)
    co2_id = None if co2 is None else co2.id
    uptakes = []
    for rxn in model.exchanges:
        if rxn.id == co2_id:
            continue
        flux = float(solution.fluxes[rxn.id])
        if flux >= -_CARBON_TOL:
            continue
        carbon = _carbon_count(rxn)
        if carbon > 0:
            uptakes.append(CarbonUptake(rxn.id, -flux, carbon))
    return tuple(sorted(uptakes, key=lambda u: -u.carbon_flux))


def _net_co2(solution, model: Model) -> float:
    for co2_id in _CO2_IDS:
        if co2_id in model.reactions:
            return float(solution.fluxes[co2_id])
    return 0.0


def _require_boundary(model: Model, reaction_id: str, role: str):
    """Fetch an exchange reaction, refusing an internal one.

    ``theoretical_yield(model, "PFK")`` used to return 17.66 mol/mol against a ceiling of
    0.1875 — a number with no meaning, because a yield is defined per unit of *exchanged*
    product. Refuse it at the door.
    """

    reaction = model.reactions.get_by_id(reaction_id)
    if not reaction.boundary:
        raise ValueError(
            f"{role} {reaction_id!r} is not a boundary (exchange) reaction: "
            f"'{reaction.reaction}'. A molar yield is defined per unit of product leaving "
            "the system, so both product and substrate must be exchange reactions; pass the "
            "corresponding EX_ reaction instead."
        )
    return reaction


def theoretical_yield(
    model: Model,
    product: str,
    substrate: str | None = None,
    *,
    condition: Condition | None = None,
) -> ProductionYield:
    """Maximum molar yield of ``product`` per ``substrate`` taken up.

    Maximizes the product exchange with the substrate uptake fixed exactly (so the yield
    denominator is deterministic). Both ``product`` and ``substrate`` must be exchange
    reactions; an internal reaction raises.

    ``carbon_ceiling`` is a carbon balance over **every** organic uptake the solution
    consumes, not the nominated substrate alone, so a co-fed model is measured against the
    carbon it actually receives. CO2 is held out of that balance and reported separately:
    ``co2_carbon_fraction`` states how much of the product's carbon net CO2 uptake supplies,
    and ``carbon_imbalance`` flags an excess that CO2 does *not* explain. A yield obtained by
    fixing CO2 also warns, because ``co2_fixed`` on its own let a 15.9%-inflated yield pass
    every guard simply by sitting below the ceiling.

    ``condition`` is the supported way to state the aeration, medium and substrate of a run;
    the applied condition is recorded in full under
    ``metadata["parameters"]["applied_condition"]``.
    """

    provenance = run_provenance(
        model,
        method="theoretical_yield",
        product=product,
        substrate=substrate,
        condition=condition.name if condition is not None else None,
    )
    with model:
        _apply_condition(model, condition=condition)
        product_rxn = _require_boundary(model, product, "product")
        substrate = substrate or _detect_substrate(model)
        substrate_rxn = _require_boundary(model, substrate, "substrate")
        substrate_uptake = abs(substrate_rxn.lower_bound)
        if substrate_uptake <= 1e-9:
            # A closed substrate (lower_bound 0) gives a 0 denominator and a NaN yield. This
            # is plausible (e.g. a genome-scale medium with glucose shut, see rule R14), so
            # fail loudly with an actionable message instead of returning a silent NaN.
            raise ValueError(
                f"substrate '{substrate}' has no uptake capacity (lower_bound "
                f"{substrate_rxn.lower_bound}); open its uptake or choose another substrate "
                f"before computing a theoretical yield."
            )
        substrate_rxn.bounds = (-substrate_uptake, -substrate_uptake)
        _record_condition(provenance, _condition_provenance(model, condition))
        run_aerobic = _is_aerobic(model)
        model.objective = product_rxn
        model.objective_direction = "max"
        solution = model.optimize()
        if solution.status != "optimal":
            raise ValueError(
                f"theoretical-yield optimization for {product!r} is {solution.status}"
            )
        product_flux = float(solution.fluxes[product])
        molar_yield = product_flux / substrate_uptake

        carbon_uptake = _carbon_uptake_fluxes(model, solution)
        product_carbon = _carbon_count(product_rxn)
        total_carbon_in = sum(u.carbon_flux for u in carbon_uptake)
        carbon_ceiling = (
            total_carbon_in / (product_carbon * substrate_uptake)
            if product_carbon > 0 and total_carbon_in > 0
            else None
        )
        result = ProductionYield(
            product=product,
            substrate=substrate,
            product_flux=product_flux,
            substrate_uptake=substrate_uptake,
            molar_yield=molar_yield,
            status=solution.status,
            aerobic=run_aerobic,
            carbon_ceiling=carbon_ceiling,
            co2_exchange=_net_co2(solution, model),
            product_carbon=product_carbon,
            carbon_uptake=carbon_uptake,
            metadata=provenance,
        )
    _warn_on_co2_contribution(result)
    return result


def _warn_on_co2_contribution(result: ProductionYield) -> None:
    """Say so when CO2 uptake is holding up the yield, and by how much.

    ``co2_fixed`` alone never fired as a guard: it is a boolean that a caller reading a
    single yield number does not see, and the only active check compared the yield with a
    ceiling it could sit below while still being CO2-inflated.
    """

    if result.co2_carbon_fraction <= 1e-9:
        return
    detail = (
        f"{result.product} yield {result.molar_yield:.4f} mol/mol {result.substrate} "
        f"draws {result.co2_carbon_fraction:.1%} of its product carbon from net CO2 uptake "
        f"({result.co2_uptake:.4f} mmol gDW-1 h-1). Close the CO2 exchange for uptake "
        "unless the condition genuinely supplies CO2."
    )
    if result.carbon_imbalance:
        detail += (
            f" The yield also exceeds the carbon ceiling {result.carbon_ceiling} by more "
            f"than CO2 uptake accounts for ({result.excess_carbon:.4f} mmol C gDW-1 h-1 "
            "unexplained) - check the model's carbon balance."
        )
    warnings.warn(detail, UserWarning, stacklevel=3)


def production_envelope(
    model: Model,
    product: str,
    *,
    objective: str | None = None,
    substrate: str | None = None,
    condition: Condition | None = None,
    points: int = 20,
) -> ProductionEnvelope:
    """Growth-vs-product envelope: the feasible product range at each growth rate.

    This is the production envelope of Burgard et al. (2003), who introduced it alongside
    OptKnock as the projection of the flux cone onto the (growth, product) plane. It is
    **not** a phenotypic phase plane: a phase plane varies two *substrate uptake* rates and
    partitions them into regions of constant shadow price, whereas this varies the product
    flux and reports the growth range at each level. "Production envelope" is the term used
    by Feist et al. (2010) for exactly this object.
    """

    if points < 2:
        raise ValueError("points must be at least 2")
    objective = objective or _objective_reaction(model)
    provenance = run_provenance(
        model,
        method="production_envelope",
        product=product,
        objective=objective,
        substrate=substrate,
        condition=condition.name if condition is not None else None,
        points=points,
    )
    with model:
        _apply_condition(model, condition=condition)
        substrate = substrate or _detect_substrate(model)
        _record_condition(provenance, _condition_provenance(model, condition))
        frame = _cobra_envelope(
            model,
            [product],
            objective=objective,
            carbon_sources=substrate,
            points=points,
        )
    envelope_points = tuple(
        EnvelopePoint(
            product_flux=float(row[product]),
            growth_min=float(row["flux_minimum"]),
            growth_max=float(row["flux_maximum"]),
        )
        for _, row in frame.iterrows()
    )
    return ProductionEnvelope(
        product=product,
        objective=objective,
        points=envelope_points,
        metadata=provenance,
    )


def _validate_scan_parameters(
    *,
    n_steps: int,
    fraction_min: float,
    fraction_max: float,
    tol: float,
    biomass_fraction: float | None = None,
) -> None:
    if n_steps < 2:
        raise ValueError("n_steps must be at least 2")
    if not 0.0 <= fraction_min < fraction_max <= 1.0:
        raise ValueError(
            "scan fractions must satisfy 0 <= fraction_min < fraction_max <= 1"
        )
    if tol < 0:
        raise ValueError("tol must be non-negative")
    if biomass_fraction is not None and not 0.0 < biomass_fraction <= 1.0:
        raise ValueError("biomass_fraction must be in (0, 1]")


def _initial_product_flux(
    model: Model,
    product: str,
    biomass: str,
    *,
    condition: Condition | None,
) -> float:
    with model:
        _apply_condition(model, condition=condition)
        model.objective = model.reactions.get_by_id(biomass)
        model.objective_direction = "max"
        solution = model.optimize()
        if solution.status != "optimal":
            raise ValueError(
                f"cannot determine initial product flux: biomass solve is {solution.status}"
            )
        return float(solution.fluxes[product])


def _scan_levels(
    initial: float,
    maximum: float,
    *,
    n_steps: int,
    fraction_min: float,
    fraction_max: float,
    tol: float,
) -> np.ndarray:
    if maximum < initial - tol:
        raise ValueError(
            f"maximum product flux ({maximum:g}) is below its growth-optimal initial flux "
            f"({initial:g}); check product direction and model bounds"
        )
    if maximum - initial <= tol:
        return np.asarray([maximum], dtype=float)
    fractions = np.linspace(fraction_min, fraction_max, n_steps)
    return initial + fractions * (maximum - initial)


#: Placeholder gene identifiers that stand for "spontaneous / no catalysing gene" rather than
#: a real gene. BiGG models (``e_coli_core``, ``iJO1366``) use ``s0001``. A reaction whose GPR
#: contains only these has no gene a strain engineer could delete, so it is not actionable even
#: though ``reaction.genes`` is non-empty.
PSEUDO_GENE_IDS = frozenset({"s0001", "spontaneous", "SPONTANEOUS"})


def _has_real_gene(rxn: Reaction) -> bool:
    """True when the reaction's GPR names at least one non-placeholder gene."""

    return any(gene.id not in PSEUDO_GENE_IDS for gene in rxn.genes)


def _actionable_reaction(model: Model, rid: str, product: str, biomass: str) -> bool:
    rxn = model.reactions.get_by_id(rid)
    return _has_real_gene(rxn) and not rxn.boundary and rid not in {product, biomass}


def _linear_slope(levels: np.ndarray, values: np.ndarray) -> float:
    if len(levels) < 2 or len(values) != len(levels) or not np.isfinite(values).all():
        return float("nan")
    if float(np.ptp(levels)) <= 1e-12:
        return 0.0
    return float(np.polyfit(levels, values, 1)[0])


def _add_linear_flux_couplings(
    model: Model, couplings: Iterable[Mapping[str, float]] | None
) -> int:
    """Add caller-supplied linear flux equalities ``sum(coeff * flux) = 0``.

    **These are not Park et al.'s grouping-reaction (GR) constraints**, and were renamed so
    the name stops implying that they are. Park derive reaction groups automatically from
    STRING genomic-context analysis and impose two constructs on them: binary on/off
    indicators tying the two members' activity (``y_v1 = y_v2``, a MILP construct), and a
    flux-scale *inequality* ``|v1 - mean| + |v2 - mean| <= delta`` with ``delta = 0.3`` on
    fluxes normalised by the carbon uptake rate. CMM takes neither: it takes whatever linear
    equalities the caller supplies, on raw fluxes, with no automated grouping, and stays a
    pure LP. Automated STRING-based grouping is tracked as future work in the roadmap.
    """

    count = 0
    for count, coefficients in enumerate(couplings or (), start=1):
        if len(coefficients) < 2:
            raise ValueError("each linear flux coupling needs at least two reactions")
        expression = 0.0
        for rid, coefficient in coefficients.items():
            if rid not in model.reactions:
                raise KeyError(
                    f"linear flux coupling references unknown reaction {rid!r}"
                )
            if not np.isfinite(coefficient) or coefficient == 0:
                raise ValueError(
                    "linear flux coupling coefficients must be finite and non-zero"
                )
            expression += coefficient * model.reactions.get_by_id(rid).flux_expression
        constraint = model.problem.Constraint(
            expression, lb=0.0, ub=0.0, name=f"_fvseof_coupling_{count}"
        )
        model.add_cons_vars([constraint])
    return count


def fseof(
    model: Model,
    product: str,
    biomass: str | None = None,
    *,
    n_steps: int = 10,
    fraction_min: float = 0.1,
    fraction_max: float = 0.9,
    condition: Condition | None = None,
    reactions: Iterable[str] | None = None,
    tol: float = 1e-3,
) -> FseofResult:
    """Scan enforced product flux and classify each reaction's flux trend (Choi 2010).

    At each enforced product level (a fraction of the theoretical maximum) the product flux
    is fixed and biomass is maximized, exactly as Choi et al. specify.

    **The selection criterion is CMM's own, not Choi et al.'s**, and is deliberately kept
    (see the deviation register). Choi select on ``|v|max > |v_initial|`` and
    ``v_max * v_min >= 0``, baselined at the growth-optimal wild-type flux. CMM instead
    requires the flux *magnitude* to rise from the first to the last scan level **and** the
    linear regression of that magnitude against enforced product to have positive slope,
    with no direction reversal anywhere in the scan; the baseline is the first scan level,
    not ``v_initial``. It does **not** test monotonicity at every step — the provenance
    string used to claim it did. CMM's rule is the stricter of the two: on anaerobic
    ``e_coli_core``/succinate Choi's additionally admits the acetate branch (``ACKr``,
    ``ACt2r``, ``EX_ac_e``, ``PTAr``), which 12 of 18 OptKnock and 8 of 41 RobustKnock
    designs for the same target *delete*.

    Classifying on magnitude rather than signed flux keeps reactions that operate in the
    reverse direction — e.g. the reductive succinate pathway — correct.
    """

    provenance = run_provenance(
        model,
        method="fseof",
        product=product,
        biomass=biomass,
        n_steps=n_steps,
        fraction_min=fraction_min,
        fraction_max=fraction_max,
        condition=condition.name if condition is not None else None,
        tol=tol,
    )
    _validate_scan_parameters(
        n_steps=n_steps,
        fraction_min=fraction_min,
        fraction_max=fraction_max,
        tol=tol,
    )
    biomass = biomass or _objective_reaction(model)
    yield_result = theoretical_yield(model, product, condition=condition)
    max_product = yield_result.product_flux
    initial_product = _initial_product_flux(
        model, product, biomass, condition=condition
    )
    levels = _scan_levels(
        initial_product,
        max_product,
        n_steps=n_steps,
        fraction_min=fraction_min,
        fraction_max=fraction_max,
        tol=tol,
    )

    rxn_ids = (
        list(reactions) if reactions is not None else [r.id for r in model.reactions]
    )
    columns: dict[float, dict[str, float]] = {}
    with model:
        _apply_condition(model, condition=condition)
        _record_condition(provenance, _condition_provenance(model, condition))
        product_rxn = model.reactions.get_by_id(product)
        model.objective = model.reactions.get_by_id(biomass)
        model.objective_direction = "max"
        for level in levels:
            product_rxn.bounds = (float(level), float(level))  # fix flux per Choi 2010
            solution = model.optimize()
            if solution.status != "optimal":
                columns[float(level)] = {rid: float("nan") for rid in rxn_ids}
                continue
            columns[float(level)] = {
                rid: float(solution.fluxes[rid]) for rid in rxn_ids
            }

    trends = pd.DataFrame(columns, index=rxn_ids)
    scan_columns = list(columns.keys())
    level_array = np.asarray(scan_columns, dtype=float)
    slopes = {
        rid: _linear_slope(
            level_array,
            np.abs(trends.loc[rid, scan_columns].to_numpy(dtype=float)),
        )
        for rid in rxn_ids
    }
    trends["slope"] = pd.Series(slopes)
    trends["classification"] = [
        _classify_trend(trends.loc[rid, scan_columns].to_numpy(dtype=float), tol)
        for rid in rxn_ids
    ]
    trends["actionable"] = [
        _actionable_reaction(model, rid, product, biomass) for rid in rxn_ids
    ]
    return FseofResult(
        product=product,
        biomass=biomass,
        enforced_levels=tuple(float(x) for x in levels),
        trends=trends,
        metadata={
            **provenance,
            "initial_product": initial_product,
            "max_product": max_product,
            "aerobic": yield_result.aerobic,
            "n_failed_levels": sum(
                not np.isfinite(trends[level].to_numpy(dtype=float)).all()
                for level in scan_columns
            ),
            # Names what is computed: endpoint difference AND positive linear-regression
            # slope of |flux|, with no direction reversal. Monotonicity is NOT tested - the
            # previous string said it was, and three amplify calls were non-monotonic.
            "criterion": (
                "endpoint_difference_and_linear_slope_of_flux_magnitude_"
                "without_direction_reversal"
            ),
            "criterion_source": "CMM (not Choi et al. 2010)",
        },
    )


def _forced_min_magnitude(minimum: float, maximum: float, tol: float = 1e-9) -> float:
    """The smallest |flux| the reaction must carry given its feasible [min, max] range."""

    if minimum > tol:
        return minimum
    if maximum < -tol:
        return -maximum
    return 0.0  # the range spans zero, so the reaction can be off


def fvseof(
    model: Model,
    product: str,
    biomass: str | None = None,
    *,
    n_steps: int = 10,
    fraction_min: float = 0.1,
    fraction_max: float = 0.9,
    biomass_fraction: float = 0.95,
    condition: Condition | None = None,
    reactions: Iterable[str] | None = None,
    linear_flux_couplings: Iterable[Mapping[str, float]] | None = None,
    tol: float = 1e-3,
) -> FvseofResult:
    """Flux Variability Scanning based on Enforced Objective Flux (Park et al. 2012).

    At each enforced product level the product flux is fixed and biomass is held by equality
    at ``biomass_fraction`` of its maximum, recomputed per step; FVA then gives every
    reaction's flux range, from which V_avg and l_sol follow. Reactions are classified into
    Park's nine types on the joint sign of the change in V_avg and the change in l_sol, and
    types 1-3 are the amplification candidates, ordered by ascending l_sol as Park prioritise
    them. ``n_steps`` defaults to 10, Park's stated minimum.

    ``linear_flux_couplings`` are **not** Park's STRING-derived grouping-reaction constraints
    — see :func:`_add_linear_flux_couplings` — and the ``robust`` flag is CMM's, not Park's.
    """

    provenance = run_provenance(
        model,
        method="fvseof",
        product=product,
        biomass=biomass,
        n_steps=n_steps,
        fraction_min=fraction_min,
        fraction_max=fraction_max,
        biomass_fraction=biomass_fraction,
        condition=condition.name if condition is not None else None,
        tol=tol,
    )
    _validate_scan_parameters(
        n_steps=n_steps,
        fraction_min=fraction_min,
        fraction_max=fraction_max,
        tol=tol,
        biomass_fraction=biomass_fraction,
    )
    biomass = biomass or _objective_reaction(model)
    yield_result = theoretical_yield(model, product, condition=condition)
    max_product = yield_result.product_flux
    initial_product = _initial_product_flux(
        model, product, biomass, condition=condition
    )
    levels = _scan_levels(
        initial_product,
        max_product,
        n_steps=n_steps,
        fraction_min=fraction_min,
        fraction_max=fraction_max,
        tol=tol,
    )

    rxn_ids = (
        list(reactions) if reactions is not None else [r.id for r in model.reactions]
    )
    mean_cols: dict[float, dict[str, float]] = {}
    forced_cols: dict[float, dict[str, float]] = {}
    capacity_cols: dict[float, dict[str, float]] = {}
    with model:
        _apply_condition(model, condition=condition)
        _record_condition(provenance, _condition_provenance(model, condition))
        product_rxn = model.reactions.get_by_id(product)
        model.objective = model.reactions.get_by_id(biomass)
        model.objective_direction = "max"
        n_couplings = _add_linear_flux_couplings(model, linear_flux_couplings)
        for level in levels:
            level = float(level)
            with model:
                product_rxn.bounds = (level, level)
                try:
                    growth_optimum = model.slim_optimize(error_value=None)
                    if growth_optimum is None or not np.isfinite(growth_optimum):
                        raise ValueError("growth optimization is infeasible")
                    target_growth = biomass_fraction * float(growth_optimum)
                    biological_objective = model.objective.expression
                    growth_constraint = model.problem.Constraint(
                        biological_objective,
                        lb=target_growth,
                        ub=target_growth,
                        name="_fvseof_growth_level",
                    )
                    model.add_cons_vars([growth_constraint])
                    model.objective = model.problem.Objective(0.0, direction="min")
                    table = _cobra_fva(
                        model,
                        reaction_list=rxn_ids,
                        fraction_of_optimum=0.0,
                        processes=1,
                    )
                except Exception:
                    mean_cols[level] = {rid: float("nan") for rid in rxn_ids}
                    forced_cols[level] = {rid: float("nan") for rid in rxn_ids}
                    capacity_cols[level] = {rid: float("nan") for rid in rxn_ids}
                    continue
            mean_cols[level] = {
                rid: 0.5 * (float(row["minimum"]) + float(row["maximum"]))
                for rid, row in table.iterrows()
            }
            forced_cols[level] = {
                rid: _forced_min_magnitude(float(row["minimum"]), float(row["maximum"]))
                for rid, row in table.iterrows()
            }
            capacity_cols[level] = {
                rid: float(row["maximum"]) - float(row["minimum"])
                for rid, row in table.iterrows()
            }

    mean = pd.DataFrame(mean_cols, index=rxn_ids)
    forced = pd.DataFrame(forced_cols, index=rxn_ids)
    capacity = pd.DataFrame(capacity_cols, index=rxn_ids)
    levels_key = list(mean_cols.keys())
    level_array = np.asarray(levels_key, dtype=float)
    park_type = pd.Series(
        [
            _park_flux_type(
                mean.loc[rid, levels_key].to_numpy(dtype=float),
                capacity.loc[rid, levels_key].to_numpy(dtype=float),
                tol,
            )
            for rid in rxn_ids
        ],
        index=rxn_ids,
        dtype=int,
    )
    classification = park_type.map(_PARK_TYPE_CLASSIFICATION)
    robust = pd.Series(
        [
            _rises_monotonically(forced.loc[rid, levels_key].to_numpy(), tol)
            for rid in rxn_ids
        ],
        index=rxn_ids,
    )
    slope = pd.Series(
        {
            rid: _linear_slope(
                level_array, np.abs(mean.loc[rid, levels_key].to_numpy(dtype=float))
            )
            for rid in rxn_ids
        }
    )
    capacity_slope = pd.Series(
        {
            rid: _linear_slope(
                level_array, capacity.loc[rid, levels_key].to_numpy(dtype=float)
            )
            for rid in rxn_ids
        }
    )
    actionable = pd.Series(
        {rid: _actionable_reaction(model, rid, product, biomass) for rid in rxn_ids},
        dtype=bool,
    )
    return FvseofResult(
        product=product,
        biomass=biomass,
        # Use the actual scan keys so enforced_levels always matches the result columns
        # (a zero-yield product collapses all levels to a single 0.0).
        enforced_levels=tuple(float(x) for x in levels_key),
        mean=mean,
        forced=forced,
        capacity=capacity,
        classification=classification,
        park_type=park_type,
        robust=robust,
        slope=slope,
        capacity_slope=capacity_slope,
        actionable=actionable,
        metadata={
            **provenance,
            "initial_product": initial_product,
            "max_product": max_product,
            "biomass_fraction": biomass_fraction,
            "n_linear_flux_couplings": n_couplings,
            "n_failed_levels": sum(
                not np.isfinite(mean[level].to_numpy(dtype=float)).all()
                for level in levels_key
            ),
            "criterion": (
                "park_nine_type_joint_sign_of_delta_Vavg_and_delta_lsol_"
                "types_1_3_amplification_ranked_by_ascending_lsol"
            ),
            "criterion_source": "Park et al. 2012",
            "robust_flag_source": "CMM (not in Park et al. 2012)",
            "coupling_constraint_source": (
                "CMM caller-supplied linear equalities (not Park et al.'s "
                "STRING-derived grouping-reaction constraints)"
            ),
        },
    )


def _rises_monotonically(values: np.ndarray, tol: float) -> bool:
    if len(values) < 2 or not np.isfinite(values).all():
        return False
    return bool(values[-1] - values[0] > tol and np.all(np.diff(values) >= -tol))


def _trend_sign(values: np.ndarray, tol: float) -> int:
    """+1 if the series rises across the scan, -1 if it falls, 0 if it does neither.

    A rise requires *both* an endpoint gain and a positive linear-regression slope, so a
    series that ends higher only because of a single noisy step is reported flat.
    """

    if len(values) < 2 or not np.isfinite(values).all():
        return 0
    start, end = float(values[0]), float(values[-1])
    slope = float(np.polyfit(np.arange(len(values), dtype=float), values, 1)[0])
    if end - start > tol and slope > tol:
        return 1
    if start - end > tol and slope < -tol:
        return -1
    return 0


def _reverses_direction(values: np.ndarray, tol: float) -> bool:
    """Whether the signed series changes sign anywhere it is non-negligible."""

    nonzero = values[np.abs(values) > tol]
    return bool(len(nonzero) > 1 and np.any(np.sign(nonzero) != np.sign(nonzero[0])))


def _classify_trend(values: np.ndarray, tol: float) -> str:
    """Classify a direction-preserving flux-magnitude trend across a scan."""

    if len(values) < 2 or not np.isfinite(values).all():
        return "none"
    # A direction reversal is not an amplification/knockdown target in FSEOF/FVSEOF.
    if _reverses_direction(values, tol):
        return "none"
    sign = _trend_sign(np.abs(values), tol)
    if sign > 0:
        return "amplify"
    if sign < 0:
        return "knockdown"
    return "none"


#: Park et al. (2012) classify reaction responses into "nine types based on combinations of
#: positive and negative changes in V_avg and l_sol", with types 1-3 positively correlated
#: with production (the amplification candidates), 4-6 negatively correlated, and 7-9
#: showing no clear correlation. The **band** (1-3 / 4-6 / 7-9) is therefore Park's and is
#: fixed by the sign of the change in V_avg. The **order within a band** is CMM's
#: convention: the paper's own type table lives in a figure that could not be resolved from
#: the text, so rather than invent their numbering CMM orders each band by the sign of the
#: change in l_sol, narrowing range first, because Park give "reactions with smaller values
#: of l_sol received higher priorities" as the tie-break among amplification candidates.
_PARK_TYPE_BY_SIGNS: dict[tuple[int, int], int] = {
    (1, -1): 1,
    (1, 0): 2,
    (1, 1): 3,
    (-1, -1): 4,
    (-1, 0): 5,
    (-1, 1): 6,
    (0, -1): 7,
    (0, 0): 8,
    (0, 1): 9,
}

_PARK_TYPE_CLASSIFICATION: dict[int, str] = {
    **{t: "amplify" for t in (1, 2, 3)},
    **{t: "knockdown" for t in (4, 5, 6)},
    **{t: "none" for t in (7, 8, 9)},
}


def _park_flux_type(vavg: np.ndarray, lsol: np.ndarray, tol: float) -> int:
    """Park's nine-type index from the joint sign of delta V_avg and delta l_sol.

    V_avg is read as a magnitude, consistent with CMM's convention everywhere else in the
    scan methods: a reaction running in reverse whose flux grows more negative is carrying
    *more* flux, not less. A reaction that reverses direction mid-scan has no clear
    correlation by construction and falls in the 7-9 band.
    """

    if len(vavg) < 2 or not np.isfinite(vavg).all() or not np.isfinite(lsol).all():
        return 8  # no trend measurable in either statistic
    vavg_sign = 0 if _reverses_direction(vavg, tol) else _trend_sign(np.abs(vavg), tol)
    return _PARK_TYPE_BY_SIGNS[(vavg_sign, _trend_sign(lsol, tol))]


def _objective_reaction(model: Model) -> str:
    for rxn in model.reactions:
        if rxn.objective_coefficient != 0:
            return rxn.id
    raise ValueError("model has no objective reaction")
