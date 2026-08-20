"""Growth-media compositions and presets.

A :class:`Medium` is a set of substrate uptake limits keyed by exchange reaction. Applying a
medium closes every other uptake and opens the listed ones (via cobra's ``model.medium`` API),
so swapping media is a single call. Exchange ids are resolved tolerantly across naming
conventions (BiGG ``EX_glc__D_e`` vs ``EX_glc_e`` vs ``EX_glc(e)``) so a preset works on many
models.

**Where the preset compositions come from.** The shared mineral/trace-element set is the
default medium distributed with the *iJO1366* reconstruction (Orth et al. 2011,
doi:10.1038/msb.2011.65) as shipped by BiGG/COBRApy (King et al. 2016,
doi:10.1093/nar/gkv1049). Applied to iJO1366, ``glucose_aerobic`` reproduces that model's
shipped default medium exactly - all 24 components matching by id and by value, plus vitamin
B12 (``EX_cbl1_e`` 0.01) which makes the match complete. **No publication defines this
composition:** Orth et al. 2011 states glucose -8 and O2 -18.5, calls those "arbitrary bounds
based on typical substrate uptake rates", and never lists the mineral exchanges at all. The
attribution is therefore to the model artefact, not to a composition table in a paper.

**These media are not M9.** The M9 laboratory recipe (Cold Spring Harbor Protocols
pdb.rec12295 + pdb.rec614) contains no trace metals and no iron; a medium restricted to the
literal M9 element set gives growth 0.0 on iJO1366, because Fe2+/Fe3+ are essential as a pair
and Mn, Zn, Cu, Co, Mo and Ni are each individually essential. The preset is *richer* than M9,
not a simplification of it. The acetate and glycerol conditions substitute the carbon source
on that same mineral base and are CMM's own choice, not published media.

**CO2 uptake is closed by default, and that is a deliberate deviation.** Both ``e_coli_core``
and iJO1366 ship with ``EX_co2_e`` open for uptake at 1000, and no published convention was
found for closing it. It is closed here because methods that maximise a *product* rather than
biomass actually exercise that free CO2 inlet: leaving it open inflates succinate yield on
``e_coli_core`` from 1.2000 to 1.3906 mol/mol below a carbon ceiling no guard then trips.
Measured cost of closing it: 0.000% of aerobic growth, at most 0.25% of anaerobic growth, and
it removes 8 of 8 measured carbon-ceiling violations. Pass ``co2_uptake=True`` for a genuinely
CO2-supplemented condition. Secretion is never restricted either way.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
import math
import warnings

from cobra import Model
import pandas as pd

#: Uptake cap for a component that is present but not intended to limit growth.
UNLIMITED_UPTAKE = 1000.0

# Mineral-salt and trace-element exchanges of the iJO1366 default medium as distributed by
# BiGG/COBRApy, each opened at UNLIMITED_UPTAKE (any naming variant present is used).
# EX_co2_e is deliberately NOT in this set - see the module docstring.
#
#   EX_nh4_e      ammonium, the nitrogen source
#   EX_pi_e       inorganic phosphate
#   EX_so4_e      sulfate, the sulfur source
#   EX_h2o_e      water
#   EX_h_e        protons (charge/pH balance)
#   EX_k_e        potassium        EX_na1_e   sodium         EX_cl_e    chloride
#   EX_ca2_e      calcium          EX_mg2_e   magnesium
#   EX_mn2_e      manganese        EX_fe2_e   ferrous iron   EX_fe3_e   ferric iron
#   EX_zn2_e      zinc             EX_cu2_e   copper         EX_cobalt2_e cobalt
#   EX_mobd_e     molybdate        EX_ni2_e   nickel
#   EX_slnt_e     selenite         EX_sel_e   selenate       EX_tungs_e tungstate
_BIGG_MINERAL_IONS: tuple[str, ...] = (
    "EX_nh4_e",
    "EX_pi_e",
    "EX_so4_e",
    "EX_h2o_e",
    "EX_h_e",
    "EX_k_e",
    "EX_na1_e",
    "EX_cl_e",
    "EX_ca2_e",
    "EX_mg2_e",
    "EX_mn2_e",
    "EX_fe2_e",
    "EX_fe3_e",
    "EX_zn2_e",
    "EX_cu2_e",
    "EX_cobalt2_e",
    "EX_mobd_e",
    "EX_ni2_e",
    "EX_slnt_e",
    "EX_sel_e",
    "EX_tungs_e",
)

#: Vitamin B12, the one component of iJO1366's shipped default medium that is capped below
#: UNLIMITED_UPTAKE. Including it makes ``glucose_aerobic`` an exact reproduction of that
#: medium; it is absent from small models such as ``e_coli_core``, where it drops harmlessly
#: (measured: growth identical with and without).
_VITAMIN_B12 = ("EX_cbl1_e", 0.01)

# Naming-convention variants for tolerant resolution.
_ALIASES: dict[str, tuple[str, ...]] = {
    "EX_glc__D_e": ("EX_glc__D_e", "EX_glc_D_e", "EX_glc_e", "EX_glc(e)"),
    "EX_o2_e": ("EX_o2_e", "EX_o2(e)"),
    "EX_ac_e": ("EX_ac_e", "EX_ac(e)"),
    "EX_succ_e": ("EX_succ_e", "EX_succ(e)"),
    "EX_glyc_e": ("EX_glyc_e", "EX_glyc(e)"),
    "EX_co2_e": ("EX_co2_e", "EX_co2(e)"),
}


def _resolve_id(model: Model, exchange_id: str) -> str | None:
    """Find the exchange reaction in the model for an id, tolerating naming variants."""

    present = model.reactions
    for candidate in _ALIASES.get(exchange_id, (exchange_id,)):
        if candidate in present:
            return candidate
    return None


def _mineral_base(*, co2_uptake: bool) -> dict[str, float]:
    """The shared mineral/trace-element base, with CO2 uptake opened only on request."""

    rates: dict[str, float] = {ion: UNLIMITED_UPTAKE for ion in _BIGG_MINERAL_IONS}
    rates[_VITAMIN_B12[0]] = _VITAMIN_B12[1]
    if co2_uptake:
        rates["EX_co2_e"] = UNLIMITED_UPTAKE
    return rates


@dataclass(frozen=True)
class MediumApplication(Mapping[str, float]):
    """What a :meth:`Medium.apply_to` call actually did to a model.

    Behaves as the mapping of applied ``{exchange_id: uptake}`` it replaces, so callers that
    only want the applied bounds keep working, while ``dropped`` and :meth:`to_provenance`
    record the components the model could not express. Recording both is the point: a run
    directory that stores only the preset key does not state the condition it ran under.
    """

    medium: str
    applied: Mapping[str, float] = field(default_factory=dict)
    dropped: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "applied", dict(self.applied))
        object.__setattr__(self, "dropped", tuple(self.dropped))

    def __getitem__(self, exchange_id: str) -> float:
        return self.applied[exchange_id]

    def __iter__(self) -> Iterator[str]:
        return iter(self.applied)

    def __len__(self) -> int:
        return len(self.applied)

    def to_provenance(self) -> dict[str, object]:
        """Provenance record of the applied condition, for ``run_provenance`` parameters."""

        return {
            "medium": self.medium,
            "applied": dict(self.applied),
            "dropped": list(self.dropped),
        }

    def to_frame(self) -> pd.DataFrame:
        """Deterministic export table: one row per declared component, applied or not."""

        records = [
            {"exchange_id": exchange_id, "uptake": rate, "status": "applied"}
            for exchange_id, rate in sorted(self.applied.items())
        ]
        records.extend(
            {"exchange_id": exchange_id, "uptake": float("nan"), "status": "dropped"}
            for exchange_id in sorted(self.dropped)
        )
        return pd.DataFrame(records, columns=["exchange_id", "uptake", "status"])


@dataclass(frozen=True)
class Medium:
    """A named medium: maximum uptake (a positive flux) per substrate exchange.

    ``required`` names the components whose absence would make this a *different*
    experiment rather than a degraded one - the carbon source, above all. Dropping one of
    them raises; dropping anything else warns and is recorded. Left as ``None`` it is
    inferred as every component capped below :data:`UNLIMITED_UPTAKE` and above zero, i.e.
    everything the medium declares as growth-limiting.
    """

    name: str
    uptake: Mapping[str, float] = field(default_factory=dict)
    required: frozenset[str] | None = None

    def __post_init__(self) -> None:
        if self.required is None:
            object.__setattr__(self, "required", _infer_required(self.uptake))
        else:
            object.__setattr__(self, "required", frozenset(self.required))

    def resolve(self, model: Model) -> dict[str, float]:
        """Map this medium onto a model's exchange ids, dropping absent components.

        A plain lookup with no policy: use :meth:`apply_to` for the tiered handling of
        components the model cannot express.
        """

        applied, _ = self._partition(model)
        return applied

    def _partition(self, model: Model) -> tuple[dict[str, float], tuple[str, ...]]:
        """Split declared components into those the model has and those it lacks."""

        applied: dict[str, float] = {}
        dropped: list[str] = []
        for exchange_id, rate in self.uptake.items():
            rate = float(rate)
            if not math.isfinite(rate) or rate < 0:
                raise ValueError(
                    f"uptake rate for {exchange_id!r} must be finite and non-negative"
                )
            actual = _resolve_id(model, exchange_id)
            if actual is None:
                dropped.append(exchange_id)
            else:
                applied[actual] = rate
        return applied, tuple(dropped)

    def apply_to(self, model: Model) -> MediumApplication:
        """Set the model's medium to this composition; report what was applied and dropped.

        Raises when a **required** component has no exchange in the model: a "glycerol
        minimal" medium on a model without a glycerol exchange is a carbon-free, infeasible
        model, and the resulting solver error names neither the medium nor the substrate.
        Missing minerals only warn, because 17 of 24 components legitimately drop on
        ``e_coli_core`` and the applied medium is then bit-identical to that model's shipped
        default.
        """

        applied, dropped = self._partition(model)
        required = self.required or frozenset()
        missing_required = tuple(sorted(set(dropped) & required))
        if missing_required:
            tried = {
                exchange_id: _ALIASES.get(exchange_id, (exchange_id,))
                for exchange_id in missing_required
            }
            raise ValueError(
                f"medium {self.name!r} requires {list(missing_required)}, which model "
                f"{model.id!r} has no exchange reaction for (tried {tried}). Applying it "
                "would run a different experiment - most often a model with no carbon "
                "source at all. Choose a medium this model can express, or add the "
                "exchange."
            )
        if self.uptake and not applied:
            raise ValueError(
                f"medium {self.name!r} has no exchange reactions present in model {model.id!r}"
            )
        if dropped:
            warnings.warn(
                f"medium {self.name!r}: {len(dropped)} declared component(s) have no "
                f"exchange in model {model.id!r} and were not applied: "
                f"{sorted(dropped)}. They are recorded in the run provenance.",
                UserWarning,
                stacklevel=2,
            )
        model.medium = applied
        return MediumApplication(medium=self.name, applied=applied, dropped=dropped)


def _infer_required(uptake: Mapping[str, float]) -> frozenset[str]:
    """Components declared as growth-limiting: capped above 0 and below unlimited."""

    required: set[str] = set()
    for exchange_id, rate in uptake.items():
        try:
            value = float(rate)
        except (TypeError, ValueError):
            continue
        if 0.0 < value < UNLIMITED_UPTAKE:
            required.add(exchange_id)
    return frozenset(required)


def glucose_minimal(
    *, aerobic: bool = True, glucose: float = 10.0, co2_uptake: bool = False
) -> Medium:
    """Glucose minimal medium on the iJO1366/BiGG mineral set (aerobic or anaerobic).

    Composition, in mmol gDW-1 h-1 as ``model.medium`` uptake caps:

    ==================  =========  ===============================================
    Component           Cap        Role
    ==================  =========  ===============================================
    ``EX_glc__D_e``     ``glucose``  carbon and energy source; the limiting substrate
    ``EX_o2_e``         1000 / 0   electron acceptor; explicitly 0 when anaerobic
    ``EX_cbl1_e``       0.01       vitamin B12
    21 mineral ions     1000       :data:`_BIGG_MINERAL_IONS`, listed at its definition
    ``EX_co2_e``        closed     opened at 1000 only when ``co2_uptake=True``
    ==================  =========  ===============================================

    Identical to the default medium distributed with iJO1366 apart from the closed CO2
    uptake; **not** the M9 laboratory recipe. See the module docstring for the attribution,
    the M9 distinction and the CO2 measurement.
    """

    uptake = _mineral_base(co2_uptake=co2_uptake)
    uptake["EX_glc__D_e"] = glucose
    # Written explicitly rather than omitted: cobra closes an unlisted exchange to uptake
    # either way, but an explicit zero makes the anaerobic condition auditable.
    uptake["EX_o2_e"] = UNLIMITED_UPTAKE if aerobic else 0.0
    return Medium(
        name=f"Glucose minimal, {_describe(aerobic=aerobic, co2_uptake=co2_uptake)}",
        uptake=uptake,
        required=frozenset({"EX_glc__D_e"}),
    )


def carbon_minimal(
    carbon_exchange: str,
    *,
    aerobic: bool = True,
    uptake: float = 10.0,
    co2_uptake: bool = False,
    label: str | None = None,
) -> Medium:
    """Minimal medium on an arbitrary carbon source over the iJO1366/BiGG mineral set.

    Same composition as :func:`glucose_minimal` with ``carbon_exchange`` substituted for
    glucose. The mineral base is sourced (module docstring); **the choice of carbon source
    and its 10 mmol gDW-1 h-1 cap are CMM's own**, carried over from the glucose default,
    not a published medium and not a measured uptake rate. ``label`` names the substrate in
    the medium's display name; it defaults to the exchange id's stem.
    """

    rates = _mineral_base(co2_uptake=co2_uptake)
    rates[carbon_exchange] = uptake
    rates["EX_o2_e"] = UNLIMITED_UPTAKE if aerobic else 0.0
    if label is None:
        label = carbon_exchange.replace("EX_", "").replace("_e", "")
    return Medium(
        name=(
            f"{label} minimal, {_describe(aerobic=aerobic, co2_uptake=co2_uptake)}, "
            "CMM-defined carbon source"
        ),
        uptake=rates,
        required=frozenset({carbon_exchange}),
    )


def _describe(*, aerobic: bool, co2_uptake: bool) -> str:
    """Condition suffix for a medium name: aeration, mineral-set provenance, CO2 policy."""

    aeration = "aerobic" if aerobic else "anaerobic"
    co2 = "CO2 uptake open" if co2_uptake else "CO2 uptake closed"
    return f"{aeration} (iJO1366/BiGG mineral set; {co2})"


PRESET_MEDIA: dict[str, Medium] = {
    "glucose_aerobic": glucose_minimal(aerobic=True),
    "glucose_anaerobic": glucose_minimal(aerobic=False),
    "acetate_aerobic": carbon_minimal("EX_ac_e", aerobic=True, label="Acetate"),
    "glycerol_aerobic": carbon_minimal("EX_glyc_e", aerobic=True, label="Glycerol"),
}


def preset_medium(name: str) -> Medium:
    """Look up a preset medium by key (see ``PRESET_MEDIA``)."""

    if name not in PRESET_MEDIA:
        raise KeyError(
            f"unknown preset medium {name!r}; available: {sorted(PRESET_MEDIA)}"
        )
    return PRESET_MEDIA[name]


def apply_medium(model: Model, medium: Medium | str) -> MediumApplication:
    """Apply a :class:`Medium` (or preset name) to the model in place."""

    if isinstance(medium, str):
        medium = preset_medium(medium)
    return medium.apply_to(model)
