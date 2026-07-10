"""Growth-media compositions and presets.

A :class:`Medium` is a set of substrate uptake limits keyed by exchange reaction. Applying a
medium closes every other uptake and opens the listed ones (via cobra's ``model.medium`` API),
so swapping media is a single call. Exchange ids are resolved tolerantly across naming
conventions (BiGG ``EX_glc__D_e`` vs ``EX_glc_e`` vs ``EX_glc(e)``) so a preset works on many
models.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
import math

from cobra import Model

# Mineral-salt ions opened in a minimal medium (any naming variant present is used).
_MINIMAL_IONS: tuple[str, ...] = (
    "EX_nh4_e",
    "EX_pi_e",
    "EX_so4_e",
    "EX_co2_e",
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

# Naming-convention variants for tolerant resolution.
_ALIASES: dict[str, tuple[str, ...]] = {
    "EX_glc__D_e": ("EX_glc__D_e", "EX_glc_D_e", "EX_glc_e", "EX_glc(e)"),
    "EX_o2_e": ("EX_o2_e", "EX_o2(e)"),
    "EX_ac_e": ("EX_ac_e", "EX_ac(e)"),
    "EX_succ_e": ("EX_succ_e", "EX_succ(e)"),
    "EX_glyc_e": ("EX_glyc_e", "EX_glyc(e)"),
}


def _resolve_id(model: Model, exchange_id: str) -> str | None:
    """Find the exchange reaction in the model for an id, tolerating naming variants."""

    present = model.reactions
    for candidate in _ALIASES.get(exchange_id, (exchange_id,)):
        if candidate in present:
            return candidate
    return None


@dataclass(frozen=True)
class Medium:
    """A named medium: maximum uptake (a positive flux) per substrate exchange."""

    name: str
    uptake: Mapping[str, float] = field(default_factory=dict)

    def resolve(self, model: Model) -> dict[str, float]:
        """Map this medium onto a model's exchange ids, dropping absent components."""

        resolved: dict[str, float] = {}
        for exchange_id, rate in self.uptake.items():
            rate = float(rate)
            if not math.isfinite(rate) or rate < 0:
                raise ValueError(
                    f"uptake rate for {exchange_id!r} must be finite and non-negative"
                )
            actual = _resolve_id(model, exchange_id)
            if actual is not None:
                resolved[actual] = rate
        return resolved

    def apply_to(self, model: Model) -> dict[str, float]:
        """Set the model's medium to this composition; returns what was applied."""

        resolved = self.resolve(model)
        if self.uptake and not resolved:
            raise ValueError(
                f"medium {self.name!r} has no exchange reactions present in model {model.id!r}"
            )
        model.medium = resolved
        return resolved


def glucose_minimal(*, aerobic: bool = True, glucose: float = 10.0) -> Medium:
    """M9-style glucose minimal medium (aerobic or anaerobic)."""

    uptake: dict[str, float] = {"EX_glc__D_e": glucose}
    for ion in _MINIMAL_IONS:
        uptake[ion] = 1000.0
    if aerobic:
        uptake["EX_o2_e"] = 1000.0
    name = f"Glucose minimal ({'aerobic' if aerobic else 'anaerobic'})"
    return Medium(name=name, uptake=uptake)


def carbon_minimal(
    carbon_exchange: str, *, aerobic: bool = True, uptake: float = 10.0
) -> Medium:
    """Minimal medium on an arbitrary carbon source exchange."""

    rates: dict[str, float] = {carbon_exchange: uptake}
    for ion in _MINIMAL_IONS:
        rates[ion] = 1000.0
    if aerobic:
        rates["EX_o2_e"] = 1000.0
    label = carbon_exchange.replace("EX_", "").replace("_e", "")
    return Medium(
        name=f"{label} minimal ({'aerobic' if aerobic else 'anaerobic'})", uptake=rates
    )


PRESET_MEDIA: dict[str, Medium] = {
    "glucose_aerobic": glucose_minimal(aerobic=True),
    "glucose_anaerobic": glucose_minimal(aerobic=False),
    "acetate_aerobic": carbon_minimal("EX_ac_e", aerobic=True),
    "glycerol_aerobic": carbon_minimal("EX_glyc_e", aerobic=True),
}


def preset_medium(name: str) -> Medium:
    """Look up a preset medium by key (see ``PRESET_MEDIA``)."""

    if name not in PRESET_MEDIA:
        raise KeyError(
            f"unknown preset medium {name!r}; available: {sorted(PRESET_MEDIA)}"
        )
    return PRESET_MEDIA[name]


def apply_medium(model: Model, medium: Medium | str) -> dict[str, float]:
    """Apply a :class:`Medium` (or preset name) to the model in place."""

    if isinstance(medium, str):
        medium = preset_medium(medium)
    return medium.apply_to(model)
