"""Solver capability detection and tiered gating.

The cobra default solver is LP-only (GLPK). Several services need richer capabilities:
QP for MOMA and revert-metabolism, MILP for OptKnock/RobustKnock, MIQP for original MTA.
Rather than failing deep inside a solve with an opaque optlang error, services call
:func:`require` up front and get a typed, actionable :class:`SolverCapabilityError`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from cobra import Configuration
from cobra.util.solver import solvers as _available_interfaces

Capability = Literal["LP", "QP", "MILP", "MIQP"]

# Capability matrix keyed by cobra/optlang short solver name. Conservative: a name is
# listed with a capability only when that interface genuinely supports it.
_CAPABILITIES: dict[str, frozenset[str]] = {
    "glpk": frozenset({"LP", "MILP"}),
    "glpk_exact": frozenset({"LP"}),
    "scipy": frozenset({"LP"}),
    "osqp": frozenset({"LP", "QP"}),
    "highs": frozenset({"LP", "MILP"}),
    "coinor_cbc": frozenset({"LP", "MILP"}),
    "cplex": frozenset({"LP", "QP", "MILP", "MIQP"}),
    "gurobi": frozenset({"LP", "QP", "MILP", "MIQP"}),
}


class SolverCapabilityError(RuntimeError):
    """Raised when the active solver cannot run a requested optimization class."""

    def __init__(self, capability: str, solver: str, *, feature: str | None = None):
        self.capability = capability
        self.solver = solver
        self.feature = feature
        capable = sorted(s for s, caps in _CAPABILITIES.items() if capability in caps)
        installed = sorted(set(available_solvers()) & set(capable))
        where = f" for {feature}" if feature else ""
        hint = (
            f"installed solvers with {capability}: {installed}"
            if installed
            else f"no installed solver provides {capability}; install one of {capable}"
        )
        super().__init__(
            f"active solver {solver!r} does not support {capability}{where}. {hint}."
        )


def _short_name(solver: object) -> str:
    """Normalize a solver module/name to its short key (e.g. 'gurobi')."""

    if solver is None:
        solver = Configuration().solver
    name = getattr(solver, "__name__", solver)
    name = str(name)
    if name.startswith("optlang."):
        name = name.split(".", 1)[1]
    if name.endswith("_interface"):
        name = name[: -len("_interface")]
    return name


def available_solvers() -> tuple[str, ...]:
    """Short names of solvers cobra can actually construct in this environment."""

    return tuple(sorted(_available_interfaces))


def active_solver(model=None) -> str:
    """Short name of the solver a model (or the global default) would use."""

    if model is not None:
        return _short_name(model.solver.interface)
    return _short_name(None)


def capabilities(solver: object | None = None) -> frozenset[str]:
    """Capability set of the given (or active) solver."""

    return _CAPABILITIES.get(_short_name(solver), frozenset({"LP"}))


def supports(capability: Capability, solver: object | None = None) -> bool:
    """Whether the given (or active) solver supports a capability."""

    return capability in capabilities(solver)


def require(capability: Capability, solver: object | None = None, *, feature: str | None = None) -> None:
    """Raise :class:`SolverCapabilityError` if the capability is unavailable."""

    if not supports(capability, solver):
        raise SolverCapabilityError(capability, _short_name(solver), feature=feature)


@dataclass(frozen=True)
class SolverStatus:
    """A snapshot of the active solver and what it can run."""

    name: str
    capabilities: tuple[str, ...]
    recommended: bool  # full LP+QP+MILP+MIQP stack (gurobi/cplex)
    available: tuple[str, ...]

    @property
    def warning(self) -> str | None:
        """A human warning when the solver cannot run the full workbench, else None."""

        missing = sorted({"QP", "MILP", "MIQP"} - set(self.capabilities))
        if not missing:
            return None
        return (
            f"Solver '{self.name}' lacks {', '.join(missing)}; "
            "MOMA/ROOM/revert (QP), OptKnock (MILP), and original-MTA (MIQP) need a "
            "commercial solver (Gurobi or CPLEX)."
        )

    def summary(self) -> str:
        tag = "recommended" if self.recommended else "limited"
        return f"{self.name} ({tag}): {', '.join(self.capabilities)}"


def solver_status(model=None) -> SolverStatus:
    """Report the active solver, its capabilities, and whether it is recommended."""

    caps = capabilities(model.solver.interface if model is not None else None)
    recommended = {"LP", "QP", "MILP", "MIQP"} <= set(caps)
    return SolverStatus(
        name=active_solver(model),
        capabilities=tuple(sorted(caps)),
        recommended=recommended,
        available=available_solvers(),
    )
