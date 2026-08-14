"""Deterministic model and software provenance for publication-facing results."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import platform
from importlib.metadata import PackageNotFoundError, version
from typing import Any

from cobra import Model

from cmm.core.solvers import active_solver

# Distribution that ships each optlang solver interface, so a record can carry the
# solver's own version and not merely its name.
_SOLVER_DISTRIBUTIONS: dict[str, str] = {
    "coinor_cbc": "mip",
    "cplex": "cplex",
    "glpk": "swiglpk",
    "glpk_exact": "swiglpk",
    "gurobi": "gurobipy",
    "highs": "highspy",
    "osqp": "osqp",
    "scipy": "scipy",
}


def _package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not-installed"


def _solver_version(solver: str) -> str:
    """Version string of the backend behind a cobra/optlang solver short name."""

    distribution = _SOLVER_DISTRIBUTIONS.get(solver)
    if distribution is None:
        return "unknown"
    return _package_version(distribution)


def _utc_timestamp() -> str:
    """Current time as an ISO-8601 UTC instant, e.g. ``2026-08-10T09:31:07.512873Z``."""

    stamp = datetime.now(timezone.utc).isoformat(timespec="microseconds")
    return stamp.replace("+00:00", "Z")


def model_fingerprint(model: Model) -> str:
    """SHA-256 of the model structure, bounds, GPRs, and objective in stable order."""

    reactions: list[dict[str, Any]] = []
    for reaction in sorted(model.reactions, key=lambda item: item.id):
        reactions.append(
            {
                "id": reaction.id,
                "bounds": [float(reaction.lower_bound), float(reaction.upper_bound)],
                "objective": float(reaction.objective_coefficient),
                "gpr": str(reaction.gene_reaction_rule),
                "metabolites": [
                    [metabolite.id, float(coefficient)]
                    for metabolite, coefficient in sorted(
                        reaction.metabolites.items(), key=lambda item: item[0].id
                    )
                ],
            }
        )
    payload = {
        "model_id": str(model.id),
        "objective_direction": str(model.objective_direction),
        "reactions": reactions,
    }
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def run_provenance(
    model: Model, *, seed: int | None = None, **parameters: object
) -> dict[str, object]:
    """Stable provenance metadata suitable for result exports and manuscript supplements.

    Carries everything a reader needs to reproduce the number: the model fingerprint, the
    UTC run time, the random seed the calling method used, the solver *and its version*,
    the platform/OS/CPU, the Python and package versions, and every call parameter.

    ``seed`` is threaded through from the calling method and is never invented here: a
    deterministic method that takes no seed records ``None`` (JSON ``null``), which states
    "this run had no seed" rather than implying an unrecorded one. It is also kept in
    ``parameters`` so records written before this field existed keep their shape.
    """

    solver = active_solver(model)
    recorded_parameters = dict(parameters)
    if seed is not None:
        recorded_parameters["seed"] = seed
    return {
        "model_id": str(model.id),
        "model_sha256": model_fingerprint(model),
        "timestamp_utc": _utc_timestamp(),
        "seed": seed,
        "solver": solver,
        "solver_version": _solver_version(solver),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor() or "unknown",
        "python": platform.python_version(),
        "cmm": _package_version("cmm"),
        "cobra": _package_version("cobra"),
        "numpy": _package_version("numpy"),
        "pandas": _package_version("pandas"),
        "scipy": _package_version("scipy"),
        "parameters": recorded_parameters,
    }
