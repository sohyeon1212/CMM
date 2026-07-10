"""Deterministic model and software provenance for publication-facing results."""

from __future__ import annotations

import hashlib
import json
import platform
from importlib.metadata import PackageNotFoundError, version
from typing import Any

from cobra import Model

from cmm.core.solvers import active_solver


def _package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not-installed"


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


def run_provenance(model: Model, **parameters: object) -> dict[str, object]:
    """Stable provenance metadata suitable for result exports and manuscript supplements."""

    return {
        "model_id": str(model.id),
        "model_sha256": model_fingerprint(model),
        "solver": active_solver(model),
        "python": platform.python_version(),
        "cmm": _package_version("cmm"),
        "cobra": _package_version("cobra"),
        "numpy": _package_version("numpy"),
        "pandas": _package_version("pandas"),
        "scipy": _package_version("scipy"),
        "parameters": dict(parameters),
    }
