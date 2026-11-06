"""What a ranked knockout changes: the MOMA flux state behind its transformation score.

An MTA/rMTA ranking says *which* knockout moves the source state furthest toward the target. A
researcher choosing one to test also needs *how*: which of the changed reactions it moves the
desired way, which it moves the wrong way, and which reactions and subsystems it blocks. That
is what decides whether a prediction is biologically plausible and what to measure.

The explanation re-solves the MOMA QP that produces ``mTS`` (Valcárcel et al. 2019; the
``moma`` component of :func:`cmm.features.revert.revert_targets`). Two properties make it the
right solution to explain:

* It is **unique**. The objective is a strictly convex sum of squares over every reaction, so
  re-solving returns the same flux vector the ranking scored. The best- and worst-case MTA
  solutions are MIQP optima and need not be unique, so explaining them could describe a
  different optimum from the one that was scored.
* It is what the score **is**. Outside Equation 9's amplification branch ``rTS = mTS``, and
  inside it ``rTS`` is ``mTS`` scaled by ``k·(bTS − wTS)``. The recomputed ``mTS`` is returned
  beside the counts so a reader can check the two agree.

Classifying a move uses the transformation score's own sign rule (any change beyond solver
tolerance), so the counts reconcile with ``mTS``. Whether each move also clears the reaction's
epsilon — MTA's threshold for a change that counts — is recorded separately.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal, cast

import pandas as pd
from cobra import Model

from cmm.core import solvers
from cmm.core.flux_state import FluxState
from cmm.core.provenance import run_provenance
from cmm.features._perturbation import (
    apply_perturbation,
    gene_perturbations,
    reaction_perturbations,
)
from cmm.features.revert import (
    DEFAULT_EPSILON,
    Epsilon,
    _epsilon_for,
    _mta_qp,
    _prepared_direction,
    _sign,
    _steady_floor,
    _target_reactions,
    _transformation_score,
)
from cmm.omics.differential import DirectionMap

MoveOutcome = Literal["desired", "opposite", "unchanged"]

_SUMMARY_COLUMNS = (
    "target_id",
    "status",
    "mTS",
    "n_changed_reactions",
    "n_desired",
    "n_opposite",
    "n_unchanged",
    "n_desired_clearing_epsilon",
    "desired_l1",
    "opposite_l1",
    "steady_l1",
    "n_blocked_reactions",
    "blocked_reactions",
    "blocked_subsystems",
)
_CHANGE_COLUMNS = (
    "target_id",
    "reaction_id",
    "reaction_name",
    "subsystem",
    "desired_direction",
    "v_ref",
    "v_knockout",
    "delta",
    "epsilon",
    "outcome",
    "clears_epsilon",
)


@dataclass(frozen=True)
class KnockoutExplanation:
    """One knockout's MOMA state, reduced to what a reader needs to judge it."""

    target_id: str
    status: str
    mts: float
    blocked_reactions: tuple[str, ...]
    blocked_subsystems: tuple[str, ...]
    #: One row per changed reaction, in the order of the changed set.
    changes: tuple[Mapping[str, object], ...]
    steady_l1: float

    def counts(self) -> dict[str, int]:
        outcomes = [str(row["outcome"]) for row in self.changes]
        return {
            "n_changed_reactions": len(outcomes),
            "n_desired": outcomes.count("desired"),
            "n_opposite": outcomes.count("opposite"),
            "n_unchanged": outcomes.count("unchanged"),
            "n_desired_clearing_epsilon": sum(
                1
                for row in self.changes
                if row["outcome"] == "desired" and row["clears_epsilon"]
            ),
        }

    def summary(self) -> dict[str, object]:
        desired = 0.0
        opposite = 0.0
        for row in self.changes:
            size = abs(float(cast(float, row["delta"])))
            if row["outcome"] == "desired":
                desired += size
            elif row["outcome"] == "opposite":
                opposite += size
        return {
            "target_id": self.target_id,
            "status": self.status,
            "mTS": self.mts,
            **self.counts(),
            "desired_l1": desired,
            "opposite_l1": opposite,
            "steady_l1": self.steady_l1,
            "n_blocked_reactions": len(self.blocked_reactions),
            "blocked_reactions": "; ".join(self.blocked_reactions),
            "blocked_subsystems": "; ".join(self.blocked_subsystems),
        }


@dataclass(frozen=True)
class KnockoutMechanism:
    """Explanations for a set of knockouts, with the provenance of the solves behind them."""

    explanations: tuple[KnockoutExplanation, ...]
    metadata: Mapping[str, object] = field(default_factory=dict)

    def summary_frame(self) -> pd.DataFrame:
        """One row per knockout: move counts, recomputed ``mTS`` and what it blocks."""

        return pd.DataFrame(
            [explanation.summary() for explanation in self.explanations],
            columns=list(_SUMMARY_COLUMNS),
        )

    def to_frame(self) -> pd.DataFrame:
        """One row per knockout and changed reaction: the flux move behind each count."""

        return pd.DataFrame(
            [
                {"target_id": explanation.target_id, **row}
                for explanation in self.explanations
                for row in explanation.changes
            ],
            columns=list(_CHANGE_COLUMNS),
        )


def _subsystem(model: Model, reaction_id: str) -> str:
    return (model.reactions.get_by_id(reaction_id).subsystem or "").strip()


def knockout_mechanism(
    model: Model,
    reference_state: FluxState,
    direction: DirectionMap,
    targets: Iterable[str],
    *,
    perturbation: Literal["gene", "reaction"] = "gene",
    epsilon: Epsilon = DEFAULT_EPSILON,
) -> KnockoutMechanism:
    """Explain each knockout by the MOMA flux state its ``mTS`` was scored on.

    ``reference_state``, ``direction`` and ``epsilon`` must be the ones the ranking used, and
    ``model`` must carry the same medium and condition: the changed set is masked against the
    unperturbed model exactly as :func:`cmm.features.revert.revert_targets` masks it, so the
    recomputed ``mTS`` matches the ranked one. Needs a QP-capable solver. A knockout whose
    solve is not optimal is kept with its status and no reaction rows.
    """

    if perturbation not in {"gene", "reaction"}:
        raise ValueError("perturbation must be 'gene' or 'reaction'")
    requested = list(dict.fromkeys(targets))
    solvers.require("QP", model.solver.interface, feature="knockout mechanism (MOMA)")
    reference_state.validate()
    provenance = run_provenance(
        model,
        method="knockout_mechanism_moma",
        perturbation=perturbation,
        n_targets=len(requested),
    )

    explanations: list[KnockoutExplanation] = []
    with model:
        target_rxns = _target_reactions(model, reference_state)
        prepared = _prepared_direction(model, reference_state, direction, target_rxns)
        changed = [rid for rid in target_rxns if prepared.get(rid, 0) != 0]
        floor = _steady_floor(epsilon)
        lookup = model.genes if perturbation == "gene" else model.reactions
        known = [target_id for target_id in requested if target_id in lookup]
        perts = (
            gene_perturbations(model, known)
            if perturbation == "gene"
            else reaction_perturbations(model, known)
        )
        by_id = {pert.target_id: pert for pert in perts}
        for target_id in requested:
            pert = by_id.get(target_id)
            if pert is None:
                explanations.append(
                    KnockoutExplanation(
                        target_id=target_id,
                        status="not_a_knockout",
                        mts=float("nan"),
                        blocked_reactions=(),
                        blocked_subsystems=(),
                        changes=(),
                        steady_l1=float("nan"),
                    )
                )
                continue
            blocked = tuple(sorted(pert.reaction_ids))
            blocked_subsystems = tuple(
                sorted({s for s in (_subsystem(model, rid) for rid in blocked) if s})
            )
            with apply_perturbation(model, pert):
                fluxes, status = _mta_qp(
                    model, reference_state, prepared, target_rxns, 0.0
                )
            if status != "optimal":
                explanations.append(
                    KnockoutExplanation(
                        target_id=target_id,
                        status=status,
                        mts=float("nan"),
                        blocked_reactions=blocked,
                        blocked_subsystems=blocked_subsystems,
                        changes=(),
                        steady_l1=float("nan"),
                    )
                )
                continue
            rows: list[Mapping[str, object]] = []
            for rid in changed:
                desired = prepared.get(rid, 0)
                v_ref = reference_state.get(rid)
                v_ko = fluxes.get(rid, 0.0)
                delta = v_ko - v_ref
                moved = _sign(delta)
                outcome: MoveOutcome = (
                    "desired"
                    if moved == desired
                    else "opposite"
                    if moved == -desired
                    else "unchanged"
                )
                threshold = _epsilon_for(epsilon, rid)
                reaction = model.reactions.get_by_id(rid)
                rows.append(
                    {
                        "reaction_id": rid,
                        "reaction_name": reaction.name or "",
                        "subsystem": (reaction.subsystem or "").strip(),
                        "desired_direction": desired,
                        "v_ref": v_ref,
                        "v_knockout": v_ko,
                        "delta": delta,
                        "epsilon": threshold,
                        "outcome": outcome,
                        "clears_epsilon": desired * delta >= threshold,
                    }
                )
            steady = sum(
                abs(fluxes.get(rid, 0.0) - reference_state.get(rid))
                for rid in target_rxns
                if prepared.get(rid, 0) == 0
            )
            explanations.append(
                KnockoutExplanation(
                    target_id=target_id,
                    status=status,
                    mts=_transformation_score(
                        fluxes,
                        reference_state,
                        prepared,
                        target_rxns,
                        steady_floor=floor,
                    ),
                    blocked_reactions=blocked,
                    blocked_subsystems=blocked_subsystems,
                    changes=tuple(rows),
                    steady_l1=steady,
                )
            )

    return KnockoutMechanism(
        explanations=tuple(explanations),
        metadata={
            **provenance,
            "solution": "moma_qp",
            "outcome_rule": "sign of v_knockout - v_ref against the desired direction "
            "(the transformation score's rule); clears_epsilon applies MTA's threshold",
            "n_changed_reactions": len(changed),
            "n_impossible_masked": prepared.metadata.get("n_impossible_masked"),
        },
    )


__all__: Sequence[str] = (
    "KnockoutExplanation",
    "KnockoutMechanism",
    "knockout_mechanism",
)
