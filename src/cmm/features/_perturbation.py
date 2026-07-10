"""Knockout perturbation enumeration and application.

Shared by batch MOMA/ROOM, gene-essentiality, and revert-metabolism: each of those needs
to enumerate gene/reaction knockouts, apply one to a model transiently, and run a solve.
A gene knockout is represented purely by the set of reactions it blocks, so application is
always "force these reactions to zero" regardless of GPR complexity.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Literal

from cobra import Model

try:  # cobra >= 0.29 renamed the helper and made it context-revertible
    from cobra.manipulation.delete import knock_out_model_genes as _knock_out_genes

    def _blocked_reaction_ids_for(
        model: Model, gene_ids: Sequence[str]
    ) -> tuple[str, ...]:
        with model:  # bounds are reverted on exit; we only read the result
            blocked = _knock_out_genes(model, list(gene_ids))
            return tuple(sorted(r.id for r in blocked))
except ImportError:  # pragma: no cover - older cobra layout
    from cobra.manipulation.delete import find_gene_knockout_reactions as _find_ko

    def _blocked_reaction_ids_for(
        model: Model, gene_ids: Sequence[str]
    ) -> tuple[str, ...]:
        genes = [model.genes.get_by_id(g) for g in gene_ids]
        blocked = _find_ko(model, genes)
        return tuple(sorted(r.id for r in blocked))


def _blocked_reaction_ids(model: Model, gene_id: str) -> tuple[str, ...]:
    return _blocked_reaction_ids_for(model, [gene_id])


def blocked_reactions_for_genes(
    model: Model, gene_ids: Iterable[str]
) -> tuple[str, ...]:
    """Reactions disabled when *all* the given genes are knocked out together (GPR-aware).

    Resolving the genes jointly (not one at a time) is what makes a multi-gene knockout
    correct: a reaction whose GPR needs two of the knocked-out genes via an ``and`` is only
    blocked when both are gone, which per-gene resolution would miss.
    """

    gene_list = [str(g) for g in gene_ids]
    if not gene_list:
        return ()
    return _blocked_reaction_ids_for(model, gene_list)


PerturbationKind = Literal["gene", "reaction"]


@dataclass(frozen=True)
class Perturbation:
    """A single knockout, expressed as the reactions it forces to zero."""

    target_id: str
    kind: PerturbationKind
    reaction_ids: tuple[str, ...]

    @property
    def is_inert(self) -> bool:
        return len(self.reaction_ids) == 0


def gene_perturbations(
    model: Model,
    genes: Iterable[str] | None = None,
    *,
    include_inert: bool = False,
) -> list[Perturbation]:
    """Enumerate gene knockouts, resolving each gene to the reactions it disables via GPR."""

    gene_objs = (
        [model.genes.get_by_id(g) for g in genes]
        if genes is not None
        else list(model.genes)
    )
    perts: list[Perturbation] = []
    for gene in gene_objs:
        reaction_ids = _blocked_reaction_ids(model, gene.id)
        pert = Perturbation(target_id=gene.id, kind="gene", reaction_ids=reaction_ids)
        if pert.is_inert and not include_inert:
            continue
        perts.append(pert)
    return perts


def grouped_gene_perturbations(
    model: Model,
    separator: str,
    genes: Iterable[str] | None = None,
    *,
    include_inert: bool = False,
) -> list[Perturbation]:
    """Group transcript-like gene IDs before knockout (COBRA rMTA ``SeparateTranscript``).

    For example, separator ``'.'`` groups ``g2.1`` and ``g2.2`` into a joint knockout named
    ``g2``. This is required to reproduce published MTA/rMTA gene-level screens.
    """

    if not separator:
        raise ValueError("separator must not be empty")
    grouped: dict[str, list[str]] = {}
    for gene in model.genes:
        base = gene.id.split(separator, 1)[0]
        grouped.setdefault(base, []).append(gene.id)
    selected = set(str(gene) for gene in genes) if genes is not None else None
    if selected is not None:
        unknown = selected - set(grouped)
        if unknown:
            raise KeyError(f"unknown grouped genes: {sorted(unknown)}")

    perturbations: list[Perturbation] = []
    for base in sorted(grouped):
        if selected is not None and base not in selected:
            continue
        reaction_ids = blocked_reactions_for_genes(model, grouped[base])
        perturbation = Perturbation(base, "gene", reaction_ids)
        if perturbation.is_inert and not include_inert:
            continue
        perturbations.append(perturbation)
    return perturbations


def reaction_perturbations(
    model: Model,
    reactions: Iterable[str] | None = None,
) -> list[Perturbation]:
    """Enumerate single-reaction knockouts."""

    reaction_ids = (
        list(reactions) if reactions is not None else [r.id for r in model.reactions]
    )
    return [
        Perturbation(target_id=rid, kind="reaction", reaction_ids=(rid,))
        for rid in reaction_ids
    ]


@contextmanager
def apply_perturbation(model: Model, pert: Perturbation):
    """Force the perturbation's reactions to zero within a reverting model context."""

    with model:
        for rid in pert.reaction_ids:
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        yield model


def run_perturbations(
    model: Model,
    perturbations: Sequence[Perturbation],
    fn: Callable[[Model, Perturbation], object],
) -> list[tuple[Perturbation, object]]:
    """Apply each perturbation transiently and collect ``fn(model, perturbation)`` results.

    The model is restored to its original bounds after every iteration, so the caller's
    model is unchanged on return.
    """

    results: list[tuple[Perturbation, object]] = []
    for pert in perturbations:
        with apply_perturbation(model, pert):
            results.append((pert, fn(model, pert)))
    return results
