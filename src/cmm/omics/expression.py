"""Single-state omics integration: map gene expression to fluxes (E-Flux2, LAD).

These methods predict an intracellular flux distribution from one transcriptomic sample by
turning gene expression into reaction-level weights through the GPR, then either scaling
reaction bounds by expression (E-Flux2, Kim et al. 2016) or fitting fluxes to
expression-derived targets (LAD, Lee et al. 2012). They complement the two-state
``differential`` module used by revert-metabolism.

References
----------
Kim MK, Lane A, Kelley JJ, Lun DS (2016) PLoS ONE 11(6):e0157101 — E-Flux2.
Lee D, Smallbone K, Dunn WB, et al. (2012) BMC Syst Biol 6:73 — LAD ("Lee-12").
"""

from __future__ import annotations

import ast
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

import numpy as np
from cobra import Model

from cmm.core import solvers
from cmm.core.flux_state import FluxState
from cmm.core.provenance import run_provenance


@dataclass(frozen=True)
class OmicsFluxResult:
    """Flux distribution predicted from expression data."""

    method: str
    status: str
    objective_value: float | None
    fluxes: dict[str, float] = field(default_factory=dict)
    detail: str = ""
    metadata: dict[str, object] = field(default_factory=dict)

    def to_flux_state(self, name: str = "omics") -> FluxState:
        return FluxState(self.fluxes, name=name, provenance="imported")


OR_RULES = ("sum", "max")


def _eval_gpr_value(
    node: ast.AST | None, expr: Mapping[str, float], *, or_rule: str = "sum"
) -> float | None:
    """Continuous GPR evaluation: AND (complex) -> min, OR (isozymes) -> ``or_rule``.

    ``or_rule="sum"`` is the published convention for both consumers of this helper:
    Kim et al. (2016) *"we took the sum of the expression values of the associated genes …
    since the total capacity of the reaction is given by the sum of the capacities of its
    isozymes"*, and Lee et al. (2012) *"the total capacity is given by the sum of its
    components"*. ``or_rule="max"`` is CMM's pre-0.4.0 behaviour and matches no source paper;
    it is retained only so an older result can be reproduced.

    A missing gene contributes ``None``; an AND with any missing subunit is ``None``, which
    makes ``gene_to_reaction_weights`` omit the reaction so its template bounds are left
    untouched — Kim et al.'s missing-data rule (*"kept as they were defined in the template
    model … so as not to constrain the model unnecessarily"*). An OR keeps the available
    alternatives.
    """

    if node is None:
        return None
    if isinstance(node, ast.Module):
        body = getattr(node, "body", None)
        if isinstance(body, list):
            body = body[0] if body else None
        return _eval_gpr_value(body, expr, or_rule=or_rule)
    if isinstance(node, ast.Expression):
        return _eval_gpr_value(node.body, expr, or_rule=or_rule)
    if isinstance(node, ast.Name):
        return expr.get(node.id)
    if isinstance(node, ast.BoolOp):
        values = [_eval_gpr_value(v, expr, or_rule=or_rule) for v in node.values]
        if isinstance(node.op, ast.And):
            if any(value is None for value in values):
                return None
            return min(value for value in values if value is not None)
        available = [v for v in values if v is not None]
        if not available:
            return None
        return sum(available) if or_rule == "sum" else max(available)
    return None


def gene_to_reaction_weights(
    model: Model,
    gene_expression: Mapping[str, float],
    *,
    or_rule: str = "sum",
) -> dict[str, float]:
    """Map gene expression to reaction weights through each reaction's GPR.

    AND (enzyme complex) takes the minimum over subunits; OR (isozymes) takes ``or_rule``,
    ``"sum"`` by default, which is what Kim et al. (2016) specify for E-Flux2 and Lee et al.
    (2012) for LAD (their worked example, ``r_0250 = 0.0976 + 0.126 = 0.224``, is a sum).
    Pass ``or_rule="max"`` only to reproduce a pre-0.4.0 CMM result; no source paper uses it.
    """

    if or_rule not in OR_RULES:
        raise ValueError(f"or_rule must be one of {OR_RULES}; got {or_rule!r}")

    lookup: dict[str, float] = {}
    for gene in model.genes:
        for key in (gene.id, gene.id.upper(), gene.id.lower()):
            if key in gene_expression:
                value = float(gene_expression[key])
                if not np.isfinite(value) or value < 0:
                    raise ValueError(
                        f"expression for gene {gene.id!r} must be finite and non-negative"
                    )
                lookup[gene.id] = value
                break

    weights: dict[str, float] = {}
    for rxn in model.reactions:
        if not rxn.genes:
            continue
        mapped_value = _eval_gpr_value(rxn.gpr, lookup, or_rule=or_rule)
        if mapped_value is not None:
            weights[rxn.id] = float(mapped_value)
    return weights


def _is_excluded(rxn) -> bool:
    return rxn.boundary or "biomass" in rxn.id.lower() or "growth" in rxn.id.lower()


def _provenance_with_method(
    provenance: Mapping[str, object], method: str
) -> dict[str, object]:
    """Copy of a provenance record whose archived ``parameters.method`` is ``method``."""

    out = dict(provenance)
    parameters = out.get("parameters")
    if isinstance(parameters, Mapping):
        out["parameters"] = {**parameters, "method": method}
    return out


#: Documented departures of CMM's E-Flux2 from Kim et al. (2016), archived with every result.
EFLUX2_DEVIATIONS: tuple[str, ...] = (
    "min_scale floor: scaled bounds are floored at min_scale, so a zero-expression "
    "reaction is never blocked; Kim et al. set the bound to g and therefore block it at g=0",
    "percentile normalisation: bounds are the template bound times "
    "clip(weight / percentile(weights), min_scale, 1); Kim et al. use the absolute "
    "expression value as the bound and apply no normalisation",
)

#: Documented departures of CMM's LAD from Lee et al. (2012), archived with every result.
LAD_DEVIATIONS: tuple[str, ...] = (
    "sigma weights: Lee et al. minimise sum((1/sigma_i)(d+ + d-)); CMM defaults to "
    "sigma_i = 1 (unweighted) because expression tables carry no per-reaction error, and "
    "exposes the weights through reaction_sigma",
    "no iterative direction assignment: Lee et al. resolve reversible-reaction ambiguity "
    "by FVA-driven iteration; CMM solves a single LP and fits |v| via the forward/backward "
    "split, so a reversible reaction's direction is left to the LP",
    "weight_threshold: reactions whose weight is below the threshold are dropped from the "
    "objective rather than driven toward zero flux; the default 0.0 drops nothing",
)


def eflux2(
    model: Model,
    reaction_weights: Mapping[str, float],
    *,
    objective_fraction: float = 1.0,
    min_scale: float = 1e-3,
    weight_threshold: float = 0.0,
    normalization_percentile: float = 100.0,
    exclude_exchange: bool = True,
    excluded_reactions: Iterable[str] | None = None,
    allow_l1_fallback: bool = False,
    gpr_or_rule: str | None = None,
) -> OmicsFluxResult:
    """E-Flux2: scale reaction bounds by normalized expression, then minimize the L2 flux.

    Scales each non-excluded reaction's bounds by ``clip(weight / denom, min_scale, 1)`` where
    ``denom`` is the ``normalization_percentile`` of the above-threshold weights; runs FBA for
    the optimum, holds the objective at ``objective_fraction`` of it, then minimizes the
    L2-norm of the flux vector (QP) for a unique distribution. The published method uses
    ``objective_fraction=1``. An explicitly requested L1 approximation is returned under the
    distinct method name ``eflux2_l1_fallback``.

    ``gpr_or_rule`` does not affect the solve — ``reaction_weights`` are already mapped. It
    records, in provenance, which GPR OR convention produced them; ``integrate_expression``
    fills it in. Left at ``None`` it is archived as ``"unspecified"``.

    Two properties of this implementation are CMM's, not Kim et al.'s, and are archived under
    ``metadata["cmm_deviations"]``:

    - **``min_scale`` floor.** Every scaled bound is floored at ``min_scale`` (default 1e-3),
      so a reaction whose enzyme is not expressed keeps a small residual capacity. Kim et al.
      set the bound to the expression value itself, so ``g = 0`` blocks the reaction outright.
      CMM's E-Flux2 therefore cannot reproduce a Kim-style "off" prediction. Set
      ``min_scale=0`` for the published behaviour, at the risk of an infeasible model.
    - **Percentile normalisation.** Kim et al. use the absolute expression value as the bound
      and explicitly leave the resulting unit mismatch unresolved. CMM divides by the
      ``normalization_percentile`` of the weights (the maximum by default) and multiplies the
      *template* bound by that ratio, which has no counterpart in the paper. Because exchange
      and biomass reactions are excluded from the scaling, the transformation is not
      scale-invariant. Below 100 the ``min(1.0, …)`` cap also erases relative expression
      differences above the chosen percentile.
    """

    if not 0.0 < objective_fraction <= 1.0:
        raise ValueError("objective_fraction must be in (0, 1]")
    if not 0.0 <= min_scale <= 1.0:
        raise ValueError("min_scale must be between 0 and 1")
    if weight_threshold < 0:
        raise ValueError("weight_threshold must be non-negative")
    if not 0.0 < normalization_percentile <= 100.0:
        raise ValueError("normalization_percentile must be in (0, 100]")

    valid: dict[str, float] = {}
    for rid, raw_weight in reaction_weights.items():
        if rid not in model.reactions:
            continue
        weight = float(raw_weight)
        if not np.isfinite(weight) or weight < 0:
            raise ValueError(
                f"reaction weight for {rid!r} must be finite and non-negative"
            )
        valid[rid] = weight
    provenance = run_provenance(
        model,
        method="eflux2",
        objective_fraction=objective_fraction,
        min_scale=min_scale,
        weight_threshold=weight_threshold,
        normalization_percentile=normalization_percentile,
        gpr_or_rule=gpr_or_rule or "unspecified",
    )
    above = {rid: w for rid, w in valid.items() if w > weight_threshold}
    if not above:
        return OmicsFluxResult(
            "eflux2",
            "no_targets",
            None,
            {},
            metadata={**provenance, "cmm_deviations": EFLUX2_DEVIATIONS},
        )

    denom = float(np.percentile(list(above.values()), normalization_percentile)) or 1.0
    excluded = set(excluded_reactions or ())

    with model:
        for rid, weight in above.items():
            rxn = model.reactions.get_by_id(rid)
            if rid in excluded or (exclude_exchange and _is_excluded(rxn)):
                continue
            ratio = min(1.0, max(min_scale, weight / denom))
            lb, ub = rxn.lower_bound, rxn.upper_bound
            rxn.bounds = (lb * ratio if lb < 0 else lb, ub * ratio if ub > 0 else ub)
        # Reactions below threshold are squeezed to the noise floor (E-Flux semantics).
        for rid in valid.keys() - above.keys():
            rxn = model.reactions.get_by_id(rid)
            if rid in excluded or (exclude_exchange and _is_excluded(rxn)):
                continue
            lb, ub = rxn.lower_bound, rxn.upper_bound
            rxn.bounds = (
                lb * min_scale if lb < 0 else lb,
                ub * min_scale if ub > 0 else ub,
            )

        optimum = model.slim_optimize()
        if optimum is None or not np.isfinite(optimum):
            return OmicsFluxResult(
                "eflux2",
                "infeasible",
                None,
                {},
                metadata={**provenance, "cmm_deviations": EFLUX2_DEVIATIONS},
            )

        prob = model.problem
        objective_expr = model.objective.expression
        biological_objective = {
            reaction.id: reaction.objective_coefficient
            for reaction in model.reactions
            if reaction.objective_coefficient != 0
        }
        # Hold the biological objective near its optimum. For a maximized objective that is a
        # lower bound; for a minimized one it is an upper bound. The min-direction band must
        # relax UPWARD from the optimum regardless of sign — dividing a *negative* minimized
        # optimum by objective_fraction tightens past the optimum and makes the QP infeasible,
        # so relax by a fraction of |optimum| instead (equals optimum/fraction for optimum>0).
        if model.objective_direction == "min":
            slack = abs(optimum) * (1.0 / objective_fraction - 1.0)
            floor = prob.Constraint(
                objective_expr, ub=optimum + slack, name="_eflux2_floor"
            )
        else:
            floor = prob.Constraint(
                objective_expr, lb=objective_fraction * optimum, name="_eflux2_floor"
            )
        model.add_cons_vars([floor])

        method = "eflux2"
        if solvers.supports("QP", model.solver.interface):
            l2 = 0.0
            for rxn in model.reactions:
                l2 = l2 + rxn.flux_expression * rxn.flux_expression
            model.objective = prob.Objective(l2, direction="min")
            solution = model.optimize()
            detail = "qp"
        else:
            if not allow_l1_fallback:
                solvers.require(
                    "QP", model.solver.interface, feature="published E-Flux2"
                )
            from cobra.flux_analysis import pfba

            solution = pfba(model)
            detail = "pfba_l1_approximation"
            method = "eflux2_l1_fallback"

        # The provenance record is built before the solver branch, so its archived method
        # parameter must be corrected here or an L1-substituted run archives "eflux2".
        provenance = _provenance_with_method(provenance, method)

        if solution.status != "optimal":
            return OmicsFluxResult(
                method,
                solution.status,
                None,
                {},
                detail,
                {
                    **provenance,
                    "unconstrained_optimum": float(optimum),
                    "cmm_deviations": EFLUX2_DEVIATIONS,
                },
            )
        fluxes = {rid: float(v) for rid, v in solution.fluxes.items()}
        achieved_objective = float(
            sum(
                coefficient * fluxes[rid]
                for rid, coefficient in biological_objective.items()
            )
        )
        return OmicsFluxResult(
            method,
            solution.status,
            achieved_objective,
            fluxes,
            detail,
            {
                **provenance,
                "unconstrained_optimum": float(optimum),
                "objective_fraction": objective_fraction,
                "normalization_percentile": normalization_percentile,
                "min_scale": min_scale,
                "objective_kind": "achieved_biological_objective",
                "cmm_deviations": EFLUX2_DEVIATIONS,
            },
        )


def lad(
    model: Model,
    reaction_weights: Mapping[str, float],
    *,
    scaling_factor: float = 1.0,
    weight_threshold: float = 0.0,
    reaction_sigma: Mapping[str, float] | None = None,
    gpr_or_rule: str | None = None,
) -> OmicsFluxResult:
    """LAD (Lee et al. 2012): fit |flux| to expression targets by least absolute deviation.

    Lee et al.'s Eq. 3 linearisation, ``min sum((1/sigma_i)(d+ + d-))`` subject to
    ``v_i = d_i + d+ - d-`` and the steady-state constraints. The residual is taken against
    the **absolute** flux: each reaction's net flux is split as ``v = f - b`` with
    ``f, b >= 0``, and ``f + b - d+ + d- == target``. COBRApy's forward/reverse variables
    already implement that split with exactly the required degeneracy rules — ``b = 0``
    when ``lb >= 0`` and ``f = 0`` when ``ub <= 0`` — so they are used directly. The problem
    stays an LP and runs on any solver.

    Fitting ``|v|`` rather than signed ``v`` matters because every target is non-negative: a
    reversible reaction operating at ``v = -5`` against a target of 5 scores deviation 0 here
    and 10 under a signed residual, which systematically pushed CMM's solutions before 0.4.0
    away from reverse-direction operation.

    Parameters
    ----------
    scaling_factor
        ``target = |weight| * scaling_factor``. CMM's extension; Lee et al. use the mapped
        value directly (``scaling_factor=1``).
    weight_threshold
        Reactions whose weight is below this are **excluded from the objective** rather than
        driven toward zero flux. Lee et al. include every reaction that has data, so the
        default is 0.0 (exclude nothing); the pre-0.4.0 default of 0.01 discarded the
        method's down-regulation signal.
    reaction_sigma
        Per-reaction ``sigma_i`` — Lee et al.'s *"error in data point i as calculated in the
        mapping of gene expression data to individual reactions"* — giving each residual the
        weight ``1/sigma_i``. **The default (``None``) departs from the paper**: it weights
        every reaction equally (``sigma_i = 1``), because a plain expression table carries no
        per-reaction error estimate. Reactions absent from the mapping keep ``sigma = 1``.
    gpr_or_rule
        Provenance only: which GPR OR convention produced ``reaction_weights``. Filled in by
        ``integrate_expression``; archived as ``"unspecified"`` when omitted.

    Departures from Lee et al. (2012) are archived under ``metadata["cmm_deviations"]``. The
    one not listed there — because it is inherent to fitting ``f + b`` in a single LP — is
    that nothing forces ``f * b = 0``, so for a reversible reaction the LP may satisfy a
    target by carrying flux in both directions at once. This is the behaviour of the
    reference formulation, whose ``delta >= (f + b) - target`` has the same property.
    """

    if scaling_factor <= 0 or not np.isfinite(scaling_factor):
        raise ValueError("scaling_factor must be finite and positive")
    if weight_threshold < 0:
        raise ValueError("weight_threshold must be non-negative")
    targets: dict[str, float] = {}
    for rid, raw_weight in reaction_weights.items():
        if rid not in model.reactions:
            continue
        weight = float(raw_weight)
        if not np.isfinite(weight) or weight < 0:
            raise ValueError(
                f"reaction weight for {rid!r} must be finite and non-negative"
            )
        if weight >= weight_threshold:
            targets[rid] = weight
    sigma: dict[str, float] = {}
    for rid, raw_sigma in (reaction_sigma or {}).items():
        value = float(raw_sigma)
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"reaction_sigma for {rid!r} must be finite and positive")
        sigma[rid] = value
    provenance = run_provenance(
        model,
        method="lad",
        scaling_factor=scaling_factor,
        weight_threshold=weight_threshold,
        sigma_weighted=bool(sigma),
        gpr_or_rule=gpr_or_rule or "unspecified",
    )
    if not targets:
        return OmicsFluxResult(
            "lad",
            "no_targets",
            None,
            {},
            metadata={**provenance, "cmm_deviations": LAD_DEVIATIONS},
        )

    with model:
        prob = model.problem
        deviations = []
        constraints = []
        terms = []
        for rid, weight in targets.items():
            rxn = model.reactions.get_by_id(rid)
            target_flux = weight * scaling_factor
            d_pos = prob.Variable(f"lad_pos_{rid}", lb=0)
            d_neg = prob.Variable(f"lad_neg_{rid}", lb=0)
            deviations.extend([d_pos, d_neg])
            coefficient = 1.0 / sigma.get(rid, 1.0)
            terms.extend([coefficient * d_pos, coefficient * d_neg])
            # |v| = f + b, with cobra's forward/reverse variables supplying f, b >= 0 and
            # collapsing to one of them for a reaction that is already irreversible.
            absolute_flux = rxn.forward_variable + rxn.reverse_variable
            constraints.append(
                prob.Constraint(
                    absolute_flux - d_pos + d_neg,
                    lb=target_flux,
                    ub=target_flux,
                    name=f"lad_{rid}",
                )
            )
        # Add variables first, then the constraints that reference them (optlang-safe order).
        model.add_cons_vars(deviations)
        model.solver.update()
        model.add_cons_vars(constraints)
        model.objective = prob.Objective(sum(terms), direction="min")

        solution = model.optimize()
        if solution.status != "optimal":
            return OmicsFluxResult(
                "lad",
                solution.status,
                None,
                {},
                metadata={
                    **provenance,
                    "objective_kind": "total_absolute_deviation",
                    "cmm_deviations": LAD_DEVIATIONS,
                },
            )
        fluxes = {rid: float(v) for rid, v in solution.fluxes.items()}
        objective = (
            float(solution.objective_value)
            if solution.objective_value is not None
            else None
        )
        return OmicsFluxResult(
            "lad",
            solution.status,
            objective,
            fluxes,
            metadata={
                **provenance,
                "objective_kind": "total_absolute_deviation",
                "scaling_factor": scaling_factor,
                "weight_threshold": weight_threshold,
                "sigma_weighted": bool(sigma),
                "cmm_deviations": LAD_DEVIATIONS,
            },
        )


def integrate_expression(
    model: Model,
    gene_expression: Mapping[str, float],
    *,
    method: str = "eflux2",
    or_rule: str = "sum",
    **kwargs,
) -> OmicsFluxResult:
    """Convenience: gene expression -> reaction weights -> predicted fluxes (E-Flux2 or LAD).

    ``or_rule`` is the GPR isozyme convention (``"sum"``, as both Kim et al. 2016 and Lee
    et al. 2012 specify; ``"max"`` reproduces pre-0.4.0 CMM). It is forwarded to the method
    so the run's provenance states which convention produced its numbers.
    """

    weights = gene_to_reaction_weights(model, gene_expression, or_rule=or_rule)
    if method == "eflux2":
        return eflux2(model, weights, gpr_or_rule=or_rule, **kwargs)
    if method == "lad":
        return lad(model, weights, gpr_or_rule=or_rule, **kwargs)
    raise ValueError(f"unknown omics method {method!r}; use 'eflux2' or 'lad'")
