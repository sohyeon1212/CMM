"""iMAT: an activity pattern for a metabolic network from discretised gene expression.

Shlomi et al. (2008) classify each reaction as highly expressed (R_H), lowly expressed (R_L) or
moderate, then solve one MILP for the flux distribution that agrees with as many of those
labels as possible::

    max  sum_{i in R_H} (y+_i + y-_i)  +  sum_{i in R_L} y+_i
    s.t. S v = 0,  v_min <= v <= v_max
         v_i + y+_i (v_min,i - eps) >= v_min,i          i in R_H   (active forward)
         v_i + y-_i (v_max,i + eps) <= v_max,i          i in R_H   (active backward)
         v_min,i (1 - y+_i) <= v_i <= v_max,i (1 - y+_i) i in R_L   (inactive)

with ``eps = 1`` and no growth objective. The MILP is solved by the MIT-licensed ``imatpy``
package (Griebel, https://github.com/correctjellyfish/imatpy), called through its public API;
CMM passes ``threshold=0`` so its low-expression constraint is exactly the paper's
``v_i = 0`` rather than imatpy's ``|v_i| <= threshold`` relaxation.

Two inputs the iMAT paper leaves to the user are fixed here from the papers that use iMAT for
MTA's reference state:

* **Discretisation** follows Valcárcel et al. (2019, rMTA, section 2.3): in each sample a gene
  is high above ``mean + delta * SD`` and low below ``mean - delta * SD`` of that sample's
  log2 expression over the model's measured genes, and a gene keeps a class only when at least
  ``beta`` of the samples agree. Unmeasured and disagreeing genes are moderate.
* **GPR mapping** uses AND -> min and OR -> max, imatpy's default and the usual reading (a
  complex is as available as its scarcest subunit; isozymes as their most available member).
  The iMAT paper's Methods text states the opposite pairing ("replacing logical 'and' and 'or'
  with 'max' and 'min', respectively"); this choice is a disclosed deviation from that text.

The returned pattern is the input to :func:`cmm.features.revert_reference.imat_reference_state`,
which samples the space consistent with it to build MTA's v_ref (Yizhak et al. 2013).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Literal

import pandas as pd
from cobra import Model

from cmm.core.provenance import run_provenance

#: Shlomi et al. (2008): "epsilon ... was set to 1".
PUBLISHED_IMAT_EPSILON = 1.0
#: rMTA requires every sample to agree on a gene's class (beta = 100%).
PUBLISHED_IMAT_BETA = 1.0
#: The paper's low-expression constraint is v = 0, i.e. imatpy's threshold at zero.
_IMAT_THRESHOLD = 0.0
#: Recorded with every result so the GPR choice is never inferred.
IMAT_GPR_RULE = "AND=min, OR=max"

ReactionActivity = Literal[
    "active_forward", "active_reverse", "inactive", "not_matched", "moderate"
]


@dataclass(frozen=True)
class ExpressionStates:
    """Genes discretised into high (1), moderate (0) and low (-1) expression."""

    states: dict[str, int]
    delta: float
    beta: float
    n_samples: int
    #: Per gene, the fraction of samples calling it high and low.
    high_fraction: dict[str, float] = field(default_factory=dict, repr=False)
    low_fraction: dict[str, float] = field(default_factory=dict, repr=False)

    def summary(self) -> dict[str, float | int]:
        values = list(self.states.values())
        return {
            "n_genes": len(values),
            "n_high": values.count(1),
            "n_moderate": values.count(0),
            "n_low": values.count(-1),
            "n_samples": self.n_samples,
            "delta": self.delta,
            "beta": self.beta,
        }

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "gene_id": gene,
                    "state": state,
                    "high_fraction": self.high_fraction.get(gene, 0.0),
                    "low_fraction": self.low_fraction.get(gene, 0.0),
                }
                for gene, state in sorted(self.states.items())
            ],
            columns=["gene_id", "state", "high_fraction", "low_fraction"],
        )


def discretize_expression(
    log2_expression: pd.DataFrame,
    *,
    delta: float,
    beta: float = PUBLISHED_IMAT_BETA,
    genes: Iterable[str] | None = None,
) -> ExpressionStates:
    """Classify genes as high/moderate/low by the rMTA rule.

    ``log2_expression`` has one row per gene and one column per sample, on the log2 scale.
    ``genes`` restricts the table first; Valcárcel et al. set the thresholds on the
    distribution of metabolic genes, so callers pass the model's measured genes. Each sample's
    mean and standard deviation (``ddof=1``) are taken over the kept genes.
    """

    if delta < 0:
        raise ValueError("delta must be non-negative")
    if not 0.0 < beta <= 1.0:
        raise ValueError("beta must be in (0, 1]")
    table = log2_expression
    if genes is not None:
        wanted = [gene for gene in dict.fromkeys(str(g) for g in genes)]
        table = table.loc[[gene for gene in wanted if gene in table.index]]
    if table.shape[0] < 2 or table.shape[1] < 1:
        raise ValueError(
            "discretisation needs at least two genes and one sample to set thresholds"
        )
    mean = table.mean(axis=0)
    sd = table.std(axis=0, ddof=1)
    high = table.gt(mean + delta * sd, axis=1)
    low = table.lt(mean - delta * sd, axis=1)
    high_fraction = high.mean(axis=1)
    low_fraction = low.mean(axis=1)
    # A small tolerance so beta = 1 with, say, 3 samples is not lost to 1/3 * 3 rounding.
    tol = 1e-12
    states = {
        str(gene): 1
        if high_fraction[gene] >= beta - tol
        else -1
        if low_fraction[gene] >= beta - tol
        else 0
        for gene in table.index
    }
    return ExpressionStates(
        states=states,
        delta=float(delta),
        beta=float(beta),
        n_samples=int(table.shape[1]),
        high_fraction={str(g): float(v) for g, v in high_fraction.items()},
        low_fraction={str(g): float(v) for g, v in low_fraction.items()},
    )


@dataclass(frozen=True)
class IMATResult:
    """One optimal iMAT solution and the activity pattern it realises."""

    status: str
    objective_value: float | None
    epsilon: float
    #: Reaction expression class from the GPR: 1 (R_H), -1 (R_L), 0 (moderate or no GPR).
    reaction_states: dict[str, int]
    #: For R_H and R_L reactions, what the solution did with them; moderate otherwise.
    activity: dict[str, ReactionActivity]
    fluxes: dict[str, float]
    metadata: dict[str, object] = field(default_factory=dict)

    @property
    def n_matched(self) -> int:
        return sum(
            1
            for state in self.activity.values()
            if state in ("active_forward", "active_reverse", "inactive")
        )

    def pattern_bounds(self, model: Model) -> dict[str, tuple[float, float]]:
        """Bounds that confine fluxes to this solution's activity pattern.

        An active R_H reaction keeps its solved direction with ``|v| >= epsilon``; an inactive
        R_L reaction is fixed at zero. Unmatched and moderate reactions keep the model's
        bounds. Fixing the direction is what keeps the space convex, and so samplable.
        """

        bounds: dict[str, tuple[float, float]] = {}
        for rid, state in self.activity.items():
            reaction = model.reactions.get_by_id(rid)
            lb, ub = reaction.lower_bound, reaction.upper_bound
            if state == "active_forward":
                bounds[rid] = (max(lb, self.epsilon), ub)
            elif state == "active_reverse":
                bounds[rid] = (lb, min(ub, -self.epsilon))
            elif state == "inactive":
                bounds[rid] = (max(lb, 0.0), min(ub, 0.0))
        return bounds

    def summary(self) -> dict[str, object]:
        values = list(self.activity.values())
        return {
            "status": self.status,
            "objective_value": self.objective_value,
            "epsilon": self.epsilon,
            "n_high_reactions": sum(1 for v in self.reaction_states.values() if v > 0),
            "n_low_reactions": sum(1 for v in self.reaction_states.values() if v < 0),
            "n_active_forward": values.count("active_forward"),
            "n_active_reverse": values.count("active_reverse"),
            "n_inactive": values.count("inactive"),
            "n_not_matched": values.count("not_matched"),
        }

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "reaction_id": rid,
                    "expression_state": self.reaction_states.get(rid, 0),
                    "activity": self.activity.get(rid, "moderate"),
                    "imat_flux": self.fluxes.get(rid, float("nan")),
                }
                for rid in sorted(self.reaction_states)
            ],
            columns=["reaction_id", "expression_state", "activity", "imat_flux"],
        )


def imat(
    model: Model,
    gene_states: Mapping[str, int],
    *,
    epsilon: float = PUBLISHED_IMAT_EPSILON,
) -> IMATResult:
    """Solve the iMAT MILP on ``model`` as constrained (medium and condition applied).

    ``gene_states`` maps gene ids to 1/0/-1. Model genes missing from it are moderate (0).
    The model is not changed. A non-optimal solve is returned with an empty pattern rather than
    raised, so the caller decides how to report it.
    """

    from cobra.core.solution import get_solution
    from imatpy.imat import add_imat_constraints, add_imat_objective_, flux_to_binary
    from imatpy.parse_gpr import IMAT_FUNC_DICT, gene_to_rxn_weights

    if epsilon <= 0:
        raise ValueError("epsilon must be positive")
    if any(value not in (-1, 0, 1) for value in gene_states.values()):
        raise ValueError("gene states must be -1, 0 or 1")
    weights = pd.Series(
        {gene.id: int(gene_states.get(gene.id, 0)) for gene in model.genes},
        dtype="int64",
    )
    # Every model gene is present, so imatpy has nothing to fill (and nothing to warn about).
    reaction_weights = gene_to_rxn_weights(
        model, weights, fn_dict=IMAT_FUNC_DICT, fill_val=0
    ).astype("int64")
    reaction_states = {str(rid): int(v) for rid, v in reaction_weights.items()}

    provenance = run_provenance(
        model,
        method="imat",
        epsilon=epsilon,
        threshold=_IMAT_THRESHOLD,
        gpr_rule=IMAT_GPR_RULE,
        implementation="imatpy",
        implementation_version=_imatpy_version(),
        n_genes_high=int((weights > 0).sum()),
        n_genes_low=int((weights < 0).sum()),
    )
    # imatpy's constraint and objective builders work on a copy, so ``model`` is untouched.
    # The status is read before any primal value: reading primals from an infeasible solve
    # raises inside the solver interface rather than reporting the status.
    problem = add_imat_constraints(
        model, reaction_weights, epsilon=epsilon, threshold=_IMAT_THRESHOLD
    )
    add_imat_objective_(problem, reaction_weights)
    problem.slim_optimize()
    status = str(problem.solver.status)
    if status != "optimal":
        return IMATResult(
            status=status,
            objective_value=None,
            epsilon=epsilon,
            reaction_states=reaction_states,
            activity={},
            fluxes={},
            metadata=provenance,
        )
    solution = get_solution(problem)
    fluxes = solution.fluxes
    tolerance = model.tolerance
    high = [rid for rid, v in reaction_states.items() if v > 0]
    low = [rid for rid, v in reaction_states.items() if v < 0]
    forward = flux_to_binary(
        fluxes[high], "forward", epsilon, _IMAT_THRESHOLD, tolerance
    )
    reverse = flux_to_binary(
        fluxes[high], "reverse", epsilon, _IMAT_THRESHOLD, tolerance
    )
    inactive = flux_to_binary(
        fluxes[low], "inactive", epsilon, _IMAT_THRESHOLD, tolerance
    )
    activity: dict[str, ReactionActivity] = {}
    for rid in high:
        activity[rid] = (
            "active_forward"
            if forward[rid]
            else "active_reverse"
            if reverse[rid]
            else "not_matched"
        )
    for rid in low:
        activity[rid] = "inactive" if inactive[rid] else "not_matched"
    return IMATResult(
        status="optimal",
        objective_value=float(solution.objective_value),
        epsilon=epsilon,
        reaction_states=reaction_states,
        activity=activity,
        fluxes={str(rid): float(v) for rid, v in fluxes.items()},
        metadata=provenance,
    )


def _imatpy_version() -> str:
    from importlib.metadata import version

    return version("imatpy")
