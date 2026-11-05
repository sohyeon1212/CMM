"""Transformation-target discovery: rank knockouts that move one metabolic state to another.

Where :mod:`cmm.workflows.production` pushes flux toward a product, this workflow asks the
inverse question — given expression for a *source* state and a *target* state, which single
knockout best moves the source toward the target? The methods are Yizhak et al.'s MTA (2013)
and Valcárcel et al.'s rMTA (2019).

Six stages, mirroring the production workflow's numbering so a reader of one recognises the
other:

``01_preflight``    inputs, gene coverage, solver capability
``02_reference``    the source flux state v_ref
``03_direction``    which reactions must change, and which way
``04_candidates``   the knockout universe
``05_transformation`` one optimisation per candidate, ranked
``06_validation``   score decomposition, epsilon sweep, and the optional MOMA baseline

**The reference state follows the papers.** The source expression is discretised by
Valcárcel et al.'s rule (per-sample mean ± ``imat_delta`` SD, ``imat_beta`` consensus), iMAT
finds one activity pattern consistent with it, and v_ref is the mean of
``reference_samples`` OptGP samples from that pattern's space
(:func:`cmm.features.revert_reference.imat_reference_state`). ``epsilon="per_reaction"``
(the default) derives each changed reaction's threshold from the same samples
(:mod:`cmm.features.revert_epsilon`). The remaining differences — the GPR pairing, COBRApy's
fixed warm-up and a single sampling chain — are recorded in the run's provenance and stated in
the report. A scalar ``epsilon`` is still accepted; it is then chosen rather than derived, and
:meth:`TransformationWorkflowConfig.suggest_epsilon` gives it a basis.
A reaction-level run additionally collapses coupled reactions to one representative each, using
full flux coupling; see :mod:`cmm.features.coupling`. A gene-level run, the default, instead
deduplicates genes sharing a blocked-reaction signature.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
import math
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Mapping, Sequence, cast

from cmm.core.condition import Condition
from cmm.core.media import Medium

if TYPE_CHECKING:  # heavy imports stay out of module import time
    import pandas as pd

    from cmm.omics.expression import ExpressionScale

SignificanceTest = Literal["ttest", "fold_change"]
TransformationMethod = Literal["mta", "rmta"]
PerturbationLevel = Literal["gene", "reaction"]
ChangedRanking = Literal["p_value", "fold_change"]
PAdjust = Literal["none", "fdr_bh"]
EpsilonMode = Literal["per_reaction"]

#: Yizhak et al. keep "the top 100-200 most differentially expressed reactions" as the changed
#: set. Each one adds a binary variable to the MIQP, so the cut is also what decides whether a
#: genome-scale run finishes.
PUBLISHED_CHANGED_SET_RANGE = (100, 200)

# Same per-replicate transform used by the desktop expression comparison.
_EXPRESSION_PSEUDOCOUNT = 1.0

#: A linear-declared table whose every value is at most this looks like log2 expression, and
#: would be logged a second time. Real linear counts, TPM or intensities exceed it.
LOG2_LOOKING_MAX = 20.0
#: A log2 value above this would be more than 2**32 linear units, which no expression measure
#: reaches; a log2-declared table that does is almost certainly linear.
LOG2_IMPLAUSIBLE_MAX = 32.0

#: What still differs from Yizhak et al. (2013) and Valcárcel et al. (2019) in building v_ref.
REFERENCE_STATE_DEVIATIONS = (
    "reaction classes use AND=min, OR=max over the GPR, the reverse of the pairing the iMAT "
    "paper's Methods text states",
    "COBRApy's OptGP sampler builds its own warm-up points (two per reaction) where "
    "optGpSampler was run with 5,000",
    "sampling is one seeded chain so the run is reproducible; optGpSampler runs several",
)

#: Config keys of the removed E-Flux2/LAD reference path, rejected with an explanation.
_REMOVED_KEYS = {
    "reference_method": "v_ref is now always built from iMAT",
    "reference_objective_fraction": "iMAT places no objective on growth",
    "epsilon_samples": "renamed to reference_samples",
    "epsilon_thinning": "renamed to reference_thinning",
    "epsilon_seed": "renamed to reference_seed",
}


@dataclass(frozen=True)
class CandidateConfig:
    """How the knockout universe is built.

    The candidate count is the denominator of any "ranked in the top *N*%" statement, so the
    construction is configuration rather than an implementation detail, and it is reported.

    Coupled-set reduction is defined on reactions: Yizhak et al. take "a member from each
    partially coupled set". It therefore applies only when ``perturbation="reaction"``. A
    gene-level run instead deduplicates genes that block the same reaction signature, which is
    the convention the production workflow already uses for its single-knockout screen.
    """

    exclude_blocked: bool = True
    exclude_essential: bool = True
    essential_growth_fraction: float = 0.20  # paper: growth reduced by more than 80%
    #: ``None`` follows the perturbation level, which is what makes the two defaults coherent:
    #: on for a reaction-level run, off for a gene-level one. Set it explicitly to override.
    collapse_coupled_sets: bool | None = None
    explicit: tuple[str, ...] | None = None  # skip construction; use these ids verbatim

    def collapse_for(self, perturbation: PerturbationLevel) -> bool:
        if self.collapse_coupled_sets is None:
            return perturbation == "reaction"
        return self.collapse_coupled_sets

    def validate(self, perturbation: PerturbationLevel) -> None:
        if not 0.0 <= self.essential_growth_fraction < 1.0:
            raise ValueError("essential_growth_fraction must be in [0, 1)")
        if self.collapse_coupled_sets is True and perturbation == "gene":
            raise ValueError(
                "collapse_coupled_sets applies to reaction-level knockouts; a gene-level run "
                "deduplicates on blocked-reaction signature instead. Set perturbation="
                "'reaction' to reproduce the published candidate set, or leave "
                "collapse_coupled_sets unset to follow the perturbation level."
            )
        if self.explicit is not None and not self.explicit:
            raise ValueError("candidates.explicit must be None or a non-empty list")


@dataclass(frozen=True)
class DirectionConfig:
    """How the two expression states become a per-reaction desired direction of change."""

    significance: SignificanceTest = "ttest"
    p_value_cutoff: float = 0.05  # paper: P < 0.05
    up_threshold: float = 1.0  # used only when significance="fold_change"
    down_threshold: float = 1.0
    #: ``"fdr_bh"`` cuts on Benjamini-Hochberg adjusted P (Valcárcel et al.'s prostate cancer
    #: filter); ``"none"`` cuts on the raw P, as Yizhak et al. do.
    p_adjust: PAdjust = "none"
    #: Also require |log2 fold change| >= this for a t-test gene to count as changed. 0 is
    #: Yizhak et al.'s test alone; Valcárcel et al. pair theirs with a 20% change (0.263).
    #: Worth setting with many replicates, where tiny changes reach significance.
    min_abs_log2fc: float = 0.0
    top_n_changed: int | None = 200
    ranking: ChangedRanking = "p_value"

    def validate(self) -> None:
        if not 0.0 < self.p_value_cutoff <= 1.0:
            raise ValueError("p_value_cutoff must be in (0, 1]")
        if self.p_adjust not in ("none", "fdr_bh"):
            raise ValueError("p_adjust must be 'none' or 'fdr_bh'")
        if self.p_adjust != "none" and self.significance != "ttest":
            raise ValueError("p_adjust applies to significance='ttest' only")
        if self.min_abs_log2fc < 0:
            raise ValueError("min_abs_log2fc must be non-negative")
        if self.min_abs_log2fc > 0 and self.significance != "ttest":
            raise ValueError(
                "min_abs_log2fc applies to significance='ttest'; a fold-change run already "
                "cuts on up_threshold/down_threshold"
            )
        if self.up_threshold < 0 or self.down_threshold < 0:
            raise ValueError("fold-change thresholds must be non-negative")
        if self.top_n_changed is not None and self.top_n_changed < 1:
            raise ValueError("top_n_changed must be None or at least 1")
        if self.ranking == "p_value" and self.significance != "ttest":
            raise ValueError(
                "ranking='p_value' needs significance='ttest'; a fold-change run has no "
                "P values to rank on"
            )


@dataclass(frozen=True)
class TransformationValidationConfig:
    """Cross-checks run after the ranking, not further predictions."""

    enabled: bool = True
    #: Yizhak et al. compare MTA against a MOMA baseline and report it as markedly inferior.
    #: Off by default: it costs one QP per candidate, and for rMTA the score decomposition
    #: (rTS against its own MOMA component, no extra solves) answers how much of the ranking
    #: is MOMA alone. Turn it on to reproduce the published comparison; it is a cross-check
    #: and negative control, not a superiority claim.
    run_moma_baseline: bool = False
    #: Empty means no sweep. The papers report epsilon sensitivity (Yizhak Supp. Fig. S6); a
    #: sweep shows whether the threshold, derived or chosen, decided the ranking.
    epsilon_sweep: tuple[float, ...] = ()

    def validate(self) -> None:
        if any(value < 0 for value in self.epsilon_sweep):
            raise ValueError("epsilon_sweep values must be non-negative")


@dataclass(frozen=True)
class TransformationWorkflowConfig:
    """Complete, serializable invocation of transformation-target discovery.

    ``source`` is the state to move *away* from and ``target`` the state to move *toward*.
    Swapping them is a different scientific question with a different answer, and nothing in
    the model can detect the mistake — the agent skill is required to confirm the direction
    rather than infer it from file names.

    ``expression_scale`` says how both expression files are scaled. Both the iMAT
    discretisation and the direction tests work on log2 values (``log2(x + 1)`` for a linear
    file, the values as given for a log2 file).

    ``imat_delta`` has no default: Valcárcel et al. chose it per data set (0.3 for their
    Alzheimer's data, 0.7 for prostate cancer and 0.7-1 for the knockout experiments), so a run
    has to state it. ``epsilon="per_reaction"`` derives each changed reaction's threshold from
    the ``reference_samples`` samples behind v_ref. A number uses that one threshold for every
    reaction.
    """

    model_path: str | Path
    source_expression_path: str | Path
    target_expression_path: str | Path
    output_dir: str | Path | None = None
    solver: str | None = None
    medium: Medium | str | None = None
    condition: Condition | None = None

    expression_scale: ExpressionScale = "linear"
    #: rMTA discretisation (Valcárcel et al. 2019, section 2.3). Required; see the docstring.
    imat_delta: float | None = None
    imat_beta: float = 1.0  # every sample must agree
    imat_epsilon: float = 1.0  # Shlomi et al. (2008)
    #: Sampling behind v_ref and ``epsilon="per_reaction"``. Both papers draw 2,000 samples;
    #: epsilon scales with 1/sqrt(n), so the count is part of the method. Valcárcel et al.
    #: used 500 steps between samples.
    reference_samples: int = 2000
    reference_thinning: int = 500
    reference_seed: int = 0

    direction: DirectionConfig = field(default_factory=DirectionConfig)
    candidates: CandidateConfig = field(default_factory=CandidateConfig)

    method: TransformationMethod = "mta"
    perturbation: PerturbationLevel = "gene"
    alpha: float = 0.66  # paper's main-text value, robust over 0.1-0.9
    epsilon: float | EpsilonMode = "per_reaction"
    parameter_k: float = 100.0  # rMTA Equation 9

    validation: TransformationValidationConfig = field(
        default_factory=TransformationValidationConfig
    )
    overwrite: bool = False

    def __post_init__(self) -> None:
        for name in ("model_path", "source_expression_path", "target_expression_path"):
            object.__setattr__(self, name, Path(getattr(self, name)))
        if self.output_dir is not None:
            object.__setattr__(self, "output_dir", Path(self.output_dir))
        self.validate()

    def validate(self) -> None:
        for name in ("model_path", "source_expression_path", "target_expression_path"):
            if not str(getattr(self, name)):
                raise ValueError(f"{name} must not be empty")
        if self.source_expression_path == self.target_expression_path:
            raise ValueError(
                "source and target expression must be different files; the transformation "
                "from a state to itself is not defined"
            )
        if self.imat_delta is None:
            raise ValueError(
                "imat_delta is required: Valcárcel et al. set it per data set (0.3 for "
                "Alzheimer's, 0.7 for prostate cancer), so the run must state it"
            )
        if self.imat_delta < 0:
            raise ValueError("imat_delta must be non-negative")
        if not 0.0 < self.imat_beta <= 1.0:
            raise ValueError("imat_beta must be in (0, 1]")
        if self.imat_epsilon <= 0:
            raise ValueError("imat_epsilon must be positive")
        if self.method not in ("mta", "rmta"):
            raise ValueError("method must be 'mta' or 'rmta'")
        if self.perturbation not in ("gene", "reaction"):
            raise ValueError("perturbation must be 'gene' or 'reaction'")
        if not 0.0 <= self.alpha <= 1.0:
            raise ValueError("alpha must be between 0 and 1")
        if self.expression_scale not in ("linear", "log2"):
            raise ValueError("expression_scale must be 'linear' or 'log2'")
        if isinstance(self.epsilon, str):
            if self.epsilon != "per_reaction":
                raise ValueError("epsilon must be a number or 'per_reaction'")
        elif self.epsilon < 0:
            raise ValueError("epsilon must be non-negative")
        if self.reference_samples < 2:
            raise ValueError("reference_samples must be at least 2")
        if self.reference_thinning < 1:
            raise ValueError("reference_thinning must be at least 1")
        if self.parameter_k <= 0:
            raise ValueError("parameter_k must be positive")
        self.direction.validate()
        self.candidates.validate(self.perturbation)
        self.validation.validate()

    @property
    def follows_published_changed_set_size(self) -> bool:
        """Whether the changed-set cut sits in the range Yizhak et al. report using."""

        low, high = PUBLISHED_CHANGED_SET_RANGE
        return self.top_n_changed is not None and low <= self.top_n_changed <= high

    @property
    def top_n_changed(self) -> int | None:
        return self.direction.top_n_changed

    @staticmethod
    def suggest_epsilon(reference_fluxes: Mapping[str, float]) -> dict[str, float]:
        """Epsilon candidates derived from the reference state's own flux magnitudes.

        For a run with a single numeric epsilon (``epsilon="per_reaction"`` derives one per
        reaction instead). A single epsilon is a *flux magnitude*, so it must be read off the
        model at hand rather than carried over from another one.

        Returns percentiles of the non-zero |v_ref| so a caller — or an agent interview — can
        offer a number with a basis instead of a default nobody chose.
        """

        import numpy as np

        magnitudes = np.array(
            [abs(v) for v in reference_fluxes.values() if abs(v) > 1e-12]
        )
        if not magnitudes.size:
            return {}
        return {
            "p10": float(np.percentile(magnitudes, 10)),
            "p25": float(np.percentile(magnitudes, 25)),
            "median": float(np.median(magnitudes)),
            "p75": float(np.percentile(magnitudes, 75)),
        }

    def to_provenance(self) -> dict[str, object]:
        """Everything a reader needs to tell this run apart from a differently-configured one."""

        from cmm.omics.imat import IMAT_GPR_RULE

        return {
            "model_path": str(self.model_path),
            "source_expression": str(self.source_expression_path),
            "target_expression": str(self.target_expression_path),
            "reference_method": "imat",
            "expression_scale": self.expression_scale,
            "imat_delta": self.imat_delta,
            "imat_beta": self.imat_beta,
            "imat_epsilon": self.imat_epsilon,
            "imat_gpr_rule": IMAT_GPR_RULE,
            "imat_discretisation_transform": "value"
            if self.expression_scale == "log2"
            else "log2(value + pseudocount)",
            "reference_samples": self.reference_samples,
            "reference_thinning": self.reference_thinning,
            "reference_seed": self.reference_seed,
            "reference_sampler": "optgp (COBRApy), one seeded chain",
            "direction_transform": "value"
            if self.expression_scale == "log2"
            else "log2(value + pseudocount)",
            "direction_pseudocount": _EXPRESSION_PSEUDOCOUNT,
            "significance": self.direction.significance,
            "p_value_cutoff": self.direction.p_value_cutoff,
            "p_adjust": self.direction.p_adjust,
            "min_abs_log2fc": self.direction.min_abs_log2fc,
            "up_threshold": self.direction.up_threshold,
            "down_threshold": self.direction.down_threshold,
            "top_n_changed": self.direction.top_n_changed,
            "changed_ranking": self.direction.ranking,
            "method": self.method,
            "perturbation": self.perturbation,
            "alpha": self.alpha,
            "epsilon": self.epsilon,
            "parameter_k": self.parameter_k,
            "candidate_exclude_blocked": self.candidates.exclude_blocked,
            "candidate_exclude_essential": self.candidates.exclude_essential,
            "candidate_essential_growth_fraction": self.candidates.essential_growth_fraction,
            "candidate_collapse_coupled_sets": self.candidates.collapse_for(
                self.perturbation
            ),
            "candidate_explicit": list(self.candidates.explicit)
            if self.candidates.explicit
            else None,
            "epsilon_sweep": list(self.validation.epsilon_sweep),
            "run_moma_baseline": self.validation.run_moma_baseline,
            # A Medium object is recorded by name here; the exchanges it actually set are
            # recorded beside it as ``medium_applied`` once the model has been conditioned.
            "medium": self.medium.name
            if isinstance(self.medium, Medium)
            else self.medium,
            "condition": self.condition.name if self.condition else None,
            # Stated on every run so a reader never has to infer them.
            "reference_state_deviations": list(REFERENCE_STATE_DEVIATIONS),
        }

    @classmethod
    def from_json(cls, path: str | Path) -> TransformationWorkflowConfig:
        """Load from a UTF-8 JSON object, resolving relative paths against the config file."""

        config_path = Path(path).expanduser().resolve()
        payload = json.loads(config_path.read_text(encoding="utf-8"))
        if not isinstance(payload, Mapping):
            raise ValueError("transformation workflow JSON must contain an object")
        values = dict(payload)
        for name in (
            "model_path",
            "source_expression_path",
            "target_expression_path",
            "output_dir",
        ):
            raw = values.get(name)
            if raw is None:
                continue
            candidate = Path(str(raw)).expanduser()
            if not candidate.is_absolute():
                candidate = config_path.parent / candidate
            values[name] = candidate.resolve()
        return cls.from_mapping(values)

    @classmethod
    def from_mapping(
        cls, payload: Mapping[str, object]
    ) -> TransformationWorkflowConfig:
        """Build a config from CLI-friendly mappings, including the nested sections."""

        from cmm.workflows.production import (
            _condition_from_payload,
            _medium_from_payload,
        )

        values = dict(payload)
        removed = sorted(key for key in _REMOVED_KEYS if key in values)
        if removed:
            raise ValueError(
                "unsupported transformation config keys: "
                + "; ".join(f"{key} ({_REMOVED_KEYS[key]})" for key in removed)
            )
        values["medium"] = _medium_from_payload(values.get("medium"))
        values["condition"] = _condition_from_payload(values.get("condition"))
        for key, factory in (
            ("direction", DirectionConfig),
            ("candidates", CandidateConfig),
            ("validation", TransformationValidationConfig),
        ):
            section = values.get(key)
            if isinstance(section, Mapping):
                section_values = dict(section)
                for tuple_field in ("explicit", "epsilon_sweep"):
                    raw = section_values.get(tuple_field)
                    if isinstance(raw, Sequence) and not isinstance(raw, (str, bytes)):
                        section_values[tuple_field] = tuple(raw)
                values[key] = factory(**section_values)  # type: ignore[arg-type]
        return cls(**values)  # type: ignore[arg-type]


# --- execution ---------------------------------------------------------------------------

_STAGE_DIRECTORIES = (
    "01_preflight",
    "02_reference",
    "03_direction",
    "04_candidates",
    "05_transformation",
    "06_validation",
    "model",
    "inputs",
    "scripts",
)


class TransformationWorkflowError(RuntimeError):
    """The transformation workflow could not run as configured."""


@dataclass(frozen=True)
class PreflightRecord:
    check: str
    status: Literal["pass", "warning"]
    value: object
    message: str


@dataclass(frozen=True)
class TransformationWorkflowResult:
    """Everything one transformation run produced, plus where it was written."""

    config: TransformationWorkflowConfig
    provenance: Mapping[str, object]
    preflight: tuple[PreflightRecord, ...]
    reference_fluxes: Mapping[str, float]
    direction: Mapping[str, int]
    candidates: tuple[str, ...]
    ranking: tuple[Mapping[str, object], ...]
    moma_baseline: tuple[Mapping[str, object], ...] = ()
    epsilon_sweep: tuple[Mapping[str, object], ...] = ()
    candidate_filtering: Mapping[str, object] = field(default_factory=dict)
    direction_summary: Mapping[str, object] = field(default_factory=dict)
    #: One row per changed reaction when ``epsilon="per_reaction"``; empty for a scalar.
    reaction_epsilons: tuple[Mapping[str, object], ...] = ()
    #: The source genes' iMAT classes (high 1, moderate 0, low -1) and sample agreement.
    gene_states: tuple[Mapping[str, object], ...] = ()
    #: One row per reaction: iMAT class and activity, pattern bounds imposed, and v_ref.
    reference_pattern: tuple[Mapping[str, object], ...] = ()
    #: rMTA only: each candidate's four scores, its Eq. 9 branch, and its rank under rTS
    #: and under the MOMA component alone. Computed from the ranking, with no extra solve.
    score_decomposition: tuple[Mapping[str, object], ...] = ()
    #: Tie structure, amplified count and rank correlations of the ranking.
    ranking_structure: Mapping[str, object] = field(default_factory=dict)
    #: How much of the model the expression covers, and how open the uptake side is.
    data_coverage: Mapping[str, object] = field(default_factory=dict)
    run_directory: Path | None = None

    def summary(self) -> dict[str, object]:
        changed = sum(1 for value in self.direction.values() if value != 0)
        best = self.ranking[0] if self.ranking else None
        return {
            "method": self.config.method,
            "perturbation": self.config.perturbation,
            "n_candidates": len(self.candidates),
            "n_scored": len(self.ranking),
            "n_changed_reactions": changed,
            "reference_method": "imat",
            "top_target": best["target_id"] if best else None,
            "top_score": best["score"] if best else None,
            # The count is the denominator of any "ranked in the top N%" reading, so it
            # travels with the result rather than being recomputed by a reader.
            "candidate_construction": dict(self.candidate_filtering),
            "direction_construction": dict(self.direction_summary),
            "ranking_structure": dict(self.ranking_structure),
            "data_coverage": dict(self.data_coverage),
        }


def _read_expression(path: Path, scale: str = "linear") -> "pd.DataFrame":
    """Read a gene-expression table: first column gene ids, the rest measurements.

    Two columns is one measurement per gene, which supports a fold-change cut only. Three or
    more are replicates, which is what the published Student's t-test needs. The shape is a
    fact about the file, so it is inspected and reported rather than asked about.
    """

    import numpy as np

    from cmm.omics.expression import describe_expression_columns, read_gene_table

    frame = read_gene_table(path, index=True)
    numeric = frame.select_dtypes("number").copy()
    if numeric.empty:
        raise TransformationWorkflowError(
            f"{path.name} has no numeric expression columns after the gene id column"
        )
    # Kept on the frame so preflight can report which columns became samples.
    numeric.attrs["columns"] = describe_expression_columns(frame)
    if numeric.isna().any().any():
        raise TransformationWorkflowError(f"{path.name} carries missing values")
    if not np.isfinite(numeric.to_numpy()).all():
        raise TransformationWorkflowError(
            f"{path.name} must contain finite expression values"
        )
    # A log2 table legitimately goes below zero; a linear one cannot.
    if scale == "linear" and (numeric < 0).any().any():
        raise TransformationWorkflowError(
            f"{path.name} has negative values, so it is not linear expression; if the file "
            "is log2-scaled, set expression_scale='log2'"
        )
    numeric.index = [str(value) for value in numeric.index]
    if numeric.index.has_duplicates:
        raise TransformationWorkflowError(f"{path.name} repeats a gene id")
    return numeric


def _check_expression_scale(check, role, table, scale):
    """Warn when a file's values contradict its declared scale.

    A linear-declared file with negative values is already refused. What still passes silently
    is a log2 file declared linear (all values small and positive), which would be logged twice,
    and a linear file declared log2. Neither is refused: the values only suggest the mistake.
    """

    low = float(table.to_numpy().min())
    high = float(table.to_numpy().max())
    value = {"declared": scale, "min": low, "max": high}
    if scale == "linear":
        ok = high > LOG2_LOOKING_MAX
        message = (
            f"values exceed {LOG2_LOOKING_MAX:g}, consistent with linear expression"
            if ok
            else f"every value lies between {low:g} and {high:g}, which looks like log2 "
            "expression; a log2 file declared linear is logged a second time. If the file "
            "is log2, set expression_scale='log2'"
        )
    else:
        ok = high <= LOG2_IMPLAUSIBLE_MAX
        message = (
            f"values stay at or below {LOG2_IMPLAUSIBLE_MAX:g}, consistent with log2 expression"
            if ok
            else f"values reach {high:g}, more than any log2 expression measure; if the file "
            "is linear, set expression_scale='linear'"
        )
    check(f"{role}_expression_scale", ok, value, message)


def _imat_state(model, expression, model_genes, config, role):
    """Discretise one state's expression and build its iMAT reference state.

    The table is put on the log2 scale the rMTA rule is stated on, restricted to the model's
    measured genes (the thresholds are set on metabolic genes), then discretised and sampled.
    """

    from cmm.features.revert_reference import IMATReferenceError, imat_reference_state
    from cmm.omics.expression import log2_expression
    from cmm.omics.imat import discretize_expression

    assert config.imat_delta is not None  # validate() requires it
    log2 = log2_expression(expression, config.expression_scale, _EXPRESSION_PSEUDOCOUNT)
    states = discretize_expression(
        log2,
        delta=config.imat_delta,
        beta=config.imat_beta,
        genes=sorted(model_genes & set(log2.index)),
    )
    try:
        state = imat_reference_state(
            model,
            states.states,
            epsilon=config.imat_epsilon,
            n_samples=config.reference_samples,
            thinning=config.reference_thinning,
            seed=config.reference_seed,
        )
    except IMATReferenceError as exc:
        raise TransformationWorkflowError(
            f"the {role} expression produced no iMAT reference state: {exc}"
        ) from exc
    return states, state


def _gene_directions(
    source: "pd.DataFrame",
    target: "pd.DataFrame",
    config: DirectionConfig,
    scale: "ExpressionScale" = "linear",
) -> "pd.DataFrame":
    """Compare expression in log2 space, whichever scale the files are stored on."""

    from cmm.omics.differential import (
        gene_directions_by_fold_change,
        gene_directions_from_replicates,
    )
    from cmm.omics.expression import log2_expression

    source_log2 = log2_expression(source, scale, _EXPRESSION_PSEUDOCOUNT)
    target_log2 = log2_expression(target, scale, _EXPRESSION_PSEUDOCOUNT)
    try:
        if config.significance == "ttest":
            return gene_directions_from_replicates(
                source_log2,
                target_log2,
                p_value_cutoff=config.p_value_cutoff,
                p_adjust=config.p_adjust,
                min_abs_log2fc=config.min_abs_log2fc,
            )
        return gene_directions_by_fold_change(
            source_log2,
            target_log2,
            up_threshold=config.up_threshold,
            down_threshold=config.down_threshold,
        )
    except ValueError as error:
        raise TransformationWorkflowError(str(error)) from error


def _build_candidates(
    model,
    config: TransformationWorkflowConfig,
) -> tuple[tuple[str, ...], dict[str, object]]:
    """The knockout universe for this config; see :func:`cmm.features.revert_candidates`."""

    from cmm.features.revert import revert_candidates

    return revert_candidates(
        model,
        perturbation=config.perturbation,
        exclude_blocked=config.candidates.exclude_blocked,
        exclude_essential=config.candidates.exclude_essential,
        essential_growth_fraction=config.candidates.essential_growth_fraction,
        collapse_coupled_sets=config.candidates.collapse_for(config.perturbation),
        explicit=config.candidates.explicit,
    )


#: What the run reads with ``cobra.io.read_sbml_model`` and archives byte-for-byte.
_SBML_SUFFIXES = (".xml", ".sbml", ".xml.gz", ".sbml.gz")


def _sbml_model_path(path: str | Path) -> Path:
    """The model file, refused up front unless it is SBML.

    The bundle archives the model as given and replays it with the SBML reader, so a COBRApy
    JSON or MATLAB model would otherwise fail deep inside the parser with no hint of the fix.
    """

    model_path = Path(path).expanduser()
    if not model_path.is_file():
        raise FileNotFoundError(f"model_path is not a file: {model_path}")
    if not model_path.name.lower().endswith(_SBML_SUFFIXES):
        raise ValueError(
            f"transformation workflow reads SBML models (.xml/.sbml); got {model_path.name}. "
            "Convert a COBRApy JSON model first: cobra.io.write_sbml_model("
            "cobra.io.load_json_model(path), 'model.xml')"
        )
    return model_path


def run_transformation_target_discovery(
    config: TransformationWorkflowConfig,
) -> TransformationWorkflowResult:
    """Run the six stages and, when ``output_dir`` is set, write a schema-v2 run bundle."""

    from cobra.io import read_sbml_model

    from cmm.core import solvers
    from cmm.core.provenance import run_provenance
    from cmm.features.revert import revert_targets
    from cmm.features.revert_epsilon import epsilons_from_samples
    from cmm.omics.differential import (
        reaction_directions,
        restrict_to_top_changed,
    )

    model = read_sbml_model(str(_sbml_model_path(config.model_path)))
    if config.solver:
        model.solver = config.solver

    preflight: list[PreflightRecord] = []

    def check(name: str, ok: bool, value: object, message: str) -> None:
        preflight.append(
            PreflightRecord(name, "pass" if ok else "warning", value, message)
        )

    # --- 01_preflight -------------------------------------------------------------------
    # The MIQP gate is a stop, not a warning: rmta_continuous is a different method and must
    # never stand in for a published one that the solver cannot run.
    solvers.require(
        "MIQP", model.solver.interface, feature=f"published {config.method.upper()}"
    )

    source_expression = _read_expression(
        Path(config.source_expression_path), config.expression_scale
    )
    target_expression = _read_expression(
        Path(config.target_expression_path), config.expression_scale
    )
    check(
        "replicates",
        min(source_expression.shape[1], target_expression.shape[1]) >= 2,
        [source_expression.shape[1], target_expression.shape[1]],
        "at least two replicates per state are needed for the published t-test",
    )
    for role, table in (("source", source_expression), ("target", target_expression)):
        _check_expression_scale(check, role, table, config.expression_scale)
        columns = table.attrs["columns"]
        notes = columns.warnings()
        check(
            f"{role}_sample_columns",
            not notes,
            {
                "samples": list(columns.samples),
                "ignored": list(columns.dropped),
                "identifier_like": list(columns.identifier_like),
            },
            "; ".join(notes)
            if notes
            else "every column after the gene id was read as one sample",
        )

    model_genes = {gene.id for gene in model.genes}
    measured = set(source_expression.index) & set(target_expression.index)
    covered = model_genes & measured
    check(
        "gene_coverage",
        bool(model_genes) and len(covered) / len(model_genes) >= 0.5,
        {"model_genes": len(model_genes), "measured_in_model": len(covered)},
        "expression must cover the model's genes; a small overlap usually means the two use "
        "different identifier systems",
    )
    if not covered:
        raise TransformationWorkflowError(
            "no model gene appears in both expression files; check the identifier system"
        )

    medium_applied: dict[str, object] | None = None
    with model:
        if config.medium is not None:
            from cmm.core.media import apply_medium

            medium_applied = apply_medium(model, config.medium).to_provenance()
        if config.condition is not None:
            config.condition.apply_to(model)
        # A bound read, not a solve: on a genome-scale model "model as loaded" can leave most
        # exchanges open for uptake, which is a scientific statement rather than a default.
        exchanges = list(model.exchanges)
        data_coverage: dict[str, object] = {
            "model_genes": len(model_genes),
            "genes_measured_in_model": len(covered),
            "gene_coverage": len(covered) / len(model_genes) if model_genes else 0.0,
            "n_exchanges": len(exchanges),
            "n_exchanges_open_for_uptake": sum(
                1 for reaction in exchanges if _permits_uptake(reaction)
            ),
        }

        growth = model.slim_optimize()
        check(
            "growth",
            growth == growth and growth > 1e-9,
            float(growth) if growth == growth else None,
            "the model must grow under the requested condition",
        )

        # --- 02_reference ---------------------------------------------------------------
        gene_states, imat_state = _imat_state(
            model, source_expression, model_genes, config, "source"
        )
        reference = imat_state.to_flux_state("source_state_imat")
        check(
            "imat_pattern",
            imat_state.imat.n_matched > 0,
            imat_state.imat.summary(),
            "iMAT must match at least one highly or lowly expressed reaction; none means "
            "the discretisation left every gene moderate",
        )

        # --- 03_direction ---------------------------------------------------------------
        genes = _gene_directions(
            source_expression,
            target_expression,
            config.direction,
            config.expression_scale,
        )
        gene_dirs = {gene: int(value) for gene, value in genes["direction"].items()}
        direction = reaction_directions(model, gene_dirs, reference=reference)
        n_labelled = len(direction.nonsteady())
        if config.direction.top_n_changed is not None:
            direction = restrict_to_top_changed(
                model,
                direction,
                genes["log2_fold_change"].to_dict(),
                config.direction.top_n_changed,
                gene_p_values=genes["p_value"].to_dict()
                if config.direction.ranking == "p_value"
                else None,
            )
        if not direction.nonsteady():
            raise TransformationWorkflowError(
                "no reaction was labelled as changed; loosen the significance threshold or "
                "check that the two expression states differ"
            )
        n_changed = len(direction.nonsteady())
        low, high = PUBLISHED_CHANGED_SET_RANGE
        # The papers tune their filters until 100-200 reactions are changed; whether a run
        # lands there is a fact about the data, so it is checked on the actual count.
        check(
            "changed_set_size",
            low <= n_changed <= high,
            n_changed,
            f"Yizhak et al. and Valcárcel et al. use {low}-{high} changed reactions; outside "
            "that range, adjust p_value_cutoff, p_adjust or top_n_changed and say so",
        )
        direction_summary = {
            "significance": config.direction.significance,
            "p_adjust": config.direction.p_adjust,
            "min_abs_log2fc": config.direction.min_abs_log2fc,
            "n_genes_compared": int(len(genes)),
            "n_genes_significant": int(genes["significant"].sum()),
            "n_reactions_labelled": n_labelled,
            "n_reactions_changed": n_changed,
            "changed_set_in_published_range": low <= n_changed <= high,
            "ranking": config.direction.ranking,
            "top_n_changed": config.direction.top_n_changed,
            "gene_states": gene_states.summary(),
            "imat": imat_state.imat.summary(),
            **dict(direction.metadata),
        }

        # --- 04_candidates --------------------------------------------------------------
        candidates, filtering = _build_candidates(model, config)
        if not candidates:
            raise TransformationWorkflowError(
                "the candidate filters removed every knockout; relax exclude_blocked or "
                "exclude_essential"
            )

        # --- 05_transformation ----------------------------------------------------------
        epsilon: float | Mapping[str, float]
        epsilon_rows: tuple[Mapping[str, object], ...] = ()
        if config.epsilon == "per_reaction":
            assert imat_state.samples is not None
            thresholds = epsilons_from_samples(
                imat_state.samples,
                sorted(direction.nonsteady()),
                metadata=imat_state.metadata,
            )
            epsilon = thresholds.epsilons
            epsilon_rows = tuple(thresholds.to_frame().to_dict("records"))
            direction_summary["epsilon_summary"] = thresholds.summary()
        else:
            epsilon = config.epsilon
        ranking = revert_targets(
            model,
            None,
            reference,
            direction,
            targets=candidates,
            method=config.method,
            alpha=config.alpha,
            epsilon=epsilon,
            parameter_k=config.parameter_k,
            perturbation=config.perturbation,
        )
        # Only rMTA computes bTS, mTS and wTS as separate quantities. For MTA the same score
        # is returned in all four slots, so emitting them as columns would present one number
        # as three independent measurements that happen to agree.
        components = ("bTS", "mTS", "wTS") if config.method == "rmta" else ()
        # The name column appears only when the model named something, so a reconstruction
        # that names nothing writes exactly the table it wrote before.
        named = any(target.target_name for target in ranking.targets)
        rows = tuple(
            {
                "target_id": target.target_id,
                **({"target_name": target.target_name} if named else {}),
                "score": float(target.score),
                "rank": index + 1,
                **{
                    key: float((target.detail or {})[key])
                    for key in components
                    if isinstance((target.detail or {}).get(key), (int, float))
                },
            }
            for index, target in enumerate(ranking.sorted().targets)
        )

        # --- 06_validation --------------------------------------------------------------
        baseline: tuple[Mapping[str, object], ...] = ()
        sweep: tuple[Mapping[str, object], ...] = ()
        if config.validation.enabled:
            if config.validation.run_moma_baseline:
                baseline = _moma_baseline(
                    model,
                    reference,
                    target_expression,
                    candidates,
                    config,
                )
            sweep = _epsilon_sweep(
                model, reference, direction, candidates, config, rows
            )

    decomposition = _score_decomposition(rows) if config.method == "rmta" else ()
    structure = _ranking_structure(ranking.metadata, rows, decomposition, baseline)

    provenance = {
        **run_provenance(model, method=f"transformation_{config.method}"),
        **config.to_provenance(),
        "reference_provenance": reference.provenance,
        "reference_state": dict(imat_state.metadata),
        "direction_provenance": dict(direction.metadata),
        "ranking_provenance": dict(ranking.metadata),
        **({"medium_applied": medium_applied} if medium_applied is not None else {}),
    }
    result = TransformationWorkflowResult(
        config=config,
        provenance=provenance,
        preflight=tuple(preflight),
        reference_fluxes=dict(reference.fluxes),
        direction=dict(direction.directions),
        candidates=candidates,
        ranking=rows,
        moma_baseline=baseline,
        epsilon_sweep=sweep,
        candidate_filtering=filtering,
        direction_summary=direction_summary,
        reaction_epsilons=epsilon_rows,
        gene_states=tuple(gene_states.to_frame().to_dict("records")),
        reference_pattern=tuple(imat_state.to_frame().to_dict("records")),
        score_decomposition=decomposition,
        ranking_structure=structure,
        data_coverage=data_coverage,
        run_directory=None,
    )
    if config.output_dir is None:
        return result
    return _export(result, genes, model)


def _write_reproduction_scripts(
    writer: Any,
    result: TransformationWorkflowResult,
    *,
    root: Path,
    model_relative: str,
    expression_relatives: Mapping[str, str],
) -> None:
    """Write the three entry points that make the directory replayable on its own.

    Months later the run directory is all that survives; a bundle that cannot say how to
    re-run, re-render or re-check itself leaves that to whoever finds it.
    """

    from dataclasses import asdict

    from cmm.workflows._bundle import _jsonable

    config_payload = cast(dict[str, object], _jsonable(asdict(result.config)))
    config_payload["model_path"] = f"../{model_relative}"
    config_payload.update(
        {field: f"../{relative}" for field, relative in expression_relatives.items()}
    )
    config_payload["output_dir"] = f"../../{root.name}__reproduced"
    config_payload["overwrite"] = False
    writer.json(
        "scripts/transformation_config.json",
        config_payload,
        stage="scripts",
        role="reproduction_config",
    )
    writer.text(
        "scripts/reproduce.py",
        """#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

from cmm.workflows.transformation import (
    TransformationWorkflowConfig,
    run_transformation_target_discovery,
)

HERE = Path(__file__).resolve().parent


def main() -> int:
    parser = argparse.ArgumentParser(description="Reproduce this CMM transformation run")
    parser.add_argument("--output-dir", type=Path, help="Override the sibling output directory")
    args = parser.parse_args()
    config = TransformationWorkflowConfig.from_json(HERE / "transformation_config.json")
    if args.output_dir is not None:
        config = replace(config, output_dir=args.output_dir.resolve(), overwrite=False)
    result = run_transformation_target_discovery(config)
    print(result.run_directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
""",
        stage="scripts",
        role="reproduce_script",
        executable=True,
    )
    writer.text(
        "scripts/render.py",
        """#!/usr/bin/env python3
from pathlib import Path

from cmm.reporting import render_transformation_report

RUN_DIR = Path(__file__).resolve().parents[1]
report = render_transformation_report(RUN_DIR)
print(report.report_standalone_html)
""",
        stage="scripts",
        role="render_script",
        executable=True,
    )
    writer.text(
        "scripts/validate.py",
        """#!/usr/bin/env python3
import json
from pathlib import Path

from cmm.reporting import validate_transformation_run

RUN_DIR = Path(__file__).resolve().parents[1]
report = validate_transformation_run(RUN_DIR)
print(json.dumps({"valid": report.valid, "issues": report.issues, "warnings": report.warnings}, indent=2))
raise SystemExit(0 if report.valid else 1)
""",
        stage="scripts",
        role="validate_script",
        executable=True,
    )


def _moma_baseline(model, reference, target_expression, candidates, config):
    """The MOMA comparison Yizhak et al. report as markedly inferior to MTA.

    It needs a target *flux* vector rather than a direction map, so the target expression is
    put through the same discretisation, iMAT and sampling that produced v_ref. Reproducing
    the contrast is what shows a ranking's signal comes from the method rather than from the
    inputs.
    """

    from cmm.features.revert_reference import IMATReferenceError
    from cmm.features.transformation import transformation_targets

    # Same construction as v_ref, so the two states differ only in their expression.
    try:
        _, imat_state = _imat_state(
            model,
            target_expression,
            {gene.id for gene in model.genes},
            config,
            "target",
        )
    except (IMATReferenceError, TransformationWorkflowError):
        return ()
    target_state = imat_state.to_flux_state("target_state_imat")
    ranking = transformation_targets(
        model,
        reference,
        target_state,
        method="moma",
        perturbation=config.perturbation,
        targets=list(candidates),
        order=2,
    )
    named = any(t.target_name for t in ranking.targets)
    return tuple(
        {
            "target_id": t.target_id,
            **({"target_name": t.target_name} if named else {}),
            "moma_score": float(t.score),
            "status": t.status,
            "rank": i + 1,
        }
        for i, t in enumerate(ranking.sorted().targets)
    )


def _permits_uptake(reaction: Any) -> bool:
    """Whether a boundary reaction can bring its metabolite into the model.

    The usual exchange is written ``A <=>`` and takes up with negative flux; a supply written
    ``-> A`` takes up with positive flux. Reading the orientation keeps both counted.
    """

    net = sum(reaction.metabolites.values())
    return (net < 0 and reaction.lower_bound < 0) or (
        net > 0 and reaction.upper_bound > 0
    )


def _score_decomposition(
    rows: Sequence[Mapping[str, object]],
) -> tuple[Mapping[str, object], ...]:
    """Split each rMTA score into its Equation 9 branch, with no extra solve.

    When the amplification branch (``bTS > 0``, ``mTS > 0``, ``wTS < 0``) does not fire,
    ``rTS = mTS``; ranking the same candidates by ``mTS`` alone therefore shows how much of
    the ordering is the MOMA component. ``rank_mTS`` breaks ties on target id, the same rule
    the ranking itself uses, so the two ranks are comparable.
    """

    scored = []
    for row in rows:
        values = [row.get(key) for key in ("bTS", "mTS", "wTS", "score")]
        if all(isinstance(v, (int, float)) and math.isfinite(v) for v in values):
            scored.append(row)
    by_moma = sorted(scored, key=lambda r: (-float(r["mTS"]), str(r["target_id"])))  # type: ignore[arg-type]
    moma_rank = {str(r["target_id"]): i + 1 for i, r in enumerate(by_moma)}
    decomposition = []
    for row in scored:
        b, m, w = (float(row[key]) for key in ("bTS", "mTS", "wTS"))  # type: ignore[arg-type]
        target = str(row["target_id"])
        rank_rts = int(row["rank"])  # type: ignore[call-overload]
        decomposition.append(
            {
                "target_id": target,
                **({"target_name": row["target_name"]} if "target_name" in row else {}),
                "bTS": b,
                "mTS": m,
                "wTS": w,
                "rTS": float(row["score"]),  # type: ignore[arg-type]
                "amplified": b > 0 and m > 0 and w < 0,
                "rank_rTS": rank_rts,
                "rank_mTS": moma_rank[target],
                "rank_delta": moma_rank[target] - rank_rts,
            }
        )
    return tuple(decomposition)


def _spearman(left: Sequence[float], right: Sequence[float]) -> float | None:
    """Spearman's rho, or ``None`` where it is undefined (fewer than 3 or a constant)."""

    from scipy import stats

    if len(left) < 3 or len(set(left)) < 2 or len(set(right)) < 2:
        return None
    rho = float(stats.spearmanr(left, right).statistic)
    return rho if math.isfinite(rho) else None


def _ranking_structure(
    metadata: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
    decomposition: Sequence[Mapping[str, object]],
    baseline: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """The facts a reader needs before quoting a top-*k* from this ranking."""

    structure: dict[str, object] = {
        "n_scored": len(rows),
        "n_distinct_scores": metadata.get("n_distinct_scores"),
        "largest_tie_block": metadata.get("largest_tie_block"),
        "n_impossible_masked": metadata.get("n_impossible_masked"),
    }
    if "n_impossible_masked_reversed" in metadata:
        structure["n_impossible_masked_reversed"] = metadata[
            "n_impossible_masked_reversed"
        ]
    if decomposition:
        structure["n_amplified"] = sum(1 for row in decomposition if row["amplified"])
        structure["spearman_rTS_vs_mTS"] = _spearman(
            [float(row["rTS"]) for row in decomposition],  # type: ignore[arg-type]
            [float(row["mTS"]) for row in decomposition],  # type: ignore[arg-type]
        )
    solved = {
        str(row["target_id"]): int(row["rank"])  # type: ignore[call-overload]
        for row in baseline
        if row.get("status") in (None, "", "optimal")
        and math.isfinite(float(row["moma_score"]))  # type: ignore[arg-type]
    }
    if solved:
        paired = [
            (int(row["rank"]), solved[str(row["target_id"])])  # type: ignore[call-overload]
            for row in rows
            if str(row["target_id"]) in solved
        ]
        structure["spearman_rank_vs_moma_baseline"] = _spearman(
            [float(a) for a, _ in paired], [float(b) for _, b in paired]
        )
    return structure


def _epsilon_sweep(model, reference, direction, candidates, config, base_rows):
    """Re-rank at each swept epsilon and report how far each top candidate moves.

    Both source papers report this sensitivity. Each swept value is applied to every changed
    reaction, so the sweep tests the ranking's dependence on the threshold's scale.
    """

    from cmm.features.revert import revert_targets

    if not config.validation.epsilon_sweep:
        return ()
    base_rank = {row["target_id"]: row["rank"] for row in base_rows}
    out: list[Mapping[str, object]] = []
    for epsilon in config.validation.epsilon_sweep:
        ranking = revert_targets(
            model,
            None,
            reference,
            direction,
            targets=list(candidates),
            method=config.method,
            alpha=config.alpha,
            epsilon=epsilon,
            parameter_k=config.parameter_k,
            perturbation=config.perturbation,
        )
        for index, target in enumerate(ranking.sorted().targets):
            out.append(
                {
                    "epsilon": epsilon,
                    "target_id": target.target_id,
                    "score": float(target.score),
                    "rank": index + 1,
                    "rank_at_configured_epsilon": base_rank.get(target.target_id),
                }
            )
    return tuple(out)


def _export(
    result: TransformationWorkflowResult,
    genes: "pd.DataFrame",
    model,
) -> TransformationWorkflowResult:
    """Write the schema-v2 run bundle: one manifest, one role per artifact."""

    import json as _json
    import shutil
    from dataclasses import asdict, replace

    import pandas as pd
    from cobra.io import write_sbml_model

    from cmm.workflows._bundle import ArtifactRecord, _ArtifactWriter, _jsonable

    config = result.config
    if config.output_dir is None:  # guarded by the caller; kept as an invariant check
        raise TransformationWorkflowError("cannot export a run with no output_dir")
    root = Path(config.output_dir).resolve()
    if root.exists() and not root.is_dir():
        raise FileExistsError(f"output path exists and is not a directory: {root}")
    if root.exists() and any(root.iterdir()) and not config.overwrite:
        raise FileExistsError(
            f"output directory is not empty: {root}; choose a new directory or set "
            "overwrite=True"
        )
    # Capture input bytes before cleanup, including when an in-place rerun reads its own
    # archived model and expression tables.
    model_path = Path(config.model_path)
    model_bytes = model_path.read_bytes() if model_path.is_file() else None
    expression_inputs = {
        field: (Path(path).suffix or ".csv", Path(path).read_bytes())
        for field, path in (
            ("source_expression_path", config.source_expression_path),
            ("target_expression_path", config.target_expression_path),
        )
    }
    if config.overwrite:
        # Remove only workflow-owned paths, so permission to rerun is never permission to
        # delete a user's unrelated files that happen to share the directory.
        for relative in (
            *_STAGE_DIRECTORIES,
            "00_config.json",
            "00_provenance.json",
            "00_summary.json",
            "00_manifest.json",
        ):
            owned = root / relative
            if owned.is_symlink() or owned.is_file():
                owned.unlink()
            elif owned.is_dir():
                shutil.rmtree(owned)
    root.mkdir(parents=True, exist_ok=True)
    for directory in _STAGE_DIRECTORIES:
        (root / directory).mkdir(parents=True, exist_ok=True)

    writer = _ArtifactWriter(root, error_type=TransformationWorkflowError)

    model_relative = f"model/{model_path.name}"
    if model_bytes is not None:
        (root / model_relative).write_bytes(model_bytes)
    else:
        write_sbml_model(model, str(root / model_relative))
    writer.existing(
        model_relative,
        stage="model",
        role="model",
        media_type="application/sbml+xml",
    )

    expression_relatives: dict[str, str] = {}
    for input_field, (suffix, payload) in expression_inputs.items():
        role = input_field.removesuffix("_path")
        relative = f"inputs/{role}{suffix}"
        (root / relative).write_bytes(payload)
        expression_relatives[input_field] = relative
        writer.existing(
            relative,
            stage="inputs",
            role=role,
            media_type="text/tab-separated-values" if suffix == ".tsv" else "text/csv",
        )

    writer.csv(
        "01_preflight/preflight.csv",
        pd.DataFrame([asdict(record) for record in result.preflight]),
        stage="01_preflight",
        role="preflight",
    )
    writer.csv(
        "02_reference/source_reference_fluxes.csv",
        pd.DataFrame(
            {
                "reaction_id": list(result.reference_fluxes),
                "flux": list(result.reference_fluxes.values()),
            }
        ),
        stage="02_reference",
        role="source_reference_fluxes",
        method="imat_sampling_mean",
    )
    writer.csv(
        "02_reference/gene_expression_states.csv",
        pd.DataFrame(
            list(result.gene_states),
            columns=["gene_id", "state", "high_fraction", "low_fraction"],
        ),
        stage="02_reference",
        role="gene_expression_states",
        method="rmta_discretisation",
    )
    writer.csv(
        "02_reference/imat_reference_pattern.csv",
        pd.DataFrame(
            list(result.reference_pattern),
            columns=[
                "reaction_id",
                "expression_state",
                "activity",
                "imat_flux",
                "pattern_lower_bound",
                "pattern_upper_bound",
                "v_ref",
            ],
        ),
        stage="02_reference",
        role="imat_reference_pattern",
        method="imat",
    )
    writer.csv(
        "03_direction/gene_differential_expression.csv",
        genes.reset_index().rename(columns={"index": "gene_id"}),
        stage="03_direction",
        role="gene_differential_expression",
        method=config.direction.significance,
    )
    writer.csv(
        "03_direction/reaction_direction_map.csv",
        pd.DataFrame(
            {
                "reaction_id": list(result.direction),
                "direction": list(result.direction.values()),
            }
        ),
        stage="03_direction",
        role="reaction_direction_map",
    )
    writer.csv(
        "04_candidates/candidates.csv",
        pd.DataFrame({"target_id": list(result.candidates)}),
        stage="04_candidates",
        role="transformation_candidates",
    )
    writer.csv(
        "05_transformation/reaction_epsilons.csv",
        pd.DataFrame(
            list(result.reaction_epsilons),
            columns=[
                "reaction_id",
                "sample_std",
                "t_test_epsilon",
                "epsilon",
                "floored",
            ],
        ),
        stage="05_transformation",
        role="reaction_epsilons",
        method="sampled_t_test",
        status="complete" if result.reaction_epsilons else "skipped",
        reason=None
        if result.reaction_epsilons
        else f"a single epsilon of {config.epsilon} was used for every reaction",
    )
    writer.csv(
        "05_transformation/transformation_ranking.csv",
        pd.DataFrame(list(result.ranking)),
        stage="05_transformation",
        role="transformation_ranking",
        method=config.method,
    )
    writer.csv(
        "06_validation/moma_baseline.csv",
        pd.DataFrame(list(result.moma_baseline)),
        stage="06_validation",
        role="moma_baseline",
        method="moma",
        status="complete" if result.moma_baseline else "skipped",
        reason=None if result.moma_baseline else "run_moma_baseline was disabled",
    )
    writer.csv(
        "06_validation/score_decomposition.csv",
        pd.DataFrame(
            list(result.score_decomposition),
            columns=None
            if result.score_decomposition
            else [
                "target_id",
                "bTS",
                "mTS",
                "wTS",
                "rTS",
                "amplified",
                "rank_rTS",
                "rank_mTS",
                "rank_delta",
            ],
        ),
        stage="06_validation",
        role="score_decomposition",
        method="rmta_equation_9",
        status="complete" if result.score_decomposition else "skipped",
        reason=None
        if result.score_decomposition
        else (
            "MTA computes one score per candidate; the bTS/mTS/wTS decomposition exists "
            "only for rMTA"
            if config.method != "rmta"
            else "no candidate has finite bTS, mTS and wTS"
        ),
    )
    writer.csv(
        "06_validation/epsilon_sensitivity.csv",
        pd.DataFrame(list(result.epsilon_sweep)),
        stage="06_validation",
        role="epsilon_sensitivity",
        status="complete" if result.epsilon_sweep else "skipped",
        reason=None
        if result.epsilon_sweep
        else "no epsilon_sweep values were configured",
    )

    exported_config = cast(dict[str, object], _jsonable(asdict(config)))
    exported_config["model_path"] = model_relative
    exported_config.update(expression_relatives)
    exported_config["output_dir"] = "."
    writer.json(
        "00_config.json", exported_config, stage="root", role="workflow_configuration"
    )
    writer.json(
        "00_provenance.json",
        {
            **result.provenance,
            "model_path": model_relative,
            **{
                field.removesuffix("_path"): relative
                for field, relative in expression_relatives.items()
            },
        },
        stage="root",
        role="provenance",
    )
    writer.json("00_summary.json", result.summary(), stage="root", role="summary")

    _write_reproduction_scripts(
        writer,
        result,
        root=root,
        model_relative=model_relative,
        expression_relatives=expression_relatives,
    )

    manifest_record = ArtifactRecord(
        path="00_manifest.json",
        stage="root",
        role="authoritative_artifact_manifest",
        media_type="application/json",
    )
    records = (*writer.records, manifest_record)
    manifest = {
        "schema_version": 2,
        "workflow": "transformation_target_discovery",
        "status": "complete"
        if all(record.status == "complete" for record in writer.records)
        else "partial",
        "artifacts": {
            record.role: {
                key: value for key, value in asdict(record).items() if value is not None
            }
            for record in records
        },
    }
    (root / "00_manifest.json").write_text(
        _json.dumps(_jsonable(manifest), indent=2, sort_keys=True, ensure_ascii=False)
        + "\n",
        encoding="utf-8",
    )
    return replace(result, run_directory=root)
