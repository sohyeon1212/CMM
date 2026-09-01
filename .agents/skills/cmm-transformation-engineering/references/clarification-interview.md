# CMM transformation clarification interview

Use this protocol only when read-only inspection leaves a decision unresolved or the user
explicitly asks to be interviewed, challenged, or guided through the setup. It resolves inputs;
the canonical CMM workflow remains responsible for numerical analysis and reporting.
Write repository instructions and artifacts in English. Conduct the actual interview in the
user's language unless the user requests another language.

## Inspect before asking

Classify every field:

- **Resolved:** explicitly supplied and internally consistent.
- **Discoverable:** uniquely determined from the model, the data files, or the environment;
  inspect it and state the resolution instead of asking.
- **Unresolved:** more than one scientifically meaningful choice remains; ask the user.
- **Blocked:** no valid choice can be established; explain what is missing.

These are **facts — inspect, never ask**:

| Fact | Where it comes from |
|---|---|
| Model fingerprint, objective, exchange inventory, current bounds | the model file |
| Whether the solver supports MIQP | `cmm.core.solvers` |
| How many replicates each expression file carries | the file's column count |
| How many of the model's genes the expression covers | intersecting the two id sets |
| Whether the model grows under an already-known condition | one LP |
| Blocked, essential and coupled-set counts | the candidate stage |

Model selection, the direction of the transformation, and the intended biological condition are
user decisions. Model bounds are evidence about a file, not evidence of the intended
experiment. Never interpret "use defaults" as permission to choose medium, substrate uptake,
aeration, the reference-state estimator, or ε.

## Resolve decisions in dependency order

Ask only the earliest unresolved decision whose prerequisites are settled. The order is forced:
ε cannot be recommended before v_ref exists, and v_ref cannot be computed before the condition
and the estimator are settled.

### 1. Exact model bytes

When the user's reference does not identify one file. "iAF1260" names a reconstruction, not a
distribution — different SBML releases of the same reconstruction use different reaction ids.

### 2. Which file is source, which is target

**Always confirm this, even when the file names appear to answer it.** It is the one decision
no inspection can check and no downstream stage can detect.

Frame it as the scientific question, not as file assignment:

> Are you asking *which perturbation would produce* state B from state A, or *which
> perturbation would revert* state A back to state B? The first makes A the source; the second
> makes B the target.

### 3. One coherent condition

Medium, substrate exchange and uptake, oxygen exchange and uptake, and other changed bounds.
Prefer complete compatible condition choices over mixing unrelated bound questions. Offer
`model-as-loaded` only as an explicit choice, after showing the inspected bounds and saying
what they do and do not establish.

### 4. Reference-state estimator

E-Flux2 or LAD. Present this as a choice with no correct answer, and state the shared caveat:

> Yizhak et al. compute the source flux state with iMAT. CMM implements no iMAT, so neither
> option reproduces the published pipeline. iMAT places no objective on growth while E-Flux2 at
> full objective fraction maximises it, so the ranking is conditioned on whichever you pick. If
> you have a reference state computed elsewhere, it can be supplied through the Python API
> instead.

E-Flux2 scales bounds by expression and needs QP; LAD fits |flux| to expression targets and
needs only LP. Neither relaxes the MIQP requirement of the search itself.

### 5. MTA or rMTA

Ask when the user has not said. Give the cost and the difference in what is claimed:

> **MTA** solves one optimisation per candidate. **rMTA** solves three — best case, MOMA, and
> the direction reversed — and demotes candidates that score well whichever way they are
> pushed. rMTA is the more conservative ranking and costs roughly three times as much; on a
> genome-scale model that is hours rather than an hour.

### 6. Epsilon

**Never offer a bare default.** Compute `TransformationWorkflowConfig.suggest_epsilon(v_ref)`
and present its percentiles, because ε is a flux magnitude and only the model at hand can say
what a meaningful one is:

> ε is how far a reaction's flux must move before the move counts as a success. It is measured
> in this model's flux units: among the reactions required to change, the median |v_ref| here
> is {median} and the upper quartile is {p75}. Yizhak et al. derive ε per data set from a
> sampled reference distribution, which CMM cannot reproduce, so this is a choice you are
> making rather than a value being recovered.

Recommend a value, say why, and offer the sweep (question 7) as the way to show the choice did
not decide the answer.

### 7. Epsilon sweep

Off by default because it multiplies runtime by the number of values. Offer it whenever the
ranking will be reported or published:

> Both source papers report an ε sensitivity analysis. Since ε is chosen here rather than
> derived, reporting how the ranking moves across a few values is the honest substitute. It
> costs one full run per value.

### 8. Perturbation level, when it matters

Default to gene level: it is what a user acts on. Ask only when the user has said they want to
compare against a published reaction-level ranking, in which case reaction level with
coupled-set reduction is the comparable construction. Do not raise reproduction as a goal the
user has not stated. State the consequence either way — the candidate count is the denominator
of any percentile claim.

## Stop conditions

These are not questions. Explain what is missing and stop.

| Condition | Why it stops the run |
|---|---|
| No MIQP-capable solver | MTA and rMTA cannot run; `rmta_continuous` is a different method |
| Expression and model share no gene ids | almost always mismatched identifier systems |
| The model does not grow under the stated condition | every downstream flux state is meaningless |
| No reaction is labelled as changed | the threshold or the data, not the method |
| Both files carry one measurement per gene and the user wants the published protocol | the t-test needs replicates |

## Confirm before running

Present the resolved definition and obtain explicit confirmation. Include, in words:

- which file is source and which is target, and the question that ordering asks;
- the condition;
- the estimator, with the iMAT caveat restated;
- the method and its cost;
- ε, and that it was chosen rather than derived.

A prompt that was complete from the start requires neither an interview nor a confirmation
round — except for the source/target direction, which is confirmed on every run.
