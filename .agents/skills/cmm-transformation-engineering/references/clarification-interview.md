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
| How many replicates each expression file carries | the file's column count; check the column names and drop any identifier column (e.g. Entrez) before the run |
| Whether the files are linear or log2 (`expression_scale`) | the values: negatives or a ~0–20 range mean log2 |
| How many reactions the direction stage labels changed | a dry run of the direction stage |
| How many of the model's genes the expression covers | intersecting the two id sets |
| Whether the model grows under an already-known condition | one LP |
| Blocked, essential and coupled-set counts | the candidate stage |

Model selection, the direction of the transformation, and the intended biological condition are
user decisions. Model bounds are evidence about a file, not evidence of the intended
experiment. Never interpret "use defaults" as permission to choose medium, substrate uptake,
aeration, or the discretisation threshold `imat_delta`. Per-reaction ε is a deterministic protocol
setting, not a question.

## Resolve decisions in dependency order

Ask only the earliest unresolved decision whose prerequisites are settled. The order is forced:
v_ref cannot be computed before the condition and `imat_delta` are settled, and the changed
set cannot be counted before v_ref exists.

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

### 4. Discretisation threshold `imat_delta`

v_ref comes from iMAT, which needs each source gene classed high, moderate or low. Valcárcel et
al. class a gene high above mean + δ·SD and low below mean − δ·SD within each sample, keeping a
class only when every sample agrees (`imat_beta` = 100%). They chose δ per data set, so it has
no default:

> δ sets how many genes iMAT treats as high or low. Valcárcel et al. used 0.3 for their
> Alzheimer's data and 0.7 for prostate cancer, where fewer samples made them stricter; their
> knockout experiments used 0.7–1. A smaller δ classes more genes and constrains the reference
> state harder. Recommended: {0.3 with many replicates per state / 0.7 with few}.

Report the resulting high/low gene counts (`gene_states` in the run summary) back to the user.
If iMAT matches no reaction at all, the δ is too wide for the data.

### 5. MTA or rMTA

Ask when the user has not said. Give the cost and the difference in what is claimed:

> **MTA** solves one optimisation per candidate. **rMTA** solves three — best case, MOMA, and
> the direction reversed — and demotes candidates that score well whichever way they are
> pushed. rMTA is the more conservative ranking and costs roughly three times as much; on a
> genome-scale model that is hours rather than an hour.

### 6. Changed-set size — only when it falls outside 100–200

Do not ask about the significance filter up front; the defaults (t-test, raw P < 0.05, top 200
by P) are Yizhak et al.'s protocol. Once v_ref exists, count the changed reactions the
direction stage produces. Running that stage alone is cheap next to the MIQP search. Ask only
when the count is outside 100–200. Give the count and a recommendation:

> With the current filter, {n} reactions are labelled changed. Both source papers aim for
> 100–200: Yizhak et al. keep the top 100–200 after a P < 0.05 t-test, and Valcárcel et al.
> loosened or tightened their filter per data set until they got there. Fewer than 100 leaves
> the MIQP with little to transform; many more slows a genome-scale run. Recommended: {e.g.
> relax the cutoff to P < 0.1 / keep the raw P instead of FDR / lower top_n_changed}.

FDR adjustment never reorders genes. It matters only when fewer genes pass than the top-N cut
would keep. With many replicates (roughly ten or more per state), recommend a fold-change floor
as well, `direction.min_abs_log2fc=0.263` (Valcárcel et al.'s 20%), because a P-ranked set
otherwise fills with changes too small to matter. Accept the count as-is if the user prefers, and state it in the report.

**Numeric ε, only on request.** Never offer a bare default. Compute
`TransformationWorkflowConfig.suggest_epsilon(v_ref)` and present its percentiles, because a
single ε is a flux magnitude that only the model at hand can calibrate.

### 7. Epsilon sweep

Off by default because it multiplies runtime by the number of values. Offer it whenever the
ranking will be reported or published:

> Both source papers report an ε sensitivity analysis. Reporting how the ranking moves across a
> few scalar values shows whether the threshold decided the answer. It costs one full run per
> value.

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
- `imat_delta` (and `imat_beta` if not 100%), with the three stated differences from the
  papers' reference state;
- the method and its cost;
- the expression scale that was inspected;
- ε: per reaction from the 2,000 samples behind v_ref (default), or the chosen number;
- the changed-set count, if question 6 was asked.

A prompt that was complete from the start requires neither an interview nor a confirmation
round — except for the source/target direction, which is confirmed on every run.
