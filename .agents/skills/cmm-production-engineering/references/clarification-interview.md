# CMM clarification interview

Use this protocol only when read-only inspection leaves a decision unresolved or the user
explicitly asks to be interviewed, challenged, or guided through the setup. It resolves inputs;
the canonical CMM workflow remains responsible for numerical analysis and reporting.
Write repository instructions and artifacts in English. Conduct the actual interview in the
user's language unless the user requests another language.

## Inspect before asking

Build a small working run definition and classify every field:

- **Resolved:** explicitly supplied and internally consistent.
- **Discoverable:** uniquely determined from the model or environment; inspect it and state the
  resolution instead of asking.
- **Unresolved:** more than one scientifically meaningful choice remains; ask the user.
- **Blocked:** no valid choice can be established from available evidence; explain what is
  missing and recommend a safe next step.

Inspect model paths and bytes, exchange and objective inventories, current bounds, growth
feasibility where the condition is already known, solver capabilities, `straindesign`, and the R
renderer. These are facts. Do not ask the user for a fingerprint, package version, installed
capability, or uniquely matched exchange.

Model selection, product identity when mapping is not unique, and the intended biological
condition are user decisions. Model bounds are evidence about a file, not evidence of the
intended experiment. Never interpret "use defaults" as permission to choose medium, substrate
uptake, aeration, or other biologically decisive bounds. Offer `model-as-loaded` only as an
explicit choice after showing the relevant inspected bounds and their limitations.

## Resolve decisions in dependency order

Ask only the earliest unresolved decision whose prerequisites are settled:

1. Exact model bytes, when the user's reference does not identify one file.
2. Product exchange, when name or synonym matching yields zero or multiple plausible exchanges.
3. One coherent condition: medium, substrate exchange and uptake, oxygen exchange and uptake,
   and other changed bounds. Prefer complete compatible condition choices over mixing unrelated
   bound questions.
4. Capability-driven scope changes, only when the requested scientific method cannot run.
5. User-specific experimental constraints such as minimum retained growth, allowed intervention
   types, or edit budget, but only when the user raised them or the requested claim depends on
   them.

Do not turn deterministic protocol settings into interview questions. Unless the user has asked
to change them, use and disclose CMM's documented method suite, ranking sizes, seeds, sampling
settings, and reporting defaults in config and provenance.

After each answer, update only the affected downstream fields. Do not repeat settled questions.
Ask one dependency-bearing question at a time; at most a few independent, simple choices may be
grouped into one message.

## Question contract

Every decision question must contain:

1. A stable positive identifier such as `C1` or `S1` and one direct question.
2. The relevant facts already inspected.
3. Why the decision changes the calculation or scientific claim.
4. Two or three compatible choices whenever possible.
5. The recommended option listed first and explicitly labeled **Recommended**, followed by the
   evidence-based rationale and the consequence of accepting it.
6. A way to provide a custom value or say that the information is unknown.

Do not phrase the question so that agreement with the recommendation reverses a yes/no answer.
Do not recommend an option merely because it is convenient to compute. Base the recommendation
on the user's stated goal, inspected model evidence, and preservation of the requested method or
claim. When no biological choice is defensible, the recommended option is to pause and obtain
the missing condition or mapping.

Use this shape, adapting the choices to inspected evidence:

```text
C1 — Which biological condition should define this run?

Observed: the model contains <substrate and oxygen facts>, but the request does not specify the
experimental condition.
Why it matters: medium and aeration change growth, product yield, and every target ranking.

A. <coherent condition A> — Recommended. <evidence-based reason>. This supports <claim>.
B. <coherent condition B>. This changes <calculation or interpretation>.
C. Provide custom medium and bounds.

Recommendation: A, because <reason tied to the stated goal and inspected evidence>.
```

For a missing capability, preserve the requested claim when feasible:

```text
S1 — How should the unavailable solver capability be handled?

A. Install or select a solver that supports the requested method — Recommended. This preserves
   the full method and its intended claim.
B. Run an explicitly narrower workflow. This omits or substitutes <method> and cannot support
   <lost claim>.
C. Stop after preflight.
```

## Exit gate

When no unresolved decision remains, show one concise summary containing:

- exact model path and inspected fingerprint;
- product exchange reaction;
- medium, substrate and uptake, oxygen and uptake, and every changed bound;
- requested method scope, solver capabilities, and any approved substitution;
- user-specified growth, intervention, or edit constraints;
- output and publication-report scope.

If an interview occurred, ask the user to confirm this summary before launching the workflow or
writing run artifacts. A clear authorization in the user's final answer can serve as that
confirmation. If the original request was complete and consistent, skip both the interview and
this additional confirmation gate. If the user explicitly requested an interview despite a
complete run definition, audit it for consequential conflicts; when none exist, say so and offer
**Proceed with the supplied definition (Recommended)** with a concise rationale instead of
manufacturing technical questions.
