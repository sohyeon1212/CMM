# Architecture decision records

Decisions whose *why* is not recoverable from the code. Each record states what was decided,
what it costs, and what a reader must not "fix" without reopening the decision.

The first four are **recorded retrospectively**: the decisions were already in force and enforced by
[`AGENTS.md`](../../AGENTS.md), [`docs/VALIDATION.md`](../VALIDATION.md), and the scenario
documents. Writing them down moves the rationale out of prose and into a place an agent or a new
contributor can find by asking "why is this like that?". The dates are the dates of the record,
not of the original decision.

| ADR | Decision | Reopen if |
|---|---|---|
| [0001](0001-solver-capability-gate.md) | A missing solver capability fails the run; it never silently downgrades the method | a substitution becomes scientifically equivalent |
| [0002](0002-reference-state-and-epsilon.md) | *Superseded by 0005.* `v_ref` came from E-Flux2/LAD; ε had no default | — |
| [0003](0003-full-coupling-not-partial.md) | Coupled sets come from the null space of S (full coupling), not the paper's partial coupling | O(*n*²) LPs become affordable |
| [0004](0004-no-synthesized-recommendations.md) | The publication report presents each method separately and proposes no strain | a validated ranking across methods exists |
| [0005](0005-imat-reference-and-per-reaction-epsilon.md) | `v_ref` is the mean of samples from one iMAT activity pattern (rMTA discretisation), and ε is derived per reaction from the same samples | a published reference implementation shows a different GPR pairing or sampling protocol |

## Writing a new one

Copy the shape of an existing record: **Status · Context · Decision · Consequences · Source**.
Keep it under a page. A record explains a choice that a reasonable person would otherwise undo;
it is not a place for API documentation — that belongs in
[`docs/agent-reference.md`](../agent-reference.md).
