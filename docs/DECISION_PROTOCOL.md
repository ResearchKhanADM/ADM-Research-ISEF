# Decision protocol — what Claude decides alone, what gets a panel, what stops

*Established 2026-08-15, session 6, by Luqmaan.*
Operative summary lives in `CLAUDE.md`; this is the reference.

**Why it exists.** Stopping to ask about decisions Claude is qualified to make was
costing more than it protected — a round trip per functional form. The fix is not
"ask less"; it is **a stated boundary**, so that what is decided alone is decided
alone *on the record*, and what stops still stops.

---

## TIER 1 — decide alone, log one line in `PROGRESS.md`, keep going

Implementation choices · functional forms where a default is defensible · naming ·
test design · refactors that keep the tests green · anything reversible in under a
day.

**Do not ask. Do not batch them — log each as one line as it happens.** A Tier 1
log line is `name — what was chosen — why, in a clause`. If it takes a paragraph,
it was not Tier 1.

## TIER 2 — convene a panel, decide, write the file, keep going

Anything architectural:

- a state-space change
- a new or removed parameter
- a changed functional form **that alters a reported number**
- a method substitution
- a conflict between two documents

**Procedure.** Spawn **two subagents with opposing mandates**, make them argue,
then decide. Write `docs/decisions/NNN-name.md` under the four headings, with
*what would reverse this* filled in **properly** — checkable conditions, not
sentiments. Flag it in the session report as **DECIDED-PENDING-REVIEW** so Luqmaan
reads the file rather than the transcript.

**The panel** lives in `.claude/agents/`:

| Agent | Mandate |
|---|---|
| `adversarial-reviewer` | assumes the proposal is wrong, argues for rejection, finds the failure mode rather than balancing |
| `methods-checker` | statistical and numerical validity; is the named method the method; does it do what it is claimed to do |
| `literature-verifier` | every empirical claim against PubMed; "no evidence found" stated explicitly, never inferred; **re-derives constants rather than inheriting them** |

Two is the minimum and they must genuinely oppose. Pick the second by what is at
stake: a method substitution wants `methods-checker`; a claim about biology wants
`literature-verifier`; a structural change wants `adversarial-reviewer` against a
case-for.

`literature-verifier`'s re-derivation rule is standing, not optional. It exists
because `n_eff = 1.34` was inherited across a change of mechanism and would have
inflated the headline number ~2.7× while everything still ran.

## TIER 3 — STOP AND ASK. Never automated.

**(a) Anything encoding Luqmaan's scientific priority rather than a technical
fact** — what he will trade off, what risk he accepts, what an endpoint means.
**(b) Anything that consumes or constrains the one wet-lab shot.**
**(c) External facts that cannot be obtained here** — the arm budget, bench
measurements, anything requiring him to talk to a person.
**(d) Anything where being wrong is expensive *and* irreversible.**
**(e) Any claim of novelty or priority.** One has already been wrong, and it would
have been made in public.

### Batching is part of the rule

**Do not stop at the first Tier 3 question.** Carry on with everything unblocked,
accumulate them, and present them together at session end as a **numbered list,
each with a recommendation and the cost of each option**. Five questions answered
once beats one question answered five times.

**If a Tier 3 question blocks everything**, stop — but say **explicitly that it is
a hard block**, and say what was tried first.

---

## The decisions folder is the interface

Luqmaan reads decision files, not transcripts. **Every Tier 2 file must stand
alone.** It carries what was at stake, what the options were, what was chosen,
what would reverse it, and **what it changes downstream**.

> **If a file cannot be understood without the conversation, it is not finished.**

This also serves the harder requirement: he must be able to explain **every line**
of this project to a judge. The decisions folder is how that happens without him
watching every keystroke. **Write for that reader** — the one who has the repo and
no memory of how it was built.

## Session end — always, in this order

1. **What was done.**
2. **Tier 2 decisions and why** — marked DECIDED-PENDING-REVIEW, with file links.
3. **Batched Tier 3 questions** — numbered, each with a recommendation and costs.
4. **What broke.**
5. **What is next.**
6. **Anything that changed a previously reported number.**

**Item 6 is not optional.** The `n_eff` correction had already propagated into
five documents before it was caught; a number that moves silently is worse than
one that was never reported.
