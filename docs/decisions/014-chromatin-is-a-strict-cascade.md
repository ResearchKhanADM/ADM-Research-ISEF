# 014 · `C` is a strict cascade — the payload cannot reach the durability endpoint

*Date:* 2026-08-15 · *Status:* **RETIRED AS MOOT, same day, by decision 015.**
The finding below is correct — `C` really is a strict cascade — but the question
it poses does not arise, because relapse turns out **not to be chromatin-limited**
at all. Neither fix is adopted; no parameter is added; the budget stays at 12.
Read `015-relapse-is-not-chromatin-limited.md` first. Preserved unedited because
the cascade property is still true and still worth knowing.
*Panel:* adversarial-reviewer mandate, arguing against the proposed handling of
the (KRAS × trametinib) plane. It found this instead, which is more important
than the question it was asked.

**Read this instead of the transcript.** It determines whether Phase 5 has a
deliverable.

---

## Question

`dC/dτ = γ·(α_C·id3/(k_w + id3) − C + ε·Hill(C, 2))`

`C` is the slow chromatin/memory state, and **`C` is what sets time-to-relapse —
the R3 deliverable.** Nothing else in the model is slow.

Look at what `dC/dτ` depends on: the input, and `C` itself. **Not `P`. Not `R`.
Not the payload.** Verified directly: across nine `(P, R)` combinations at fixed
`C` and input, `dC/dτ` takes exactly **one** distinct value; with payload dosed at
`u_P = u_R = 10` it takes **the same** value.

So `C` is a strict cascade — the input drives it, and it drives the rest of the
model, but nothing flows back. **The payload has no channel to the durability
endpoint at all.**

## Why this is fatal to a phase, not merely untidy

Phase 5's deliverable is a **dose × redosing-interval map**, three regions,
undershoot / success / toxic, with time-to-relapse as the outcome. If relapse
timing is set by `C`, and `C` cannot see dose or interval, **that map is flat by
construction.** It would render as stripes, and it would be read — correctly — as
*"dose doesn't matter"*, which is not a finding about biology but an artefact of
a missing term.

`PHASE2_PARAMETER_BUDGET.md` §4.5 flagged active erasure as optional: *"If Gate B
or Phase 5 needs active erasure, add it deliberately and displace something."*
That was too casual. It is not an optional refinement; it is the only route by
which the intervention reaches the endpoint the project is named after.

## The two fixes, and why the choice is not mine

**Option A — add active erasure.** A term `−η·P·E_free·R` in `dC/dτ`: the restored
acinar complex actively erases metaplasia chromatin. Costs **one group, 12 → 13**
(recoverable to 12 by lumping `ε` and `α_C`, which §4.5 already licenses if
profile likelihood cannot separate them — though see decision 013, which removes
the mechanism that would have decided that).

- **For:** it is the mechanism most people would assume is already there, it
  restores a dose-dependent relapse time, and it makes Phase 5's map meaningful.
- **Against:** it is an *unmeasured* edge added to make a deliverable work, which
  is uncomfortably close to fitting the model to the figure. There is no
  measurement of acinar-complex-driven erasure at metaplasia loci in this system.

**Option B — keep the cascade, restate Phase 5's deliverable as reachability.**
The payload still matters, just not through `C`: it holds `P` and `R` up while
`C` decays passively, so what dose and schedule determine is **whether the acinar
state survives long enough to reach the far side of the memory**, not how fast the
memory clears. The deliverable becomes *"which schedules reach and hold the acinar
basin"* rather than *"how relapse time varies with dose"*.

- **For:** adds no unmeasured parameter, and is arguably the more honest reading
  of what a transient payload can do.
- **Against:** it changes what Phase 5 promises, and the three-region dose ×
  interval map — a named figure in v3 Part 5 — becomes a different figure.

**Why this is Tier 3 and not mine to settle.** Option A adds an unmeasured
mechanism to preserve a promised deliverable. Option B keeps the model honest and
changes what an endpoint means. That is a trade between scientific conservatism
and deliverable scope — **Luqmaan's call under Tier 3(a)**, and it constrains the
wet-lab arms under 3(b). Batched for this session's report.

**Recommendation: B, with A pre-registered as a stated extension.** Restate the
endpoint as reachability now, and register in advance that if the bench shows
chromatin marks clearing faster in reverted cells than passive decay predicts,
`η` enters with a measurement behind it rather than to rescue a figure. That
ordering — honest model first, mechanism added when evidence arrives — is the one
this project has taken every other time.

## What would reverse this

1. **Evidence of complex-driven erasure.** Any measurement that metaplasia-locus
   marks clear faster in cells that have re-expressed the acinar program than in
   cells that have not. That makes `η` a measured edge, and Option A becomes
   correct rather than convenient. Falvo 2023's retained H3K4me1 at metaplasia
   genes is the natural place to look.
2. **`C` turns out not to set relapse timing.** A separate panel finding (recorded
   in the session log) measured `d ln(time-to-relapse)/d ln γ ≈ 0.000` at high ERK
   — relapse there is a fast `(P, R)` collapse driven by free E-protein crashing,
   not a chromatin clock. **If relapse is generically not chromatin-limited, this
   entire decision is moot and `C`'s role in the model needs re-examining before
   either option is implemented.** Check this first — it is cheaper than both
   fixes and it could remove the question.
3. **Phase 5 is cut under compression.** It is third in the cut order (after
   Phase 6 and Phase 4's simplex interior). If it goes, so does the pressure.

## Downstream — what this changes

- **Phase 5's deliverable is in question until this is answered.** Do not build
  the dose × interval map until the endpoint is restated or `η` is added.
- **`ε = 0.5` makes `C` a filter, not a memory.** Bistability of the `C`
  subsystem needs `ε > 1.5396` (exact: the reciprocal of `max d/dC[C²/(1+C²)] =
  3√3/8`). The default ships below it, and `Params.eps` was commented "memory, not
  filter". Corrected in code, with the threshold exported as
  `core.EPS_MEMORY_THRESHOLD`. **A lagged filter still delays relapse, so this is
  not necessarily wrong — but the word "memory" must not appear on a poster
  without checking the value**, and any prior range for `ε` should straddle 1.54
  so the question is tested rather than assumed.
- **Gate C (the KRAS-history effect) is at structural risk from the same
  cascade.** History enters as `C(0)`, and `C(0)` is *forgotten* during a long
  drug hold — the panel measured identical outcomes for `C(0)` from 0.05 to 6.0
  at a hold of 600 scaled units. Gate C therefore survives only if the drug hold
  is short relative to `1/γ`. **That ratio depends on the unmeasured PTF1A
  half-life (bench item 9) and must be pre-registered before Phase 6**, or Gate C
  fails for structural reasons that have nothing to do with the biology it is
  testing.
