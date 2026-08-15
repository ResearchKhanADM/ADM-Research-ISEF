# 015 · Relapse is not chromatin-limited — what R3 actually is

*Date:* 2026-08-15 · *Status:* **DECIDED-PENDING-REVIEW** (Tier 2)
*Evidence:* `results/relapse_mechanism/` — three independent tests, run before
answering the question they were attached to.

**Read this instead of the transcript.** It retires decision 014's question,
**and it contradicts the R3 reframing Luqmaan proposed in the same message.**

---

## Question

Decision 014 found that `dC/dτ` depends only on the exogenous input — not on `P`,
`R`, or the payload — so the payload has no channel to the chromatin state. Two
fixes were on the table: add active erasure (one parameter), or restate Phase 5's
deliverable as reachability.

**Both assume `C` sets time-to-relapse.** That assumption was never tested. If it
is false, neither fix is needed and the question dissolves.

## What was measured

`scripts/run_relapse_mechanism.py`, at placeholder parameters. All times in
`1/δ_P`; **no conversion to hours is possible or attempted.**

**1 · Sensitivity of relapse timing to `γ`.** If `C` were the clock,
`t_relapse ∝ 1/γ` and the log-log slope would be −1.

| post-withdrawal ERK | `d ln t_relapse / d ln γ` |
|---|---|
| 0.75 | −0.253 |
| 0.90 | +0.066 |
| 1.20 | +0.002 |
| 2.00 | +0.001 |
| 3.00 | +0.001 |

**2 · Freeze `C`** (`dC/dτ := 0`) — removes the mechanism rather than perturbing it.

| ERK | free | frozen | ratio |
|---|---|---|---|
| 0.75 | 110.1 | **never relapses** | — |
| 0.90 | 4.6 | 4.8 | 1.036 |
| 1.20 – 3.00 | 2.1 | 2.1 | 1.000 |

**3 · Vary how much memory is written before withdrawal.** `C(0)` from 0.1 to 4.0:
relapse time **4.6, 4.6, 4.6, 4.6** at ERK 0.9 and **2.1, 2.1, 2.1, 2.1** at ERK
2.0 — identical to the digit.

Zero solver failures across all three.

## Decision

**Relapse is not chromatin-limited. Decision 014's question is retired, not
answered** — active erasure is not added, and Phase 5's endpoint is not restated
as reachability-through-`C`, because neither addresses what actually sets the
outcome.

**What sets durability is where the post-withdrawal ERK drive sits relative to the
upper fold.** Below it (ERK < 0.7123 at these parameters) the acinar state is an
attractor and reversion persists indefinitely. Above it, the state falls over, and
*how fast* is set by how far above the fold it sits — not by the payload, not by
the schedule, not by the memory.

**`C` is not useless, and the exception is precise.** In a narrow band immediately
above the fold — ERK ≈ 0.72 to 0.90 — freezing `C` prevents relapse entirely, and
`γ` has its only non-zero slope (−0.25). There, chromatin is what carries the
state back across. That band is ~25% wide in ERK and sits just outside the
persistence window. **`C` is a boundary-layer effect, not the clock.**

## ⚠ This contradicts the proposed R3 reframing

Luqmaan proposed, conditional on this check:

> *"R3 becomes: durability is set by drug-hold duration, not by payload
> composition."*

**The first half does not survive.** Drug-hold duration acts on relapse *only*
through `C` at withdrawal — and test 3 shows relapse timing is identical across a
40-fold range of `C` at withdrawal. **Drug-hold duration does not set durability
either.** The logic was right — "the only lever on `C` is how long the input is
held" — but it inherited the untested premise that `C` is the lever on durability.

**The second half survives and is strengthened.** Payload composition does not set
durability. What does:

> **R3 — Durability is a threshold property of the post-withdrawal KRAS drive,
> not a graded property of dose or schedule. The payload buys the *crossing*;
> whether it *sticks* is decided by whether the drive sits inside the persistence
> window, and nothing the payload or the schedule does moves that window.**

This is sharper than either version, and it is *more* falsifiable: it predicts
durability is **all-or-nothing in KRAS level and flat in payload dose**, which one
experiment can refute. A graded dose-response for durable reversal would kill it
outright.

It is also honest about what the model cannot promise. "Hold trametinib for N
before withdrawal" is not a protocol this model supports — and offering it would
have been a recommendation with no mechanism behind it.

## What would reverse this

1. **The narrow band widens.** The `C`-matters band is ERK ≈ 0.72–0.90 at
   placeholder parameters. If the prior-predictive sweep shows it covering a large
   fraction of plausible parameter space, `C` is not a boundary layer and this
   decision needs revisiting. **Check this in the Sobol/prior-predictive run** —
   it is one extra output, not a separate study.
2. **A bench result showing graded durability in dose.** If durable reversal
   fraction rises smoothly with payload dose rather than switching, the
   threshold claim is wrong and `C` — or something else slow — is doing more than
   this model says.
3. **`ε` above the memory threshold changes the picture.** The default `ε = 0.5`
   makes `C` a lagged filter, not a bistable memory (`EPS_MEMORY_THRESHOLD` =
   1.5396). **A genuinely bistable `C` could hold its own state after withdrawal
   and re-suppress the acinar program**, which is exactly the mechanism this
   decision says is absent. **This is the single most likely reverser and it is
   cheap to test** — re-run the three tests at `ε = 2.0`. Recorded here rather
   than run now because it belongs inside the pre-registered sweep over `ε`,
   whose range must straddle 1.54.

## Downstream

- **Decision 014 is retired as moot.** No new parameter; budget stays at 12.
- **v3's R3 wording changes**, and Part 5's `fig05_durability` — a dose × interval
  three-region map — is not a figure this model can produce. The persistence
  window (`fig02`) already carries R3.
- **Phase 5's ordering arms are unaffected** — they test the ERK→ID3 edge, not
  durability.
- **Every time axis stays in `1/δ_P`** until bench item 9 lands.
