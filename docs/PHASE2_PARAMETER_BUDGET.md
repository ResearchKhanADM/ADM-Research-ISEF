# Phase 2 · parameter budget for the minimal core

*Date:* 2026-08-15 · **Status: implemented in `src/core.py`.**
Target from v3 Part 2: **3–4 states, ~9–12 parameters, nondimensionalized.**
Delivered: **3 states, 12 free groups.** Previous system: 11 states, 61 parameters.

> ## ⚠ REVISED after writing the right-hand side
>
> The first version of this file, written before `src/core.py`, said **11 groups**.
> **The correct number is 12**, and the revision is recorded rather than quietly
> patched because the count is the gate.
>
> Two parameters were missed by counting on paper:
>
> 1. **`k_w`, the chromatin write half-max.** The pERK input scale had already
>    been spent defining the ID3 titration level, so a second ERK-driven process
>    needs its own half-max relative to that scale. It cannot be absorbed.
> 2. **`b_P`, the basal ignition gain** — and this one was not a counting slip but
>    a **missing mechanism**. See §4.1: without it the model predicts trametinib
>    alone can never revert anything, contradicting Collins 2014.
>
> `n_C` is then pinned at 1 on the structural argument in §4.2, which is what
> brings the count back inside target. **The margin is not comfortable — 12 of
> 12 — so any new term must displace an existing one.**

---

## 1 · Why profile likelihood is the deliverable, not a side-quest

The three parameters chosen for profiling map **one-to-one onto the three
headline results**:

| Profiled parameter | What it controls | Headline result |
|---|---|---|
| **`a_P`** autoregulatory gain | does the loop close | **R2 composition** |
| **`γ = δ_C/δ_P`** memory timescale | how long it holds | **R3 durability** |
| **`κ = K_d/E_tot`** binding regime | threshold sharpness | **R1 formulation** |

That is not a coincidence to bury in a methods section. It means **profile
likelihood is the uncertainty bar on each deliverable**, not an identifiability
exercise run alongside them. Write it that way in the figure captions: the
profile on `κ` *is* the error bar on the co-formulation gap.

It also answers the standing attack cleanly. *"None of your parameters are
measured"* — correct, and **each unmeasured parameter has a named consequence and
a computed interval**. That is a stronger position than a fitted point estimate
with no interval at all.

---

## 2 · States — 3, plus prescribed inputs

| Symbol | Meaning | Why a state |
|---|---|---|
| `P` | PTF1A activity | autoregulatory; the loop's driver |
| `R` | RBPJL | **no P-independent production — that zero is the bootstrap claim** |
| `C` | chromatin/memory at metaplasia loci | the slow variable; sets time-to-relapse, which *is* R3 |

**Inputs — prescribed in time, never states:** `erk(τ)` with a swept withdrawal
rebound profile (decision 002 amendment); `u_P(τ)`, `u_R(τ)` as analytic pulses
(decision 011).
**Algebraic — never integrated:** `ID3 = f(erk)`; `E_free` from the exact binding
solution (decision 006 amendment).

---

## 3 · The 12 groups

| # | Group | Form | Role |
|---|---|---|---|
| 1 | `a_P` | `α_P/(δ_P·K_P)` | autoregulatory gain — **bistability lives here.** Profiled → R2 |
| 2 | `b_P` | `β_ign/(δ_P·K_P)` | basal ignition gain (ERK-suppressed) — see §4.1 |
| 3 | `c_rep` | `K_C/θ_C` | strength of repression by memory |
| 4 | `n_P` | — | Hill exponent, **scanned 1–4** |
| 5 | `a_R` | `α_R/(δ_R·K_R)` | *Rbpjl* gain (constant RBPJ pool absorbed here) |
| 6 | `n_R` | — | Hill exponent, **scanned 1–4** |
| 7 | `ρ` | `δ_R/δ_P` | RBPJL vs PTF1A turnover |
| 8 | `γ` | `δ_C/δ_P` | **the durability knob.** Profiled → R3 |
| 9 | `α_C` | `α_C/(δ_C·θ_C)` | chromatin write gain |
| 10 | `k_w` | `κ_I·K_wC/E_tot` | ERK half-max for writing, in ID3 units |
| 11 | `ε` | `ε_C/(δ_C·θ_C)` | self-reinforcement — memory rather than filter |
| 12 | `κ` | `K_d/E_tot` | **the binding regime.** Profiled → R1 |

**Pinned, not free:** `n_C = 1` (§4.2) and the self-reinforcement exponent = 2
(same double-counting argument).
**Absorbed into scales:** `δ_P`, `K_P`, `K_R`, `θ_C`, `E_tot`, `κ_I`, `IC50`, and
the constant RBPJ pool.
**Not counted — they parameterize inputs, and are swept or measured, not fitted:**
the rebound profile (rise, overshoot amplitude, peak time — bench item 7) and the
pulse shape (dose, interval, translation and decay rates).

`κ = K_d/E_tot` surfacing as a single group is a **consistency check that
passed**: `prereg/id3_kd_prior_justification.md` stated that the
nondimensionalization should produce exactly this ratio and that the
nondimensionalization was wrong if it did not.

---

## 4 · What writing the right-hand side changed

### 4.1 · The two complexes are not one complex — and collapsing them broke the model

The first draft wrote *Rbpjl* production against the **PTF1-L** complex
(`P·E_free·R`). That makes RBPJL production require RBPJL, so **`R = 0` becomes
absorbing**: no amount of MEK inhibition can restore it, and the model predicts
**trametinib alone never reverts anything**. That contradicts Collins 2014 head-on
and would have destroyed Phase 7's trametinib-only positive control — the arm that
makes the one-shot experiment un-failable.

The fix is the biology: *Rbpjl* is driven by **PTF1-J** (PTF1A + E-protein +
**RBPJ**, broadly expressed and not the bottleneck, so its constant pool is
absorbed into `a_R` at no parameter cost), while the *Ptf1a* autoregulatory
enhancer needs **PTF1-L** (with **RBPJL**). That asymmetry is the developmental
handoff, and it is what makes the loop a loop.

Stated precisely, because the two are easy to conflate and only one is defensible:

> The bootstrap claim is **"nothing but PTF1A makes RBPJL"** — not **"nothing but
> RBPJL makes RBPJL"**. The second is stronger, the biology does not support it,
> and it passes every guard test.

**A second missing piece followed from the same check.** `b_P` — ERK-suppressed
basal ignition — is required for the same reason: without it `P = 0` is absorbing
too, and only the payload could ever restore anything. Its ERK sensitivity is tied
to the ID3 scale rather than given its own half-max, a declared simplification
since both are direct responses to the same "ERK is high" state and nothing
measured separates their EC50s. Below `b_P ≈ 0.4` in the placeholder set, the
model still says trametinib alone cannot revert — so this is a **regime the
literature rules out**, not merely a parameter value.

### 4.2 · `n_C = 1`, on a structural argument

`C`'s memory already comes from its own self-reinforcement term. Making the
repression of `P` by `C` cooperative **as well** encodes the same switching
physics twice, and the two mechanisms then cannot be distinguished — a small
identifiability trap of exactly the kind this whole reduction exists to avoid.

So `n_C = 1` is adopted **because of the double-counting argument**, not because
it empirically made no difference. Same reasoning pins the self-reinforcement
exponent at 2: it is the minimal value that lets the C-subsystem hold a memory at
all, and scanning both would parameterize one phenomenon twice.

### 4.3 · One `K_d`, with a range sweep — not a single-point check

Two independent reasons one suffices: Langlands 1997 gives the rank order (Id3
binds E-proteins tightly, class B factors weakly) and PTF1A is class B; and
mechanistically the E47 arm is load-bearing even in Dufresne, because **E47
carries nuclear import** — that is the pathway the mislocalization phenotype runs
through.

**But PTF1A was not in the Langlands panel, so this is an extrapolation, and a
single sensitivity point cannot tell you whether the extrapolation matters.** The
check therefore sweeps the **full range from the Langlands rank order (E47 tight,
PTF1A weak) to parity (both equal)**. Stable across that entire span → one `K_d`
is justified *and the reason is stated*. Moves anywhere in the span → the second
`K_d` returns, **and the location of the flip is itself the finding.** Recorded in
decision 006.

### 4.4 · The `n_eff` prefactor is 0.5, not 1.34 — and that one mattered

`n_eff` is the sharpness of the bootstrap threshold, and **Phase 3 convolves the
per-cell LNP dose distribution against it to produce the co-formulation gap —
the project's headline number.**

The derivation's `1.34·√(E_tot/K_d)` was measured on the deleted model's *ternary
complex* under *two-target* titration, where ID3 taxes PTF1A and E-protein through
two multiplied factors (which is also why its loose limit is 2). The core titrates
**one** target and takes the slope on `E_free`: measured prefactor **0.5**, loose
limit **1**.

**Carrying 1.34 across that change of mechanism would have overstated threshold
sharpness by ~2.7× and inflated R1 by the same factor** — silently, because the
model still runs and the figure still renders. Caught by
`tests/test_core.py::test_n_eff_agrees_with_the_measured_slope_of_the_exact_solution`,
which permanently ties the formula to the function it describes.

Generalisation worth keeping: **a constant measured on one observable of one
mechanism is not a property of "the model".** When the mechanism changes,
re-measure rather than re-cite.

**Open for Phase 3:** the quantity actually convolved is the *complex* threshold
`P·E_free·R`, whose sharpness compounds `n_P`/`n_R` on top of the binding term.
Phase 3 must define sharpness on the quantity it convolves and **measure it
there**, not reuse either constant.

### 4.5 · Still open

- **`ε` and `α_C` are the least-constrained pair.** If profile likelihood cannot
  separate them, lump them and report **11 groups**, naming which two were lumped
  and why. A lumped combination reported honestly is stronger than two numbers
  reported separately when the data cannot tell them apart.
- **`C` is erased passively**, at rate `γ`, with no active erasure by the acinar
  complex. That keeps time-to-relapse a property of `γ` alone, which is what makes
  `γ` readable as the durability knob. If Gate B or Phase 5 needs active erasure,
  add it deliberately and displace something.

---

## 5 · `A`, the observable — a Phase 0 item, and the sharpest unresolved thing here

`A` (acinar output) is **not** a state. Decision 007's warning is unresolved and
gets sharper under this core, not softer:

> Collins reports a **fraction of cells** determined by a staining threshold. `P`
> is a concentration. **A deterministic single-cell model cannot produce that
> number at all**, because there is no population and nothing to threshold.

Phase 0's *"define the reversal score"* must resolve it in three parts:

**(a) A minimal observation model.** Spread on one or two parameters, propagated,
then thresholded. **Not** a full stochastic layer.

**(b) The staining threshold as a sampled nuisance parameter**, never a fixed pick.

**(c) Phase 6 targets framed as *timing* comparisons wherever possible**, not
absolute-percentage matches. The KRAS-history prediction is fundamentally
*"long-history cells take longer"* — a time difference, far less sensitive to
where the cutoff sits. If the timing conclusion survives the threshold prior,
*"our result does not depend on the staining cutoff"* is a sentence worth having.

**Convenient overlap: the same per-cell spread machinery is what Phase 3 needs for
LNP dose heterogeneity. Build it once.**

---

## 6 · Held against v3's targets

| | Target | Delivered |
|---|---|---|
| States | 3–4 | **3** |
| Parameters | ~9–12 | **12 free groups** (10 fitted + 2 scanned) |
| Nondimensionalized | required | yes |
| `E_free` | exact, not linear | exact closed form; `κ` is its regime diagnostic |
| pERK | input with rebound profile | yes; profile parameters swept or measured, not fitted |
| Profile likelihood | 3 key parameters | `a_P` → R2, `γ` → R3, `κ` → R1 |

**11 states → 3. 61 parameters → 12.**
