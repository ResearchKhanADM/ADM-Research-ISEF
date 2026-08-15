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

## 1 · The parameter → deliverable mapping — **CLAIM WITHDRAWN, 2026-08-15**

> ## ⚠ THIS SECTION PREVIOUSLY ASSERTED A ONE-TO-ONE MAPPING. IT IS NOT TRUE.
>
> It read: *`a_P` → R2, `γ` → R3, `κ` → R1, one-to-one; therefore profile
> likelihood is the uncertainty bar on each deliverable; write it that way in the
> figure captions.* **Do not write that.** A methods panel measured it and the
> R1 leg is false in the model as shipped. Two separate problems, recorded
> together because they were found together — see `docs/decisions/013`.

**Problem 1 — the mapping is not one-to-one, and κ is not even dominant for R1.**
R1's co-formulation gap is set by the *bootstrap threshold in delivered dose*,
which depends on the whole complex `P·E_free·R` — not on the binding term alone.
A one-at-a-time diagnostic around the placeholder set found `a_R` (**not
profiled**) moving the gap far more than `κ` does over `κ`'s own declared range,
and `n_P` (**a scanned exponent, not a fitted parameter**) able to abolish the
conversion entirely, while `κ` is nearly flat across its first three decades.

The irony is instructive: **`core.py`'s own `n_eff` docstring and §4.4 below
already say the convolved quantity is `P·E_free·R` and that its sharpness
compounds `n_P`/`n_R`.** The 1.34→0.5 correction was found by exactly that
reasoning, and then this section was written contradicting it. *A claim can be
refuted by a file in the same repository.*

**Problem 2 — a parameter interval is not a deliverable interval anyway.** R1 is
a percentage gap, R2 a marginal-value curve, R3 a window in days. An interval on
`κ` becomes an interval on R1 only by pushing it through the Phase 3 computation,
a nonlinear map depending on the other eleven groups. Even a *correct* interval
on `κ` would not be R1's error bar without that propagation.

**What replaces it.** Global sensitivity analysis (Sobol) of each deliverable
against every group, which settles the mapping by measurement instead of
assertion — pre-registered, then run. Until it has run, the defensible phrasing is
**"the dominant single contributor to"**, and only where the sensitivity analysis
supports even that. See `prereg/2026-08-15_phase2_sensitivity_*`.

**And the profiling itself is separately misnamed** — there is no data, so there
is no likelihood. See §1a.

## 1a · "Profile likelihood" is not available, and the gate depends on it

Profile likelihood needs a likelihood — a density for **observed data** — plus
re-optimisation of the nuisance parameters at each fixed value of the profiled
one, plus a Wilks threshold. **This project has no data.** Building a Gaussian
penalty around two published timings with tolerances chosen by the modeller is
not a likelihood; it is the **same error as v3 §0.7's "stability selection"**
(a data-based method's name applied to a data-free operation on parameter space),
and it belongs in that table rather than in a caption.

Measured consequence, from the panel: the resulting "interval" is **the prior box
itself** for `a_P` and `γ` — 100% of the declared range — and where it is not,
**halving the invented tolerance changed it from 100% to 67% of the box**.
Nothing about the published observations changed. The width is a function of a
number the modeller typed.

**There is a real result in there, and it is the diagnosis rather than the
number:** a flat profile is the textbook signature of **structural
non-identifiability**. *"These parameters are not identifiable from the available
published observations"* is correct, correctly named, and useful.

**⚠ Gate B is currently unpassable as worded.** v3 defines it as *"key parameters
surviving profile likelihood"* — a gate defined in terms of a method that cannot
be run, which fails for reasons unrelated to the model. Rewording it is a **Tier 3
question** (it changes what an endpoint means) and is in this session's batch.

## 1b · The model has no clock — and this reaches Phase 5

Time is `τ = t·δ_P`. **`δ_P` — the PTF1A protein turnover rate — was absorbed by
the nondimensionalisation, is not one of the 12 groups, and is not on Bench
Handshake #1.** CLAUDE.md's own limitations say there is no measured PTF1A
half-life.

So any observation quoted in **days** requires reintroducing `δ_P` as a 13th
unknown, and because it multiplies every predicted time identically, *k*
day-valued timings supply only *k−1* dimensionless constraints. **Collins' two
timings therefore give one constraint — the ratio ~7/3 — and one of them is
consumed setting a clock nobody has measured.**

**This is not only a Phase 2 problem. Phase 5's deliverable is a redosing interval
in hours or days, and it has no units without `δ_P`.** v3 says bench items 3 and 4
give those axes units; they do not — mRNA half-life constrains a *product* of
dimensionless groups, not the protein clock. **Bench Handshake #1 gains item 9:
PTF1A protein half-life by cycloheximide chase.**

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
