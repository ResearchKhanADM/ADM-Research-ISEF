# Phase 2 · parameter budget for the minimal core

*Date:* 2026-08-15 · **Status: proposed. Report before the RHS is written.**
Target from v3 Part 2: **3–4 states, ~9–12 parameters, nondimensionalized.**
Current implementation for comparison: **11 states, 61 parameters.**

This file exists because the parameter count *is* the identifiability budget. On
the old model it was discovered after the fact; here it is decided first, and the
right-hand side is written to fit it rather than the reverse.

---

## 1 · States — 3, plus two prescribed inputs

| Symbol | Meaning | Why it is a state |
|---|---|---|
| `P` | PTF1A activity | autoregulatory; the loop's driver |
| `R` | RBPJL | **no P-independent production term — that zero is the bootstrap claim** |
| `C` | chromatin/memory at metaplasia loci | the slow variable; it sets time-to-relapse, which *is* R3 |

**Inputs, prescribed functions of time, not states:**

- `k(t)` — pERK. Trametinib sets it while present; **on withdrawal it follows a
  rebound profile** (decision 002 amendment), not a step.
- `u_P(t)`, `u_R(t)` — delivered mRNA, as analytic pulses (decision 011 survives).

**Algebraic, not integrated:**

- `ID3 = f(k)` — a saturating function of pERK.
- `E_free` — from the exact binding equilibrium (decision 006 amendment), which
  reduces to `E_total − ID3` in the tight-binding limit.

---

## 2 · Dimensional parameters, before nondimensionalization

Grouped by the equation they enter. **19 dimensional parameters.**

**`dP/dt` — autoregulation, gated by complex availability, opposed by memory (6)**

| # | Parameter | Role |
|---|---|---|
| 1 | `α_P` | max autoregulatory production rate |
| 2 | `K_P` | half-max of the P·E·R complex term |
| 3 | `n_P` | Hill exponent, **sampled 1–4, never fixed** (decision 004) |
| 4 | `K_C` | half-max of repression by `C` |
| 5 | `n_C` | Hill exponent for that repression |
| 6 | `δ_P` | PTF1A turnover |

**`dR/dt` — driven only by P (4)**

| # | Parameter | Role |
|---|---|---|
| 7 | `α_R` | max *Rbpjl* production |
| 8 | `K_R` | half-max |
| 9 | `n_R` | Hill exponent |
| 10 | `δ_R` | RBPJL turnover |

*There is no 11th entry here, and its absence is the model's central claim.*

**`dC/dt` — slow memory, written by ERK, opposed by the acinar complex (5)**

| # | Parameter | Role |
|---|---|---|
| 11 | `α_C` | write rate |
| 12 | `K_wC` | ERK half-max for writing |
| 13 | `δ_C` | erasure rate — **with `δ_C` the slowest rate in the model, this is what sets time-to-relapse** |
| 14 | `ε_C` | self-reinforcement strength (memory rather than filter) |
| 15 | `θ_C` | self-reinforcement half-max |

**Titration and binding (3)**

| # | Parameter | Role |
|---|---|---|
| 16 | `E_tot` | total E-protein pool |
| 17 | `K_d` | ID3 dissociation constant (one, not two — see §4) |
| 18 | `κ_I` | gain from pERK to ID3 |

**Drug input (1)**

| # | Parameter | Role |
|---|---|---|
| 19 | `IC50` | trametinib → pERK |

*Not counted here, because they parameterize inputs rather than the system:* the
pERK rebound profile (rise, overshoot, settle) and the mRNA pulse shape
(dose, interval, translation and decay rates). They are **swept or measured**, and
the rebound profile is Bench Handshake item 7. Listing them as model parameters
would inflate the count with quantities the bench supplies.

---

## 3 · After nondimensionalization — **11 groups**

Scale time by `1/δ_P`, `P` by `K_P`, `R` by `K_R`, `C` by `θ_C`, ID3 by `E_tot`.
Six dimensional parameters collapse into scales:

| # | Group | Form | What it means |
|---|---|---|---|
| 1 | `a_P` | `α_P/(δ_P·K_P)` | autoregulatory gain — **bistability lives here** |
| 2 | `a_R` | `α_R/(δ_R·K_R)` | *Rbpjl* gain |
| 3 | `ρ` | `δ_R/δ_P` | RBPJL vs PTF1A turnover ratio |
| 4 | `γ` | `δ_C/δ_P` | **memory vs protein timescale — the durability knob** |
| 5 | `c` | `K_C/θ_C` | how strongly memory represses the acinar program |
| 6 | `ε` | `ε_C/δ_C` | self-reinforcement relative to erasure |
| 7 | `α_C'` | `α_C/(δ_C·θ_C)` | write gain |
| 8 | `n_P` | — | sampled 1–4 |
| 9 | `n_R` | — | sampled 1–4 |
| 10 | `n_C` | — | sampled 1–4 |
| 11 | `κ` | `K_d/E_tot` | **the binding regime — sets threshold sharpness, and therefore Phase 3** |

**11 groups. Inside the 9–12 target.**

Three of them are Hill exponents that are scanned rather than fitted, so the
effective fitting dimension is **8**. Profile likelihood runs on the three that
carry the results: **`a_P`** (does the loop close), **`γ`** (how long it holds),
**`κ`** (how sharp the threshold is).

`κ = K_d/E_tot` is the same dimensionless ratio the old derivation identified as
decisive, which is a consistency check on the nondimensionalization rather than a
coincidence — `prereg/id3_kd_prior_justification.md` predicted it would surface as
a group, and it does.

---

## 4 · Four judgement calls, flagged rather than buried

1. **One `K_d`, not two.** The old model sampled `K_IE` and `K_IP` independently on
   the strength of Langlands 1997's rank-ordering (Id3 binds E-proteins tightly,
   class B factors weakly). That asymmetry is real but it costs a parameter and it
   was only ever supported by a rank order in yeast against *different* partners.
   **Proposal: one `K_d`, with the asymmetry as a declared sensitivity check.**
   Reverses if the check moves the answer.
2. **`n_C` may be removable.** If repression by `C` is not sharp, `n_C = 1` and the
   count drops to 10. Test it before spending the parameter.
3. **`A` (acinar output) is not a state.** It was the validation observable in the
   old model. Under v3, reversal score is a Phase 0 decision and `P` is its natural
   proxy — but **decision 007's warning still stands**: Collins reports *% of cells*
   and `P` is a concentration. That map is a Phase 0 item, not a state.
4. **`ε_C` and `θ_C` are the least constrained pair here.** They exist to make `C` a
   memory rather than a filter. If profile likelihood cannot separate them, collapse
   the self-reinforcement to a single lumped term and report 10 groups.

---

## 5 · Held against v3's targets

| | Target | This proposal |
|---|---|---|
| States | 3–4 | **3** |
| Parameters | ~9–12 | **11 groups** (8 fitted + 3 scanned) |
| Nondimensionalized | required | yes, ~19 → 11 |
| `E_free` | exact, not linear | exact; `κ` is its regime diagnostic |
| pERK | input with rebound profile | yes — profile parameters are swept/measured, not fitted |
| Profile likelihood | 3 key parameters | `a_P`, `γ`, `κ` |

**Reduction: 11 states → 3, 61 parameters → 11.**
