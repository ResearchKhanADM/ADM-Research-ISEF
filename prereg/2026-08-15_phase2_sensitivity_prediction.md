# Pre-registration · prediction — Phase 2 global sensitivity analysis

**Predictions written by Luqmaan, 2026-08-15, before the sweep ran.**
Ranges: `2026-08-15_phase2_sensitivity_ranges.yaml`.
**This file is committed and pushed before any sample is drawn.** The pushed
timestamp is the proof.

Method: Sobol variance decomposition (Saltelli sampling, SALib), first-order `S1`
and total-order `ST`, for each deliverable. Design detail in §"P5/P6" below —
**read that before running**, it changes what is computed.

---

## Why this runs

`docs/PHASE2_PARAMETER_BUDGET.md` §1 asserted that `a_P`, `γ` and `κ` map
one-to-one onto R2, R3 and R1. A one-at-a-time diagnostic contradicted the R1 leg.
**The claim is withdrawn.** This replaces assertion with measurement.

---

## The predictions

**P1 · Persistence-window width (fold separation).** First-order Sobol dominated by
`b_P` and `a_P`: **`S1(b_P) + S1(a_P) > 0.6`**.

**P2 · `γ` contributes `S1 < 0.05` to every durability output.** This follows
directly from the relapse-mechanism check (decision 015): freezing `C` changed
relapse timing by ≤3.6%, and a 40-fold change in written memory changed it not at
all.

> **⚠ STOP CONDITION. If `γ` appears anywhere with `S1 > 0.15`, that is
> INCONSISTENT with test 3 of decision 015. Treat it as a BUG HYPOTHESIS FIRST, a
> finding second — stop and report before running anything else.**

**P3 · `ε` shows a STEP, not a slope** — total-order sensitivity near zero when
sampled below 1.54, and **`ST > 0.3` above it**.

> *"This is the prediction I most want tested. If `ε` is flat across the whole
> registered range, `C` never becomes bistable in the plausible box and R3 is
> unconditional. If it steps, R3 is conditional and the conditional must go on the
> poster."*

**P4 · `a_P` and `b_P` interact:** `ST − S1 > 0.1` for each, because they jointly
set loop closure. First-order indices alone will under-attribute both.

**P5 · For the co-formulation gap: the `n_P` scan and `a_R` outweigh `κ`.**
**Declared WEAK — a consistency check, not a test**, because it was already
measured by a one-at-a-time diagnostic. **It must not be counted as a hit.**
*(Operationalisation changed from the original wording — see below.)*

**P6 · Hill exponents are scanned, not sampled, so they must not appear in the
Sobol indices at all.** If they do, it is an implementation error.

---

## ⚠ P5 and P6 contradict each other — resolved here, not silently

**The conflict.** P5 as originally written predicted `ST(n_P) > ST(κ)`, which
requires `n_P` to *have* a Sobol index. P6 says the Hill exponents must not appear
in the indices at all. Both cannot hold: computing `ST(n_P)` means sampling `n_P`
as an uncertain input, which is exactly what P6 forbids.

**Resolution — a stratified design, which honours P6 and preserves P5's content.**

- Sobol runs over the **10 continuous groups only**. `n_P` and `n_R` never enter
  the index computation, so **P6 holds by construction** and remains a live check
  on the implementation: if either appears, it is a bug.
- The sweep is **repeated at each point of the `(n_P, n_R)` scan**. The exponents
  are then a stratum label, not an uncertain input — which is what "scanned, not
  sampled" means operationally.
- **P5 becomes a between-strata comparison:** does the co-formulation gap move
  more *across the `n_P` strata* than `ST(κ)` accounts for *within* a stratum?
  That preserves what P5 was actually claiming — the scanned exponent matters more
  than the binding regime — while keeping the exponents out of the indices.

**Why this is the right resolution and not a dodge.** Averaging over a parameter
inside a variance decomposition asserts it is uncertain. `n_P` is not uncertain in
that sense; it is unknown-and-scanned, and decision 004 fixed that deliberately
("we did not use a sequence model to get a Hill coefficient, we scan it").
Sampling it would quietly convert a declared scan into a prior — the same class of
move as calling a tolerance an uncertainty.

**Cost:** the sweep runs once per stratum. With a 3×3 scan that is 9× the samples.
Acceptable; it is cheap relative to what a mislabelled index would cost.

**This changed Luqmaan's wording, so it is flagged rather than absorbed.** If the
between-strata form is not what P5 meant, say so and it will be rerun.

---

## What counts as a hit

- **P1, P2, P3, P4** — real predictions, counted.
- **P5** — consistency check, **not counted**.
- **P6** — an implementation check, **not counted** as a scientific prediction.

Reporting rule: report **all six outcomes**, hit or miss, including P2's stop
condition if it fires. A pre-registration that only reports its hits is not one.

## Convergence — the standing rule

Failures logged individually and reported **as a function of every swept
parameter, binned, never as a scalar**. Bound 0.05.

Two specific warnings from prior diagnostics, both already observed:

- failure rate swung ~13× across `κ` quintiles and was flat across `ρ`;
- failures correlated hard with `γ` — 25.8% in the slowest quintile vs 0.0% in the
  fastest — traced to a **fixed integration horizon while `1/γ` reaches 500**.
  **The horizon must scale with `1/γ`.** A fixed horizon manufactures failures in
  exactly the slow-memory regime P3 exists to probe.

If failures correlate with any swept parameter, **every result is reported as
conditional on convergence and labelled so on the figure.**

## Sign-off

- [x] Ranges reviewed and confirmed
- [x] Predictions written by Luqmaan
- [x] P5/P6 conflict flagged and resolved in writing, before running
- [ ] Committed **and pushed** before the sweep runs
