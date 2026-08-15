# 006 · Two-target competitive titration, solved at equilibrium

*Date:* 2026-08-15 · *Status:* **accepted, and PROMOTED 2026-08-15** — see the
amendment below. This is now load-bearing for the v3 core and for Phase 3.

---

## AMENDMENT — this is no longer a topology discriminator; it is the model

**What changed.** The T1-vs-T2 topology competition this decision was written to
serve is cut (decision 012). For about a day that made this file look like dead
work. It is the opposite.

**v3's `E_free = E_total − k·ID3` is the tight-binding limit of the equilibrium
solved here.** When `Kd ≪` the protein totals, each ID3 molecule takes one E
molecule 1:1 and the exact solution reduces to exactly that linear subtraction.
So the linear form is not an independent modelling choice competing with this
one — it is **this decision's mechanism, in its ultrasensitive limit**, and this
file is its justification.

**Consequence 1 — ship the exact form, not the linear one with a floor.** The
linear expression goes negative once `k·ID3 > E_total`. That negativity is not a
numerical nuisance to be clipped; it is the approximation announcing it has left
its domain. A `max(0, ·)` hack hides the domain violation instead of fixing it.
Use the exact solution from `docs/derivations/binding_polynomial.md`, and state
in the code comment that it reduces to `E_total − k·ID3` in the tight regime,
with

```
n_eff ≈ 0.5·√(E_tot/Kd)           tight limit, MEASURED on the shipped 1:1 form
n_eff → 1                          loose limit
```

as the check on which regime a given sample is actually in. No invalid region, no
hack.

> **Correction, 2026-08-15 — the constant below is not transferable, and this is
> the amendment's most important line.** The table further down measures
> `n_eff ≈ 1.34·√(E_tot/Kd)` and `n_eff → 2`. Both are correct **for what they
> measured**: the log-log slope of the *ternary complex* `C_L` under *two-target*
> titration, where ID3 taxes PTF1A and E-protein through two multiplied factors —
> which is exactly why the loose limit is 2 rather than 1.
>
> The Phase 2 core titrates **one** target and takes the diagnostic on `E_free`
> itself. Measured against the shipped closed form: prefactor **0.5**, loose limit
> **1**. Carrying 1.34 across the change of mechanism overstates threshold
> sharpness by ~2.7×, and since Phase 3's co-formulation gap grows with sharpness,
> it would have inflated **R1 — the project's headline number** — by that factor,
> silently. Found on the first integration of `src/core.py`, not by reading.
>
> Generalisation worth keeping: **a constant measured on one observable of one
> mechanism is not a property of "the model".** When the mechanism changes,
> re-measure rather than re-cite.

**Consequence 2 — the ultrasensitivity result now feeds Phase 3, the headline.**
Phase 3 convolves the per-cell LNP dose distribution against the P–R bootstrap
threshold to get converted fraction. **How sharp that threshold is determines the
answer**, and `n_eff` is exactly that sharpness. A soft threshold and an
ultrasensitive one give different converted-fraction curves and therefore a
different co-formulation gap — which is the number the whole project reports. The
table in the original decision below is no longer a statement about model
selection; it is a statement about the shape of the headline result.

**Consequence 3 — the Kd prior is still decisive, for a new reason.** It no
longer decides whether a topology competition can discriminate. It decides how
sharp the threshold is, and therefore the size of R1. `prereg/id3_kd_prior_justification.md`
survives for that reason; its Stage 2 framing does not.

**What does NOT survive:** T1-vs-T2 as competing topologies, the Q-value
comparison, and `first_order_sequestration` as a rival mechanism to be beaten.
The first-order sink is now simply a wrong limit, not a candidate.

### One `K_d`, and how the extrapolation gets checked

The core titrates **one** target with a closed form, where the deleted model
sampled `K_IE` and `K_IP` independently. Two independent reasons that suffices:

1. **Rank order.** Langlands 1997 (PMID 9242638) found all three Ids bind
   E-proteins with high affinity while Id3 interacts *weakly* with class B
   factors. **PTF1A is class B**, so the E47 arm should dominate.
2. **Mechanism.** The E47 arm is load-bearing even in Dufresne 2010, because
   **E47 carries nuclear import of the PTF1 complex** — that is the pathway the
   mislocalization phenotype actually runs through.

**But PTF1A was not in the Langlands panel.** That is an extrapolation, and **a
single sensitivity point cannot tell you whether an extrapolation matters** — it
tells you about one point. So the check is a **range sweep across the full span
from the Langlands rank order (E47 tight, PTF1A weak) to parity (both equal)**:

- **stable across the entire span** → one `K_d` is justified, and the sweep is the
  reason, stated rather than asserted;
- **moves anywhere in the span** → the second `K_d` comes back, **and the location
  of the flip is itself the finding** — it says how much asymmetry the conclusion
  can tolerate, which is more informative than either answer alone.

`src/supplementary/binding.py` is retained to run exactly this check; it is the
only implementation of the two-target mechanism.

### The prefactor does not transfer — see the correction box above

`n_eff` for the core is `0.5/√κ`, measured. The 1.34 below belongs to the ternary
complex under two-target titration. Reusing it would inflate Phase 3's
co-formulation gap — **R1, the headline number** — by ~2.7×.

---

*Original decision, 2026-08-15, preserved below.*

## Question

ID3 traps both E47 and PTF1A (Dufresne 2010, PMID 20830706, in AR4-2J). How is
that represented — and is the equilibrium solved algebraically or are complex
formation kinetics integrated explicitly?

## Positions considered

**Position A — explicit mass-action ODEs for every complex.** Strictly
unreduced, which is what Step 1 asked for: you cannot eliminate what you have
not written. No equilibrium assumption to justify.

Against: it adds four fast states (`ID3·E`, `ID3·P`, `C_L`, `C_J`) whose
association rates are seconds-to-minutes against transcriptional dynamics of
hours-to-weeks. That is a stiffness ratio of 10⁴–10⁵ *inside the fast block
alone*, on top of the 10⁴–10⁵ the model already has. And **no `k_on` has ever
been measured for any ID3 interaction** (Part 8) — the panel additionally
established that `k_on` is structurally non-identifiable here, because binding
is quasi-steady-state on transcriptional timescales. Integrating rates nobody
can constrain, to resolve dynamics nobody can observe, buys nothing.

**Position B — one shared equilibrium, solved for free species.** Three
conserved pools, four complexes, solved simultaneously each RHS call.

*Where they actually disagree:* on whether "unreduced" means *no algebraic
elimination anywhere*, or *no elimination that destroys structure*.

## Decision

**Position B — equilibrium, one polynomial, both targets.**

```
E_tot = E_free + [ID3·E] + C_L + C_J
P_n   = P_free + [ID3·P] + C_L + C_J
I     = I_free + [ID3·E] + [ID3·P]
R     = R_free + C_L
```

**Why this elimination is safe when the acetylation one was not.** Eliminating
`H` collapsed `Hill(C_L, H)` into `Hill(C_L, C_L)` — the variable vanished
*algebraically* and took its own mechanism with it. Here the eliminated species
retain a nontrivial nonlinear dependence on the slow states; nothing collapses.
**Elimination is dangerous when it destroys structure, not when it is fast.**
That distinction is the whole content of the rule, and it is worth stating
plainly because "we eliminated a fast variable" sounds identical in both cases.

Ternary complexes use a **single-step assembly convention** — `C_L = P·E·R/Kd_L`
— rather than tracking dimer intermediates, on the same grounds: the
intermediates are faster still and their rates are unmeasured.

**Measured consequence, and it constrains Stage 2.** Effective Hill coefficient
`n_eff = max |d ln C_L / d ln I|`:

| regime | titration | first-order sink (T2) | ratio |
|---|---|---|---|
| `Kd = 1.0` (loose) | 2.07 | 1.94 | 1.07× |
| `Kd = 0.01` (tight) | 13.40 | 2.00 | 6.70× |
| `Kd = 0.001` (tight) | 41.66 | 2.00 | 20.83× |

The first-order sink sits at exactly 2.0 throughout — its analytic slope, since
ID3 taxes `P` and `E` through two multiplied factors. **Titration only beats it
when `Kd` << protein totals.** If Stage 2's box does not reach that regime, T1
and T2 are indistinguishable and the competition silently fails to discriminate.
No `Kd` has ever been measured for any ID3 interaction, so this prior is a free
choice that determines the answer, and it must be pre-registered deliberately
rather than defaulted.

The full algebra, every step, is Step 2 → `docs/derivations/binding_polynomial.md`.

## What would reverse this

1. **If the effective `Kd` turns out to depend on assembly order** — i.e. the
   single-step ternary convention is not a good approximation of sequential
   assembly — switch to explicit sequential binding with a stated intermediate.
   Testable by comparing the two conventions at fixed overall `Kd`.
2. **If the equilibrium solve dominates runtime** in the Stage 2 sweep (>~50% of
   wall clock), replace the numeric solve with Step 2's closed form, or fall
   back to explicit fast kinetics with a stiff solver.
3. **If a measured `k_on` appears for any ID3 pair** and the association
   timescale turns out to be comparable to transcription rather than far faster,
   the quasi-equilibrium assumption fails and Position A becomes correct.
