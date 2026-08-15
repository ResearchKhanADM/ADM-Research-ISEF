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
n_eff ≈ 1.34·√(E_tot/Kd)          derivation §6 — the validity diagnostic
n_eff → 2                          derivation §5 — the loose limit
```

as the check on which regime a given sample is actually in. No invalid region, no
hack.

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
