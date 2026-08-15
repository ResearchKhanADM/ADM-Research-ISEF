# 006 · Two-target competitive titration, solved at equilibrium

*Date:* 2026-08-15 · *Status:* accepted

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
