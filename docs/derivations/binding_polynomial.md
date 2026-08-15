# The three-complex competitive binding polynomial

*Stage 0, Step 2. Every step shown.*

This is the term that decides whether the model can be bistable at all, so it is
written to be followable at a whiteboard rather than only in code. The
implementation in `src/binding.py` follows this document line for line, and the
test suite checks the two against each other.

---

## 1 · Species and complexes

**Free species** (what we solve for): `P` free PTF1A · `E` free E-protein ·
`I` free ID3 · `R` free RBPJL.
**Constant pool:** `J` = RBPJ, broadly expressed and not the bottleneck.

**Four complexes**, each written with its dissociation constant:

| complex | composition | expression |
|---|---|---|
| ID3·E | ID3 + E-protein | `[ID3·E] = I·E / K_IE` |
| ID3·P | ID3 + PTF1A | `[ID3·P] = I·P / K_IP` |
| `C_L` | PTF1A + E + **RBPJL** — adult, Notch-refractory | `C_L = P·E·R / K_L` |
| `C_J` | PTF1A + E + **RBPJ** — immature, NICD-vulnerable | `C_J = P·E·J / K_J` |

ID3 traps **both** partners of an obligate heterodimer — Dufresne 2010
(PMID 20830706), in AR4-2J. That is the two-target structure, and it is what
makes this sharper than titrating either alone.

Ternary complexes use a **single-step assembly convention**: `C_L = P·E·R/K_L`
rather than tracking dimer intermediates. Justified in decision 006 — the
intermediates are faster still and no `k_on` has ever been measured for any of
these interactions.

---

## 2 · Conservation laws, and the fact that makes this tractable

Four conserved pools. Write each total as free + everything it is bound into:

```
P_tot = P + [ID3·P] + C_L + C_J
E_tot = E + [ID3·E] + C_L + C_J
I_tot = I + [ID3·E] + [ID3·P]
R_tot = R + C_L
```

Now substitute the complex expressions. **The key observation: every species
appears exactly linearly in every complex it belongs to, so each equation
factorises.**

```
P_tot = P + I·P/K_IP + P·E·R/K_L + P·E·J/K_J
      = P·[ 1 + I/K_IP + E·( R/K_L + J/K_J ) ]

E_tot = E + I·E/K_IE + P·E·R/K_L + P·E·J/K_J
      = E·[ 1 + I/K_IE + P·( R/K_L + J/K_J ) ]

I_tot = I + I·E/K_IE + I·P/K_IP
      = I·[ 1 + E/K_IE + P/K_IP ]

R_tot = R + P·E·R/K_L
      = R·[ 1 + P·E/K_L ]
```

Each bracket is a **binding partition function** — the factor by which a free
species' total exceeds its free concentration.

Define the shared **complex-forming capacity**

```
Φ ≡ R/K_L + J/K_J
```

which appears in both the `P` and `E` equations. `Φ` is the only place the two
ternary complexes enter, which is worth noticing: PTF1-L and PTF1-J compete for
the *same* `P·E` dimer pool.

---

## 3 · Eliminating `R`

The `R_tot` equation is linear in `R` once `P·E` is known:

```
R = R_tot / ( 1 + P·E/K_L )
```

Substituting into `Φ` and simplifying:

```
Φ(P,E) = R_tot / ( K_L·(1 + P·E/K_L) ) + J/K_J
       = R_tot / ( K_L + P·E ) + J/K_J
```

`R` is now gone. **Four unknowns have become three.**

---

## 4 · The reduced system

Divide each conservation law by its partition function:

```
    P = P_tot / [ 1 + I/K_IP + E·Φ(P,E) ]        (I)
    E = E_tot / [ 1 + I/K_IE + P·Φ(P,E) ]       (II)
    I = I_tot / [ 1 + E/K_IE + P/K_IP ]         (III)
```

Three equations, three unknowns. This is the binding polynomial system.

### Why this form and not the residual form

The obvious implementation solves

```
    P + [ID3·P] + C_L + C_J − P_tot = 0
```

and that is what the first implementation did. **It is numerically wrong in the
regime we care about.** When binding is tight, free `P` falls to ~10⁻¹² while
`P_tot` is ~1, so `P` is recovered as the difference between two nearly equal
numbers — catastrophic cancellation. The solver could not reach better than
~10⁻⁶ relative accuracy and failed to converge outright at `K_d ≤ 0.01`, which
is *precisely the regime where the mechanism lives* (§6).

Form (I)–(III) computes small `P` **directly as a ratio**. No subtraction of
nearly-equal quantities occurs anywhere. Residuals in log space,

```
    f₁ = ln P − ln P_tot + ln(1 + I/K_IP + E·Φ)      and similarly f₂, f₃
```

are all O(1), well-scaled, and keep the unknowns positive automatically.

**Measured effect.** Newton on the log-ratio residuals, seeded by 15 cheap
fixed-point sweeps of (I)–(III):

| regime | function evals | max residual |
|---|---|---|
| `K_d = 1` (loose) | 7 | 0 |
| `K_d = 10⁻²` | 25 | 4.4×10⁻¹⁶ |
| `K_d = 10⁻³` | 30 | 4.4×10⁻¹⁶ |
| `K_d = 10⁻⁵` | 19 | 8.9×10⁻¹⁶ |
| `K_d = 10⁻⁸` | 21 | 0 |

Machine precision everywhere, including where the previous solver failed. The
1e-6 tolerance is retired; the tolerance is now 1e-12.

---

## 5 · Loose-binding limit — recovers the first-order sink

Take `K_IE, K_IP ≫` the protein totals, and complex formation weak (`E·Φ ≪ 1`).

From (III), the denominator → 1, so **`I ≈ I_tot`**: ID3 is not appreciably
consumed by binding. Then (I) and (II) decouple:

```
    P ≈ P_tot / ( 1 + I_tot/K_IP )
    E ≈ E_tot / ( 1 + I_tot/K_IE )
```

and therefore

```
    C_L = P·E·R/K_L  ∝  P_tot·E_tot / [ (1 + I_tot/K_IP)(1 + I_tot/K_IE) ]
```

At large `I_tot` both brackets grow linearly, so

```
    C_L ∝ I_tot⁻²        ⟹   d ln C_L / d ln I_tot → −2
```

**This is exactly the first-order sink.** Two multiplied linear taxes, log-log
slope −2, independent of everything else. Verified numerically: asymptotic slope
**−1.990**.

> **So T1 and T2 are the same model in the loose-binding limit.** Not similar —
> the same. This is the analytic explanation for the numerical result in
> `test_titration_is_ultrasensitive_only_in_the_tight_binding_regime`.

---

## 6 · Tight-binding limit — threshold-linear, and where the sharpness comes from

Take `K_IE, K_IP ≪` the totals. Now ID3 binds essentially stoichiometrically,
and (III) matters: ID3 **is** consumed.

Consider the single-target intuition first, then the two-target case.

**Below threshold** (`I_tot` < available `E_tot + P_tot`): almost every ID3
molecule is bound, so free `E` is depleted roughly one-for-one,

```
    E ≈ E_tot − (ID3 bound to E)         a LINEAR decline, not a hyperbolic one
```

**Above threshold** (`I_tot` exceeds what the targets can absorb): the targets
are exhausted, free `E` collapses toward

```
    E ≈ K_IE · E_tot / ( I_tot − E_tot )
```

**The crossover width is set by `K_d`, not by the totals.** The transition
happens over a range of `I_tot` of order `K_IE`, while its *position* is at
`I_tot ~ E_tot`. So the response is steep in proportion to how much smaller
`K_IE` is than `E_tot`.

The standard result for molecular titration is that the effective Hill
coefficient scales as the square root of that separation. Measured here, with
`n_eff ≡ max |d ln C_L / d ln I_tot|`:

| `E_tot/K_d` | `n_eff` | `√(E_tot/K_d)` | ratio |
|---|---|---|---|
| 10 | 4.25 | 3.16 | 1.34 |
| 10² | 13.45 | 10.0 | 1.34 |
| 10³ | 43.47 | 31.6 | 1.37 |
| 10⁴ | 133.26 | 100.0 | 1.33 |
| 10⁵ | 333.75 | 316.2 | 1.06 |

```
    n_eff  ≈  1.34 · √( E_tot / K_d )
```

Clean square-root scaling over four decades. (The drift at 10⁵ is the sampling
grid, not the physics.) **The prefactor above 1 is the two-target effect:**
`C_L ∝ P·E`, and *both* factors are being titrated, so both collapse together.

**No cooperativity was assumed anywhere.** There is no Hill coefficient in this
derivation — the exponents `n` in the enhancer terms are a separate matter
(decision 004). The sharpness here is purely stoichiometric. That is the whole
point of molecular titration as an ultrasensitivity generator, and it is why
this term can produce switching in a model whose Hill exponents are sampled low.

---

## 7 · The consequence for Stage 2

Combining §5 and §6:

```
    K_d ≫ totals   →   n_eff → 2,  identical to a first-order sink
    K_d ≪ totals   →   n_eff ≈ 1.34·√(totals/K_d),  unbounded as K_d falls
```

**The discriminating variable is the dimensionless ratio `K_d/E_tot`, not `K_d`
and `E_tot` separately.** It should therefore emerge as one of the dimensionless
groups under nondimensionalisation, and it is sampled log-uniformly across 4–5
decades deliberately bracketing the transition at `K_d ~ totals`.

Stage 2 must report **discrimination power as a function of that ratio**, not a
single Q-value comparison:

> *"T1 and T2 are distinguishable only when `K_d/E_tot < X`. Within the
> pre-registered plausible range, `Y`% of samples fall in that regime. The
> topology competition has power only there, and here is what it concludes."*

`K_IE` and `K_IP` are sampled **independently**. Langlands 1997 (PMID 9242638,
*J Biol Chem* 272(32):19785–93) reports all three Ids binding E-proteins with
**high affinity** but Id3 interacting **weakly** with all four class B MRFs —
and PTF1A is a class B bHLH. PTF1A itself was not among the factors tested, so
this is an extrapolation, but it argues against tying the two constants to a
shared prior. Justification for the prior widths: `prereg/` (see the ID3 `K_d`
prior document).

---

## 8 · What the code does

`src/binding.py::solve_binding` implements (I)–(III) exactly:

1. 15 fixed-point sweeps of (I)–(III) to land in the basin — cheap, no Jacobian.
2. Newton (`scipy.optimize.root`, `hybr`) on the log-ratio residuals `f₁,f₂,f₃`.
3. Recover `R` from §3, then all four complexes from §1.
4. **Every solve outcome is logged**, converged or not (see `CLAUDE.md`'s
   standing rule on convergence accounting). Silently dropping failures would
   deplete the sample set precisely in the tight-binding regime — biasing the
   Q-value comparison *against* the discriminating regime and producing a
   confident wrong negative.
