# 011 · mRNA pulse forcing — analytic, and why Stage 6 survives it

*Date:* 2026-08-15 · *Status:* accepted

## Question

§1.3 insists the payload is a **pulse of known shape** — uptake, translation
ramp, first-order decay — and that this is what makes "duration" a real quantity
rather than a free variable. §3.2 wrote the inputs as bare additive `u(t)` with
no shape at all. What shape, and does it cost states?

## Positions considered

**Position A — three explicit mRNA states.** `dm_i/dt = −delta_m·m_i + dosing`,
with protein production proportional to `m_i`. Mechanistically transparent, and
redosing is just repeated impulses on a state.

Against: three more states on a model already carrying a state-count schedule
risk, to represent a first-order decay whose analytic solution is exact. The
states buy no dynamics that the closed form lacks.

**Position B — analytic forcing.** `u_i(t)` computed directly as a difference of
exponentials, normalised so each pulse delivers exactly `dose`.

*Where they actually disagree:* only on state count — the trajectories are
identical, because a linear ODE forced by impulses *has* this closed form.

## Decision

**Position B.**

```
u_i(t) = dose · Σ_k [exp(−delta_m·(t−t_k)) − exp(−k_translate·(t−t_k))] · 1{t ≥ t_k} / area
area   = 1/delta_m − 1/k_translate            requires delta_m < k_translate
```

Parameterised by **(dose per pulse, redosing interval)** — which *is* Stage 5's
axis, so this decision also defines that axis in units a bench can execute.

**The normalisation is not cosmetic.** Each pulse integrates to exactly `dose`,
which is what lets the necessity analysis compare payload subsets at matched
**total** delivered dose. Without it, changing an mRNA half-life would silently
change how much protein was delivered, and every subset comparison would be
confounded by the thing it is trying to control for.

`delta_m < k_translate` is enforced with an exception rather than a comment: if
decay outruns translation the pulse inverts and becomes negative, which would
mean the payload *removes* protein — physically impossible, and the constraint
that inputs are non-negative is load-bearing for Stage 5's control problem.

### Stage 6 precondition — verified symbolically, not assumed

Stage 6's closed form requires the exogenous supply to enter as a **constant
vector field in state space**, so that `[g₁,g₂] = 0` and the entire ordering
effect comes from the nonlinearity of `f`. Checked with SymPy on a deliberately
nonlinear `f`:

- `[g₁, g₂] = 0` with constant directions ✓
- with a **time-varying coefficient** `u(t)`: `∂g/∂x = 0` still holds, so
  constancy in `x` survives and `[g₁(t), g₂(t)] = 0` ✓
- `[f_A, f_B] = J(x)·(g₁ − g₂)` exactly ✓ — confirming the ordering effect is
  carried entirely by the Jacobian, i.e. by the nonlinearity of `f`

**A prescribed time-varying coefficient multiplying a constant direction
preserves the property.** Stage 6's closed-form result survives the pulse
forcing.

**One caveat the algebra exposes, recorded so Stage 6 does not trip on it.** The
expansion `Δ(s) = s²·[f_A,f_B](x₀) + O(s³)` assumes **autonomous** vector
fields. A pulse makes the flow non-autonomous, so the closed form applies to
**constant-amplitude holds** — which is what `Φ_A^s` denotes anyway. The pulse
shape therefore belongs to Stage 6's *simulation* arm (the isodose phase sweep),
not to the closed-form arm. Both arms are in the plan; they must not be
conflated.

## What would reverse this

1. **If measured LNP uptake in Matrigel is not first-order** — pre-flight item 6
   in Part 7 — the difference-of-exponentials shape is wrong and should be
   replaced by the measured profile. The interface takes `u_i(t)`, so this is a
   one-function change by construction.
2. **If translation saturates at high dose** (ribosome limitation), the input
   stops being proportional to delivered mRNA, `g_i` stops being a constant
   direction, and **Stage 6's closed form genuinely dies** — it would then need
   the state-dependent version. This is the one failure mode that costs a stage.
3. **If redosing intervals shorter than the translation ramp are of interest**,
   pulse overlap makes the closed form awkward and Position A becomes simpler.
