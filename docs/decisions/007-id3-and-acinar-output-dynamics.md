# 007 · `dI/dt` and `dA/dt` — the two equations the plan implies but never writes

*Date:* 2026-08-15 · *Status:* accepted

## Question

Neither ID3 nor acinar output has an equation anywhere in §3.2, yet both are
load-bearing: ID3 is the titrator that decides bistability, and `A` is the
observable every validation target in Part 6 is expressed in.

## `dI/dt` — ID3

```
dI/dt = beta_I0 + beta_I·hill_activate(K_eff, k_I, n_I) − delta_I·I
```

**Positions considered.** *ID3 as a fast variable slaved to ERK* — §3.1 calls it
"intermediate", and if it is fast it can be substituted as `I ≈ I_ss(K_eff)`,
removing a state. Against: that hard-codes ID3 as a pure ERK readout with no
dynamics of its own, which forecloses the possibility that ID3's *relaxation
time* is what sets the lag between MEK inhibition and complex recovery — one of
the candidate explanations for the 3-day/3-week asymmetry that Stage 2 exists to
adjudicate.

**Decision.** Keep it as a state with a **sampled** `delta_I`, and let the
reduction step decide. It is the leading QSS candidate — eliminating it is how
the slow count reaches 6 — but that elimination must be *earned by timescale
separation measured in the fitted ensemble*, not assumed now. A basal term
`beta_I0` is included because ID3 is not zero in untreated acinar cells.

**What would reverse this.** If Stage 0's reduction shows `delta_I` is fast
across the whole plausible range, eliminate it and drop the state. If instead
`delta_I` lands in a stiff FIM direction, ID3 dynamics are carrying real weight
and the elimination must not happen at any stage.

## `dA/dt` — acinar output, and the observable map

```
dA/dt = k_A · output_gain · hill_activate(C_L, k_AC, n_AC) − delta_A·A
```

`output_gain` is 1 except in T6b, where NR5A2 acts here rather than on the
enhancer (decision 010).

**The part that matters more than the form.** Collins 2014 reports *"amylase-
positive acinar cells represented approximately 30% of the epithelial cells"* —
a **percentage of cells**, not a concentration. `A` is a continuous
concentration-like quantity. Comparing them directly is a category error, and it
would silently contaminate every validation target in Part 6.

**Decision:** `A` requires an explicit, stated map to *"% amylase-positive"*
before any comparison to published numbers. The map is a thresholded
single-cell readout, not a rescaling of the population mean — which is the same
lesson as PU.1/GATA1 (Stage 7): a population average and a fraction-of-cells are
different observables, and validating on the wrong one is how a beautiful
bistable model gets refuted later. The map is deferred to Step 3/4 where the
scaling is set, and is flagged there as a required deliverable.

**What would reverse this.** If the threshold that reproduces Collins' 30% is
not stable across the parameter ensemble — i.e. the reported percentage depends
mostly on where the threshold is put rather than on the dynamics — then `A` is
not a usable validation observable and the targets must be restated in terms of
something measurable in a dish.
