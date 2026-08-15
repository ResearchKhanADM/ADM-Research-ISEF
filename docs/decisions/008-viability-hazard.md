# 008 · `dS/dt`, `γ(capacity)`, and the U-shaped hazard

*Date:* 2026-08-15 · *Status:* accepted

## Question

§3.2 said both *"death when `S` exceeds tolerance"* and *"the curve is
U-shaped."* Those are different objects. What is the viability axis actually
made of — given it is the y-axis of the headline reversal–viability figure?

## Positions considered

**Position A — threshold on `S`.** Simple, and it is what the plan literally
said. A trajectory either crosses `S_crit` or it does not.

Against, decisively: a one-sided threshold on a cargo/capacity ratio **is the
`U_crit` construct the expert panel killed**, in new clothing. The panel's
objection was that secretory capacity is co-induced with cargo by the same
differentiation program, so at steady state the ratio is invariant and the
ceiling does not exist — and in every published experiment where full
redifferentiation occurred, the cells survived. A threshold also cannot
represent the U-shape at all: it has no mechanism by which *too little* PTF1A
kills.

**Position B — an integrated hazard rising at both ends.**

*Where they actually disagree:* on whether viability is an **amplitude** limit
or a **rate** limit. The panel's finding was that the defensible version is a
kinetic mismatch during the transient — cargo ramps faster than capacity — which
is a rate statement, and a threshold cannot express it.

## Decision

```
dS/dt        = k_cargo·A − gamma(capacity)·S
h(S, P_free) = h_max_cargo·hill_activate(S, s_crit, nu_s)      ← cargo outruns capacity
             + h_max_chop ·hill_repress (P_free, p_crit, mu_p) ← CHOP under PTF1A loss
survival(t)  = exp(−∫₀ᵗ h dτ)
```

`h_low` is grounded: CHOP-dependent apoptosis under PTF1A loss (Sakikubo 2018,
PMID 30361559), and blocking dedifferentiation under stress *increases* death
(Backx 2021, PMID 33762742).

**The two arms are additive, not multiplicative.** The cell can die of either
cause independently; a product would make each arm require the other to be
non-zero, which would mean a cell with catastrophic cargo overload but healthy
PTF1A cannot die. That is wrong and it would quietly flatter the model.

`gamma(capacity) = gamma0 + gamma1·MIST1` in **T5**, and `gamma0` alone
elsewhere. Folding capacity into `gamma(A)` was proposed and **rejected**: it
makes T5 unimplementable, silently reduces Stage 2 from five topologies to four,
and rigs Stage 3's ranking of MIST1 as u₃ by making it undoseable.

**Three reasons this is not merely a smoothed threshold:**
1. It is genuinely **U-shaped**, which the plan asserts and Position A cannot express.
2. It **integrates along the trajectory**, so a brief excursion costs survival
   without automatically killing — the rate-limit statement, not an amplitude one.
3. It yields a **continuous 0–1 number**. A binary alive/dead flag cannot be
   plotted against a continuous reversal axis, and the headline figure needs
   exactly that.

Implemented as a `cumhaz` state so survival is read straight off the trajectory
rather than post-processed — which also keeps trajectories integrable to
completion instead of terminating at a stopping condition.

## What would reverse this

1. **If the reversal–viability plane comes out degenerate** — every intervention
   landing on one line, so viability adds no information beyond reversal — the
   hazard's two arms are not separating anything and the axis needs rethinking
   before it goes on a board.
2. **If HO-2 cannot be reproduced** (trametinib and KRAS-extinction at opposite
   corners), the hazard is too simple: it has no representation of *how many
   survival branches an intervention severs*, which is the stated mechanism of
   the headline claim. That would be the signal to add branch structure.
3. **If Stage 3 shows `h_max_chop` is completely unidentifiable** and the low-
   PTF1A arm never fires anywhere in the ensemble, the U-shape is decorative and
   should be reported as such rather than claimed.
