# 002 · `W` (phospho-MEK) as an 11th state, protected from elimination

*Date:* 2026-08-15
*Status:* accepted

## Question

§1.2 predicts that trametinib shows a longer relapse delay than PD325901 **at
matched pERK suppression**. §3.2 specified `K_eff = K·f_act(v)·f_cat(v)`.

Can that prediction be derived from that equation — and if not, what is the
minimum change, given that the obvious fix adds a *fast* state to a model whose
entire Stage 0 premise is that fast states get eliminated?

## Positions considered

**Position A — keep ten states; drop or weaken the prediction.**

The model's reduction logic is its main defence against the schedule risk that
kills projects like this: reduce to 5–6 slow states, then continuation is
tractable. Adding a state that is *fast by turnover* and then exempting it from
the reduction is exactly the kind of special case that makes a model look
gerrymandered. A reader who accepts "fast variables are eliminated, that is why
`H` was cut" will immediately ask why `W` is different, and "because otherwise my
prediction disappears" is a bad answer.

There is also a quantitative objection, and it is strong. MEK phosphorylation
turns over in **minutes**. The model runs on **hours to weeks**. Within minutes
of washout both drugs converge to the same `W`. Collins measures relapse at
**7 days**. A minutes-long divergence in a days-long process contributes almost
nothing: the impulse `∫ΔK_eff dt` is small, and slow states integrate impulse,
not peak. On its face the asymmetry is invisible, so the state buys nothing.

Cheaper and more honest: model the two drugs as differing in potency, drop the
withdrawal claim, and keep ten states.

**Position B — add `W`, and protect it.**

The prediction is not derivable from the static product, and this is provable
rather than arguable. At matched pERK the two drugs produce identical `K_eff(t)`.
The ten states see the drug *only* through `K_eff`. Identical input to identical
dynamics gives identical trajectories, so both drugs arrive at withdrawal in the
same state and relapse identically. **A static product has nowhere to store the
difference between the drugs.** Either a state is added or the prediction is
false as stated — those are the only options.

Position A's timescale objection is correct but is an argument about *magnitude*,
not about *existence*, and it is answerable — see the decision.

*Where they actually disagree:* both accept the algebra, and both accept that `W`
is fast by turnover. They disagree on whether a **transient** can have a
**permanent** consequence in this system.

## Decision

**Add `W`. Protect it from elimination. Restate the prediction conditionally.**

Three parts, and the third is what makes the first two legitimate.

**1. The structure.** `f_act` and `f_cat` act at different points, which is what
requires the state:

```
dW/dt  = k_on · RAF_drive(K_eff) · f_act(v)  −  k_off · W
K_eff  = K · W · f_cat(v)
```

**2. `RAF_drive` is strictly decreasing in `K_eff`.** Phospho-MEK accumulates
under a catalytic inhibitor *because* falling ERK relieves ERK-mediated negative
feedback on RAF, so RAF drive rises as ERK falls. **An increasing implementation
inverts the mechanism and predicts the opposite result while still running
cleanly** — the worst class of bug, because nothing fails. The sign is asserted
in code and covered by a unit test that constructs an inverted `RAF_drive` and
requires the test to fail.

**3. Why the timescale objection does not kill it — and what it costs.**
Position A is right that the divergence lasts minutes. It is wrong that this
makes it invisible, because **the overshoot does not need to persist. It needs to
be large enough to carry the slow state back across the separatrix.** Once
crossed, relapse proceeds on the slow manifold regardless of what `W` does
afterwards. In a bistable system a transient can have a permanent consequence;
that is what bistability *is*.

But this is a demanding condition, not a free one, and the honest form of the
claim says so:

> *"Withdrawal asymmetry occurs **if and only if** the transient phospho-MEK
> overshoot is sufficient to return the state across the separatrix. We report
> the fraction of the plausible parameter ensemble in which that condition is
> met, and the overshoot magnitude required."*

The deciding quantity is the **impulse** `∫ΔK_eff dt`, not the peak — so `τ_W =
1/k_off` is **sampled across minutes-to-hours**, not fixed. Pure phospho-MEK
turnover is minutes; the feedback-relief components (DUSP, SPRY) are
transcriptional and run to hours. Fixing `τ_W` would presuppose the answer to the
question the stage is asking. **If the condition is met in 5% of the ensemble,
the result is 5%**, and that is still a publishable, falsifiable statement with a
named mechanism — considerably better than an unconditional claim that happens to
be true for unexamined reasons.

**Why elimination is forbidden.** QSS on `W` gives `W_ss ∝
RAF_drive(K_eff)·f_act(v)/k_off`; substituting into `K_eff = K·W·f_cat(v)`
recovers a static (implicit) relation between `K_eff` and the drug — i.e. the
static product, and with it the proof that the prediction vanishes. **Eliminating
`W` is not an approximation of this model; it is a different model that makes a
different prediction.** This is the same failure mode that killed `H` in the
pre-panel draft, and it is being consciously refused here rather than walked into.

`W` is therefore exempt from the Stage 0 fast-variable sweep, flagged in code and
in `CLAUDE.md`, and the exemption must be argued in the writeup where a reader
will find it — otherwise it reads as an inconsistency, and Position A's
"gerrymandered" objection lands.

## What would reverse this

1. **The primary reversal condition, checked at Stage 1.** If the two-parameter
   bifurcation shows the **ADM attractor sits far from the separatrix relative to
   the achievable overshoot** — i.e. the impulse from any physically plausible
   `W` transient cannot return the state across it anywhere in the ensemble —
   then `W` earns nothing. Cut back to Position A: ten states, drug difference
   modelled as potency, withdrawal claim dropped. **Report that as a result**;
   it is a clean negative and it is cheap to state.
2. **Ensemble fraction near zero.** If the §1.2 condition is met in <1% of the
   ensemble, the prediction is technically true and practically useless. Keep the
   state only if it costs nothing downstream; do not put it on the board.
3. **`τ_W` measured and short.** If a direct measurement in AR42J puts phospho-MEK
   recovery firmly at the minutes end with no slower component, the impulse
   argument weakens and condition 1 should be re-checked before `W` is defended.

Process note: as with 001, both positions were argued directly rather than by
subagent panel, per this session's tool policy.
