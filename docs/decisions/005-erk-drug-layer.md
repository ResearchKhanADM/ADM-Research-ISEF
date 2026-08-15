# 005 · `g(K_eff)`, `f_cat`, `f_act` — the ERK/drug layer

*Date:* 2026-08-15 · *Status:* accepted
*Companion to 002, which decides `W` and `RAF_drive`'s sign.*

## Question

What forms do the ignition promoter `g` and the two drug terms take, and — the
part that actually bites — **what argument does `g` receive?** §3.2 wrote
`g(K,v)` while also defining `K_eff` as carrying the drug.

## Positions considered

**Position A — keep `g(K, v)` as written.** It is what the plan says, and it
makes the drug's action on ignition explicit at the point of use rather than
hidden inside `K_eff`.

**Position B — `g(K_eff)` only.** `K_eff` already contains the drug through `W`
and `f_cat`. Passing `v` again means trametinib acts on ignition **twice**: once
by lowering `K_eff`, and again through `g`'s second argument. The model still
runs; it merely responds to the drug about twice as strongly as it should, and
every downstream dose recommendation is wrong by a factor nobody can see.

*Where they actually disagree:* on whether explicitness at the call site is
worth the risk of double-counting. It is not, because the double-count is
invisible and the explicitness is recoverable by a comment.

## Decision

**`g(K_eff)`, one argument.** `K_eff` is the single channel through which the
drug reaches every downstream term, and **no other term may take `v` directly.**
Enforced by convention and stated at the top of `functional_forms.g_ignition`.

Forms, all fractional-residual so that the no-drug case is exactly 1:

```
g(K_eff)  = 1/(1 + (K_eff/k_ign)^q)      DECREASING — cut #1, ERK shuts ignition down
f_cat(v)  = 1/(1 + v/ic50_cat)           MEK catalytic inhibition, both drugs
f_act(v)  = 1/(1 + v/ic50_act)           RAF→MEK phosphorylation, trametinib only
          ≡ 1                            when is_trametinib=False → recovers PD325901
```

The `f_act ≡ 1` branch is not a convenience. It is what preserves the direct
comparison to Collins 2014, who used PD325901 — and it is tested, because a
silently-diverging PD325901 arm would invalidate every validation target in
Part 6 at once.

Hill-type saturation rather than exponential decay for `g`: the ignition
promoter is a transcriptional response to ERK, and a saturating form is the
standard minimal choice. `q` is sampled, as `n` is (decision 004).

## What would reverse this

1. **If a measured trametinib dose–response in AR42J is inconsistent with a
   single `ic50_cat`** — pre-flight item 5 in Part 7 — the one-channel structure
   may need a second term, and this should be revisited before Stage 1 fixes the
   trametinib axis.
2. **If Stage 3's FIM shows `ic50_cat` and `ic50_act` are not separately
   identifiable**, the two-term drug model is not earning its extra parameter
   from the data, and the honest report is that the trametinib/PD325901
   distinction is a structural assumption rather than a fitted difference.
3. **If ERK is found to act on ignition through a route independent of MEK
   output**, the single-channel rule breaks and `g` genuinely needs a second
   argument — but it would then need to be a different variable, not `v`.
