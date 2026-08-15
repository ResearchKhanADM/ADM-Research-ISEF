# 004 · Hill forms, and how `M` enters the acinar enhancer

*Date:* 2026-08-15 · *Status:* accepted

## Question

Two things. (1) What functional form do the enhancer and *Rbpjl* activation
terms take, and is the exponent fixed or sampled? (2) `Hill(C_L, M; n)` is
written as a two-argument function — how does `M` actually enter, given that
constraint 3 places `M` at **metaplasia** loci while the 2.3-kb autoregulatory
enhancer is an **acinar** locus?

## Positions considered

**On the exponent.** *Fix `n` from stoichiometry* — PTF1-L is an obligate trimer
with two conserved PTF1 sites, which arguably justifies `n = 2`. Against:
site count is not Hill coefficient. The Hill exponent is a phenomenological
summary of cooperativity, not a subunit count, and no Hill coefficient for PTF1
enhancer occupancy has ever been measured (Part 8). Fixing it asserts precisely
the quantity we claim to be uncertain about — and it is the exact move the panel
killed when it rejected using a sequence model to obtain `n`.

**On `M`'s entry.** *Direct gating* — write `M` inside the enhancer's binding
expression, as the plan's `Hill(C_L, M; n)` literally suggests. Against: that
places the metaplasia chromatin state physically at the acinar enhancer, which
contradicts constraint 3 and the evidence behind it. Falvo 2023's retained
H3K4me1 is at **metaplasia** genes; the published memory primes the cell to
*leave* the acinar state, not to return.

*Where they actually disagree:* not on whether `M` opposes the acinar program —
both agree it does — but on whether that opposition is **local** (same locus) or
**programmatic** (one program represses the other).

## Decision

**Separable activation × repression, exponents sampled.**

```
enhancer = hill_activate(C_L, k_auto, n_auto) × hill_repress(M, k_M_rep, n_M_rep)
```

`n` is **sampled over 1–4**, never fixed — this is the standing answer to *"you
used a sequence model to get a Hill coefficient"*: we did not, we scan it.

`M` enters **multiplicatively as repression, not as a co-argument of the
activation term.** This keeps `M` at metaplasia loci (constraint 3) while
producing the same qualitative effect: the metaplasia program represses the
acinar program. Separability is the minimal assumption — there is no evidence
for any specific interaction form between the two, and inventing a coupled
binding expression would be asserting mechanism we do not have.

## What would reverse this

1. **If Stage 2 shows the Q-value ranking depends on separability** — test by
   swapping in a coupled form (`M` competing within the activation expression)
   and re-running T4. If the ranking flips, separability was load-bearing and
   must be argued from evidence rather than parsimony.
2. **If sampled `n` turns out to be a stiff direction in Stage 3's FIM
   spectrum** while the conclusions depend on it, the "we scan it" defence
   weakens and the exponent needs an independent constraint.
3. **If evidence appears for metaplasia-chromatin marks at the *Ptf1a*
   autoregulatory enhancer itself**, direct gating becomes correct and this
   decision inverts.
