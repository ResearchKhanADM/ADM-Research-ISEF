# Prior justification — ID3 dissociation constants `K_IE` and `K_IP`

*Written 2026-08-15, before any sweep. To be incorporated verbatim into the
Stage 2 pre-registration.*

**This document justifies the WIDTH and WEIGHTING of a prior. It does not pick a
value, and no sweep may fix one.**

---

## Why this prior is not a routine choice

`K_d/E_tot` is the single parameter that determines whether Stage 2 can
discriminate T1 from T2 at all. From `docs/derivations/binding_polynomial.md`:

```
K_d ≫ totals  →  n_eff → 2          identical to a first-order sink (§5)
K_d ≪ totals  →  n_eff ≈ 1.34·√(totals/K_d)   unbounded as K_d falls (§6)
```

In the loose limit **T1 and T2 are the same model** — not similar, the same.
So a prior that sits entirely in the loose regime guarantees a null result, and
the null would be a property of the prior rather than of the biology.

**No `K_d`, `k_on` or `k_off` has ever been measured for ID3·E47, ID3·PTF1A, or
PTF1A·RBPJ/RBPJL** (master plan Part 8, after targeted search). The prior is
therefore entirely a choice, and it is a choice that determines the answer.

## What is sampled

- **The dimensionless ratio `K_d/E_tot` is the quantity that matters**, not
  `K_d` and `E_tot` separately. Nondimensionalisation should surface it as one
  of the dimensionless groups; if it does not, that is a signal the
  nondimensionalisation is wrong.
- **Log-uniform, spanning at least 4–5 decades**, deliberately bracketing the
  transition at `K_d ~ totals`. The prior must contain both regimes, or the
  discrimination-power analysis has nothing to vary over.
- **`K_IE` and `K_IP` are sampled INDEPENDENTLY** — see the asymmetry below.

## Evidence bearing on the prior — two strands, both soft

Neither strand justifies fixing a value. Both justify **weighting** the prior
toward tighter binding, and both must be stated as the qualitative evidence they
are.

**1. Direct, qualitative — and it argues for an asymmetry.**
Langlands, Yin, Anand & Prochownik, *J Biol Chem* 272(32):19785–93 (1997),
PMID 9242638. A quantitative yeast two-hybrid assay across Id1/Id2/Id3 against
class A E-proteins and class B factors. **All three Ids bound E-proteins with
"high affinity"** — but **Id3 interacted *weakly* with all four class B MRFs**
(MyoD, myogenin, Myf-5, MRF4). Rank-order only; **no dissociation constants are
reported**, so nothing here sets a number.

**PTF1A is a class B bHLH.** It was *not* among the factors tested — the class B
panel was myogenic and hematopoietic — so applying this to PTF1A is an
extrapolation and must be labelled as one. But it argues clearly against tying
`K_IE` and `K_IP` to a shared prior: the evidence that exists points to
**ID3·E-protein tighter than ID3·(class B factor)**. `K_IP`'s prior is therefore
centred looser than `K_IE`'s, with both spanning the full range.

**2. Functional, indirect — and the stronger of the two.**
Dufresne 2010 (PMID 20830706) observes that ID3 **does** effectively mislocalise
PTF1A to the cytoplasm at physiological expression levels in AR4-2J — the
wet-lab line — and that silencing Id3 reverses it. **A loose titrator could not
produce that phenotype.** Sequestration sufficient to relocate the master
regulator implies binding tight enough to compete with complex formation at
native concentrations.

This is the more useful strand precisely because it is a statement about the
system of interest at physiological concentrations, rather than an affinity
rank measured in yeast against different partners.

## What Stage 2 must report

Not a Q-value comparison at one prior. **Discrimination power as a function of
the binding regime:**

> *"T1 and T2 are distinguishable only when `K_d/E_tot < X`. Within the
> pre-registered plausible range, `Y`% of samples fall in that regime. The
> topology competition has power only there, and here is what it concludes."*

This is a candidate headline figure — **discrimination power vs binding regime**
— and it is a genuine methodological result: most modelling papers never state
their discrimination power at all.

## Convergence caveat that must accompany any Q-value

Per the `CLAUDE.md` standing rule: failure rates are reported as a function of
`K_d/E_tot`, because dropped failures would deplete the sample set precisely in
the discriminating regime and bias the comparison toward a false null. Current
measured rate after the Step 2 solver rewrite: **0 failures in 4000 solves
across 12 decades of `K_d/E_tot`**. Note also that **T2 is closed-form and
cannot fail**, so any nonzero T1 failure rate is an asymmetry that must be
disclosed rather than averaged away.
