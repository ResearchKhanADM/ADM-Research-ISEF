# 013 · "Profile likelihood" is not available with no data — what replaces it

*Date:* 2026-08-15 · *Status:* **DECIDED-PENDING-REVIEW** (Tier 2)
*Panel:* methods-checker mandate (against, with measurements) + advocate mandate
(for). Both argued directly rather than through the registered agent files, which
were created this session and are not yet loaded into the agent registry.

**Read this instead of the transcript.** It changes wording in v3 Part 0.7, the
Gate B criterion, and `PHASE2_PARAMETER_BUDGET.md` §1.

---

## Question

v3 Phase 2 requires *"identifiability on the three parameters that matter"*, Gate
B is worded *"key parameters surviving profile likelihood"*, and
`PHASE2_PARAMETER_BUDGET.md` §1 asserted that the profile on each of `a_P`, `γ`,
`κ` **is** the uncertainty on R2, R3, R1 respectively — to be written into figure
captions in those words.

**There is no dataset.** Nothing quantitative has been fitted. The available
observations are qualitative: Collins 2014 (PMID 24315826) reports reversion under
MEK inhibition in ~3 days and relapse within ~7 days of withdrawal.

Can anything called profile likelihood be run and reported, and if not, what
replaces it?

## Positions considered

**For — run something profile-likelihood-shaped now.** Deferring all uncertainty
reporting to a wet lab that may return a negative leaves the poster with no
uncertainty statement at all, which is weaker than a correctly-hedged one. Gate B
is worded around it. Phase 7's discrimination-power argument needs some notion of
parameter uncertainty. And the computation is cheap.

**Against — it is a misnomer of a class this project has already been caught by.**
Profile likelihood requires a likelihood (a density for *observed data*),
re-optimisation of nuisance parameters at each fixed value of the profiled one,
and a Wilks threshold. With no data there is no likelihood; a Gaussian penalty
built around two published timings with tolerances chosen by the modeller is a
penalty, not a likelihood. **This is the same shape as v3 §0.7 row 1**, where
"stability selection" was used for resampling a parameter ensemble — a
data-based method's name applied to a data-free operation on parameter space.

*Where they actually disagree:* not on whether the computation runs — both agree
it does — but on whether its output is an **uncertainty** or an **artefact of the
inputs the modeller chose**.

## Decision

**The name is withdrawn. Three things replace it.**

The against case was settled by measurement rather than argument:

- The profile is **flat**. A crude sample already drives the objective to ~0 in
  nearly every bin of `a_P` and `γ`, and a real optimiser would flatten it
  further.
- The resulting "interval" is **100% of the prior box** for `a_P` and `γ`.
  Reporting it means reporting the numbers typed into the ranges file.
- Where it is not the whole box, **the tolerance is the answer**: halving the
  chosen σ moved `κ`'s interval from 100% to 67% of the box. Nothing about
  Collins changed.

**What is reported instead:**

1. **The fold locus from continuation** — exact for the model, not calibrated,
   and already computed by Gate B. *"The model is bistable only for `a_P` above
   this fold, whose location as a function of the other groups is this"* is a
   sharp, correctly-named statement about a parameter that matters.
2. **Constraint-filtered prior-predictive intervals on the DELIVERABLES**, not on
   the parameters — sample the declared box, keep vectors reproducing Collins
   qualitatively, report the interval of R1/R2/R3 over survivors, **plus the
   acceptance fraction, plus how the interval moves when the box moves.** This is
   v3 §0.7 row 3's own prescription, currently applied only to Phase 6; applying
   it to Phase 2 needs no new terminology.
3. **The flat profile itself, correctly named.** In Raue's framework a flat
   profile is the textbook signature of **structural non-identifiability**.
   *"These parameters are not identifiable from the available published
   observations"* is true, useful, and costs nothing to say.

**Reserve the name for the wet-lab timecourse, and pre-register it now.** When
the bench returns a timecourse with replicates, a genuine likelihood exists and
profile likelihood is the right method under its right name. Pre-registering it
before the design is locked is worth more than running a fake version today — and
it is exactly the "computation that designed an experiment" framing.

**Also decided: the one-to-one mapping claim is withdrawn** (see §"Downstream").

## What would reverse this

1. **A quantitative dataset appears** — the wet-lab timecourse, or a published
   one with replicates and error bars. Then a likelihood exists, `σ` comes from
   the data rather than from the modeller, and profile likelihood is correct and
   should be run under its own name. This is the expected path, not a remote one.
2. **A source supplies real dispersion for the Collins timings.** If those figures
   yield a mean ± SEM with a stated *n*, (A) becomes a genuine — if tiny —
   likelihood and the label becomes defensible for those two observations. The
   Wilks threshold would still not be, with 13 unknowns and 2 constraints.
3. **The prior-predictive intervals turn out to be narrow.** If constraint
   filtering leaves a tight interval on a deliverable, that is a real result and
   this decision's replacement is doing better than the thing it replaced.

## Downstream — what this changes

- **`PHASE2_PARAMETER_BUDGET.md` §1: the one-to-one mapping is withdrawn.**
  Measured: `a_R` — *not profiled* — moves R1's co-formulation gap far more than
  `κ` does across `κ`'s own range, and `n_P` — *a scanned exponent* — can abolish
  the conversion entirely, while `κ` is nearly flat over its first three decades.
  A Sobol analysis replaces the assertion; template pre-registered at
  `prereg/2026-08-15_phase2_sensitivity_*`. Until it runs, the strongest
  defensible phrasing is **"the dominant single contributor to"**, and only where
  the analysis supports it.
- **A parameter interval was never a deliverable interval anyway.** R1 is a
  percentage gap, R3 a window in days. Converting an interval on `κ` into one on
  R1 requires pushing it through the Phase 3 computation — a nonlinear map
  depending on the other eleven groups.
- **v3 Part 0.7 gains a fourth row**, so the error is logged in the project's own
  table of misused terms rather than only here.
- **⚠ Gate B's wording is now unpassable and must change.** It is defined as
  *"key parameters surviving profile likelihood"* — a gate defined in terms of a
  method that cannot be run fails for reasons unrelated to the model. Rewording
  it changes what an endpoint means, so it is a **Tier 3 question**, batched for
  Luqmaan this session. The structural half of Gate B (two stable states, a
  saddle, an identifiable separatrix) **is met** — see `results/gate_b/`.
- **The model has no clock**, and it reaches Phase 5. Time is `τ = t·δ_P`, and
  `δ_P` — PTF1A protein turnover — is absorbed by the nondimensionalisation and
  unmeasured. Any day-valued observation needs it reintroduced as a 13th unknown,
  so *k* timings give only *k−1* dimensionless constraints: **Collins' two give
  one, the ratio ~7/3.** More damaging, **Phase 5's redosing interval has no
  units without it** — mRNA half-life (bench item 4) constrains a product of
  dimensionless groups, not the protein clock. **Bench Handshake item 9 added:
  PTF1A protein half-life by cycloheximide chase.**
