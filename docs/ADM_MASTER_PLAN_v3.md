# Making ADM reversal *stick*
## Master plan v3 — supersedes all earlier versions

> ## ⚠ READ FIRST
>
> **This replaces `ADM_INSILICO_MASTER_PLAN.md`, `INSILICO_PLAN_v2.md`, and `INSILICO_PLAN_8months.md`.** They have been deleted. If any file on this machine still describes a seven-stage pipeline, a five-way topology competition, a 2^k subset enumeration, or a Lie-bracket ordering analysis, it is obsolete.
>
> A three-reviewer panel found four problems that the previous plan could not survive. All four are recorded in Part 0 with the evidence. Read Part 0 before Part 1 — it explains why the project changed shape.

> ## AMENDMENTS — v3.1, 2026-08-15
>
> Four changes, made after auditing v3 against the existing repo. Each is recorded in `docs/decisions/012-durability-framing-architecture-change.md`; the affected sections below are edited in place rather than appended to, so this file stays the single plan of record.
>
> 1. **§2 Phase 2 — pERK is still an input, but no longer a step.** On withdrawal it follows a prescribed, swept **rebound profile**. v3's spine is durability *after withdrawal*, and pERK rebound is a withdrawal-specific mechanism; deleting it silently would have deleted the thing most relevant to the new endpoint. This converts an unmeasured state into a **measured input** — a pERK western timecourse after washout — which is the move this plan makes everywhere else. **Bench Handshake #1 gains that timecourse.**
> 2. **§1.4 — the ID3 limitation was too broad and is rewritten.** Dufresne 2010 (PMID 20830706) *does* document the two-target node, in AR4-2J. The real gap is one edge: **ERK→ID3**. **Bench Handshake #1 gains an ID3 western ± trametinib**, which is cheap and decides whether Phase 5 means anything.
> 3. **§2 Phase 2 — `E_free` uses the exact binding solution, not the linear approximation.** `E_total − k·ID3` is the *tight-binding limit* of the exact form and goes negative outside its domain. The derivation already exists (`docs/derivations/binding_polynomial.md`); it becomes this plan's justification rather than dead work.
> 4. **Viability is a bench-measured floor, not a modelled hazard.** The U-shaped hazard is retired (decision 008 amendment). One thing is kept as an **output flag, not a term**: if the model predicts `P` below the CHOP-apoptosis threshold, say so on the output.

---

# PART 0 — WHY THIS PLAN REPLACED THE LAST ONE

## 0.1 "Reversal" was already solved in 2014 — and the gap is what happens next

Collins et al., *Gastroenterology* 2014 (PMID 24315826): MEK inhibition reverts established PanIN in **~3 days**, ~50% amylase-positive by day 3, cleaved-caspase-3 rare throughout. A generic MEK inhibitor. In vivo. No transcription-factor payload.

Then, from the same paper's full text: after 5 days of inhibition they **stopped the drug** while leaving KRAS on for 7 more days —

> *"Analysis of the tissues showed that ADM and PanIN formation had resumed."*

Collins names the gap in her own discussion: *"MAPK inhibition might result in quiescence of the redifferentiated cells, but release of the inhibition might lead to rapid relapse."*

**So a project framed as "reverse ADM" is eleven years late, and a judge will say so.** The unoccupied claim is the one Collins pointed at and did not take:

> **Trametinib reverts. It does not stick. A payload that re-closes the PTF1A↔RBPJL loop should convert a drug-dependent reversion into a self-sustaining one. This project finds the minimum payload, ratio, formulation and schedule that does it.**

Everything below is built on that sentence.

## 0.2 The "discovery" framing does not survive contact with the literature

Three reviewers found this independently, three different ways:

1. **The answer is already in an abstract.** Jiang et al., *Gastro Hep Advances* 2023 (PMID 37425649), primary human acinar cells, 3D Matrigel, 14 donors, day 0 vs day 6, verbatim: *"transcription factors with reduced activity (PTF1A, RBPJL, and BHLHA15)."* BHLHA15 **is** MIST1. A regulon-activity screen run correctly returns exactly the assumed payload with the third slot filled.
2. **The screen is circular by construction.** Regulon activity is scored from target-gene expression. The PTF1A regulon's targets are the digestive-enzyme genes — the same genes that *define* the acinar cluster. The regulon best separating acinar from metaplastic is PTF1A by construction. This is a **different and harder** circularity than the ODE-representability bias previously logged.
3. **It is not a screen.** No cells are perturbed. It is computational prioritization of public data. Applying "stability selection" to bootstraps of one dataset measures robustness to resampling, not to biology.

**Consequence:** the identity claim is dead. The composition claim is alive and unpublished.

> Not *"I discovered which factors."* → **"I independently recovered the published regulon loss, then asked what the published work did not: how many components, in what ratio, under a fixed delivery budget."**

## 0.3 The 10× kinetic asymmetry — the discriminator for the whole topology competition — is not measured

The previous plan's Stage 2 (4 weeks, called "the scientific core") was discriminated by "MEKi reverts in 3 days, forced PTF1A takes 3 weeks."

**Krah 2019's DOX timepoints are 24 hours, 3 weeks, and 6 weeks. There is nothing in between.** The paper does not measure that PTF1A takes three weeks. It measures that it takes *more than one day and no more than twenty-one*.

And Collins 2014 gets a ~2× spread within one paper: iKras\* reverts by day 3, KC shows no acinar cells until after day 5, same drug, pERK suppressed equally in both.

So the asymmetry is bounded on the evidence at **(1×, 21×)**, from two papers differing in species-model, lesion stage, readout modality, and — decisively — temporal sampling density. The number 10 was chosen from that interval.

**A model-selection step resting on an unsourced number cannot be the scientific core.** The five-way competition is cut. What replaces it is a smaller model whose parameters are actually identifiable (§2.2).

## 0.4 The ordering-novelty claim was false

The previous plan claimed no published precedent for formalizing intervention ordering in cell-fate control. **Letsou & Cai, *PLoS Comput Biol* 2016 (PMID 27560383), "Noncommutative Biology: Sequential Regulation of Complex Networks"** formalizes exactly this and derives scaling laws. Order-dependent induction navigating cells to different terminal states is also established (PMID 28397688).

Being caught on a false priority claim costs more than the correct result earns. The ordering question survives as **one number in plain English with a correct citation** — not as a formalism.

## 0.5 What the panel found that is *better* than what was cut

**Delivery is missing from the plan, and it is where the strongest available result is hiding.**

Per-cell LNP dose is a **distribution**, not a number. Rees et al., *Nat Commun* 2019 (PMID 31138801): the statistical distribution of nanoparticle dose per endosome is independent of administered dose; what scales with dose is the *number* of NP-containing endosomes — so per-cell dose behaves as a compound-Poisson process and heterogeneity is intrinsic. Dobrowolski et al., *Nat Nanotechnol* 2022 (PMID 35768613) measured LNP-delivered barcode, mRNA, protein output and transcriptome at single-cell resolution and found cell subtypes with systematically high and low uptake.

Now combine that with this project's own mechanism, which is **stoichiometric and thresholded**: PTF1A needs an E-protein partner, RBPJL cannot bootstrap, so the pair is required together in the *same cell*.

**Co-formulated in one particle**, per-cell doses of two mRNAs are near-perfectly correlated → fraction of cells with both above threshold ≈ the marginal, ~30%.
**In separate particles**, doses are approximately independent → ≈ 0.30 × 0.30 ≈ **9%**.

**For an obligate stoichiometric partnership, co-formulation is not an optimization — it is a requirement, and the model can prove it and quantify the gap.** Neither a screen nor a 13-state ODE would ever generate that result. It costs ~2–3 weeks. It is the reason to keep a model at all.

## 0.6 2^k subset enumeration is the wrong method

A fixed total mRNA mass budget makes this a **mixture experiment**: components whose proportions sum to a constant live on a simplex. **Subset enumeration visits only the vertices of that simplex** — {PTF1A alone}, {PTF1A + RBPJL at 50:50} — and structurally cannot express 70:30, which is precisely where a stoichiometric-titration mechanism lives.

The standard tool is a **mixture-amount design** with a Scheffé polynomial: proportions on the simplex, crossed with total amount (Cornell, *Experiments with Mixtures*, 3rd ed.).

Also note: submodular optimization is *inappropriate* here. The (1−1/e) greedy guarantee requires diminishing returns. This project's central premise — a bootstrap threshold where PTF1A alone does nothing and PTF1A+RBPJL crosses it — is a claim of **super**modularity. Using submodular optimization would assume away the phenomenon under study.

## 0.7 Corrected terminology — misused statistical terms a reviewer will catch

| Was called | Actually is | Fix |
|---|---|---|
| "stability selection" | Resampling a parameter ensemble | Stability selection resamples **data**. Do it across ≥2 **independent datasets** and report per-dataset frequency, or call it a sensitivity analysis. |
| "selection / confirmation split" | Train/test across halves of the same generative model | Controls parameter overfitting, says nothing about **model misspecification**, which is the dominant error source here. Rename to dataset-level cross-validation and say what it does not cover. |
| "held-out prediction" | Constrained on 2 observations with ~30 parameters | With 30 parameters, 2 constraints remove almost nothing. Report the prediction as an **interval over all parameter vectors passing the filter**. If the interval spans everything, there is no prediction. |

---

# PART 1 — THE PROJECT

## 1.1 The system

| | |
|---|---|
| **Cells** | AR42J, ATCC parental (CRL-1492) — **not** B13 |
| **ADM induction** | KRAS G12D, inducible (Tet-On) |
| **Brake** | Trametinib |
| **Payload** | mRNA in LNPs — composition, count and ratio are **outputs** |
| **Format** | Growth-factor-reduced Matrigel, 3D |
| **Wet lab** | Same cycle. **One shot, no iteration.** |
| **Primary endpoint** | **Durable reversal after payload clearance with trametinib withdrawn** |
| **Constraint** | Viability ≥ a floor set before results are seen |

## 1.2 The mechanism, in one paragraph

PTF1A works as a trimer: PTF1A + an E-protein (E47/HEB) + RBPJL. PTF1-L drives **both** its own 2.3-kb autoregulatory enhancer **and** its partner *Rbpjl*. *Rbpjl* has **no PTF1A-independent promoter**, so once the loop opens it cannot re-close on its own — a bootstrap failure. KRAS→ERK cuts the loop four ways: suppresses the PTF1-independent ignition promoter; starves *Rbpjl*; induces ID3, which stoichiometrically titrates the E-protein that PTF1A requires and that carries nuclear import; and drives TRIP12-mediated PTF1A degradation (K312 human = K311 rat = K309 mouse, same residue, species numbering offset). Chromatin at metaplasia loci retains a slow memory that primes relapse.

**Trametinib lifts the ERK-dependent cuts. It cannot rebuild RBPJL, because nothing but the loop makes RBPJL. That is why reversion is drug-dependent, and that is what the payload is for.**

## 1.3 The three results this project is built around

**R1 — Formulation.** Co-formulation versus separate particles, quantified across the plausible uptake-CV range. For an obligate stoichiometric pair, this is a requirement, not a preference.

**R2 — Composition.** How many components, at what ratio, under a fixed total mRNA mass. Stated as a number with a marginal-value curve: *"two components; the third buys under X% of durable reversal for 33% of the mass budget."*

**R3 — Durability.** The drug-free persistence window — the region of parameter and dose space where the reverted state survives trametinib withdrawal. **This is what the bifurcation diagram is for.** In the old plan the hysteresis wedge was a decorative figure with no job; here it is the answer.

Everything else in the plan is in service of these three.

## 1.4 Known limitations — state these before anyone asks

- **The prioritization is a positive control, not a discovery.** Say it on the poster.
- **A regulon screen is blind to post-translational mechanisms.** TCF3/E47 is not transcriptionally lost — it is titrated by ID3 — so no threshold will make it appear. It enters by **declared mechanism**, registered before the screen runs. That is not an excuse; it is a correct statement about two instruments with disjoint blind spots, and it is the strongest argument that both arms are needed.
- **The ID3 titration node is documented in this cell line. The load-bearing gap is one edge: ERK→ID3.** *(Corrected 2026-08-15 — the previous text said targeted search returned no evidence in the pancreatic context. That was too broad, and it would have had this project disclaiming published work on the poster. Verified against PubMed.)*

  **Documented** — Dufresne et al. 2010, *Int J Cancer* 129(2):295–306, PMID 20830706, doi 10.1002/ijc.25668, *"Id3 modulates cellular localization of bHLH Ptf1-p48 protein"*, **in AR4-2J**: Id3 binds **both** E47 and Ptf1-p48; a proliferative signal raises Id3 and raises Id3/E47 **and** Id3/Ptf1-p48 while E47/Ptf1-p48 **falls**; **silencing Id3 reverses the cytoplasmic mislocalization**; and the pattern — Id3 overexpressed, Ptf1-p48 absent or mislocalized — holds in **human and murine preneoplastic lesions**. E47 carries nuclear import of the PTF1 complex, which is the stated mechanism of the mislocalization.

  **Not documented** — any **Kd** for either interaction. Unchanged, and it is why binding affinities are sampled across decades rather than fixed.

  **Assumed, and this is the load-bearing gap** — **that ERK drives ID3.** Dufresne's driver is **gastrin**. The ERK→ID3 edge is an inference, not a measurement.

  **Why it matters where it matters:** Phase 5's ordering prediction lives on exactly that edge. If trametinib does not lower ID3 in these cells, trametinib-first loses its mechanistic basis. That is a **one-western experiment** — ID3 by western, ± trametinib, in our cells — and it is on Bench Handshake #1 as high-value. Net: the prediction is **better supported than v3 originally claimed**, and the residual uncertainty is one edge rather than the whole node. Say it that way.

- **Viability is measured at the bench, not modelled.** The old U-shaped death hazard is retired (decision 008 amendment): its parameters were never measurable, ER-stress dynamics were the least-constrained part of the old model, and the bench measures viability directly. Modelling something you can measure, with parameters you cannot, is backwards. **Kept as an output flag, not a term:** if the model predicts `P` below the CHOP-apoptosis threshold, surface it as a warning. The low-PTF1A arm is real biology (Sakikubo 2018, PMID 30361559); it does not need to be a differential equation.
- **PTF1A is pleiotropic.** It reprograms fibroblasts into tripotential neural stem cells, and is dosage-sensitive in both pancreas and cerebellum (PMID 30470852). "Preserving viability" does not cover "did not create a neural-program-expressing pancreatic cell." Lineage fidelity is a separate endpoint.
- **Converted cells often fail to silence the starting program** (CellNet, PMID 25126793). The candidate list is all *turn acinar back on* and needs an **ADM-repressor axis**.
- **AR42J requires dexamethasone to express amylase at a differentiated level at all.** Baseline PTF1A/RBPJL and the achievable dynamic range must be measured, or the experiment's falsification power is unknown.
- **AR42J Kras genotype is unverified.** Azaserine-induced rat pancreatic tumours are classically Kras-mutant. Sequence Kras codons 12/13 and 61 in week 1.
- **Reference data is mouse/human; the bench is rat, and a transformed line.** Cross-species transfer must be argued, not assumed.

---

# PART 2 — THE PIPELINE

Eight phases, ~18 weeks of work. Dependencies and gates stated. Phase 1 runs parallel to Phase 2.

## PHASE 0 · Decision spec — 1 week · do this first

**Do.** Write one page and lock it: cell system and passage; **how many experimental arms/wells actually exist**; assays and readouts; the viability floor X; the durability timepoint; and the precise definition of the reversal score.

**Why first.** The previous plan optimized before knowing how many wells exist — which is how you end up with a 256-point Pareto front and twelve wells. Every downstream design decision is a function of the arm budget.

**The arm budget is a blocking EXTERNAL input, not a modelling choice** *(amended v3.1)*. It depends on the wet lab and is not known yet. Do not guess it and do not let it block Phases 2–3, which do not depend on it. Instead:

- put an explicit **placeholder range in the decision spec — 12 / 24 / 48 wells** — and
- build Phase 4's mixture-amount design so it can be **instantiated at any tier**, and
- **state what changes at each tier**, so the number can be argued for rather than accepted.

A design that degrades gracefully across tiers is worth more than one that is optimal at a well count nobody has confirmed.

**Also issue Bench Handshake #1 now** — these have long lead times, and Phase 3 cannot run without them. **All items are blocking.**

| # | Item | Blocks | Note |
|---|---|---|---|
| 1 | *Kras* genotype — sequence codons 12/13, 61 | the forcing term | highest-stakes unknown in the project |
| 2 | Trametinib IC50 in *these* cells | Phase 2 drug axis, Phase 5 | a published IC50 from another line sets the wrong axis origin |
| 3 | **LNP transfection efficiency AND its cell-to-cell CV** | **Phase 3 ★** | the CV *is* the x-axis of the headline figure; a mean without a spread cannot produce the result |
| 4 | mRNA half-life in AR42J / Matrigel | Phase 5 redosing axis | 3D is not 2D; a plastic-dish number is a different number |
| 5 | Baseline PTF1A/RBPJL by qPCR ± dex | Phase 2 ICs, Gate A | without dynamic range, a null is unreadable: "payload failed" or "no room to move" |
| 6 | ADM stability to ≥14–21 days | Phase 5 durability endpoint | no contrast if the state does not hold unassisted |
| 7 | **pERK western timecourse after trametinib washout** | **Phase 2 input, Phase 5 ★** | *added v3.1.* Supplies the rebound profile directly instead of sampling its shape — converts an unmeasured state into a measured input |
| 8 | **ID3 by western, ± trametinib** | **Phase 5 — decides whether it means anything** | *added v3.1.* **High value, one western.** Tests the ERK→ID3 edge, the single assumption Phase 5's ordering prediction rests on. If trametinib does not lower ID3 here, trametinib-first loses its mechanistic basis |

Items **3** and **8** are the ones to chase hardest: 3 is the only input to the strongest result and no literature value substitutes for it; 8 is cheap and decides the fate of a whole phase.

## PHASE 1 · Candidate generation — 2 weeks · parallel with Phase 2

**Do.** Regulon activity (pySCENIC-style) across **≥2 independent public datasets**, reported per dataset with cross-dataset frequency — that is what makes a stability claim mean anything. Add an **ADM-repressor axis** to the candidate classes, not just activators. Cross-reference every candidate against the Joung TF Atlas (PMID 36608654) and the human TFome hit list (PMID 33257861) for known potency and known toxicity — free prior information, one day of work.

**Name it honestly: prioritization, not discovery.** Report it as a positive control that had to work.

**CellOracle: one run, as a declared negative control.** Pre-register the prediction that it will fail to recover PTF1A/RBPJL because the GRN is fit in a cluster where those genes have no variance. When it fails, that is a demonstration of method understanding. Framed as "orthogonal validation," it is a self-inflicted wound.

**GATE A (end of week 2):** does the pipeline recover PTF1A / RBPJL / BHLHA15? Cheap, fast, decisive. If no, the pipeline is broken and you know in month one.

**→ Layer 1 payload lock. Order PTF1A + RBPJL mRNA now.** These two are not in doubt — published *and* mechanically forced by the bootstrap argument. The bench's longest-lead item (IVT synthesis, formulation, QC) starts in week 2 instead of week 11.

## PHASE 2 · Minimal mechanistic core — 4 weeks

**Do.** Build **3–4 states, ~9–12 parameters**, not thirteen states and thirty parameters:

- `P` — PTF1A activity, autoregulatory, requiring an E-protein partner
- `R` — RBPJL, produced **only** as a function of P. No P-independent term. That zero *is* the bootstrap claim.
- `C` — slow chromatin/memory state at metaplasia loci
- **Algebraic, not differential: `E_free` from the exact binding solution**, with `ID3 = f(pERK)`. *(Amended — see below. `E_free = E_total − k·ID3` is the tight-binding limit of that solution, not an independent modelling choice, and it goes negative once `k·ID3 > E_total`.)*
- **pERK is an input, not a state — and on withdrawal it is not a step.** *(Amended — see below.)*

**`E_free`: use the exact solution, cite the derivation.** Linear subtraction is what two-target titration reduces to when `Kd ≪` the protein totals: in that limit each ID3 molecule takes one E molecule 1:1, so `E_free → E_total − k·ID3`. That makes it the **ultrasensitive** limit, not a soft approximation — but it is only valid inside its domain, and a `max(0, ·)` floor is a hack that hides the domain violation rather than fixing it. Use the exact form from `docs/derivations/binding_polynomial.md` (already derived, previously written for a cut stage) and state in the code comment that it reduces to the linear form in the tight regime, with

```
n_eff ≈ 0.5·√(E_tot/Kd)         the diagnostic for when the approximation holds
```

as the check. No invalid region, no floor, and the derivation earns its keep.

> **Correction, 2026-08-15, found on first integration of the core.** This line originally read `1.34·√(E_tot/Kd)`, carried over from `docs/derivations/binding_polynomial.md` §6. **That constant does not apply to the shipped mechanism.** It was derived for the deleted 11-state model's *ternary complex* under *two-target* titration (ID3 taxing both PTF1A and E-protein, so the log-log slope compounded). The Phase 2 core titrates one target, and the measured prefactor on `E_free` is **0.5**. Carrying 1.34 across that change of mechanism would have overstated threshold sharpness by ~2.7×.

**Why that matters downstream, and why it was worth catching:** Phase 3 convolves the per-cell dose distribution against this threshold, and a sharper threshold produces a **larger co-formulation gap — the headline number of the entire project**. An inherited constant would have inflated R1 by ~2.7× while the model still ran and the figure still rendered.

**Open item for Phase 3:** the quantity actually convolved is the *complex* threshold `P·E_free·R`, not `E_free` alone, and its sharpness compounds `n_P`/`n_R` on top of the binding term. **Phase 3 must define sharpness on the quantity it convolves and measure it there**, rather than reusing either constant.

**pERK on withdrawal: a swept rebound profile, not a step.** Trametinib sets pERK while present. When it is withdrawn, pERK does not jump instantly to baseline; it recovers on its own timescale, with a possible overshoot from relief of ERK-mediated negative feedback on RAF. Because the primary endpoint of this whole plan is *what happens after withdrawal*, that recovery shape is not a detail — it sets how hard the system is pushed back toward the metaplastic basin at exactly the moment durability is being tested.

Model it as a **prescribed, parameterized input curve** — rise time, overshoot amplitude, settling time — **swept across its plausible shape**, not as a state. This buys the phenomenon without a phospho-MEK state, without a coupled algebraic loop, and without spending identifiability on parameters nobody has measured.

**And it is measurable.** A pERK western timecourse after trametinib washout in AR42J hands the model the actual curve, converting an unmeasured state into a measured input. That is on Bench Handshake #1.

**Explicitly retired:** the trametinib-vs-PD325901 withdrawal-asymmetry prediction. It was a free-prize prediction about an experiment nobody is running. See decision 002's amendment — retired, not merely set aside, because an unretired claim is one someone builds on later.

**Nondimensionalize.** Thirty parameters with no nondimensionalization is conspicuous to anyone with the relevant PhD. Collapsing to ~8 dimensionless groups takes days and directly blunts the "none of your parameters are measured" attack, because *ratios* and *timescale ratios* are often constrained even when absolute rates are not.

**Use profile likelihood on the three key parameters, not an FIM eigenspectrum.** Sloppiness was a symptom of over-parameterization that no longer exists. Profile likelihood runs in an afternoon and is more convincing.

**Deliverables:** the P–R bootstrap threshold as a function of free E47; the **two-parameter bifurcation in (KRAS × trametinib) reframed as the drug-free persistence window**; identifiability on the three parameters that matter; and the signed ordering prediction (below).

**GATE B (end of week 6):** two stable states plus a saddle, with an identifiable separatrix, in a model whose key parameters survive profile likelihood.

## PHASE 3 · Delivery layer — 3 weeks · ★ THE HEADLINE RESULT

**Do.**
1. Model per-cell LNP dose as a distribution calibrated to Rees et al. (PMID 31138801) — compound-Poisson / lognormal, parameterized by the uptake CV the bench measures.
2. Convolve with the Phase-2 threshold → **converted fraction as a function of mean dose**, which is not the same shape as the single-cell dose–response and is a first-order effect on the headline endpoint.
3. **The decisive computation:** double-above-threshold fraction under **correlated** (co-formulated) versus **independent** (separate particles) per-cell dose, swept across the plausible uptake-CV range.
4. Convert the mass budget into a **protein-stoichiometry** budget. Translation output per unit mRNA mass varies with ORF length and UTR; a budget in mass with a mechanism in molecules is unphysical without a stated conversion.

**Deliverable:** a formulation recommendation, a mass ratio, and a total mass — each with its sensitivity to uptake CV. **One figure, one recommendation.**

Also specify **route and tropism**. Standard LNPs traffic to liver systemically. If the intent is intraductal, intrapancreatic or ex vivo, say so.

## PHASE 4 · Composition by mixture design — 2 weeks

**Do.** First, cut the candidate list to **≤4–5 components, and justify the cut with arithmetic rather than preference**: Phase 3 lets you state quantitatively that you cannot deliver 8 mRNAs at above-threshold per-cell dose within a fixed mass budget. That argument alone dissolves the 2^k problem.

Then run a **mixture-amount design** — proportions on the simplex crossed with total amount, analyzed with a Scheffé polynomial. If the assay only permits present/absent, use a **definitive screening design** instead (~2m+1 runs for m factors, main effects unconfounded with two-factor interactions).

**Test synergy against a stated null** — Bliss independence or Loewe additivity. Without an additivity null, "three components beat two" is not evidence of interaction.

**Deliverable:** the marginal-value curve — best achievable durable reversal as a function of payload size, under fixed total mass — with the knee marked and the sentence written: *"k components capture X% of achievable durable reversal."* **Report the shape, not the absolute numbers.** Draw the ensemble, not a single curve; a single curve from an unidentifiable model is a lie by graphic design.

**→ Layer 2 payload lock (end of week 11).** Composition, count, mass split, trametinib arm.

## PHASE 5 · Durability and schedule — 2 weeks

**Do.** The primary endpoint. Time-to-relapse after payload clearance **with trametinib withdrawn**, set by the slow chromatin state `C`. Then the dose × redosing-interval map, three regions (undershoot / success / toxic), **using the bench's LNP efficiency and mRNA half-life so the axes have units.**

**Deliverable:** the persistence window, the redosing interval, and — the sentence that carries the project — how long the payload must hold before withdrawal no longer relapses.

**Ordering, reduced to its useful core.** Not n! and not a Lie-bracket formalism. A **two-block question**: does the ERK-lowering (permissive) block precede the TF (instructive) block? Three arms — before / simultaneous / after — at most two gap lengths. Six wells.

The model makes a **signed, falsifiable prediction**: PTF1A mRNA delivered while ERK is high arrives in a cell where its obligate partner is sequestered by ID3 — wasted payload. **Trametinib-first should strictly dominate.** Report as one number in plain English — *"doing it backwards costs N-fold"* — citing Letsou & Cai (PMID 27560383) for the general formalism. No priority claim.

**Be precise about what a null here would mean** *(amended v3.1)*. The titration node itself is documented in AR4-2J (§1.4, Dufresne 2010) — a null does **not** refute it. The assumption this arm actually tests is the **ERK→ID3 edge**: that trametinib lowers ID3 in these cells at all. **Check that first with Bench Handshake item 8** — one ID3 western, ± trametinib. If trametinib does not lower ID3 here, trametinib-first has no mechanistic basis and the ordering arms should not be run; six wells are better spent elsewhere. That is a cheap pre-flight, not a post-hoc caveat.

## PHASE 6 · Held-out prediction — 2 weeks

**Do.** Constrain the model on WT caerulein recovery kinetics and the Krah 2015 PTF1A dose ladder. Then predict, tuning nothing:

1. **The KRAS-history effect** — iKras\* reverts by day 3; KC shows no acinar cells until after day 5, with pERK suppressed equally in both. Same drug, same dose, different history.
2. **The viability dissociation** — MEK inhibition gives reversal with cCasp3 rare (Collins 2014); KRAS extinction at the same timepoint gives reversal with *"extensive cell death"* (Collins 2012, PMID 22232209). Two interventions, same pathway, opposite corners.

**Report both as intervals over all parameter vectors passing the constraint filter, never as point estimates.** With few constraints and many parameters, a point prediction is meaningless and an interval that spans everything is an honest null.

**And report one thing the model gets wrong.** This is the highest score-per-week component in the project because it is the only part that *could have been false*.

**GATE C:** does the history effect fall out of a single parameter?

## PHASE 7 · Pre-registration and discrimination power — 2 weeks

**Do.** Simulate the planned wet-lab design under the Phase 2–3 model and verify the arms can actually **separate the competing hypotheses**. This is expected-value-of-information logic done cheaply.

**Design to discriminate, not to confirm.** With one shot and no iteration, the highest-value allocation is not triplicates of the predicted best — it is arms that separate surviving hypotheses:

- PTF1A alone / RBPJL alone / PTF1A+RBPJL at two ratios → separates the bootstrap-threshold story from the additive story
- Trametinib-first vs simultaneous → tests the titration node
- **Trametinib alone** → the internal positive control for relapse. If it reverts and then relapses, you have validated the *system* even if the payload does nothing. **The experiment then cannot fail completely.**
- eGFP-mRNA LNP at matched mass → LNPs are immunostimulatory and alter the transcriptome (PMID 35768613); this control is mandatory
- Vehicle

**Spend one arm on a predicted negative** — a composition the model says is dominated. A one-shot experiment designed only to succeed wastes the shot if it fails.

**Pre-specified decision rule, written before any data exist:** maximize durable reversal subject to viability ≥ X. Report the **therapeutic index** — the ratio of viability-limiting dose to efficacious dose. **No Pareto front in the deliverable.** A front is preference-free by construction; a one-shot experiment consumes a decision, not a set.

## PHASE 8 · Figures, writeup, buffer — 3 weeks

---

# PART 3 — TIMELINE

| Phase | Weeks | Depends on | Gate | Locks |
|---|---|---|---|---|
| 0 · Decision spec | 1 | — | — | Arm budget, endpoints, viability floor |
| 1 · Candidate generation *(parallel)* | 2 | — | **A: recovers PTF1A/RBPJL/BHLHA15?** | **Layer 1 payload → order mRNA** |
| 2 · Minimal core | 4 | 0 | **B: bistable + identifiable?** | Persistence window |
| 3 · Delivery layer ★ | 3 | 2, bench CV | — | **Formulation + ratio** |
| 4 · Mixture design | 2 | 1, 2, 3 | — | **Layer 2 payload: count + split** |
| 5 · Durability + schedule | 2 | 3, 4, bench PK | — | **Redosing interval + order** |
| 6 · Held-out prediction | 2 | 2 | **C: history effect from one parameter?** | Credibility |
| 7 · Pre-registration + power | 2 | all | — | **Locked protocol** |
| 8 · Figures, writeup, buffer | 3 | — | — | — |
| **Total** | **~18–21** | | | |

**Never cut:** Phase 0, Phase 2, Phase 3, Phase 7.
**Cut order under compression:** Phase 6 → Phase 4's simplex interior (keep vertices only) → Phase 5's ordering arms.

---

# PART 4 — WHAT WAS CUT, AND WHY

Do not re-propose these. Each has a specific reason.

| Cut | Reason |
|---|---|
| **Five-way topology competition** | Its discriminator (the 10× kinetic asymmetry) is not measured — Krah 2019 has no timepoint between 24 h and 3 weeks. With ~30 free parameters, Q-value under a shared sampling box ranks prior volume, not evidence. |
| **13-state ODE** | Cannot rank a candidate it does not contain. A "modular framework" where any TF is inserted with a declared mechanism adds ~3 unidentifiable parameters per candidate — the least defensible option available and the most sophisticated-looking. |
| **FIM sloppiness analysis** | A symptom of over-parameterization now removed. Replaced by profile likelihood on three parameters. |
| **Lie-bracket ordering formalism** | Precedent exists (Letsou & Cai 2016, PMID 27560383) — the novelty claim was false. On a 3-state system the sign follows from the titration term by inspection. Survives as one number in plain English. |
| **2^k subset enumeration** | The mass budget makes this a mixture problem; enumeration visits only simplex vertices and cannot express the ratios where titration lives. |
| **Pareto front as deliverable** | Preference-free by construction, and **a one-shot experiment consumes a decision, not a set.** The two axes are also coupled through MEK — trametinib buys redifferentiation by blocking proliferation (PMID 28090569) — so much of the "front" is one dose axis re-plotted in two dimensions. **This is an explicit reversal, not an omission:** the previous plan and decision 001 both named the reversal-vs-viability Pareto question "a headline result and the stated objective". **Replaced by constrained optimization** — maximize durable reversal subject to viability ≥ X — **plus a therapeutic index**, the ratio of viability-limiting dose to efficacious dose. |
| **Trametinib-vs-PD325901 withdrawal asymmetry** | *Retired v3.1.* A free-prize prediction about an experiment nobody is running. It required a phospho-MEK state (`W`) to have anywhere to live, and pERK is an input here. Retired rather than set aside — see decision 002's amendment. |
| **U-shaped viability hazard as a differential equation** | *Retired v3.1.* Its parameters were never measurable, ER-stress dynamics were the least-constrained part of the old model, and the bench measures viability directly. Modelling what you can measure, with parameters you cannot, is backwards. Survives as an **output flag**: warn when predicted `P` falls below the CHOP threshold. See decision 008's amendment. |
| **Gillespie bimodality vs public scRNA-seq** | Replaced by delivery-induced per-cell heterogeneity, which is the bimodality that *changes the recommendation* and is falsifiable at the bench. |
| **CellOracle as validation** | Predicted to fail on the known true positive for structural reasons. Keep as a declared negative control; presenting it as orthogonal validation hands a hostile judge a loaded question. |
| **Enformer/Borzoi, AlphaFold→k_on, Perturb-seq, ipTM linker scans, one-sided U_crit, structural/Kalman controllability, full-dimension HJ reachability, MPC, all-atom MD** | Killed by earlier panels. Reasons in git history. |

---

# PART 5 — FIGURES

**General-purpose figures. No IEEE or conference-template geometry — sizing is handled downstream.**

**The rule that makes this work: figure modules never compute science.** A figure module loads `results/`, does display arithmetic, and draws. It never runs an ODE solve or a sample. If a figure needs a number not in `results/`, the fix is a stage that writes it. Consequence: figure rebuild is always under a minute, figures cannot drift from the analysis, and a slow stage never blocks a figure.

```
figures/
  _style.py          one house style, one scale knob
  _provenance.py     stamp_run(), save_figure(), is_stale()
  figNN_<topic>.py   one module per figure
  out/               TRACKED: final PDFs + _source.csv + .prov.json
build/figures/       GITIGNORED working renders and proofs
make_figures.py      the single command
```

**Style.** Sans-serif; a single `SCALE` constant so everything resizes together; Okabe–Ito palette capped at four nominal colours (`#0072B2` acinar, `#D55E00` metaplastic, `#009E73` intervention/success, `#6E6E6E` toxic — **grey, not red**, because green/red is the textbook deuteranopia collision); `pdf.fonttype 42`, `svg.fonttype "none"`, `savefig.bbox "standard"` **not** `"tight"` (tight silently resizes and breaks width checks); no gridlines, no top/right spines, axis labels never bold. Use an `rc_context` manager, never a bare rcParams mutation — one process imports every module and a leaked mutation is the classic "it looked right alone" bug.

**Global grammar, fixed once:** a **solid line is a model**; an **open marker with a dark edge is data**. Never connect experimental points with a line.

**Provenance.** Every stage calls `stamp_run()` writing git commit + dirty flag, seed, versions, and SHA-256 of every output. Every figure declares inputs; `save_figure()` hashes them and writes `<slug>.prov.json` and `<slug>_source.csv`. **`paper` and `poster` profiles refuse to render from a dirty tree.** `make_figures.py --check` fails if anything is stale.

**The core figures:**

| Slug | Phase | What it shows |
|---|---|---|
| `fig01_loop_schematic` | — | The reduced 3–4 node core. Not all thirteen states. |
| `fig02_persistence_window` ★ | 2 | (KRAS × trametinib), bistable wedge shaded + hatched, labelled **"reversion persists at zero drug"**. The headline. |
| `fig03_formulation` ★ | 3 | Double-above-threshold fraction vs uptake CV, co-formulated vs separate particles. Two curves, one gap, one recommendation. |
| `fig04_marginal_value` ★ | 4 | Best durable reversal vs payload size under fixed mass. Knee marked. Ensemble drawn, not a single curve. |
| `fig05_durability` | 5 | Time-to-relapse after clearance; dose × interval three-region map. |
| `fig06_heldout` | 6 | Prediction interval vs published numbers, with the prediction box **dated and placed left of the outcome box**. |
| `fig07_prioritization` | 1 | Ranked regulon activity across ≥2 datasets, labelled as a positive control. |
| `figS0x_*` | — | Profile likelihood, QSSA error, CellOracle negative control. |

**Traps to avoid:** no 3-D cusp surface; no hairball of all states; no continuous colormap under a three-class map; no fake error bars on deterministic output; no p-value where there is no sampling; no truncated percentage axis; no bimodality claimed from a KDE shape; no over-precise numbers (two significant figures unless you can defend more).

---

# PART 6 — HOW TO PRESENT IT

**Category: Computational Biology & Bioinformatics (CBIO).** Judges are assigned by subcategory and CBIO explicitly covers mathematical modelling and simulation; TMED's category text never mentions computational work.

*Evidence worth knowing, including the part that cuts against a pure-mechanism framing:* across ISEF 2024 and 2025 full award lists, **zero CBIO or TMED First Awards went to a mechanistic dynamical-systems model.** Winners were biomarker discovery, virtual screens, and named deliverables — PT150, CtBP2, OncoNote, a nanosphere cocktail. Every one leaves the judge holding a *thing*. The confound is real (absence of mechanistic winners is confounded with absence of mechanistic entrants), but the safe conclusion is: **lead with a validated prediction wrapped in a mechanism, delivering a named, concrete recommendation.** That is exactly what R1–R3 produce.

**Frame the wet lab as computation-that-designed-an-experiment, never as validation.** If the result is negative and you claimed validation, your conclusion is deleted and you have nothing for the remaining eight minutes. If you claimed design, a negative result is *data about the model* — a finding. The recovery answer is strong: *"Correct, we don't know if it's true. Here is the experiment that would tell us, here is what I pre-committed to concluding from each outcome, and here is what happened when I ran it once."*

**Poster structure that survives a failed experiment:**
- A **dated prediction box, placed left of the results box**, saying "registered [date], before any wet-lab work." Physically separating prediction from outcome with a date is the highest-credibility move available and almost nobody does it.
- A **pre-specified interpretation table**: three rows saying what each possible outcome means. Written *before* the run, it kills every post-hoc-rationalization question.
- **Effect size and pre-registered direction. Never a p-value.** One shot, transformed line — "directionally consistent, n=1, not powered for inference, here's the effect size" is a *stronger* answer than a p-value, and a good judge knows it.
- **Bound what the system can falsify.** Ordering and relative timing: yes. Absolute dose: no. Saying this unprompted is worth more than any positive result.

**The one sentence:**

> **Reversing the earliest step of pancreatic cancer already works — and already fails, because it relapses within a week of stopping the drug. I show the relapse is a property of a broken two-node transcriptional loop, that closing it with a minimal co-formulated mRNA payload makes reversal self-sustaining, and I predict the composition, ratio and redosing interval that does it.**

---

## Sources for Part 0

Collins 2014, *Gastroenterology* 146:822–834, PMID 24315826 · Collins 2012, *J Clin Invest*, PMID 22232209 · Jiang 2023, *Gastro Hep Adv* 2:532–543, PMID 37425649 · Krah 2015, *eLife*, PMID 26151762 · Krah 2019, *Dev Cell*, PMID 31422917 · Letsou & Cai 2016, *PLoS Comput Biol*, PMID 27560383 · Rees 2019, *Nat Commun*, PMID 31138801 · Dobrowolski 2022, *Nat Nanotechnol*, PMID 35768613 · Joung 2023, *Cell*, PMID 36608654 · Ng 2021, *Nat Biotechnol*, PMID 33257861 · Cahan 2014, *Cell*, PMID 25126793 · Jin & Xiang 2019, *Cell Mol Life Sci*, PMID 30470852 · Halbrook 2017, *Cell Mol Gastroenterol Hepatol*, PMID 28090569 · Kamimoto 2023, *Nature*, PMID 36755098 · Ahlmann-Eltze 2025, *Nat Methods*, PMID 40759747 · Masui 2008, *Mol Cell Biol*, PMID 18606784.

**No evidence found** after targeted search: any measured Kd for an ID3–E-protein interaction in pancreatic acinar cells; any measured PTF1A half-life; any prior ODE model of ADM; any head-to-head kinetic comparison establishing the MEKi-versus-PTF1A asymmetry; any application of multi-objective Bayesian optimization to TF cocktail selection.
