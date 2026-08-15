# 012 · Reframe from reversal-discovery to durability-and-delivery

*Date:* 2026-08-15 · *Status:* accepted
*Supersedes in whole or in part:* 001 (screen arm), 002 (`W`), 003 (composable
topologies), 006 (two-target titration as a Stage 2 discriminator), 010 (NR5A2
placement), and the Stage-6 precondition section of 011.

**Process note.** `CLAUDE.md` rule 4 calls for a 2-subagent adversarial panel on
architectural decisions. This one was already decided by a three-reviewer panel
whose findings are recorded in Part 0 of `docs/ADM_MASTER_PLAN_v3.md`; this file
records that outcome rather than re-litigating it. Luqmaan directed that no panel
be convened.

---

## Question

The project was framed as *"reverse KRAS-driven ADM while keeping the cells
alive"*, executed as a seven-stage pipeline whose scientific core was a five-way
topology competition over an 11–13 state ODE, with payload identity to be
recovered by an in-silico perturbation screen.

Three reviewers, independently, found that framing unsurvivable. **Is the project
rebuilt around a different claim, and if so what happens to the model, the screen,
and the code already written against the old plan?**

## Positions considered

**Position A — keep the reversal framing, patch the weak points.**

The seven stages were sequenced deliberately, ~2,000 lines of decision documents
and ~800 lines of tested code already exist, and the model is genuinely more
complete than anything published on ADM — there *is* no prior ODE model of ADM. A
five-way topology competition is real model selection and it is the kind of thing
that reads as serious computational work. The specific objections are patchable:
widen the kinetic-asymmetry interval and report it honestly, rename the screen from
"discovery" to "prioritization", state the circularity as a limitation. Rebuilding
costs weeks against a schedule with a wet lab in the same cycle, and throws away
work that is correct.

**Position B — rebuild around durability and delivery.**

Each of the three objections is structural, not cosmetic:

1. **The claim was occupied.** Collins 2014 (PMID 24315826) reverts established
   PanIN in ~3 days with a generic MEK inhibitor, in vivo, no TF payload — *and*
   reports that lesions resumed when the drug was withdrawn. A judge who knows the
   paper deletes the project's headline in one question. No amount of framing fixes
   being eleven years late.
2. **The discriminator does not exist.** Stage 2 — four weeks, called "the
   scientific core" — was discriminated entirely by a 10× kinetic asymmetry between
   MEKi reversion and forced-PTF1A reversion. **Krah 2019 has no timepoint between
   24 hours and 3 weeks.** The measurement bounds the ratio at (1×, 21×); 10 was
   chosen from inside that interval. Collins 2014 alone spans ~2× between iKras\*
   and KC at equal pERK suppression. A model-selection step resting on a number
   nobody measured cannot be the core, and with ~30 free parameters a Q-value under
   a shared sampling box ranks prior volume rather than evidence.
3. **The identity claim was already published and the screen was circular.** Jiang
   2023 (PMID 37425649) names *"transcription factors with reduced activity (PTF1A,
   RBPJL, and BHLHA15)"* — BHLHA15 **is** MIST1, which fills the third slot. And
   regulon activity is scored from target-gene expression, where the PTF1A regulon's
   targets *are* the digestive-enzyme genes that define the acinar cluster: the
   regulon best separating acinar from metaplastic is PTF1A by construction. That is
   a harder circularity than the ODE-representability one already logged in 001.

Against this, the delivery layer — absent from the old plan entirely — contains the
strongest available result, and it is one neither a screen nor a 13-state ODE could
produce. Per-cell LNP dose is a **distribution**, not a number: Rees 2019 (PMID
31138801) shows the dose-per-endosome distribution is independent of administered
dose, so per-cell dose is compound-Poisson and heterogeneity is intrinsic;
Dobrowolski 2022 (PMID 35768613) measures systematically high- and low-uptake
subtypes. Combine that with this project's own mechanism, which is **stoichiometric
and thresholded** — PTF1A needs an E-protein partner and RBPJL cannot bootstrap, so
the pair is required together *in the same cell* — and co-formulation stops being an
optimization and becomes a requirement, with a quantifiable gap (~30% vs ~9% at a
30% marginal).

*Where they actually disagree:* not on whether the old model is well built — both
agree it is — but on **whether a well-built model of an occupied question beats a
smaller model of an unoccupied one.**

## Decision

**Position B. Rebuild around durability and delivery.**

The claim becomes: *trametinib reverts, it does not stick; a payload that re-closes
the PTF1A↔RBPJL loop should convert a drug-dependent reversion into a
self-sustaining one; this project finds the minimum payload, ratio, formulation and
schedule that does it.* Everything is now in service of R1 formulation, R2
composition, R3 durability.

Consequences, stated explicitly so nothing is reversed by silence:

| Artefact | Fate | Why |
|---|---|---|
| Seven-stage pipeline | **replaced** by eight phases (v3 Part 2) | — |
| `src/topology.py`, five-way competition, decision 003 | **dead** | the discriminator is unmeasured |
| 11–13 state ODE, `src/model.py` | **dead** as written; rebuilt at 3–4 states | cannot rank a candidate it does not contain |
| `W` state, decision 002 | **dead** — pERK is an input, not a state | see below |
| `src/binding.py`, `docs/derivations/binding_polynomial.md`, decision 006 | **demoted to supplementary**, not deleted | see below |
| `payload_subsets()` (2^k) | **dead** → mixture-amount design on the simplex | mass budget makes this a mixture problem |
| Stage 6 Lie-bracket precondition in 011 | **dead**; pulse forcing itself survives | precedent exists (PMID 27560383) |
| Screen as second arm, decision 001 | **demoted** to prioritization + a declared CellOracle negative control | identity is published |
| Pareto front as deliverable | **dead** → maximize durable reversal s.t. viability ≥ X, report therapeutic index | preference-free by construction |
| U-shaped hazard, decision 008 | **retired** — viability becomes a bench-measured floor + a CHOP output flag | see resolution 3 below |

**Two of these deserve their own sentence, because they are being reversed by
architecture change rather than by their own stated reversal conditions.**

**`W` (decision 002).** 002's primary reversal condition was a Stage 1 finding that
the ADM attractor sits too far from the separatrix for any plausible phospho-MEK
impulse to matter. **That test was never run.** `W` is dying instead because v3
specifies pERK as an *input*, which leaves no place for it. The cost is real and
should be stated rather than absorbed: the trametinib-vs-PD325901 withdrawal
asymmetry was 002's entire justification, and dropping `W` drops that prediction.
Note this is a *different* asymmetry from the 10× MEKi-vs-PTF1A one that Part 0
kills — Part 0 does not address it. **The honest position is that the drug-identity
prediction is set aside for schedule reasons, not refuted.**

**The binding polynomial (decision 006).** Its *purpose* — discriminating T1 from
T2 — is gone with the topology competition. Its *content* is not. The derivation's
result that `n_eff` scales as `√(E_tot/Kd)` in the tight regime and saturates in the
loose one is a statement about **how sharp the P–R bootstrap threshold is**, and v3
Phase 3 step 2 convolves the per-cell dose distribution against exactly that
threshold. A soft threshold and an ultrasensitive one give different
converted-fraction curves and therefore a different co-formulation gap — which is
the headline number. **Keep the derivation as supplementary and re-read it when
Phase 3 sets the threshold's sharpness.** v3 specifies only `E_free = E_total −
k·ID3`, a linear form; whether that is sufficient for Phase 3 is an open question,
not a settled one.

## Resolutions — the six open items, decided 2026-08-15

The audit above left six things unresolved. All six are now decided by Luqmaan.
v3 is amended in place (v3.1); the affected decision files carry their own
amendments.

**1 · `W` — pERK stays an input, but gains a withdrawal rebound profile.**
Neither "add `W` back" nor "drop it". v3's spine is durability *after
withdrawal*, and pERK rebound is a withdrawal-specific mechanism — deleting it
silently would have deleted the thing most relevant to the new endpoint. So pERK
stays an input and stops being a *step*: on withdrawal it follows a prescribed,
**swept rebound profile**. No state, no coupled algebraic loop, no identifiability
cost. The gain is that it converts an unmeasured state into a **measured input** —
a pERK western timecourse after washout, now Bench Handshake item 7. **The
trametinib-vs-PD325901 comparison is retired outright**, not set aside: it was a
free-prize prediction about an experiment nobody is running, and an unretired
claim is one someone builds on later. Decision 002 amended.

**2 · v3 Part 1.4 rewritten, narrower and stronger.** See the correction section
below, which is what prompted it. The limitation as written would have had this
project disclaiming published work on a poster. The residual gap is **one edge —
ERK→ID3 — not the whole node**, and it is testable with a single western
(Bench Handshake item 8). Phase 5's ordering prediction is therefore *better*
supported than v3 claimed. Say it that way.

**3 · Viability becomes a bench-measured floor; the hazard is retired.** Its six
parameters were never measurable, ER-stress dynamics were the least-constrained
part of the old model, and the bench measures viability directly. Modelling what
you can measure with parameters you cannot is backwards. **Kept as a flag, not a
term:** warn when predicted `P` falls below the CHOP-apoptosis threshold — the
low-PTF1A arm is real biology (Sakikubo 2018) that does not need a differential
equation. Decision 008 retired.

**4 · `E_free` uses the exact binding solution.** The observation that
`E_total − k·ID3` is the *tight-binding limit* of decision 006's equilibrium is a
save, not a problem: it makes the derivation v3's justification rather than dead
work. **Do not ship the linear form with a floor hack** — the negativity at
`k·ID3 > E_total` is the approximation announcing it has left its domain, and
clipping hides that rather than fixing it. Ship the exact form; comment that it
reduces to the linear one in the tight regime, with `n_eff` as the validity
diagnostic. Decision 006 promoted. **Follow-up correction (2026-08-15, on first
integration):** the diagnostic's prefactor is **0.5**, not the derivation's 1.34 —
that value belongs to the deleted two-target ternary complex, and reusing it would
have overstated threshold sharpness ~2.7× and inflated R1 with it. See 006.

**5 · The guard test is the flagship.** *"`dR/dt` has no P-independent term"* is
the bootstrap claim expressed as a zero, and violating it would run clean while
predicting the opposite — the same failure class the retired `raf_drive` test
guarded. It gets the same treatment: a **guard-on-the-guard** that constructs the
violating implementation and *requires* it to fail, because an assertion passing
for both the correct and incorrect version tests nothing. Plus a **dynamic**
companion: `R` must not rise from zero while `P` is held at zero, which catches an
accidental basal term that a structural check would miss.

**6 · The Pareto cut is an explicit reversal, recorded as one.** The previous plan
and decision 001 both called the reversal-vs-viability Pareto question "a headline
result and the stated objective". It is cut deliberately: **a front is
preference-free by construction, and a one-shot experiment consumes a decision,
not a set.** Replaced by constrained optimization — maximize durable reversal
subject to viability ≥ X — plus a therapeutic index. Recording it as a decision
rather than letting it read as an omission.

**On the arm budget.** Not a modelling choice — a **blocking external input** that
depends on the wet lab and is not known yet. It does not block Phases 2–3. Phase 0
carries an explicit placeholder range (12 / 24 / 48 wells), Phase 4's
mixture-amount design must be instantiable at any tier, and the deliverable
includes **what changes at each tier**, so the number can be argued for rather
than accepted.

## A correction to v3 itself

**v3 Part 1.4 says: *"The ID3→E47 titration node in pancreatic ADM is a modelling
hypothesis, not a documented mechanism. Targeted PubMed search returned no evidence
for it in this specific context."* That is too strong, and it should not go on a
poster as written.**

Checked against PubMed this session. Dufresne et al. 2010, *Int J Cancer*
129(2):295–306, PMID 20830706, [DOI 10.1002/ijc.25668](https://doi.org/10.1002/ijc.25668),
is titled *"Id3 modulates cellular localization of bHLH Ptf1-p48 protein"* and
reports, in **AR4-2J** — the same rat acinar line as the bench — that a proliferative
signal (gastrin) *"leads to increases in Id3 protein expression and levels of
Id3/E47 and Id3/Ptf1-p48 interactions, and a decrease in the level of E47/Ptf1-p48
interaction"*, that *"Id3 silencing reversed the cytoplasmic mislocalization"*, and
that Id3 is overexpressed with Ptf1-p48 mislocalized in human and murine pancreatic
preneoplastic lesions.

That is a documented two-target titration, in the right cell type, with a
loss-of-function reversal and an in-vivo correlate. **What is genuinely undocumented
is narrower and should be stated in exactly this form:**

- no **Kd** has been measured for ID3·E47 or ID3·PTF1A — correct, and unchanged;
- the driver in Dufresne is **gastrin**, not KRAS/ERK. **The ERK→ID3 edge is the
  modelling assumption**, not the titration node itself.

Stating a limitation that is not true is not conservative — it is a different error
in the same family, and a judge who knows Dufresne will ask why the supporting
paper was missed. The Phase 5 ordering prediction (trametinib-first should strictly
dominate) rests on the titration node and is therefore **better supported than v3
claims**, while remaining conditional on the ERK→ID3 edge.

## What would reverse this

1. **A measurement that pins the MEKi-vs-forced-PTF1A kinetic ratio.** This is the
   specific and most likely reverser. The five-way topology competition was cut for
   exactly one reason: its discriminator is unmeasured, bounded only at (1×, 21×) by
   Krah 2019's 24 h / 3 weeks / 6 weeks sampling. **If a paper, a dataset, or a bench
   experiment supplies a timepoint between 24 hours and 3 weeks — or any head-to-head
   kinetic comparison at matched pERK — the competition becomes runnable and should
   be reconsidered on its merits.** The code to run it is in git history at
   `4afa014` and does not need rewriting from scratch. Reconsider, do not
   auto-restore: the ~30-parameter Q-value objection (prior volume, not evidence)
   is independent of the discriminator and would still have to be answered.
2. **The delivery result collapses.** If the bench's measured uptake CV is small
   enough that correlated and independent per-cell doses give near-identical
   double-above-threshold fractions, R1 evaporates and the strongest reason to
   prefer this framing goes with it. Check this *first*, at Phase 3, before Phase 4
   or 5 is built on top of it.
3. **Gate B fails.** If the 3–4 state core cannot produce two stable states plus a
   saddle with an identifiable separatrix, there is no persistence window, R3 is
   gone, and the project is a delivery-and-composition study with no dynamics. That
   is still a project, but it is a different one and should be renamed rather than
   quietly continued.
4. **Phase 3 needs threshold sharpness the linear `E_free` form cannot express.**
   Then decision 006's polynomial comes back — as a Phase 3 input, not as a topology
   discriminator.
