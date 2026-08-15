# Reversing KRAS-driven acinar-to-ductal metaplasia without killing the cell
## In-silico master plan

> ## ⚠ STATUS — READ FIRST
>
> **This is the only current plan.** It supersedes every earlier draft.
>
> If you find `INSILICO_PLAN_8months.md` or `INSILICO_PLAN_v2.md` anywhere on this machine: **do not read them, do not build from them, do not merge them in.** `INSILICO_PLAN_8months.md` is the pre-panel draft and it specifies a model that a five-specialist review demolished — an 11-state system with fast acetylation as a state variable, a one-sided toxicity ceiling, an AlphaFold-derived rate constant, and a Perturb-seq module built on cell lines that do not express the genes. Building from it wastes weeks and reintroduces errors that were already found and fixed.
>
> **This document does not contain a fully-closed ODE system, and that is deliberate.** §3 specifies the architecture, the state list, the load-bearing equation forms, and the reasoning behind each. **Writing the complete right-hand side is Stage 0's job** — see §3.4 for exactly what is specified here versus what you must derive. Do not go looking for the "real" equations in another file. There isn't one. Derive them.

---

# PART 1 · THE SYSTEM

| | |
|---|---|
| **Cells** | AR42J, **ATCC parental (CRL-1492)** — *not* the B13 subclone |
| **ADM induction** | **KRAS G12D**, inducible (Tet-On), stably integrated |
| **Brake** | **Trametinib** (MEK1/2 inhibitor) |
| **Payload** | **3 synthetic mRNAs**, LNP-delivered — see §1.3 |
| **Format** | Growth-factor-reduced Matrigel, 3D, daily media change |
| **Wet lab** | Same cycle as this modeling. **One shot, no iteration.** |
| **Objective** | Reversal **and** viability. Reversal > viability, but both is the goal. |

## 1.1 Why KRAS induction is the right driver

Moving from dex + TGF-α/EGF to inducible KRAS G12D fixes three structural problems at once:

1. **The model system now matches every validation dataset.** Collins 2014 is iKras\* (inducible KRAS G12D). Krah 2015 and Krah 2019 are KRAS-driven. A dex+EGF wet lab would use a *different induction route* than every paper the model is validated against.

2. **The forcing input becomes an experimental dial.** In the ODE, ERK activity `K` was an imposed forcing function. With Tet-On KRAS it is a knob in the dish: dox on, off, titrated. That makes **hysteresis directly testable** — up-sweep and down-sweep the KRAS dose and look for a gap. One sweep direction cannot distinguish a switch from a steep graded response; two can. This is the single most important experiment the model can design.

3. **The history experiment becomes doable in vitro.** Collins' strongest datum is that iKras\* and KC mice revert on *different* timescales under identical MEK inhibition with pERK suppressed equally in both — differing only in how long KRAS had been on. With inducible KRAS you reproduce that in a dish: short KRAS exposure vs long, then identical trametinib.

It also **removes a live hazard.** In AR42J-B13, dex *inhibits* amylase and drives hepatocyte transdifferentiation, and dex + EGF specifically enriches the ductal population via C/EBPβ (Al-Adsani 2010, PMID 21048969) — a hepatic third attractor the state space cannot represent.

## 1.2 Trametinib is not interchangeable with PD325901 — model the difference

Collins used PD325901. Trametinib is the clinical drug and the right choice, but it is mechanistically distinct in a way that matters to the equations:

- **PD325901** — allosteric MEK kinase inhibitor; blocks MEK's catalytic activity.
- **Trametinib** — additionally inhibits **RAF-mediated phosphorylation of MEK**, blocking MEK *activation*, not just MEK output. This suppresses the adaptive feedback reactivation that limits other MEK inhibitors.

**Model these as two terms, not one.** A pure activity inhibitor lets the RAF→MEK feedback wind up during treatment, so ERK rebounds fast on withdrawal. An activation inhibitor does not.

> **Free prediction — conditional form.**
>
> *"Withdrawal asymmetry between a MEK-catalytic inhibitor and a MEK-activation inhibitor occurs **if and only if** the transient phospho-MEK overshoot is sufficient to return the state across the separatrix. We report the fraction of the plausible parameter ensemble in which that condition is met, and the overshoot magnitude required."*
>
> **Why it must be stated conditionally.** The earlier version — *"trametinib shows a longer relapse delay than PD325901 at matched pERK suppression"* — is not derivable from a static `K_eff = K·f_act(v)·f_cat(v)`. At matched pERK the states see an identical `K_eff`, evolve identically, and reach withdrawal at the same point; a static product has nowhere to store the difference. The asymmetry requires a **state** (`W`, the phospho-MEK pool, §3.2).
>
> **And `W` is fast** — MEK phosphorylation turns over in minutes against a model running hours to weeks — so within minutes of washout both drugs converge to the same `W`, while Collins measures relapse at 7 days. The naive reading is that the asymmetry is invisible.
>
> **It survives only through bistability.** The overshoot does not need to *persist*; it needs to be large enough to carry the slow state back across the separatrix. Once crossed, relapse proceeds regardless of what `W` does afterwards — **a transient in a bistable system can have a permanent consequence.** That is a real mechanism and a demanding condition, and the honest deliverable is the ensemble fraction, not the claim. **If the condition holds in 5% of the ensemble, report 5% and say what it means.**
>
> Note the quantity that decides it is the **impulse** `∫ΔK_eff dt`, not the peak: a minutes-long spike delivers little impulse to a process running on days. Whether it is enough depends strongly on `τ_W`, which is why `τ_W` is sampled across minutes-to-hours rather than fixed (§3.2).

## 1.3 The payload — three mRNAs, and why mRNA sharpens the control problem

| | mRNA | Role | Why it is in the payload |
|---|---|---|---|
| **u₁** | **RBPJL** | Rebuilds the mature partner | *Rbpjl* has **no PTF1A-independent promoter.** Its only driver is PTF1 itself. It cannot bootstrap. This is the deepest hole in the loop and the reason ERK release alone is slow. |
| **u₂** | **PTF1A** | The master factor | Re-expression reverts established lesions with cells surviving (Krah 2019). Supplied alone it is slow, and it is a sequestration target. |
| **u₃** | **E47 / TCF3** *(leading candidate)* | Restores the limiting subunit | ID3 traps **both** E47 and PTF1A, and **E47 carries nuclear import** of the complex (Dufresne 2010, PMID 20830706). A stoichiometric trap is not relieved by adding the trapped species — it is relieved by saturating the trap. You cannot deliver an shRNA as an mRNA, so the mRNA-native answer to ID3 is **flooding the E-protein pool**, not knocking ID3 down. |

**The third slot is a model output, not an assumption.** Stage 3 ranks candidates for u₃ — E47/TCF3, NR5A2, MIST1/BHLHA15 — by how much each moves the reversal boundary per unit dose across the parameter ensemble. If the model nominates NR5A2 over E47, that is a result, not a correction. Keep the slot parameterized in code.

> ### ⚠ The circularity limit on what this model may claim
>
> RBPJL, PTF1A and the E-protein are **state variables**. A model whose states are the payload species cannot be asked *"which molecules should be supplied?"* — it can only nominate what was written into it, and the answer is fixed by the model's construction rather than by the biology. Ranking a **pre-specified** three-candidate list for u₃ narrows this but does not escape it, because the list was curated by the same person who chose the states.
>
> **The mechanistic model is therefore asked only questions it can answer without circularity:**
>
> 1. **Necessity — which components are actually required.** Four interventions (trametinib, u₁, u₂, u₃) give **16 subsets**. Run all 16 across the parameter ensemble; report which components are necessary, which are redundant, and the marginal contribution of each. **"Three mRNAs is over-engineering" is a valid and useful outcome** — it is cheaper at the bench and it is a real finding, not a failure.
> 2. **Dose, schedule and order** — Stages 5 and 6, unchanged.
> 3. **Whether the reversal-optimal payload and the viability-optimal payload are the same payload.** This is a Pareto question, not a single optimum. If the two differ, that is a headline result and it is exactly the stated objective — reversal **and** viability, not either alone. If they coincide, say so; it is a cleaner recommendation for the bench.
>
> **Molecule identity comes from the other arm** — the unbiased screen in **Stage 3B**, which starts from a candidate set the modeller did not curate and is capable of ranking PTF1A, RBPJL and NR5A2 low. See `docs/decisions/001-two-arm-payload-derivation.md`.
>
> **Fair-comparison rule for the necessity analysis.** Dropping a component changes total delivered material, so a naive subset comparison rediscovers *"more protein is better."* Report each subset **twice**: once at matched per-component dose, and once at matched **total** delivered dose (∫Σuᵢdt held constant, redistributed over the surviving components). A component is only "necessary" if it survives the matched-total comparison. This is the same discipline Stage 6 imposes on ordering, for the same reason.

**Why mRNA makes the control problem sharper, not softer.** With a viral vector, expression is a step of unknown height and indefinite duration. With mRNA the input is a **pulse of known shape**: uptake, translation ramp, first-order decay with a measurable half-life. Therefore:

- **"Duration" is not a free variable.** It is set by mRNA half-life plus redosing interval. The therapeutic window is not an abstract (dose × duration) plane — it is **(dose per pulse) × (redosing interval)**, a schedule someone can execute.
- Inputs are bounded and non-negative by construction: `uᵢ(t) ≥ 0`. You can add mRNA; you cannot subtract it. Exactly the constraint structure Stage 5 assumes.
- Nothing integrates, so relapse after payload clearance is a **real prediction the model must make** rather than an artifact of constitutive expression.

**Construct note:** if any two are bicistronic, keep **RBPJL translated ahead of PTF1A** so a residual 2A tag does not land on PTF1A's C-terminal RBPJ-contact motif (W301, inside the disordered 305–328 region). Three separate mRNAs makes this moot — an argument for three separate mRNAs.

---

# PART 2 · THE PATHWAY COLLAPSE

The full map is 228 nodes and 465 edges (`adm_nodes.csv`, `adm_edges.csv`). Every arm — MAPK, PI3K, RAL/RAC, TGF-β, NOTCH, WNT, HH, HIPPO, inflammatory, ECM — converges on **one broken positive-feedback loop**.

## 2.1 The healthy loop

The functional unit is a **trimer**:

```
PTF1-L  =  PTF1A  +  E-protein (E47/HEB)  +  RBPJL     ← adult, Notch-refractory
PTF1-J  =  PTF1A  +  E-protein (E47/HEB)  +  RBPJ      ← immature, NICD-vulnerable
```

PTF1-L drives **two** targets that matter (Masui 2008, PMID 18606784):

1. **Its own gene**, via a 2.3-kb enhancer carrying two conserved PTF1 sites.
2. **Its own partner**, *Rbpjl*.

That is a **dual autoregulatory loop** — self-sustaining while closed. A separate, PTF1-independent **13.4-kb proximal "ignition" promoter** is the only route by which PTF1A can be made when the loop is open, and it is what ERK suppresses.

## 2.2 What KRAS G12D does

Gly12 sits within van der Waals distance of **both** Arg789 of p120GAP and Gln61 (Scheffzek 1997, PMID 9219684). The substitution is **steric, not electrostatic** — even alanine at position 12 disturbs the transition state. The arginine finger cannot be positioned, GAP-accelerated hydrolysis fails, RAS stays GTP-loaded. Constitutive ERK.

## 2.3 The four cuts

ERK does not break the loop in one place. It cuts it four ways, and that redundancy is why single-agent reversal is slow:

1. **Ignition suppressed** — ERK shuts down the 13.4-kb PTF1-independent promoter. No new PTF1A from scratch.
2. **RBPJL starved** — with PTF1 output down, *Rbpjl* transcription falls. Having no independent promoter, it **cannot recover on its own**. The loop cannot re-close even if PTF1A returns.
3. **ID3 induced** — downstream of growth-factor/ERK signalling, ID3 rises and traps E47 *and* PTF1A in the cytoplasm. Subunits exist but are not in the nucleus.
4. **PTF1A degraded** — TRIP12 targets a degron at **K312 (human) = K311 (rat) = K309 (mouse)** — the same residue, second lysine of the conserved `PEDPRKLNSK` motif. The offset is a species-numbering artifact: human residue 311 is a serine.

Chromatin at acinar enhancers loses H3K27ac while **retaining H3K4me1** — "decommissioned but primed."

> **Correction that matters:** the retained H3K4me1 Falvo 2023 reports is at **metaplasia** genes. The published epigenetic memory primes the cell to *leave* the acinar state, not return to it. Oriented correctly, that memory explains **relapse** — a better use of it.

## 2.4 The one-line collapse

> **ADM is locked because a self-sustaining transcriptional loop was opened, and the one component that cannot bootstrap itself — RBPJL — has no way back in.**

This is why the therapy is trametinib **plus** payload rather than either alone. Trametinib lifts ERK suppression on ignition and lowers ID3; the mRNA supplies what cannot rebuild itself fast enough. Neither is sufficient. That is a mechanistic combination rationale, not "let's try both."

## 2.5 The observation that discriminates every model

- Remove ERK → **~50% acinar in 3 days** (Collins 2014).
- Force PTF1A directly → **zero acinar at 24 h; 3 weeks to redifferentiate** (Krah 2019).

**Forcing the master TF is an order of magnitude slower than removing the upstream kinase.** Any model in which PTF1A is the sole bottleneck predicts the opposite. Not a fitting problem — a statement about architecture, and the observable Stage 2 uses to choose between architectures.

---

# PART 3 · THE MODEL

## 3.1 States

Write the full system, then **nondimensionalize and adiabatically eliminate the fast variables** down to 5–6 slow states. Not simplification for its own sake: it is what makes continuation tractable, identifiability interpretable, and reachability computable.

**Not every state exists in every topology.** The model is **composable**: a CORE right-hand side plus optional terms switched on by config flag. See §3.5 — this is what keeps the slow count at 6 despite the state list growing.

### CORE states — present in every topology

| Symbol | Meaning | Timescale |
|---|---|---|
| `P_n` | PTF1A, nuclear **free** (ID3-bound and complexed pools are derived) | slow |
| `R` | RBPJL protein | **slow — the bottleneck** |
| `E_tot` | **total** E-protein (E47/HEB) — synthesis, degradation, and the `u₃=E47` input | **slow** |
| `I` | ID3 | intermediate → leading QSS candidate at reduction |
| `M` | slow self-reinforcing chromatin at **metaplasia** loci | **slow — the memory** |
| `A` | acinar output (amylase / CPA1 proxy) | slow |
| `S` | secretory cargo / capacity ratio | slow |
| `W` | phospho-MEK pool / RAF→MEK activation drive | **fast — PROTECTED FROM ELIMINATION, §3.2** |

### VARIANT states — exist only in the topologies that test them

| Symbol | Meaning | Present in |
|---|---|---|
| `MIST1` | secretory capacity arm (BHLHA15) | **T5 only** |
| `NR5A2` | enhancer co-activator (or acinar-output co-activator — competing variants) | NR5A2 topology variants only |

### Derived, never integrated — recovered algebraically from the binding polynomial

`P_c` (ID3·PTF1A) · `E_free` · `C_L` (PTF1-L) · `C_J` (PTF1-J)

> **`P_c` is no longer a state and its sequestration did not disappear.** PTF1A partitions into free and ID3-bound through the **same competitive equilibrium** that partitions E — one polynomial, two targets (§3.2). Deleting `−k_seq·I·P_n` without moving sequestration into the equilibrium would swap one error for its opposite.

**Slow-state arithmetic, since state count is the schedule risk.** CORE is 8 states. At reduction: `C_L`, `C_J`, `E_free`, `P_c` are already algebraic; `I` is the QSS candidate; `W` is retained by exemption. That leaves **6 slow states — `P_n`, `R`, `E_tot`, `M`, `A`, `S` — plus `W`.** Within the stated 5–6 target, with `W` as the one declared exception. Variant states are absent from Stage 1's continuation, which runs on CORE with the payload at zero.

**Forcing input:** `K` = ERK activity, set by KRAS-G12D induction level.
**Control inputs:** `u₁` RBPJL mRNA, `u₂` PTF1A mRNA, `u₃` third mRNA, `v` trametinib.

## 3.2 The equation forms that carry the argument

**PTF1A — two sources, one PTF1-independent**
```
dP_n/dt = α_ign·g(K_eff) + α_auto·Hill(C_L, M; n) − δ_P·P_n + u₂(t)
```
`α_ign·g(K_eff)` is the 13.4-kb ignition, suppressed by ERK, released by trametinib — note it takes `K_eff`, which already carries the drug through `W` and `f_cat`, so writing `g(K,v)` here would double-count it. `α_auto·Hill(·)` is the 2.3-kb autoregulatory enhancer.

> ⚠ **The term `− k_seq·I·P_n` was here and has been deleted. Do not restore it.** It is a first-order sink, which §3.4's constraint 2 forbids and which this section's own titration paragraph calls out as generating no ultrasensitivity. **Sequestration has not been removed — it has been moved into the binding polynomial below**, where PTF1A and E partition through one shared competitive equilibrium. `P_n` here is *free* nuclear PTF1A; the ID3-bound and complexed pools are derived, not integrated.

**RBPJL — the asymmetry that makes the system interesting**
```
dR/dt = β·Hill(C_L; m) − δ_R·R + u₁(t)
```
**No ignition term.** *Rbpjl* has no PTF1-independent promoter. This asymmetry is the structural reason RBPJL is the deepest hole, and it falls out of the equations rather than being asserted.

**ID3 sequestration — TWO-TARGET competitive titration, not a first-order sink**

ID3 traps **both** partners of an obligate heterodimer. Verified: Dufresne 2010 (*Int J Cancer* 129(2):295–306, PMID 20830706) reports gastrin raising Id3 and increasing **both** Id3/E47 **and** Id3/Ptf1-p48 interactions while *decreasing* E47/Ptf1-p48 — in AR4-2J, the wet-lab line. Silencing Id3 reversed the mislocalisation.

One conserved equilibrium, three complexes, solved simultaneously:

```
E_tot   = E_free + [ID3·E] + (E bound in C_L and C_J)
P_tot,n = P_n    + [ID3·P] + (PTF1A bound in C_L and C_J)
I_tot   = I_free + [ID3·E] + [ID3·P]
```

**Why this is the term that decides whether the model can be bistable.** A first-order sink `−k·I·P_n` generates **no** ultrasensitivity. Molecular titration of a stoichiometric partner is a classical ultrasensitivity generator and can produce switching without high Hill coefficients — and a titrator that sequesters **both** members of an obligate heterodimer gives a **sharper** threshold than one sequestering either alone, because both routes to complex formation are shut simultaneously.

**The full algebra, every step, is in `docs/derivations/binding_polynomial.md` — not only in code.** This term is load-bearing for the central claim and must be followable line by line.

**Chromatin as slow, self-reinforcing memory at metaplasia loci**
```
dM/dt = k_w·φ(K, C_L) − k_e·M + ε·M²/(θ² + M²)
```
> **Do not use `dH/dt = k_write·C_L − k_erase·H`.** Acetylation turnover at p300 sites is **< 30 min** (Weinert 2018, PMID 29804834), so H is adiabatically slaved and `Hill(C_L, H)` collapses to `Hill(C_L, C_L)` — the variable algebraically vanishes. The memory must be slow, self-reinforcing, and located at **metaplasia** loci where the evidence puts it. **This is what sets relapse timing after payload clearance.**

**Trametinib — two terms acting at different points, and one protected state**

> ⚠ **The static product `K_eff = K·f_act(v)·f_cat(v)` was here and has been deleted. Do not restore it.** It cannot produce the §1.2 prediction: at matched pERK the two drugs give identical `K_eff`, so every state evolves identically and withdrawal is identical. The two terms must act at *different points in the cascade*, which requires a state.

```
dW/dt  = k_on · RAF_drive(K_eff) · f_act(v)  −  k_off · W        ← f_act blocks FORMATION
K_eff  = K · W · f_cat(v)                                         ← f_cat blocks OUTPUT
```

`W` is the phospho-MEK pool. `f_cat` = MEK catalytic inhibition, shared with PD325901. `f_act` = inhibition of RAF-mediated MEK phosphorylation, trametinib-specific. **Setting `f_act ≡ 1` still recovers PD325901**, so the direct comparison to Collins is preserved.

**`RAF_drive` must be a *decreasing* function of `K_eff`.** Phospho-MEK accumulates under a catalytic inhibitor *because* falling ERK output relieves ERK-mediated negative feedback on RAF, so RAF drive rises as ERK falls. **Implementing it as increasing inverts the mechanism and the model predicts the opposite result** — silently, since the code still runs. The sign is asserted in the code and covered by a unit test that fails on an inverted implementation.

**`W` is PROTECTED FROM ELIMINATION.** It is fast by turnover (minutes) and the Stage 0 reduction would ordinarily eliminate it — but QSS on `W` substitutes `W_ss ∝ RAF_drive(K_eff)·f_act(v)/k_off` straight back into `K_eff` and **recovers the static product**, destroying the prediction it was added to support. It is the one place in this model where a fast variable is deliberately retained. It is exempt from the fast-variable sweep, and the writeup must say why. Reasoning: `docs/decisions/002-w-state-protected-from-elimination.md`.

**`τ_W = 1/k_off` is sampled, not fixed** — minutes (pure phospho-MEK turnover) to hours (the feedback-relief components, DUSP/SPRY, are transcriptional). The overshoot's *duration* sets the impulse `∫ΔK_eff dt`, which is what decides whether the separatrix is crossed, so `τ_W` is one of the parameters that determines the §1.2 ensemble fraction. Fixing it would presuppose the answer.

**Viability — a U-shaped hazard integrated along the trajectory, not a ceiling**
```
dS/dt   = φ(A) − γ(capacity)·S
h(S,P_n) = h_high(S) + h_low(P_n)                     ← rises at BOTH ends
survival(t) = exp( −∫₀ᵗ h(S(τ), P_n(τ)) dτ )
```

> ⚠ **"death when `S` exceeds tolerance" was here and has been deleted. Do not restore it.** A one-sided threshold on `S` is the `U_crit` construct the panel killed, in new clothing.

`h_high(S)` — cargo outrunning capacity. `h_low(P_n)` — CHOP-dependent apoptosis under PTF1A loss (Sakikubo 2018, PMID 30361559; Backx 2021, PMID 33762742).

**Three reasons the hazard is right and not merely a smoothed threshold:**
1. It is **genuinely U-shaped**, which the plan asserts and the deleted form contradicted.
2. It **integrates along the trajectory**, so a brief excursion costs survival without automatically killing — matching the panel's finding that the risk is a **rate** mismatch during the transient, not an amplitude ceiling.
3. It yields a **continuous 0–1 viability number**, which is what the y-axis of the reversal–viability figure actually needs. A binary alive/dead flag cannot be plotted against a continuous reversal axis.

`γ(capacity)` is carried by the **`MIST1` state in T5** and by a constant elsewhere — see §3.5. Folding capacity into `γ(A)` was considered and rejected: it makes T5 unimplementable and rigs Stage 3's ranking of MIST1 as u₃.
> **Do not use a one-sided `U_crit`.** Secretory *capacity* — chaperones, ER membrane, XBP1, vesicle trafficking — is co-induced by the same differentiation program via MIST1 (Jakubison 2018, PMID 29719936). At steady state the cargo/capacity ratio is invariant and a ceiling does not exist. In every published experiment where full redifferentiation occurred, the cells survived.

The defensible version is a **kinetic mismatch during the transient**: cargo ramps faster than capacity, so the risk is a *rate* limit, not an *amplitude* limit. And the curve is **U-shaped** — too little PTF1A also kills, via CHOP-dependent apoptosis (Sakikubo 2018, PMID 30361559), and blocking dedifferentiation under stress increases death (Backx 2021, PMID 33762742).

## 3.3 The headline claim

> **Published interventions that reverse ADM occupy two distinct positions in a reversal–viability plane, and the difference is predictable from how many survival branches the intervention severs. The model predicts where a trametinib + mRNA combination lands, on what schedule — and which components of it are actually necessary.**

*(Deliberately "a combination" rather than "the 3-mRNA combination": the necessity analysis in §1.3 is allowed to conclude that fewer components suffice, and the headline must not presuppose the answer to a question the project is asking.)*

Better than "reversal and viability are two boundaries of one window" because it **explains an existing published discrepancy** (§6, HO-2) instead of asserting a geometry.

## 3.5 Composable topologies — one right-hand side, not five files

**The problem this solves.** With `W`, `E_tot`, `MIST1` and `NR5A2` the state list grows from 10 to ~13. Stage 0 targets **5–6 slow states**, and continuation on a large system is the single biggest schedule risk in this project. Adding four states without addressing that is how the plan quietly becomes unrunnable.

**The resolution: not every state exists in every topology.**

- **CORE** — the loop itself: `P_n`, `R`, `E_tot`, `I`, `M`, `A`, `S`, `W`. This is what **Stage 1 continues on**, with the payload at zero. Reduces to **6 slow + `W`**.
- **VARIANT** — `MIST1` exists only in **T5**. `NR5A2` exists only in the variants that test it. Switched on for the stages that need them, absent otherwise.

**Implementation rule, binding on all later stages:** a **base right-hand side plus optional terms selected by a config flag** — *never* five copied files that drift apart. Stage 2 must swap topologies with **one argument**; Stage 3's sensitivity analysis must run over whichever state set is active; Stage 5's necessity analysis must toggle each of the four interventions independently and swap `u₃`'s identity.

A topology is therefore a **configuration object**, not a source file. The five Stage 2 candidates (T1–T5) differ by which optional terms are active, and the sampler, integrator and scoring code are shared by construction — which is exactly what the Ma et al. 2009 Q-value methodology requires, since an unequal comparison across topologies is not model selection.

Reasoning and reversal conditions: `docs/decisions/003-composable-topology-architecture.md`.

## 3.4 ⚠ SPECIFIED HERE vs. DERIVED IN STAGE 0

This is the division of labor. Read it before deciding something is "missing."

**Specified in this document — treat as given, do not re-litigate:**
- The full state list and which states are fast vs slow (§3.1)
- Which terms exist in each equation and **why** each is there (§3.2)
- The four structural constraints that make the model non-generic:
  1. `dR/dt` has **no** ignition term
  2. ID3 acts by **titration**, not as a first-order sink
  3. chromatin memory is **slow and at metaplasia loci**
  4. trametinib is **two terms**, not one
- Viability is a **U-shaped ratio**, never a one-sided ceiling
- Which inputs are actuatable, and that all are non-negative and bounded

**Derived in Stage 0 — your job, and the reason Stage 0 exists:**
- Explicit functional forms for `Hill(·)`, `g(K,v)`, `φ(A)`, `f_act`, `f_cat` — choose them, justify each choice in a `docs/decisions/` file
- The complete closed right-hand side as runnable code
- The binding polynomial for E-protein titration, solved explicitly
- Nondimensionalization: which groups, which scales — expect ~30 parameters to collapse to ~15–18 dimensionless groups
- Quasi-steady-state elimination of `P_c`, `E`, `C_L`, `C_J` — and the algebra showing it is valid
- The parameter table: every parameter, its plausible range, and its source (or "unmeasured — sampled")
- Initial conditions for the ADM and acinar states
- Numerical settings: stiff solver choice (BDF/Radau), tolerances, scaling

**If a functional form is genuinely ambiguous, that is a `docs/decisions/` panel, not a blocker.** Pick one, document what would reverse the choice, and move. Stage 2 exists precisely to test whether such choices matter.

---

# PART 4 · THE PIPELINE

Seven stages. One argument. Each stage changes what the next stage *does*, not just what number it receives.

## STAGE 0 · Reduce and certify — 2 weeks

**Do.** Everything in §3.4's "derived" list. Then: multi-start Newton for all fixed points; classify by eigenvalue; locate the saddle; trace the separatrix by backward integration of its stable manifold.

**Proves.** The system has the structure being claimed — two stable states with the right identity, plus a saddle.

**Buys.** Everything downstream, plus tractability. At 5–6 states you can hand-roll pseudo-arclength continuation in ~200 lines of numpy. PyDSTool is unmaintained and painful on Windows; AUTO-07p is a Fortran build with a custom input format. **Continuation on the full-dimension system is where 6 weeks becomes 5 months.** Reduce first.

**Never cut this.** Assuming bistability from a cartoon and discovering monostability at week 20 is the classic failure.

## STAGE 1 · Two-parameter bifurcation — 2 weeks · ★ MONEY FIGURE

**Do.** Continuation in **KRAS dose × trametinib dose**. Locate saddle-node curves, cusp points, the hysteresis region. Only possible because KRAS is an inducible dial.

**Proves.** Whether a bistable region exists, where its boundaries lie, how wide the hysteresis gap is — in units both axes can be set to in a dish.

**Buys.** The figure judges remember, and **the experimental design for the wet lab**. The gap width in KRAS-dose units *is* the up-sweep/down-sweep experiment: which doses, how many points to resolve the gap.

## STAGE 2 · Topology competition — 4 weeks · ★ SCIENTIFIC CORE

**Do.** Five candidate architectures, all from the literature, all sampled from the **same** Latin hypercube box with the **same** code. Report the **Q-value** for each — fraction of sampled parameter sets achieving the target function.

| | Topology | Source |
|---|---|---|
| **T1** | ID3 as stoichiometric titration of the E-protein pool | Dufresne 2010 |
| **T2** | ID3 as a first-order sink on PTF1A | the naive version |
| **T3** | RBPJ→RBPJL handoff as the memory | Masui 2010, PMC2902682 |
| **T4** | Slow self-reinforcing chromatin at metaplasia loci | Falvo 2023 |
| **T5** | MIST1 parallel arm carrying secretory capacity | Jakubison 2018 |
| **T6a/T6b** | **NR5A2 as enhancer co-activator** (on `α_auto`) vs **NR5A2 as acinar-output co-activator** (on `A`) | Holmstrom 2011 — see §3.2 caveat |

**T6a vs T6b exists because NR5A2's placement is an assumption, not a finding.** Holmstrom 2011 (PMID 21852532, GSE34295) shows LRH-1/NR5A2 and PTF1-L co-regulating an **exocrine transcriptional network** — co-occupancy on digestive-enzyme genes. It does **not** demonstrate binding at the *Ptf1a* autoregulatory enhancer. Rather than pick one and hide the assumption, both placements enter the competition and the Q-values decide. This also gives Stage 3 a genuine comparison for the u₃ slot: **E47 helps by relieving titration, NR5A2 helps by boosting transcription — two different mechanisms competing for one payload slot.**

**Discriminator is §2.5:** MEKi reverts in 3 days, forced PTF1A takes 3 weeks. Topologies where PTF1A is the sole bottleneck **cannot** produce a 10× asymmetry. Topologies where the rate-limiting step is RBPJL accumulation, chromatin, or E-protein availability **can**.

**Proves.** Which architecture the published data supports.

**Buys.** This converts the weakest claim into the strongest. *"X% of parameter space is bistable"* for one topology is a property of your prior box — halve the box and the number changes. The **same** number across five topologies under an identical box is **model selection**, a published methodology: **Ma, Trusina, El-Samad, Lim & Tang, *Cell* 2009, PMID 19703401.** Same sampler, same integrator, swap the right-hand side.

**Mandatory companion:** prior-sensitivity. Halve and double the box; show the *ranking* is invariant even though absolute Q-values move.

> **Free prize.** Krah 2019 and Jakubison 2018 disagree on whether Panc1 responds to PTF1A, and Krah's printed explanation is that Jakubison *"may have achieved higher level expression of Ptf1a, overcoming an inhibitory threshold."* A published, unresolved, **threshold-shaped** discrepancy. If topology selection identifies which architecture predicts a sharp PTF1A dose threshold, and at what fold-change, the model reconciles two papers.

## STAGE 3 · Identifiability, sloppiness, third mRNA — 3 weeks

**Do (a).** Fisher Information Matrix from trajectory sensitivities alone — `JᵀJ` over the observables, eigendecompose, plot the spectrum. **No data required.**

**Do (b).** Rank u₃ candidates — E47/TCF3, NR5A2, MIST1/BHLHA15 — by boundary shift per unit dose **across the surviving ensemble**, not at one parameter set.

**Proves.** (a) You know what the model can and cannot know. (b) Payload composition is derived, not asserted.

**Buys.** (a) pre-empts the hardest question. *"Thirty unmeasured parameters, three observations — why isn't this curve-fitting?"* becomes: *"the model is sloppy, here is the spectrum, here are the stiff directions carrying every prediction, and here is why the topology ranking is robust to the sloppy ones."* Gutenkunst et al. (PMID 17922568) — whose own prescription is **"focus on predictions rather than parameters."**

**Do not compute AIC/BIC.** With n ≈ 10 and k ≈ 30 it is illegitimate, and someone will ask.

## STAGE 3B · In-silico perturbation screen — 2–3 weeks · ★ THE UNBIASED ARM

**Why this stage exists.** The mechanistic model's states *are* the payload species, so it cannot derive payload identity without circularity (§1.3). This arm supplies identity from data the modeller did not curate.

**Do.** **CellOracle** (Kamimoto, Stringa, Hoffmann, Jindal, Solnica-Krezel & Morris, *Nature* **614**(7949):742–751, 2023, PMID 36755098) on pancreatic scRNA-seq containing **both** acinar and metaplastic populations. Infer the GRN, simulate overexpression across the **full TF repertoire present in the network** — not a three-candidate shortlist — score the shift of the ADM/metaplastic cluster toward the acinar cluster, and produce a **ranked list**.

**Dataset.** Start with **GSE207938** — smallest, cleanest, <1 GB, and already required for Stage 7. **Download once, use twice.** Repeat on **GSE172380** as the stability check (below).

**Precedent for this exact task shape.** Kamimoto, Adil, Jindal, Hoffmann, Kong, Yang & Morris, *Stem Cell Reports* **18**(1):97–112, 2022, PMID 36584685 — in-silico perturbation identified the factors that rescue a **failing** fibroblast→iEP conversion, nominating **Fos** with **Yap1**. Note what makes it a good precedent: the nominated factors were *not* the field's prior expectation. A screen that can only confirm is not a screen.

### Why this is not the Perturb-seq module that was killed

Reasonable objection, and it must be answered in the writeup before a judge raises it. The killed approach was **vector arithmetic on Replogle Perturb-seq**, and the reason it was killed was **cell context**: RBPJL is "Not detected" in every Human Protein Atlas cell line, the K562 library targeted only K562-expressed genes, and the RPE1 library is common-essential only — so the genes of interest **structurally could not be in the library**. That is a property of the *data source*, not of in-silico perturbation as an idea.

CellOracle on a pancreatic dataset changes the data source: acinar and metaplastic cells are present, and the network therefore contains the relevant factors. **The method is different too** — GRN inference plus signal propagation, rather than treating measured perturbation profiles as displacement vectors to be added.

**This is not a free pass.** CellOracle has its own limits, stated here so they are not discovered late:

- **It is local and first-order.** CellOracle propagates a perturbation through a linearised GRN and projects the result onto the existing embedding. It predicts the *direction of movement in the neighbourhood of observed cells*; it does not and cannot establish that a cell crosses a separatrix into another basin. **That is precisely why the two arms are complementary rather than redundant: the screen nominates, the dynamical model tests whether the nomination can actually complete the transition.** State it this way round — it converts a limitation into the division of labour.
- **Species.** The screen is mouse; the wet lab is rat. This is the same bridge every validation dataset in this project already crosses (Collins and Krah are mouse), so it introduces no new inconsistency — but it is not zero, and it should be said aloud rather than left for someone to notice.

### ⚠ Pre-flight — check these before trusting any output

1. **Is RBPJL perturbable?** CellOracle's edges come from motif scanning, so a TF with no motif in the database has **no outgoing edges** and its simulated perturbation returns ≈ 0 — *silently*, and indistinguishably from "RBPJL doesn't matter." RBPJL is obscure enough that this is a live risk, and it is the one component the entire mechanistic argument rests on. **Confirm RBPJL is present as a regulator with outgoing edges before running anything.** If it is not, report that limitation explicitly and do not let a near-zero RBPJL score be read as a biological result.
2. **Is the base GRN actually pancreatic?** If the dataset has no matched scATAC, CellOracle falls back to a generic base GRN. Check that the resulting network contains PTF1A→acinar-target edges. If it does not, the screen is measuring the motif database, not the data.
3. **Is the ranking stable?** Repeat on GSE172380. **Top-20 overlap below ~50% means the ranking is inference noise and no result may be reported from it.**

### Pre-registration — non-negotiable for this stage

Before the screen runs, commit **and push**:
- `prereg/<date>_celloracle_screen_ranges.yaml` — dataset, preprocessing, GRN settings, TF set, scoring metric.
- `prereg/<date>_celloracle_screen_prediction.md` — what the screen is predicted to nominate, and **explicitly where PTF1A, RBPJL and NR5A2 are expected to rank.**

**Whatever it returns gets reported, including if it ranks those three low.** A screen whose unwelcome outcomes would not have been published is not evidence, and the pushed timestamp is what makes that commitment checkable. The prediction may legitimately be informed by Stage 3's ranking — say so in the file.

**Proves.** That payload identity was derived by a method capable of returning a different answer.

**Buys.** The answer to *"your model only nominates what you put in it."* Convergence between arms is worth more than either alone, since the two have unrelated failure modes; divergence is more interesting still and arrives early enough to act on.

## STAGE 4 · Held-out prediction — 2 weeks · ★ THE HEADLINE RESULT

**Constrain on:** WT caerulein recovery (~7 days absent oncogenic KRAS) + the PTF1A dose ladder (§6, HO-3).

**Then predict, nothing else tuned:**
1. **The KRAS-history effect** — long-exposure cells revert markedly slower than short-exposure under *identical* trametinib, from a single history parameter. (HO-1)
2. **The viability dissociation** — trametinib and KRAS-extinction at opposite corners of the reversal × viability plane. (HO-2)

**Pre-register both** — commit and push before running the comparison.

**Proves.** The model predicts things it was not built on.

**Buys.** Reproducing what you fit is a consistency check. Predicting what you held out is validation. **And report one thing the model gets wrong** — *"predicted the history effect to within a day but overpredicted recovery in the apoptosis arm, which tells me my survival branch is too simple"* demonstrates understanding rather than operation.

## STAGE 5 · The dosing schedule — 3 weeks

**Do.** Sweep **(dose per mRNA pulse) × (redosing interval)** — the real actuatable axes. Pure simulation, no optimizer. Classify every cell:

- **FAIL–undershoot** — never leaves the ADM basin
- **SUCCESS** — reaches the acinar basin with stress within tolerance throughout
- **FAIL–toxic** — would reach the target, but the stress axis crosses first

**State it correctly.** This is a 2D slice, in input space, through the **capture basin of the acinar target viable in the safe set** — Aubin's viability theory (Aubin 1991; *SIAM J Control Optim* 40(3):853–881, 2002), equivalently the **backward reach-avoid set** (Margellos & Lygeros, IEEE TAC 56(8):1849, 2011).

> **Get this right or be corrected in public:** the **viability kernel** is the set from which you can remain safe *forever*. That is not what you want — you want to *reach* the acinar basin. The object is the **capture basin**.

Write verbatim: *"Exact Hamilton–Jacobi computation is intractable at this dimension, since grid-based solvers scale as O(Mⁿ); we therefore compute a sampled inner approximation, which is conservative by construction."*

**Proves.** The efficacy and toxicity boundaries have **different slopes**.

**Buys.** A falsifiable structural claim at zero extra compute. If dose × duration were all that mattered, both boundaries would be hyperbolas `D·τ = const`. Deviation from that is the signature of the slow chromatin integrator (a leaky accumulator weights *sustained* over *peak*) and of ID3 titration saturating at high dose.

Also report **minimum total mRNA along the efficacy boundary**, and: *"trametinib shifts the efficacy boundary left by X-fold while leaving the toxicity boundary unmoved, widening the window by Y."*

## STAGE 6 · Intervention ordering via Lie bracket — 1 week · ★ THE NOVEL BIT

**The question.** Four inputs — trametinib, RBPJL, PTF1A, third mRNA. Does order matter?

**Framing warning.** "Order matters" is trivially true for nonlinear systems; flows do not commute. The non-trivial questions are the **sign**, the **magnitude**, and whether the sign is **invariant across the ensemble**.

**The formalization.**
```
Δ(s) = Φ_B^s ∘ Φ_A^s (x₀) − Φ_A^s ∘ Φ_B^s (x₀) = s²·[f_A, f_B](x₀) + O(s³)
```
Order matters to leading order **iff the Lie bracket of the closed-loop vector fields is nonzero at the ADM fixed point**, growing as s². Since exogenous mRNA enters additively, the supply vector fields are constant and their mutual bracket vanishes — **the entire ordering effect comes from the nonlinearity of f alone**, governed by the failure of the `P_n` and `R` Jacobian columns to commute.

**Five lines of CasADi or SymPy.** Build f symbolically, `jacobian(f, x)`, extract two columns, form the bracket, evaluate at the ADM fixed point, report `‖[f_A,f_B]‖ / ‖f‖`. A closed-form answer **before** a single expensive simulation.

**Mechanistic prediction, stated up front.** RBPJL cannot bootstrap — endogenous R only rises after C_L rises — but exogenous u₁ bypasses that. **RBPJL-first pre-loads the pool**, so arriving PTF1A converts to complex as fast as it arrives. **PTF1A-first wastes protein** into the ID3 and degradation sinks while waiting for the slow R loop. Trametinib-first lowers ID3 before either arrives, reducing the sequestration tax on both.

> **Predicted optimal sequence: trametinib → RBPJL → PTF1A**, with the effect collapsing if RBPJL turnover is fast. That conditional makes it a prediction rather than an observation.

**Simulation design — non-negotiable:** hold total delivered dose constant. Fix ∫u₁dt and ∫u₂dt; vary only the phase offset Δ. Otherwise "RBPJL-first is better" reduces to "more protein delivered." Report **ordering gain** M_min(0)/M_min(Δ\*), optimal offset in hours, and ensemble fraction preserving sign(Δ\*). Gain of 1.05× → say so and drop the claim. 2–5× → real result.

**Then close the loop:** show `‖[f_A,f_B](x₀;θ)‖` **correlates with measured ordering gain across the ensemble.**

**Precedent.** Lee, Ye, Gardino, Heijink, Sorger, MacBeath & Yaffe, *Cell* 149:780 (2012), PMID 22579283 — **time-staggered** EGFR inhibition then genotoxic chemo dramatically sensitizes TNBC; simultaneous does not. Also Behar et al., *Cell* 155:448 (2013), PMID 24120141.

**No published precedent found for Lie-bracket formalization of cell-fate intervention ordering.** If that holds it is claimable, and it costs an afternoon.

## STAGE 7 · Single-cell falsification — 2 weeks

**Do.** Gillespie / chemical-Langevin on the reduced network with log-normal per-cell LNP dose spread. Predict the **distribution**, not the mean.

| Accession | What | Why |
|---|---|---|
| **GSE314765** | >300,000 cells, **7 KRAS timepoints (14 h → 12 wk)** | Only densely time-resolved epithelium-enriched ADM atlas; RAP/CP arms give a built-in replication cohort |
| **GSE207938** | FACS mKate2⁺ epithelium, N1→N2→K1–K6, **3 .h5ad, <1 GB** | Best effort-to-value ratio on a laptop (Burdziak, *Science* 2023, PMID 37167403) |
| **GSE172380** | >13,000 **EYFP lineage-traced acinar** cells | Cleanest ADM ground truth — you know they started acinar (Ma 2021, PMID 34695382) |
| **GSE141017** | 6 timepoints, **per-barcode metaplastic identities** | Six published metaplastic states = a concrete discrete hypothesis (Schlesinger 2020, PMID 32908137) |

**Proves.** Bistability at the level of the observable that can falsify it — **bimodality**, not population means.

**Buys.** The PU.1/GATA1 lesson applied correctly. That toggle was a widely accepted bistable TF model built on **population-averaged** data. Single-cell imaging killed it twice — Hoppe 2016 (PMID 27411635) found dynamics *"incompatible"* with stochastic switching preceding lineage choice; Strasser 2018 (PMID 30002371) titled their paper *"…refutes the PU.1/GATA1 toggle switch paradigm."* The failure was **validating on the wrong observable.** Raise it yourself and you look like a scientist.

**The formal discreteness test is the contribution, not the trajectory.** Decipher (Nazaret 2025, PMID 40702544) already models trajectories on GSE207938 — pseudotime alone is not novel.

## OPTIONAL · Sequence module — 2 weeks, parallel, cut first

**Not Enformer, not Borzoi.** No rat head in Enformer, Borzoi, or AlphaGenome; Enformer scores **r = 0.137 on enhancers** vs 0.81 on promoters (Karollus 2023, PMID 36973806); and in-silico mutagenesis is single-nucleotide **additive**, so it structurally cannot represent cooperativity — which is what a Hill coefficient is.

**Instead:** a **bipartite PTF1 scanner** on the mouse locus — `E-box(CANNTG) + N₅ + TC-box`, spacer length an explicit scannable parameter. JASPAR's six Ptf1a matrices encode only the E-box half; **none encodes the TC-box or spacer**, so naive FIMO returns noise. ~50–100 lines of Biopython.

**Validate three ways:** (a) recovers both known sites in the 2.3-kb enhancer; (b) composite hits enriched in **GSE86262** PTF1A ChIP-seq peaks vs shuffled; (c) scanning N₃/N₄/N₆/N₇ shows a sharp optimum at **N₅** — reproducing Beres 2006 (PMID 16354684) from data.

**Rat bridge:** rn7 has **no conservation tracks at all**. Use **mm39 `cons35way`** (rat is an aligned species), then liftOver to rn7 via the verified `chainMm39`/`netMm39` tracks.

**Added task — settle the NR5A2 placement assumption qualitatively.** T6a assumes NR5A2 acts at the *Ptf1a* autoregulatory enhancer, which Holmstrom 2011 does not show. **GSE34295 can partially settle it:** *"is there LRH-1 signal at the Ptf1a locus, yes or no"* is answerable from that data. **Qualitative only** — it is 2 samples, single-replicate, 2011-era ChIP-seq with no biological replication, so it may support presence/absence and nothing more. A yes/no here does not set a parameter; it tells T6a whether its central assumption has any footing. Do not report a number from it.

---

# PART 5 · TIMELINE

**Spine estimate:** 240–320 focused hours for the ODE core — model coding 30–40 h, numerical debugging 40–60 h (stiffness ratio ~10⁴–10⁵ across minutes→weeks; needs BDF/Radau and careful scaling), continuation 40 h *if reduced first and hand-rolled*, bistability detection 25–40 h, sweep 40–60 h, analysis 60–80 h. **6–8 weeks full-time.**

| Stage | Weeks | Cut order |
|---|---|---|
| 0 · Reduce and certify | 2 | **never** |
| 1 · Two-parameter bifurcation ★ | 2 | never |
| 2 · Topology competition ★ | 4 | never |
| 3 · Identifiability + third mRNA | 3 | 5th |
| 3B · In-silico perturbation screen ★ | 2–3 | **2nd** |
| 4 · Held-out prediction ★ | 2 | never |
| 5 · Dosing schedule | 3 | degrade, don't drop |
| 6 · Ordering / Lie bracket ★ | 1 | 4th |
| 7 · Single-cell falsification | 2 | 3rd |
| Sequence module (parallel) | 2 | **1st** |
| Figures, writeup, buffer | 3 | — |
| **Total** | **~24–25** | |

**Stage 3B placement.** After Stage 3, deliberately. Running it first would make Stage 3's third-mRNA ranking a downstream confirmation of the screen; running it second keeps the two arms **independent**, so agreement between them is evidence rather than construction. It shares datasets with Stage 7 — **download once, use twice**. In the cut order it goes **2nd**, after the sequence module and before Stage 7: it is more expendable than the falsification arm, and it is never cut before Stage 4.

**Because the wet lab runs in the same cycle**, the stages the bench actually needs are **1** (which KRAS × trametinib doses to test), **5** (mRNA dose and interval), and **6** (the order). Stages 3 and 7 strengthen the writeup but do not change what gets pipetted. If the lab date compresses the schedule, run **0 → 1 → 2 → 4 → 5 → 6** and hand over a protocol.

**Gates — decide now, not later.**
- **End of Stage 2** — does *any* topology reproduce the 3-day/3-week asymmetry? If none, the loop architecture as described cannot account for the kinetics and something is missing. A genuine result, arriving in month 2.
- **End of Stage 4** — does the KRAS-history effect fall out of one parameter? If not, report the mismatch and name which observation breaks it.

**Neither failure ends the project.** What ends it is discovering either one in month 7.

**Order note:** run Stage 5's sweep before any optimizer. Pure simulation, presentable figure in week 1, and it immediately reveals whether the window is non-empty. Never spend three weeks on an optimizer for a problem with no feasible point.

---

# PART 6 · VALIDATION TARGETS, WITH NUMBERS

Collins 2014 is ***Gastroenterology* 146(3):822–834, PMID 24315826** — not JCI, not Cancer Res. Dosing: PD325901 10 mg/kg q12h, begun 5 weeks after PanIN establishment.

### HO-1 · The history effect — the best held-out test

Same drug, same dose, same readout, **different KRAS history**:

- **iKras\***: 1.5 d — PanINs prevalent, pERK already dramatically reduced. 2 d — *"amylase-positive acinar cells represented approximately 30% of the epithelial cells."* 3 d — *"increased to almost 50%."* 5 d — *"very few PanIN lesions remained."* 7 d — only acinar + ADM.
- **KC**: *"acinar cells did not arise until after 5 days of treatment, and PanIN lesions still persisted at this time point"* — and decisively, *"the delay in tissue repair of KC pancreata was not caused by insufficient inhibition of MAPK signaling because pERK1/2 levels were reduced at all time points."*

**With inducible KRAS in AR42J this is directly reproducible in vitro.**

### HO-2 · The viability dissociation — the headline claim, already published

- **MEK inhibition** (Collins 2014): cCasp3 *"rare positive staining at all time points"*; *"the death of PanIN cells was unlikely to explain their elimination from the tissue over time."* → **reversal, low death.**
- **KRAS extinction** (Collins 2012, PMID 22232209): *"KrasG12D was required for tumor cell survival"*; removal at this stage *"resulted in extensive cell death."* → **reversal, high death.**

If the model cannot reproduce this dissociation, the headline claim is dead on arrival. **Figure 1, not an afterthought.**

### HO-3 · The PTF1A dose ladder (Krah 2015, PMID 26151762)

- *Ptf1a* cKO + Kras vs Kras alone: **>15-fold increase** in Alcian-blue⁺ PanIN burden per cm².
- *Ptf1a* **heterozygous** on KC at 1 month: *"increased PanINs at this early stage"* when KC PanIN formation is normally minimal.
- *Ptf1a*-het **KPC**: PDAC *"much earlier,"* log-rank **p < 0.01**; liver mets **3/9 vs 0/9**.
- **Threshold signature:** *"we were surprised that a moderate level of acinar cell recombination (~25%) failed to produce an overt, short-term phenotype"* — while a **4.5× higher** tamoxifen dose produced rapid ADM. A dose that does nothing and a higher dose that does a lot is exactly the nonlinearity a bistable model predicts.
- **Constraint on over-claiming:** even extensive deletion gives *"rapid but incomplete acinar-ductal metaplasia."* Your model must not predict 100%.

### HO-4 · Reversal is partial and hybrid

Krah 2019: *"nearly all PanIN-associated acinar cells in [tetO-Ptf1a] mice maintained strong Sox9 expression"*; *"a subset of morphologically normal acinar cells did remain pERK+."* Schlesinger 2020 resolves **six** metaplastic states. Del Poggetto 2021 (PMID 34529467): after *complete* resolution of transient inflammation, cells retain *"an enduring adaptive response associated with sustained transcriptional and epigenetic reprogramming."* Lo 2025 (PMID 40156071): differential methylation **persists after expression returns to normal.**

**The model must permit reversal on the transcriptional axis while SOX9 and pERK remain elevated** — and it does, because `M` relaxes far slower than the TFs. That also predicts **relapse risk stays elevated long after phenotypic reversal looks complete.**

### HO-5 · Baseline recovery kinetics

- Collins 2014 controls: 2 d post-caerulein, *"tissue-wide acinar-ductal metaplasia"*; 1 week, *"the tissue had recovered."* Spontaneous full reversal absent oncogenic KRAS ≈ 7 days.
- Krah 2019 controls: *"fully recovered at 3 weeks following induction of pancreatitis."*

**Architecture:** constrain on HO-5 + HO-3. **Predict HO-1 and HO-2.** Use HO-4 as the qualitative structural test.

---

# PART 7 · PRE-FLIGHT — verify before the wet lab

1. **⚠ Sequence KRAS in your AR42J stock. Highest-stakes unknown in the project.** AR42J derives from an azaserine-induced rat pancreatic tumour, and azaserine-induced rat acinar tumours are classically Kras-mutant. **If AR42J already carries an activating Kras mutation, inducing G12D on top is close to meaningless** — the "KRAS-induced ADM" framing collapses and the design pivots to **KRAS dose titration on a mutant background** (still clean and defensible, but a different experiment with a different forcing term). *Unverified — PubMed queries returned zero. Resolve before anything else.*
2. **Confirm parental AR42J (ATCC CRL-1492), not B13.** Pre-register transferrin / CK7 counterscreens against hepatic transdifferentiation.
3. **Does ADM hold ≥14–21 days in your protocol?** Krah's redifferentiation took 3 weeks.
4. **ID3 western + PTF1A localization IF.** Determines the third mRNA empirically. PTF1A cytoplasmic + ID3 high → E-protein limiting, u₃ = E47. PTF1A absent, ID3 normal → transcription is off rather than protein trapped.
5. **Trametinib dose–response in AR42J** — local IC50 for pERK, not a literature value from another line.
6. **LNP transfection efficiency and mRNA half-life in your cells, in Matrigel.** Stage 5's x-axis has no units without this.

**The experiment the model designs is a dose titration with up- and down-sweeps, not a single dose held for a duration.** Hysteresis requires both directions. Readout must be **single-cell and bimodality-capable**.

---

# PART 8 · KNOWN GAPS — state aloud, do not paper over

Explicit **"no evidence found"** after targeted search:

- Experimental hysteresis or bistability of acinar/ductal identity — never tested by anyone.
- Any dose–response for the PTF1A autoregulatory loop.
- PTF1A protein or mRNA half-life.
- Any Hill coefficient for PTF1 enhancer occupancy.
- Any measured Kd, k_on, or k_off for ID3·E47, ID3·PTF1A, or PTF1A·RBPJ/RBPJL. Only relative yeast-two-hybrid (Langlands 1997, PMID 9242638) and co-IP.
- Any quantification of a secretory-capacity ceiling in acinar cells.
- **Any prior ODE or dynamical-systems model of ADM.** The field is open.
- AR42J Kras genotype (see §7.1).
- No public AR42J transcriptomics with TGF-α or EGF; nothing annotated as ADM. All eight GEO series containing "AR42J" were enumerated.
- Human ADM single-cell data does not exist publicly.

## Ruled out — do not re-propose

| Cut | Reason |
|---|---|
| AlphaFold/Boltz → `k_on` | AF outputs regress model accuracy, never affinity. And `k_on` is structurally non-identifiable here — binding is quasi-steady-state on transcriptional timescales |
| Enformer/Borzoi → Hill coefficient | ISM is additive, cannot represent cooperativity; r = 0.137 on enhancers; no rat head |
| ipTM linker-length scan | ipTM normalizes over chain length — the curve is arithmetic, not biophysics |
| Perturb-seq cocktail derivation | RBPJL is "Not detected" in every HPA cell line; the libraries structurally cannot contain the genes. **Still ruled out** — but see the note below, this does *not* rule out in-silico perturbation on a pancreatic dataset |
| Fast acetylation as a state | Turnover < 30 min → adiabatically slaved → the variable algebraically vanishes |
| One-sided `U_crit` | Capacity is co-induced with cargo; no ceiling exists at steady state |
| Structural / Kalman controllability | Linear theorem, generic answer, and *local* while reversal is *global* |
| Full-dimension HJ reachability, MPC, all-atom MD | O(Mⁿ) intractable; no sensor for feedback; needs µs sampling |

> **Perturb-seq vs Stage 3B — do not confuse these.** Perturb-seq cocktail derivation was killed on **cell context**: the K562/RPE1 libraries could not contain the genes, because those lines do not express them. That is a fact about the *data source*, not a verdict on in-silico perturbation as a method. Stage 3B changes the data source to a pancreatic dataset where acinar and metaplastic cells — and therefore the relevant factors — are present, and changes the method to GRN inference plus signal propagation rather than treating measured perturbation profiles as displacement vectors. The original objection does not transfer. Stage 3B's own distinct limitations are stated in full in that section, including the RBPJL-motif risk, which is the analogous failure mode and must be checked before the screen is trusted.

---

## THE ONE-LINE SUMMARY

**Five architectures are consistent with the biology of KRAS-driven ADM. One reproduces the 3-day-versus-3-week kinetic asymmetry between MEK inhibition and forced PTF1A. That architecture — whose parameters are demonstrably sloppy — nonetheless predicts the KRAS-history effect and the viability dissociation without being fit to either, and it outputs a trametinib-plus-three-mRNA schedule: which dose, at what interval, in what order.**
