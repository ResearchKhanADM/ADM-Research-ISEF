# ADM durability — in-silico project

**Goal: make ADM reversal *stick*.** Reversal itself was solved in 2014. The wet
lab runs in the same cycle and gets **one shot, no iteration** — this model exists
to design it, and must produce answers in time to inform it.

Plan of record: `docs/ADM_MASTER_PLAN_v3.md` (**Part 0 first, not last** — it
explains why the project changed shape). Nothing else is the plan. Running log:
`docs/PROGRESS.md`, read the last entry every session.
Repo: https://github.com/ResearchKhanADM/ADM-Research-ISEF (public)

## The claim

Collins 2014 (PMID 24315826) showed MEK inhibition reverts established PanIN in
~3 days — and that ADM/PanIN **resumed** when the drug was stopped with KRAS left
on. "Reverse ADM" is eleven years late. The unoccupied claim is durability:
**trametinib reverts, it does not stick.** *Rbpjl* has no PTF1A-independent
promoter, so the loop cannot re-close on its own — that bootstrap failure is why
reversion is drug-dependent, and it is what the payload is for. Payload
**identity is not the claim**; it is published (Jiang 2023, PMID 37425649: PTF1A,
RBPJL, BHLHA15/MIST1). The claim is **composition, ratio, formulation and schedule
under a fixed delivery budget** — three results everything serves: **R1
formulation ★** (co-formulated vs separate particles across the plausible uptake-CV
range; for an obligate stoichiometric pair a requirement, not a preference) · **R2
composition** (how many components at what ratio under fixed total mRNA mass — a
marginal-value curve, not a subset table) · **R3 durability** (the drug-free
persistence window — **what the bifurcation diagram is for**, the answer not
decoration).

## Layout

**Committed:** `src/core.py` the 3-state model · `src/supplementary/` kept, not part
of the core · `scripts/` entry points · `prereg/` ranges + predictions written
*before* runs · `figures/` `_style.py`, `_provenance.py`, one module per figure,
`out/` holding final PDFs + `_source.csv` + `.prov.json` · `docs/decisions/` one file
per architectural decision · `tests/golden/` pre-rewrite fixture (**arbitrary
parameters — report no number from it**). **Gitignored:** `data/` `results/`
`logs/` `build/`.

## Phases — Part 2. Phase 1 runs parallel to Phase 2.

| # | Phase | Wk | Gate / lock |
|---|---|---|---|
| 0 | **Decision spec** — cell system, assays, viability floor X, durability timepoint, reversal score. **Issue Bench Handshake #1 — all 8 items blocking.** Arm budget is an unknown **external** input: placeholder 12/24/48 wells, design instantiable at any tier, and say what changes at each. | 1 | Arm budget |
| 1 | **Candidate generation** — regulon activity across **≥2 independent datasets**, per-dataset frequency. Add an **ADM-repressor axis**. Cross-ref Joung TF Atlas (PMID 36608654) + human TFome (PMID 33257861). | 2 | **GATE A** → Layer 1 lock, **order PTF1A + RBPJL mRNA** |
| 2 | **Minimal mechanistic core** — 3–4 states, ~9–12 parameters, **nondimensionalized**. | 4 | **GATE B** |
| 3 | **Delivery layer** ★ **the headline** — per-cell LNP dose as a distribution; correlated vs independent double-above-threshold fraction. | 3 | Formulation + ratio |
| 4 | **Mixture-amount design** — simplex × total amount, Scheffé polynomial. | 2 | Layer 2 lock: count + mass split |
| 5 | **Durability + schedule** — time-to-relapse after clearance **with trametinib withdrawn**; dose × redosing map; the two-block ordering question. **Pre-flight: bench item 8 (ID3 ± trametinib) before running the ordering arms.** | 2 | Redosing interval + order |
| 6 | **Held-out prediction** — KRAS-history effect; viability dissociation. | 2 | **GATE C** |
| 7 | **Pre-registration + discrimination power** — design to *separate hypotheses*, not to confirm. | 2 | Locked protocol |
| 8 | **Figures, writeup, buffer** | 3 | — |

**Gates — stop and report. A failed gate is a real result, not a reason to push on.**
**A** (wk 2) does the pipeline recover PTF1A / RBPJL / BHLHA15? — if no, the pipeline
is broken in month one. **B** (wk 6) two stable states plus a saddle, identifiable
separatrix, key parameters surviving profile likelihood. **C** does the KRAS-history
effect fall out of a **single** parameter? **Never cut:** 0, 2, 3, 7. **Cut order
under compression:** 6 → Phase 4's simplex interior → Phase 5's ordering arms.

## The Phase 2 model — 3–4 states, not thirteen

- `P` PTF1A, autoregulatory via PTF1-L (**RBPJL**), needs an E-protein partner,
  plus ERK-suppressed basal ignition — **without which `P = 0` is absorbing too**
- `R` RBPJL, from PTF1-J (**RBPJ**, constant pool). **No P-independent term — that
  zero *is* the bootstrap claim.** *Nothing but PTF1A makes RBPJL* — **not** *nothing
  but RBPJL makes RBPJL*: that is stronger, unsupported, makes `R = 0` absorbing, and
  contradicts Collins. It also passes every guard test.
- `C` slow chromatin/memory at metaplasia loci — sets time-to-relapse
- `E_free` from the **exact** binding solution, `ID3 = f(pERK)` — algebraic. The
  linear `E_total − k·ID3` is its tight limit, negative outside its domain; **no
  floor hack.** Sharpness **`n_eff ≈ 0.5/√κ`, measured on the shipped form** — the
  derivation's 1.34 is the *deleted* two-target complex; re-citing it overstates
  sharpness ~2.7× and inflates R1.
- **pERK is an input, not a state — and on withdrawal it is not a step.** It follows
  a prescribed **rebound profile**, swept over its shape, measured by bench item 7.
  Withdrawal *is* the endpoint, so the recovery shape is mechanism (002 amendment).
- **Viability is measured at the bench, not modelled** (008 retired). One survivor,
  a flag not a term: warn when predicted `P` drops below the CHOP threshold.
- Profile likelihood on `a_P`→R2, `γ`→R3, `κ`→R1 — **it is the uncertainty on each
  deliverable**, not an identifiability side-quest. **Not** an FIM eigenspectrum.
- **Flagship guard test:** *`dR/dt` has no P-independent term* — as a
  guard-on-the-guard (construct the violating version, require it to fail), plus a
  dynamic companion: `R` must not rise from zero while `P` is held at zero.

## Hard rules

**1. Pre-registration before any sweep or fit.** Ranges to
`prereg/<date>_<name>_ranges.yaml`, prediction to `prereg/<date>_<name>_prediction.md`,
**commit and push**, then run — the pushed timestamp is the proof. **I write the
templates; Luqmaan writes the content.** Asked to skip this: refuse and say why.

**2. Nothing large or secret in git.** Fetch datasets by accession with an idempotent
script — the script is the artifact. Public repo: no tokens, no `.env`.

**3. Long jobs run detached, never as a Claude Code background process.** Bash caps
near 10 min; session-owned processes die on exit. Checkpoint to `results/<run_id>/`,
`--resume`, log to `logs/`, launch detached — designed up front.

**4. Architectural decisions get a 2-subagent adversarial panel**, then
`docs/decisions/NNN-name.md`: question, both positions, decision, **what would
reverse it**. Routine coding choices don't get a panel — just decide.

**5. The student must be able to explain every line.** Comment *why*, not *what*.
Stop and explain unfamiliar techniques (continuation, profile likelihood, Scheffé
polynomials, compound-Poisson dose, Bliss independence) in chat. Clear over clever.
**Push back when he's wrong.**

**6. Append to `docs/PROGRESS.md` every session** — done / broke / next / open
questions. Push. **7. A retired claim is retired, not paused** — 002 and 008 were
retired, not shelved; an unretired claim is one someone builds on later.

## ⚠ STANDING RULE — never silently drop a failed solve

**Silently dropped failures are the easiest way to produce a confident wrong
answer**, because failures *correlate with swept parameters* — drop them and the
sample set is depleted precisely in the regime the sweep was built to probe, which
looks like a clean negative. Everywhere a solve happens (continuation, profile
likelihood, Phase 3's Monte Carlo): **log every outcome**; **report failure rate as
a function of every swept parameter, never as a scalar**; **if it correlates, say
so** and treat results as *conditional on convergence*; **keep a test sweeping the
hard regime**, bound <1%. The core's binding step is closed-form and exempt.

## ⚠ STANDING RULE — figure modules never compute science

**A figure module loads `results/`, does display arithmetic, and draws. It never
runs an ODE solve or a sample.** If a figure needs a number not in `results/`, the
fix is a stage that writes it. Consequence: rebuild is always under a minute,
figures cannot drift from the analysis, and a slow stage never blocks a figure.
**Produce figures as you go, every session — not in Phase 8:** `python make_figures.py`.

- **General-purpose geometry. No IEEE or conference-template sizing** — one `SCALE`
  in `figures/_style.py`, scaled downstream by Luqmaan.
- Okabe–Ito, **four nominal colours**: `#0072B2` acinar · `#D55E00` metaplastic ·
  `#009E73` intervention · `#6E6E6E` toxic — **grey, not red** (deuteranopia).
- `rc_context`, never a bare rcParams mutation — one process imports every module.
  `savefig.bbox="standard"` **not** `"tight"` (tight silently resizes, breaking width
  checks). `pdf.fonttype 42`, `svg.fonttype "none"`.
- **Global grammar: a solid line is a model; an open marker with a dark edge is data.
  Never connect experimental points with a line.**
- Every stage calls `stamp_run()`; figures declare inputs and `save_figure()` hashes
  them. `paper`/`poster` **refuse a dirty tree**; `--check` fails on staleness.
- **Traps:** no 3-D cusp surface · no hairball of all states · no continuous colormap
  under a three-class map · no fake error bars on deterministic output · no p-value
  without sampling · no truncated percentage axis · no bimodality from a KDE shape ·
  two significant figures unless defensible.

## Cut — do not re-propose. Reasons in Part 4 and decision 012.

**Five-way topology competition** (10× discriminator unmeasured — Krah 2019 has no
timepoint between 24 h and 3 wk; bound is (1×, 21×)) · **13-state ODE** · **FIM
sloppiness** · **Lie-bracket ordering formalism** (precedent: Letsou & Cai 2016,
PMID 27560383 — the novelty claim was false) · **2^k subset enumeration** (mass
budget makes this a mixture problem; enumeration visits only vertices) ·
**submodular/greedy optimization** (assumes diminishing returns; the bootstrap
threshold claims *super*modularity) · **Pareto front as deliverable** (explicit
reversal — a front is preference-free and a one-shot experiment consumes a
decision; replaced by constrained optimization + therapeutic index) · **U-shaped
viability hazard as an ODE** · **trametinib-vs-PD325901 asymmetry** · **Gillespie
bimodality vs scRNA-seq** · **CellOracle as validation** (keep as a *declared
negative control*) · **Enformer/Borzoi · AlphaFold→k_on · Perturb-seq · one-sided
U_crit · structural/Kalman controllability · full-dim HJ reachability · MPC ·
all-atom MD**. To revisit one, argue it in `docs/decisions/`. Don't just start.

## Standing limitations — state these before anyone asks

- **A regulon screen is blind to post-translational mechanisms.** TCF3/E47 is not
  transcriptionally lost, it is titrated — so no threshold surfaces it; it enters
  by **declared mechanism, registered before the screen runs.**
- **ID3 titration is DOCUMENTED in this cell line. The gap is one edge: ERK→ID3.**
  Dufresne 2010 (PMID 20830706), *in AR4-2J*: Id3 binds both E47 and Ptf1-p48, those
  interactions rise while E47/Ptf1-p48 falls, silencing Id3 reverses the
  mislocalization, pattern holds in human/murine lesions. Its driver is **gastrin**
  — that KRAS/ERK drives ID3 is the inference, testable in one western (bench item
  8). **No Kd** measured. Don't disclaim the node; v3.1 Part 1.4 is the wording.
- **PTF1A is pleiotropic** and dosage-sensitive (PMID 30470852). "Preserving
  viability" does not cover "did not make a neural-program-expressing cell" —
  **lineage fidelity is a separate endpoint.**
- **Converted cells often fail to silence the starting program** (CellNet, PMID
  25126793) — hence the ADM-repressor axis. **The prioritization is a positive
  control, not a discovery.**
- **AR42J needs dexamethasone** to express amylase at all; without baseline and
  dynamic range, falsification power is unknown. **Its *Kras* genotype is
  unverified.** Both week-1 blocking. **Reference data is mouse/human; the bench is
  rat, and transformed** — argue the transfer, don't assume it. **No measured PTF1A
  half-life, no Hill coefficient, no Kd for any ID3 interaction, no prior ODE model
  of ADM** — sampled with stated priors, never fitted points.

## Framing — the rest is v3 Part 6

Category **CBIO**. **Frame the wet lab as computation-that-designed-an-experiment,
never as validation** — claim validation, get a negative, conclusion deleted.
**Effect size and pre-registered direction, never a p-value.** Bound what it can
falsify: ordering and timing yes, dose no.

## Environment

`venv\Scripts\activate`, `pip install -r requirements.txt`. Continuation is
hand-rolled pseudo-arclength on 3 states; no PyDSTool or AUTO-07p without a writeup.
