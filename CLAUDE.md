# ADM durability — in-silico project

**Goal: make ADM reversal *stick*.** Reversal itself was solved in 2014. The wet
lab runs in the same cycle and gets **one shot, no iteration** — this model exists
to design it, and must produce answers in time to inform it.

Plan of record: `docs/ADM_MASTER_PLAN_v3.md`. **Read Part 0 first, not last** — it
explains why the project changed shape. Nothing else is the plan.
Running log: `docs/PROGRESS.md`. Read the last entry at the start of every session.
Repo: https://github.com/ResearchKhanADM/ADM-Research-ISEF (public)

## The claim

Collins 2014 (PMID 24315826) showed MEK inhibition reverts established PanIN in
~3 days — and that ADM/PanIN **resumed** when the drug was stopped with KRAS left
on. "Reverse ADM" is eleven years late. The unoccupied claim is durability:
**trametinib reverts, it does not stick.** *Rbpjl* has no PTF1A-independent
promoter, so the loop cannot re-close on its own — that bootstrap failure is why
reversion is drug-dependent, and it is what the payload is for. Payload
**identity is not the claim**; it is already published (Jiang 2023, PMID 37425649:
PTF1A, RBPJL, BHLHA15/MIST1). The claim is **composition, ratio, formulation and
schedule under a fixed delivery budget.** Three results everything serves:

- **R1 · Formulation** ★ — co-formulated vs separate particles across the plausible
  uptake-CV range. For an obligate stoichiometric pair, a requirement not a preference.
- **R2 · Composition** — how many components, at what ratio, under fixed total
  mRNA mass. A marginal-value curve, not a subset table.
- **R3 · Durability** — the drug-free persistence window. **This is what the
  bifurcation diagram is for.** It is the answer, not decoration.

## Layout

| Path | Contents | Committed? |
|---|---|---|
| `src/` | model + analysis code | yes |
| `scripts/` | runnable entry points (sweeps, downloads) | yes |
| `prereg/` | parameter ranges + predictions, written *before* runs | yes |
| `figures/` | `_style.py`, `_provenance.py`, one module per figure; `out/` holds final PDFs + `_source.csv` + `.prov.json` | yes |
| `docs/decisions/` | one file per architectural decision | yes |
| `data/` `results/` `logs/` `build/` | datasets, sweep output, working renders | **no — gitignored** |

## Phases — Part 2. Phase 1 runs parallel to Phase 2.

| # | Phase | Wk | Gate / lock |
|---|---|---|---|
| 0 | **Decision spec** — cell system, **how many wells actually exist**, assays, viability floor X, durability timepoint, reversal score. Issue Bench Handshake #1. | 1 | Arm budget |
| 1 | **Candidate generation** — regulon activity across **≥2 independent datasets**, per-dataset frequency. Add an **ADM-repressor axis**. Cross-ref Joung TF Atlas (PMID 36608654) + human TFome (PMID 33257861). | 2 | **GATE A** → Layer 1 lock, **order PTF1A + RBPJL mRNA** |
| 2 | **Minimal mechanistic core** — 3–4 states, ~9–12 parameters, **nondimensionalized**. | 4 | **GATE B** |
| 3 | **Delivery layer** ★ **the headline** — per-cell LNP dose as a distribution; correlated vs independent double-above-threshold fraction. | 3 | Formulation + ratio |
| 4 | **Mixture-amount design** — simplex × total amount, Scheffé polynomial. | 2 | Layer 2 lock: count + mass split |
| 5 | **Durability + schedule** — time-to-relapse after clearance **with trametinib withdrawn**; dose × redosing map; the two-block ordering question. | 2 | Redosing interval + order |
| 6 | **Held-out prediction** — KRAS-history effect; viability dissociation. | 2 | **GATE C** |
| 7 | **Pre-registration + discrimination power** — design to *separate hypotheses*, not to confirm. | 2 | Locked protocol |
| 8 | **Figures, writeup, buffer** | 3 | — |

**Gates — stop and report. A failed gate is a real result, not a reason to push on.**
**A** (wk 2) does the pipeline recover PTF1A / RBPJL / BHLHA15? — cheap, fast,
decisive, and if no the pipeline is broken in month one. **B** (wk 6) two stable
states plus a saddle, identifiable separatrix, key parameters surviving profile
likelihood. **C** does the KRAS-history effect fall out of a **single** parameter?
**Never cut:** 0, 2, 3, 7. **Cut order under compression:** 6 → Phase 4's simplex
interior (vertices only) → Phase 5's ordering arms.

## The Phase 2 model — 3–4 states, not thirteen

- `P` PTF1A activity, autoregulatory, requiring an E-protein partner
- `R` RBPJL, produced **only** as a function of `P`. **No P-independent term. That
  zero *is* the bootstrap claim.**
- `C` slow chromatin/memory at metaplasia loci — sets time-to-relapse
- `E_free = E_total − k·ID3`, `ID3 = f(pERK)` — **algebraic, not differential**
- **pERK is an input, not a state.** Trametinib sets it.
- Profile likelihood on the three key parameters. **Not** an FIM eigenspectrum —
  sloppiness was a symptom of over-parameterization that no longer exists.

## Hard rules

**1. Pre-registration before any sweep or fit.** Ranges to
`prereg/<date>_<name>_ranges.yaml`, expected outcome to
`prereg/<date>_<name>_prediction.md`, **commit and push**, then run — the pushed
timestamp is the proof. **I write the templates; Luqmaan writes the content —
never fabricate his predictions.** Asked to skip this: refuse and say why.

**2. Nothing large or secret in git.** Datasets run to tens of GB; fetch by accession
with an idempotent script — the script is the artifact. Public repo: no tokens, no `.env`.

**3. Long jobs run detached, never as a Claude Code background process.** Bash caps
near 10 min and session-owned processes die on exit. Checkpoint to `results/<run_id>/`
every N samples, `--resume`, log to `logs/`, launch detached — designed up front.

**4. Architectural decisions get a 2-subagent adversarial panel**, then
`docs/decisions/NNN-name.md`: question, both positions, decision, **what would
reverse it**. Routine coding choices do not get a panel — just decide.

**5. The student must be able to explain every line.** Comment *why*, not *what*.
Stop and explain unfamiliar techniques (continuation, profile likelihood, Scheffé
polynomials, compound-Poisson dose, Bliss independence) in chat. Clear over clever.
**Push back when he's wrong.**

**6. Append to `docs/PROGRESS.md` every session** — done / broke / next / open
questions. Push.

## ⚠ STANDING RULE — never silently drop a failed solve

**Silently dropped failures are the single easiest way for this project to produce
a confident wrong answer**, because failures *correlate with swept parameters*.
Drop them and the surviving sample set is depleted precisely in the regime the
sweep was built to probe — which looks like a clean negative result. So everywhere
a solve happens (continuation, root-finding, profile likelihood, Phase 3's Monte
Carlo): **log every outcome**; **report failure rate as a function of every swept
parameter, never as a scalar**; **if it correlates with a swept parameter, say so
in the writeup** and treat results as *conditional on convergence*; and **keep a
test that sweeps deliberately in the hard regime** asserting a stated bound (<1%).

## ⚠ STANDING RULE — figure modules never compute science

**A figure module loads `results/`, does display arithmetic, and draws. It never
runs an ODE solve or a sample.** If a figure needs a number not in `results/`, the
fix is a stage that writes it. Consequence: rebuild is always under a minute,
figures cannot drift from the analysis, and a slow stage never blocks a figure.
**Produce figures as you go, every session — not in Phase 8.** One command:
`python make_figures.py`.

- **General-purpose geometry. No IEEE or conference-template sizing** — one `SCALE`
  constant in `figures/_style.py`, scaled downstream by Luqmaan.
- Okabe–Ito, **four nominal colours**: `#0072B2` acinar · `#D55E00` metaplastic ·
  `#009E73` intervention/success · `#6E6E6E` toxic — **grey, not red** (green/red is
  the textbook deuteranopia collision).
- `rc_context`, never a bare rcParams mutation — one process imports every module.
- `savefig.bbox="standard"`, **not** `"tight"` (tight silently resizes and breaks
  width checks). `pdf.fonttype 42`, `svg.fonttype "none"`.
- **Global grammar: a solid line is a model; an open marker with a dark edge is
  data. Never connect experimental points with a line.**
- Every stage calls `stamp_run()`; every figure declares inputs and `save_figure()`
  hashes them. `paper`/`poster` **refuse to render from a dirty tree**.
  `make_figures.py --check` fails on staleness.
- **Traps:** no 3-D cusp surface · no hairball of all states · no continuous colormap
  under a three-class map · no fake error bars on deterministic output · no p-value
  where there is no sampling · no truncated percentage axis · no bimodality from a
  KDE shape · two significant figures unless you can defend more.

## Cut — do not re-propose. Reasons in Part 4 and decision 012.

**Five-way topology competition** (its 10× discriminator is not measured — Krah
2019 has no timepoint between 24 h and 3 weeks; the bound is (1×, 21×)) ·
**13-state ODE** · **FIM sloppiness analysis** · **Lie-bracket ordering formalism**
(precedent exists — Letsou & Cai 2016, PMID 27560383 — so the novelty claim was
false) · **2^k subset enumeration** (the mass budget makes this a mixture problem;
enumeration visits only simplex vertices) · **submodular/greedy optimization**
(assumes diminishing returns; the bootstrap threshold claims *super*modularity) ·
**Pareto front as deliverable** · **Gillespie bimodality vs public scRNA-seq** ·
**CellOracle as validation** (keep as a *declared negative control*, pre-registered
to fail) · **Enformer/Borzoi · AlphaFold→k_on · Perturb-seq · one-sided U_crit ·
structural/Kalman controllability · full-dimension HJ reachability · MPC · all-atom
MD**. To revisit one, argue it in `docs/decisions/`. Don't just start doing it.

## Standing limitations — state these before anyone asks

- **The prioritization is a positive control, not a discovery.** Say it on the poster.
- **A regulon screen is blind to post-translational mechanisms.** TCF3/E47 is not
  transcriptionally lost — it is titrated — so no threshold surfaces it; it enters
  by **declared mechanism, registered before the screen runs.**
- **ID3→E47/PTF1A titration is documented in this cell line; the ERK→ID3 link is
  the assumption.** Dufresne 2010 (PMID 20830706), *in AR4-2J*: gastrin raises Id3,
  raises Id3/E47 *and* Id3/Ptf1-p48, lowers E47/Ptf1-p48, and Id3 silencing reverses
  the mislocalization. **v3 Part 1.4 understates this — see decision 012.**
  Unmeasured: any **Kd**, and whether **KRAS/ERK** rather than gastrin drives ID3.
- **PTF1A is pleiotropic** and dosage-sensitive (PMID 30470852). "Preserving
  viability" does not cover "did not make a neural-program-expressing cell" —
  **lineage fidelity is a separate endpoint.**
- **Converted cells often fail to silence the starting program** (CellNet, PMID
  25126793) — the candidate list is all *turn acinar back on*, hence the
  ADM-repressor axis.
- **AR42J needs dexamethasone** to express amylase at a differentiated level at all;
  baseline PTF1A/RBPJL and dynamic range must be measured or falsification power is
  unknown.
- **AR42J Kras genotype is unverified** — azaserine-induced rat tumours are
  classically Kras-mutant. Week-1 blocking item.
- **Reference data is mouse/human; the bench is rat, and a transformed line.**
  Cross-species transfer must be argued, not assumed.
- **No measured PTF1A half-life, no Hill coefficient, no Kd for any ID3
  interaction, no prior ODE model of ADM anywhere** — sampled with stated priors,
  never fitted point values.

## Framing — Part 6

Category **CBIO**. **Frame the wet lab as computation-that-designed-an-experiment,
never as validation** — if it comes back negative and you claimed validation, your
conclusion is deleted. Dated prediction box **left of** the outcome box;
pre-specified interpretation table; **effect size and pre-registered direction,
never a p-value**. Bound what it can falsify: ordering and timing yes, dose no.

## Environment

`venv\Scripts\activate` then `pip install -r requirements.txt`. Stack:
`numpy scipy matplotlib pandas sympy casadi SALib h5py anndata scanpy`.
Continuation is hand-rolled pseudo-arclength on a 3–4 state system — far smaller
than the old plan's, so the schedule risk that motivated the reduction is largely
gone. Do not adopt PyDSTool or AUTO-07p without a decision writeup.
