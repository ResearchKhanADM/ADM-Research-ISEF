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
**identity is not the claim**; it is published (Jiang 2023, PMID 37425649). The
claim is **composition, ratio, formulation and schedule under a fixed delivery
budget** — three results everything serves: **R1 formulation ★** (co-formulated vs
separate particles across the uptake-CV range; for an obligate stoichiometric pair
a requirement, not a preference) · **R2 composition** (count and ratio under fixed
total mRNA mass — a marginal-value curve, not a subset table) · **R3 durability**
(the drug-free persistence window — **what the bifurcation diagram is for**).

## Layout

**Committed:** `src/core.py` the 3-state model · `src/supplementary/` kept, not part
of the core · `scripts/` entry points · `prereg/` ranges + predictions written
*before* runs · `figures/` `_style.py`, `_provenance.py`, one module per figure,
`out/` holding final PDFs + `_source.csv` + `.prov.json` · `docs/decisions/` one file
per architectural decision · `tests/golden/` pre-rewrite fixture (**arbitrary
parameters — report no number from it**). **Gitignored:** `data/` `results/`
`logs/` `build/`.

## Phases — Part 2. Phase 1 runs parallel to Phase 2.

**0** decision spec, 1 wk — endpoints, viability floor X, reversal score; **issue
Bench Handshake #1, all 8 items blocking**; arm budget is an unknown **external**
input (placeholder 12/24/48 wells, design instantiable at any tier, say what
changes at each) · **1** candidate generation, 2 wk, **≥2 independent datasets** +
an ADM-repressor axis → **GATE A**, Layer 1 lock, order PTF1A + RBPJL mRNA ·
**2** minimal core, 4 wk → **GATE B** · **3 ★ delivery layer, 3 wk, the headline** →
formulation + ratio · **4** mixture-amount design, 2 wk → Layer 2 lock, count +
mass split · **5** durability + schedule, 2 wk (**pre-flight bench item 8 before
the ordering arms**) → redosing interval + order · **6** held-out prediction, 2 wk
→ **GATE C** · **7** prereg + discrimination power, 2 wk → locked protocol ·
**8** figures, writeup, buffer, 3 wk.

**Gates — stop and report. A failed gate is a real result, not a reason to push on.**
**A** does the pipeline recover PTF1A / RBPJL / BHLHA15? — if no, it is broken in
month one. **B** two stable states plus a saddle, identifiable separatrix, key
parameters surviving profile likelihood. **C** does the KRAS-history effect fall out
of a **single** parameter? **Never cut:** 0, 2, 3, 7. **Cut order under
compression:** 6 → Phase 4's simplex interior → Phase 5's ordering arms.

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

## ⚠ IF A CHECK COSTS LESS THAN THE ROUND-TRIP, RUN THE CHECK FIRST

**If a computation could make a Tier 3 question moot, it is Tier 1 — regardless of
what it is attached to.** Run it, then ask the question that survives. A question
Luqmaan answers that turns out to be moot costs both of us more than the check
did. Never write "cheaper than either fix, recommend running this first" and then
wait; that sentence *is* the instruction to run it.

## ⚠ TIME IS DIMENSIONLESS UNTIL BENCH ITEM 9 LANDS

Time is `τ = t·δ_P`, and `δ_P` — PTF1A protein turnover — is **absorbed by the
nondimensionalisation and unmeasured**. So **the model has no clock.**

**Report every timing in units of `1/δ_P`. Label every time axis `1/δ_P`. No
placeholder conversion to hours or days, not even in drafts.** Converting requires
a measured PTF1A protein half-life — **Bench Handshake item 9, now equal in
priority to item 8**: it is what makes the schedule deliverable quotable at all.
Item 4 (mRNA half-life) does not cover it — that constrains a product of
dimensionless groups, not the protein clock.

## ⚠ DECISION TIERS — reference: `docs/DECISION_PROTOCOL.md`

**T1 · decide alone, log one line in `PROGRESS.md`, keep going.** Implementation,
functional forms with a defensible default, naming, test design, refactors that
keep tests green, anything reversible in under a day. **Don't ask, don't batch —
just log.**

**T2 · convene a panel, decide, write the file, keep going.** Anything
architectural: state-space change · new or removed parameter · a changed
functional form **that alters a reported number** · method substitution ·
a conflict between two documents · **rewording a gate or a deliverable when the
correct wording, the reasoning and the cost are all in hand**. Spawn **2 subagents
with opposing mandates**
from `.claude/agents/` (`adversarial-reviewer`, `methods-checker`,
`literature-verifier` — the last **re-derives constants rather than inheriting
them**, standing instruction). Then `docs/decisions/NNN-name.md`, four headings,
*what would reverse this* checkable. Report it as **DECIDED-PENDING-REVIEW**.

**T3 · STOP AND ASK. Never automated.** (a) anything encoding Luqmaan's scientific
priority rather than a technical fact — what he trades off, what risk he accepts,
what an endpoint means · (b) anything consuming or constraining the one wet-lab
shot · (c) external facts unobtainable here (arm budget, bench measurements,
anything needing a person) · (d) wrong is expensive **and** irreversible · (e) any
claim of novelty or priority.
**BATCH THEM.** Don't stop at the first — carry on with everything unblocked and
present them numbered at session end, each with a recommendation and the cost of
each option. If one blocks *everything*, stop and **say it is a hard block**, and
what you tried first.

**The decisions folder is the interface** — Luqmaan reads files, not transcripts.
Every T2 file stands alone: stake, options, choice, reversal condition, downstream
effect. **If it needs the conversation to be understood, it is not finished.**

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

**4. The student must be able to explain every line.** Comment *why*, not *what*.
Stop and explain unfamiliar techniques (continuation, profile likelihood, Scheffé
polynomials, compound-Poisson dose, Bliss independence) in chat. Clear over clever.
**Push back when he's wrong.** **A retired claim is retired, not paused** — 002 and
008 were retired, not shelved; an unretired claim is one someone builds on later.

**5. Session end, always, in order:** what was done · **T2 decisions and why**
(DECIDED-PENDING-REVIEW, with links) · **batched T3 questions**, numbered, with
recommendations and costs · what broke · what is next · **anything that changed a
previously reported number** — not optional; the `n_eff` correction had reached
five documents before it was caught. Append to `docs/PROGRESS.md` and push.

**6. When Luqmaan hands you a framing, check it against the code before adopting
it.** The `a_P`→R2 / `γ`→R3 / `κ`→R1 mapping was contradicted by `core.py`'s own
`n_eff` docstring, in this repository, before it was written down. His framings
get the same verification as anyone's — that is what he asks for.

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
- **Global grammar: a solid line is a model; an open marker with a dark edge is data.
  Never connect experimental points with a line.**
- Okabe–Ito, **four nominal colours**, `#6E6E6E` toxic — **grey, not red**
  (deuteranopia). `rc_context` never a bare rcParams mutation.
  `savefig.bbox="standard"` **not** `"tight"`.
- Every stage calls `stamp_run()`; figures declare inputs and `save_figure()` hashes
  them. `paper`/`poster` **refuse a dirty tree**; `--check` fails on staleness.
- **House style, exact palette, and the trap list live in `figures/README.md`**,
  next to the code they govern. Read it before writing a figure module.

## Cut — do not re-propose. Reasons in Part 4 and decision 012.

**Five-way topology competition** · **13-state ODE** · **FIM sloppiness** ·
**Lie-bracket ordering formalism** · **2^k subset enumeration** ·
**submodular/greedy optimization** · **Pareto front as deliverable** · **U-shaped
viability hazard as an ODE** · **trametinib-vs-PD325901 asymmetry** · **Gillespie
bimodality vs scRNA-seq** · **CellOracle as validation** (survives only as a
*declared negative control*) · **Enformer/Borzoi · AlphaFold→k_on · Perturb-seq ·
one-sided U_crit · structural/Kalman controllability · full-dim HJ reachability ·
MPC · all-atom MD**.
**Each has a specific reason in v3 Part 4 — read it before arguing, they are not
interchangeable.** To revisit one, argue it in `docs/decisions/`. Don't just start.

## Standing limitations — say these first. Full wording: v3.1 Part 1.4

- **ID3 titration is DOCUMENTED in this cell line — the gap is ONE EDGE, ERK→ID3.**
  Dufresne 2010 (PMID 20830706) in AR4-2J; its driver is **gastrin**, so KRAS/ERK→ID3
  is the inference, testable in one western (bench item 8). **No Kd** measured.
  **Don't disclaim the node itself** — that error was already made once.
- **The prioritization is a positive control, not a discovery.** And a regulon
  screen is blind to post-translational mechanisms: TCF3/E47 is titrated, not
  transcriptionally lost, so it enters by **declared mechanism registered before
  the screen runs**. Converted cells often fail to silence the starting program
  (CellNet, PMID 25126793) — hence the ADM-repressor axis.
- **PTF1A is pleiotropic** and dosage-sensitive (PMID 30470852) — **lineage fidelity
  is a separate endpoint** from viability.
- **AR42J needs dexamethasone** to express amylase at all, and **its *Kras* genotype
  is unverified**; both week-1 blocking, and without dynamic range falsification
  power is unknown. **Reference data is mouse/human, the bench is rat and
  transformed** — argue the transfer. **No measured PTF1A half-life, Hill
  coefficient, ID3 Kd, or prior ODE model of ADM** — sampled priors, never fitted
  points.

## Framing — the rest is v3 Part 6

Category **CBIO**. **Frame the wet lab as computation-that-designed-an-experiment,
never as validation** — claim validation, get a negative, conclusion deleted.
**Effect size and pre-registered direction, never a p-value.** Bound what it can
falsify: ordering and timing yes, dose no.

## Environment

`venv\Scripts\activate`, `pip install -r requirements.txt`. Continuation is
hand-rolled pseudo-arclength on 3 states; no PyDSTool or AUTO-07p without a writeup.
