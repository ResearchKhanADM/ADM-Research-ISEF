# ADM Reversal — in-silico project

**Goal: reverse KRAS-driven acinar-to-ductal metaplasia while keeping the cells alive.**
The wet-lab experiment runs in the same cycle and gets **one shot, no iteration** — this model exists to design it, and must produce answers in time to inform it.

Full plan: `docs/ADM_INSILICO_MASTER_PLAN.md`. **Read it before proposing anything.**
Running log: `docs/PROGRESS.md`. Read the last entry at the start of every session.

Repo: https://github.com/ResearchKhanADM/ADM-Research-ISEF (public)

## Layout

| Path | Contents | Committed? |
|---|---|---|
| `src/` | model + analysis code | yes |
| `scripts/` | runnable entry points (sweeps, downloads) | yes |
| `prereg/` | parameter ranges + predictions, written *before* runs | yes |
| `figures/` | board figures | yes |
| `docs/decisions/` | one file per architectural decision | yes |
| `data/` `results/` `logs/` | datasets, sweep output, checkpoints | **no — gitignored** |

## Hard rules

**1. Pre-registration before any sweep or fit.** Write ranges to `prereg/<date>_<name>_ranges.yaml` and the expected outcome to `prereg/<date>_<name>_prediction.md`, then **commit and push**, then run. The value of the held-out prediction depends on proving it was recorded first, and a pushed commit timestamp is that proof. If asked to run a sweep without this, refuse and say why.

**2. Nothing large or secret in git.** Datasets run to tens of GB. `scripts/download_data.py` fetches by accession and is idempotent — the script is the artifact. Repo is public: no tokens, no `.env`.

**3. Long jobs run detached, never as a Claude Code background process.** Bash caps near 10 minutes and session-owned background processes are killed on exit. Sweeps must checkpoint to `results/<run_id>/` every N samples, support `--resume`, log to `logs/`, and launch detached. Design checkpointing before the first long run.

**4. Architectural decisions get a 2-subagent adversarial panel**, then a writeup in `docs/decisions/NNN-name.md` stating the question, both positions, the decision, and **what would reverse it**. Routine coding choices do not get a panel — just decide.

**5. The student must be able to explain every line.** Comment *why*, not *what*. Stop and explain unfamiliar techniques (continuation, Sobol indices, Lie brackets, reach-avoid sets) in chat. Prefer clear over clever. Push back when he's wrong.

**6. Append to `docs/PROGRESS.md` at the end of every session** — done / broke / next / open questions. Push.

## ⚠ `W` IS PROTECTED FROM ELIMINATION

The model has **11 states**, not 10. `W` (phospho-MEK pool) is **fast by turnover — minutes — and must NOT be eliminated in the Stage 0 fast-variable sweep.**

QSS on `W` substitutes it straight back into `K_eff` and recovers the static product `K_eff = K·f_act·f_cat`, which provably cannot produce the §1.2 withdrawal-asymmetry prediction: at matched pERK both drugs give identical `K_eff`, so all states evolve identically and withdrawal is identical. **Eliminating `W` is not an approximation of this model — it is a different model making a different prediction.**

Two things that must not drift:
- **`RAF_drive` is strictly DECREASING in `K_eff`.** Falling ERK relieves negative feedback on RAF, so drive rises as ERK falls. An increasing implementation inverts the mechanism, predicts the opposite, and still runs clean. Asserted in code, covered by a unit test.
- **`τ_W` is sampled across minutes–hours**, never fixed. It sets the impulse `∫ΔK_eff dt`, which decides whether the separatrix is crossed.

The §1.2 prediction is **conditional**: asymmetry occurs iff the transient overshoot suffices to cross the separatrix. Report the ensemble fraction. 5% is a result; say 5%.

This is the one deliberate fast-variable retention in the model. It will read as an inconsistency to anyone who does not find the justification, so it must be argued in the writeup. Reasoning: `docs/decisions/002-w-state-protected-from-elimination.md`.

## What the mechanistic model may and may not claim

The payload species **are** state variables, so the ODE cannot derive payload *identity* without circularity. It is never asked "which molecules." It is asked:

- **necessity** — all 16 subsets of {trametinib, u₁, u₂, u₃} across the ensemble; which components are required, which are redundant. Compare at matched **total** dose as well as matched per-component dose, or the result is just "more protein is better." "Three mRNAs is over-engineering" is a valid finding.
- **dose, schedule, order** — Stages 5 and 6.
- **whether the reversal-optimal and viability-optimal payloads are the same payload** — a Pareto question. If they differ, that is a headline result and it is the stated objective.

Identity comes from the other arm, Stage 3B. Reasoning in `docs/decisions/001-two-arm-payload-derivation.md`.

## Already ruled out — do not re-propose

Killed by expert panel with reasons in the master plan: AlphaFold/Boltz → binding rate constants · Enformer/Borzoi (no rat support, wrong output class) · Perturb-seq vector arithmetic for cocktail derivation · fast histone acetylation as a standalone state · one-sided toxicity ceiling · structural/Kalman controllability · full-dimension Hamilton–Jacobi reachability, MPC, all-atom MD.

**Perturb-seq ≠ Stage 3B.** Perturb-seq was killed on *cell context* — K562/RPE1 do not express RBPJL, so the libraries could not contain it. That objection is about the data source and does not transfer to CellOracle on a pancreatic dataset, which has the factors present and uses a different method. Stage 3B carries its own separate limitations; they are listed in the master plan.

To revisit one, argue it in `docs/decisions/` against the plan. Don't just start doing it.

## Stage order — sequential, each depends on the last

0. **Reduce and certify** — nondimensionalize, eliminate fast variables to 5–6 slow states, classify all fixed points, locate the separatrix. *Never skip.*
1. **Two-parameter bifurcation** — KRAS dose × trametinib dose.
2. **Topology competition** — 5 architectures, one sampling box, compare Q-values. **GATE: does any topology reproduce the 3-day (MEKi) vs 3-week (forced PTF1A) asymmetry?**
3. **Identifiability + third-mRNA selection** — FIM eigenspectrum; rank E47/TCF3 vs NR5A2 vs MIST1.
3B. **In-silico perturbation screen** — CellOracle (PMID 36755098) on GSE207938; overexpress across the full TF repertoire, rank by ADM→acinar shift. **Pre-registration is non-negotiable here**, including predicted ranks for PTF1A/RBPJL/NR5A2. Whatever it returns gets reported, including a low rank for those three. Pre-flight: confirm RBPJL has outgoing edges in the base GRN, or its score is meaningless rather than negative.
4. **Held-out prediction** — **GATE: does the KRAS-history effect fall out of one parameter?**
5. **Dosing schedule** — (dose per mRNA pulse) × (redosing interval), three-region classification.
6. **Ordering** — Lie bracket, then isodose phase sweep.
7. **Single-cell falsification** — Gillespie vs GSE314765 / GSE207938 / GSE172380 / GSE141017.

**Stop at each gate and report before continuing.** A failed gate is a real result, not a reason to push on.

## Environment

```
venv\Scripts\activate                  # Windows
pip install -r requirements.txt
```
Core stack: `numpy scipy matplotlib pandas sympy casadi SALib h5py anndata scanpy`.
Continuation is hand-rolled pseudo-arclength on the *reduced* system (~200 lines). Do not adopt PyDSTool (unmaintained, poor on Windows) or AUTO-07p (Fortran build) without a decision writeup — continuation on the full-dimension system is the single biggest schedule risk in the project.

## Standing unknowns

- **AR42J Kras genotype is unverified.** The line came from an azaserine-induced rat tumour and those are classically Kras-mutant. If it is already mutant, "KRAS-induced ADM" changes to "KRAS dose titration on a mutant background" and the forcing term changes. Highest-priority open question.
- No measured PTF1A half-life, no measured Hill coefficient, no Kd for any ID3 interaction, no prior ODE model of ADM anywhere. Treat these as sampled parameters with stated priors, never as fitted point values.
