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

## Already ruled out — do not re-propose

Killed by expert panel with reasons in the master plan: AlphaFold/Boltz → binding rate constants · Enformer/Borzoi (no rat support, wrong output class) · Perturb-seq vector arithmetic for cocktail derivation · fast histone acetylation as a standalone state · one-sided toxicity ceiling · structural/Kalman controllability · full-dimension Hamilton–Jacobi reachability, MPC, all-atom MD.

To revisit one, argue it in `docs/decisions/` against the plan. Don't just start doing it.

## Stage order — sequential, each depends on the last

0. **Reduce and certify** — nondimensionalize, eliminate fast variables to 5–6 slow states, classify all fixed points, locate the separatrix. *Never skip.*
1. **Two-parameter bifurcation** — KRAS dose × trametinib dose.
2. **Topology competition** — 5 architectures, one sampling box, compare Q-values. **GATE: does any topology reproduce the 3-day (MEKi) vs 3-week (forced PTF1A) asymmetry?**
3. **Identifiability + third-mRNA selection** — FIM eigenspectrum; rank E47/TCF3 vs NR5A2 vs MIST1.
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
