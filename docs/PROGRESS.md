# Progress log

Append one entry per session. Newest at the bottom, so `git log` and this file
read in the same direction. Every entry: **done / broke / next / open questions.**

---

## 2026-08-14 · Session 1 — setup only, no science

### Done

- Project created at `C:\Users\luqma\adm-insil-des`. No spaces in the path, by design.
- Git initialised, `core.autocrlf true`, remote `origin` →
  `https://github.com/khanluqmaanresearch-byte/ADM-Research.git`, branch `main`.
- Directory structure built per the kickoff spec. `data/`, `results/`, `logs/`
  exist locally but are gitignored, so **a fresh clone will not have them** —
  anything that writes there must `mkdir` first.
- `.gitignore` written to spec. Checked before the first push: no `.env`, no
  keys, no credentials, no datasets staged. The repo is public.
- venv on **Python 3.13.9** (the standalone install at
  `AppData\Local\Programs\Python\Python313`, deliberately not the miniconda base
  — mixing conda and venv on Windows is a debugging tax nobody needs in month 4).
- All eleven required packages installed **and imported**. CasADi and scanpy both
  had native Windows wheels; no substitutions, no workarounds.
- `scripts/verify_env.py` — repeatable import check, exits 1 on failure. It exists
  because "pip said success" and "the import works" are different facts on
  Windows, where CasADi ships compiled binaries that can fail on a missing DLL.
- `requirements.txt` (floors) + `requirements-lock.txt` (exact, verified). The
  lockfile is the one that matters in five months.
- Master plan copied to `docs/ADM_INSILICO_MASTER_PLAN.md`; `CLAUDE.md` at root.
- `docs/decisions/000-TEMPLATE.md` seeded.

### Verified versions

python 3.13.9 · numpy 2.5.2 · scipy 1.18.0 · matplotlib 3.11.1 · pandas 3.0.5 ·
sympy 1.14.0 · casadi 3.7.2 · SALib 1.5.2 · h5py 3.16.0 · anndata 0.13.2 ·
scanpy 1.12.3 · pyyaml 6.0.3

### Broke

- Nothing blocking. Two things worth recording:
  - The GitHub repo already carried an `Initial commit` with a one-line stub
    README. Built **on top of it** rather than force-pushing over it — the
    history is the evidence trail and overwriting it would have cost the
    timestamp on the first commit for no gain.
  - `verify_env.py` first reported SALib as version "unknown" and raised
    FutureWarnings from anndata/scanpy, because it read `module.__version__`.
    Rewritten to use `importlib.metadata`. Fixed, not suppressed.
- **The push did not complete.** `git push` failed with
  `fatal: User cancelled dialog` → `could not read Username for
  'https://github.com'`. The credential helper is Git Credential Manager, which
  needs to open an interactive auth window, and the session running this had no
  TTY. **Nothing is wrong with the repo** — the commit is intact and `main` is
  simply ahead of `origin/main` by one. It needs one interactive
  `git push -u origin main` from a normal terminal; GCM then caches the
  credential and later sessions push without prompting.
  This matters beyond convenience: **pre-registration is enforced by a pushed
  commit timestamp**, so pushing must be known to work before the first sweep,
  not discovered to be broken at that moment.

### Next

- **Stage 0 · Reduce and certify.** Per master plan §3.4 the deliverables are:
  explicit functional forms for `Hill(·)`, `g(K,v)`, `φ(A)`, `f_act`, `f_cat`;
  the E-protein binding polynomial solved; nondimensionalisation (~30 parameters
  → ~15–18 groups); QSS elimination of `P_c`, `E`, `C_L`, `C_J` **with the
  algebra showing it is valid**; the parameter table with ranges and sources;
  ADM and acinar initial conditions; solver settings.
- Ambiguous functional forms go to `docs/decisions/`, not to a blocker.
- No sweep may run until its `prereg/` ranges and prediction are committed **and
  pushed**. That is not a formality — the headline validation in Stage 4 is a
  held-out prediction, and its value is the pushed timestamp.

### Open questions — need Luqmaan

1. **AR42J Kras genotype. Highest-stakes unknown in the project** (master plan
   §7.1). AR42J came from an azaserine-induced rat tumour and those are
   classically Kras-mutant. If the stock is already mutant, inducing G12D on top
   is close to meaningless: the forcing term changes and the project becomes
   "KRAS dose titration on a mutant background". This changes Stage 0's forcing
   input, so it wants resolving early — not before Stage 0 starts, but before
   Stage 1 fixes the KRAS axis. Resolve by sequencing the stock.
2. **Two superseded plan drafts are still on this machine** —
   `INSILICO_PLAN_8months.md` and `INSILICO_PLAN_v2.md`, in the Claude outputs
   folder. The master plan says explicitly not to read or merge them. Worth
   deleting or moving them so a future session cannot pick one up by accident.
3. **No `gh` CLI installed.** Not needed for anything so far; push works over
   HTTPS. Flagging only so it is not a surprise later.

---

## 2026-08-15 · Session 2 — Stage 0, part 1: write the unreduced system

### Housekeeping

- **The GitHub remote moved.** Old:
  `khanluqmaanresearch-byte/ADM-Research`. New:
  **`ResearchKhanADM/ADM-Research-ISEF`**. Luqmaan ran `git remote set-url`,
  corrected `user.name` / `user.email`, and force-pushed before this session.
  Verified here: `git fetch` then compared SHAs directly rather than trusting
  the tracking ref — local and remote `main` are both `40bdaa4`. In sync, and
  both session-1 commits survived the force-push. The session-1 log entry above
  still names the old URL **on purpose**: it records what was true that day, and
  a log that gets edited to match the present is not a log.
- Repo URL updated in `CLAUDE.md` and `README.md`.
- PubMed and bioRxiv MCP servers confirmed live by real calls, not just by
  loading schemas. PubMed round-tripped PMID 24315826 and returned the exact
  citation the master plan gives for Collins 2014 — *Gastroenterology*
  146(3):822–834.e7, doi 10.1053/j.gastro.2013.11.052. That is one independent
  check on the plan's most load-bearing reference.

### Plan change — two-arm payload derivation (Luqmaan's call, before any code)

Caught a real problem in the plan before Step 1 was written: **the payload
species are state variables, so asking the ODE "which molecules?" is circular.**
Fixed by splitting the question across two arms rather than by rewording it.

- **Mechanistic arm reframed** to three non-circular questions — necessity
  (16 subsets of the four interventions), dose/schedule/order, and whether the
  reversal-optimal and viability-optimal payloads are the *same* payload.
- **New Stage 3B** — CellOracle in-silico perturbation screen (PMID 36755098) on
  GSE207938, full TF repertoire, ranked output, pre-registered. Placed after
  Stage 3 to keep the two arms independent; shares datasets with Stage 7;
  2nd in cut order.
- Both PMIDs verified against PubMed before being written into the plan:
  36755098 = *Nature* 614(7949):742–751 (2023); 36584685 = *Stem Cell Reports*
  18(1):97–112 (2022), which rescued a failing conversion by nominating Fos and
  Yap1 — factors the field did not expect, which is what makes it a real
  precedent.
- Reasoning and reversal conditions: `docs/decisions/001-two-arm-payload-derivation.md`.

Three things added that were not requested, because they are where this gets
attacked:

1. **The RBPJL-motif risk.** CellOracle propagates through motif-derived edges;
   a TF with no motif has no outgoing edges and scores ≈ 0 *silently*, looking
   like "RBPJL doesn't matter" rather than "RBPJL is not representable." RBPJL is
   the one component the whole argument rests on. Pre-flight check, not a
   post-hoc caveat.
2. **A fair-dose rule for the necessity analysis.** Dropping a component changes
   total delivered material, so subsets must be compared at matched *total* dose
   as well as matched per-component dose — otherwise the analysis rediscovers
   "more protein is better." Same discipline Stage 6 already imposes on ordering.
3. **An explicit Perturb-seq vs Stage 3B note** in both the plan's ruled-out
   table and `CLAUDE.md`. Without it, the obvious reading is that a killed module
   was smuggled back in.

### Two plan-level errors caught by reading before coding

Both found by checking the plan against itself, and both would have been
expensive later.

1. **The trametinib model could not produce its own free prediction.** §3.2's
   static `K_eff = K·f_act(v)·f_cat(v)` gives both drugs an identical `K_eff` at
   matched pERK, so every state evolves identically and withdrawal is identical.
   A static product has nowhere to store the difference. Fixed by adding `W`
   (phospho-MEK) as an 11th state, with `f_act` blocking **formation** and
   `f_cat` blocking **output**. Luqmaan's two corrections on top: `RAF_drive`
   must be *decreasing* in `K_eff` (an increasing version inverts the mechanism
   and still runs clean), and `W` is *fast*, so the prediction survives only
   through bistability — the overshoot need not persist, it need only carry the
   state across the separatrix. Prediction restated conditionally, ensemble
   fraction to be reported. Decision 002.
2. **§3.2 violated its own constraint 2**, carrying `−k_seq·I·P_n` two
   paragraphs above the text explaining why a first-order sink is wrong. My
   proposed fix — delete the term — was itself wrong, and Luqmaan caught it:
   ID3 genuinely does sequester PTF1A, so sequestration had to move **into the
   shared equilibrium**, not disappear. Deleting it would have been the opposite
   error.

### Stage 0 · Step 1 — the full unreduced system

`src/functional_forms.py`, `src/topology.py`, `src/binding.py`, `src/model.py`,
plus 15 passing invariant tests. All seven topologies assemble and evaluate.
Nothing reduced — nondimensionalisation and QSS are next session.

**Two findings from actually writing it:**

1. **Titration is ultrasensitive only in the tight-binding regime, and this
   constrains Stage 2's sampling box.** Effective Hill coefficient
   `n_eff = max|d ln C_L / d ln I|`:

   | regime | titration | first-order (T2) | ratio |
   |---|---|---|---|
   | `Kd=1.0` loose | 2.07 | 1.94 | 1.07× |
   | `Kd=0.01` tight | 13.40 | 2.00 | 6.70× |
   | `Kd=0.001` tight | 41.66 | 2.00 | 20.83× |

   The first-order sink sits at exactly 2.0 throughout — its analytic slope.
   **If the `Kd` prior does not reach `Kd` << protein totals, T1 and T2 are
   indistinguishable and the topology competition silently fails to
   discriminate**, looking like "the data cannot separate these architectures"
   when the truth is "the box never sampled where they differ." No `Kd` has ever
   been measured for any ID3 interaction, so this prior is a free choice that
   determines the answer. **Must be pre-registered deliberately, not defaulted.**

2. **I measured ultrasensitivity wrongly first.** My initial test used global
   fold-change across a fixed ID3 range and ranked the first-order sink as
   *sharper*. That metric is wrong: titration is threshold-linear — steep near
   the threshold, flat away from it — so a global fold-change averages away the
   thing being measured. Corrected to the local effective Hill coefficient, and
   the reasoning is written into the test so it is not repeated.

Also fixed a real solver bug the test exposed: the log-space root find failed to
converge in the tight-binding regime, i.e. exactly where the mechanism lives.
Replaced with bounded multi-start least-squares on totals-scaled residuals.
Tolerance is 1e-6 relative rather than machine-tight, because deep tight binding
requires catastrophic cancellation to do better and four-decade priors cannot
justify chasing it — Step 2's explicit polynomial removes the cancellation
instead of out-iterating it.

**Decision files written:** 001 two-arm payload · 002 `W` protected · 003
composable topologies · 004 Hill forms · 005 ERK/drug layer · 006 two-target
titration · 007 `dI/dt` and `dA/dt` · 008 viability hazard · 010 NR5A2
placement · 011 pulse forcing.

**Stage 6 precondition verified symbolically rather than assumed:** `[g₁,g₂]=0`
holds, a time-varying coefficient preserves `∂g/∂x = 0`, and
`[f_A,f_B] = J·(g₁−g₂)` exactly. Caveat recorded: the `s²` expansion assumes
autonomous fields, so the closed form applies to constant-amplitude holds and
the pulse shape belongs to Stage 6's simulation arm.

### Next

- **Step 2** — solve the binding polynomial explicitly, every step, to
  `docs/derivations/binding_polynomial.md`. It also removes the numerical
  cancellation problem found above.
- **Step 3** — parameter table. Expect the "unmeasured — sampled" column to be
  long (Part 8). **The ID3 `Kd` prior is now known to be decisive** and needs
  its width justified explicitly.
- **Step 4** — initial conditions and solver settings; report the observed
  stiffness ratio.
- Then nondimensionalisation and QSS — with `W` **exempt**.

### Open questions

1. **AR42J Kras genotype** — unchanged, still the highest-stakes unknown.
2. **`A` → "% amylase-positive cells" map is not yet defined.** Collins reports a
   *fraction of cells*; `A` is a concentration. Comparing them directly is a
   category error that would contaminate every Part 6 validation target. Needs
   settling in Step 3/4 as a thresholded single-cell readout, not a rescaled
   population mean — the same wrong-observable failure that killed PU.1/GATA1.
3. **Luqmaan's CORE list omitted `I` (ID3).** I included it — the unreduced
   system needs the titrator as a state — and flagged it as the leading QSS
   candidate, which is how the slow count still reaches 6. Flagging in case the
   omission was deliberate.
