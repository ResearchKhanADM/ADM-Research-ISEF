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

---

## 2026-08-15 · Session 3 — architecture change absorbed; documents and cleanup only

No science code, no modelling, by instruction. The model rewrite is next session.

### 🚨 BLOCKING — Week-1 bench items. Phase 3 cannot run without these.

**These are blocking, not "nice to have", and they have long lead times.** Phase 3
is the headline result and every one of its axes has units the bench must supply;
Phase 5's dose × interval map is unplottable without PK. Issue **Bench Handshake
#1 now** (v3 Phase 0), in week 1, not when Phase 3 starts.

| # | Item | Blocks | Why it cannot be assumed |
|---|---|---|---|
| 1 | **AR42J *Kras* genotype** — sequence codons 12/13 and 61 | the entire forcing term | Azaserine-induced rat tumours are classically Kras-mutant. If the stock is already mutant, "KRAS-induced ADM" becomes "KRAS dose titration on a mutant background" and the input changes shape. **Highest-stakes unknown in the project, unchanged since session 1.** |
| 2 | **Trametinib IC50 *in these cells*** | Phase 2's drug axis; Phase 5 | A published IC50 from another line sets the wrong axis origin, and every dose recommendation inherits the error invisibly. |
| 3 | **LNP transfection efficiency AND its cell-to-cell CV** | **Phase 3 ★, directly** | The CV *is* the x-axis of the headline figure. Efficiency alone is not enough — the whole co-formulation result is a statement about the *spread*, so a mean without a CV cannot produce it. |
| 4 | **mRNA half-life in AR42J / Matrigel** | Phase 5's redosing interval | The interval axis has no units without it. Matrigel-specific: 3D culture is not 2D and a plastic-dish number is a different number. |
| 5 | **Baseline PTF1A / RBPJL by qPCR, ± dexamethasone** | Phase 2 initial conditions; Gate A interpretation | AR42J needs dex to express amylase at a differentiated level at all. Without baseline and achievable dynamic range, **the experiment's falsification power is unknown** — a null could mean "payload failed" or "there was no room to move". |
| 6 | **ADM stability to ≥14–21 days** | Phase 5's durability endpoint | If the ADM state does not hold that long unassisted, "durable reversal" has no measurable contrast. |

Item **3** is the one to chase hardest: it is the only input to the project's
strongest result, and it is the one no literature value can substitute for.

### Done

- **Read `ADM_MASTER_PLAN_v3.md` end to end** and audited it against the repo.
  Conflicts reported to Luqmaan before anything was touched; the material ones are
  recorded in decision 012 rather than left in chat.
- **Deleted the superseded plans**, from the repo and from the Claude outputs
  folder on this machine: `docs/ADM_INSILICO_MASTER_PLAN.md`,
  `INSILICO_PLAN_v2.md`, `INSILICO_PLAN_8months.md`. Git history preserves them.
  This closes session 1's open question 2 — they had been read by mistake once and
  nearly built from.
  *(The two AppData paths that appeared to hold separate copies —
  `Roaming\Claude\...` and `Local\Packages\Claude_pzs8sxrjxfjjc\LocalCache\Roaming\Claude\...`
  — are the **same directory** behind the packaged-app VFS redirect. Verified with
  a probe file rather than assumed, so one delete genuinely covered both.)*
- **v3 copied into `docs/ADM_MASTER_PLAN_v3.md`** as the single plan of record,
  hash-verified identical before anything was deleted.
- **`CLAUDE.md` rewritten** against v3 — eight phases, three gates, the cut list,
  the standing limitations, and the Part 5 figure rule. 198 lines.
- **`CLAUDE_CODE_KICKOFF_PROMPT.md` fixed** (it lives in the outputs folder, not
  the repo): banner saying it is historical, the seven-stage work order replaced
  by the eight-phase table, and the one-line "how to start a session now".
  Also banner-marked the **stale `CLAUDE.md` copy** sitting in the same folder —
  it named the old GitHub URL and the deleted plan.
- **Decision 012 written** — the architecture change, what it supersedes, and what
  would reverse it.
- **`README.md` rewritten** — it pointed at the deleted plan and described five
  candidate architectures.
- **Superseded banners** on `docs/STAGE0_PLAN_GAPS.md` (every `§` in it points
  into a deleted file) and `prereg/id3_kd_prior_justification.md` (its Stage 2
  destination is gone, but its physics feeds Phase 3).
- **Figure infrastructure built** per v3 Part 5, and **verified end to end**:
  `figures/_style.py` (one `SCALE`, four Okabe–Ito roles, model-vs-data grammar,
  `savefig.bbox="standard"`), `figures/_provenance.py` (`stamp_run`,
  `save_figure`, `is_stale`), `make_figures.py`, `figures/README.md`, `build/`
  gitignored. Confirmed by running: `--selftest` renders a proof sheet with a real
  `.prov.json` and `_source.csv`; `profile="paper"` **refuses a dirty tree**;
  `is_stale` correctly reports a never-rendered figure. No figure modules yet —
  there are no `results/` to draw from, and a figure module that computes its own
  numbers is the one thing the architecture forbids.

### Broke

- Nothing broke. One thing was **found wrong, in v3 itself** — see open question 1.
- Worth recording as a near-miss rather than a break: the plan-of-record file and
  three superseded drafts were all sitting in the same folder, and the folder is
  reachable by any session. That is how the wrong plan got read the first time.

### Next — the model rewrite

- **Phase 2 core: 3–4 states, ~9–12 parameters, nondimensionalized.** `P`, `R`,
  `C`, with `E_free = E_total − k·ID3` and `ID3 = f(pERK)` algebraic, and **pERK an
  input, not a state.**
- Before writing it, decide what happens to `src/` (audit is in this session's
  report — nothing has been moved or deleted). `binding.py` and the derivation
  should be **kept as supplementary**, not deleted: they quantify how *sharp* the
  bootstrap threshold is, and Phase 3 convolves the LNP dose distribution against
  exactly that threshold.
- Profile likelihood on three parameters. **Not** an FIM eigenspectrum.
- Produce figures as the stages produce results, not in Phase 8.

### Open questions — need Luqmaan

1. **v3 Part 1.4 is wrong about the ID3→E47 titration node, and it should not go
   on a poster as written.** It says targeted search returned no evidence in the
   pancreatic context. **Dufresne 2010 (PMID 20830706) is exactly that evidence,
   in AR4-2J — the bench line.** Verified against PubMed this session: gastrin
   raises Id3, raises Id3/E47 **and** Id3/Ptf1-p48, lowers E47/Ptf1-p48; Id3
   silencing reverses the mislocalization; Id3 is overexpressed with Ptf1-p48
   mislocalized in human and murine preneoplastic lesions. What *is* undocumented
   is narrower and should be stated in that narrower form: **no Kd**, and the
   driver in Dufresne is **gastrin, not KRAS/ERK — the ERK→ID3 edge is the
   assumption.** This makes Phase 5's ordering prediction better supported than v3
   claims. **Confirm you want v3 Part 1.4 amended**; I have not edited the plan.
2. **`W` and the trametinib-vs-PD325901 withdrawal asymmetry are being dropped by
   architecture, not by evidence.** Decision 002's own reversal condition — a
   Stage 1 finding that the attractor sits too far from the separatrix — was never
   tested. v3 makes pERK an input, which leaves nowhere for `W` to live. Confirm
   the prediction is **set aside for schedule reasons, not refuted**; that is how
   decision 012 records it.
3. **Does viability stay modelled, or become a bench measurement plus a floor?**
   v3 keeps a viability floor X (Phase 0) and the therapeutic index (Phase 7) but
   **never mentions the U-shaped hazard** — neither keeping nor cutting it. At 3–4
   states there is no `S` (secretory cargo) state, so the hazard's high-`S` arm has
   nothing to attach to. Decision 008 is therefore in limbo. **Needs a call at
   Phase 0**, because it defines an endpoint.
4. **What is the arm budget?** v3 Phase 0's first question is *how many
   experimental wells actually exist*. Everything in Phases 4, 5 and 7 is a
   function of that number, and the previous plan's mistake was optimizing before
   knowing it.

---

## 2026-08-15 · Session 4 — six calls absorbed; golden pinned; parameter budget

All four open questions from session 3 answered by Luqmaan, plus two he raised.
Documents and fixtures done; **the RHS is deliberately not written yet** — the
parameter budget is a gate.

### Done

- **Golden trajectory pinned before anything is deleted** —
  `scripts/pin_golden_trajectory.py` → `tests/golden/`. Two initial conditions
  (acinar, ADM) integrated on the 11-state system for 1500 h, with a manifest,
  hashes and a `stamp_run()`. **The old implementation has two attractors**, and
  the ADM one sits at `R ≈ 0` — the bootstrap claim showing up in the fixture
  rather than in a paragraph. That is the qualitative property the rewrite gets
  compared against; it will not match numerically and is not meant to.
  **The parameter set is arbitrary and unmeasured** — labelled as such in the
  script, the manifest and `CLAUDE.md`. No number from it may ever be reported.
- **v3 amended in place to v3.1** with calls 1–4, plus an amendment banner at the
  top so the edits are visible rather than silent. Bench Handshake #1 is now a
  **table of 8 blocking items** with what each one blocks.
- **Decision files amended, originals preserved underneath** — deleting a decision
  makes its reversal unauditable:
  - **002** — `W` removed; the phenomenon becomes a **swept pERK rebound profile**;
    trametinib-vs-PD325901 **retired outright**.
  - **006** — **promoted.** It is no longer a topology discriminator; it is the
    justification for v3's `E_free`, and its `n_eff` is Phase 3's threshold sharpness.
  - **008** — **retired.** Viability is a bench-measured floor; the CHOP arm
    survives as an output flag, not a term.
  - **012** — the six calls recorded as resolutions, including the Pareto cut as an
    **explicit reversal** rather than an omission.
- **`CLAUDE.md` updated and still under 200 lines** (199): exact `E_free` with no
  floor hack, the pERK rebound profile, viability retired, the flagship guard test,
  the arm-budget tiering, and hard rule 7 — *a retired claim is retired, not paused*.
- **`docs/PHASE2_PARAMETER_BUDGET.md` written** — the gate. See below.

### The number, before the RHS is written

**11 states → 3. 61 parameters → 11 dimensionless groups**, of which 3 are Hill
exponents that are scanned rather than fitted, so the effective fitting dimension
is **8**. Inside v3's 9–12 target.

Profile likelihood goes on the three that carry the results: **`a_P`** (does the
loop close), **`γ = δ_C/δ_P`** (how long it holds — the durability knob), **`κ =
K_d/E_tot`** (threshold sharpness, which feeds Phase 3).

`κ` surfacing as a group is a **consistency check that passed**:
`prereg/id3_kd_prior_justification.md` predicted the nondimensionalization would
produce `K_d/E_tot` as a single group, and said that if it did not, the
nondimensionalization was wrong. It did.

Four judgement calls are flagged in that file rather than buried — one `K_d`
instead of two, `n_C` possibly removable (→ 10 groups), `A` not a state, and
`ε_C`/`θ_C` as the least-constrained pair.

### Broke

- Nothing. One sequencing note: **`topology.py` and `model.py` are still present.**
  `model.py` imports `topology.py`, so deleting either alone breaks the other, the
  test suite, and the golden-trajectory script that was just built. They come out
  **in the same commit as the new core**, not before it. 20 tests still pass.

### Next

1. Luqmaan signs off on the 11 groups (or cuts `n_C` to 10).
2. Write `src/core.py` — the 3-state RHS — with the flagship guard test **first**:
   *`dR/dt` has no P-independent term*, as a guard-on-the-guard, plus the dynamic
   companion (`R` must not rise from zero while `P` is held at zero).
3. Same commit: delete `topology.py`, `model.py`, `payload_subsets()` and the tests
   that die with them; move `binding.py` and its derivation to supplementary.
4. Nondimensionalize in code, not just on paper. Then Gate B: two stable states plus
   a saddle, identifiable separatrix.

### Open questions

1. **Sign-off on the parameter budget** — 11 groups, or 10 with `n_C = 1`?
2. **One `K_d` or two?** Proposal is one, with the Langlands 1997 asymmetry as a
   declared sensitivity check rather than a second parameter. It is a real
   asymmetry but it rests on a rank order measured in yeast against different
   partners, and it costs a parameter out of a budget of eleven.
3. **Arm budget** — still an external unknown. Placeholder 12/24/48 is now written
   into v3 Phase 0, and Phase 4's design must be instantiable at each tier.

---

## 2026-08-15 · Session 5 — the Phase 2 core exists

Guard test written **before** the right-hand side, as instructed. The old model is
deleted in the same commit as its replacement. 38 tests pass.

### The pre-registration discipline paid, and it is worth being able to point at

`prereg/id3_kd_prior_justification.md` was written for **Stage 2 of a plan that no
longer exists**. It stated that nondimensionalization *should* surface
`K_d/E_tot` as a single dimensionless group, and that **if it did not, the
nondimensionalization was wrong.** Two plan rewrites later, on a model with a
different state space, a different mechanism and a different purpose — it does,
as group 12, and it is one of the three profiled parameters.

A structural prediction, recorded in advance, that survived a full architecture
change and came true. That is small, but it is evidence the pre-registration rule
is doing real work rather than being ceremony, and it is the kind of thing that is
much more convincing when you can point at the commit timestamp.

### Two things the right-hand side caught that the parameter table could not

**1 · The two complexes are not one complex.** The first draft wrote *Rbpjl*
production against the PTF1-**L** complex — i.e. RBPJL production requiring RBPJL.
That makes `R = 0` **absorbing**, so no amount of MEK inhibition can restore it and
**the model predicts trametinib alone never reverts anything.** That contradicts
Collins 2014 head-on and would have destroyed Phase 7's trametinib-only positive
control — the arm that makes the one-shot experiment un-failable.

The fix is the biology: *Rbpjl* is driven by PTF1-**J** (with **RBPJ**, broadly
expressed, constant pool absorbed into `a_R` at no parameter cost); the *Ptf1a*
enhancer needs PTF1-**L** (with **RBPJL**). Stated precisely, because the two are
one word apart and only one is defensible:

> The claim is **"nothing but PTF1A makes RBPJL"** — **not** "nothing but RBPJL
> makes RBPJL". The second is stronger, unsupported, **and it passes every guard
> test**.

A missing basal ignition term `b_P` fell out of the same check, for the same
reason. Below `b_P ≈ 0.4` the model is in a regime the literature already rules
out, so `default_params` sits above it — a sanity floor, not a calibration.

**2 · The `n_eff` prefactor is 0.5, not 1.34 — and this one had a price tag.**
The 1.34 was measured on the *deleted* model's ternary complex under *two-target*
titration. The core titrates one target and takes the slope on `E_free`: measured
**0.5**, loose limit 1 rather than 2. **Phase 3 convolves the per-cell LNP dose
distribution against this threshold, so carrying 1.34 over would have inflated the
co-formulation gap — R1, the headline number — by ~2.7×, silently.** The model
would still run; the figure would still render.

It had propagated into five documents (CLAUDE.md, v3, decisions 006 and 012, the
Kd prereg). All corrected, each with the reason rather than a silent edit.
Generalisation now written into 006: **a constant measured on one observable of
one mechanism is not a property of "the model" — when the mechanism changes,
re-measure rather than re-cite.**

### Done

- **`tests/test_bootstrap_guard.py` first, then `src/core.py`.** Seven guard tests:
  the structural claim, an exact-zero check at the origin, **two
  guards-on-the-guard** (a naive `+1e-6` basal term *and* a subtler
  offset-inside-the-Hill version, both of which the guard must catch), the dynamic
  companion (`R` cannot rise from zero with `P` clamped at zero, ERK swept from
  fully suppressed to high), a collapse test, and a positive control that the
  payload can still raise `R` — without which a model where `R` is simply inert
  would pass everything.
- **Budget: 12 free groups**, not the 11 signed off. Two were missed by counting
  on paper: `k_w` (the pERK scale was already spent on ID3, so a second
  ERK-driven process needs its own half-max) and `b_P` (a missing *mechanism*, per
  above). `n_C = 1` on the double-counting argument brings it back inside v3's
  9–12. **Margin is 12 of 12 — any new term must displace an existing one**, and
  `tests/test_core.py` asserts the count so it cannot creep.
- **Qualitative structure confirmed against the golden fixture.** Two attractors,
  metaplastic branch at `R ≈ 0` — the same signature as the pre-rewrite fixture,
  across a complete change of state space. The bistable region is a **bounded
  window**: monostable *acinar* below it (trametinib alone reverts — Collins),
  monostable metaplastic above it. Withdrawal relapses. That is the R3 structure.
- **Deletions, same commit as the replacement:** `model.py`, `topology.py`,
  `functional_forms.py`, `payload_subsets()`, `pin_golden_trajectory.py` (it
  imported the deleted modules; the fixture and its full parameter set are
  committed, and git history holds the script). `binding.py` → `src/supplementary/`,
  with its surviving tests.
- Calls 1–4 written into `docs/PHASE2_PARAMETER_BUDGET.md` and decision 006,
  including the **range-sweep** design for the one-`K_d` check.

### Broke

- Three tests failed on first run, all of them usefully: the parameter count (13
  fields, 12 free), the `n_eff` constant (above), and the bistability test, which
  had been written against low ERK before `b_P` moved the low-ERK regime to
  monostable acinar. All three were wrong assertions, not wrong code.
- **Nothing was tuned to make a result appear.** The payload-rescue question — can
  a payload prevent relapse — was explored, found to depend on drug-hold duration
  relative to `γ`, and **left alone**: that is a Phase 5 result and it needs
  pre-registration before it is swept, not a parameter set chosen until it looks
  right. Recorded here so the restraint is auditable.

### Next

1. **Gate B.** Continuation on `(KRAS × trametinib)`, the persistence window as a
   bounded wedge, saddle and separatrix identified properly rather than by
   settling from two initial conditions.
2. **Profile likelihood on `a_P`, `γ`, `κ`** — and report each as the uncertainty
   on R2, R3, R1 respectively, in those words.
3. **The one-`K_d` range sweep**, Langlands rank order → parity.
4. `fig01_loop_schematic` and `fig02_persistence_window` as soon as Gate B writes
   to `results/`.

### Open questions

1. **`ε` / `α_C` separability** — decided by profile likelihood, not in advance. If
   they lump, report 11 groups and name the combination.
2. **The observation model for `A`** (Phase 0). Three parts agreed: minimal spread
   propagated then thresholded; staining threshold as a **sampled nuisance
   parameter**; Phase 6 targets framed as **timing** rather than absolute
   percentages. **Build it once — Phase 3 needs the same per-cell spread machinery
   for LNP dose heterogeneity.**
3. **Arm budget** — unchanged, external, blocking Phase 4 only.
4. **Bench item 8** (ID3 western ± trametinib) — Luqmaan chasing this week. It
   tests the ERK→ID3 edge, the one assumption Phase 5's ordering prediction rests
   on.

---

## 2026-08-15 · Session 6 — decision protocol; Gate B met; four defects found by panel

First session under the tier protocol. Four panellists ran; three reported. What
they found is more important than what they were asked.

### Tier 1 decisions (decided alone, logged)

- **CLAUDE.md held at ~200 lines, not under it.** The tier protocol is ~30 lines
  of new load-bearing content added after the 200 cap was set; the alternative was
  deleting operating rules. Full protocol in `docs/DECISION_PROTOCOL.md`, operative
  summary in CLAUDE.md. Duplicated reference content (cut-list reasons, limitation
  wording, figure palette detail) now points at v3 and `figures/README.md`.
- **Panel agents run as inlined mandates this session.** `.claude/agents/*.md` were
  written but the registry does not pick up new definitions mid-session; the files
  are correct for future sessions.
- **`refine_fold` added** — fold locations are the persistence-window edges, and
  grid resolution put the fold state out by ~3%.
- **Separatrix bounded to the physical orthant** — below zero `core.rhs` clamps, so
  the integrator was tracing clamped dynamics, not the model.
- **`continue_branch` now Newton-corrects its seed** and raises if it cannot.
- **fig02 renders at `paper` profile** into `figures/out/` from a clean tree.

### Gate B — STRUCTURE MET

  two folds     ERK = 0.0588 and 0.7123 (refined)
  classes       252 stable, 82 saddle
  saddle        2-D stable manifold in 3-D — so "separatrix" is the right noun
  separatrix    48/48 rays, all terminating at the domain edge
  convergence   333 corrector calls, 0 failures

Continuation validated against the cusp normal form, whose folds are known
exactly, plus a guard-on-the-guard that builds the naive sweep and requires it to
miss the middle branch. `fig02_persistence_window` built from `results/` only.

**The other half of Gate B — "key parameters surviving profile likelihood" —
cannot be run at all. See decision 013 and T3 question 1.**

### Tier 2 — DECIDED-PENDING-REVIEW

- **[013 · profile likelihood is not available](decisions/013-profile-likelihood-is-not-available.md)**
  No data, so no likelihood. Measured: the "interval" is 100% of the prior box for
  `a_P` and `γ`, and halving the invented tolerance moved `κ`'s from 100% to 67% —
  the width is a number the modeller typed. Same shape as v3 §0.7 row 1. Replaced
  by the fold locus (exact), constraint-filtered prior-predictive intervals **on
  the deliverables**, and the flat profile reported correctly as *structural
  non-identifiability*. Name reserved for the wet-lab timecourse and pre-registered
  now.
- **[014 · `C` is a strict cascade](decisions/014-chromatin-is-a-strict-cascade.md)**
  `dC/dτ` depends on the input and on `C` — **not on `P`, `R`, or the payload**.
  Verified: identical across nine `(P,R)` combinations and identical at payload 10.
  So the payload has no channel to the durability endpoint and **Phase 5's dose ×
  interval map is flat by construction**. Two fixes; the choice is T3 question 2.

### Broke — four defects, three of them mine

1. **The one-to-one mapping `a_P`→R2, `γ`→R3, `κ`→R1 is false for R1.** `a_R`
   (not profiled) and `n_P` (a scanned exponent) both move the co-formulation gap
   more than `κ` does over `κ`'s own range. **`core.py`'s own `n_eff` docstring
   already said the convolved quantity is `P·E_free·R` — the claim was refuted by a
   file in the same repository.** Withdrawn; Sobol analysis pre-registered.
2. **`b_P` shipped 2% above a saddle-node.** Measured critical value **0.4903**;
   the default was 0.5 and the docs said the boundary was "roughly 0.4" — at 0.40
   trametinib alone leaves `R = 0.18`, i.e. not reverted. `b_P` is *fitted*, so a
   2% nudge would have crossed into a regime the literature forbids with every
   figure still rendering. Default raised to 0.6; `B_P_CRITICAL` exported; two
   regression tests.
3. **`ε = 0.5` makes `C` a filter, not a memory** — bistability needs
   `ε > 3√3/8 = 1.5396`, and `Params.eps` was commented "memory, not filter".
   Corrected; threshold exported as `EPS_MEMORY_THRESHOLD`.
4. **`ID3` is literally the `erk` input**, while the module docstring promised
   "a saturating function of erk". No such function existed. Corrected — and the
   linearity *is* the ERK→ID3 edge, so bench item 8 now has an explicit slot: it
   measures this input map.

### Changed a previously reported number

- **`b_P` critical: "roughly 0.4" → 0.4903 measured.** Default 0.5 → 0.6.
- **The `a_P`/`γ`/`κ` → R2/R3/R1 mapping is withdrawn**, in
  `PHASE2_PARAMETER_BUDGET.md` §1 and everywhere it was to be quoted. It was never
  reported outside the repo.
- **Bench Handshake gains item 9: PTF1A protein half-life.** Time is `τ = t·δ_P`
  and `δ_P` is absorbed by the nondimensionalisation and unmeasured, so **the model
  has no clock**: Collins' two day-valued timings collapse to one dimensionless
  constraint, and **Phase 5's redosing interval has no units**. Item 4 (mRNA
  half-life) does not cover it.

### Next

1. Luqmaan answers the T3 batch (below) — questions 1 and 2 gate Phases 2 and 5.
2. Fill in `prereg/2026-08-15_phase2_sensitivity_prediction.md`, push, **then** run
   the Sobol analysis. Blocked until then by hard rule 1.
3. Re-check whether relapse is chromatin-limited at all before implementing either
   fix in 014 — a panel measured `d ln(relapse)/d ln γ ≈ 0` at high ERK, which
   would make the question moot.
4. `fig01_loop_schematic`; `figS0x` for the fold locus.

### Open questions — batched for Luqmaan, see the session report
