# READ FIRST — live state
Updated 2026-08-15 · commit `1ab0ed3`+ · session 9

Read this before `CLAUDE.md`. It exists so someone with **zero context** can
resume in one read. Then: `CLAUDE.md` → last `docs/PROGRESS.md` entry →
`docs/decisions/INDEX.md`.

## WHERE I AM

**Phase 2 complete; Gate B structurally met.** The 3-state core (`src/core.py`,
12 free groups) is built, tested, and independently cross-checked two ways. The
persistence window, the saddle and the basin boundary are computed and figured
(`fig02`). **The pre-registered Sobol sensitivity sweep is running detached right
now** — it is the thing that settles which parameters actually drive each
deliverable, after the previously-asserted one-to-one mapping was withdrawn.
Nothing is blocked on me; two things are blocked on Luqmaan (below).

## RUNNING DETACHED

| Job | Launched | Expected | Output | Quarantines |
|---|---|---|---|---|
| *(none — Sobol completed)* | — | — | `results/sobol/` | — |

**Sobol completed. P2 did NOT fire** (`S1(γ) = 0.0000`), **P6 holds**, failure rate
**0.0000%** overall and in every quintile of every parameter. Scored: P2 HIT,
P6 HIT, **P1 MISS**, **P4 partial**, **P3 partial — the step is real**.

## LIVE NUMBERS

| Quantity | Value | Units | Source | Status |
|---|---|---|---|---|
| Persistence window | `[0.14007, 0.71786]` | ERK (ID3 units) | `results/gate_b/summary.json` | **LIVE** |
| Window width | `0.5778` | ERK | same | **LIVE** |
| Same, by independent method | `[0.1401647, 0.7178549]` | ERK | `src/equilibria.py` scalar reduction | **LIVE** — agrees to 2e-7 (upper), 1e-4 (lower) |
| `B_P_CRITICAL` | `0.4903` | — | `src/core.py`, bisection | **LIVE** — below it, trametinib alone cannot revert |
| `EPS_MEMORY_THRESHOLD` | `1.5396` = `8/(3√3)` | — | `src/core.py`, exact | **LIVE** |
| `n_eff` prefactor | `0.5` (not 1.34) | — | `src/core.py`, measured | **LIVE** |
| Parameter budget | 12 free groups (10 fitted + 2 scanned) | — | `docs/PHASE2_PARAMETER_BUDGET.md` | **LIVE** |
| Separatrix angular coverage | largest gap `185°` vs `7.5°` uniform | deg | `results/gate_b/summary.json` | **LIVE** — means the 48 rays **collapsed**; what is computed is the *in-plane* 1-D separatrix, not the 2-D `W^s` |
| Fold loci vs `a_P`, `b_P` | see `results/fold_loci/` | ERK | rerun completed post-fix | **LIVE** |
| Sobol ST, window width | `κ` 0.663 · `a_P` 0.365 · `c_rep` 0.222 · `α_C` 0.215 | — | `results/sobol/summary.json` | **LIVE** |
| Sobol S1(`γ`), S1(`ρ`) | `0.0000` exactly | — | same | **LIVE** — a *theorem*, not a measurement: both multiply a whole RHS row, so neither moves an equilibrium |
| `C` multistable in the box | `1.73%` (954/55,296) | — | `results/sobol/p3_verdict.json` | **LIVE** — **R3's caveat STANDS** |
| `ε` in every multistable sample | `> 1.7919` | — | same | **LIVE** — all above `EPS_MEMORY_THRESHOLD`, **zero violations in 55,296** |
| Window width when `C` multistable | `4.57` vs `1.44` | ERK | same | **LIVE** — 3.2× wider |

**All times are in `1/δ_P`.** There is no clock until bench item 9 lands. No
figure gets an hours axis; no placeholder conversion, ever, including drafts.

## CHANGED OR RETRACTED — append-only

1. **`n_eff` prefactor 1.34 → 0.5.** The 1.34 was measured on the *deleted*
   two-target ternary complex; this core titrates one target. Reusing it would
   have inflated R1 ~2.7×. Corrected in 5 documents. *(session 6)*
2. **`b_P` critical "roughly 0.4" → 0.4903 measured**; default 0.5 → 0.6. The old
   default was **2% clear of a saddle-node** on a *fitted* parameter. *(session 6)*
3. **The `a_P`→R2 / `γ`→R3 / `κ`→R1 one-to-one mapping — WITHDRAWN.** `a_R` (not
   profiled) and `n_P` (a scanned exponent) move R1 more than `κ` does over its
   own range. Refuted by `core.py`'s own `n_eff` docstring, in-repo, before it was
   written down. *(session 6–7)*
4. **"Profile likelihood" — WITHDRAWN as a misnomer.** No data ⟹ no likelihood.
   The "interval" was 100% of the prior box, and its width was a modeller-chosen
   tolerance. Logged in v3 §0.7 as row 4. *(decision 013)*
5. **R3 restated as a threshold claim.** Relapse is not chromatin-limited, so
   "durability is set by drug-hold duration" was **also** wrong. *(decision 015)*
6. **`fig05_durability` STRUCK** — not producible; relapse timing is flat in both
   dose and drug-hold. Replacement must show the threshold in post-withdrawal
   drive, not a graded dose response. *(session 7)*
7. **Persistence window `[0.0588, 0.7123]` → `[0.1401, 0.7179]`.** Cause: a
   finite-difference step crossing the `max(y,0)` clamp returned **exactly half**
   the true Jacobian entry on the metaplastic branch. Signs survived, so Gate B's
   structure was never wrong; eigenvalues and timescales off that branch were 2×
   out. *(session 8)*
8. **Separatrix relabelled** — the 48-ray output is the **in-plane** separatrix at
   `C = C*` (1-D), not the 2-D `W^s`. *(session 8)*
9. **R2 merged into R1** — see below. *(session 9)*

**No document still quotes a superseded value**; the repo was swept after each.

## THE RESULTS — current wording

**R1 (formulation + composition, merged) — reversal is threshold-limited, not
ratio-limited.** What matters is getting both components above threshold **in the
same cell**, which is why co-formulation gives `p` rather than `p²`. Above
threshold, **composition is free**: the ratio does not move the attractor, and the
model is explicit that it cannot.
*Evidence:* every PTF1A:RBPJL mass split from 2% to 98%, at totals 2/10/40, lands
within **1% of maximum**, with the final attractor identical to 4 dp. And `p ≥ p²`
always, so the co-formulation direction cannot invert — a panel claim that it
could was **checked and refuted analytically**.
*What would reverse it:* a bench result showing durable reversal depends on the
mix ratio at fixed total mass.

**R2 (persistence window) — durability is a threshold property of the
post-withdrawal KRAS drive**, not a graded property of dose or schedule. The
payload buys the crossing; whether it sticks is decided by whether the drive sits
inside the window, and nothing the payload or schedule does moves that window.
*What would reverse it:* **`ε > 1.5396`.** If the plausible range for `ε` sits
above that boundary, `C` becomes a bistable memory that can hold its own state
after withdrawal and re-suppress the acinar program, **and R2 inverts.** The
running sweep straddles the boundary. *This is the single named, numbered thing
that would make us wrong — say it unprompted.*

## DECISIONS INDEX

Full files in `docs/decisions/`; one-line index in `docs/decisions/INDEX.md`
(read that, not the folder). Live and load-bearing: **012** durability framing ·
**013** profile likelihood unavailable · **015** relapse not chromatin-limited ·
**004** Hill forms · **006** two-target titration (amended) · **011** pulse
forcing. Retired or moot: **002, 003, 005, 008, 010, 014**.

## BLOCKED ON HUMAN

1. **Arm budget** (12 / 24 / 48 wells?) — blocks Phase 4's mixture design.
   *Workaround in progress:* the design is being built instantiable at any tier,
   with what changes at each stated.
2. **Bench item 8 — ID3 western ± trametinib** — tests the ERK→ID3 edge, the one
   assumption Phase 5's ordering prediction rests on. Cheap, one western.
   *No workaround exists;* if trametinib does not lower ID3, the ordering arms
   should not be run.
3. **Bench item 9 — PTF1A protein half-life** (cycloheximide chase) — **gives
   every timing result units.** Until it lands, all times stay in `1/δ_P`.
   *No workaround exists.* Equal priority to item 8.
4. **Future pre-registration content** — predictions are Luqmaan's beliefs and are
   never generated here. When a sweep needs one, the template is written and this
   list gains an entry; downstream work that does not depend on the outcome
   continues meanwhile.

## REVIEW QUEUE — decided without asking, most consequential first

| # | Decision | Read |
|---|---|---|
| 1 | **R2 merged into R1**; three results become two. Ratio is not quotable, so it is not quoted. | this file + v3 Part 1.3 |
| 2 | **"Profile likelihood" withdrawn**, replaced by fold locus + prior-predictive intervals on deliverables. | `decisions/013` |
| 3 | **R3/R2 restated as a threshold claim**; contradicted a framing Luqmaan proposed. | `decisions/015` |
| 4 | **`fig05_durability` struck.** | `decisions/015` downstream |
| 5 | **P5/P6 conflict in the prereg resolved** by a stratified design; P5's operationalisation changed. | `prereg/2026-08-15_..._prediction.md` |
| 6 | Jacobian clamp bug fixed; window number changed. | `PROGRESS` session 8 |
| 7 | Separatrix relabelled to the in-plane 1-D manifold. | `PROGRESS` session 8 |
| 8 | `b_P` default raised 0.5 → 0.6 for saddle-node margin. | `src/core.py` |

Nothing here is flagged **CONSUMES-THE-SHOT**. When the wet-lab protocol, arm
allocation or go/no-go thresholds are drafted, they will be, and they will not be
treated as settled.

## CUT LIST — do not re-propose

Five-way topology competition · 13-state ODE · FIM sloppiness · Lie-bracket
ordering formalism · 2^k subset enumeration · submodular/greedy optimization ·
Pareto front as deliverable · U-shaped viability hazard as an ODE ·
trametinib-vs-PD325901 asymmetry · Gillespie bimodality vs scRNA-seq · CellOracle
as validation (survives only as a declared negative control) · Enformer/Borzoi ·
AlphaFold→k_on · Perturb-seq · one-sided U_crit · structural/Kalman
controllability · full-dim HJ reachability · MPC · all-atom MD.
**Each has a specific reason in v3 Part 4 — they are not interchangeable.**
