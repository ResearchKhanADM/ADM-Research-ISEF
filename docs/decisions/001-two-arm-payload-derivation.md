# 001 · Two arms for payload derivation — mechanistic model plus unbiased screen

*Date:* 2026-08-15
*Status:* accepted

## Question

The intended payload is RBPJL + PTF1A + NR5A2 (+ trametinib). Three of those are
state variables in the ODE. If the mechanistic model is then asked *"which
molecules should be supplied?"*, it can only nominate species that were written
into it — the answer is determined by the model's construction, not by the
biology.

Should the project (a) keep a single mechanistic arm and answer the charge by
restating what the model actually claims, or (b) add a second, data-driven arm
that starts from a candidate set the modeller did not curate?

This is architectural, not cosmetic. It changes what the model is *for*, adds a
stage, adds a dataset dependency, and changes what goes on the board.

## Positions considered

**Position A — one arm. The circularity charge is answerable, not structural.**

The model was never asked to discover molecules. It is asked whether a
*specified* payload can push a cell across a separatrix, at what dose, on what
schedule, in what order, and at what cost in viability. Those are exactly the
questions a mechanistic model exists to answer, and none of them are circular.
Molecule identity came from the literature — RBPJL from Masui 2010 and the
no-independent-promoter asymmetry, E-protein from Dufresne 2010, PTF1A from Krah
2019 — and citing a source is not the same as assuming a conclusion.

Adding a second arm is not free. It introduces a second toolchain, a second set
of failure modes (GRN inference quality, motif-database coverage, species
mismatch), 2–3 weeks against a schedule that already has a wet lab in the same
cycle, and the possibility of manufacturing a contradiction with the first arm
that then has to be explained. The cheapest correct fix is to state plainly
which claims are literature-derived and which are model-derived.

**Position B — two arms. Only an arm that can return a bad answer is evidence.**

"The literature chose the molecules" does not survive contact with a hostile
question, because the literature selection and the model structure were produced
by the same person reading the same papers. The state variables *are* the payload
species. That loop is closed regardless of how carefully it is described.

Falsifiability requires that the method could have produced an unwelcome result.
An unbiased screen over the full TF repertoire present in a pancreatic dataset
can rank PTF1A, RBPJL and NR5A2 *low*. That possibility is the entire value. If
the screen independently nominates the same factors, the convergence is worth far
more than either arm alone — two methods with different failure modes agreeing.
If it nominates something else, that is a better result still, and it arrives
early enough to act on.

The cost is bounded because the task shape has a published precedent: CellOracle
was used to rescue a *failing* lineage conversion by in-silico perturbation,
nominating Fos and Yap1 (Kamimoto 2022, PMID 36584685) — factors that were not
the field's prior expectation. That is the same shape as this problem.

*Where they actually disagree:* both positions agree the mechanistic model
legitimately answers dose, schedule, ordering and viability. Neither disputes
that. They disagree on one point only — **whether payload *identity* can be
defended by a method that was structurally incapable of nominating anything
else.**

## Decision

**Two arms.**

1. **Mechanistic arm — reframed.** It is no longer asked "which molecules." It
   is asked three things it can answer without circularity:
   - **Necessity.** Four interventions → 16 subsets. Run all 16 across the
     parameter ensemble. Report which components are necessary, which are
     redundant, and each one's marginal contribution. *"Three mRNAs is
     over-engineering"* is a valid and useful outcome.
   - **Dose, schedule, order** — Stages 5 and 6, unchanged.
   - **Whether the reversal-optimal payload and the viability-optimal payload
     are the same payload.** If they differ, that is a headline result and it is
     precisely the stated objective, which is reversal *and* viability rather
     than either alone.
2. **Screen arm — new Stage 3B.** CellOracle (Kamimoto 2023, PMID 36755098) on
   pancreatic scRNA-seq containing both acinar and metaplastic populations.
   Simulate overexpression across the full TF repertoire, score the shift of the
   metaplastic cluster toward acinar, return a **ranked list**, pre-registered.

The reason that was decisive: **the necessity and Pareto questions are things the
mechanistic model can answer and the screen cannot, and identity is a thing the
screen can answer and the model cannot.** The two arms are not redundant and not
in competition — they answer disjoint questions. Position A was right that the
model is not circular *for its own questions*; it was wrong that identity is one
of them.

Numbered **3B** rather than renumbering Stages 4–7, because those numbers are
referenced from `CLAUDE.md`, the progress log, and prereg filenames, and a
renumbering five months long is a source of silent confusion for no benefit.

## What would reverse this

Concrete, and checked **before** the screen is trusted — these are pre-flight,
not post-hoc:

1. **RBPJL has no motif in the base GRN.** CellOracle propagates a perturbation
   through edges derived from motif scanning. A TF with no motif has no outgoing
   edges and its simulated perturbation returns approximately zero — *silently*,
   looking like "RBPJL doesn't matter" rather than "RBPJL is not representable."
   RBPJL is obscure enough that this is a live risk, and it is the single
   component the whole mechanistic argument rests on. **If RBPJL cannot be
   perturbed, the screen cannot speak to the most important component**, its
   answer to the circularity charge collapses, and it degrades to ranking
   PTF1A vs NR5A2 vs others. Revert to Position A plus honest framing, and say so
   on the board.
2. **Rank instability across datasets.** Run the screen on GSE207938 and repeat
   on GSE172380. If top-20 overlap is below ~50%, the ranking is measuring
   inference noise, not biology, and no result may be reported from it.
3. **The base GRN is generic rather than pancreatic.** If the dataset has no
   matched scATAC and the fallback base GRN yields a network with no
   PTF1A→acinar-target edges, the screen is measuring the motif database rather
   than the data. Check for those edges before running any simulation.

Note on process: `CLAUDE.md` rule 4 asks for a 2-subagent adversarial panel on
architectural decisions. This session's tool policy does not permit spawning
subagents, so both positions above were argued directly rather than by panel.
The reversal conditions are unaffected, but the record should show how the file
was produced.
