# 003 · Composable topologies — one right-hand side, optional terms behind flags

*Date:* 2026-08-15
*Status:* accepted

## Question

Resolving B1–B7 grew the state list from 10 to roughly 13 (`W`, `E_tot`,
`MIST1`, `NR5A2`). Stage 0 targets **5–6 slow states**, and continuation on a
large system is the single biggest schedule risk in the project.

How is the model structured so that the state list can grow to cover five
competing topologies without the continuation problem becoming unrunnable?

## Positions considered

**Position A — one maximal model containing every state.**

Simplest to reason about: there is one right-hand side, one parameter table, one
integrator path, and no configuration logic that can silently disagree with
itself. Topologies are then expressed by setting the parameters of unused terms
to zero.

It is also the honest default. A model that changes shape depending on a flag is
harder to describe in a writeup, and "which version produced this figure?" becomes
a real question that has to be answered every time.

**Position B — five separate model files, one per topology.**

Maximum clarity per file. Anyone can read `topology_T3.py` top to bottom and see
the whole system without tracing conditionals. No risk that a flag combination
produces a system nobody intended.

**Position C — a CORE right-hand side plus optional terms selected by config.**

Topologies are configuration objects. The sampler, integrator, scoring code and
parameter machinery are shared by construction.

*Where they actually disagree:* all three produce the same trajectories for any
single topology. They differ entirely in what happens **across** topologies and
over five months of edits.

## Decision

**Position C.**

Two reasons, and the first is fatal to the alternatives.

**1. Stage 2 is a controlled comparison, and A and B both break the control.**
The Q-value methodology (Ma, Trusina, El-Samad, Lim & Tang, *Cell* 2009,
PMID 19703401) compares the fraction of an identical sampling box that achieves
the target function, across topologies, **with the same code**. Its entire
validity rests on nothing differing except the right-hand side.

- **Position B fails outright.** Five files drift. A tolerance tightened in three
  of them, a sampler bug fixed in one — and the Q-value comparison silently stops
  being a comparison. The failure is invisible: every file still runs.
- **Position A fails more subtly.** Zeroing a term is not the same as the term
  being absent. A zeroed state still enters the Jacobian, still contributes a
  (zero) eigenvalue, still costs continuation dimensions, and still appears in
  the FIM eigenspectrum in Stage 3 as a spurious sloppy direction. **The model
  selection would be run at a dimension no candidate topology actually has.**

**2. It is the only option that keeps the slow count at 6.** CORE is `P_n`, `R`,
`E_tot`, `I`, `M`, `A`, `S`, `W`. At reduction, `P_c`/`E_free`/`C_L`/`C_J` are
already algebraic, `I` is the QSS candidate, and `W` is retained by exemption
(002) — leaving **6 slow states plus `W`**. `MIST1` and `NR5A2` are simply not
present when they are not being tested, so Stage 1's continuation never sees
them. Under Position A they would be present always, and the schedule risk that
motivated the whole reduction returns.

**Position A's legitimate objection is answered, not dismissed.** *"Which version
produced this figure?"* is a real question. Every result therefore records its
topology configuration alongside the output, and the config is part of the run
identifier — not a note in a log.

**Binding consequences for later stages:**
- Stage 2 swaps topology with **one argument**.
- Stage 3 runs sensitivity over **whichever states are active**, not a fixed list.
- Stage 5's necessity analysis toggles all four interventions independently, and
  `u₃`'s identity is swappable — built in from the first line, not retrofitted.
- T6a/T6b (NR5A2 placement) exist *because* composition makes an extra topology
  nearly free. Under Position B, testing an assumption would have cost a sixth
  file, and it simply would not have been tested.

## What would reverse this

1. **A flag combination produces a system nobody intended, and it reaches a
   reported result.** If configuration logic ever silently yields an unintended
   right-hand side, the clarity argument for Position B wins. Mitigation, and the
   thing to check first: every topology validates its own state set and term list
   at construction, and refuses to run an undeclared combination.
2. **Continuation still fails at 6 slow states + `W`.** If Stage 1 cannot be made
   to converge at that dimension, composition was not the binding constraint and
   the reduction itself needs revisiting — not the code architecture.
3. **The config layer costs more time than five files would have.** If more than
   roughly a week goes into configuration plumbing before Stage 1 produces a
   figure, that is a signal the abstraction was built too early and too generally.
   Collapse to the topologies actually needed.
