# Decisions — index

**Read this at session start. Open a full file only when it bears on the task in
hand; never load the folder.** Kept current in the same commit as any new or
retired decision — an index that drifts is worse than none.

*Status:* **live** = in force · **amended** = superseded in part, read the
amendment at the top · **retired** = the claim is withdrawn, not paused ·
**moot** = correct, but the question it posed does not arise.

| # | Title | One sentence | Status |
|---|---|---|---|
| 000 | TEMPLATE | Four headings: question, positions, decision, what would reverse this. | live |
| 001 | Two-arm payload derivation | Payload identity cannot come from a model whose states *are* the payload; a second data-driven arm was added. | **amended** — identity is published (Jiang 2023); the screen is demoted to a declared negative control (012) |
| 002 | `W` protected from elimination | Phospho-MEK added as a state so the trametinib-vs-PD325901 withdrawal asymmetry could exist. | **amended** — `W` removed, pERK is an input with a swept **rebound profile**; the drug-identity prediction is **retired outright** |
| 003 | Composable topology architecture | One right-hand side plus config flags, so five topologies could be compared under identical code. | **retired** — the topology competition is cut (012) |
| 004 | Hill functional forms | Separable activation × repression; exponents **sampled 1–4, never fixed**. | live |
| 005 | ERK / drug layer | `g(K_eff)` takes one argument so the drug cannot be double-counted. | **retired** with `W` (002) |
| 006 | Two-target competitive titration | ID3 titrates **both** E47 and PTF1A; solved at equilibrium, not integrated. | **amended** — promoted: it is now the justification for the core's exact `E_free`, but its `n_eff` prefactor **1.34 does not transfer** (the core's is 0.5) |
| 007 | `dI/dt` and `dA/dt` | ID3 and acinar-output equations, plus the warning that `A` is a concentration while Collins reports a **fraction of cells**. | **amended** — `A` is no longer a state; the observable warning **still stands** and is a Phase 0 item |
| 008 | U-shaped viability hazard | Death as an integrated hazard rising at both ends, never a threshold. | **retired** — viability is a bench-measured floor; only a CHOP output **flag** survives |
| 009 | *(merged into 008)* | — | retired |
| 010 | NR5A2 placement | Enhancer vs acinar-output co-activator, tested as competing topologies. | **retired** with the topology competition (012) |
| 011 | mRNA pulse forcing | Analytic difference-of-exponentials, normalised so each pulse delivers exactly `dose`. | **live** — the normalisation matters *more* now (Phase 4 compares at matched total mass); only its Lie-bracket section is dead |
| 012 | Durability framing — the architecture change | Reversal was solved in 2014 and relapses; the project is rebuilt around durability, delivery and composition. | **live** — the plan of record |
| 013 | Profile likelihood is not available | No data means no likelihood; the "interval" was the prior box, and its width was a modeller-chosen tolerance. | **live** |
| 014 | `C` is a strict cascade | `dC/dτ` sees neither the states nor the payload, so the payload has no channel to durability. | **moot** — the cascade is real, but relapse is not chromatin-limited (015) |
| 015 | Relapse is not chromatin-limited | Measured three ways; **R3 becomes a threshold claim** in post-withdrawal drive, flat in dose and in drug-hold. | **live** — conditional on `ε < 1.54`; see its headline caveat |

## Standing reversal conditions worth knowing without opening a file

- **015 inverts if `ε > 1.5396`** — a bistable `C` could hold its own state after
  withdrawal and re-suppress the acinar program. The pre-registered sweep
  straddles the boundary. **This is the single finding that would reverse the
  current headline conclusion.**
- **012 reopens** if a measurement pins the MEKi-vs-forced-PTF1A kinetic ratio;
  the topology competition becomes runnable.
- **013 resolves** the day a quantitative timecourse exists — profile likelihood
  is pre-registered for it, under its own name.
