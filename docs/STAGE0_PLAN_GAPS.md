# Stage 0 · complete gap audit of master plan §3.1 and §3.2

*Date:* 2026-08-15
*Purpose:* one exhaustive pass, so all plan changes can be absorbed at once
rather than one per session.

**Method.** §3.2 writes six equations: `dP_n/dt`, `dR/dt`, the E-titration
relation, `dM/dt`, the trametinib/`W` block, and `dS/dt`. §3.1 lists **eleven**
states. Every state without an equation is listed below, together with every
place where §3.2 contradicts itself or contradicts a §3.4 constraint.

**Categories.**
**[A] — I resolve it.** §3.4 assigns functional forms to Stage 0. Each gets a
`docs/decisions/` file. No sign-off needed.
**[B] — needs Luqmaan's decision.** Changes the state list, contradicts a
declared constraint, or changes what the model claims.

---

## [B] — needs a decision before Step 1 is written

### B1 · §3.2's `dP_n/dt` violates constraint 2, in its own section

`dP_n/dt` contains **`− k_seq·I·P_n`** — a first-order sink. Two paragraphs
later the same section says: *"A first-order sink `−k·I·P_n` generates no
ultrasensitivity"*, and §3.4 lists "ID3 acts by titration, not a first-order
sink" as non-negotiable.

The equation and the constraint are in direct contradiction. §3.4 declares the
constraint authoritative, so **the term is deleted** and ID3 enters only through
the binding polynomial. Flagging rather than silently fixing, because the
contradiction is in the plan and will be read by others.

### B2 · ID3 titrates **two** targets, not one — the polynomial is competitive

§3.2's relation names only `(ID3·E complexes)`. But §1.3 and §2.3 both say ID3
traps **E47 *and* PTF1A**, and this is confirmed: Dufresne 2010 (*Int J Cancer*
129(2):295–306, PMID 20830706, doi 10.1002/ijc.25668) reports gastrin raising Id3
and increasing **both** Id3/E47 and Id3/Ptf1-p48 interactions while decreasing
E47/Ptf1-p48 — in AR4-2J, the wet-lab line.

So the binding polynomial is a **two-target competitive titration**: ID3, E and
PTF1A compete across ID3·E, ID3·PTF1A, and E·PTF1A(·RBPJL/RBPJ). This raises the
polynomial's degree and is materially more algebra in Step 2 — and it is also
*good* for the model, since two-target titration is a stronger ultrasensitivity
generator than one. Confirming this is the intended reading.

### B3 · Total vs free E-protein are conflated — a state is missing

§3.1 lists `E` as "free E-protein, fast → eliminate". Two different quantities
are hiding here:
- **free E** — set by binding equilibrium, genuinely fast, correctly eliminated;
- **total E** — set by synthesis and degradation, genuinely **slow**, and *the
  thing `u₃ = E47` adds to*.

You cannot eliminate the pool the payload doses. **Proposing `E_tot` as a slow
state, with `E_free` recovered algebraically from the polynomial.** This is a
state-list change, hence [B].

### B4 · `γ(MIST1-driven capacity)` — MIST1 has no state and no equation

`dS/dt` depends on a MIST1-driven capacity that exists nowhere in §3.1. MIST1 is
also a u₃ candidate, so it must be addressable by an input.

**Proposing `γ = γ(A)`** — capacity co-induced with the acinar program, which is
what Jakubison 2018 supports — rather than adding a MIST1 state. Cheaper, and it
preserves the §3.2 claim that capacity and cargo are co-induced. **The cost: if
MIST1 is only a function of `A`, then "MIST1 as u₃" cannot be dosed
independently of the acinar program, which weakens Stage 3's ability to rank it.**
That trade is yours to make.

### B5 · NR5A2 is not represented anywhere, but Stage 3 must rank it

Stage 3 ranks u₃ candidates E47/TCF3 vs NR5A2 vs MIST1. E47 enters via `E_tot`
(B3); MIST1 via capacity (B4); **NR5A2 enters nowhere.** It has no state, no
term, and no route into any equation.

A ranking exercise cannot rank a candidate that has no representation — it will
score zero by construction, exactly the silent-null failure mode flagged for
RBPJL in Stage 3B. Needs either a representation (NR5A2 as a co-activator on the
acinar enhancer is the usual reading) or explicit removal from the candidate
list. **This is also the payload you said you ideally want to test**, so it
should not be the one that is unrepresentable.

### B6 · The viability axis is described but not defined — and it is self-contradictory

§3.2 says both:
- *"death when `S` exceeds tolerance"* — a **one-sided** threshold on `S`;
- *"the curve is **U-shaped** — too little PTF1A also kills, via CHOP-dependent
  apoptosis"*.

Those are different objects, and a one-sided threshold on `S` is close to the
`U_crit` construct the panel already killed. Half the headline claim rides on
this axis, so it needs an explicit definition: a death hazard depending on
**both** the transient cargo/capacity mismatch **and** low differentiation.
Proposing a hazard function rather than an absorbing threshold, so trajectories
stay integrable and "viability" is a reported quantity rather than a stop
condition. Needs your call because it defines the y-axis of the headline figure.

### B7 · mRNA inputs have no shape — Stage 5's axis depends on it

§1.3 is emphatic that the inputs are **pulses** — uptake, translation ramp,
first-order decay — and that this is what makes duration a real quantity rather
than a free variable. §3.2 writes them as bare additive `u₂(t)`, `u₁(t)`.

**Proposing analytic pulse forcing** (sum of decaying exponentials with a
translation ramp), not three extra mRNA states — same shape, no state inflation,
and it keeps the reduced system small. Flagging because Stage 5 sweeps *(dose per
pulse) × (redosing interval)* and that axis is undefined until this is fixed.

---

## [A] — I resolve these, each with a decision file

| # | Gap | Resolution |
|---|---|---|
| A1 | **`dI/dt` missing** (ID3). Timescale "intermediate". | ERK-driven production, first-order decay. Timescale matters: if `I` is fast it slaves to `K_eff` and titration loses its own dynamics. Sampled, not fixed. |
| A2 | **`dA/dt` missing** (acinar output). | Driven by `C_L`, first-order decay. Also needs an explicit map from `A` to *"% amylase-positive cells"* — Collins' numbers are percentages of cells, not concentrations, and validation is meaningless without the map. |
| A3 | **`dP_c/dt` missing.** | Follows from B1/B2: `P_c` is the ID3·PTF1A complex, set by the polynomial, not by a separate sink term. |
| A4 | **`RBPJ` has no state or parameter**, yet `C_J = PTF1A + E + RBPJ`. | Constant pool parameter — RBPJ is broadly expressed and not the bottleneck. Declared explicitly rather than assumed. |
| A5 | **PTF1A mass balance ambiguous** — is `P_n` free or total nuclear? | Free. Complex formation debits it. Conservation stated explicitly in code. |
| A6 | **`C_L`/`C_J` are ternary complexes**; assembly is not one bimolecular step. | Sequential assembly with a stated intermediate; affects the polynomial's degree. |
| A7 | **`K` has no map from KRAS/dox dose.** | Explicit dose→`K` map. Stage 1 sweeps this axis, so it needs units the bench can set. |
| A8 | **`Hill(C_L, M; n)` — `M` sits at metaplasia loci but gates the acinar enhancer.** | Constraint 3 puts `M` at metaplasia loci, so it cannot directly gate the 2.3-kb acinar enhancer. `M` acts **indirectly**: the metaplasia program represses the acinar program. Same effect, correct location. |
| A9 | **`φ(K, C_L)` in `dM/dt` undefined; `g(K,v)` now redundant with `K_eff`.** | With `W` added, `g` takes `K_eff` only. Removes a double-count of the drug. |
| A10 | **Death: state, absorbing condition, or post-hoc hazard?** | Post-hoc hazard (see B6), so trajectories integrate to completion. |

---

## Not a gap, but a code requirement

The necessity analysis (§1.3) needs **all four interventions independently
toggleable**, and `u₃`'s identity **swappable**, from the first line of code. If
that is retrofitted later it will be retrofitted badly. Building it in now.

---

## Status

Nothing in `src/` yet. Step 1 is not written and should not be written until
B1–B7 are settled, because five of the seven change either the state list or an
equation that everything downstream inherits.
