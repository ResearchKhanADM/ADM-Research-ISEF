# Pre-registration · prediction — Phase 2 global sensitivity analysis

**STATUS: TEMPLATE. EMPTY ON PURPOSE. Claude wrote the questions; Luqmaan writes
the answers.**

Hard rule 1: *I write the templates; Luqmaan writes the content — never fabricate
his predictions.* Every answer below is blank because a prediction written by the
model that is about to be tested is not a prediction. **The sweep does not run
until this file is filled in, committed and pushed.** The pushed timestamp is the
only thing that will make the result mean anything in five months.

Ranges: `2026-08-15_phase2_sensitivity_ranges.yaml` — **confirm or edit those
first**, since a prediction against an unconfirmed box is not registered against
anything.

---

## Why this analysis is being run

`docs/PHASE2_PARAMETER_BUDGET.md` §1 claimed the three profiled parameters map
one-to-one onto the three headline results: `a_P` → R2, `γ` → R3, `κ` → R1. A
one-at-a-time diagnostic contradicted the R1 leg — `a_R`, which is not profiled,
and `n_P`, which is a scanned exponent rather than a fitted parameter, both moved
the co-formulation gap more than `κ` did across `κ`'s own declared range.

**The claim is withdrawn.** This analysis replaces it with a measurement: Sobol
first-order and total-order indices for every deliverable against all 12 groups.

The prediction matters more than usual here, because a sensitivity analysis run
*after* a claim has been contradicted is exactly the situation where post-hoc
rationalisation is easiest. Registering the expectation first is what makes the
outcome informative either way.

---

## 1 · Which parameter dominates each deliverable?

For each, name the group you expect to carry the **largest total-order Sobol
index**, and — this is the useful part — **name your second choice**, because
being wrong about first place while right about the top two is a different result
from being wrong about both.

| Deliverable | Expected dominant | Expected second | Confidence (low / med / high) |
|---|---|---|---|
| `persistence_window_width` (R3) | | | |
| `coformulation_gap` (R1) | | | |
| `threshold_dose` (R1) | | | |
| `time_to_relapse` (R3) | | | |

## 2 · Does the withdrawn mapping survive in weakened form?

The strong claim (one-to-one) is dead. The weak claim would be *"`a_P`, `γ` and
`κ` are each the dominant single contributor to R2, R3 and R1 respectively."*

- Do you expect the weak claim to hold? ☐ yes ☐ no ☐ partly — which legs?

- **If it fails, what do you want the project to report instead?** Write this
  *now*, before the numbers exist:

## 3 · Interactions

Total-order minus first-order measures how much of a parameter's influence runs
through interactions with others.

- Which parameter do you expect to have the **largest gap** between total and
  first order — i.e. to matter mostly *through* other parameters?

- Do you expect the deliverables to be dominated by main effects or by
  interactions? A model dominated by interactions is harder to summarise and
  changes how Phase 4's mixture design should be read.

## 4 · The scanned exponents

`n_P` and `n_R` are scanned 1–4 rather than fitted, on the grounds that a Hill
coefficient is a phenomenological summary nobody has measured (decision 004).

- Do you expect `n_P` to rank in the top three for **any** deliverable?

- **If it does, does the scan-don't-fit position still hold?** A parameter that
  dominates a headline result is hard to describe as a nuisance being scanned
  over. State now what you would do:

## 5 · The kill condition

Name a result that would make you **stop and rethink Phase 2** rather than
proceed to Phase 3. Examples of the shape (do not just tick one):

- every deliverable dominated by one parameter, so the model is effectively 1-D;
- no parameter with a total-order index above ~0.1 for R1, so the co-formulation
  gap is not controlled by anything the model contains;
- `κ` ranking last for R1, which would mean the binding treatment — decisions 006,
  the exact `E_free`, the `n_eff` correction — buys nothing for the headline.

Your kill condition:

## 6 · Convergence

A prior diagnostic found solve-failure rates swinging ~13× across `κ` quintiles
and flat across `ρ`.

- Do you expect failures to correlate with `κ` again here? ☐ yes ☐ no

- If they do and the rate exceeds the 0.05 bound in the loose-binding regime,
  what should happen — fix the solver, shrink the box, or report the results as
  conditional on convergence and say so on the figure?

---

## Sign-off

- [ ] Ranges reviewed and confirmed (or edited) in the `_ranges.yaml`
- [ ] Predictions above written by Luqmaan
- [ ] Committed **and pushed** before the sweep runs

Signed: ______________________  Date: ______________
