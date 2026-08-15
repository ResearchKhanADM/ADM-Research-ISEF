"""Two-target competitive titration — the term that decides bistability.

ID3 traps BOTH partners of an obligate heterodimer. Verified: Dufresne 2010
(Int J Cancer 129(2):295-306, PMID 20830706) reports gastrin raising Id3 and
increasing both Id3/E47 AND Id3/Ptf1-p48 interactions while decreasing
E47/Ptf1-p48 — in AR4-2J, which is the wet-lab line. Silencing Id3 reversed the
cytoplasmic mislocalisation.

WHY THIS MATTERS MORE THAN ANY OTHER TERM. A first-order sink (-k*I*P) generates
no ultrasensitivity at all. Molecular titration of a stoichiometric partner is a
classical ultrasensitivity generator and can produce switching without large
Hill coefficients. A titrator that sequesters BOTH members of an obligate
heterodimer is sharper still, because both routes to complex formation shut
simultaneously. If the model turns out not to be bistable anywhere in plausible
parameter space, this is the first place to look.

This module solves the equilibrium numerically. The explicit algebra, every step
shown, is Step 2 -> docs/derivations/binding_polynomial.md. Code and derivation
must agree; the test suite checks them against each other.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares


def _residuals(log_free, p_tot, e_tot, i_tot, r_tot, kd, rbpj):
    """Conservation residuals in log space.

    Log-parameterised because every unknown is a concentration and must stay
    strictly positive. A linear-space solver will happily step negative, return
    a mathematically valid root that is physically meaningless, and the ODE will
    then integrate nonsense without raising anything.
    """
    p, e, i, r = np.exp(log_free)

    # Pairwise ID3 complexes — the two-target trap.
    ie = i * e / kd["IE"]
    ip = i * p / kd["IP"]

    # Ternary PTF1 complexes. Single-step convention: assembly intermediates are
    # far faster than anything else here, so tracking them buys nothing. The
    # sequential-assembly alternative is noted in decision 006 as the reversal
    # condition if the effective Kd turns out to depend on assembly order.
    c_l = p * e * r / kd["L"]          # PTF1-L: PTF1A + E + RBPJL  (adult)
    c_j = p * e * rbpj / kd["J"]       # PTF1-J: PTF1A + E + RBPJ   (immature)

    # Conservation. Each state variable is a TOTAL pool; free species are what
    # we solve for. This is the structural fix from B1/B3: sequestration lives
    # in the equilibrium, not in a sink term bolted onto dP_n/dt.
    #
    # Residuals are scaled by their own total. Unscaled, a pool at 1e-3 and a
    # pool at 1e3 contribute incomparably to the cost and the solver optimises
    # the large one while leaving the small one badly wrong — which in this
    # model is usually the free-PTF1A pool that everything downstream reads.
    return np.array([
        (p + ip + c_l + c_j - p_tot) / max(p_tot, 1e-30),
        (e + ie + c_l + c_j - e_tot) / max(e_tot, 1e-30),
        (i + ie + ip - i_tot) / max(i_tot, 1e-30),
        (r + c_l - r_tot) / max(r_tot, 1e-30),
    ])


def solve_binding(p_tot, e_tot, i_tot, r_tot, kd, rbpj, guess=None):
    """Return dict of free species and complexes at binding equilibrium.

    Quasi-equilibrium is justified by timescale separation: protein-protein
    association here is seconds-to-minutes against transcriptional dynamics of
    hours-to-weeks — three or more orders of magnitude. This is the ONE place
    the model uses equilibrium rather than explicit kinetics, and it is safe for
    exactly the reason the fast-acetylation variable was not: eliminating H
    collapsed Hill(C_L,H) into Hill(C_L,C_L) and the variable vanished
    algebraically, whereas here the eliminated species retain a nontrivial
    nonlinear dependence on the slow states. Elimination is only dangerous when
    it destroys structure.
    """
    totals = np.array([p_tot, e_tot, i_tot, r_tot])
    args = (p_tot, e_tot, i_tot, r_tot, kd, rbpj)

    # Multiple starts, tried in order. A single naive guess converges fine at
    # loose binding and FAILS in the tight-binding regime (Kd << totals) — which
    # is precisely the regime where titration generates ultrasensitivity. A
    # solver that only works where the mechanism is absent would let Stage 2
    # conclude "titration isn't sharp" when it had simply never sampled there.
    starts = []
    if guess is not None:
        starts.append(np.asarray(guess, dtype=float))
    starts.append(np.log(np.maximum(totals * 0.5, 1e-12)))
    # Tight-binding limit: the titrator consumes its targets almost completely,
    # so free species sit far below their totals. Seed several decades down.
    for depth in (1e-2, 1e-4, 1e-6):
        starts.append(np.log(np.maximum(totals * depth, 1e-15)))

    # Bounded least-squares, not a plain root find. Free concentration can never
    # exceed its own total, and imposing that bound is what keeps the solver in
    # the physical region — an unbounded log-space search overflows exp() and
    # wanders off in the tight-binding regime.
    upper = np.log(np.maximum(totals, 1e-30))
    lower = np.full(4, np.log(1e-30))

    best = None
    for start in starts:
        start = np.clip(start, lower + 1e-9, upper - 1e-9)
        sol = least_squares(
            _residuals, start, args=args,
            bounds=(lower, upper), xtol=1e-14, ftol=1e-14, gtol=1e-14,
        )
        if best is None or sol.cost < best.cost:
            best = sol
        # Relative mass-balance tolerance. 1e-6 rather than machine-tight on
        # purpose: in deep tight binding free PTF1A falls to ~1e-12 against a
        # total of 1, so recovering it to 1e-10 demands catastrophic
        # cancellation. One part per million is far beyond what parameters with
        # four-decade priors can justify, and chasing further is false
        # precision. Step 2's explicit polynomial is the real fix — it removes
        # the cancellation instead of out-iterating it.
        if np.max(np.abs(sol.fun)) < 1e-6:
            break
    if best is None or np.max(np.abs(best.fun)) >= 1e-6:
        raise RuntimeError(
            "binding equilibrium failed to converge from any start "
            f"(p={p_tot:.3g} e={e_tot:.3g} i={i_tot:.3g} r={r_tot:.3g}, kd={kd}); "
            f"best relative residual {np.max(np.abs(best.fun)):.3e}"
        )
    p, e, i, r = np.exp(best.x)
    best_x = best.x
    return {
        "P_free": p, "E_free": e, "I_free": i, "R_free": r,
        "ID3_P": i * p / kd["IP"],      # this IS P_c — derived, never integrated
        "ID3_E": i * e / kd["IE"],
        "C_L": p * e * r / kd["L"],
        "C_J": p * e * rbpj / kd["J"],
        "_guess": best_x,              # warm-start the next RHS call
    }


def first_order_sequestration(p_tot, e_tot, i_tot, r_tot, kd, rbpj, guess=None):
    """T2: ID3 as a first-order sink. The NAIVE version, kept deliberately.

    This exists to be beaten. It is topology T2 in the Stage 2 competition, and
    the master plan's prediction is that it cannot produce the required
    ultrasensitivity. Implementing it honestly — rather than strawmanning it —
    is what makes the Q-value comparison meaningful.
    """
    # No competition: complexes form from whatever is left after a linear tax.
    p = p_tot / (1.0 + i_tot / kd["IP"])
    e = e_tot / (1.0 + i_tot / kd["IE"])
    r = r_tot
    return {
        "P_free": p, "E_free": e, "I_free": i_tot, "R_free": r,
        "ID3_P": p_tot - p,
        "ID3_E": e_tot - e,
        "C_L": p * e * r / kd["L"],
        "C_J": p * e * rbpj / kd["J"],
        "_guess": guess,
    }
