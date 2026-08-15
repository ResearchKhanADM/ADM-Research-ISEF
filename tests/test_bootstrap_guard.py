"""THE FLAGSHIP GUARD. `dR/dt` has no P-independent term.

This file was written **before** `src/core.py`, deliberately, because the claim it
guards is the one the whole project rests on:

    *Rbpjl* has no PTF1A-independent promoter. Once the PTF1A-RBPJL loop opens it
    cannot re-close on its own. That bootstrap failure is why MEK-inhibitor
    reversion is drug-dependent, and it is what the payload is for.

In the model that claim is **a zero** — the absence of a basal production term in
`dR/dt`. An absence is the most fragile thing a codebase can carry: adding
`+ beta_R0` looks like a reasonable leak term, makes the model *better behaved*
numerically, and silently converts the project's central claim into its opposite.
Nothing crashes. Every figure still renders. The conclusion inverts.

Two tests, and the second exists because the first is not enough on its own.

  1. **Structural.** With `P = 0`, `dR/dt` must have no positive contribution.
  2. **Guard-on-the-guard.** Construct the violating implementation and *require
     test 1 to fail on it*. An assertion that passes for both the correct and the
     incorrect version tests nothing — this is the same discipline that was
     applied to `raf_drive`'s sign in the retired 11-state model, and it is the
     part people leave out.
  3. **Dynamic companion.** `R` must not rise from zero while `P` is held at zero,
     integrated over a long horizon. A structural check inspects one derivative;
     this catches an accidental basal term that only shows up once the system is
     actually run — e.g. one hidden inside an input, or one that only activates
     away from the initial condition.

Run:  venv\\Scripts\\python.exe -m pytest tests/test_bootstrap_guard.py -v
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from src import core


# ---------------------------------------------------------------------------
# The guard itself, factored out so it can be pointed at a deliberately broken
# implementation as well as the real one. That reuse IS the guard-on-the-guard:
# the same code must pass on `core.rhs` and fail on the violating version.
# ---------------------------------------------------------------------------


def _assert_no_p_independent_r_production(rhs_fn, n=2000, seed=0):
    """Fail unless `dR/dt <= 0` everywhere `P = 0` and no payload is delivered.

    The logic: with P at zero the PTF1-L complex cannot form, so the *only*
    admissible terms in `dR/dt` are decay (negative) and the delivered payload
    (switched off here). Any strictly positive value is a P-independent source.

    Swept randomly over parameters, states and ERK level rather than checked at
    one point, because a basal term could easily be small enough to hide at a
    single arbitrary parameter set.
    """
    rng = np.random.default_rng(seed)
    idx_R = core.STATES.index("R")

    worst = -np.inf
    worst_at = None
    for _ in range(n):
        p = core.Params(
            a_P=10 ** rng.uniform(-1, 1.5), b_P=10 ** rng.uniform(-3, 0),
            c_rep=10 ** rng.uniform(-1, 1), n_P=rng.uniform(1, 4),
            a_R=10 ** rng.uniform(-1, 1.5), n_R=rng.uniform(1, 4),
            rho=10 ** rng.uniform(-1, 1), gamma=10 ** rng.uniform(-3, 0),
            alpha_C=10 ** rng.uniform(-1, 1), k_w=10 ** rng.uniform(-1, 1),
            eps=10 ** rng.uniform(-2, 0.5), kappa=10 ** rng.uniform(-4, 2),
            n_C=1.0,
        )
        # P is pinned at zero; R and C roam, ERK roams, payload is OFF.
        y = np.array([0.0, 10 ** rng.uniform(-6, 1), 10 ** rng.uniform(-6, 1)])
        erk = 10 ** rng.uniform(-3, 1)

        dy = rhs_fn(0.0, y, p, core.Inputs.constant(erk=erk))
        if dy[idx_R] > worst:
            worst, worst_at = dy[idx_R], (p, y, erk)

    assert worst <= 0.0, (
        "dR/dt is POSITIVE while P = 0 with no payload delivered. That is a "
        "P-independent source of RBPJL, which is precisely the term this model "
        "claims does not exist. The bootstrap argument — and with it the reason "
        f"reversion is drug-dependent — is void.\n  max dR/dt = {worst!r}\n"
        f"  at {worst_at!r}"
    )


# ---------------------------------------------------------------------------
# 1 · Structural
# ---------------------------------------------------------------------------


def test_dR_dt_has_no_P_independent_term():
    """THE claim, as a zero. See the module docstring."""
    _assert_no_p_independent_r_production(core.rhs)


def test_R_production_is_exactly_zero_at_the_origin():
    """At `P = R = 0` with no payload, `dR/dt` must be *exactly* 0.0.

    Stronger than `<= 0` and worth stating separately: decay contributes nothing
    at `R = 0`, so anything other than exact zero is a source term, however
    small. `== 0.0` rather than `approx`, on purpose — a basal term of 1e-14 is
    still a basal term, and floating-point noise cannot enter here because every
    admissible contribution is an exact product with a zero factor.
    """
    p = core.default_params()
    idx_R = core.STATES.index("R")
    for erk in (0.0, 0.1, 1.0, 10.0, 1e3):
        for c in (0.0, 0.5, 5.0):
            dy = core.rhs(0.0, np.array([0.0, 0.0, c]), p,
                          core.Inputs.constant(erk=erk))
            assert dy[idx_R] == 0.0, (
                f"dR/dt = {dy[idx_R]!r} at P = R = 0 (erk={erk}, C={c}); "
                f"expected exactly 0.0"
            )


# ---------------------------------------------------------------------------
# 2 · Guard-on-the-guard — the part people leave out
# ---------------------------------------------------------------------------


def _violating_rhs(t, y, p, inp):
    """`core.rhs` plus the bug: a small P-independent basal RBPJL source.

    This is written to look *reasonable*, because that is the point. A modeller
    adding `+ beta_R0` to keep `R` off the boundary, or to represent "leaky
    transcription", would produce exactly this and would not think of it as
    changing a claim. The magnitude is deliberately tiny — 1e-6 — so this also
    checks the guard is not merely catching gross errors.
    """
    dy = np.array(core.rhs(t, y, p, inp), dtype=float)
    dy[core.STATES.index("R")] += 1e-6
    return dy


def test_the_guard_would_actually_catch_a_violation():
    """The guard must FAIL on the violating implementation.

    Without this, `test_dR_dt_has_no_P_independent_term` could be passing for a
    reason unrelated to the model — a sign convention, an index mix-up, a
    tolerance so loose that nothing can trip it. A monotonicity or
    non-positivity assertion that holds for both the correct and the incorrect
    version is not evidence about the model; it is decoration.
    """
    with pytest.raises(AssertionError, match="P-independent source"):
        _assert_no_p_independent_r_production(_violating_rhs, n=50)


def test_the_guard_catches_a_violation_hidden_in_the_complex_term():
    """A subtler violation: RBPJL production that does not vanish with P.

    The realistic version of this bug is not an added constant — it is a Hill
    function written with the wrong argument, or an offset inside it, so that
    `H(0) > 0`. The result is a *P-dependent-looking* term that nonetheless fires
    at `P = 0`, which is much harder to spot by reading the code than `+ beta_R0`.
    """
    def offset_hill_rhs(t, y, p, inp):
        dy = np.array(core.rhs(t, y, p, inp), dtype=float)
        # H(x + 0.1) instead of H(x): looks like the real term, but H(0) != 0.
        dy[core.STATES.index("R")] += p.rho * p.a_R * core.hill(0.1, p.n_R)
        return dy

    with pytest.raises(AssertionError, match="P-independent source"):
        _assert_no_p_independent_r_production(offset_hill_rhs, n=50)


# ---------------------------------------------------------------------------
# 3 · Dynamic companion
# ---------------------------------------------------------------------------


def _integrate_with_P_clamped_to_zero(p, inp, r0, c0, t_end):
    """Integrate (R, C) with P held at zero throughout.

    P is clamped both in the state passed to the right-hand side and in the
    returned derivative, so nothing can raise it — this asks "if PTF1A never
    comes back, can RBPJL return on its own?", which is the bootstrap claim
    stated as an experiment rather than as an algebraic property.
    """
    idx_P = core.STATES.index("P")

    def clamped(t, y):
        y = np.asarray(y, dtype=float).copy()
        y[idx_P] = 0.0
        dy = np.array(core.rhs(t, y, p, inp), dtype=float)
        dy[idx_P] = 0.0
        return dy

    return solve_ivp(clamped, (0.0, t_end), [0.0, r0, c0],
                     method="LSODA", rtol=1e-10, atol=1e-12,
                     dense_output=False, t_eval=np.linspace(0, t_end, 200))


def test_R_cannot_rise_from_zero_while_P_is_held_at_zero():
    """The claim, run rather than inspected.

    Long horizon, ERK swept from fully suppressed (as under trametinib) to high.
    **Trametinib being on must not be enough** — that is the entire point of the
    project. If `R` recovers here, the model says MEK inhibition alone re-closes
    the loop, and the payload has no reason to exist.
    """
    p = core.default_params()
    for erk in (0.0, 1e-3, 0.1, 1.0, 10.0):
        sol = _integrate_with_P_clamped_to_zero(
            p, core.Inputs.constant(erk=erk), r0=0.0, c0=0.5, t_end=5_000.0)
        assert sol.success, f"integration failed at erk={erk}: {sol.message}"
        r = sol.y[core.STATES.index("R")]
        assert np.max(r) <= 1e-9, (
            f"R rose to {np.max(r):.3e} from zero with P clamped at zero "
            f"(erk={erk}). RBPJL bootstrapped itself; the loop re-closed without "
            f"PTF1A. The model no longer makes the claim the project is about."
        )


def test_R_decays_to_zero_from_a_high_start_while_P_is_held_at_zero():
    """The complement, and it is not implied by the test above.

    Starting `R` at zero and staying there is consistent with `dR/dt = 0`
    everywhere — a model where `R` is simply frozen would pass. This requires
    the loop to actually *collapse*: with PTF1A gone, RBPJL must decay away
    rather than persist, which is what makes the metaplastic state an attractor
    the payload has to push the cell out of.
    """
    p = core.default_params()
    sol = _integrate_with_P_clamped_to_zero(
        p, core.Inputs.constant(erk=1.0), r0=2.0, c0=0.5, t_end=5_000.0)
    assert sol.success
    r_final = sol.y[core.STATES.index("R")][-1]
    assert r_final < 1e-6, (
        f"R settled at {r_final:.3e} with P clamped at zero; RBPJL must decay "
        f"away once its only driver is gone"
    )


def test_payload_is_the_only_route_back_for_R():
    """Positive control on the guard: `u_R` MUST be able to raise `R`.

    Every test above asserts that something cannot happen. Taken alone they
    would all pass on a model where `R` is inert — which would be a broken model
    that satisfies the letter of the bootstrap claim while making the payload
    useless. This checks the delivered payload still works, so the suite
    distinguishes "no P-independent term" from "no term at all".
    """
    p = core.default_params()
    inp = core.Inputs.constant(erk=1.0, u_R=0.5)
    dy = core.rhs(0.0, np.array([0.0, 0.0, 0.5]), p, inp)
    assert dy[core.STATES.index("R")] > 0.0, (
        "delivered RBPJL mRNA must raise R; if it cannot, the payload has no "
        "route into the model and every downstream result is vacuous"
    )
