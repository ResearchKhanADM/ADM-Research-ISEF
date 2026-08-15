"""Continuation is validated against a problem with a known answer, first.

The whole point of pseudo-arclength continuation here is that it turns folds
where a parameter sweep cannot. That claim has to be checked on something whose
folds are known analytically, because on the real model there is nothing to
check against — a wrong continuation produces a smooth, plausible curve.

The test problem is the cusp normal form

    f(x, p) = p + x - x^3

whose equilibrium branch is `p = x^3 - x`, with folds where `dp/dx = 3x^2 - 1 = 0`,
i.e. `x = +/- 1/sqrt(3)` and `p = -/+ 2/(3*sqrt(3)) ~= -/+ 0.3849`. Both numbers are
exact, so the continuation can be graded rather than eyeballed.
"""

from __future__ import annotations

import numpy as np
import pytest

from src import continuation as cont

FOLD_P = 2.0 / (3.0 * np.sqrt(3.0))          # ~0.3849
FOLD_X = 1.0 / np.sqrt(3.0)                  # ~0.5774


def _cusp(x, p):
    return np.array([p + x[0] - x[0] ** 3])


# ---------------------------------------------------------------------------
# Against the analytic answer
# ---------------------------------------------------------------------------


def test_continuation_turns_both_folds_of_the_cusp():
    """The load-bearing capability, graded against exact fold locations."""
    led = cont.ContinuationLedger()
    X, P = cont.continue_branch(_cusp, [-1.5], -1.0, ds=0.02, ds_max=0.05,
                                p_bounds=(-2.0, 2.0), ledger=led, max_steps=4000)

    folds = cont.find_folds(P)
    assert len(folds) == 2, f"expected 2 folds, found {len(folds)}"
    got = sorted(float(P[i]) for i in folds)
    assert got[0] == pytest.approx(-FOLD_P, abs=2e-3)
    assert got[1] == pytest.approx(+FOLD_P, abs=2e-3)


def test_refined_folds_beat_the_grid_and_hit_the_analytic_state():
    """Grid resolution is not good enough for a reported number.

    `find_folds` returns the nearest sample, so its state is out by up to one
    arclength step — measured at ~3% here, on a quantity (the edge of the
    persistence window) that goes on a poster. `refine_fold` fits a parabola in
    arclength and takes its vertex. This test requires the refinement to be
    *strictly better* than the grid value, not merely close, so the extra code
    has to keep earning its place.
    """
    X, P = cont.continue_branch(_cusp, [-1.5], -1.0, ds=0.02, ds_max=0.05,
                                p_bounds=(-2.0, 2.0), max_steps=4000)
    folds = cont.find_folds(P)
    assert len(folds) == 2

    for i in folds:
        p_ref, x_ref = cont.refine_fold(X, P, i)
        err_grid = abs(abs(float(X[i][0])) - FOLD_X)
        err_ref = abs(abs(float(x_ref[0])) - FOLD_X)
        assert err_ref < err_grid, (
            f"refinement did not improve the fold state: grid err {err_grid:.4g}, "
            f"refined err {err_ref:.4g}"
        )
        assert abs(p_ref) == pytest.approx(FOLD_P, abs=2e-3)
        assert abs(float(x_ref[0])) == pytest.approx(FOLD_X, abs=6e-3)


def test_the_branch_is_non_monotone_in_the_parameter():
    """The S-curve, stated as a property.

    If `P` came back monotone the continuation would have missed the middle
    branch entirely — the failure mode of a naive sweep — while still returning a
    perfectly smooth curve.
    """
    X, P = cont.continue_branch(_cusp, [-1.5], -1.0, ds=0.02,
                                p_bounds=(-2.0, 2.0), max_steps=4000)
    assert np.any(np.diff(P) < 0) and np.any(np.diff(P) > 0)


def test_every_branch_point_actually_solves_the_system():
    """A curve that is smooth but not a solution set is the worst outcome here."""
    X, P = cont.continue_branch(_cusp, [-1.5], -1.0, ds=0.02,
                                p_bounds=(-2.0, 2.0), max_steps=4000)
    resid = np.array([abs(float(_cusp(x, p)[0])) for x, p in zip(X, P)])
    assert resid.max() < 1e-8, f"max residual {resid.max():.3e}"


def test_a_naive_parameter_sweep_cannot_do_this():
    """Guard-on-the-guard: show the thing continuation replaces actually fails.

    A sweep steps `p` and Newton-solves in `x` alone. Past the fold there is no
    nearby solution on that sheet, so it either fails or jumps branches. This
    constructs the naive version and requires it to miss the middle branch —
    otherwise the tests above would be passing for reasons unrelated to
    continuation, and the module would not have earned its complexity.
    """
    from scipy.optimize import fsolve

    xs, ps = [], []
    x = -1.5
    for p in np.linspace(-1.0, 1.0, 400):
        sol, _, ier, _ = fsolve(lambda z: _cusp(z, p), [x], full_output=True)
        if ier == 1:
            x = float(sol[0])
            xs.append(x)
            ps.append(p)
    xs = np.array(xs)

    # The middle (unstable) branch lives at |x| < 1/sqrt(3). A sweep starting on
    # the lower sheet jumps straight to the upper one at the fold and never
    # visits it.
    visited_middle = np.sum(np.abs(xs) < FOLD_X * 0.9)
    assert visited_middle <= 2, (
        f"the naive sweep visited the middle branch {visited_middle} times; if a "
        f"sweep can do that, pseudo-arclength is not buying what this module "
        f"claims it buys"
    )


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


def test_classify_separates_stable_saddle_and_unstable():
    def linear(x, p):
        return np.array([-x[0], p * x[1]])

    lab, _ = cont.classify(linear, np.zeros(2), -1.0)
    assert lab == "stable"
    lab, ev = cont.classify(linear, np.zeros(2), +1.0)
    assert lab == "saddle"
    assert cont.stable_manifold_dim(ev) == 1


def test_stable_manifold_dimension_is_counted_not_inferred():
    """In 3-D a saddle can have a 1-D or 2-D stable manifold, and the separatrix
    construction depends on which. The label alone does not distinguish them, so
    the dimension is counted from the eigenvalues."""
    def sys3(x, p):
        return np.array([-x[0], -x[1], x[2]])

    lab, ev = cont.classify(sys3, np.zeros(3), 0.0)
    assert lab == "saddle"
    assert cont.stable_manifold_dim(ev) == 2


# ---------------------------------------------------------------------------
# The standing rule: failures are counted, and counted AGAINST the parameter
# ---------------------------------------------------------------------------


def test_ledger_reports_failure_rate_as_a_function_of_the_parameter():
    """A scalar rate is explicitly not enough.

    Continuation failures concentrate at folds, and folds are the reported
    quantity. A 2% scalar rate that is 0% everywhere except a wall at the fold is
    a completely different result from 2% spread evenly, and only the binned form
    can tell them apart.
    """
    led = cont.ContinuationLedger()
    for i in range(100):
        led.record(ok=(i % 10 != 0), param=float(i) / 100.0, iters=3,
                   residual=1e-12)
    assert led.attempts == 100 and led.failures == 10
    binned = led.rate_vs_param(bins=5)
    assert len(binned) == 5
    assert sum(b["n"] for b in binned) == 100
    assert all("failure_rate" in b for b in binned)


def test_hard_regime_failure_rate_is_bounded():
    """MANDATORY (CLAUDE.md standing rule) — sweep deliberately where it is hard.

    Small `ds_min` and a tight tolerance across the full cusp, folds included.
    The bound is on the fraction of corrector calls that fail, and it is asserted
    rather than reported so a regression cannot pass quietly.
    """
    led = cont.ContinuationLedger()
    cont.continue_branch(_cusp, [-1.5], -1.0, ds=0.05, ds_max=0.08,
                         p_bounds=(-2.0, 2.0), ledger=led, max_steps=4000)
    assert led.attempts > 50
    assert led.failure_rate < 0.05, (
        f"corrector failure rate {led.failure_rate:.2%} across the cusp including "
        f"both folds; {led.summary()}"
    )


def test_jacobian_is_not_halved_by_the_state_clamp():
    """Regression: the clamp trap that returned exactly half the true slope.

    `core.rhs` clamps states with `max(y, 0)`. On the metaplastic branch `R*`
    falls far below a fixed finite-difference step (6e-11 at ERK 5.0), so the
    `x - h` evaluation landed in the clamped region, the perturbation was only
    half-applied, and the central difference returned **exactly half** the true
    slope — `J[1,1] = -0.5003` against a true `-rho = -1.0`.

    Stability *labels* survived, because the sign stayed negative, so Gate B's
    structural result was never wrong. Every eigenvalue and timescale quoted off
    that branch was wrong by 2x. It ran clean and it rendered — the class of bug
    this project keeps finding.
    """
    from scipy.integrate import solve_ivp

    from src import core

    p = core.default_params()

    def f(x, erk):
        return core.rhs(0.0, np.asarray(x, float), p,
                        core.Inputs.constant(erk=float(erk)))

    for erk in (1.6, 3.0, 5.0):
        s = solve_ivp(core.rhs, (0, 20_000), [0.0, 0.0, 3.0],
                      args=(p, core.Inputs.constant(erk=erk)),
                      method="LSODA", rtol=1e-12, atol=1e-14)
        assert s.success
        J = cont.jacobian_x(f, s.y[:, -1], erk)
        # dR/dtau = rho*(a_R*hill(...) - R), so d/dR = -rho exactly.
        assert J[1, 1] == pytest.approx(-p.rho, rel=1e-4), (
            f"J[1,1] = {J[1,1]:.6f} at ERK={erk} with R* = {s.y[1, -1]:.2e}; "
            f"expected {-p.rho}. A value near {-p.rho / 2} means the clamp is "
            f"halving the finite difference again.")


def test_jacobian_handles_a_coordinate_exactly_at_the_boundary():
    """At `x_i = 0` a symmetric step is impossible; forward must be used."""
    def g(x, p):
        return np.array([-2.0 * x[0] + p, -3.0 * x[1]])

    J = cont.jacobian_x(g, np.array([0.0, 1.0]), 0.0)
    assert J[0, 0] == pytest.approx(-2.0, rel=1e-5)
    assert J[1, 1] == pytest.approx(-3.0, rel=1e-5)
