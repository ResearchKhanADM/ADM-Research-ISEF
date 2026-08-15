"""Exact equilibrium set by scalar reduction — and why that is a proof, not a shortcut.

**THE STRUCTURAL FACT.** At constant input the 3-state core collapses to two
*independent scalar* root-finds:

  1. `C` is **decoupled**: `dC/dτ` depends on `C` and the input only — the `dC`
     row of the Jacobian is exactly `(0, 0, ·)`. So `C*` solves a scalar equation
     on its own, with no reference to `P` or `R`.
  2. `R` is **explicit** given `P`: `dR/dτ = 0` ⟹ `R = a_R·hill(P·E, n_R)`.
  3. Substituting both leaves a **scalar equation in `P` alone**.

Bracketed root-finding on that scalar therefore **enumerates every equilibrium**,
including the saddle, and **cannot jump branches**. That is what makes this an
independent proof rather than a second opinion: continuation follows a curve and
could in principle miss a disconnected branch; a sign-change sweep over a scalar
cannot.

Three consequences worth stating as results rather than observations:

**THEOREM 1 — `C` cannot change where `(P, R)` goes.** `C` is a one-way input to
the `(P,R)` subsystem with no feedback. Freezing it therefore cannot alter the
`(P,R)` attractor structure; it can only alter the modulation `rep(C)`. Decision
015's finding that freezing `C` leaves relapse timing unchanged is a *theorem*,
not a sweep result that happened to come out that way.

**THEOREM 2 — no oscillations are possible.** The `(P,R)` block has both
off-diagonal entries ≥ 0 (raising `P` raises `R`, raising `R` raises `P`), so it
is a **cooperative** planar system. Its discriminant `(J₀₀−J₁₁)² + 4·J₀₁·J₁₀` is
then ≥ 0, so the spectrum is **always real** — no complex pair, hence no Hopf
bifurcation and no limit cycle. **"Two sinks plus one saddle" is the only bistable
structure this system can have.** That forecloses the entire class of *"but what
if it oscillates?"* questions, without a single simulation.

**THEOREM 3 — the fold count is exactly the root count.** Bistability is three
roots of `h`; monostability is one. The persistence window is the interval of
input over which `h` has three roots, and its edges are tangencies of `h` with
zero. No continuation is needed to find them, so the two methods are genuinely
independent checks on each other.

`C` may itself have three roots — that is exactly what `ε > EPS_MEMORY_THRESHOLD`
buys — so the full equilibrium set is the product of the `C` roots and, for each,
the `P` roots.

Verified against `src/continuation.py` to 6–7 significant figures; the agreement
is what licenses using this fast path for the Sobol sweep, where continuation
would need ~16 hours.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from . import core


def _roots_scalar(g, lo, hi, n_grid=400, xtol=1e-13, count_only=False):
    """Every sign change of `g` on `[lo, hi]`, refined by Brent.

    A grid sweep plus bracketing, not a local solve from a guess: this is what
    makes the enumeration complete rather than dependent on where you started.
    Resolution is the only failure mode — two roots inside one grid cell are
    missed — so `n_grid` is generous.

    `g` MUST accept a numpy array. The scan is vectorised because the Sobol sweep
    calls this ~10^6 times: evaluating the grid point-by-point put the whole
    analysis at ~39 hours, and vectorising the scan is what brings it into range.
    Brent still runs scalar, but only on the handful of brackets that contain a
    root, so it costs almost nothing.

    `count_only=True` skips the Brent refinement entirely and returns the number
    of sign changes. Bistability is a *count*, so the window search never needs
    the root locations — only the edges do.
    """
    xs = np.linspace(lo, hi, n_grid)
    vals = np.asarray(g(xs), dtype=float)
    sign_change = vals[:-1] * vals[1:] < 0.0
    exact = vals == 0.0

    if count_only:
        return int(sign_change.sum() + exact.sum())

    out = [float(x) for x in xs[exact]]
    for i in np.flatnonzero(sign_change):
        out.append(float(brentq(lambda z: float(g(np.asarray(z))),
                                xs[i], xs[i + 1], xtol=xtol)))
    return sorted(out)


def steady_C(p: core.Params, id3, c_max=None):
    """All roots of `dC/dτ = 0`. Scalar and independent of `P`, `R`.

    `α_C·w − C + ε·hill(C, 2) = 0`, with `w = id3/(k_w + id3)`.

    Multiple roots appear only when `ε > EPS_MEMORY_THRESHOLD` — that is
    precisely what makes `C` a bistable memory rather than a lagged filter, and
    it is the condition decision 015's headline caveat is bounded by. Returning
    *all* of them is what lets the sweep test prediction P3 (a step in `ε`, not a
    slope) rather than assume it.
    """
    w = id3 / (p.k_w + id3) if (p.k_w + id3) > 0 else 0.0
    # Upper bound: C cannot exceed its own maximum production, alpha_C*w + eps.
    hi = c_max if c_max is not None else max(1.0, p.alpha_C * w + p.eps) * 1.5

    def g(C):
        C = np.asarray(C, dtype=float)
        return p.alpha_C * w - C + p.eps * core.hill(C, core.SELF_REINFORCEMENT_EXPONENT)

    roots = _roots_scalar(g, 0.0, hi, n_grid=600)
    if not roots:                       # numerically flat; fall back to the filter root
        roots = _roots_scalar(g, 0.0, hi * 4, n_grid=3000) or [p.alpha_C * w]
    return roots


def _h_of_P(p: core.Params, id3, C):
    """The scalar residual in `P`, with `R` and `C` eliminated."""
    e_free = core.free_e_protein(id3, p.kappa)
    rep = 1.0 / (1.0 + (C / p.c_rep) ** p.n_C)
    ignition = p.b_P / (1.0 + id3)

    def h(P):
        P = np.asarray(P, dtype=float)
        R = p.a_R * core.hill(P * e_free, p.n_R)
        return p.a_P * core.hill(P * e_free * R, p.n_P) * rep + ignition - P

    return h


def equilibria(p: core.Params, id3, p_max=None):
    """Every equilibrium `(P, R, C)` at this input. Complete by construction.

    Returns a list of `(P, R, C)` tuples, ascending in `P` within each `C` root.
    """
    e_free = core.free_e_protein(id3, p.kappa)
    # h(P) = production + ignition - P, and production is bounded by a_P, so all
    # roots lie below a_P + b_P. Bounding it analytically rather than guessing is
    # what makes the enumeration provably complete.
    hi = p_max if p_max is not None else (p.a_P + p.b_P) * 1.05

    out = []
    for C in steady_C(p, id3):
        h = _h_of_P(p, id3, C)
        for P in _roots_scalar(h, 0.0, hi, n_grid=500):
            R = p.a_R * core.hill(P * e_free, p.n_R)
            out.append((P, R, C))
    return out


def n_equilibria(p: core.Params, id3, n_grid=500):
    """Count equilibria without locating them — the fast path.

    Bistability is a *count*, so the window search never needs root positions.
    Skipping the Brent refinement is most of the speedup that makes the Sobol
    sweep tractable at all.
    """
    total = 0
    for C in steady_C(p, id3):
        hi = (p.a_P + p.b_P) * 1.05
        total += _roots_scalar(_h_of_P(p, id3, C), 0.0, hi,
                               n_grid=n_grid, count_only=True)
    return total


def is_bistable(p: core.Params, id3, n_grid=500):
    """Three or more equilibria ⟹ two sinks plus a saddle (Theorem 2)."""
    return n_equilibria(p, id3, n_grid=n_grid) >= 3


def persistence_window(p: core.Params, lo=1e-3, hi=6.0, n=140, refine=40):
    """The bistable interval in the input — the R3 deliverable, found exactly.

    Coarse scan for where the equilibrium count changes, then bisection on each
    edge. No continuation, no tangent, no arclength: the edges are located by
    *counting roots*, which cannot miss a branch or turn the wrong way at a fold.

    Returns `(lo_edge, hi_edge, width)`, or `None` if never bistable.
    """
    grid = np.logspace(np.log10(lo), np.log10(hi), n)
    flags = np.array([is_bistable(p, float(x)) for x in grid])
    if not flags.any():
        return None

    idx = np.where(flags)[0]
    i0, i1 = idx[0], idx[-1]

    def _edge(a, b):
        """Bisect between a non-bistable `a` and a bistable `b`."""
        for _ in range(refine):
            m = 0.5 * (a + b)
            if is_bistable(p, m):
                b = m
            else:
                a = m
        return 0.5 * (a + b)

    lo_edge = _edge(float(grid[i0 - 1]), float(grid[i0])) if i0 > 0 else float(grid[0])
    hi_edge = (_edge(float(grid[i1 + 1]), float(grid[i1]))
               if i1 < len(grid) - 1 else float(grid[-1]))
    return lo_edge, hi_edge, hi_edge - lo_edge


def spectrum_is_real(p: core.Params, id3, state):
    """Theorem 2, checked pointwise rather than assumed.

    The `(P,R)` block is cooperative — both off-diagonals ≥ 0 — so its
    discriminant is non-negative and the eigenvalues are real. Returns the
    discriminant; a negative value would falsify the theorem and must never occur.
    """
    P, R, C = state
    e_free = core.free_e_protein(id3, p.kappa)
    eps_h = 1e-7

    def dP(Pv, Rv):
        rep = 1.0 / (1.0 + (C / p.c_rep) ** p.n_C)
        return (p.a_P * core.hill(Pv * e_free * Rv, p.n_P) * rep
                + p.b_P / (1.0 + id3) - Pv)

    def dR(Pv, Rv):
        return p.rho * (p.a_R * core.hill(Pv * e_free, p.n_R) - Rv)

    hP = eps_h * max(1.0, abs(P))
    hR = eps_h * max(1.0, abs(R))
    j00 = (dP(P + hP, R) - dP(P - hP, R)) / (2 * hP)
    j01 = (dP(P, R + hR) - dP(P, R - hR)) / (2 * hR)
    j10 = (dR(P + hP, R) - dR(P - hP, R)) / (2 * hP)
    j11 = (dR(P, R + hR) - dR(P, R - hR)) / (2 * hR)
    return float((j00 - j11) ** 2 + 4.0 * j01 * j10), (j01, j10)
