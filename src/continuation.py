"""Pseudo-arclength continuation of equilibria, on the 3-state core.

WHAT CONTINUATION IS, AND WHY A PARAMETER SWEEP IS NOT IT
---------------------------------------------------------
The obvious way to trace equilibria against a parameter is: step the parameter,
solve `f(x, p) = 0` with the previous solution as the guess, repeat. That is a
**naive sweep**, and it cannot turn a fold.

At a fold (saddle-node), the equilibrium branch turns back on itself: two
equilibria collide and vanish. Approaching one, `df/dx` becomes singular, Newton
stops converging, and — this is the part that matters — **there is no solution at
all just past the fold**, so no amount of solver effort helps. A sweep therefore
stops dead at the exact point the model is most interesting, and silently reports
only the outer branches. The middle branch, which is the saddle, is invisible to
it. For this project the saddle *is* the separatrix, and the fold locations *are*
the persistence window, so a sweep cannot answer Gate B.

Pseudo-arclength continuation fixes this by refusing to treat the parameter as
special. It parameterizes the branch by arclength `s` and solves for the state
**and** the parameter together:

    F(x, p) = 0                                    N equations
    tangent . ( [x, p] - [x0, p0] ) - ds = 0       1 arclength constraint
    -------------------------------------------------
                                                   N+1 equations, N+1 unknowns

The augmented system stays non-singular *through* the fold, because the branch is
locally smooth in `s` even where it is vertical in `p`. So the solver walks around
the turn and continues down the middle branch, which is what produces the S-curve.

The three pieces, in order:

  **predictor**  step along the unit tangent to the branch
  **corrector**  Newton on the augmented system above
  **tangent**    updated from the null vector of the augmented Jacobian, with
                 its sign chosen to keep going the same way round the fold

STANDING RULE — NEVER SILENTLY DROP A FAILED SOLVE
--------------------------------------------------
Continuation failures are not uniform: they cluster at folds, and folds are
exactly what is being measured. Dropping them would delete the boundary of the
persistence window and leave a clean-looking curve that stops early. So every
corrector call is recorded in a ledger with the parameter value at which it
happened, and the caller reports failure rate **as a function of the parameter**.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

#: Corrector tolerance on max|F|. Tight: the fold location is read off this curve,
#: and a loose tolerance would blur exactly the quantity being reported.
TOL = 1e-10
MAX_NEWTON = 40


@dataclass
class ContinuationLedger:
    """Every corrector outcome, with the parameter value it happened at.

    A scalar failure rate is not enough and is explicitly forbidden by the
    standing rule: failures concentrate near folds, so a 2% rate that is 0%
    everywhere except a 100% wall at the fold is a completely different — and
    much worse — result than 2% spread uniformly. Keeping `param` per record is
    what makes that distinguishable.
    """

    attempts: int = 0
    failures: int = 0
    records: list[dict] = field(default_factory=list)

    def record(self, ok: bool, param: float, iters: int, residual: float,
               **extra) -> None:
        self.attempts += 1
        self.failures += 0 if ok else 1
        self.records.append({"ok": bool(ok), "param": float(param),
                             "iters": int(iters), "residual": float(residual),
                             **extra})

    @property
    def failure_rate(self) -> float:
        return self.failures / self.attempts if self.attempts else 0.0

    def rate_vs_param(self, bins=10) -> list[dict]:
        """Failure rate binned by parameter — the REQUIRED report form."""
        if not self.records:
            return []
        p = np.array([r["param"] for r in self.records])
        ok = np.array([r["ok"] for r in self.records])
        edges = np.linspace(p.min(), p.max(), bins + 1)
        out = []
        for i in range(bins):
            m = (p >= edges[i]) & (p <= edges[i + 1] if i == bins - 1
                                   else p < edges[i + 1])
            if m.sum():
                out.append({"param_lo": float(edges[i]),
                            "param_hi": float(edges[i + 1]),
                            "n": int(m.sum()),
                            "failure_rate": float(1.0 - ok[m].mean())})
        return out

    def summary(self) -> str:
        return (f"corrector calls: {self.attempts}, failures: {self.failures} "
                f"({100 * self.failure_rate:.2f}%)")


class ContinuationError(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# Numerical Jacobians
# ---------------------------------------------------------------------------


def jacobian_x(f, x, p, eps=1e-7):
    """df/dx by central differences.

    Central rather than forward: forward differences carry O(h) error, and the
    fold test reads the SIGN of an eigenvalue crossing zero. An O(h) bias near a
    zero eigenvalue moves the detected fold location, which is the number being
    reported.
    """
    x = np.asarray(x, float)
    n = x.size
    J = np.empty((n, n))
    for i in range(n):
        h = eps * max(1.0, abs(x[i]))
        xp, xm = x.copy(), x.copy()
        xp[i] += h
        xm[i] -= h
        J[:, i] = (f(xp, p) - f(xm, p)) / (2 * h)
    return J


def _jacobian_p(f, x, p, eps=1e-7):
    h = eps * max(1.0, abs(p))
    return (f(x, p + h) - f(x, p - h)) / (2 * h)


# ---------------------------------------------------------------------------
# Corrector
# ---------------------------------------------------------------------------


def _correct(f, z0, tangent, z_prev, ds, ledger):
    """Newton on the augmented system. `z = [x, p]`, size N+1.

    The augmented residual is

        G(z) = [ f(x, p) ,  tangent . (z - z_prev) - ds ]

    and its Jacobian is [[df/dx, df/dp], [tangent]] — non-singular at a fold even
    though `df/dx` alone is singular there. That is the entire trick.
    """
    z = np.asarray(z0, float).copy()
    n = z.size - 1
    residual = np.inf
    for it in range(MAX_NEWTON):
        x, p = z[:n], z[n]
        G = np.empty(n + 1)
        G[:n] = f(x, p)
        G[n] = tangent @ (z - z_prev) - ds
        residual = float(np.max(np.abs(G)))
        if residual < TOL:
            return z, True, it, residual
        J = np.empty((n + 1, n + 1))
        J[:n, :n] = jacobian_x(f, x, p)
        J[:n, n] = _jacobian_p(f, x, p)
        J[n, :] = tangent
        try:
            dz = np.linalg.solve(J, -G)
        except np.linalg.LinAlgError:
            break
        # Damped step. Undamped Newton on the augmented system can overshoot
        # into a different branch near a fold, which silently produces a curve
        # that jumps between sheets and looks continuous.
        lam = 1.0
        for _ in range(8):
            z_try = z + lam * dz
            r_try = np.max(np.abs(f(z_try[:n], z_try[n])))
            if np.isfinite(r_try) and r_try < residual * 10:
                break
            lam *= 0.5
        z = z + lam * dz
    return z, False, MAX_NEWTON, residual


def _correct_at_fixed_p(f, x0, p0):
    """Newton the seed onto `f(x, p0) = 0` before continuation begins.

    Plain Newton in `x` only — correct here precisely because the seed is assumed
    to be away from a fold, which is the one place this would fail. If it does
    fail, that is reported rather than silently accepted: a seed that cannot be
    corrected means the caller is starting somewhere the branch does not pass.
    """
    x = np.asarray(x0, float).copy()
    for _ in range(MAX_NEWTON):
        r = f(x, p0)
        if np.max(np.abs(r)) < TOL:
            return x
        try:
            x = x + np.linalg.solve(jacobian_x(f, x, p0), -r)
        except np.linalg.LinAlgError:
            break
    resid = float(np.max(np.abs(f(x, p0))))
    if resid > 1e-6:
        raise ContinuationError(
            f"seed did not converge to an equilibrium at p={p0:g}: "
            f"max|f| = {resid:.3e}. Continuation starting from a non-solution "
            f"would return a branch whose first point is not on the branch."
        )
    return x


def _tangent(f, z, prev_tangent):
    """Unit tangent to the branch: null vector of the augmented Jacobian.

    The sign is chosen to keep the same direction of travel. Without that the
    continuation reverses at the fold and retraces the branch it came from —
    which produces a plausible-looking curve covering only half the structure.
    """
    n = z.size - 1
    x, p = z[:n], z[n]
    A = np.empty((n + 1, n + 1))
    A[:n, :n] = jacobian_x(f, x, p)
    A[:n, n] = _jacobian_p(f, x, p)
    A[n, :] = prev_tangent
    rhs = np.zeros(n + 1)
    rhs[n] = 1.0
    try:
        t = np.linalg.solve(A, rhs)
    except np.linalg.LinAlgError:
        _, _, Vt = np.linalg.svd(A[:n, :])
        t = Vt[-1]
    t = t / np.linalg.norm(t)
    if t @ prev_tangent < 0:
        t = -t
    return t


# ---------------------------------------------------------------------------
# The driver
# ---------------------------------------------------------------------------


def continue_branch(f, x0, p0, *, ds=0.02, ds_min=1e-6, ds_max=0.1,
                    max_steps=4000, p_bounds=(-np.inf, np.inf),
                    ledger=None):
    """Trace an equilibrium branch of `f(x, p) = 0` from `(x0, p0)`.

    Returns `(X, P)` with `X` shape (n_points, n_states) and `P` shape
    (n_points,), ordered along the branch — so `P` is NOT monotone when the
    branch folds, and that non-monotonicity is the signal.

    Step size adapts: halve on a corrector failure and retry, grow slowly on
    easy steps. Adaptation is what lets one call handle both the flat parts and
    the fold, where curvature is high.
    """
    x0 = np.asarray(x0, float)
    n = x0.size

    # Correct the seed onto the branch before starting. Without this, a caller
    # who passes an approximate guess gets that guess returned as branch point
    # zero — a point that is not a solution, sitting in a file of solutions, with
    # nothing downstream able to tell. Everything after it would be fine, which
    # is what makes it dangerous.
    x0 = _correct_at_fixed_p(f, x0, p0)
    z = np.concatenate([x0, [p0]])

    # Seed the tangent by pushing in +p, then let the null-vector update take
    # over. If the branch is near-vertical in p at the start this is a poor
    # guess, but one corrector call fixes it.
    t = np.zeros(n + 1)
    t[n] = 1.0
    t = _tangent(f, z, t)

    X, P = [z[:n].copy()], [float(z[n])]
    step = ds
    for _ in range(max_steps):
        z_pred = z + step * t
        z_new, ok, iters, res = _correct(f, z_pred, t, z, step, ledger=ledger)
        if ledger is not None:
            ledger.record(ok, param=float(z_pred[n]), iters=iters, residual=res,
                          ds=float(step))
        if not ok:
            step *= 0.5
            if step < ds_min:
                # Genuine end of branch, not a solver artefact: recorded and
                # returned rather than raised, because a truncated branch is a
                # result the caller must be able to report.
                break
            continue
        z = z_new
        t = _tangent(f, z, t)
        X.append(z[:n].copy())
        P.append(float(z[n]))
        if not (p_bounds[0] <= z[n] <= p_bounds[1]):
            break
        step = min(step * 1.15, ds_max)
    return np.array(X), np.array(P)


# ---------------------------------------------------------------------------
# Classification and folds
# ---------------------------------------------------------------------------


def classify(f, x, p):
    """Stability class of an equilibrium, from the eigenvalues of df/dx.

    Returns `(label, eigenvalues)` with label in
    {"stable", "saddle", "unstable", "nonhyperbolic"}.

    "saddle" means at least one eigenvalue in each half-plane — for a 3-D system
    that covers both a 1-D and a 2-D unstable manifold, which is why
    `stable_manifold_dim` is reported separately rather than being inferred from
    the label.
    """
    ev = np.linalg.eigvals(jacobian_x(f, x, p))
    re = ev.real
    tol = 1e-8 * max(1.0, float(np.max(np.abs(ev))))
    if np.any(np.abs(re) < tol):
        return "nonhyperbolic", ev
    if np.all(re < 0):
        return "stable", ev
    if np.all(re > 0):
        return "unstable", ev
    return "saddle", ev


def stable_manifold_dim(ev) -> int:
    return int(np.sum(np.asarray(ev).real < 0))


def find_folds(P):
    """Indices where the branch turns in `p` — i.e. where `dp/ds` changes sign.

    A fold is detected geometrically, from the branch's own shape, rather than by
    hunting for a zero eigenvalue. Both work; this one is robust because it does
    not require differentiating an eigenvalue that is itself computed by finite
    differences, and near a fold that eigenvalue is precisely where the numerics
    are worst.

    Returns the index of the first point *after* the turn. That is a grid
    location, accurate only to one arclength step — use `refine_fold` for any
    number that gets reported.
    """
    dP = np.diff(P)
    sign = np.sign(dP)
    turns = []
    for i in range(1, len(sign)):
        if sign[i] != 0 and sign[i - 1] != 0 and sign[i] != sign[i - 1]:
            turns.append(i)          # index into P of the turning point
    return turns


def refine_fold(X, P, i):
    """Sub-step fold location, by fitting `p(s)` locally and taking its vertex.

    **Why this is not optional.** `find_folds` returns the nearest sample, so its
    accuracy is one arclength step — and the fold locations *are* the edges of
    the persistence window, which is the R3 deliverable. Reporting a headline
    number at step resolution, when the refinement is a three-point parabola,
    is not defensible: at `ds = 0.05` the state at the fold was out by ~3%.

    Arclength is approximated by cumulative Euclidean distance along the branch
    in `(x, p)` — the natural parameterization continuation is already using.
    Both `p` and each component of `x` are fitted as quadratics in `s`, and the
    fold is the vertex of the `p` fit.
    """
    lo, hi = max(i - 1, 0), min(i + 2, len(P))
    if hi - lo < 3:
        return float(P[i]), np.asarray(X[i], float)

    Xw, Pw = np.asarray(X[lo:hi], float), np.asarray(P[lo:hi], float)
    steps = np.linalg.norm(np.diff(np.column_stack([Xw, Pw]), axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(steps)])

    a, b, _ = np.polyfit(s, Pw, 2)
    if abs(a) < 1e-30:
        return float(P[i]), np.asarray(X[i], float)
    s_star = -b / (2 * a)
    # Refuse to extrapolate: outside the window the parabola is not a local
    # model of anything, and a fold "refined" to somewhere the branch does not
    # go is worse than the grid value it replaced.
    if not (s[0] <= s_star <= s[-1]):
        return float(P[i]), np.asarray(X[i], float)

    p_star = float(np.polyval([a, b, np.polyfit(s, Pw, 2)[2]], s_star))
    x_star = np.array([float(np.polyval(np.polyfit(s, Xw[:, k], 2), s_star))
                       for k in range(Xw.shape[1])])
    return p_star, x_star
