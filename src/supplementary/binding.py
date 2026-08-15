"""Three-complex competitive binding — the term that decides bistability.

ID3 traps BOTH partners of an obligate heterodimer. Verified: Dufresne 2010
(Int J Cancer 129(2):295-306, PMID 20830706) reports gastrin raising Id3 and
increasing both Id3/E47 AND Id3/Ptf1-p48 interactions while decreasing
E47/Ptf1-p48 — in AR4-2J, which is the wet-lab line.

FULL ALGEBRA: docs/derivations/binding_polynomial.md. This module follows that
document step for step, and the test suite checks the two against each other.

The short version. Each conservation law factorises, because every species
appears linearly in every complex it belongs to:

    P_tot = P*[1 + I/K_IP + E*Phi]
    E_tot = E*[1 + I/K_IE + P*Phi]
    I_tot = I*[1 + E/K_IE + P/K_IP]
    R_tot = R*[1 + P*E/K_L]      ->  Phi = R_tot/(K_L + P*E) + J/K_J

so we solve P = P_tot/(partition function), NOT the residual
P + bound - P_tot = 0. That distinction is not stylistic. In the tight-binding
regime free P falls to ~1e-12 against a total of ~1, so the residual form
recovers P by catastrophic cancellation. The first implementation did exactly
that: it could not do better than 1e-6 relative and failed to converge outright
below Kd ~ 0.01 — which is precisely where the mechanism lives. The ratio form
computes small P directly and reaches machine precision everywhere.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import root

# Machine-precision residuals are achievable in the ratio form, so this is tight
# on purpose: a loose tolerance would hide exactly the failures the convergence
# ledger exists to count.
TOL = 1e-12
_FIXED_POINT_SWEEPS = 15


@dataclass
class ConvergenceLedger:
    """Counts every binding solve, converged or not.

    STANDING RULE (CLAUDE.md): never silently drop a failed solve. At 1e5
    samples across five topologies, convergence failures correlate with Kd,
    because conditioning worsens as binding tightens. Dropping them would
    deplete the surviving sample set precisely in the regime where T1 and T2
    differ — biasing the Q-value comparison AGAINST the discriminating regime,
    and looking like a clean negative result.
    """

    attempts: int = 0
    failures: int = 0
    failed_conditions: list[dict] = field(default_factory=list)

    @property
    def failure_rate(self) -> float:
        return self.failures / self.attempts if self.attempts else 0.0

    def record(self, ok: bool, **conditions) -> None:
        self.attempts += 1
        if not ok:
            self.failures += 1
            # Keep the conditions, not just the count: the required report is
            # failure rate AS A FUNCTION OF Kd/E_tot and every swept parameter.
            self.failed_conditions.append(conditions)

    def summary(self) -> str:
        return (f"binding solves: {self.attempts}, failures: {self.failures} "
                f"({100 * self.failure_rate:.3f}%)")


class BindingSolveError(RuntimeError):
    """Raised on non-convergence. Never caught-and-ignored inside a sweep:
    catch at the sweep level, record in the ledger, report the rate."""


def _phi(p, e, r_tot, kd, rbpj):
    """Shared complex-forming capacity, after eliminating R (derivation §3).

    Both ternary complexes enter only here, which encodes that PTF1-L and PTF1-J
    compete for the same P*E dimer pool.
    """
    return r_tot / (kd["L"] + p * e) + rbpj / kd["J"]


def _fixed_point(p_tot, e_tot, i_tot, r_tot, kd, rbpj, sweeps=_FIXED_POINT_SWEEPS):
    """Cheap sweeps of the ratio form to land in Newton's basin (derivation §4).

    Converges from any positive start, but slowly when binding is tight (~6000
    iterations at Kd=1e-5), so it is a seed and never used alone.
    """
    p, e, i = p_tot, e_tot, i_tot
    for _ in range(sweeps):
        phi = _phi(p, e, r_tot, kd, rbpj)
        p = p_tot / (1.0 + i / kd["IP"] + e * phi)
        e = e_tot / (1.0 + i / kd["IE"] + p * phi)
        i = i_tot / (1.0 + e / kd["IE"] + p / kd["IP"])
    return np.log([p, e, i])


def _log_residuals(z, p_tot, e_tot, i_tot, r_tot, kd, rbpj):
    """Log-space residuals of the ratio form. All terms O(1), no cancellation.

    z is clipped to the physical bound (a free concentration cannot exceed its
    own total) before exponentiating. This does not move the root — the
    constraint holds at any true solution — but it stops Newton's rejected trial
    steps from overflowing exp() and spraying warnings that would mask a real
    numerical problem later.
    """
    z = np.minimum(z, np.log([p_tot, e_tot, i_tot]))
    p, e, i = np.exp(z)
    phi = _phi(p, e, r_tot, kd, rbpj)
    return np.array([
        z[0] - np.log(p_tot) + np.log1p(i / kd["IP"] + e * phi),
        z[1] - np.log(e_tot) + np.log1p(i / kd["IE"] + p * phi),
        z[2] - np.log(i_tot) + np.log1p(e / kd["IE"] + p / kd["IP"]),
    ])


def solve_binding(p_tot, e_tot, i_tot, r_tot, kd, rbpj, guess=None, ledger=None):
    """Free species and complexes at binding equilibrium.

    Quasi-equilibrium is justified by timescale separation: association is
    seconds-to-minutes against transcription at hours-to-weeks. Safe here for
    the reason the fast-acetylation elimination was not: eliminating H collapsed
    Hill(C_L,H) into Hill(C_L,C_L) and the variable vanished algebraically,
    whereas these species retain a nontrivial nonlinear dependence on the slow
    states. Elimination is dangerous when it destroys structure, not when it is
    merely fast.
    """
    args = (p_tot, e_tot, i_tot, r_tot, kd, rbpj)
    z0 = np.asarray(guess, float) if guess is not None else _fixed_point(*args)

    sol = root(_log_residuals, z0, args=args, method="hybr", tol=TOL)
    ok = bool(sol.success) and np.max(np.abs(sol.fun)) < TOL

    if not ok and guess is not None:
        # A stale warm start can land outside the basin. Retry cold before
        # calling it a failure, or the ledger records solver bookkeeping as
        # though it were a property of the parameter set.
        z0 = _fixed_point(*args)
        sol = root(_log_residuals, z0, args=args, method="hybr", tol=TOL)
        ok = bool(sol.success) and np.max(np.abs(sol.fun)) < TOL

    if not ok:
        # ESCALATION LADDER. Newton is fast but not globally convergent, and its
        # failures concentrate in the tight-binding regime — exactly where T1
        # and T2 differ. A bare Newton failed on ~4% of tight-regime samples,
        # which is the precise mechanism that would bias Q-values against the
        # discriminating regime. The fixed-point map on the ratio form is slow
        # but converges from any positive start, so it is the safety net.
        z0 = _fixed_point(*args, sweeps=4000)
        sol = root(_log_residuals, z0, args=args, method="hybr", tol=TOL)
        ok = bool(sol.success) and np.max(np.abs(sol.fun)) < TOL

    if not ok:
        # Last resort: pure fixed point, no Newton. Guaranteed to converge in
        # this system, just slowly. Cost is irrelevant on the <<1% of samples
        # that reach here; a dropped sample is not.
        z_fp = _fixed_point(*args, sweeps=200_000)
        if np.max(np.abs(_log_residuals(z_fp, *args))) < TOL:
            sol = type(sol)(x=z_fp, fun=_log_residuals(z_fp, *args), success=True)
            ok = True

    if ledger is not None:
        ledger.record(
            ok,
            # Kd/E_tot is the discriminating dimensionless ratio (derivation §7),
            # so it is what the failure rate must be reported against.
            kd_over_etot=kd["IE"] / e_tot if e_tot > 0 else np.inf,
            kd_IE=kd["IE"], kd_IP=kd["IP"],
            p_tot=p_tot, e_tot=e_tot, i_tot=i_tot, r_tot=r_tot,
        )
    if not ok:
        raise BindingSolveError(
            f"binding equilibrium did not converge (p={p_tot:.3g} e={e_tot:.3g} "
            f"i={i_tot:.3g} r={r_tot:.3g}, Kd_IE={kd['IE']:.3g}); "
            f"max|residual|={np.max(np.abs(sol.fun)):.3e}"
        )

    p, e, i = np.exp(sol.x)
    r = r_tot / (1.0 + p * e / kd["L"])          # derivation §3
    return {
        "P_free": p, "E_free": e, "I_free": i, "R_free": r,
        "ID3_P": i * p / kd["IP"],               # this IS P_c — derived, never integrated
        "ID3_E": i * e / kd["IE"],
        "C_L": p * e * r / kd["L"],
        "C_J": p * e * rbpj / kd["J"],
        "_guess": sol.x,                         # warm-start the next RHS call
    }


def first_order_sequestration(p_tot, e_tot, i_tot, r_tot, kd, rbpj,
                              guess=None, ledger=None):
    """T2: ID3 as a first-order sink. The naive version, kept deliberately.

    This exists to be beaten, and implementing it honestly rather than as a
    strawman is what makes the Q-value comparison meaningful.

    Derivation §5 shows something worth knowing: in the LOOSE-binding limit the
    full titration model reduces to exactly this — T1 and T2 are then not merely
    similar but the same model, both with log-log slope -2. They separate only
    when Kd << totals.

    Closed form, so it cannot fail to converge. That is itself an asymmetry the
    convergence ledger must account for: if T1 loses samples to non-convergence
    where T2 cannot, the Q-values are not comparable until that is reported.
    """
    p = p_tot / (1.0 + i_tot / kd["IP"])
    e = e_tot / (1.0 + i_tot / kd["IE"])
    r = r_tot
    if ledger is not None:
        ledger.record(True, kd_over_etot=kd["IE"] / e_tot if e_tot > 0 else np.inf,
                      kd_IE=kd["IE"], kd_IP=kd["IP"], p_tot=p_tot, e_tot=e_tot,
                      i_tot=i_tot, r_tot=r_tot)
    return {
        "P_free": p, "E_free": e, "I_free": i_tot, "R_free": r,
        "ID3_P": p_tot - p, "ID3_E": e_tot - e,
        "C_L": p * e * r / kd["L"],
        "C_J": p * e * rbpj / kd["J"],
        "_guess": guess,
    }
