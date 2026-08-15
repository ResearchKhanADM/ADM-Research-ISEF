"""GATE B — bistability, the saddle, the separatrix, and the persistence window.

Gate B asks: *two stable states plus a saddle, with an identifiable separatrix, in
a model whose key parameters survive profile likelihood.* This script does the
first three; profile likelihood is a separate stage.

WHAT IT COMPUTES

  1. **The equilibrium branch** in the ERK input, by pseudo-arclength continuation
     (`src/continuation.py`). Continuation rather than a parameter sweep because a
     sweep cannot turn a fold, and the folds are the answer — see that module's
     docstring.
  2. **Stability class of every point** on the branch, from Jacobian eigenvalues.
  3. **The two folds**, which bracket the bistable window. **That window is the
     drug-free persistence window** — the R3 deliverable — expressed in the one
     scalar the drug and the oncogene both act through.
  4. **The separatrix**: the saddle's 2-D stable manifold, grown by integrating
     `-f` outward from a small circle in the saddle's stable eigenplane. Points
     on this surface flow *into* the saddle forward in time, so backward
     integration traces the surface itself. In 3-D it is a genuine dividing
     surface, not a curve.
  5. **Convergence failures as a function of ERK**, never as a scalar — the
     standing rule. Continuation failures cluster at folds, and folds are exactly
     what is being measured, so a scalar rate would hide the only failures that
     would matter.

NOTHING HERE IS A CALIBRATED RESULT. It runs at `core.default_params()`, which is
a placeholder set. What Gate B tests is whether the STRUCTURE exists and is
computable, not where the window sits in real units.

Run:  venv\\Scripts\\python.exe scripts/run_gate_b.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src import continuation as cont, core            # noqa: E402
from figures import _provenance as prov               # noqa: E402

OUT = ROOT / "results" / "gate_b"
ERK_MAX = 5.0
SEED = 20260815


def _rhs_of_erk(p):
    """f(x, erk) with everything else fixed — the continuation's right-hand side."""
    def f(x, erk):
        return core.rhs(0.0, np.asarray(x, float), p,
                        core.Inputs.constant(erk=float(erk)))
    return f


def _seed_equilibrium(p, erk, y0):
    """Settle to an equilibrium to start continuation from."""
    s = solve_ivp(core.rhs, (0, 20_000), y0, args=(p, core.Inputs.constant(erk=erk)),
                  method="LSODA", rtol=1e-12, atol=1e-14)
    if not s.success:
        raise RuntimeError(f"seed integration failed: {s.message}")
    return s.y[:, -1]


def trace_branch(p, ledger):
    """The full S-curve, from the acinar branch at low ERK around both folds."""
    f = _rhs_of_erk(p)
    x0 = _seed_equilibrium(p, 0.05, [3.0, 3.0, 0.02])
    resid = float(np.max(np.abs(f(x0, 0.05))))
    if resid > 1e-8:
        raise RuntimeError(f"seed is not an equilibrium: |f| = {resid:.3e}")
    X, P = cont.continue_branch(f, x0, 0.05, ds=0.01, ds_max=0.04,
                                p_bounds=(0.0, ERK_MAX), ledger=ledger,
                                max_steps=8000)
    labels, sdims = [], []
    for x, erk in zip(X, P):
        lab, ev = cont.classify(f, x, erk)
        labels.append(lab)
        sdims.append(cont.stable_manifold_dim(ev))
    return X, P, labels, sdims


def separatrix(p, x_saddle, erk, n_rays=48, radius=1e-4, t_back=400.0,
               domain=None):
    """The saddle's stable manifold, grown backwards from its stable eigenplane.

    The saddle here has a 2-D stable manifold in a 3-D state space, so `W^s` is a
    **surface** and it is what separates the two basins. Construction: take the
    two stable eigenvectors, lay a small circle of radius `radius` in the plane
    they span, and integrate `-f` from each point. Backward time carries them out
    along `W^s`, sweeping the surface.

    `radius` is small on purpose. The eigenplane is only tangent to `W^s` at the
    saddle, so seeding far out starts on the wrong surface — a systematic error
    that produces a smooth, plausible, wrong separatrix.

    **`domain` is not cosmetic.** In backward time the two stable directions
    become expanding (eigenvalues negate), which is exactly what sweeps out the
    manifold — and also means the trajectory grows exponentially without bound.
    `W^s` extends to infinity; the first version of this ran to ~1e24 and
    produced a figure whose axes were meaningless. So each ray terminates when it
    leaves a bounding box drawn around the states the model actually occupies.
    The manifold is not being truncated arbitrarily: it is being reported over
    the region where the model is defined, which is the only region a basin
    boundary means anything in.

    Every ray's integration outcome is returned so failures can be counted rather
    than dropped.
    """
    f = _rhs_of_erk(p)
    J = cont.jacobian_x(f, x_saddle, erk)
    ev, evec = np.linalg.eig(J)
    stable = np.where(ev.real < 0)[0]
    if stable.size != 2:
        raise RuntimeError(
            f"expected a 2-D stable manifold at the saddle, found dim={stable.size}; "
            f"eigenvalues {ev}"
        )
    v1, v2 = np.real(evec[:, stable[0]]), np.real(evec[:, stable[1]])
    v1 /= np.linalg.norm(v1)
    v2 = v2 - (v2 @ v1) * v1
    v2 /= np.linalg.norm(v2)

    def backward(t, y):
        return -core.rhs(t, y, p, core.Inputs.constant(erk=erk))

    hi = np.asarray(domain, float)

    def leave_domain(t, y):
        """Terminal event: 0 at the edge of the PHYSICAL box `[0, hi]`.

        Two-sided, and the lower bound matters as much as the upper. `core.rhs`
        clamps negative states to zero, so below the origin the integrator is
        following the *clamped* dynamics, which are not a smooth extension of the
        model. A manifold traced through there is an artefact of the clamp, not
        geometry — and it looked entirely plausible in the first draft of the
        figure, extending to RBPJL = -7.8.
        """
        y = np.asarray(y, float)
        return float(min(np.min(hi - y), np.min(y + 1e-3)))

    leave_domain.terminal = True
    leave_domain.direction = -1

    rays, outcomes = [], []
    for th in np.linspace(0, 2 * np.pi, n_rays, endpoint=False):
        y0 = x_saddle + radius * (np.cos(th) * v1 + np.sin(th) * v2)
        s = solve_ivp(backward, (0, t_back), y0, method="LSODA",
                      rtol=1e-9, atol=1e-12, events=leave_domain, max_step=1.0)
        left = s.status == 1          # 1 == a terminal event fired
        outcomes.append({"theta": float(th), "ok": bool(s.success),
                         "left_domain": bool(left), "t_end": float(s.t[-1]),
                         "message": str(s.message)})
        if s.success and s.y.shape[1] > 1:
            rays.append(s.y.T)

    # COVERAGE, not just termination. "48/48 rays terminated" says nothing about
    # whether they SPAN the manifold. In backward time the fast stable direction
    # expands ~165x faster than the slow one, so a uniform seed circle is
    # stretched into a sliver and the rays collapse onto the two branches of the
    # IN-PLANE 1-D stable manifold. Measured here as the largest angular gap
    # between adjacent ray endpoints, projected back onto the eigenplane: uniform
    # coverage would give 360/n_rays.
    gaps = None
    if rays:
        ang = []
        for r in rays:
            d = r[-1] - x_saddle
            ang.append(np.arctan2(d @ v2, d @ v1))
        a = np.sort(np.asarray(ang))
        gaps = float(np.max(np.diff(np.concatenate([a, [a[0] + 2 * np.pi]]))))
    return rays, outcomes, (v1, v2), ev, gaps


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    p = core.default_params()
    ledger = cont.ContinuationLedger()

    print("tracing equilibrium branch (pseudo-arclength continuation in ERK)...")
    X, P, labels, sdims = trace_branch(p, ledger)
    print(f"  {len(P)} points, ERK in [{P.min():.4f}, {P.max():.4f}]")
    print(f"  {ledger.summary()}")

    # Refined, not grid-resolution: these two numbers ARE the edges of the
    # persistence window, i.e. the R3 deliverable, so they are reported to
    # sub-step accuracy rather than to the nearest continuation sample.
    folds = cont.find_folds(P)
    refined = [cont.refine_fold(X, P, i) for i in folds]
    fold_erk = sorted(float(e) for e, _ in refined)
    fold_grid = sorted(float(P[i]) for i in folds)
    print(f"  folds at ERK = {[round(e, 5) for e in fold_erk]} "
          f"(grid: {[round(e, 5) for e in fold_grid]})")
    if len(fold_erk) != 2:
        print(f"  !! expected 2 folds, found {len(fold_erk)} — GATE B NOT MET")

    counts = {k: labels.count(k) for k in set(labels)}
    print(f"  stability classes: {counts}")

    # --- the separatrix, at the midpoint of the bistable window
    sep_rays, sep_outcomes, basis, sep_ev, sep_gap = [], [], None, None, None
    erk_mid = None
    if len(fold_erk) == 2:
        erk_mid = 0.5 * (fold_erk[0] + fold_erk[1])
        saddle_idx = [i for i, (lab, e) in enumerate(zip(labels, P))
                      if lab == "saddle" and abs(e - erk_mid) < 0.02]
        if saddle_idx:
            i = saddle_idx[len(saddle_idx) // 2]
            # Bounding box: twice the largest equilibrium value the branch
            # reaches in each state, so the manifold is reported over the region
            # the model actually occupies rather than out to numerical infinity.
            box = 2.0 * np.max(np.abs(X), axis=0)
            print(f"growing separatrix from the saddle at ERK={P[i]:.4f} "
                  f"(domain box {np.round(box, 3)}) ...")
            sep_rays, sep_outcomes, basis, sep_ev, sep_gap = separatrix(
                p, X[i], P[i], domain=box)
            n_ok = sum(o["ok"] for o in sep_outcomes)
            n_left = sum(o["left_domain"] for o in sep_outcomes)
            print(f"  {n_ok}/{len(sep_outcomes)} rays integrated, "
                  f"{n_left} reached the domain edge")
            print(f"  angular coverage: largest gap {np.degrees(sep_gap):.1f} deg "
                  f"(uniform would be {360/len(sep_outcomes):.1f}) -> "
                  f"{'COLLAPSED onto the in-plane manifold' if np.degrees(sep_gap) > 45 else 'spread'}")
        else:
            print("  !! no saddle found near the window midpoint")

    # ---------------------------------------------------------------- write
    np.savetxt(OUT / "branch.csv",
               np.column_stack([P, X, sdims]), delimiter=",",
               header="erk,P,R,C,stable_manifold_dim", comments="", fmt="%.10g")
    (OUT / "branch_labels.csv").write_text(
        "erk,label\n" + "\n".join(f"{e:.10g},{l}" for e, l in zip(P, labels)),
        encoding="utf-8")

    if sep_rays:
        rows = ["ray,step,P,R,C"]
        for k, r in enumerate(sep_rays):
            for j, pt in enumerate(r):
                rows.append(f"{k},{j},{pt[0]:.10g},{pt[1]:.10g},{pt[2]:.10g}")
        (OUT / "separatrix.csv").write_text("\n".join(rows), encoding="utf-8")

    # Failure rate AS A FUNCTION OF the swept parameter — the standing rule. A
    # scalar here would hide the only failures that could matter, since they
    # concentrate at the folds.
    (OUT / "convergence.json").write_text(json.dumps({
        "corrector_calls": ledger.attempts,
        "failures": ledger.failures,
        "failure_rate_scalar_DO_NOT_REPORT_ALONE": ledger.failure_rate,
        "failure_rate_vs_erk": ledger.rate_vs_param(bins=12),
        "separatrix_ray_outcomes": sep_outcomes,
        "note": "Continuation failures concentrate at folds, which are exactly "
                "what is measured. Report the binned rate, never the scalar.",
    }, indent=2), encoding="utf-8")

    summary = {
        "gate": "B",
        "parameters": "core.default_params() — PLACEHOLDER, not calibrated",
        "n_branch_points": len(P),
        "erk_range": [float(P.min()), float(P.max())],
        "folds_erk": fold_erk,
        "folds_erk_grid_resolution": fold_grid,
        "persistence_window_erk": (
            [fold_erk[0], fold_erk[1]] if len(fold_erk) == 2 else None),
        "persistence_window_width_erk": (
            fold_erk[1] - fold_erk[0] if len(fold_erk) == 2 else None),
        "stability_counts": counts,
        "saddle_present": "saddle" in counts,
        "saddle_stable_manifold_dim": 2,
        "separatrix_rays": len(sep_rays),
        "separatrix_erk": erk_mid,
        "separatrix_max_angular_gap_deg": (float(np.degrees(sep_gap))
                                           if sep_gap is not None else None),
        "separatrix_WHAT_IS_ACTUALLY_COMPUTED": (
            "The IN-PLANE separatrix at C = C*: the 1-D stable manifold of the "
            "saddle within the invariant plane {C = C*}. NOT the full 2-D W^s. "
            "C is decoupled from (P,R) — the dC row of the Jacobian is (0,0,.) — "
            "so {C = C*} is invariant, and in backward time the fast stable "
            "direction expands ~165x faster than the slow one, collapsing a "
            "uniform seed circle onto the in-plane branches. The angular gap "
            "above measures that collapse; a large gap means coverage of the "
            "2-D manifold was NOT achieved and must not be claimed."),
        "gate_b_structure_met": bool(len(fold_erk) == 2 and "saddle" in counts
                                     and len(sep_rays) > 0),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2),
                                      encoding="utf-8")
    prov.stamp_run(OUT, seed=SEED,
                   params={"note": "placeholder parameter set; no number from "
                                   "this run is a calibrated result"},
                   note="Gate B: continuation, folds, saddle, separatrix")

    print("\n" + ("GATE B STRUCTURE MET" if summary["gate_b_structure_met"]
                  else "GATE B NOT MET"))
    print(f"  persistence window in ERK: {summary['persistence_window_erk']}")
    print(f"  wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
