"""What actually limits time-to-relapse — the chromatin clock, or a fast collapse?

**Why this runs before the Phase 5 endpoint is rewritten.** Decision 014 offers
two fixes for the fact that the payload cannot reach `C`. Both are expensive, and
**both are moot if relapse is not chromatin-limited in the first place** — if the
acinar state falls over in a fast `(P, R)` collapse driven by free E-protein
crashing, then `C` is not the clock, `γ` is not the durability knob, and the
question decision 014 poses does not arise.

A panel measured `d ln(t_relapse)/d ln γ ≈ 0` at high ERK. This checks that
properly, three ways, because one slope at one operating point is not enough to
retire a state.

  1. **Sensitivity of relapse timing to `γ`**, across decades, at several
     post-withdrawal ERK levels. If `C` is the clock, `t_relapse ∝ 1/γ` and the
     log-log slope is −1. If the slope is ~0, `C` is a spectator.
  2. **Freeze `C`.** Re-run with `dC/dτ` forced to zero. If relapse timing is
     unchanged, `C`'s *dynamics* contribute nothing — the strongest form of the
     test, because it removes the mechanism rather than perturbing it.
  3. **Vary `C` at withdrawal.** If relapse timing does not depend on how much
     memory is written, `C` is not setting it.

TIME IS DIMENSIONLESS. Every duration here is in units of `1/δ_P`, and `δ_P` —
PTF1A protein turnover — is unmeasured (Bench Handshake item 9). **No conversion
to hours or days appears anywhere in this file, including in drafts.** Quoting a
time in days would require a number nobody has measured.

Placeholder parameters. The *mechanism* conclusion is the output; no duration
here is a result.
"""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src import core                                    # noqa: E402
from figures import _provenance as prov                 # noqa: E402

OUT = ROOT / "results" / "relapse_mechanism"
HOLD = 400.0            # drug-hold duration, in 1/delta_P
ERK_DRUG = 0.05
R_RELAPSE = 0.5         # RBPJL below this counts as relapsed
T_MAX = 20_000.0


def _revert(p, c0):
    """Settle under drug from a metaplastic start; return the state at withdrawal."""
    s = solve_ivp(core.rhs, (0, HOLD), [0.0, 0.0, c0],
                  args=(p, core.Inputs.constant(erk=ERK_DRUG)),
                  method="LSODA", rtol=1e-10, atol=1e-12)
    return s.y[:, -1], bool(s.success)


def _time_to_relapse(p, y0, erk_high, freeze_C=False):
    """Time from withdrawal until RBPJL falls below threshold, in 1/delta_P.

    Returns `None` if it never relapses inside `T_MAX` — reported as a
    non-relapse rather than silently coerced to the horizon, which would look
    like a very long but finite relapse and bias any slope computed from it.
    """
    inp = core.Inputs.constant(erk=erk_high)

    def f(t, y):
        dy = core.rhs(t, y, p, inp)
        if freeze_C:
            dy = np.array([dy[0], dy[1], 0.0])
        return dy

    def crossed(t, y):
        return y[1] - R_RELAPSE

    crossed.terminal = True
    crossed.direction = -1

    s = solve_ivp(f, (0, T_MAX), y0, method="LSODA", rtol=1e-10, atol=1e-12,
                  events=crossed, max_step=50.0)
    if not s.success:
        return None, False
    if s.status == 1 and len(s.t_events[0]):
        return float(s.t_events[0][0]), True
    return None, True


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    base = core.default_params()
    erk_levels = [0.75, 0.9, 1.2, 2.0, 3.0]
    gammas = np.logspace(-2.7, -0.3, 9)          # ~0.002 to 0.5
    failures = []

    # ---------------------------------------------- 1 · sensitivity to gamma
    rows = []
    for erk_high in erk_levels:
        ts, gs = [], []
        for g in gammas:
            p = replace(base, gamma=float(g))
            y_w, ok0 = _revert(p, c0=1.5)
            if not ok0 or y_w[1] < R_RELAPSE:
                failures.append({"stage": "revert", "gamma": float(g),
                                 "erk_high": erk_high, "reason": "did not revert"})
                continue
            t, ok = _time_to_relapse(p, y_w, erk_high)
            if not ok:
                failures.append({"stage": "relapse", "gamma": float(g),
                                 "erk_high": erk_high, "reason": "solver"})
                continue
            rows.append({"erk_high": erk_high, "gamma": float(g),
                         "t_relapse": t, "relapsed": t is not None})
            if t is not None:
                ts.append(t)
                gs.append(g)
        if len(ts) >= 3:
            slope = float(np.polyfit(np.log(gs), np.log(ts), 1)[0])
        else:
            slope = float("nan")
        print(f"  erk_high={erk_high:4.2f}  n={len(ts):2d}  "
              f"d ln(t_relapse)/d ln gamma = {slope:+.3f}")
        rows.append({"erk_high": erk_high, "gamma": None, "t_relapse": None,
                     "relapsed": None, "slope": slope})

    slopes = {r["erk_high"]: r["slope"] for r in rows if r.get("slope") is not None}

    # ---------------------------------------------- 2 · freeze C
    print("\nfreezing C (dC/dtau := 0):")
    frozen = []
    for erk_high in erk_levels:
        y_w, _ = _revert(base, c0=1.5)
        t_free, _ = _time_to_relapse(base, y_w, erk_high, freeze_C=False)
        t_frozen, _ = _time_to_relapse(base, y_w, erk_high, freeze_C=True)
        frozen.append({"erk_high": erk_high, "t_free": t_free,
                       "t_frozen": t_frozen})
        both = t_free is not None and t_frozen is not None
        ratio = (t_frozen / t_free) if both else None
        print(f"  erk_high={erk_high:4.2f}  free={_fmt(t_free)}  "
              f"frozen={_fmt(t_frozen)}  ratio={'%.3f' % ratio if ratio else '—'}")

    # ---------------------------------------------- 3 · vary C at withdrawal
    print("\nvarying C written before withdrawal:")
    c_rows = []
    for erk_high in (0.9, 2.0):
        for c0 in (0.1, 0.5, 1.5, 4.0):
            y_w, _ = _revert(base, c0=c0)
            t, _ = _time_to_relapse(base, y_w, erk_high)
            c_rows.append({"erk_high": erk_high, "c0": c0,
                           "C_at_withdrawal": float(y_w[2]), "t_relapse": t})
        got = [r["t_relapse"] for r in c_rows if r["erk_high"] == erk_high]
        print(f"  erk_high={erk_high:4.2f}  t_relapse over C(0) 0.1->4.0: "
              f"{[_fmt(x) for x in got]}")

    # ---------------------------------------------- verdict
    finite = [s for s in slopes.values() if np.isfinite(s)]
    chromatin_limited = bool(finite and np.median(np.abs(finite)) > 0.5)
    verdict = ("CHROMATIN-LIMITED" if chromatin_limited
               else "NOT chromatin-limited — relapse is a fast (P,R) collapse")

    summary = {
        "question": "is time-to-relapse set by the chromatin state C?",
        "units": "ALL times in 1/delta_P. delta_P is unmeasured (bench item 9). "
                 "No conversion to hours or days is possible or attempted.",
        "parameters": "core.default_params() — PLACEHOLDER",
        "drug_hold": HOLD,
        "relapse_threshold_R": R_RELAPSE,
        "dlnt_dlngamma": slopes,
        "median_abs_slope": float(np.median(np.abs(finite))) if finite else None,
        "freeze_C": frozen,
        "vary_C_at_withdrawal": c_rows,
        "chromatin_limited": chromatin_limited,
        "verdict": verdict,
        "solver_failures": failures,
        "failure_count": len(failures),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2),
                                      encoding="utf-8")
    (OUT / "gamma_sweep.csv").write_text(
        "erk_high,gamma,t_relapse_over_delta_P,relapsed\n" + "\n".join(
            f"{r['erk_high']},{r['gamma']},{r['t_relapse']},{r['relapsed']}"
            for r in rows if r.get("gamma") is not None), encoding="utf-8")
    prov.stamp_run(OUT, params={"note": "placeholder set; mechanism is the "
                                        "output, no duration here is a result"},
                   note="does C set time-to-relapse?")

    print(f"\nVERDICT: {verdict}")
    print(f"  median |d ln t / d ln gamma| = "
          f"{summary['median_abs_slope']:.3f}" if finite else "  no slopes")
    print(f"  solver failures: {len(failures)}")
    print(f"  wrote {OUT}")
    return 0


def _fmt(x):
    return "no relapse" if x is None else f"{x:.1f}"


if __name__ == "__main__":
    raise SystemExit(main())
