"""Gate B, second half · the fold loci bounding the persistence window.

Gate B (as reworded) asks for *"the fold loci bounding the bistable window"*. A
single bifurcation diagram gives two fold points at one parameter set; the
**locus** is how those two points move as a second parameter varies, and it is
what makes the persistence window a statement about the model rather than about
one arbitrary parameter vector.

**This is exact, not calibrated.** No sampling, no tolerance, no prior — it is
continuation, so it is a property of the equations. That is precisely why it
survives decision 013: it is the one quantitative statement about `a_P` available
with no data at all.

Method: for each value of the second parameter, re-trace the equilibrium branch in
ERK and refine both folds. Every corrector outcome goes to a ledger, and the
failure rate is reported against **both** axes — folds are where continuation is
hardest, and the fold locations are the reported number, so a scalar rate would
average away the only failures that could matter.

Times do not appear here. Where they do elsewhere, they are in `1/δ_P`.
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

from src import continuation as cont, core            # noqa: E402
from figures import _provenance as prov               # noqa: E402

OUT = ROOT / "results" / "fold_loci"
ERK_MAX = 6.0


def _folds_at(p, ledger):
    """Both fold locations in ERK, refined, for one parameter set.

    Returns `(lo, hi)` or `None` if the branch does not fold twice — which is a
    result, not an error: it means the model is monostable everywhere for that
    parameter set, and the persistence window has closed.
    """
    def f(x, erk):
        return core.rhs(0.0, np.asarray(x, float), p,
                        core.Inputs.constant(erk=float(erk)))

    s = solve_ivp(core.rhs, (0, 20_000), [3.0, 3.0, 0.02],
                  args=(p, core.Inputs.constant(erk=0.02)),
                  method="LSODA", rtol=1e-12, atol=1e-14)
    if not s.success:
        return None
    try:
        X, P = cont.continue_branch(f, s.y[:, -1], 0.02, ds=0.01, ds_max=0.04,
                                    p_bounds=(0.0, ERK_MAX), ledger=ledger,
                                    max_steps=8000)
    except cont.ContinuationError:
        return None
    idx = cont.find_folds(P)
    if len(idx) != 2:
        return None
    vals = sorted(float(cont.refine_fold(X, P, i)[0]) for i in idx)
    return vals[0], vals[1]


def locus(param_name, values, ledger):
    rows = []
    base = core.default_params()
    for v in values:
        p = replace(base, **{param_name: float(v)})
        got = _folds_at(p, ledger)
        rows.append({param_name: float(v),
                     "fold_lo": got[0] if got else None,
                     "fold_hi": got[1] if got else None,
                     "window_width": (got[1] - got[0]) if got else 0.0,
                     "bistable": got is not None})
        w = f"{rows[-1]['window_width']:.4f}" if got else "CLOSED"
        print(f"  {param_name}={v:8.4g}   window = {w}")
    return rows


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    ledger = cont.ContinuationLedger()
    base = core.default_params()

    print(f"fold locus vs a_P (b_P fixed at {base.b_P}):")
    a_P_vals = np.logspace(np.log10(1.0), np.log10(40.0), 22)
    rows_aP = locus("a_P", a_P_vals, ledger)

    print(f"\nfold locus vs b_P (a_P fixed at {base.a_P}); "
          f"B_P_CRITICAL = {core.B_P_CRITICAL}:")
    b_P_vals = np.logspace(np.log10(0.05), np.log10(2.0), 22)
    rows_bP = locus("b_P", b_P_vals, ledger)

    print(f"\n{ledger.summary()}")

    for name, rows in (("a_P", rows_aP), ("b_P", rows_bP)):
        lines = [f"{name},fold_lo,fold_hi,window_width,bistable"]
        lines += [f"{r[name]},{r['fold_lo']},{r['fold_hi']},"
                  f"{r['window_width']},{r['bistable']}" for r in rows]
        (OUT / f"locus_{name}.csv").write_text("\n".join(lines), encoding="utf-8")

    def _edge(rows, name):
        """Where the window closes — the exact bound this locus delivers."""
        bis = [r for r in rows if r["bistable"]]
        if not bis:
            return None
        return {"min": min(r[name] for r in bis), "max": max(r[name] for r in bis)}

    summary = {
        "what": "fold loci bounding the persistence window — EXACT, from "
                "continuation. No sampling, no tolerance, no prior.",
        "parameters": "core.default_params() with one group varied — PLACEHOLDER",
        "bistable_range_a_P": _edge(rows_aP, "a_P"),
        "bistable_range_b_P": _edge(rows_bP, "b_P"),
        "b_P_critical_for_trametinib_alone": core.B_P_CRITICAL,
        "n_bistable_a_P": sum(r["bistable"] for r in rows_aP),
        "n_bistable_b_P": sum(r["bistable"] for r in rows_bP),
        "corrector_calls": ledger.attempts,
        "failures": ledger.failures,
        "failure_rate_vs_erk": ledger.rate_vs_param(bins=12),
        "note": "Failure rate is reported against ERK, binned. Folds are where "
                "continuation is hardest and the fold locations are the reported "
                "number, so a scalar rate would hide the failures that matter.",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2),
                                      encoding="utf-8")
    prov.stamp_run(OUT, params={"note": "placeholder base set; the LOCUS SHAPE is "
                                        "the result, not the absolute values"},
                   note="Gate B second half: fold loci")

    print(f"\nbistable for a_P in {summary['bistable_range_a_P']}")
    print(f"bistable for b_P in {summary['bistable_range_b_P']}")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
