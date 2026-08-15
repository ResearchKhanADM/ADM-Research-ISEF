"""Pin a golden trajectory from the 11-state implementation before it is deleted.

WHY THIS EXISTS, AND WHAT IT IS NOT.

The Phase 2 rewrite drops from 11 states / 61 parameters to 3-4 states / ~9-12,
so the new core **cannot** reproduce these numbers and is not expected to. What
survives a change of state space is *qualitative*: does the core still have two
attractors, is the metaplastic one still reachable from an acinar start under
high KRAS, and are the timescales still hours-to-weeks rather than seconds or
years. This file records the old answer so that comparison is possible at all.
An hour now; impossible once `model.py` is gone.

**THE PARAMETER SET BELOW IS ARBITRARY.** No parameter table was ever written
(old Stage 0 Step 3 never ran), and none of these values is measured. They are
order-of-magnitude choices on a nondimensional-ish scale, fixed with a seed so
the record is reproducible. Nothing here is a scientific claim, no number from
this file may be reported, and it must never be cited as a model prediction.
It is a regression fixture, in the same sense as a checksum.

Run:  venv\\Scripts\\python.exe scripts/pin_golden_trajectory.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src import model                                    # noqa: E402
from src.topology import CORE_UNTREATED                  # noqa: E402
from figures import _provenance as prov                  # noqa: E402

GOLDEN = ROOT / "tests" / "golden"

# Time in hours. The old model's slowest process is chromatin memory, and the
# validation targets were days-to-weeks, so 1500 h (~9 weeks) is long enough for
# the slow manifold to be visibly settled rather than still in transient.
T_END = 1500.0
N_OUT = 300


def parameters() -> dict:
    """One fixed, arbitrary parameter set. See the module docstring."""
    return {
        # --- ERK / drug layer. v=0 throughout (payload and drug both off), so
        # f_cat and f_act are identically 1 and only K_kras and W matter.
        "K_kras": 1.0, "ic50_cat": 1.0, "ic50_act": 1.0,
        "k_fb": 1.0, "r_fb": 2.0,          # raf_drive: DECREASING in K_eff
        "k_on_W": 1.0, "k_off_W": 2.0,     # tau_W ~ 30 min, the fast end

        # --- binding. Tight relative to totals (~1), which is the regime where
        # two-target titration is ultrasensitive and the model can be bistable.
        "kd": {"IE": 0.02, "IP": 0.05, "L": 0.5, "J": 5.0},
        "rbpj_total": 0.5,

        # --- PTF1A
        "alpha_ign": 0.30, "k_ign": 0.30, "q_ign": 3.0,
        "alpha_auto": 2.00, "k_auto": 0.25, "n_auto": 3.0,
        "k_M_rep": 0.40, "n_M_rep": 3.0,
        "delta_P": 0.30,

        # --- RBPJL. NO ignition term in the RHS; that zero is the bootstrap claim.
        "beta_R": 1.20, "k_R": 0.25, "m_R": 3.0, "delta_R": 0.25,

        # --- E-protein pool
        "beta_E": 0.30, "delta_E": 0.20,

        # --- ID3
        "beta_I0": 0.02, "beta_I": 1.00, "k_I": 0.30, "n_I": 2.0, "delta_I": 0.50,

        # --- chromatin memory at metaplasia loci
        "k_w": 0.05, "k_wM": 0.30, "a_wM": 3.0,
        "k_CM": 0.30, "b_CM": 3.0,
        "k_e": 0.02, "eps_M": 0.03, "theta_M": 0.50,

        # --- acinar output and secretory load
        "k_A": 1.00, "k_AC": 0.25, "n_AC": 2.0, "delta_A": 0.20,
        "gamma0": 0.50, "k_cargo": 0.30,

        # --- viability hazard (retired by decision 008 amendment; kept here
        # only because the old RHS integrates cumhaz as a state)
        "h_max_cargo": 0.02, "s_crit": 2.0, "nu_s": 4.0,
        "h_max_chop": 0.02, "p_crit": 0.05, "mu_p": 4.0,
    }


def initial_conditions(topo) -> dict:
    """Two starts: differentiated acinar, and metaplastic. Order = topo.states."""
    idx = topo.index
    n = len(topo.states)

    acinar = np.zeros(n)
    acinar[idx("P_n")] = 1.20      # loop closed
    acinar[idx("R")] = 1.00
    acinar[idx("E_tot")] = 1.50
    acinar[idx("I")] = 0.05        # titrator low
    acinar[idx("M")] = 0.02        # metaplasia chromatin near-blank
    acinar[idx("A")] = 1.50
    acinar[idx("S")] = 0.50
    acinar[idx("W")] = 0.20

    adm = np.zeros(n)
    adm[idx("P_n")] = 0.05         # loop open
    adm[idx("R")] = 0.02           # RBPJL is the deepest hole
    adm[idx("E_tot")] = 1.50
    adm[idx("I")] = 1.20           # titrator high
    adm[idx("M")] = 0.80           # metaplasia chromatin written
    adm[idx("A")] = 0.05
    adm[idx("S")] = 0.02
    adm[idx("W")] = 0.50

    return {"acinar": acinar, "adm": adm}


def integrate(y0, p, topo):
    cache: dict = {}
    sol = solve_ivp(
        model.rhs_for_solver, (0.0, T_END), y0, args=(p, topo, cache),
        method="LSODA",            # stiff: binding is fast, chromatin is slow
        t_eval=np.linspace(0.0, T_END, N_OUT),
        rtol=1e-8, atol=1e-10,
    )
    if not sol.success:
        # STANDING RULE: never silently drop a failed solve. A golden fixture
        # built from a half-integrated trajectory is worse than none.
        raise RuntimeError(f"golden integration failed: {sol.message}")
    return sol


def main() -> int:
    GOLDEN.mkdir(parents=True, exist_ok=True)
    topo = CORE_UNTREATED                    # payload and drug both off
    p = parameters()
    states = list(topo.states)

    written, summary = [], {}
    for name, y0 in initial_conditions(topo).items():
        sol = integrate(y0, p, topo)
        path = GOLDEN / f"golden_{name}.csv"
        header = "t," + ",".join(states)
        np.savetxt(path, np.column_stack([sol.t, sol.y.T]), delimiter=",",
                   header=header, comments="", fmt="%.10g")
        written.append(path)
        summary[name] = {s: float(sol.y[i, -1]) for i, s in enumerate(states)}
        print(f"{name:7s} -> {path.name}")
        print("          final " + "  ".join(
            f"{s}={summary[name][s]:.4g}" for s in ("P_n", "R", "M", "A")))

    # The qualitative property the rewrite should be compared against: do the
    # two starts land in DIFFERENT places? Reported as measured, never tuned to.
    sep = max(abs(summary["acinar"][s] - summary["adm"][s]) for s in states)
    bistable = sep > 1e-3
    print(f"\nendpoints differ by {sep:.4g} -> "
          f"{'TWO attractors' if bistable else 'ONE attractor (monostable)'}")

    manifest = {
        "purpose": "regression fixture for the Phase 2 rewrite; NOT a result",
        "warning": "parameter set is arbitrary and unmeasured; report no number "
                   "from this file",
        "implementation": "11-state composable-topology model, deleted after this",
        "topology": topo.name,
        "states": states,
        "t_end_hours": T_END,
        "parameters": p,
        "final_states": summary,
        "endpoint_separation": sep,
        "two_attractors": bool(bistable),
        "files": {f.name: prov.sha256(f) for f in written},
    }
    (GOLDEN / "manifest.json").write_text(json.dumps(manifest, indent=2),
                                          encoding="utf-8")
    prov.stamp_run(GOLDEN, seed=None, params={"note": "arbitrary fixture set"},
                   outputs=written + [GOLDEN / "manifest.json"],
                   note="golden trajectory pinned before the Phase 2 rewrite")
    print(f"wrote {GOLDEN / 'manifest.json'} and _run.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
