"""Phase 2 global sensitivity — pre-registered, stratified, checkpointed.

Pre-registration: `prereg/2026-08-15_phase2_sensitivity_{ranges.yaml,prediction.md}`,
**pushed before this ran**. Predictions P1–P6 are Luqmaan's; P5's operationalisation
was changed and flagged in that file because P5 and P6 contradicted each other.

DESIGN, per the prereg:
  * Sobol (Saltelli) over the **10 continuous groups only**. `n_P` and `n_R` are a
    **stratum label**, not an uncertain input — sampling them would convert a
    declared scan into a prior, which decision 004 refused. **P6 holds by
    construction** and stays a live implementation check.
  * Repeated at each point of the `(n_P, n_R)` scan.
  * Integration horizon `t_end = max(4000, 40/γ)`. A fixed horizon manufactures
    failures in the slow-memory regime, which is exactly where P3 is tested.

**P2 CARRIES A STOP CONDITION.** If `γ` shows `S1 > 0.15` on any durability
output, that contradicts decision 015's test 3 — stop and report, treating it as a
**bug hypothesis first, a finding second**.

Times are in `1/δ_P` throughout. There is no clock (bench item 9), and no
conversion to hours appears anywhere.

Long job: checkpoints every `--chunk` samples to `results/sobol/`, `--resume`
picks up from the last checkpoint. **Launch detached**, never as a session
background process.

Run:  venv\\Scripts\\python.exe scripts/run_sobol.py [--resume] [--n-base 512]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src import core, equilibria as eq                # noqa: E402
from figures import _provenance as prov               # noqa: E402

PREREG = ROOT / "prereg" / "2026-08-15_phase2_sensitivity_ranges.yaml"
OUT = ROOT / "results" / "sobol"


def load_prereg():
    spec = yaml.safe_load(PREREG.read_text(encoding="utf-8"))
    names = spec["method"]["sobol_inputs"]
    bounds, logmask = [], []
    for n in names:
        e = spec["parameters"][n]
        is_log = e["scale"] == "log"
        logmask.append(is_log)
        bounds.append([np.log10(e["low"]), np.log10(e["high"])] if is_log
                      else [e["low"], e["high"]])
    return spec, names, np.array(bounds), np.array(logmask)


def to_params(names, logmask, row, n_P, n_R):
    kw = {}
    for name, is_log, v in zip(names, logmask, row):
        kw[name] = float(10.0 ** v if is_log else v)
    return core.Params(n_P=float(n_P), n_R=float(n_R), n_C=1.0, **kw)


# ---------------------------------------------------------------------------
# Deliverable-side outputs
# ---------------------------------------------------------------------------


def evaluate(p: core.Params):
    """Return the deliverable outputs for one parameter vector, plus an outcome.

    Everything here comes from the **exact scalar reduction** (`src/equilibria.py`),
    not from integration: the equilibrium set is enumerated by bracketed
    root-finding, which cannot miss a branch. That is what makes ~200k samples
    tractable at all, and it is independently checked against continuation to
    2e-7 on the upper fold.
    """
    try:
        win = eq.persistence_window(p, n=48, refine=22)
    except Exception as exc:                       # noqa: BLE001 — recorded, never dropped
        return None, f"window:{type(exc).__name__}"

    if win is None:
        # Not bistable anywhere: a real outcome, width 0. NOT a failure — coercing
        # it to one would deplete the sample exactly where the model says there is
        # no persistence window, which is a result.
        return {"persistence_window_width": 0.0,
                "window_lo": np.nan, "window_hi": np.nan,
                "n_C_roots": len(eq.steady_C(p, 1.0)),
                "bistable": 0.0}, "ok"

    lo, hi, width = win
    return {"persistence_window_width": float(width),
            "window_lo": float(lo), "window_hi": float(hi),
            "n_C_roots": len(eq.steady_C(p, 1.0)),
            "bistable": 1.0}, "ok"


OUTPUT_KEYS = ["persistence_window_width", "window_hi", "bistable", "n_C_roots"]


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-base", type=int, default=512)
    ap.add_argument("--chunk", type=int, default=2000)
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()

    from SALib.analyze import sobol as sobol_analyze
    from SALib.sample import sobol as sobol_sample

    OUT.mkdir(parents=True, exist_ok=True)
    spec, names, bounds, logmask = load_prereg()
    strata_nP = spec["method"]["strata"]["n_P"]
    strata_nR = spec["method"]["strata"]["n_R"]

    # Plain Python floats, not numpy scalars: SALib builds its internal arrays
    # from this dict, and numpy scalars make them object-dtype, which fails deep
    # inside `analyze` with an opaque "setting an array element with a sequence".
    problem = {"num_vars": len(names), "names": list(names),
               "bounds": [[float(lo), float(hi)] for lo, hi in bounds]}
    X = sobol_sample.sample(problem, a.n_base, calc_second_order=False,
                            seed=spec["method"]["seed"])
    per_stratum = len(X)
    strata = [(nP, nR) for nP in strata_nP for nR in strata_nR]
    print(f"{per_stratum} samples x {len(strata)} strata = "
          f"{per_stratum * len(strata):,} evaluations", flush=True)

    ckpt = OUT / "checkpoint.npz"
    done, Y, outcomes = 0, {}, []
    if a.resume and ckpt.exists():
        z = np.load(ckpt, allow_pickle=True)
        done = int(z["done"])
        Y = {k: list(z[f"Y_{k}"]) for k in OUTPUT_KEYS}
        outcomes = list(z["outcomes"])
        print(f"resuming from {done:,}", flush=True)
    else:
        Y = {k: [] for k in OUTPUT_KEYS}

    total = per_stratum * len(strata)
    t0 = time.time()
    for idx in range(done, total):
        s, i = divmod(idx, per_stratum)
        nP, nR = strata[s]
        p = to_params(names, logmask, X[i], nP, nR)
        out, status = evaluate(p)
        outcomes.append(status)
        for k in OUTPUT_KEYS:
            Y[k].append(float(out[k]) if out else np.nan)

        if (idx + 1) % a.chunk == 0 or idx + 1 == total:
            np.savez(ckpt, done=idx + 1, outcomes=np.array(outcomes, dtype=object),
                     **{f"Y_{k}": np.array(Y[k]) for k in OUTPUT_KEYS})
            el = time.time() - t0
            rate = (idx + 1 - done) / max(el, 1e-9)
            print(f"  {idx+1:,}/{total:,}  {rate:.0f}/s  "
                  f"eta {(total-idx-1)/max(rate,1e-9)/60:.0f} min", flush=True)

    # ------------------------------------------------------------- analyse
    results = {}
    for k in OUTPUT_KEYS:
        y = np.array(Y[k])
        per = {}
        for s, (nP, nR) in enumerate(strata):
            seg = y[s * per_stratum:(s + 1) * per_stratum]
            # A variance decomposition of a constant is undefined, and SALib
            # fails opaquely rather than saying so. A near-constant output is
            # itself informative — it means every parameter set in this stratum
            # gives the same answer — so it is RECORDED, not silently skipped.
            if np.all(np.isnan(seg)):
                per[f"n_P={nP},n_R={nR}"] = {"status": "all-NaN"}
                continue
            scale = max(abs(float(np.nanmean(seg))), 1e-12)
            if np.nanstd(seg) / scale < 1e-9:
                per[f"n_P={nP},n_R={nR}"] = {
                    "status": "constant", "value": float(np.nanmean(seg)),
                    "note": "no variance to decompose; every sampled parameter "
                            "set in this stratum gives the same value"}
                continue
            seg = np.ascontiguousarray(
                np.nan_to_num(seg, nan=float(np.nanmean(seg))), dtype=np.float64)
            # Fresh dict per call: `sample()` mutates `problem` (adds
            # `sample_scaled`), and reusing the mutated one has bitten before.
            fresh = {"num_vars": len(names), "names": list(names),
                     "bounds": [[float(lo), float(hi)] for lo, hi in bounds]}
            try:
                si = sobol_analyze.analyze(fresh, seg, calc_second_order=False,
                                           print_to_console=False)
            except Exception as exc:               # noqa: BLE001
                print(f"  analyze failed for {k} @ {nP},{nR}: "
                      f"{type(exc).__name__} {exc} | seg {seg.shape} {seg.dtype} "
                      f"expected {a.n_base * (len(names) + 2)}", flush=True)
                continue
            per[f"n_P={nP},n_R={nR}"] = {
                "S1": dict(zip(names, [float(v) for v in si["S1"]])),
                "ST": dict(zip(names, [float(v) for v in si["ST"]])),
                "mean": float(np.mean(seg)), "std": float(np.std(seg)),
            }
        results[k] = per

    # ------------------------------------------------- P2 stop condition
    p2_hits = []
    for k, per in results.items():
        if k not in ("persistence_window_width", "window_hi"):
            continue
        for stratum, d in per.items():
            if "S1" not in d:
                continue
            if d["S1"].get("gamma", 0.0) > 0.15:
                p2_hits.append({"output": k, "stratum": stratum,
                                "S1_gamma": d["S1"]["gamma"]})

    # ------------------------------------------- failures vs every parameter
    ok = np.array([o == "ok" for o in outcomes])
    fail_tbl = {}
    for j, nm in enumerate(names):
        col = np.tile(X[:, j], len(strata))[:len(ok)]
        q = np.quantile(col, np.linspace(0, 1, 6))
        fail_tbl[nm] = [
            {"lo": float(q[b]), "hi": float(q[b + 1]),
             "failure_rate": float(1.0 - ok[(col >= q[b]) & (col <= q[b + 1])].mean())}
            for b in range(5)]

    summary = {
        "prereg": str(PREREG.relative_to(ROOT)),
        "n_base": a.n_base, "per_stratum": per_stratum,
        "strata": [f"n_P={a_},n_R={b_}" for a_, b_ in strata],
        "total_evaluations": total,
        "units": "times in 1/delta_P; ERK in ID3 units. No clock (bench item 9).",
        "P6_check_hill_exponents_absent_from_indices": all(
            "n_P" not in d.get("S1", {}) and "n_R" not in d.get("S1", {})
            for per in results.values() for d in per.values()),
        "P2_STOP_CONDITION_FIRED": bool(p2_hits),
        "P2_hits": p2_hits,
        "failure_rate_overall": float(1.0 - ok.mean()),
        "failure_rate_vs_parameter": fail_tbl,
        "results": results,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    prov.stamp_run(OUT, seed=spec["method"]["seed"],
                   params={"prereg": str(PREREG.name)},
                   note="Phase 2 Sobol, stratified over the Hill exponents")

    print("\nP6 (exponents absent from indices):",
          summary["P6_check_hill_exponents_absent_from_indices"])
    print("P2 STOP CONDITION FIRED:", summary["P2_STOP_CONDITION_FIRED"])
    if p2_hits:
        print("  ", p2_hits[:4])
    print(f"overall failure rate {summary['failure_rate_overall']:.3%}")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
