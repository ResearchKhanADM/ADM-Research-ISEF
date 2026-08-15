"""Build every figure. The single command.

    python make_figures.py                 # all figures, draft profile -> build/
    python make_figures.py fig02 fig03     # just these
    python make_figures.py --profile paper # -> figures/out/, refuses a dirty tree
    python make_figures.py --check         # render nothing; fail if anything is stale
    python make_figures.py --selftest      # prove the style/provenance machinery works
    python make_figures.py --list          # what modules exist

Each figure lives in `figures/figNN_<topic>.py` and exposes:

    INPUTS = [...]        # results/ files it reads (may be empty for a schematic)
    def build(profile="draft"): ...

**THE LOAD-BEARING RULE: a figure module never computes science.** It loads
`results/`, does display arithmetic, and draws. If a figure needs a number that
is not in `results/`, the fix is a stage that writes it — not a solve in the
figure. That is what keeps this command under a minute, keeps figures from
drifting out of step with the analysis, and stops a slow stage from blocking a
figure. A figure module that imports `src.model` is a bug.

Missing inputs are a SKIP, never a crash: early in the project most stages have
not run, and a build that dies on the first absent result file is a build nobody
runs. Skips are reported explicitly so "it built fine" cannot mean "it built
nothing".
"""

from __future__ import annotations

import argparse
import importlib
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from figures import _provenance as prov  # noqa: E402


def discover() -> list[tuple[str, str]]:
    """(slug, module name) for every figNN_*.py, in filename order."""
    out = []
    for p in sorted((ROOT / "figures").glob("fig*_*.py")):
        out.append((p.stem.split("_")[0], f"figures.{p.stem}"))
    return out


def _missing(mod) -> list[str]:
    return [str(p) for p in getattr(mod, "INPUTS", []) if not Path(p).exists()]


def cmd_list() -> int:
    mods = discover()
    if not mods:
        print("no figure modules yet — figures/figNN_<topic>.py")
        return 0
    for slug, name in mods:
        mod = importlib.import_module(name)
        ins = getattr(mod, "INPUTS", [])
        print(f"{slug:8s} {name:34s} inputs: {len(ins)}"
              f"{'  [MISSING]' if _missing(mod) else ''}")
    return 0


def cmd_check(profile: str) -> int:
    """Fail if any figure is stale. This is the pre-poster gate."""
    mods = discover()
    if not mods:
        print("nothing to check — no figure modules yet")
        return 0
    bad = 0
    for slug, _ in mods:
        stale, why = prov.is_stale(slug, profile=profile)
        print(f"{'STALE' if stale else 'ok   '}  {slug:8s} {why}")
        bad += stale
    if bad:
        print(f"\n{bad} figure(s) stale at profile {profile!r}. "
              f"Rebuild before this goes anywhere.")
    return 1 if bad else 0


def cmd_build(only, profile: str) -> int:
    mods = [m for m in discover() if not only or m[0] in only]
    if not mods:
        print("no figure modules to build yet.")
        print("Infrastructure is in place; add figures/figNN_<topic>.py as "
              "stages start writing results/.")
        return 0

    built = skipped = failed = 0
    for slug, name in mods:
        mod = importlib.import_module(name)
        missing = _missing(mod)
        if missing:
            # Expected for most of the project. Say exactly what is absent so
            # this reads as "the stage has not run yet", not as a failure.
            print(f"SKIP  {slug}: missing {len(missing)} input(s) -> {missing[0]}")
            skipped += 1
            continue
        try:
            mod.build(profile=profile)
            print(f"BUILT {slug}  [{profile}]")
            built += 1
        except prov.DirtyTreeError as e:
            print(f"FAIL  {slug}: {e}")
            failed += 1
        except Exception:
            print(f"FAIL  {slug}:")
            traceback.print_exc()
            failed += 1

    print(f"\n{built} built, {skipped} skipped, {failed} failed.")
    return 1 if failed else 0


def cmd_selftest() -> int:
    """Render a style proof sheet — no project data, no science.

    This exists so the style and provenance machinery can be verified before any
    stage has produced a result. It draws the four nominal colours and the
    model-vs-data grammar on synthetic numbers, and writes a real provenance
    sidecar, exercising the whole path end to end.
    """
    import numpy as np
    from figures import _style as st

    with st.house_style():
        fig, ax = st.figure(width=5.4)
        x = np.linspace(0, 10, 200)
        for i, (name, c) in enumerate(
                (("acinar", st.ACINAR), ("metaplastic", st.METAPLASTIC),
                 ("intervention", st.INTERVENTION), ("toxic", st.TOXIC))):
            ax.plot(x, np.sin(x - 0.5 * i) - i, label=f"{name} (model)",
                    **st.model_kw(color=c))
        xd = np.arange(1, 10, 1.5)
        ax.plot(xd, np.sin(xd) + 0.35, label="data (open marker, no line)",
                **st.data_kw())
        ax.set_xlabel("x (arbitrary)")
        ax.set_ylabel("y (arbitrary)")
        ax.set_title("style proof sheet — synthetic data, no science")
        ax.legend(loc="lower left", ncol=1)
        st.label_panel(ax, "A")
        prov.save_figure(
            fig, "fig00_style_selftest", inputs=(), profile="draft",
            caption="Style/provenance self-test. Synthetic numbers only.",
            source=[{"x": float(a), "y": float(b)} for a, b in zip(xd, np.sin(xd) + 0.35)],
        )

    out = prov.BUILD / "fig00_style_selftest.pdf"
    print(f"selftest OK -> {out}")
    print(f"             {prov.BUILD / 'fig00_style_selftest.prov.json'}")
    print(f"             {prov.BUILD / 'fig00_style_selftest_source.csv'}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("only", nargs="*", help="figure slugs, e.g. fig02 fig03")
    ap.add_argument("--profile", default="draft", choices=sorted(prov.PROFILES),
                    help="draft -> build/figures/; paper/poster -> figures/out/ "
                         "and refuse a dirty tree")
    ap.add_argument("--check", action="store_true",
                    help="render nothing; exit 1 if any figure is stale")
    ap.add_argument("--selftest", action="store_true",
                    help="render a style proof sheet (no project data)")
    ap.add_argument("--list", action="store_true", help="list figure modules")
    a = ap.parse_args()

    if a.list:
        return cmd_list()
    if a.selftest:
        return cmd_selftest()
    if a.check:
        return cmd_check("paper" if a.profile == "draft" else a.profile)
    return cmd_build(set(a.only), a.profile)


if __name__ == "__main__":
    raise SystemExit(main())
