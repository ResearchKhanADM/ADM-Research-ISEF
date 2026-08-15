"""Provenance: what produced this figure, from what, at what commit.

The problem this solves is specific. Five months from now there will be a PDF on
a board and a question — *"is this the current model?"* — and the honest answer
has to be checkable rather than remembered. So:

  * every **stage** that writes to `results/` calls `stamp_run()`, recording the
    git commit, whether the tree was dirty, the seed, library versions, and a
    SHA-256 of every output file;
  * every **figure** declares the result files it read, and `save_figure()`
    hashes them and writes `<slug>.prov.json` beside the PDF;
  * `is_stale()` compares the recorded hashes against the files on disk, so
    `make_figures.py --check` can fail loudly instead of shipping a stale figure.

Why hashes and not timestamps: a `--resume`d sweep rewrites files with new
mtimes and identical contents, and a figure rebuilt from a partially-overwritten
checkpoint has a *newer* mtime than the data it disagrees with. Content hashes
do not have either failure mode.

**The `paper` and `poster` profiles refuse to render from a dirty tree.** That
is deliberately annoying: a figure whose provenance says "commit abc123" but
which was actually built from uncommitted edits is worse than no provenance,
because it is confidently wrong. Draft rendering to `build/` is unrestricted.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIG_OUT = ROOT / "figures" / "out"       # TRACKED — final PDFs + provenance
BUILD = ROOT / "build" / "figures"        # GITIGNORED — drafts and proofs

#: Rendering profiles. `draft` is for working; the other two go in front of
#: people and therefore require a clean tree.
PROFILES = {
    "draft": {"require_clean": False, "dir": BUILD},
    "paper": {"require_clean": True, "dir": FIG_OUT},
    "poster": {"require_clean": True, "dir": FIG_OUT},
}


class DirtyTreeError(RuntimeError):
    """Raised when a `paper`/`poster` render is attempted from a dirty tree."""


# ---------------------------------------------------------------------------
# git and environment
# ---------------------------------------------------------------------------

def _git(*args: str) -> str:
    try:
        return subprocess.run(("git", *args), cwd=ROOT, capture_output=True,
                              text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        # Not a git repo, or git unavailable. Record the absence rather than
        # crashing — a figure built outside git is a fact worth writing down.
        return ""


def git_state() -> dict:
    status = _git("status", "--porcelain")
    return {
        "commit": _git("rev-parse", "HEAD") or None,
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD") or None,
        # `dirty` counts *tracked* modifications only. Untracked scratch files
        # are not a provenance problem; edited source is.
        "dirty": bool([ln for ln in status.splitlines() if not ln.startswith("??")]),
    }


def env_state() -> dict:
    """Versions of the libraries that can change a number or a rendering."""
    versions = {}
    for name in ("numpy", "scipy", "matplotlib", "pandas", "sympy"):
        try:
            from importlib.metadata import version
            versions[name] = version(name)
        except Exception:
            versions[name] = "unknown"
    return {"python": sys.version.split()[0], "platform": platform.platform(),
            "packages": versions}


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# ---------------------------------------------------------------------------
# Stage side — call this from anything that writes results/
# ---------------------------------------------------------------------------

def stamp_run(run_dir, *, seed=None, params=None, outputs=None, note=None) -> Path:
    """Write `_run.json` into a results directory.

    Call at the END of a stage, once its outputs exist, so the hashes describe
    files that are actually complete. A stamp written up front would describe
    intent rather than result.
    """
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    outs = list(outputs) if outputs is not None else [
        p for p in sorted(run_dir.iterdir()) if p.is_file() and p.name != "_run.json"
    ]
    stamp = {
        "written": _now(),
        "git": git_state(),
        "env": env_state(),
        "seed": seed,
        "params": params,
        "note": note,
        "outputs": {str(Path(p).relative_to(ROOT)): sha256(p) for p in outs},
    }
    path = run_dir / "_run.json"
    path.write_text(json.dumps(stamp, indent=2), encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Figure side
# ---------------------------------------------------------------------------

def save_figure(fig, slug, *, inputs=(), profile="draft", source=None,
                caption=None, formats=("pdf",)):
    """Save a figure with its provenance sidecar.

    `inputs`  — every `results/` file this figure read. Declaring nothing is
                allowed (a pure schematic has no inputs) but declaring the wrong
                thing defeats staleness checking, so keep it honest.
    `source`  — the numbers actually plotted, as a list of dict rows or a pandas
                DataFrame. Written to `<slug>_source.csv` so any number on the
                figure can be traced without re-running anything. This is the
                single cheapest credibility move available.
    """
    if profile not in PROFILES:
        raise ValueError(f"unknown profile {profile!r}; expected {sorted(PROFILES)}")
    cfg = PROFILES[profile]
    git = git_state()
    if cfg["require_clean"] and git["dirty"]:
        raise DirtyTreeError(
            f"profile {profile!r} refuses to render from a dirty tree. Commit "
            f"first, or render with profile='draft' into build/figures/. A "
            f"figure stamped with a commit it was not built from is worse than "
            f"an unstamped one."
        )

    outdir = Path(cfg["dir"])
    outdir.mkdir(parents=True, exist_ok=True)

    written = []
    for ext in formats:
        p = outdir / f"{slug}.{ext}"
        fig.savefig(p)             # bbox/pad come from the rc block, not here
        written.append(p)

    src_path = None
    if source is not None:
        src_path = outdir / f"{slug}_source.csv"
        _write_source_csv(src_path, source)

    prov = {
        "slug": slug,
        "profile": profile,
        "written": _now(),
        "git": git,
        "env": env_state(),
        "caption": caption,
        "inputs": {str(Path(p)): (sha256(p) if Path(p).exists() else None)
                   for p in inputs},
        "outputs": {p.name: sha256(p) for p in written},
        "source_csv": src_path.name if src_path else None,
    }
    (outdir / f"{slug}.prov.json").write_text(json.dumps(prov, indent=2),
                                              encoding="utf-8")
    return written


def _write_source_csv(path, source) -> None:
    """Write plotted numbers to CSV. Accepts a DataFrame or a list of dicts."""
    if hasattr(source, "to_csv"):
        source.to_csv(path, index=False)
        return
    import csv
    rows = list(source)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def is_stale(slug, *, profile="paper") -> tuple[bool, str]:
    """Has anything this figure was built from changed since it was rendered?

    Returns `(stale, reason)`. A missing figure counts as stale, because the
    caller wants "is this figure current and present", not "has it drifted".
    """
    outdir = Path(PROFILES[profile]["dir"])
    prov_path = outdir / f"{slug}.prov.json"
    if not prov_path.exists():
        return True, "no provenance record — never rendered at this profile"
    prov = json.loads(prov_path.read_text(encoding="utf-8"))

    for name, recorded in prov.get("outputs", {}).items():
        p = outdir / name
        if not p.exists():
            return True, f"output {name} is missing"
        if sha256(p) != recorded:
            return True, f"output {name} was modified after rendering"

    for path, recorded in prov.get("inputs", {}).items():
        p = Path(path)
        if recorded is None:
            return True, f"input {path} did not exist when rendered"
        if not p.exists():
            return True, f"input {path} has been deleted"
        if sha256(p) != recorded:
            return True, f"input {path} changed since rendering"

    head = git_state()["commit"]
    if head and prov.get("git", {}).get("commit") != head:
        # Not necessarily wrong — the commit may not have touched anything this
        # figure reads — so this is reported as staleness, not as an error.
        return True, "rendered at a different commit"
    return False, "current"


def relpath(p) -> str:
    """Repo-relative path, for readable provenance records."""
    try:
        return str(Path(p).resolve().relative_to(ROOT))
    except ValueError:
        return str(p)


__all__ = [
    "ROOT", "RESULTS", "FIG_OUT", "BUILD", "PROFILES", "DirtyTreeError",
    "git_state", "env_state", "sha256", "stamp_run", "save_figure", "is_stale",
    "relpath",
]
