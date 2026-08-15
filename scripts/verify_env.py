"""Verify the scientific stack actually imports, not just that pip claimed success.

Why this exists as a committed script rather than a one-off command: pip
reporting "Successfully installed casadi" and `import casadi` working are
different facts on Windows. CasADi ships compiled binaries and can install
cleanly then fail at import on a missing runtime DLL. This project runs for
months across possible machine changes, so the check needs to be repeatable.

Run:  venv\\Scripts\\python.exe scripts\\verify_env.py
Exit code is 1 if anything failed, so this can gate a longer job.
"""

import importlib
import importlib.metadata
import platform
import sys

# (import name, pip name) — they differ for SALib/PyYAML, which is exactly the
# kind of mismatch that makes a naive check pass while the import fails.
PACKAGES = [
    ("numpy", "numpy"),
    ("scipy", "scipy"),
    ("matplotlib", "matplotlib"),
    ("pandas", "pandas"),
    ("sympy", "sympy"),
    ("casadi", "casadi"),
    ("SALib", "SALib"),
    ("h5py", "h5py"),
    ("anndata", "anndata"),
    ("scanpy", "scanpy"),
    ("yaml", "pyyaml"),
]


def main() -> int:
    print(f"python   {sys.version.split()[0]}  ({platform.system()} {platform.release()})")
    print(f"exe      {sys.executable}\n")

    failures = []
    for import_name, pip_name in PACKAGES:
        try:
            module = importlib.import_module(import_name)
        except Exception as exc:  # noqa: BLE001 — we want the reason, whatever it is
            failures.append((pip_name, exc))
            print(f"  FAIL   {pip_name:<12} {type(exc).__name__}: {exc}")
            continue
        # Ask the package metadata, not module.__version__: anndata and scanpy
        # now deprecate the attribute, and SALib never defined one.
        del module
        try:
            version = importlib.metadata.version(pip_name)
        except importlib.metadata.PackageNotFoundError:
            version = "unknown"
        print(f"  ok     {pip_name:<12} {version}")

    if failures:
        print(f"\n{len(failures)} import(s) failed. Do not work around this silently —")
        print("a substituted library changes the numerics and nobody will remember why.")
        return 1

    print("\nAll imports succeeded.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
