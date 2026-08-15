"""Supplementary implementations — NOT part of the Phase 2 core.

Code here is kept because it may be needed later or because it documents a
mechanism the core simplifies. It is not imported by `src/core.py` and nothing in
`results/` currently depends on it.

`binding.py` — the full three-species competitive equilibrium (ID3 against both
E-protein and PTF1A, two ternary complexes). The core titrates one target with a
closed form; this is the only implementation of the two-target mechanism, and
Phase 3 may need it if threshold sharpness turns out to depend on ID3 taxing
PTF1A as well. See `docs/decisions/006-two-target-titration.md`.
"""
