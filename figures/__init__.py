"""Figure package.

Figure modules live here as `figNN_<topic>.py`, one per figure, each exposing
`INPUTS` (the `results/` files it reads) and `build(profile="draft")`.

**Figure modules never compute science.** They load `results/`, do display
arithmetic, and draw. See `make_figures.py` and CLAUDE.md for why.
"""
