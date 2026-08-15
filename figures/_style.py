"""One house style for every figure in the project.

General-purpose geometry. **No IEEE or conference-template sizing** — Luqmaan
scales figures downstream, so there is exactly one knob here (`SCALE`) and
everything is expressed relative to it. A figure that hard-codes a width in
inches to fit some template will be wrong the moment the template changes.

THE RULE THIS MODULE EXISTS TO ENFORCE (CLAUDE.md standing rule):
**figure modules never compute science.** This module provides paint, not data.
If you find yourself importing `src.model` into a figure, stop — the number
belongs in `results/`, written by a stage.

Use `with house_style():` around every render. NEVER mutate `plt.rcParams`
directly: `make_figures.py` imports every figure module into one process, so a
leaked mutation from fig02 silently restyles fig05. That is the classic "it
looked right when I ran it alone" bug, and it is invisible in review.
"""

from __future__ import annotations

from contextlib import contextmanager

import matplotlib as mpl
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# The one knob
# ---------------------------------------------------------------------------

#: Global size multiplier. Everything below is a multiple of this, so changing
#: it resizes fonts, line widths, markers and figure dimensions *together* —
#: which is the only way a figure stays visually consistent when rescaled.
SCALE = 1.0

#: Base figure width in inches at SCALE=1. A neutral single-column-ish default,
#: chosen because it looks right on screen and in a slide; it carries no
#: conference geometry and is expected to be overridden per figure.
BASE_WIDTH = 5.0

#: Golden-ish default aspect. Figures with a natural aspect (phase planes,
#: square heatmaps) should override rather than distort.
BASE_ASPECT = 0.68


# ---------------------------------------------------------------------------
# Colour — Okabe–Ito, capped at four nominal roles
# ---------------------------------------------------------------------------

# Four roles, fixed project-wide, so a colour means the same thing in every
# figure. TOXIC IS GREY, NOT RED, on purpose: green/red is the textbook
# deuteranopia collision and "success vs toxic" is exactly the comparison a
# judge must be able to make at a glance from two metres away.
ACINAR = "#0072B2"        # blue   — the target state
METAPLASTIC = "#D55E00"   # orange — the ADM state
INTERVENTION = "#009E73"  # green  — intervention / success
TOXIC = "#6E6E6E"         # grey   — toxic / failed / excluded

#: Ink for axes, text and data-marker edges. Not pure black — pure black against
#: white is harsher than it needs to be in print.
INK = "#1A1A1A"
#: For de-emphasised annotation (reference lines, shading edges).
MUTED = "#8C8C8C"

#: Ordered cycle for nominal categories. Capped at four *by design*: if a figure
#: needs a fifth colour, the figure is showing too much and should be split.
CATEGORICAL = (ACINAR, METAPLASTIC, INTERVENTION, TOXIC)

#: Named roles, so figure code says what it means rather than naming a hex.
ROLE = {
    "acinar": ACINAR,
    "metaplastic": METAPLASTIC,
    "intervention": INTERVENTION,
    "success": INTERVENTION,
    "toxic": TOXIC,
    "ink": INK,
    "muted": MUTED,
}


# ---------------------------------------------------------------------------
# Global grammar — model vs data
# ---------------------------------------------------------------------------

# Fixed once, project-wide, so a reader learns it from one figure and it holds
# for all of them:
#
#     a SOLID LINE is a MODEL
#     an OPEN MARKER WITH A DARK EDGE is DATA
#
# and experimental points are NEVER connected by a line — a connecting line
# asserts an interpolation the experiment did not measure. Use these helpers
# rather than re-specifying kwargs, so the grammar cannot drift figure to figure.

def model_kw(color=ACINAR, **over):
    """Style for a model curve: solid line, no markers."""
    kw = dict(color=color, linestyle="-", linewidth=1.6 * SCALE, marker="",
              solid_capstyle="round", zorder=2)
    kw.update(over)
    return kw


def data_kw(color=INK, **over):
    """Style for measured points: open marker, dark edge, NO connecting line.

    `linestyle="none"` is not a default to be overridden casually — connecting
    experimental points is the single most common way a figure claims more than
    the experiment measured.
    """
    kw = dict(linestyle="none", marker="o", markersize=5.0 * SCALE,
              markerfacecolor="white", markeredgecolor=color,
              markeredgewidth=1.2 * SCALE, zorder=3)
    kw.update(over)
    return kw


def ensemble_kw(color=ACINAR, **over):
    """Style for one member of a drawn ensemble.

    v3 Part 5: *"Draw the ensemble, not a single curve; a single curve from an
    unidentifiable model is a lie by graphic design."* Thin and translucent, so
    N curves read as a band rather than as N claims.
    """
    kw = dict(color=color, linestyle="-", linewidth=0.6 * SCALE, alpha=0.18,
              marker="", zorder=1)
    kw.update(over)
    return kw


# ---------------------------------------------------------------------------
# The rc block
# ---------------------------------------------------------------------------

def _rc() -> dict:
    s = SCALE
    return {
        # --- text. Sans-serif; DejaVu ships with matplotlib, so a figure built
        # on this machine renders identically on a machine without Helvetica.
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 9.0 * s,
        "axes.titlesize": 9.5 * s,
        "axes.labelsize": 9.0 * s,
        "xtick.labelsize": 8.0 * s,
        "ytick.labelsize": 8.0 * s,
        "legend.fontsize": 8.0 * s,
        # Axis labels are NEVER bold. Bold labels compete with the data for
        # attention and read as shouting in a figure that has real content.
        "axes.labelweight": "normal",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",

        # --- spines and grid. No gridlines, no top/right spines: both add ink
        # that carries no information.
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.linewidth": 0.8 * s,
        "xtick.major.width": 0.8 * s,
        "ytick.major.width": 0.8 * s,
        "xtick.major.size": 3.0 * s,
        "ytick.major.size": 3.0 * s,
        "xtick.direction": "out",
        "ytick.direction": "out",

        # --- lines and colour cycle
        "lines.linewidth": 1.6 * s,
        "lines.markersize": 5.0 * s,
        "axes.prop_cycle": mpl.cycler(color=list(CATEGORICAL)),

        # --- legend: no frame, no shadow
        "legend.frameon": False,
        "legend.handlelength": 1.8,
        "legend.borderpad": 0.2,

        # --- figure
        "figure.figsize": (BASE_WIDTH * s, BASE_WIDTH * BASE_ASPECT * s),
        "figure.dpi": 110,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 300,

        # --- export. fonttype 42 embeds TrueType so text stays selectable and
        # editable in Illustrator; svg.fonttype "none" keeps SVG text as text.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",

        # --- THE ONE THAT BITES. bbox="tight" silently changes the output size
        # to fit whatever happens to be drawn, so two figures asked for the same
        # width come out different widths and every width check is meaningless.
        "savefig.bbox": "standard",
        "savefig.pad_inches": 0.02,
    }


@contextmanager
def house_style(scale: float | None = None):
    """Apply the house style for the duration of one figure.

    Always a context manager, never a bare rcParams mutation — see the module
    docstring for why (one process imports every figure module).
    """
    global SCALE
    prev = SCALE
    if scale is not None:
        SCALE = scale
    try:
        with mpl.rc_context(_rc()):
            yield
    finally:
        SCALE = prev


def figure(width=None, aspect=BASE_ASPECT, **kw):
    """A figure sized off SCALE. Call inside `house_style()`."""
    w = (width if width is not None else BASE_WIDTH) * SCALE
    return plt.subplots(figsize=(w, w * aspect), **kw)


def label_panel(ax, letter, dx=-0.02, dy=1.02):
    """Panel letter in the axes' top-left. Plain, not bold, not parenthesised."""
    ax.text(dx, dy, letter, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=10.0 * SCALE, color=INK)
