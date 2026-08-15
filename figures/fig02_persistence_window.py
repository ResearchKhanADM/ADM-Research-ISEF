"""fig02 — the drug-free persistence window. ★ the R3 deliverable.

**This module computes no science.** It reads `results/gate_b/`, does display
arithmetic (picking branch segments, shading between two numbers), and draws. The
folds, the stability classes and the separatrix were all computed by
`scripts/run_gate_b.py` and are read from CSV.

Panel A is the bifurcation diagram: RBPJL at equilibrium against the ERK input,
with both folds marked and the bistable interval between them shaded. **That
shaded interval is the answer**, not decoration — it is the range of oncogenic
drive over which a reverted cell stays reverted with no drug present.

Panel B is the separatrix: the saddle's 2-D stable manifold, the surface that
decides which basin a cell falls into. In the old plan the hysteresis wedge was a
decorative figure with no job; here it is the deliverable, and Panel B is what
makes "the reverted state is a real attractor with a real basin" visible rather
than asserted.

**This is NOT a two-parameter bifurcation diagram and must not be captioned as
one.** For equilibria the (KRAS × trametinib) plane is degenerate: both act
through a single scalar, so iso-contours are hyperbolae and the second axis
carries no independent information. A genuine two-parameter figure does exist —
the withdrawal protocol has independent drug-on and drug-off levels — but it is a
**two-protocol operating window**. Presenting the degenerate version as the money
figure would be dismantled by anyone who noticed.

**NO TIME AXIS APPEARS HERE.** When one does, it is labelled `1/δ_P`: the model
has no clock until Bench Handshake item 9 (PTF1A protein half-life) lands, and no
placeholder conversion to hours is permitted, including in drafts.

CONVENTION NOTE. Project grammar is *solid line = model, open marker = data*.
Within a model, stable and unstable branches are distinguished solid vs dashed —
the universal bifurcation-diagram convention, and a distinction between two kinds
of model output rather than between model and data. Nothing here is data; there
is no data yet.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import _provenance as prov, _style as st

GATE_B = prov.RESULTS / "gate_b"
INPUTS = [GATE_B / "branch.csv", GATE_B / "branch_labels.csv",
          GATE_B / "separatrix.csv", GATE_B / "summary.json"]


def _load():
    import json

    branch = np.genfromtxt(INPUTS[0], delimiter=",", names=True)
    labels = np.array([ln.split(",")[1]
                       for ln in INPUTS[1].read_text(encoding="utf-8")
                       .strip().splitlines()[1:]])
    sep = np.genfromtxt(INPUTS[2], delimiter=",", names=True)
    summary = json.loads(INPUTS[3].read_text(encoding="utf-8"))
    return branch, labels, sep, summary


def build(profile="draft"):
    branch, labels, sep, summary = _load()
    erk, R = branch["erk"], branch["R"]
    lo, hi = summary["persistence_window_erk"]

    with st.house_style():
        fig, (axA, axB) = st.figure(width=9.0, aspect=0.42, ncols=2)

        # ------------------------------------------------ A · bifurcation
        # Shade first so the branches draw over it.
        axA.axvspan(lo, hi, color=st.INTERVENTION, alpha=0.13, lw=0)

        # Split the branch into contiguous runs of one stability class, so a
        # solid segment never bridges a fold. Plotting the whole array at once
        # would draw a line straight across the turn and imply equilibria that
        # do not exist.
        start = 0
        for i in range(1, len(labels) + 1):
            if i == len(labels) or labels[i] != labels[start]:
                seg = slice(start, i + 1 if i < len(labels) else i)
                stable = labels[start] == "stable"
                # Colour by which state the branch represents: the high-RBPJL
                # sheet is acinar, the collapsed one metaplastic.
                colour = (st.ACINAR if np.nanmean(R[seg]) > 0.5
                          else st.METAPLASTIC)
                axA.plot(erk[seg], R[seg],
                         **st.model_kw(color=colour if stable else st.MUTED,
                                       linestyle="-" if stable else (0, (4, 2)),
                                       linewidth=(1.9 if stable else 1.3) * st.SCALE))
                start = i

        for x in (lo, hi):
            axA.plot([x], [np.interp(x, erk[np.argsort(erk)],
                                     R[np.argsort(erk)])],
                     marker="o", ms=4.5 * st.SCALE, color=st.INK,
                     linestyle="none", zorder=5)
            axA.axvline(x, color=st.MUTED, lw=0.7 * st.SCALE,
                        linestyle=(0, (2, 3)), zorder=0)

        axA.set_xlim(0, 1.4)
        # Honest axis: for EQUILIBRIA, KRAS and trametinib act through this one
        # scalar, so a (KRAS x trametinib) plane would be degenerate — iso-
        # contours are hyperbolae and a "two-parameter bifurcation" there is one
        # parameter re-plotted. The genuine second axis is the withdrawal
        # protocol (drug-on vs drug-off levels), which is a two-protocol
        # operating window, not a bifurcation diagram.
        axA.set_xlabel("ERK drive — KRAS and trametinib act only through this scalar")
        axA.set_ylabel("RBPJL at equilibrium")
        axA.set_title("A · reversion persists at zero drug")
        axA.text(0.5 * (lo + hi), axA.get_ylim()[1] * 0.93,
                 f"persistence window\n{lo:.2g} – {hi:.2g}",
                 ha="center", va="top", fontsize=8 * st.SCALE, color=st.INK)

        # ------------------------------------------------ B · separatrix
        # Projection onto (P, R). The manifold is a surface in 3-D; this is a
        # projection of it, said plainly in the axis label rather than implied.
        for k in np.unique(sep["ray"]):
            m = sep["ray"] == k
            axB.plot(sep["P"][m], sep["R"][m],
                     **st.model_kw(color=st.MUTED, linewidth=0.7 * st.SCALE,
                                   alpha=0.75))

        i_sad = np.where(labels == "saddle")[0]
        near = i_sad[np.argmin(np.abs(erk[i_sad] - summary["separatrix_erk"]))]
        axB.plot([branch["P"][near]], [branch["R"][near]], marker="o",
                 ms=6 * st.SCALE, markerfacecolor="white",
                 markeredgecolor=st.INK, markeredgewidth=1.4 * st.SCALE,
                 linestyle="none", zorder=5)
        axB.annotate("saddle", (branch["P"][near], branch["R"][near]),
                     textcoords="offset points", xytext=(8, 6),
                     fontsize=8 * st.SCALE, color=st.INK)

        # One marker per attractor, not every branch sample inside the ERK
        # tolerance — a cluster of near-identical points reads as a data cloud
        # under this project's grammar, which is the one thing a filled marker
        # must never imply here.
        at_erk = np.where((np.abs(erk - summary["separatrix_erk"]) < 0.02)
                          & (labels == "stable"))[0]
        for hi_branch, colour, name in ((True, st.ACINAR, "acinar"),
                                        (False, st.METAPLASTIC, "metaplastic")):
            cand = [i for i in at_erk if (branch["R"][i] > 0.5) == hi_branch]
            if not cand:
                continue
            i = cand[int(np.argmin(np.abs(erk[cand] - summary["separatrix_erk"])))]
            axB.plot([branch["P"][i]], [branch["R"][i]], marker="o",
                     ms=7 * st.SCALE, color=colour, linestyle="none", zorder=6)
            axB.annotate(name, (branch["P"][i], branch["R"][i]),
                         textcoords="offset points", xytext=(0, -14),
                         ha="center", fontsize=8 * st.SCALE, color=colour)

        axB.set_xlabel("PTF1A  (projection of a 3-D surface)")
        axB.set_ylabel("RBPJL")
        axB.set_title(f"B · separatrix at ERK = {summary['separatrix_erk']:.2g}")
        # The 48 rays project onto one curve because the manifold is almost flat
        # in C at this ERK — C moves on 1/gamma, ~50x slower than the proteins.
        # Said in the figure rather than left for a reader to wonder about.
        cspread = float(np.ptp(sep["C"]))
        axB.text(0.97, 0.04, f"48 rays; manifold flat in C (spread {cspread:.0e})",
                 transform=axB.transAxes, ha="right", va="bottom",
                 fontsize=7 * st.SCALE, color=st.MUTED)

        fig.tight_layout()

        source = [{"erk": float(e), "R": float(r), "P": float(pp),
                   "C": float(c), "label": str(l)}
                  for e, r, pp, c, l in zip(erk, R, branch["P"], branch["C"],
                                            labels)]
        prov.save_figure(
            fig, "fig02_persistence_window", inputs=INPUTS, profile=profile,
            caption=(
                "Drug-free persistence window. Equilibrium RBPJL vs ERK drive; "
                "solid = stable, dashed = saddle. Shaded interval is bistable — "
                "a reverted cell stays reverted with no drug. Panel B projects "
                "the saddle's 2-D stable manifold. NOT a two-parameter "
                "bifurcation diagram: KRAS and trametinib act through one scalar. "
                "PLACEHOLDER PARAMETERS: the structure is the result, the numbers "
                "are not."),
            source=source,
        )
    return True
