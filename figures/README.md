# Figures

General-purpose geometry. **No IEEE or conference-template sizing** — one `SCALE`
constant in `_style.py`, scaled downstream.

```bash
python make_figures.py            # all, draft profile -> build/figures/
python make_figures.py fig02      # one
python make_figures.py --profile paper   # -> figures/out/, refuses a dirty tree
python make_figures.py --check    # exit 1 if anything is stale
python make_figures.py --selftest # prove the machinery works, no project data
```

## The rule

**A figure module never computes science.** It loads `results/`, does display
arithmetic, and draws. It never runs an ODE solve or a sample. If a figure needs
a number that is not in `results/`, the fix is a stage that writes it.

Consequences: rebuild is always under a minute, figures cannot drift from the
analysis, and a slow stage never blocks a figure. A figure module that imports
`src.model` is a bug.

## Writing one

```python
# figures/fig03_formulation.py
from pathlib import Path
from . import _style as st, _provenance as prov

INPUTS = [prov.RESULTS / "phase3_delivery" / "converted_fraction.csv"]

def build(profile="draft"):
    import pandas as pd
    df = pd.read_csv(INPUTS[0])            # load — never compute
    with st.house_style():
        fig, ax = st.figure()
        ax.plot(df.cv, df.coformulated, **st.model_kw(st.INTERVENTION))
        ax.plot(df.cv, df.separate,      **st.model_kw(st.TOXIC))
        ax.set_xlabel("per-cell uptake CV")
        ax.set_ylabel("fraction of cells above threshold for both mRNAs")
        prov.save_figure(fig, "fig03_formulation", inputs=INPUTS,
                         profile=profile, source=df)
```

`INPUTS` is what makes `--check` work: `save_figure()` hashes those files, and
`is_stale()` compares the hashes later. Missing inputs are a **skip**, not a
crash — early on, most stages have not run.

## Grammar, fixed project-wide

- **A solid line is a model.** `st.model_kw(...)`
- **An open marker with a dark edge is data.** `st.data_kw(...)` — and
  **experimental points are never connected by a line.**
- Ensembles are drawn as thin translucent members (`st.ensemble_kw`), not as a
  single representative curve.
- Four nominal colours only: `#0072B2` acinar · `#D55E00` metaplastic ·
  `#009E73` intervention/success · `#6E6E6E` toxic — **grey, not red**.
- Always `with st.house_style():`, never a bare `rcParams` mutation — one process
  imports every module.

## Core figures (v3 Part 5)

| Slug | Phase | Shows |
|---|---|---|
| `fig01_loop_schematic` | — | the reduced 3–4 node core |
| `fig02_persistence_window` ★ | 2 | (KRAS × trametinib), bistable wedge, "reversion persists at zero drug" |
| `fig03_formulation` ★ | 3 | double-above-threshold fraction vs uptake CV, co-formulated vs separate |
| `fig04_marginal_value` ★ | 4 | durable reversal vs payload size at fixed mass; knee marked; ensemble drawn |
| `fig05_durability` | 5 | time-to-relapse after clearance; dose × interval three-region map |
| `fig06_heldout` | 6 | prediction interval vs published numbers; dated box **left of** outcome |
| `fig07_prioritization` | 1 | ranked regulon activity across ≥2 datasets, labelled positive control |
| `figS0x_*` | — | profile likelihood, QSSA error, CellOracle negative control |

## Traps

No 3-D cusp surface · no hairball of all states · no continuous colormap under a
three-class map · no fake error bars on deterministic output · no p-value where
there is no sampling · no truncated percentage axis · no bimodality claimed from
a KDE shape · two significant figures unless you can defend more.
