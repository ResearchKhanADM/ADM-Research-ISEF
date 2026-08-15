# ADM-Research-ISEF — making reversal of acinar-to-ductal metaplasia *stick*

<https://github.com/ResearchKhanADM/ADM-Research-ISEF>

Acinar-to-ductal metaplasia (ADM) is the earliest and still-reversible step from a
healthy pancreas toward pancreatic cancer. **Reversing it is not the open problem:**
MEK inhibition reverts established lesions in about three days (Collins 2014, PMID
24315826) — and the same paper reports that lesions come back once the drug is
withdrawn. The reversal is real and it is drug-dependent.

The reason is structural. Oncogenic KRAS opens a transcriptional loop built on
PTF1A and its partner RBPJL, and *Rbpjl* has **no PTF1A-independent promoter** — so
once the loop is open it cannot re-close on its own. A MEK inhibitor lifts the
KRAS-dependent cuts but cannot rebuild RBPJL, because nothing except the loop makes
RBPJL. That bootstrap failure is what a delivered mRNA payload is for.

This repository holds a small dynamical model of that loop plus a delivery layer,
built to answer three questions the published work leaves open:

1. **Formulation** — per-cell LNP dose is a *distribution*, not a number. For an
   obligate stoichiometric pair that must arrive in the *same cell*, co-formulating
   both mRNAs in one particle should beat separate particles by far more than
   intuition suggests. How much, across the plausible uptake-CV range?
2. **Composition** — how many components, at what ratio, under a fixed total mRNA
   mass budget.
3. **Durability** — the drug-free persistence window: where in dose and parameter
   space does the reverted state survive trametinib withdrawal?

The output is an experimental design — composition, ratio, formulation, dose and
redosing interval. A wet-lab experiment in AR42J rat acinar cells runs in the same
cycle and gets **one shot, no iteration**, so the modelling exists to aim it. It is
framed as computation that designed an experiment, never as validation of one.

Plan of record: [`docs/ADM_MASTER_PLAN_v3.md`](docs/ADM_MASTER_PLAN_v3.md) — read
Part 0 first; it explains why the project changed shape in August 2026.
Running log: [`docs/PROGRESS.md`](docs/PROGRESS.md).
Working agreement: [`CLAUDE.md`](CLAUDE.md).

## Setup

```
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
```

`data/`, `results/`, `logs/` and `build/` are gitignored and are not created by a
clone — they are recreated locally on first run. Datasets are fetched by accession
with a script; the script is the artifact, never the data.

## Layout

| Path | Contents |
|---|---|
| `src/` | model and analysis code |
| `scripts/` | runnable entry points (sweeps, downloads) |
| `prereg/` | parameter ranges and predictions, committed **before** any sweep runs |
| `figures/` | house style, provenance, one module per figure; `out/` holds the final PDFs |
| `docs/decisions/` | one file per architectural decision, each stating what would reverse it |

Figures build with `python make_figures.py`. **Figure modules never compute
science** — they read `results/` and draw, so a rebuild is always under a minute
and a figure cannot drift out of step with the analysis.

## A note on the history

The August 2026 commits contain a larger model (11–13 states, a five-way topology
competition) built against a superseded plan. It was cut, not abandoned quietly:
its discriminating measurement turns out not to exist in the literature. The
reasoning, and the specific finding that would bring it back, are in
[`docs/decisions/012-durability-framing-architecture-change.md`](docs/decisions/012-durability-framing-architecture-change.md).
