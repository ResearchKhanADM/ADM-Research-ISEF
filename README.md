# ADM-Research — reversing KRAS-driven acinar-to-ductal metaplasia in silico

Acinar-to-ductal metaplasia (ADM) is the earliest and, in principle, still
reversible step on the path from a healthy pancreas to pancreatic cancer.
Oncogenic KRAS drives acinar cells out of their identity by opening a
self-sustaining transcriptional loop built on PTF1A and its partner RBPJL.
RBPJL has no PTF1A-independent promoter, so once the loop is open it cannot
restart itself — which is why lifting the KRAS signal alone reverses ADM only
slowly.

This repository holds a dynamical-systems model of that loop. It asks a
question with two sides: what combination of MEK inhibition and delivered mRNA
returns a cell to the acinar state, and does so **without killing it**. The
model is reduced to 5–6 slow states, certified for bistability, compared across
five candidate architectures, and validated by predicting published
observations it was never fit to.

The output is an experimental design: which doses, at what interval, in what
order. A wet-lab experiment in AR42J rat acinar cells runs in the same cycle and
gets one shot, so the modeling exists to aim it.

Plan: [`docs/ADM_INSILICO_MASTER_PLAN.md`](docs/ADM_INSILICO_MASTER_PLAN.md).
Running log: [`docs/PROGRESS.md`](docs/PROGRESS.md).

## Setup

```
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
```

`data/`, `results/`, and `logs/` are gitignored and are not created by a clone —
they are recreated locally on first run. Datasets are fetched by accession with
`scripts/download_data.py`; the script is the artifact, never the data.

## Layout

| Path | Contents |
|---|---|
| `src/` | model and analysis code |
| `scripts/` | runnable entry points (sweeps, downloads) |
| `prereg/` | parameter ranges and predictions, committed *before* any sweep runs |
| `figures/` | figures |
| `docs/decisions/` | one file per architectural decision, with what would reverse it |
