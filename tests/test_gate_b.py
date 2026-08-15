"""Gate B's structural claims, asserted against what the stage actually wrote.

Gate B is *"two stable states plus a saddle, with an identifiable separatrix."*
That is a claim about structure, and a claim that goes in a report should be
checkable from the artefacts rather than from a print statement someone
remembers seeing.

These tests read `results/gate_b/` and **skip** if it is absent — `results/` is
gitignored, so a fresh clone has no run to check. Skipping rather than failing is
deliberate: a red suite on a clean checkout trains people to ignore red suites.
Run `venv\\Scripts\\python.exe scripts/run_gate_b.py` first.

Nothing here asserts a *value*. The window edges depend on placeholder
parameters; only the structure is being tested.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
GATE_B = ROOT / "results" / "gate_b"

pytestmark = pytest.mark.skipif(
    not (GATE_B / "summary.json").exists(),
    reason="no Gate B run present; run scripts/run_gate_b.py",
)


@pytest.fixture(scope="module")
def summary():
    return json.loads((GATE_B / "summary.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def branch():
    return np.genfromtxt(GATE_B / "branch.csv", delimiter=",", names=True)


@pytest.fixture(scope="module")
def labels():
    lines = (GATE_B / "branch_labels.csv").read_text(encoding="utf-8")
    return np.array([ln.split(",")[1] for ln in lines.strip().splitlines()[1:]])


# ---------------------------------------------------------------------------
# The gate itself
# ---------------------------------------------------------------------------


def test_two_stable_states_plus_a_saddle(labels):
    """The literal wording of Gate B."""
    assert "saddle" in set(labels), "no saddle on the branch — Gate B fails"
    assert (labels == "stable").sum() > 0
    assert (labels == "saddle").sum() > 0


def test_exactly_two_folds_bracket_the_window(summary):
    """A bistable *window*, not a half-line.

    Two folds means the bistable region is bounded on both sides: below it
    trametinib alone reverts, above it KRAS holds metaplasia regardless. One fold
    would mean the reverted state persists at arbitrarily high oncogenic drive,
    which would be a different — and much less believable — claim.
    """
    assert len(summary["folds_erk"]) == 2
    lo, hi = summary["persistence_window_erk"]
    assert lo < hi
    assert summary["persistence_window_width_erk"] > 0


def test_the_branch_is_non_monotone_in_erk(branch):
    """The S-curve. A monotone branch means continuation missed the middle sheet,
    which is the failure mode a naive parameter sweep has."""
    erk = branch["erk"]
    assert np.any(np.diff(erk) < 0) and np.any(np.diff(erk) > 0)


def test_the_saddle_has_a_two_dimensional_stable_manifold(summary, branch, labels):
    """Why "separatrix" is a well-defined word here.

    In a 3-D state space a 2-D stable manifold is a *surface*, and a surface can
    separate two basins. Had the saddle carried a 1-D stable manifold instead,
    the object would be a curve and could not divide a 3-D space — "the
    separatrix" would then be the wrong noun, and Gate B's wording would need
    changing rather than the answer being reported.
    """
    assert summary["saddle_stable_manifold_dim"] == 2
    dims = branch["stable_manifold_dim"][labels == "saddle"]
    assert np.all(dims == 2)


def test_stable_branches_have_full_dimensional_stable_manifolds(branch, labels):
    dims = branch["stable_manifold_dim"][labels == "stable"]
    assert np.all(dims == 3)


def test_the_two_stable_branches_are_the_acinar_and_metaplastic_states(
        branch, labels, summary):
    """The high-RBPJL sheet and the collapsed one, both inside the window.

    `R ~ 0` on the metaplastic branch is the bootstrap claim showing up in the
    equilibrium structure — the same signature as the pre-rewrite golden fixture,
    across a complete change of state space.
    """
    lo, hi = summary["persistence_window_erk"]
    mid = 0.5 * (lo + hi)
    inside = (np.abs(branch["erk"] - mid) < 0.05) & (labels == "stable")
    R = branch["R"][inside]
    assert R.max() > 0.5, "no acinar branch inside the window"
    assert R.min() < 0.1, "no collapsed branch inside the window"


# ---------------------------------------------------------------------------
# The separatrix artefact
# ---------------------------------------------------------------------------


def test_separatrix_stays_in_the_physical_orthant():
    """Negative concentrations are an artefact of the clamp in `core.rhs`.

    Below zero the integrator follows clamped dynamics, which are not a smooth
    extension of the model, so a manifold traced there is a numerical fiction. It
    looked entirely plausible before this was bounded — RBPJL reached -7.8.
    """
    sep = np.genfromtxt(GATE_B / "separatrix.csv", delimiter=",", names=True)
    for k in ("P", "R", "C"):
        assert sep[k].min() > -1e-2, f"{k} went to {sep[k].min():.3g}"


def test_separatrix_rays_all_terminated_at_the_domain_edge():
    """Every ray should reach the boundary, not stop early or run away.

    A ray that stopped without leaving the domain either failed or hit the time
    limit, and a manifold assembled from short rays understates the basin
    boundary while looking fine.
    """
    conv = json.loads((GATE_B / "convergence.json").read_text(encoding="utf-8"))
    outcomes = conv["separatrix_ray_outcomes"]
    assert outcomes, "no rays recorded"
    assert all(o["ok"] for o in outcomes)
    assert all(o["left_domain"] for o in outcomes)


# ---------------------------------------------------------------------------
# The standing rule
# ---------------------------------------------------------------------------


def test_convergence_is_reported_against_the_parameter_not_as_a_scalar():
    """The standing rule, checked on the artefact rather than trusted.

    Continuation failures concentrate at folds, and the folds are the reported
    numbers. A scalar rate would average away the only failures that could
    change the answer, so the binned form must be present in the output.
    """
    conv = json.loads((GATE_B / "convergence.json").read_text(encoding="utf-8"))
    assert "failure_rate_vs_erk" in conv
    binned = conv["failure_rate_vs_erk"]
    assert len(binned) >= 5
    assert all({"param_lo", "param_hi", "n", "failure_rate"} <= set(b)
               for b in binned)
    assert sum(b["n"] for b in binned) == conv["corrector_calls"]


def test_failure_rate_is_below_the_stated_bound():
    conv = json.loads((GATE_B / "convergence.json").read_text(encoding="utf-8"))
    worst = max((b["failure_rate"] for b in conv["failure_rate_vs_erk"]),
                default=0.0)
    assert worst < 0.05, (
        f"worst per-bin corrector failure rate {worst:.1%}; results would be "
        f"conditional on convergence and must be reported that way")


def test_the_run_is_stamped_and_labelled_as_placeholder(summary):
    """Provenance, and the label that stops a placeholder number being quoted."""
    assert (GATE_B / "_run.json").exists()
    assert "PLACEHOLDER" in summary["parameters"].upper()
