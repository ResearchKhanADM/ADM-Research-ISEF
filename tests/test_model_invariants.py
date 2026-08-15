"""Invariants that fail silently if broken — which is why they are tested.

None of these check that the model is *right*. They check that it has not been
quietly turned into a different model. Every one corresponds to a specific way
this project could produce a clean-looking figure that means the opposite of
what it claims.

Run:  venv\\Scripts\\python.exe -m pytest tests -v
"""

from __future__ import annotations

import numpy as np
import pytest

from src import binding, functional_forms as ff
from src.topology import CORE_STATES, Topology, payload_subsets


# ---------------------------------------------------------------------------
# The sign that carries the mechanism
# ---------------------------------------------------------------------------


def test_raf_drive_is_strictly_decreasing():
    """RAF drive MUST fall as ERK rises.

    Phospho-MEK accumulates under a catalytic inhibitor *because* falling ERK
    relieves ERK-mediated negative feedback on RAF. An increasing implementation
    inverts the mechanism, predicts the opposite withdrawal asymmetry, and runs
    perfectly cleanly. This test is the only thing standing between that bug and
    a figure on the board.
    """
    k_eff = np.linspace(0.01, 10.0, 200)
    drive = ff.raf_drive(k_eff, k_fb=1.0, r=2.0)
    assert np.all(np.diff(drive) < 0), "RAF_drive must be strictly decreasing in K_eff"
    assert drive[0] > drive[-1]


def test_inverted_raf_drive_would_be_caught():
    """Guard on the guard: confirm the check above actually discriminates.

    A monotonicity assertion that passes for both an increasing and a decreasing
    function tests nothing. This constructs the inverted implementation and
    requires it to fail.
    """
    k_eff = np.linspace(0.01, 10.0, 200)
    inverted = ff.hill_activate(k_eff, 1.0, 2.0)   # increasing — the bug
    assert not np.all(np.diff(inverted) < 0)


def test_ignition_is_suppressed_by_erk():
    """Cut #1: ERK shuts down the PTF1-independent ignition promoter."""
    k_eff = np.linspace(0.01, 10.0, 100)
    g = ff.g_ignition(k_eff, k_ign=1.0, q=2.0)
    assert np.all(np.diff(g) < 0)


# ---------------------------------------------------------------------------
# Structural constraints from master plan §3.4 — the four non-negotiables
# ---------------------------------------------------------------------------


def test_pd325901_recovers_when_f_act_disabled():
    """f_act ≡ 1 must recover PD325901 exactly, or the Collins comparison dies."""
    for v in (0.0, 0.5, 5.0):
        assert ff.f_act(v, ic50_act=1.0, is_trametinib=False) == 1.0
    # ...and trametinib must actually do something, or the two drugs are one.
    assert ff.f_act(1.0, ic50_act=1.0, is_trametinib=True) < 1.0


def test_two_target_titration_traps_both_species():
    """ID3 must reduce free PTF1A *and* free E — Dufresne 2010, PMID 20830706.

    A titrator acting on only one partner is a different (weaker) mechanism and
    would not give the sharp threshold the bistability argument needs.
    """
    kd = {"IE": 1.0, "IP": 1.0, "L": 1.0, "J": 10.0}
    lo = binding.solve_binding(1.0, 1.0, 0.01, 1.0, kd, rbpj=0.1)
    hi = binding.solve_binding(1.0, 1.0, 5.00, 1.0, kd, rbpj=0.1)
    assert hi["P_free"] < lo["P_free"], "ID3 must sequester PTF1A"
    assert hi["E_free"] < lo["E_free"], "ID3 must sequester E-protein"
    assert hi["C_L"] < lo["C_L"], "and therefore suppress the PTF1-L complex"


def test_binding_conserves_mass():
    """Every total pool must be accounted for. A conservation leak here would
    look like spontaneous protein synthesis and drift the whole trajectory."""
    kd = {"IE": 2.0, "IP": 3.0, "L": 1.5, "J": 8.0}
    p_tot, e_tot, i_tot, r_tot, rbpj = 1.3, 0.9, 0.7, 1.1, 0.2
    b = binding.solve_binding(p_tot, e_tot, i_tot, r_tot, kd, rbpj)
    assert b["P_free"] + b["ID3_P"] + b["C_L"] + b["C_J"] == pytest.approx(p_tot, rel=1e-6)
    assert b["E_free"] + b["ID3_E"] + b["C_L"] + b["C_J"] == pytest.approx(e_tot, rel=1e-6)
    assert b["I_free"] + b["ID3_E"] + b["ID3_P"] == pytest.approx(i_tot, rel=1e-6)
    assert b["R_free"] + b["C_L"] == pytest.approx(r_tot, rel=1e-6)


def _effective_hill(fn, kd, i_grid):
    """n_eff = max |d ln C_L / d ln I| — the standard measure of ultrasensitivity.

    NOT a global fold-change across a fixed ID3 range. Titration produces a
    threshold-linear response: it is very steep NEAR the threshold and flat
    away from it, so a global fold-change averages the steep part away and can
    rank a first-order sink as 'sharper'. That mistake was made here once; this
    docstring exists so it is not made again.
    """
    c_l = np.array([fn(1.0, 1.0, i, 1.0, kd, 0.1)["C_L"] for i in i_grid])
    return np.max(np.abs(np.gradient(np.log(np.maximum(c_l, 1e-300)), np.log(i_grid))))


def test_titration_is_ultrasensitive_only_in_the_tight_binding_regime():
    """The T1-vs-T2 discriminator, and the constraint it puts on Stage 2's box.

    Two-target titration beats a first-order sink ONLY when Kd << the protein
    totals. Measured here:

        Kd = 1.0   (loose)   titration n_eff ~ 2.1   vs first-order 1.9   (1.07x)
        Kd = 0.01  (tight)   titration n_eff ~ 13.4  vs first-order 2.0   (6.7x)
        Kd = 0.001 (tight)   titration n_eff ~ 41.7  vs first-order 2.0   (20.8x)

    The first-order sink sits at exactly 2.0 regardless, which is the analytic
    answer: ID3 taxes P and E through two multiplied factors, so the log-log
    slope is 2 everywhere.

    CONSEQUENCE FOR STAGE 2 — this is the reason the test is written this way.
    If the sampling box's Kd prior does not reach into Kd << totals, T1 and T2
    are numerically indistinguishable and the topology competition silently
    fails to discriminate. That would look like "the data cannot tell these
    architectures apart" when the truth is "the box never sampled the regime
    where they differ". There is no measured Kd for ANY ID3 interaction
    (master plan Part 8), so this prior is entirely a choice — and this test
    says the choice determines the answer.
    """
    i_grid = np.logspace(-2, 1.5, 200)

    loose = {"IE": 1.0, "IP": 1.0, "L": 1.0, "J": 10.0}
    n_tit_loose = _effective_hill(binding.solve_binding, loose, i_grid)
    n_fo_loose = _effective_hill(binding.first_order_sequestration, loose, i_grid)
    assert n_tit_loose / n_fo_loose < 1.5, (
        "at loose binding the two mechanisms should be nearly indistinguishable"
    )

    tight = {"IE": 0.01, "IP": 0.01, "L": 1.0, "J": 10.0}
    n_tit_tight = _effective_hill(binding.solve_binding, tight, i_grid)
    n_fo_tight = _effective_hill(binding.first_order_sequestration, tight, i_grid)
    assert n_tit_tight > 5.0, "titration must be strongly ultrasensitive when tight"
    assert n_tit_tight / n_fo_tight > 3.0, (
        "two-target titration must clearly beat a first-order sink in the "
        "tight-binding regime; if this fails, the bistability argument in §3.2 "
        "needs re-examining before anything downstream is trusted"
    )
    # The first-order sink's analytic slope is exactly 2 (two multiplied factors).
    assert n_fo_tight == pytest.approx(2.0, abs=0.05)


# ---------------------------------------------------------------------------
# Composable topology — decision 003's first reversal condition
# ---------------------------------------------------------------------------


def test_w_is_present_in_every_topology():
    """W is PROTECTED FROM ELIMINATION. If it ever drops out of a state set,
    that topology silently reverts to the static-product model."""
    assert "W" in CORE_STATES
    from src.topology import TOPOLOGIES
    for name, topo in TOPOLOGIES.items():
        assert "W" in topo.states, f"{name} lost W"


def test_dosing_an_absent_species_is_refused():
    """The silent-null trap. Dosing a species the topology does not contain
    completes normally and scores ~0, which reads as 'it doesn't help' rather
    than 'it wasn't there'. Same failure mode flagged for RBPJL in Stage 3B."""
    with pytest.raises(ValueError, match="requires mist1_arm"):
        Topology("bad", u3_identity="MIST1", mist1_arm=False)
    with pytest.raises(ValueError, match="requires nr5a2_mode"):
        Topology("bad", u3_identity="NR5A2", nr5a2_mode="absent")


def test_undeclared_configuration_is_refused():
    with pytest.raises(ValueError):
        Topology("bad", id3_mode="magic")


def test_necessity_analysis_has_all_16_subsets():
    subsets = payload_subsets()
    assert len(subsets) == 16
    assert len({tuple(sorted(s.items())) for s in subsets}) == 16, "subsets must be distinct"


# ---------------------------------------------------------------------------
# Pulse forcing — Stage 5's axis and Stage 6's precondition
# ---------------------------------------------------------------------------


def test_pulse_delivers_the_stated_dose():
    """Each pulse must integrate to `dose`, or the matched-total-dose arm of the
    necessity analysis is confounded: changing a half-life would silently change
    how much was delivered."""
    dose, interval = 2.0, 1e9          # one effectively isolated pulse
    t = np.linspace(0, 400, 400_000)
    u = np.array([ff.mrna_pulse(ti, dose, interval, 1.0, 0.05, 1) for ti in t])
    assert np.trapezoid(u, t) == pytest.approx(dose, rel=1e-3)


def test_pulse_is_nonnegative():
    """Inputs are bounded and non-negative by construction: you can add mRNA,
    you cannot subtract it. Stage 5's control problem assumes this."""
    t = np.linspace(0, 200, 5000)
    u = np.array([ff.mrna_pulse(ti, 1.0, 24.0, 1.0, 0.05, 5) for ti in t])
    assert np.all(u >= -1e-12)


def test_pulse_rejects_inverted_timescales():
    with pytest.raises(ValueError, match="decay slower than translation"):
        ff.mrna_pulse(1.0, dose=1.0, interval=10.0,
                      k_translate=0.01, delta_m=1.0, n_pulses=1)


# ---------------------------------------------------------------------------
# Viability
# ---------------------------------------------------------------------------


def test_hazard_is_u_shaped():
    """Rises at BOTH ends: high S (cargo outruns capacity) and low free PTF1A
    (CHOP-dependent apoptosis). A monotone hazard is the U_crit ceiling again."""
    p = dict(h_max_cargo=1.0, s_crit=1.0, nu_s=4.0,
             h_max_chop=1.0, p_crit=1.0, mu_p=4.0)
    p_mid, s_mid = 1e3, 0.0            # healthy: high PTF1A, no cargo backlog
    assert ff.hazard(5.0, p_mid, p) > ff.hazard(0.0, p_mid, p), "high S must kill"
    assert ff.hazard(s_mid, 1e-3, p) > ff.hazard(s_mid, p_mid, p), "low PTF1A must kill"
