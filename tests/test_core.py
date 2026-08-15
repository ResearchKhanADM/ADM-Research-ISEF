"""Invariants of the Phase 2 core.

None of these check that the model is *right*. They check it has not been quietly
turned into a different model. The bootstrap claim has its own file
(`test_bootstrap_guard.py`), written before the right-hand side.

**No test here asserts a scientific result.** `default_params()` is a placeholder
set; anything that looked like a finding would be a finding about those
placeholders. The qualitative-structure tests at the bottom assert only that the
core *can* do what Gate B will ask of it, and they say so.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from src import core


# ---------------------------------------------------------------------------
# Binding — the exact solution and its two limits
# ---------------------------------------------------------------------------


def test_free_e_recovers_the_linear_tight_binding_limit():
    """kappa -> 0 must give `E_free = max(0, 1 - id3)` — floor included, free.

    This is the whole argument of decision 006's amendment: v3's linear form is
    this expression's tight limit, so the exact solution reproduces it *and*
    reproduces the flooring at `id3 > 1` without a `max(0, .)` hack. If this
    fails, the justification for using the exact form has evaporated.
    """
    for id3 in (0.0, 0.25, 0.5, 0.99, 1.01, 2.0, 10.0):
        got = core.free_e_protein(id3, kappa=1e-10)
        assert got == pytest.approx(max(0.0, 1.0 - id3), abs=1e-5)


def test_free_e_recovers_the_loose_binding_limit():
    """kappa -> inf must give `1/(1 + id3/kappa)` — the first-order form."""
    kappa = 1e8
    for id3 in (0.1, 1.0, 10.0, 100.0):
        assert core.free_e_protein(id3, kappa) == pytest.approx(
            1.0 / (1.0 + id3 / kappa), rel=1e-8)


def test_free_e_never_leaves_the_unit_interval():
    """A free fraction outside [0, 1] is unphysical.

    The linear approximation violates this for `id3 > 1`; that violation is the
    reason the exact form is used. Swept hard across both regimes, because the
    negativity is exactly what a floor hack would have hidden.
    """
    rng = np.random.default_rng(0)
    id3 = 10 ** rng.uniform(-6, 4, size=4000)
    kappa = 10 ** rng.uniform(-8, 8, size=4000)
    e = core.free_e_protein(id3, kappa)
    assert np.all(e >= 0.0) and np.all(e <= 1.0)


def test_free_e_is_monotone_decreasing_in_id3():
    """More titrator, less free partner. An increasing implementation would
    invert the titration mechanism and still run cleanly."""
    id3 = np.logspace(-4, 3, 500)
    for kappa in (1e-4, 1e-2, 1.0, 1e2):
        e = core.free_e_protein(id3, kappa)
        assert np.all(np.diff(e) < 0)


def test_free_e_is_numerically_stable_at_loose_binding():
    """The catastrophic-cancellation guard.

    `(A - sqrt(A^2 - 4*id3))/2` and `2*id3/(A + sqrt(A^2 - 4*id3))` are
    algebraically identical; only the second survives large `kappa`. This checks
    the implementation kept the stable one by comparing against the analytic
    loose-limit expansion, which the unstable form fails by orders of magnitude.
    """
    kappa, id3 = 1e12, 1.0
    assert core.free_e_protein(id3, kappa) == pytest.approx(
        1.0 - id3 / kappa, rel=1e-6)


def test_n_eff_reports_the_binding_regime():
    """The sharpness diagnostic: ~0.5/sqrt(kappa) when tight, 1 when loose.

    Not decoration — Phase 3 convolves the per-cell dose distribution against a
    threshold whose sharpness is this number, so it is an input to the headline
    co-formulation gap.
    """
    assert core.n_eff(1e-2) == pytest.approx(5.0, rel=0.01)
    assert core.n_eff(1e-4) == pytest.approx(50.0, rel=0.01)
    assert core.n_eff(1e3) == pytest.approx(1.0)          # floors at the loose limit


def test_n_eff_agrees_with_the_measured_slope_of_the_exact_solution():
    """Cross-check the diagnostic against the function it actually describes.

    **This test is why the prefactor is 0.5 and not 1.34.** The docs carried
    1.34 from a derivation written for the deleted 11-state model's *ternary
    complex* under *two-target* titration; this core titrates one target and
    takes the slope on `E_free`. Measuring `max |d ln E_free / d ln id3|` on the
    shipped closed form gives 0.5/sqrt(kappa), and the discrepancy was caught
    here rather than in a figure.

    Kept as a permanent tie between formula and implementation: if either drifts,
    Phase 3's threshold sharpness is wrong and nothing downstream would notice,
    because the model would still run and the figure would still render.
    """
    for kappa in (1e-2, 1e-3, 1e-4):
        id3 = np.logspace(-2.5, 2.5, 200_000)
        e = np.maximum(core.free_e_protein(id3, kappa), 1e-300)
        measured = np.max(np.abs(np.gradient(np.log(e), np.log(id3))))
        assert measured == pytest.approx(core.n_eff(kappa), rel=0.02)


# ---------------------------------------------------------------------------
# Right-hand side: shape, signs, inputs
# ---------------------------------------------------------------------------


def test_state_vector_is_three_states():
    """v3 Phase 2 asks for 3-4 states. A silent fourth would mean the reduction
    drifted back toward the model that was just deleted."""
    assert core.STATES == ("P", "R", "C")
    dy = core.rhs(0.0, np.array([1.0, 1.0, 0.5]), core.default_params(),
                  core.Inputs.constant())
    assert dy.shape == (3,)


def test_parameter_budget_is_twelve_free_groups():
    """The identifiability budget, asserted rather than trusted.

    Thirteen fields, of which `n_C` is **pinned at 1 by default** on the
    double-counting argument in its docstring — so twelve are free. v3 Part 2
    asks for ~9-12. This is the gate that stops the count creeping back up one
    convenient term at a time, which is how the previous model reached 61.
    """
    from dataclasses import MISSING, fields

    f = fields(core.Params)
    assert len(f) == 13
    pinned = [x for x in f if x.default is not MISSING]
    assert [x.name for x in pinned] == ["n_C"]
    assert pinned[0].default == 1.0
    assert len(f) - len(pinned) == 12


def test_erk_suppresses_basal_ignition():
    """Cut #1: ERK shuts down the PTF1-independent ignition promoter.

    Checked with the loop broken (`R = 0`) so ignition is the only P source and
    nothing else can mask the sign.
    """
    p = core.default_params()
    y = np.array([0.0, 0.0, 0.0])
    dP = [core.rhs(0.0, y, p, core.Inputs.constant(erk=e))[0]
          for e in (0.0, 0.5, 2.0, 10.0)]
    assert np.all(np.diff(dP) < 0), "rising ERK must suppress PTF1A ignition"


def test_ignition_is_present_at_all():
    """Without basal ignition, `P = 0` is absorbing and trametinib alone can
    never revert anything — contradicting Collins 2014 and destroying Phase 7's
    internal positive control. `b_P = 0` is the bug this guards."""
    p = core.default_params()
    dP = core.rhs(0.0, np.array([0.0, 0.0, 0.0]), p,
                  core.Inputs.constant(erk=0.0))[0]
    assert dP > 0.0


def test_memory_represses_the_acinar_program():
    p = core.default_params()
    dP = [core.rhs(0.0, np.array([1.0, 1.0, c]), p, core.Inputs.constant(erk=0.1))[0]
          for c in (0.0, 0.5, 2.0, 10.0)]
    assert np.all(np.diff(dP) < 0)


def test_erk_writes_chromatin_memory():
    p = core.default_params()
    dC = [core.rhs(0.0, np.array([0.0, 0.0, 0.0]), p, core.Inputs.constant(erk=e))[2]
          for e in (0.0, 0.5, 2.0, 10.0)]
    assert dC[0] == 0.0, "no ERK, no writing"
    assert np.all(np.diff(dC) > 0)


def test_chromatin_is_the_slow_variable():
    """`gamma = delta_C/delta_P` must be able to be << 1, and `C` must actually
    move on that timescale. `C` sets time-to-relapse, which *is* R3 — if it
    relaxes as fast as the proteins, there is no memory and no durability
    question to ask."""
    p = core.default_params()
    assert p.gamma < 0.1
    inp = core.Inputs.constant(erk=2.0)
    s = solve_ivp(core.rhs, (0, 1.0 / p.gamma), [0.0, 0.0, 0.0], args=(p, inp),
                  method="LSODA", rtol=1e-9, atol=1e-11)
    assert 0.1 < s.y[2, -1] < 5.0, "C should be mid-transit after ~1/gamma"


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def test_erk_withdrawal_is_not_a_step():
    """Decision 002's amendment, as a test.

    A step would be the thing that decision explicitly rejected: the whole
    endpoint is what happens *after* withdrawal, so the recovery shape is
    mechanism. This requires a genuine ramp, not a jump.
    """
    erk = core.erk_withdrawal_profile(erk_high=3.0, erk_drug=0.05, t_withdraw=100.0,
                                      tau_rise=5.0, overshoot=0.0)
    assert erk(99.0) == pytest.approx(0.05)
    assert erk(100.0) == pytest.approx(0.05)
    mid = erk(105.0)
    assert 0.05 < mid < 3.0, "recovery must be gradual, not a step"
    assert erk(1000.0) == pytest.approx(3.0, rel=1e-6)


def test_erk_overshoot_peaks_where_it_says_it_does():
    """The overshoot is parameterized so `tau_overshoot` is the peak time and
    `overshoot` is its readable amplitude — because bench item 7 is meant to
    replace this curve with measured points, and a parameterization nobody can
    read off a western blot would not be replaceable."""
    erk = core.erk_withdrawal_profile(erk_high=2.0, erk_drug=0.0, t_withdraw=0.0,
                                      tau_rise=1e-6, overshoot=0.5,
                                      tau_overshoot=7.0)
    t = np.linspace(0.01, 60, 6000)
    v = np.array([erk(x) for x in t])
    assert t[np.argmax(v)] == pytest.approx(7.0, rel=0.05)
    assert v.max() == pytest.approx(2.0 * 1.5, rel=0.02)


def test_pulse_delivers_exactly_the_stated_dose():
    """Each pulse integrates to `dose`.

    Load-bearing for Phase 4: the mixture-amount design compares compositions at
    matched *total* mRNA mass, and without this normalisation changing a
    half-life would silently change how much was delivered — confounding every
    comparison by the quantity it is trying to control for.
    """
    u = core.mrna_pulse(dose=2.0, interval=1e9, k_translate=1.0, delta_m=0.05,
                        n_pulses=1)
    t = np.linspace(0, 400, 200_000)
    assert np.trapezoid([u(x) for x in t], t) == pytest.approx(2.0, rel=1e-3)


def test_pulse_rejects_inverted_timescales():
    """Decay faster than translation inverts the pulse, which would mean the
    payload *removes* protein."""
    with pytest.raises(ValueError, match="decay slower than translation"):
        core.mrna_pulse(1.0, 10.0, k_translate=0.01, delta_m=1.0, n_pulses=1)


def test_chop_flag_is_a_flag_not_a_term():
    """Decision 008 retired: viability is measured at the bench. The CHOP arm
    survives as a warning on the output only — it must not appear in the RHS."""
    assert core.chop_flag(np.array([1.0, 0.5, 0.01])) is True
    assert core.chop_flag(np.array([1.0, 0.5])) is False
    p = core.default_params()
    dy_lo = core.rhs(0.0, np.array([1e-6, 1.0, 0.5]), p, core.Inputs.constant())
    assert np.all(np.isfinite(dy_lo)), "low P must not be a stopping condition"


# ---------------------------------------------------------------------------
# Qualitative structure — what Gate B will test properly
# ---------------------------------------------------------------------------
# These assert the core CAN exhibit the required structure at a placeholder
# parameter set. They are not measurements of the persistence window, and no
# number from them may be reported. Gate B does this properly, with continuation
# and pre-registered ranges.


def _settle(p, erk, y0, t_end=8000.0):
    s = solve_ivp(core.rhs, (0, t_end), y0, args=(p, core.Inputs.constant(erk=erk)),
                  method="LSODA", rtol=1e-9, atol=1e-11)
    assert s.success, s.message
    return s.y[:, -1]


def test_core_admits_two_attractors_at_intermediate_erk():
    """Gate B's precondition: two stable states must exist somewhere.

    Same qualitative signature as the pre-rewrite golden fixture — and the
    metaplastic attractor again sits at `R ~ 0`, across a complete change of
    state space. That is the bootstrap claim surviving the rewrite.

    Note the ERK level: **intermediate**, not low. At low ERK this core is
    monostable *acinar* — trametinib alone reverts, which is what Collins 2014
    requires and what `default_params`' `b_P` was chosen to respect. Bistability
    lives in a window between that and the high-ERK metaplastic-only regime, and
    that window is the drug-free persistence window R3 is about.
    """
    p = core.default_params()
    acinar = _settle(p, 0.4, [3.0, 3.0, 0.02])
    adm = _settle(p, 0.4, [0.0, 0.0, 5.0])
    assert np.max(np.abs(acinar - adm)) > 1e-3, "expected two distinct endpoints"
    assert acinar[1] > 0.5, "acinar branch must sustain RBPJL"
    assert adm[1] < 0.1, "metaplastic branch must have collapsed RBPJL"


def test_the_bistable_window_is_bounded_on_both_sides():
    """A wedge, not a half-plane — which is what makes it a *window*.

    Below it, trametinib alone reverts (Collins 2014). Above it, KRAS holds
    metaplasia and nothing sticks. R3 asks where the upper edge sits, and a model
    whose bistable region ran to zero ERK or to infinity would not pose that
    question at all.
    """
    p = core.default_params()

    def bistable(erk):
        a = _settle(p, erk, [3.0, 3.0, 0.02])
        d = _settle(p, erk, [0.0, 0.0, 5.0])
        return np.max(np.abs(a - d)) > 1e-3

    assert not bistable(0.05), "low ERK should be monostable acinar"
    assert bistable(0.4), "an intermediate window must exist"
    assert not bistable(3.0), "high ERK should be monostable metaplastic"


def test_high_erk_collapses_the_acinar_state():
    """KRAS-driven ADM: from a healthy start, high ERK must drive the cell out.
    Without this the model cannot induce the disease it is about."""
    p = core.default_params()
    assert _settle(p, 3.0, [3.0, 3.0, 0.02])[1] < 1e-3


def test_withdrawal_relapses_without_a_payload():
    """Collins 2014's actual observation: revert under drug, relapse after it.

    This is the phenomenon the whole project exists to address, so the core must
    reproduce it before anything is built on top. Note what is NOT asserted here:
    that a payload prevents the relapse. That is a Phase 5 result, it depends on
    the drug-hold duration relative to `gamma`, and it must be pre-registered
    before it is swept — not tuned into a test.
    """
    p = core.default_params()
    erk = core.erk_withdrawal_profile(erk_high=3.0, erk_drug=0.05,
                                      t_withdraw=600.0, tau_rise=5.0,
                                      overshoot=0.3, tau_overshoot=3.0)
    s = solve_ivp(core.rhs, (0, 4000.0), [0.0, 0.0, 1.5],
                  args=(p, core.Inputs(erk=erk)), method="LSODA",
                  rtol=1e-9, atol=1e-11, t_eval=np.linspace(0, 4000, 2000))
    assert s.success
    R = s.y[1]
    at_withdrawal = R[np.searchsorted(s.t, 600.0)]
    assert at_withdrawal > 0.5, "drug should have reverted the cell first"
    assert R[-1] < 1e-3, "and it should relapse once the drug is withdrawn"


# ---------------------------------------------------------------------------
# Regressions on defects found by the Tier 2 panels (session 6)
# ---------------------------------------------------------------------------


def test_b_P_default_clears_the_saddle_node_with_margin():
    """`b_P` sits near a boundary the literature forbids crossing.

    Below `B_P_CRITICAL` the model says trametinib ALONE never reverts anything,
    contradicting Collins 2014 and destroying Phase 7's trametinib-only positive
    control. The default previously shipped at 0.5 against a measured boundary of
    0.4903 — **2% of margin on a fitted parameter**, so any sweep nudging it down
    a couple of percent crossed into the forbidden regime while every figure
    still rendered.

    The docs also said the boundary was "roughly 0.4", which was wrong: at 0.40
    trametinib alone leaves R ~ 0.18 against a high branch of ~4.
    """
    p = core.default_params()
    assert p.b_P > core.B_P_CRITICAL * 1.15, (
        f"b_P = {p.b_P} is within 15% of the critical value "
        f"{core.B_P_CRITICAL}; a fit or sweep could cross it silently")


def test_below_b_P_critical_trametinib_alone_cannot_revert():
    """Guard-on-the-guard: confirm the boundary is real and in the stated place."""
    from dataclasses import replace
    p = core.default_params()
    below = _settle(replace(p, b_P=core.B_P_CRITICAL * 0.9), 0.05, [0.0, 0.0, 1.5])
    above = _settle(replace(p, b_P=core.B_P_CRITICAL * 1.1), 0.05, [0.0, 0.0, 1.5])
    assert below[1] < 0.5, "below the critical b_P, trametinib alone must fail"
    assert above[1] > 0.5, "above it, trametinib alone must revert"


def test_eps_default_is_a_filter_not_a_memory_and_the_code_says_so():
    """`C` is only a bistable memory above `EPS_MEMORY_THRESHOLD` = 3*sqrt(3)/8.

    The default is below it, so `C` is a lagged filter. That may be acceptable —
    a filter still delays relapse — but the model must not be *described* as
    carrying a memory while running here. This test exists so the discrepancy is
    visible rather than buried in a comment.
    """
    assert core.EPS_MEMORY_THRESHOLD == pytest.approx(1.5396, abs=1e-3)
    assert core.default_params().eps < core.EPS_MEMORY_THRESHOLD


def test_chromatin_is_currently_a_strict_cascade():
    """DOCUMENTS A KNOWN LIMITATION — decision 014, not yet resolved.

    `dC/dtau` depends on the input and on `C`, and on nothing else: not `P`, not
    `R`, not the payload. So the payload has no channel to the durability
    endpoint, and Phase 5's dose x interval map would be flat by construction.

    This test asserts the limitation *as it currently stands*, so that whichever
    fix is chosen (active erasure, or restating Phase 5's endpoint as
    reachability) **this test must be updated deliberately** rather than the
    change slipping through unnoticed.
    """
    p = core.default_params()
    base = core.rhs(0.0, np.array([1.0, 1.0, 0.5]), p,
                    core.Inputs.constant(erk=1.0))[2]
    for P_ in (0.0, 5.0):
        for R_ in (0.0, 4.0):
            got = core.rhs(0.0, np.array([P_, R_, 0.5]), p,
                           core.Inputs.constant(erk=1.0))[2]
            assert got == pytest.approx(base, rel=1e-12)
    dosed = core.rhs(0.0, np.array([1.0, 1.0, 0.5]), p,
                     core.Inputs.constant(erk=1.0, u_P=10.0, u_R=10.0))[2]
    assert dosed == pytest.approx(base, rel=1e-12), (
        "if this now differs, active erasure was added — update decision 014")
