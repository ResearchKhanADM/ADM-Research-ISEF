"""Tests for the supplementary three-species binding solver.

`src/supplementary/binding.py` is **not part of the Phase 2 core.** The core
titrates one target and uses a closed form (`core.free_e_protein`); this module
solves the full competitive equilibrium — ID3 against both E-protein and PTF1A,
with two ternary complexes — numerically.

It is kept, and tested, for one reason: it is the only implementation of the
*two-target* mechanism, and Phase 3 may need it. Phase 3 convolves the per-cell
LNP dose distribution against the bootstrap threshold, and if that threshold's
sharpness turns out to depend on ID3 taxing PTF1A as well as E, the core's
one-target simplification is not good enough and this comes back as the binding
step. Deleting it would mean re-deriving 283 lines of algebra to find out.

What was **dropped** from the previous version of this file, and why:

  * everything referencing `Topology`, `payload_subsets()` or the T1/T2 Q-value
    comparison — the five-way topology competition is cut (decision 012);
  * `test_pd325901_recovers_when_f_act_disabled` and the `raf_drive` sign tests —
    they guarded the `W` state, retired with decision 002's amendment. Their
    *discipline* (a guard-on-the-guard that constructs the violating version and
    requires it to fail) moved to `tests/test_bootstrap_guard.py`, which is where
    the project's central claim now lives.

Run:  venv\\Scripts\\python.exe -m pytest tests/test_supplementary_binding.py -v
"""

from __future__ import annotations

import numpy as np
import pytest

from src.supplementary import binding


def _effective_hill(fn, kd, i_grid):
    """n_eff = max |d ln C_L / d ln I| — ultrasensitivity of the TERNARY COMPLEX.

    NOT a global fold-change across a fixed ID3 range. Titration is
    threshold-linear: very steep near the threshold, flat away from it, so a
    global fold-change averages the steep part away and can rank a first-order
    sink as 'sharper'. That mistake was made here once; this docstring exists so
    it is not made again.

    Note the observable. This measures the slope of `C_L`, the ternary complex,
    under two-target titration — which is why its tight-limit prefactor is 1.34
    and its loose limit is 2. The Phase 2 core measures the slope of `E_free`
    under one-target titration and gets 0.5 and 1. **The constants are not
    interchangeable**; carrying 1.34 into the core would overstate threshold
    sharpness ~2.7x and inflate the Phase 3 co-formulation gap with it.
    """
    c_l = np.array([fn(1.0, 1.0, i, 1.0, kd, 0.1)["C_L"] for i in i_grid])
    return np.max(np.abs(np.gradient(np.log(np.maximum(c_l, 1e-300)),
                                     np.log(i_grid))))


# ---------------------------------------------------------------------------
# The mechanism
# ---------------------------------------------------------------------------


def test_two_target_titration_traps_both_species():
    """ID3 must reduce free PTF1A *and* free E — Dufresne 2010, PMID 20830706.

    This is the property the core does NOT have, and the reason this module is
    kept: a titrator acting on only one partner is a weaker mechanism and gives a
    softer threshold.
    """
    kd = {"IE": 1.0, "IP": 1.0, "L": 1.0, "J": 10.0}
    lo = binding.solve_binding(1.0, 1.0, 0.01, 1.0, kd, rbpj=0.1)
    hi = binding.solve_binding(1.0, 1.0, 5.00, 1.0, kd, rbpj=0.1)
    assert hi["P_free"] < lo["P_free"], "ID3 must sequester PTF1A"
    assert hi["E_free"] < lo["E_free"], "ID3 must sequester E-protein"
    assert hi["C_L"] < lo["C_L"], "and therefore suppress the PTF1-L complex"


def test_binding_conserves_mass():
    """A conservation leak would look like spontaneous protein synthesis."""
    kd = {"IE": 2.0, "IP": 3.0, "L": 1.5, "J": 8.0}
    p_tot, e_tot, i_tot, r_tot, rbpj = 1.3, 0.9, 0.7, 1.1, 0.2
    b = binding.solve_binding(p_tot, e_tot, i_tot, r_tot, kd, rbpj)
    assert b["P_free"] + b["ID3_P"] + b["C_L"] + b["C_J"] == pytest.approx(p_tot, rel=1e-11)
    assert b["E_free"] + b["ID3_E"] + b["C_L"] + b["C_J"] == pytest.approx(e_tot, rel=1e-11)
    assert b["I_free"] + b["ID3_E"] + b["ID3_P"] == pytest.approx(i_tot, rel=1e-11)
    assert b["R_free"] + b["C_L"] == pytest.approx(r_tot, rel=1e-11)


def test_binding_matches_the_derivation():
    """Code and `docs/derivations/binding_polynomial.md` must agree.

    Recomputes the conservation laws from the derivation's factorised form
    (§2-§4) and checks the solver satisfies them. If the derivation is edited
    without the code, or vice versa, this fails — which matters more now that the
    derivation is the stated justification for the core's algebraic `E_free`.
    """
    kd = {"IE": 0.02, "IP": 0.5, "L": 1.5, "J": 8.0}
    p_tot, e_tot, i_tot, r_tot, rbpj = 1.3, 0.9, 0.7, 1.1, 0.2
    b = binding.solve_binding(p_tot, e_tot, i_tot, r_tot, kd, rbpj)
    P, E, I, R = b["P_free"], b["E_free"], b["I_free"], b["R_free"]

    phi = r_tot / (kd["L"] + P * E) + rbpj / kd["J"]          # derivation §3
    assert P * (1 + I / kd["IP"] + E * phi) == pytest.approx(p_tot, rel=1e-11)
    assert E * (1 + I / kd["IE"] + P * phi) == pytest.approx(e_tot, rel=1e-11)
    assert I * (1 + E / kd["IE"] + P / kd["IP"]) == pytest.approx(i_tot, rel=1e-11)
    assert R * (1 + P * E / kd["L"]) == pytest.approx(r_tot, rel=1e-11)


# ---------------------------------------------------------------------------
# The two limits — these are what tell Phase 3 which regime it is in
# ---------------------------------------------------------------------------


def test_loose_limit_has_log_log_slope_minus_two():
    """Derivation §5. Two multiplied factors (ID3 taxes P and E), hence -2.

    Contrast with the core's one-target form, whose loose limit is -1. That
    factor of two between the mechanisms is exactly why their `n_eff` constants
    cannot be swapped.
    """
    kd = {"IE": 50.0, "IP": 50.0, "L": 1.0, "J": 10.0}
    i_grid = np.logspace(2.5, 4.0, 60)
    c_l = np.array([binding.solve_binding(1.0, 1.0, i, 1.0, kd, 0.1)["C_L"]
                    for i in i_grid])
    slope = np.gradient(np.log(c_l), np.log(i_grid))[-1]
    assert slope == pytest.approx(-2.0, abs=0.05)


def test_tight_limit_scales_as_sqrt_of_total_over_kd():
    """Derivation §6: n_eff ~= 1.34*sqrt(E_tot/Kd) for the ternary complex.

    The scaling — `sqrt(totals/Kd)`, unbounded as Kd falls — is what transfers to
    the core. **The prefactor does not.** See the module docstring.
    """
    i_grid = np.logspace(-1.2, 1.2, 400)
    for kd_val in (1e-2, 1e-3):
        kd = {"IE": kd_val, "IP": kd_val, "L": 1.0, "J": 10.0}
        n = _effective_hill(binding.solve_binding, kd, i_grid)
        predicted = 1.34 * np.sqrt(1.0 / kd_val)
        assert n == pytest.approx(predicted, rel=0.15), (
            f"tight-binding scaling broke at Kd={kd_val}: n_eff={n:.2f}, "
            f"expected ~{predicted:.2f}"
        )


def test_two_target_is_sharper_than_the_core_one_target_form():
    """Quantifies what the core's simplification costs, at matched Kd.

    This is the number that decides whether the core's one-target `E_free` is
    good enough for Phase 3. If two-target titration is much sharper, the
    converted-fraction curve and therefore the co-formulation gap change, and the
    core must borrow this module for its binding step. Recorded as a test so the
    comparison is a measurement rather than an assumption.
    """
    from src import core

    kd_val = 1e-2
    kd = {"IE": kd_val, "IP": kd_val, "L": 1.0, "J": 10.0}
    two_target = _effective_hill(binding.solve_binding, kd,
                                 np.logspace(-1.2, 1.2, 400))
    one_target = core.n_eff(kd_val)
    assert two_target > one_target, (
        "two-target titration must be the sharper mechanism; if it is not, the "
        "core's simplification is not a simplification and this module's reason "
        "for existing is gone"
    )
    # ~13.4 vs ~5.0 at Kd/E_tot = 1e-2. Asserted loosely, as a regression guard
    # on the ratio rather than a claim about either number.
    assert 2.0 < two_target / one_target < 4.0


# ---------------------------------------------------------------------------
# Convergence accounting — the standing rule
# ---------------------------------------------------------------------------


def test_convergence_failure_rate_is_bounded_in_the_tight_regime():
    """MANDATORY while this solver is in use (CLAUDE.md standing rule).

    Convergence failures correlate with Kd, so silently dropping them would
    deplete the sample set precisely in the tight regime — the regime that
    matters — and produce a confident wrong answer that looks like a clean null.

    Note the core does **not** need this: its binding step is closed-form and
    cannot fail to converge. The standing rule still applies to continuation and
    to Phase 3's Monte Carlo.
    """
    rng = np.random.default_rng(0)
    ledger = binding.ConvergenceLedger()
    n = 400
    for _ in range(n):
        kd_val = 10 ** rng.uniform(-6, -1)
        kd = {"IE": kd_val, "IP": 10 ** rng.uniform(-6, 0),
              "L": 10 ** rng.uniform(-1, 1), "J": 10 ** rng.uniform(0, 2)}
        totals = 10 ** rng.uniform(-1, 1, size=4)
        try:
            binding.solve_binding(*totals, kd, rbpj=0.1, ledger=ledger)
        except binding.BindingSolveError:
            pass          # counted by the ledger, never silently dropped
    assert ledger.attempts == n, "every attempt must be recorded"
    assert ledger.failure_rate < 0.01, (
        f"tight-regime failure rate {ledger.failure_rate:.3%} exceeds the 1% "
        f"bound. {ledger.summary()}"
    )


def test_ledger_records_failures_rather_than_hiding_them():
    """Guard on the guard: a ledger that never records a failure proves nothing."""
    ledger = binding.ConvergenceLedger()
    ledger.record(True, kd_over_etot=1.0)
    ledger.record(False, kd_over_etot=1e-6)
    assert ledger.attempts == 2 and ledger.failures == 1
    assert ledger.failure_rate == 0.5
    assert ledger.failed_conditions[0]["kd_over_etot"] == 1e-6
