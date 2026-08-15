"""Explicit functional forms for the ADM model.

Master plan §3.4 specifies *which terms exist and why*; choosing the forms is
Stage 0's job. Every choice here has a `docs/decisions/` file. They live in one
module rather than inline in the right-hand side for a reason: Stage 2 compares
five topologies under an identical sampling box, and that comparison is only
valid if the shared machinery is literally shared. A Hill function copied into
two topologies is a Hill function that can drift.

Sign conventions that carry mechanism are asserted here, not just commented.
An inverted sign in this file produces a model that runs cleanly and predicts
the opposite result, which is the worst class of bug available to us.
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Hill functions
# ---------------------------------------------------------------------------


def hill_activate(x, k, n):
    """Activating Hill:  x^n / (k^n + x^n).  Rises 0 → 1.

    Used for the 2.3-kb autoregulatory enhancer (activated by the PTF1-L
    complex) and for Rbpjl transcription. `n` is *sampled over 1–4*, never
    fixed: the master plan's answer to "where did the Hill coefficient come
    from" is that PTF1-L is an obligate trimer with two conserved sites, which
    justifies cooperativity but not a specific value. Fixing n would be
    asserting the thing we are supposed to be uncertain about.
    """
    # np.maximum guards x=0 with fractional n, which would otherwise be nan.
    x = np.maximum(x, 0.0)
    xn = x**n
    return xn / (k**n + xn)


def hill_repress(x, k, n):
    """Repressing Hill:  1 / (1 + (x/k)^n).  Falls 1 → 0.

    Separable from activation on purpose. There is no evidence for a specific
    interaction form between the metaplasia program and the acinar enhancer, so
    the minimal assumption is that repression multiplies activation rather than
    competing within one binding expression.
    """
    x = np.maximum(x, 0.0)
    return 1.0 / (1.0 + (x / k) ** n)


# ---------------------------------------------------------------------------
# ERK / drug layer  — the four cuts enter here
# ---------------------------------------------------------------------------


def g_ignition(k_eff, k_ign, q):
    """The 13.4-kb PTF1-independent 'ignition' promoter, suppressed by ERK.

    DECREASING in k_eff. This is cut #1 of the four in §2.3: ERK shuts down the
    only route by which PTF1A can be made while the loop is open.

    Takes k_eff, NOT (K, v). k_eff already carries the drug through W and f_cat,
    so passing the drug in again here would double-count it — an easy and
    invisible error, since the model would still run and merely respond to
    trametinib twice as strongly as it should.
    """
    return hill_repress(k_eff, k_ign, q)


def f_cat(v, ic50_cat):
    """MEK *catalytic* inhibition — blocks MEK's output. Shared by both drugs.

    Fractional residual activity, so f_cat(0) = 1 (no drug, full activity).
    """
    return 1.0 / (1.0 + v / ic50_cat)


def f_act(v, ic50_act, is_trametinib):
    """Inhibition of RAF-mediated MEK *phosphorylation* — blocks MEK activation.

    Trametinib-specific. Setting is_trametinib=False returns 1.0 identically,
    which recovers PD325901 — and that is what preserves the direct comparison
    to Collins 2014, who used PD325901. The two drugs are not a potency
    difference in this model; they act at different points in the cascade.
    """
    if not is_trametinib:
        return 1.0
    return 1.0 / (1.0 + v / ic50_act)


def raf_drive(k_eff, k_fb, r):
    """RAF→MEK activating drive. MUST BE STRICTLY DECREASING IN k_eff.

    This is the mechanism, and its sign is the whole point. Phospho-MEK
    accumulates under a catalytic inhibitor *because* falling ERK output
    relieves ERK-mediated negative feedback on RAF, so RAF drive rises as ERK
    falls. Implement it as increasing and the model predicts the opposite
    result while running perfectly cleanly.

    Covered by tests/test_model_invariants.py::test_raf_drive_is_decreasing,
    which is written to fail if this is ever inverted.
    """
    return hill_repress(k_eff, k_fb, r)


# ---------------------------------------------------------------------------
# Viability — a U-shaped hazard, never a threshold
# ---------------------------------------------------------------------------


def hazard(s, p_free, params):
    """Instantaneous death hazard h(S, P_n). Rises at BOTH ends.

    survival(t) = exp(-∫h dt) is accumulated as a state, so this returns a rate.

    Two arms, deliberately additive rather than multiplicative — the cell can
    die of either cause independently, and a product would make each arm
    require the other to be non-zero:

      h_high(S)      cargo outrunning secretory capacity during the transient
      h_low(P_free)  CHOP-dependent apoptosis under PTF1A loss
                     (Sakikubo 2018 PMID 30361559; Backx 2021 PMID 33762742)

    Why a hazard and not "death when S > S_crit": a one-sided threshold on S is
    the U_crit construct the expert panel killed. The hazard integrates along
    the trajectory, so a brief excursion costs survival without automatically
    killing the cell — which matches the panel's actual finding, that the risk
    is a RATE mismatch during the transient rather than an amplitude ceiling.
    It also yields a continuous 0–1 number, which is what the y-axis of the
    reversal–viability figure needs; a binary alive/dead flag cannot be plotted
    against a continuous reversal axis.
    """
    h_high = params["h_max_cargo"] * hill_activate(s, params["s_crit"], params["nu_s"])
    h_low = params["h_max_chop"] * hill_repress(p_free, params["p_crit"], params["mu_p"])
    return h_high + h_low


# ---------------------------------------------------------------------------
# mRNA pulse forcing — analytic, not extra states
# ---------------------------------------------------------------------------


def mrna_pulse(t, dose, interval, k_translate, delta_m, n_pulses):
    """Delivered mRNA input u_i(t): repeated pulses, difference of exponentials.

    Shape: uptake/translation ramp (rate k_translate) then first-order decay
    (rate delta_m). §1.3 insists the input is a pulse of known shape rather than
    a step of indefinite duration — that is what makes "duration" a real
    quantity set by half-life plus redosing interval, and it is what makes
    Stage 5's axes (dose per pulse) × (redosing interval) actuatable.

    Normalised so each pulse delivers exactly `dose` (∫ over one pulse = dose),
    which is what lets Stage 5's necessity analysis compare subsets at matched
    TOTAL delivered dose. Without that normalisation, changing the half-life
    would silently change the amount delivered and every comparison would be
    confounded.

    STAGE 6 PRECONDITION — verified symbolically, see decision 011. The input
    enters as (scalar coefficient in t) × (constant direction in state space),
    so ∂u/∂x = 0 and the supply vector fields still commute: [g1,g2] = 0. The
    entire ordering effect therefore still comes from the nonlinearity of f,
    and Stage 6's closed form survives. Caveat recorded there: the s² expansion
    assumes autonomous fields, so the closed form applies to constant-amplitude
    holds; the pulse shape belongs to Stage 6's simulation arm.
    """
    if delta_m >= k_translate:
        raise ValueError(
            "delta_m must be < k_translate: decay slower than translation. "
            "Otherwise the 'ramp' outruns the decay and the pulse inverts."
        )
    # ∫₀^∞ (e^{-at} - e^{-bt}) dt = 1/a - 1/b, with a=delta_m < b=k_translate.
    area = 1.0 / delta_m - 1.0 / k_translate
    total = 0.0
    for k in range(n_pulses):
        dt = t - k * interval
        if dt >= 0.0:
            total += np.exp(-delta_m * dt) - np.exp(-k_translate * dt)
    return dose * total / area
