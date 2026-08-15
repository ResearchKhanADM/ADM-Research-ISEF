"""The Phase 2 minimal core — 3 states, 12 dimensionless groups.

Replaces the 11-state / 61-parameter system (deleted in the same commit; see
`docs/decisions/012-durability-framing-architecture-change.md`). Everything here
is nondimensional: time in units of 1/delta_P, P in units of K_P, R in units of
K_R, C in units of theta_C, and E-protein / ID3 in units of E_tot.

STATES
    P   PTF1A activity — autoregulatory, requires an E-protein partner
    R   RBPJL           — produced ONLY as a function of P
    C   chromatin/memory at metaplasia loci — sets time-to-relapse

INPUTS, prescribed functions of time, never states
    erk       **the ERK-driven ID3 level, in units of E_tot** — see the naming
              note below. Trametinib sets it; on withdrawal it follows a REBOUND
              PROFILE, not a step (decision 002 amendment)
    u_P, u_R  delivered mRNA, as analytic pulses (decision 011)

ALGEBRAIC, never integrated
    e_free    the EXACT binding solution (decision 006 amendment) — not
              `E_tot - k*ID3`, which is only its tight-binding limit

NAMING NOTE — `erk` IS `ID3`, AND THAT IS DELIBERATE
    An earlier version of this docstring promised "`id3`, a saturating function
    of `erk`". **No such function exists and none is intended.** The input is
    taken as ID3 directly, proportional to pERK, with the constant of
    proportionality absorbed into the E_tot scale — which is exactly what makes
    `kappa = K_d/E_tot` a single group rather than two.

    So the field is called `erk` but carries an ID3 level, and every figure axis
    must say so. The linearity is not a hidden assumption: it IS the ERK->ID3
    edge, the one flagged in v3 Part 1.4 as the load-bearing gap, and
    **Bench Handshake item 8 (ID3 western +/- trametinib) measures precisely this
    input map.** If that western comes back non-linear or flat, this line is
    where the model changes.

READ THIS BEFORE EDITING
  * `dR/dt` HAS NO P-INDEPENDENT TERM. That zero is the project's central claim:
    *Rbpjl* has no PTF1A-independent promoter, so the loop cannot re-close on its
    own. Adding a basal term looks harmless, improves numerical behaviour, and
    inverts the conclusion without breaking anything. Guarded by
    `tests/test_bootstrap_guard.py`, which was written before this file.
  * `e_free` uses the exact quadratic. Do NOT substitute `1 - id3` with a
    `max(0, .)` floor — the negativity that floor would hide is the
    approximation announcing it has left its domain, not a numerical nuisance.
  * Viability is NOT modelled (decision 008 retired). `chop_flag()` is a warning
    on the output, not a term in the right-hand side.

The three profiled parameters map one-to-one onto the three headline results:
    a_P    does the loop close    -> R2 composition
    gamma  how long it holds      -> R3 durability
    kappa  threshold sharpness    -> R1 formulation
So profile likelihood is not an identifiability side-quest — it is the
uncertainty bar on each deliverable. See `docs/PHASE2_PARAMETER_BUDGET.md`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

#: Fixed state order. Load-bearing — every index downstream depends on it.
STATES: tuple[str, ...] = ("P", "R", "C")

#: Exponent of the chromatin self-reinforcement term, FIXED at 2 rather than
#: sampled. Two reasons, and the second is the important one: 2 is the minimal
#: value that lets the C-subsystem hold a memory at all, and sampling it as well
#: as `n_C` would parameterize the same switching physics twice — the
#: double-counting trap that argues for `n_C = 1` (see Params.n_C).
SELF_REINFORCEMENT_EXPONENT = 2.0

#: `eps` above which the C-subsystem is genuinely BISTABLE — i.e. a memory that
#: holds a written state on its own — rather than a low-pass filter with a lag.
#: Exact: `max_C d/dC [C^2/(1+C^2)] = 3*sqrt(3)/8 = 0.6495`, so the threshold is
#: its reciprocal, 1.5396.
#:
#: **`default_params()` ships `eps = 0.5`, which is BELOW this.** At the default,
#: `C` is a filter, not a memory. That is not necessarily wrong — a lagged filter
#: still delays relapse — but the word "memory" must not be used for it on a
#: poster, and any prior range for `eps` should straddle this threshold so the
#: question is tested rather than assumed.
EPS_MEMORY_THRESHOLD = 8.0 / (3.0 * np.sqrt(3.0))     # = 1.5396

#: Measured boundary of `b_P` below which trametinib alone cannot revert, at the
#: default parameter set (bisection). See `default_params`. Any sweep or fit that
#: allows `b_P` below this is exploring a regime the literature has ruled out,
#: and should say so rather than average over it.
B_P_CRITICAL = 0.4903


# ---------------------------------------------------------------------------
# Parameters — 12 dimensionless groups
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Params:
    """The whole model. 12 groups: 10 fitted, 2 scanned.

    Frozen because a parameter set that mutates mid-run stops describing the
    result it produced.
    """

    # --- P: autoregulation, basal ignition, repression by memory
    a_P: float          # autoregulatory gain — BISTABILITY LIVES HERE. Profiled.
    b_P: float          # basal (PTF1-independent ignition) gain, ERK-suppressed
    c_rep: float        # how strongly memory represses the acinar program
    n_P: float          # Hill exponent, SCANNED 1-4, never fitted (decision 004)

    # --- R: driven only by P. There is no fourth parameter here, and its
    #     absence is the bootstrap claim.
    a_R: float          # Rbpjl gain
    n_R: float          # Hill exponent, SCANNED 1-4
    rho: float          # delta_R / delta_P — RBPJL vs PTF1A turnover

    # --- C: slow memory
    gamma: float        # delta_C / delta_P — THE DURABILITY KNOB. Profiled.
    alpha_C: float      # write gain
    k_w: float          # ERK half-max for writing, in ID3 units
    eps: float          # self-reinforcement strength. Below EPS_MEMORY_THRESHOLD
                        # (1.5396) the C-subsystem is a lagged FILTER, not a
                        # bistable memory — and the default is below it. Do not
                        # call C a "memory" without checking this value.

    # --- binding
    kappa: float        # K_d / E_tot — THE BINDING REGIME. Profiled.

    #: Repression exponent. Defaults to 1 on a structural argument, not an
    #: empirical one: C's switching already comes from its own self-reinforcement
    #: term, so making the repression of P by C cooperative *as well* encodes the
    #: same physics twice and the two mechanisms cannot then be distinguished.
    #: That is a small identifiability trap of exactly the kind this reduction
    #: exists to avoid. Set it to a sampled value only to run that check.
    n_C: float = 1.0

    def __post_init__(self) -> None:
        if self.kappa <= 0:
            raise ValueError("kappa = K_d/E_tot must be > 0")
        for name in ("n_P", "n_R", "n_C"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be > 0")


def default_params() -> Params:
    """A mid-range set for smoke tests and guard tests.

    NOT a calibration and NOT a prediction. Every value is an order-of-magnitude
    placeholder until the Phase 0 bench items land; report no number obtained
    from these.

    **`b_P` sits close to a saddle-node and that is a hazard, not a detail.**
    Below `B_P_CRITICAL` the model says MEK inhibition *alone* can never restore
    RBPJL — contradicting Collins 2014 and destroying Phase 7's trametinib-only
    positive control. The measured boundary is **0.4903** (bisection, at this
    parameter set); an earlier version of this docstring said "roughly 0.4",
    which is **wrong** — at `b_P = 0.40` trametinib alone leaves `R = 0.18`
    against a high branch of 3.98, i.e. not reverted.

    The default is set to **0.6**, ~22% above the boundary, rather than the 0.5
    it previously shipped at — 0.5 is only **2% clear**, and `b_P` is one of the
    ten *fitted* groups, so any fit or sweep that nudges it down by a couple of
    percent silently crosses into the forbidden regime while every figure still
    renders. See `B_P_CRITICAL`.
    """
    return Params(
        a_P=6.0, b_P=0.6, c_rep=1.0, n_P=3.0,
        a_R=4.0, n_R=3.0, rho=1.0,
        gamma=0.02, alpha_C=1.5, k_w=1.0, eps=0.5,
        kappa=0.01, n_C=1.0,
    )


# ---------------------------------------------------------------------------
# Inputs — prescribed in time, never states
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Inputs:
    """pERK and the two mRNA pulses, as functions of scaled time.

    Callables rather than parameters, because the point of decision 002's
    amendment is that pERK on withdrawal is a *shape*, not a level. Shape
    parameters are swept or measured (bench item 7) — they are not fitted, and
    they are deliberately not in `Params`.
    """

    erk: object          # callable tau -> pERK-driven ID3 level, in E_tot units
    u_P: object = None   # callable tau -> delivered PTF1A mRNA
    u_R: object = None   # callable tau -> delivered RBPJL mRNA

    @staticmethod
    def constant(erk=1.0, u_P=0.0, u_R=0.0) -> "Inputs":
        return Inputs(erk=lambda tau: erk, u_P=lambda tau: u_P,
                      u_R=lambda tau: u_R)

    def eval(self, tau) -> tuple[float, float, float]:
        return (float(self.erk(tau)),
                0.0 if self.u_P is None else float(self.u_P(tau)),
                0.0 if self.u_R is None else float(self.u_R(tau)))


def erk_withdrawal_profile(erk_high, erk_drug, t_withdraw,
                           tau_rise, overshoot=0.0, tau_overshoot=1.0):
    """pERK across trametinib treatment and withdrawal. **Not a step.**

    Why this exists at all: the primary endpoint of the whole project is what
    happens *after* withdrawal, so the recovery shape is mechanism rather than
    detail. It sets how hard the system is pushed back toward the metaplastic
    basin at exactly the moment durability is being tested.

    Shape, three parameters, all swept or measured:
      `tau_rise`       recovery timescale toward `erk_high`
      `overshoot`      relative amplitude of the transient above `erk_high`,
                       from relief of ERK-mediated negative feedback on RAF
      `tau_overshoot`  when that transient peaks

    The overshoot is an alpha function `(s/tau)*exp(1 - s/tau)`, which peaks at
    exactly `s = tau_overshoot` with amplitude 1 — so `overshoot` is readable
    directly off a western blot rather than being an abstract coefficient. That
    matters because bench item 7 is meant to *replace* this curve with measured
    points, and a parameterization nobody can read off a gel would not be
    replaceable.

    In a bistable system a transient can have a permanent consequence: the
    overshoot need not persist, it need only carry the slow state back across
    the separatrix. That was the surviving insight of decision 002, and it now
    lives here instead of in a phospho-MEK state.
    """
    if tau_rise <= 0 or tau_overshoot <= 0:
        raise ValueError("timescales must be positive")

    def erk(tau):
        if tau < t_withdraw:
            return erk_drug
        s = tau - t_withdraw
        recovery = erk_drug + (erk_high - erk_drug) * (1.0 - np.exp(-s / tau_rise))
        spike = (overshoot * (erk_high - erk_drug)
                 * (s / tau_overshoot) * np.exp(1.0 - s / tau_overshoot))
        return recovery + spike

    return erk


def mrna_pulse(dose, interval, k_translate, delta_m, n_pulses):
    """Delivered mRNA as repeated pulses: difference of exponentials.

    Carried over from decision 011, which survives the rewrite unchanged and in
    fact matters *more* now: Phase 4 is a mixture-amount design under a fixed
    total mRNA mass, and comparing compositions at matched total mass is only
    meaningful if each pulse delivers exactly `dose`. Without the normalisation,
    changing a half-life would silently change how much was delivered and every
    composition comparison would be confounded by the thing it controls for.
    """
    if delta_m >= k_translate:
        raise ValueError(
            "delta_m must be < k_translate: decay slower than translation, or "
            "the pulse inverts and the payload would *remove* protein"
        )
    area = 1.0 / delta_m - 1.0 / k_translate

    def u(tau):
        total = 0.0
        for k in range(n_pulses):
            s = tau - k * interval
            if s >= 0.0:
                total += np.exp(-delta_m * s) - np.exp(-k_translate * s)
        return dose * total / area

    return u


# ---------------------------------------------------------------------------
# Algebra — binding and Hill terms
# ---------------------------------------------------------------------------


def hill(x, n):
    """x^n / (1 + x^n). Rises 0 -> 1, half-max at x = 1 by construction.

    Half-max is 1 because every concentration here is already scaled by its own
    half-max — that is what nondimensionalization bought, and it is why 19
    dimensional parameters became 12 groups.
    """
    x = np.maximum(x, 0.0)
    xn = x ** n
    return xn / (1.0 + xn)


def free_e_protein(id3, kappa):
    """Free E-protein fraction at binding equilibrium. **Exact, closed form.**

    1:1 titration of E by ID3 with dissociation constant `K_d`, in units of
    `E_tot` (so `E_tot = 1` and `kappa = K_d/E_tot`):

        E_free + ID3_free <-> complex,   complex = E_free*ID3_free/K_d
        1 = E_free + complex,  id3 = ID3_free + complex

    which is a quadratic with the standard root

        complex = 2*id3 / (A + sqrt(A^2 - 4*id3)),   A = 1 + id3 + kappa
        E_free  = 1 - complex

    **Why this form of the root and not `(A - sqrt(A^2 - 4*id3))/2`.** The two
    are algebraically identical, but for large `kappa` the second subtracts two
    nearly equal large numbers and loses all precision — catastrophic
    cancellation in exactly the loose-binding regime. The form above has no
    subtraction of comparable quantities anywhere.

    **Why exact rather than `1 - id3` (decision 006 amendment).** The linear form
    is this expression's *tight-binding limit*, not an independent choice, and it
    goes negative once `id3 > 1`. Clipping with `max(0, .)` hides a domain
    violation rather than fixing one. Check the two limits directly:

        kappa -> 0     complex -> min(1, id3),  E_free -> max(0, 1 - id3)
        kappa -> inf   complex -> id3/kappa,    E_free -> 1/(1 + id3/kappa)

    so the exact form reproduces the linear tight limit *including* its floor,
    without ever leaving its domain.

    Note what this removes: the 11-state model needed a numerical root find with
    an escalation ladder and a convergence ledger here. A closed form cannot
    fail to converge, so the standing "never silently drop a failed solve" rule
    no longer has anything to catch at this step. It still applies to
    continuation and to Phase 3's Monte Carlo.
    """
    id3 = np.maximum(id3, 0.0)
    a = 1.0 + id3 + kappa
    disc = np.maximum(a * a - 4.0 * id3, 0.0)   # >= 0 analytically; guards roundoff
    complexed = 2.0 * id3 / (a + np.sqrt(disc))
    return 1.0 - complexed


#: Prefactor in the tight-binding ultrasensitivity law for THIS mechanism.
#: **Measured on the shipped closed form, not inherited.** See `n_eff`.
_N_EFF_PREFACTOR = 0.5


def n_eff(kappa):
    """Sharpness of `free_e_protein` — `max |d ln E_free / d ln id3|`.

    `n_eff ~= 0.5/sqrt(kappa)` in the tight regime, floored at 1 when binding is
    loose (`E_free` then falls as `1/(1 + id3/kappa)`, log-log slope -1).

    **The 0.5 was measured against the function shipped here, and it is NOT the
    1.34 from `docs/derivations/binding_polynomial.md` §6.** That constant was
    derived for the deleted 11-state model's *ternary complex* `C_L` under
    *two-target* titration — ID3 taxing both PTF1A and E-protein, so its log-log
    slope compounded. This core titrates one target (E) and the diagnostic is on
    `E_free` itself, giving 0.5. Carrying 1.34 across that change of mechanism
    would overstate the threshold sharpness by ~2.7x.

    **Why that would have mattered, quantitatively.** Phase 3 convolves the
    per-cell LNP dose distribution against the bootstrap threshold; a sharper
    threshold produces a larger co-formulation gap, which is the headline number.
    An inherited constant would have inflated R1 by a factor nobody would have
    seen, because the model would still run and the figure would still render.

    **Open for Phase 3:** the quantity actually convolved is the *complex*
    threshold `P*E_free*R`, not `E_free` alone, and its sharpness compounds the
    Hill exponents `n_P`/`n_R` on top of this. Phase 3 must define sharpness on
    the quantity it convolves and measure it there, rather than reusing this.
    Flagged here so the substitution is a decision rather than an accident.
    """
    return float(np.maximum(1.0, _N_EFF_PREFACTOR / np.sqrt(kappa)))


def chop_flag(p_value, threshold=0.05):
    """True when predicted PTF1A is low enough to risk CHOP-dependent apoptosis.

    A **warning on the output, not a term in the right-hand side**. Decision 008
    is retired: viability is measured at the bench, and modelling something you
    can measure using six parameters you cannot is backwards. But the low-PTF1A
    arm is real biology (Sakikubo 2018, PMID 30361559; Backx 2021, PMID
    33762742), so a recommendation that quietly relies on driving PTF1A through
    the floor should say so out loud.
    """
    return bool(np.any(np.asarray(p_value) < threshold))


# ---------------------------------------------------------------------------
# The right-hand side
# ---------------------------------------------------------------------------


def rhs(tau, y, p: Params, inp: Inputs):
    """dy/dtau for (P, R, C). Nondimensional throughout.

        E       = E_free(id3)                    exact binding, see above
        PTF1_J  = P * E                          PTF1A + E-protein + RBPJ
        PTF1_L  = P * E * R                      PTF1A + E-protein + RBPJL

        dP/dtau = a_P*Hill(PTF1_L, n_P) / (1 + (C/c_rep)^n_C)
                  + b_P/(1 + id3)                ERK-suppressed basal ignition
                  - P + u_P

        dR/dtau = rho*(a_R*Hill(PTF1_J, n_R) - R) + u_R       NO P-INDEPENDENT TERM

        dC/dtau = gamma*(alpha_C*id3/(k_w + id3) - C + eps*Hill(C, 2))

    **The two complexes are different, and collapsing them breaks the model.**
    *Rbpjl* is driven by the PTF1-**J** complex, which uses **RBPJ** — broadly
    expressed, not the bottleneck, so its constant pool is absorbed into `a_R`
    and costs no parameter. The 2.3-kb *Ptf1a* autoregulatory enhancer needs
    PTF1-**L**, which uses **RBPJL**. That asymmetry is the developmental handoff
    (Masui 2007) and it is what makes the loop a loop.

    Writing `dR/dtau` against `PTF1_L` instead — i.e. making RBPJL production
    require RBPJL — was the first draft here, and it is wrong in a way that
    passes every guard test: `R = 0` becomes **absorbing**, so no amount of MEK
    inhibition can ever restore RBPJL and the model predicts trametinib alone
    never reverts anything. That contradicts Collins 2014 head-on and would
    destroy Phase 7's internal positive control. The bootstrap claim is *"nothing
    but PTF1A makes RBPJL"*, **not** *"nothing but RBPJL makes RBPJL"* — the
    second is a stronger claim the biology does not support.

    Notes on three further choices that are easy to get wrong:

    **Basal ignition `b_P/(1 + id3)` is required, and its absence would be a
    silent bug.** Without it, PTF1A can only be made by the loop, so `P = 0` is
    absorbing and the payload becomes the *only* route back. The model would then
    predict that trametinib alone can never revert anything — contradicting
    Collins 2014, and destroying Phase 7's internal positive control (the
    trametinib-alone arm is what makes the one-shot experiment un-failable). Its
    ERK sensitivity is tied to the ID3 scale rather than given its own half-max:
    a declared simplification, since both are direct responses to the same "ERK
    is high" state and nothing measured separates their EC50s.

    **Repression by C multiplies rather than competing inside the Hill term.**
    C sits at *metaplasia* loci and cannot gate the acinar enhancer directly, so
    the coupling is programmatic, not local — decision 004, which survives.

    **C is erased passively, at rate gamma.** No active erasure by the acinar
    complex. That keeps time-to-relapse a property of `gamma` alone, which is
    what makes `gamma` readable as the durability knob. If Gate B or Phase 5
    needs active erasure, add it deliberately and say why.
    """
    P = max(float(y[0]), 0.0)
    R = max(float(y[1]), 0.0)
    C = max(float(y[2]), 0.0)

    id3, u_P, u_R = inp.eval(tau)

    e_free = free_e_protein(id3, p.kappa)
    ptf1_j = P * e_free                 # with RBPJ — constant pool, absorbed into a_R
    ptf1_l = P * e_free * R             # with RBPJL — the autoregulatory complex

    # --- P. Self-sustaining production needs PTF1-L, hence needs R: that is the
    # loop. Basal ignition is the only P source when the loop is open.
    memory_repression = 1.0 / (1.0 + (C / p.c_rep) ** p.n_C)
    ignition = p.b_P / (1.0 + id3)      # ERK high -> ignition shut down (cut #1)
    dP = (p.a_P * hill(ptf1_l, p.n_P) * memory_repression
          + ignition
          - P
          + u_P)

    # --- R.  NO P-INDEPENDENT TERM. `hill` is exactly 0 at ptf1_j = 0, so with
    # P = 0 the only contributions are decay and the delivered payload. Note it
    # is driven by PTF1-J, NOT PTF1-L — see the docstring; using PTF1-L here
    # makes R = 0 absorbing and contradicts Collins 2014.
    # tests/test_bootstrap_guard.py exists to keep both properties.
    dR = p.rho * (p.a_R * hill(ptf1_j, p.n_R) - R) + u_R

    # --- C
    write = p.alpha_C * id3 / (p.k_w + id3)
    dC = p.gamma * (write - C + p.eps * hill(C, SELF_REINFORCEMENT_EXPONENT))

    return np.array([dP, dR, dC])


def rhs_for_solver(tau, y, p: Params, inp: Inputs):
    """scipy.integrate signature. Present so callers never index STATES wrongly."""
    return rhs(tau, y, p, inp)
