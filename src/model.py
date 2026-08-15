"""The full unreduced ADM system — Stage 0, Step 1.

Nothing is reduced here. Nondimensionalisation and quasi-steady-state
elimination are the next session; you cannot eliminate what you have not
written. The one exception is the binding equilibrium (see binding.py), which
was agreed as the treatment for ID3 sequestration rather than being a reduction
performed here.

State vector layout is defined by topology.Topology.states and is fixed.

READ THIS BEFORE EDITING:
  * W is PROTECTED FROM ELIMINATION (decision 002). QSS on W recovers the static
    product K_eff = K*f_act*f_cat, which provably cannot produce the §1.2
    withdrawal-asymmetry prediction.
  * RAF_drive is strictly DECREASING in K_eff. Inverting it flips the mechanism
    and the model still runs.
  * dR/dt has NO ignition term. Rbpjl has no PTF1A-independent promoter; that
    asymmetry is the structural reason RBPJL is the deepest hole in the loop,
    and it is constraint #1 of the four non-negotiables in §3.4.
"""

from __future__ import annotations

import numpy as np

from . import binding, functional_forms as ff
from .topology import Topology


def rhs(t, y, p, topo: Topology, warm=None):
    """Right-hand side of the full system.

    `p` is a flat parameter dict; `topo` selects which optional terms exist.
    `warm` carries the previous binding solution to warm-start the equilibrium
    solve — the RHS is called many times per step and a cold solve each time is
    the difference between minutes and hours.
    """
    idx = topo.index
    P_n   = max(y[idx("P_n")],   0.0)
    R     = max(y[idx("R")],     0.0)
    E_tot = max(y[idx("E_tot")], 0.0)
    I     = max(y[idx("I")],     0.0)
    M     = max(y[idx("M")],     0.0)
    A     = max(y[idx("A")],     0.0)
    S     = max(y[idx("S")],     0.0)
    W     = max(y[idx("W")],     0.0)

    MIST1 = max(y[idx("MIST1")], 0.0) if topo.mist1_arm else 0.0
    NR5A2 = max(y[idx("NR5A2")], 0.0) if topo.nr5a2_mode != "absent" else 0.0

    # ---------------------------------------------------------------- inputs
    v = _trametinib(t, p) if topo.use_trametinib else 0.0
    u1 = _pulse(t, p, "u1") if topo.use_u1_rbpjl else 0.0   # RBPJL mRNA
    u2 = _pulse(t, p, "u2") if topo.use_u2_ptf1a else 0.0   # PTF1A mRNA
    u3 = _pulse(t, p, "u3") if topo.use_u3 else 0.0         # third mRNA

    # u3 routes to a DIFFERENT pool depending on identity. This is what makes
    # the Stage 3 ranking a real comparison rather than a relabelling: E47
    # helps by relieving titration, NR5A2 by boosting transcription, MIST1 by
    # raising secretory capacity. Three mechanisms, one payload slot.
    u3_to_E = u3 if topo.u3_identity == "E47" else 0.0
    u3_to_NR5A2 = u3 if topo.u3_identity == "NR5A2" else 0.0
    u3_to_MIST1 = u3 if topo.u3_identity == "MIST1" else 0.0

    # ------------------------------------------------- ERK / drug layer (W)
    # K_eff is what every downstream term sees. It carries the drug, so no
    # other term may take v directly — that would double-count the drug.
    K_eff = p["K_kras"] * W * ff.f_cat(v, p["ic50_cat"])

    drive = ff.raf_drive(K_eff, p["k_fb"], p["r_fb"])       # DECREASING in K_eff
    dW = (p["k_on_W"] * drive * ff.f_act(v, p["ic50_act"], topo.is_trametinib)
          - p["k_off_W"] * W)

    # --------------------------------------------------- binding equilibrium
    solver = (binding.solve_binding if topo.id3_mode == "titration"
              else binding.first_order_sequestration)
    b = solver(P_n, E_tot, I, R, p["kd"], p["rbpj_total"], guess=warm)
    C_L, C_J, P_free = b["C_L"], b["C_J"], b["P_free"]

    # T3: RBPJ->RBPJL handoff as the memory. The immature complex competes with
    # the adult one, so a cell stuck with PTF1-J cannot express the adult
    # program even with PTF1A present.
    handoff = ff.hill_repress(C_J, p["k_handoff"], p["n_handoff"]) if topo.rbpj_handoff else 1.0

    # ------------------------------------------------------------ PTF1A (P_n)
    # Two sources: the PTF1-independent ignition promoter (ERK-suppressed) and
    # the 2.3-kb autoregulatory enhancer (needs the complex, opposed by the
    # metaplasia program).
    enhancer = ff.hill_activate(C_L, p["k_auto"], p["n_auto"])
    if topo.chromatin_memory:
        # M sits at METAPLASIA loci, so it cannot gate the acinar enhancer
        # directly. It acts indirectly: the metaplasia program represses the
        # acinar program. Same effect, correct location (constraint #3).
        enhancer *= ff.hill_repress(M, p["k_M_rep"], p["n_M_rep"])
    if topo.nr5a2_mode == "enhancer":
        # T6a. ASSUMPTION under test — Holmstrom 2011 shows an exocrine network,
        # not Ptf1a enhancer binding. T6b is the alternative. See decision 010.
        enhancer *= (1.0 + p["kappa_nr5a2"] * NR5A2)

    dP_n = (p["alpha_ign"] * ff.g_ignition(K_eff, p["k_ign"], p["q_ign"])
            + p["alpha_auto"] * enhancer * handoff
            - p["delta_P"] * P_n
            + u2)
    # NOTE: no -k_seq*I*P_n term. Sequestration is in the binding equilibrium.

    # ------------------------------------------------------------ RBPJL (R)
    # NO IGNITION TERM. Constraint #1. Rbpjl's only driver is PTF1 itself, so
    # it cannot bootstrap — the loop cannot re-close even if PTF1A returns.
    dR = (p["beta_R"] * ff.hill_activate(C_L, p["k_R"], p["m_R"])
          - p["delta_R"] * R
          + u1)

    # ------------------------------------------------------- E-protein (E_tot)
    # TOTAL pool: slow synthesis/degradation plus the u3=E47 input. Free E is
    # derived from the polynomial, never integrated. You cannot eliminate the
    # pool the payload doses (B3).
    dE_tot = p["beta_E"] - p["delta_E"] * E_tot + u3_to_E

    # ---------------------------------------------------------------- ID3 (I)
    # ERK-driven. Intermediate timescale, and the leading QSS candidate at the
    # reduction step — which is how the slow count comes down to 6.
    dI = (p["beta_I0"]
          + p["beta_I"] * ff.hill_activate(K_eff, p["k_I"], p["n_I"])
          - p["delta_I"] * I)

    # ------------------------------------------------- chromatin memory (M)
    # Slow, self-reinforcing, at METAPLASIA loci. Written by ERK, opposed by
    # the acinar complex, plus a self-reinforcement term that is what makes it
    # a memory rather than a filter — and what sets relapse timing after the
    # payload clears.
    if topo.chromatin_memory:
        write = (ff.hill_activate(K_eff, p["k_wM"], p["a_wM"])
                 * ff.hill_repress(C_L, p["k_CM"], p["b_CM"]))
        dM = (p["k_w"] * write - p["k_e"] * M
              + p["eps_M"] * ff.hill_activate(M, p["theta_M"], 2))
    else:
        dM = -p["k_e"] * M

    # ------------------------------------------------------ acinar output (A)
    # The validation observable. Collins reports "% amylase-positive cells", so
    # A needs an explicit map to that (see observables.py) or the comparison to
    # published numbers is meaningless.
    output_gain = 1.0
    if topo.nr5a2_mode == "output":
        output_gain += p["kappa_nr5a2"] * NR5A2      # T6b alternative placement
    dA = p["k_A"] * output_gain * ff.hill_activate(C_L, p["k_AC"], p["n_AC"]) - p["delta_A"] * A

    # --------------------------------------------- secretory load / capacity
    # Capacity is co-induced with the differentiation program, which is why a
    # one-sided ceiling does not exist at steady state. In T5 that co-induction
    # is carried by an explicit MIST1 state so it can be varied independently;
    # elsewhere it is a constant.
    gamma = p["gamma0"] + (p["gamma1"] * MIST1 if topo.mist1_arm else 0.0)
    dS = p["k_cargo"] * A - gamma * S

    # ------------------------------------------------------- variant states
    dMIST1 = None
    if topo.mist1_arm:
        dMIST1 = (p["beta_M1"] * ff.hill_activate(C_L, p["k_M1"], p["n_M1"])
                  - p["delta_M1"] * MIST1 + u3_to_MIST1)
    dNR5A2 = None
    if topo.nr5a2_mode != "absent":
        dNR5A2 = (p["beta_N"] * ff.hill_activate(C_L, p["k_N"], p["n_N"])
                  - p["delta_N"] * NR5A2 + u3_to_NR5A2)

    # ------------------------------------------------------------- viability
    # Accumulated hazard, so survival = exp(-cumhaz) is read straight off the
    # trajectory. U-shaped: high S (cargo outruns capacity) AND low free PTF1A
    # (CHOP-dependent apoptosis). Never a threshold — that was U_crit.
    dcumhaz = ff.hazard(S, P_free, p)

    out = np.zeros_like(y)
    out[idx("P_n")] = dP_n
    out[idx("R")] = dR
    out[idx("E_tot")] = dE_tot
    out[idx("I")] = dI
    out[idx("M")] = dM
    out[idx("A")] = dA
    out[idx("S")] = dS
    out[idx("W")] = dW
    out[idx("cumhaz")] = dcumhaz
    if dMIST1 is not None:
        out[idx("MIST1")] = dMIST1
    if dNR5A2 is not None:
        out[idx("NR5A2")] = dNR5A2
    return out, b


def rhs_for_solver(t, y, p, topo, cache):
    """scipy.integrate wrapper: returns only dy/dt, keeps the binding warm start."""
    dy, b = rhs(t, y, p, topo, warm=cache.get("guess"))
    cache["guess"] = b.get("_guess")
    return dy


def _trametinib(t, p):
    """Drug exposure. Constant hold for now; Stage 5 replaces this with a
    schedule. Kept as a function so that swap costs nothing."""
    return p.get("v_dose", 0.0)


def _pulse(t, p, which):
    return ff.mrna_pulse(
        t,
        dose=p[f"{which}_dose"],
        interval=p[f"{which}_interval"],
        k_translate=p[f"{which}_k_translate"],
        delta_m=p[f"{which}_delta_m"],
        n_pulses=p[f"{which}_n_pulses"],
    )
