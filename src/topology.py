"""Topology configuration — a topology is a config object, never a source file.

Decision 003. Stage 2 compares five architectures by Q-value under an identical
sampling box (Ma et al., Cell 2009, PMID 19703401), and that comparison is only
model selection if *nothing* differs except the right-hand side. Five copied
files drift silently. A single maximal model with unused terms zeroed is also
wrong, more subtly: a zeroed state still enters the Jacobian, still costs
continuation dimensions, and still shows up in Stage 3's FIM spectrum as a
spurious sloppy direction — so model selection would run at a dimension no
candidate topology actually has.

Hence: one base right-hand side, optional terms selected here, and the active
state vector assembled per topology.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# CORE states are present in every topology. Order is fixed and load-bearing:
# it defines the layout of the state vector everywhere downstream, so it must
# never be reordered casually.
CORE_STATES: tuple[str, ...] = (
    "P_n",      # total nuclear PTF1A (free/ID3-bound/complexed all derived)
    "R",        # RBPJL protein — the bottleneck, no independent promoter
    "E_tot",    # TOTAL E-protein; free E is derived. u3=E47 doses this pool.
    "I",        # ID3 — the titrator
    "M",        # slow self-reinforcing chromatin at METAPLASIA loci
    "A",        # acinar output (amylase / CPA1 proxy)
    "S",        # secretory cargo / capacity ratio
    "W",        # phospho-MEK pool — PROTECTED FROM ELIMINATION (decision 002)
    "cumhaz",   # ∫h dt, so survival = exp(-cumhaz) is read off directly
)

# Variant states exist only in the topology that tests them. Keeping them out
# of CORE is what holds the reduced system at 6 slow states + W.
VARIANT_STATES: dict[str, str] = {
    "MIST1": "T5 only — secretory capacity arm (Jakubison 2018)",
    "NR5A2": "T6a/T6b only — enhancer vs acinar-output co-activator",
}


@dataclass(frozen=True)
class Topology:
    """Which optional terms are active. Frozen: a topology must not mutate
    mid-run, or the recorded configuration stops describing the result."""

    name: str

    # --- T1 vs T2: how ID3 acts. This is the term that decides whether the
    # model can be bistable at all, so it is the primary axis of competition.
    id3_mode: str = "titration"      # "titration" (T1) | "first_order" (T2)

    # --- T3: RBPJ→RBPJL handoff carries the memory
    rbpj_handoff: bool = False

    # --- T4: slow self-reinforcing chromatin carries the memory
    chromatin_memory: bool = True

    # --- T5: MIST1 parallel arm carries secretory capacity
    mist1_arm: bool = False

    # --- T6a/T6b: NR5A2 placement. This is an ASSUMPTION under test, not a
    # finding — Holmstrom 2011 shows an exocrine network, not Ptf1a enhancer
    # binding. Both placements compete; the Q-values decide. See decision 010.
    nr5a2_mode: str = "absent"       # "absent" | "enhancer" (T6a) | "output" (T6b)

    # --- payload composition. All four independently toggleable from the first
    # line of code, because the necessity analysis needs all 16 subsets and
    # retrofitting this later would retrofit it badly.
    use_trametinib: bool = True
    use_u1_rbpjl: bool = True
    use_u2_ptf1a: bool = True
    use_u3: bool = True
    u3_identity: str = "E47"         # "E47" | "NR5A2" | "MIST1"

    # --- drug identity. f_act ≡ 1 recovers PD325901 and preserves the direct
    # comparison to Collins 2014.
    is_trametinib: bool = True

    def __post_init__(self) -> None:
        # Refuse undeclared combinations rather than silently producing a
        # right-hand side nobody intended — decision 003's first reversal
        # condition. Validation at construction, not at first bad figure.
        if self.id3_mode not in {"titration", "first_order"}:
            raise ValueError(f"unknown id3_mode {self.id3_mode!r}")
        if self.nr5a2_mode not in {"absent", "enhancer", "output"}:
            raise ValueError(f"unknown nr5a2_mode {self.nr5a2_mode!r}")
        if self.u3_identity not in {"E47", "NR5A2", "MIST1"}:
            raise ValueError(f"unknown u3_identity {self.u3_identity!r}")
        # A payload cannot dose a species the topology does not contain. This
        # is exactly the silent-null failure flagged for RBPJL in Stage 3B:
        # the run would complete and score ~0, reading as "MIST1 doesn't help"
        # rather than "MIST1 was not in the model".
        if self.use_u3 and self.u3_identity == "MIST1" and not self.mist1_arm:
            raise ValueError(
                "u3_identity='MIST1' requires mist1_arm=True. Dosing an absent "
                "species scores ~0 and reads as a biological result."
            )
        if self.use_u3 and self.u3_identity == "NR5A2" and self.nr5a2_mode == "absent":
            raise ValueError(
                "u3_identity='NR5A2' requires nr5a2_mode in {'enhancer','output'}. "
                "Same silent-null trap as above."
            )

    @property
    def states(self) -> tuple[str, ...]:
        """Active state vector for this topology, in fixed order."""
        extra: list[str] = []
        if self.mist1_arm:
            extra.append("MIST1")
        if self.nr5a2_mode != "absent":
            extra.append("NR5A2")
        return CORE_STATES + tuple(extra)

    def index(self, name: str) -> int:
        return self.states.index(name)


# The five Stage 2 candidates plus the two NR5A2 placement variants. Stage 2
# swaps between these with one argument, which is the point.
TOPOLOGIES: dict[str, Topology] = {
    "T1": Topology("T1_id3_titration", id3_mode="titration"),
    "T2": Topology("T2_id3_first_order", id3_mode="first_order"),
    "T3": Topology("T3_rbpj_handoff", rbpj_handoff=True),
    "T4": Topology("T4_chromatin_memory", chromatin_memory=True),
    "T5": Topology("T5_mist1_arm", mist1_arm=True),
    "T6a": Topology("T6a_nr5a2_enhancer", nr5a2_mode="enhancer"),
    "T6b": Topology("T6b_nr5a2_output", nr5a2_mode="output"),
}

# CORE with the payload at zero — what Stage 1's continuation runs on.
CORE_UNTREATED = Topology(
    "CORE_untreated",
    use_trametinib=False,
    use_u1_rbpjl=False,
    use_u2_ptf1a=False,
    use_u3=False,
)


def payload_subsets() -> list[dict[str, bool]]:
    """All 16 on/off combinations of the four interventions.

    §1.3's necessity analysis. Reported twice — at matched per-component dose
    AND at matched total delivered dose — because dropping a component changes
    total material, and without the matched-total arm the analysis just
    rediscovers "more protein is better".
    """
    keys = ("use_trametinib", "use_u1_rbpjl", "use_u2_ptf1a", "use_u3")
    subsets = []
    for mask in range(16):
        subsets.append({k: bool(mask >> i & 1) for i, k in enumerate(keys)})
    return subsets
