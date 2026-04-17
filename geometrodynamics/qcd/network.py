"""
Network data structures for the QCD flux-tube model.

A HadronicNetwork consists of Nodes (quarks / junction points),
Branches (flux tubes), and Junctions (where branches meet with
Kirchhoff-like current conservation).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import pi
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from geometrodynamics.qcd.constants import (
    MU0_CONF_SQ,
    MUJ_CONF_SQ,
    TopologyKind,
    TopologyKindBranch,
    BranchRole,
)
from geometrodynamics.qcd.color import (
    ANTI_OF,
    color_charge_vector,
    is_singlet,
)
from geometrodynamics.qcd.bridge import BridgeField


@dataclass
class Node:
    node_id: int
    endpoint: bool
    color: Optional[str] = None


@dataclass
class Branch:
    branch_id: int
    node_left: int
    node_right: int
    length: float
    v_speed: float = 1.0
    lambda_nl: float = 0.02
    color_pair: Tuple[str, str] = ("r", "r̄")
    orientation: int = +1
    tau_delay: float = 0.0
    topology_kind: TopologyKindBranch = "open"
    attach_end: str = "left"
    U_potential: Callable = field(default=lambda s, L: 0.0)
    role: BranchRole = "quark_arm"
    rep_src: str = "fund"
    rep_dst: str = "anti"
    flux_tag: str = ""
    bridge: Optional[BridgeField] = field(default=None, repr=False)
    active: bool = True
    psi: Optional[np.ndarray] = field(default=None, repr=False)
    psi_old: Optional[np.ndarray] = field(default=None, repr=False)

    # Loop-current state for attached_loop branches
    I_loop: float = 0.0
    Phi_loop: float = 0.0
    kappa_loop: float = 0.0
    loop_couple: float = 0.18
    kappa_mobius: float = 0.22
    sat_scale: float = 1.0

    @property
    def periodic(self) -> bool:
        return self.topology_kind == "periodic"

    @property
    def mobius(self) -> bool:
        return self.topology_kind == "mobius"

    @property
    def attached_ring(self) -> bool:
        return self.topology_kind == "attached_loop"

    @property
    def L_loop_inductance(self) -> float:
        if self.attached_ring and self.v_speed > 0:
            return self.length / (2.0 * pi * self.v_speed)
        return float("inf")


@dataclass
class Junction:
    junction_id: int
    node_id: int
    branch_ids: Tuple[int, ...]
    side_signs: Tuple[int, ...]
    J_nonlinear: float = 0.05


@dataclass
class HadronicNetwork:
    topology: TopologyKind
    nodes: Dict[int, Node] = field(default_factory=dict)
    branches: Dict[int, Branch] = field(default_factory=dict)
    junctions: Dict[int, Junction] = field(default_factory=dict)
    N_grid: int = 200
    dt: float = 0.005

    def initialize_fields(
        self, psi0: Optional[Dict[int, np.ndarray]] = None
    ) -> None:
        for bid, b in self.branches.items():
            N = self.N_grid
            if psi0 and bid in psi0:
                p = psi0[bid].copy()
            else:
                p = np.zeros(N)
            b.psi = p.copy()
            b.psi_old = p.copy()

    def endpoint_colors(self) -> List[str]:
        return [n.color for n in self.nodes.values() if n.endpoint and n.color]

    def is_color_singlet(self) -> bool:
        return is_singlet(self.endpoint_colors())

    def branch_color_charge(self, bid: int) -> np.ndarray:
        b = self.branches[bid]
        return color_charge_vector([b.color_pair[0]]) + color_charge_vector(
            [b.color_pair[1]]
        )

    def junction_color_charge(self, jid: int) -> np.ndarray:
        junc = self.junctions[jid]
        q = np.zeros(3)
        for bid, sgn in zip(junc.branch_ids, junc.side_signs):
            b = self.branches[bid]
            if b.periodic:
                continue
            c = b.color_pair[1] if sgn == +1 else b.color_pair[0]
            q += color_charge_vector([c])
        return q

    def confinement_mass_sq(self, bid: int) -> float:
        b = self.branches[bid]
        if b.periodic:
            return 0.0
        dQ = self.branch_color_charge(bid)
        mu2 = MU0_CONF_SQ * float(np.dot(dQ, dQ))
        for junc in self.junctions.values():
            if bid in junc.branch_ids:
                dQj = self.junction_color_charge(junc.junction_id)
                mu2 += MUJ_CONF_SQ * float(np.dot(dQj, dQj))
                break
        return mu2

    def nucleate_pair_on_branch(self, bid: int) -> Tuple[bool, str]:
        """Split branch at mid-point into two singlet daughter mesons."""
        if bid not in self.branches:
            return False, "branch not found"
        parent = self.branches[bid]
        if not parent.active:
            return False, "branch already inactive"

        N = self.N_grid
        L = parent.length
        L2 = L / 2.0
        c_q = parent.color_pair[0]
        c_qbar = parent.color_pair[1]
        c_new = ANTI_OF.get(c_q, "r̄")
        c_new_q = ANTI_OF.get(c_new, "r")

        left_colors = [c_q, c_new]
        right_colors = [c_new_q, c_qbar]
        if not is_singlet(left_colors) or not is_singlet(right_colors):
            return False, f"daughter colors not singlet: {left_colors}, {right_colors}"

        max_nid = max(self.nodes.keys()) if self.nodes else -1
        max_bid = max(self.branches.keys()) if self.branches else -1
        nid_qnew = max_nid + 1
        nid_qbnew = max_nid + 2
        bid_L = max_bid + 1
        bid_R = max_bid + 2

        self.nodes[nid_qnew] = Node(nid_qnew, True, c_new_q)
        self.nodes[nid_qbnew] = Node(nid_qbnew, True, c_new)

        psi_now = parent.psi if parent.psi is not None else np.zeros(N)
        psi_L = psi_now[: N // 2 + 1]
        psi_R = psi_now[N // 2 :]

        def _pad(arr, tgt):
            if len(arr) == tgt:
                return arr.copy()
            return np.interp(
                np.linspace(0, 1, tgt), np.linspace(0, 1, len(arr)), arr
            )

        self.branches[bid_L] = Branch(
            branch_id=bid_L,
            node_left=parent.node_left,
            node_right=nid_qbnew,
            length=L2,
            v_speed=parent.v_speed,
            lambda_nl=parent.lambda_nl,
            color_pair=(c_q, c_new),
            orientation=parent.orientation,
            tau_delay=L2 / parent.v_speed,
            topology_kind="open",
            U_potential=parent.U_potential,
            role="quark_arm",
            rep_src=parent.rep_src,
            rep_dst="anti",
            flux_tag=f"{c_q}{c_new}_L",
            bridge=None,
            psi=_pad(psi_L, N),
            psi_old=_pad(psi_L, N),
        )
        self.branches[bid_R] = Branch(
            branch_id=bid_R,
            node_left=nid_qnew,
            node_right=parent.node_right,
            length=L2,
            v_speed=parent.v_speed,
            lambda_nl=parent.lambda_nl,
            color_pair=(c_new_q, c_qbar),
            orientation=parent.orientation,
            tau_delay=L2 / parent.v_speed,
            topology_kind="open",
            U_potential=parent.U_potential,
            role="quark_arm",
            rep_src="fund",
            rep_dst=parent.rep_dst,
            flux_tag=f"{c_new_q}{c_qbar}_R",
            bridge=None,
            psi=_pad(psi_R, N),
            psi_old=_pad(psi_R, N),
        )
        parent.active = False
        return True, "ok"
