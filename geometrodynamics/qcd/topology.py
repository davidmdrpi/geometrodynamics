"""
Topology constructors for hadronic networks.

Each function builds a specific topology (meson, baryon, glueball,
hybrid, Möbius, tetraquark) as a HadronicNetwork ready for simulation.
"""

from math import pi

import numpy as np

from geometrodynamics.qcd.constants import SAT_HYBRID_SCALE
from geometrodynamics.qcd.network import (
    Node,
    Branch,
    Junction,
    HadronicNetwork,
)
from geometrodynamics.qcd.bridge import (
    BridgeField,
    make_cornell_branch_potential,
)

_U0 = lambda s, L: 0.0


def _br(branch_id, nl, nr, L, v, lam, cpair, ori, tau, tkind, U, role,
         rs, rd, ftag, bridge=None, **extra):
    """Branch constructor shorthand."""
    return Branch(
        branch_id=branch_id, node_left=nl, node_right=nr, length=L,
        v_speed=v, lambda_nl=lam, color_pair=cpair, orientation=ori,
        tau_delay=tau, topology_kind=tkind, U_potential=U, role=role,
        rep_src=rs, rep_dst=rd, flux_tag=ftag, bridge=bridge, **extra,
    )


def make_meson_tube(L, v=1.0, N=200, dt=0.005, bridge_g=0.0):
    net = HadronicNetwork("meson", N_grid=N, dt=dt)
    net.nodes = {0: Node(0, True, "r"), 1: Node(1, True, "r̄")}
    bf = BridgeField(g=bridge_g) if bridge_g > 0 else None
    U_cornell = make_cornell_branch_potential(L)
    net.branches = {
        0: _br(0, 0, 1, L, v, 0.02, ("r", "r̄"), +1, L / v, "open",
               U_cornell, "quark_arm", "fund", "anti", "rr̄", bf),
    }
    net.junctions = {}
    net.initialize_fields()
    return net


def make_baryon_y_network(L1, L2, L3, v=1.0, N=200, dt=0.005):
    net = HadronicNetwork("baryon", N_grid=N, dt=dt)
    net.nodes = {
        0: Node(0, True, "r"), 1: Node(1, True, "g"),
        2: Node(2, True, "b"), 3: Node(3, False, None),
    }
    net.branches = {
        0: _br(0, 0, 3, L1, v, 0.02, ("r", "g"), +1, L1 / v, "open",
               _U0, "quark_arm", "fund", "fund", "rg"),
        1: _br(1, 1, 3, L2, v, 0.02, ("g", "b"), +1, L2 / v, "open",
               _U0, "quark_arm", "fund", "fund", "gb"),
        2: _br(2, 2, 3, L3, v, 0.02, ("b", "r"), -1, L3 / v, "open",
               _U0, "quark_arm", "fund", "fund", "br"),
    }
    net.junctions = {0: Junction(0, 3, (0, 1, 2), (+1, +1, +1), 0.05)}
    net.initialize_fields()
    return net


def make_glueball_ring(C, v=1.0, N=200, dt=0.005):
    net = HadronicNetwork("glueball", N_grid=N, dt=dt)
    net.nodes = {0: Node(0, False, None)}
    net.branches = {
        0: _br(0, 0, 0, C, v, 0.01, ("r", "r̄"), +1, C / (2 * v),
               "periodic", _U0, "ring", "fund", "anti", "ring"),
    }
    net.junctions = {}
    net.initialize_fields()
    return net


def make_mobius_tube(L, v=1.0, N=200, dt=0.005):
    net = HadronicNetwork("mobius", N_grid=N, dt=dt)
    net.nodes = {0: Node(0, True, "r"), 1: Node(1, False, None)}
    net.branches = {
        0: _br(0, 0, 1, L, v, 0.02, ("r", "r̄"), -1, 2 * L / v,
               "mobius", _U0, "mobius_arm", "fund", "anti", "mobius"),
    }
    net.junctions = {}
    net.initialize_fields()
    return net


def make_hybrid_excitation(
    meson_length, ring_circumference, v=1.0, N=200, dt=0.005,
):
    """Hybrid: two meson arms + attached_loop ring at junction."""
    La = meson_length / 2.0
    omega_loop = 2.0 * pi * v / ring_circumference
    L_loop_val = ring_circumference / (2.0 * pi * v)
    kappa_default = 0.60 * L_loop_val * omega_loop ** 2

    amp_seed = 0.05
    I_loop_seed = amp_seed / (L_loop_val + 1e-30)

    net = HadronicNetwork("hybrid", N_grid=N, dt=dt)
    net.nodes = {
        0: Node(0, True, "r"), 1: Node(1, True, "r̄"),
        2: Node(2, False, None),
    }
    net.branches = {
        0: _br(0, 0, 2, La, v, 0.02, ("r", "r̄"), +1, La / v, "open",
               _U0, "quark_arm", "fund", "anti", "rr̄_L",
               sat_scale=SAT_HYBRID_SCALE),
        1: _br(1, 1, 2, La, v, 0.02, ("r̄", "r"), -1, La / v, "open",
               _U0, "quark_arm", "anti", "fund", "rr̄_R",
               sat_scale=SAT_HYBRID_SCALE),
        2: _br(2, 2, 2, ring_circumference, v, 0.01, ("r", "r̄"), +1,
               ring_circumference / (2 * v), "attached_loop", _U0, "ring",
               "fund", "anti", "ring_hyb",
               kappa_loop=kappa_default, I_loop=I_loop_seed,
               loop_couple=0.40, sat_scale=SAT_HYBRID_SCALE),
    }
    net.junctions = {0: Junction(0, 2, (0, 1, 2), (+1, +1, +1), 0.03)}
    net.initialize_fields()
    return net


def make_tetraquark_double_y(
    L_arm=0.5, L_center=1.0, v=1.0, N=200, dt=0.005, bridge_g=0.0,
):
    """Tetraquark: two Y-junctions connected by a central arm."""
    net = HadronicNetwork("tetraquark", N_grid=N, dt=dt)
    net.nodes = {
        0: Node(0, True, "r"), 1: Node(1, True, "g"),
        2: Node(2, True, "r̄"), 3: Node(3, True, "ḡ"),
        4: Node(4, False, None), 5: Node(5, False, None),
    }
    bf = BridgeField(g=bridge_g) if bridge_g > 0 else None
    U_arm = make_cornell_branch_potential(L_arm)
    U_center = make_cornell_branch_potential(L_center)
    net.branches = {
        0: _br(0, 0, 4, L_arm, v, 0.02, ("r", "r̄"), +1, L_arm / v,
               "open", U_arm, "quark_arm", "fund", "anti", "r_left"),
        1: _br(1, 1, 4, L_arm, v, 0.02, ("g", "ḡ"), +1, L_arm / v,
               "open", U_arm, "quark_arm", "fund", "anti", "g_left"),
        2: _br(2, 4, 5, L_center, v, 0.02, ("r", "r̄"), +1, L_center / v,
               "open", U_center, "junction_bridge", "fund", "anti",
               "center", bf),
        3: _br(3, 5, 2, L_arm, v, 0.02, ("r̄", "r"), -1, L_arm / v,
               "open", U_arm, "quark_arm", "anti", "fund", "r_right"),
        4: _br(4, 5, 3, L_arm, v, 0.02, ("ḡ", "g"), -1, L_arm / v,
               "open", U_arm, "quark_arm", "anti", "fund", "g_right"),
    }
    net.junctions = {
        0: Junction(0, 4, (0, 1, 2), (+1, +1, -1), 0.03),
        1: Junction(1, 5, (2, 3, 4), (+1, +1, +1), 0.03),
    }
    net.initialize_fields()
    return net


def make_mobius_baryon_y_network(L1, L2, L3, v=1.0, N=200, dt=0.005):
    """Y-network with one Möbius arm and two ordinary arms."""
    net = HadronicNetwork("mobius_baryon", N_grid=N, dt=dt)
    net.nodes = {
        0: Node(0, True, "r"), 1: Node(1, True, "g"),
        2: Node(2, True, "b"), 3: Node(3, False, None),
    }
    net.branches = {
        0: _br(0, 0, 3, L1, v, 0.02, ("r", "g"), +1, L1 / v, "open",
               _U0, "quark_arm", "fund", "fund", "rg"),
        1: _br(1, 1, 3, L2, v, 0.02, ("g", "b"), -1, 2 * L2 / v,
               "mobius", _U0, "mobius_arm", "fund", "fund", "gb_mob"),
        2: _br(2, 2, 3, L3, v, 0.02, ("b", "r"), -1, L3 / v, "open",
               _U0, "quark_arm", "fund", "fund", "br"),
    }
    net.junctions = {0: Junction(0, 3, (0, 1, 2), (+1, +1, +1), 0.05)}
    net.initialize_fields()
    return net


def make_mobius_baryon_v12(L1=0.5, L2=0.5, L3=0.5, v=1.0, N=80, dt=0.003):
    """Möbius baryon with dedicated Kirchhoff coupling."""
    net = HadronicNetwork("mobius_baryon", N_grid=N, dt=dt)
    net.nodes = {
        0: Node(0, True, "r"), 1: Node(1, True, "g"),
        2: Node(2, True, "b"), 3: Node(3, False, None),
    }
    net.branches = {
        0: _br(0, 0, 3, L1, v, 0.02, ("r", "g"), +1, L1 / v, "open",
               _U0, "quark_arm", "fund", "fund", "rg"),
        1: _br(1, 1, 3, L2, v, 0.02, ("g", "b"), -1, 2 * L2 / v,
               "mobius", _U0, "mobius_arm", "fund", "fund", "gb_mob",
               kappa_mobius=0.30),
        2: _br(2, 2, 3, L3, v, 0.02, ("b", "r"), -1, L3 / v, "open",
               _U0, "quark_arm", "fund", "fund", "br"),
    }
    net.junctions = {0: Junction(0, 3, (0, 1, 2), (+1, +1, +1), 0.01)}
    net.initialize_fields()
    return net
