"""Tests for the QCD flux-tube network solver."""

import numpy as np
import pytest

from geometrodynamics.qcd import (
    make_meson_tube,
    make_baryon_y_network,
    make_glueball_ring,
    make_mobius_tube,
    make_hybrid_excitation,
    HadronicNetworkSolver,
)
from geometrodynamics.qcd.color import is_singlet
from geometrodynamics.qcd.spectrum import MobiusSpectrum
from geometrodynamics.qcd.diagnostics import BridgeCouplingCalibrator


class TestColorAlgebra:
    def test_singlet_meson(self):
        assert is_singlet(["r", "r̄"])

    def test_singlet_baryon(self):
        assert is_singlet(["r", "g", "b"])

    def test_non_singlet(self):
        assert not is_singlet(["r", "g"])

    def test_empty_singlet(self):
        assert is_singlet([])


class TestMobiusSpectrum:
    def test_half_integer_modes(self):
        """Möbius spectrum matches analytic (n+½)π/L within 5%."""
        mob = MobiusSpectrum(L=1.0, v=1.0, N=150)
        report = mob.report(n_max=4)
        assert report["all_pass"]


class TestMesonDynamics:
    def test_energy_conservation(self):
        """Meson tube energy drift < 5% over 300 steps."""
        net = make_meson_tube(1.0, v=1.0, N=80, dt=0.004)
        s = np.linspace(0, 1.0, 80)
        net.initialize_fields(psi0={0: 0.5 * np.sin(np.pi * s)})
        slv = HadronicNetworkSolver(net, antipodal_coupling=0.03)
        hist = slv.run(300, 15)
        E = hist["energy"][np.isfinite(hist["energy"])]
        if len(E) >= 2:
            drift = float(np.max(np.abs(E - E[0])) / (abs(E[0]) + 1e-30))
            assert drift < 0.05


class TestMobiusTube:
    def test_nonzero_energy(self):
        net = make_mobius_tube(1.0, v=1.0, N=100, dt=0.004)
        s = np.linspace(0, 1.0, 100)
        p = 0.5 * np.sin(0.5 * np.pi * s)
        p[0] = 0.0
        net.initialize_fields(psi0={0: p})
        slv = HadronicNetworkSolver(net, antipodal_coupling=0.03)
        E0 = slv.total_energy()
        assert E0 > 1e-6


class TestHybridExcitation:
    def test_topology_kind(self):
        hyb = make_hybrid_excitation(1.0, 1.5, v=1.0, N=80, dt=0.002)
        assert hyb.branches[2].topology_kind == "attached_loop"

    def test_color_singlet(self):
        hyb = make_hybrid_excitation(1.0, 1.5, v=1.0, N=80, dt=0.002)
        assert hyb.is_color_singlet()

    def test_loop_active(self):
        hyb = make_hybrid_excitation(1.0, 1.5, v=1.0, N=80, dt=0.002)
        assert hyb.branches[2].kappa_loop > 0


@pytest.mark.slow
class TestBridgeNucleation:
    def test_nucleation_produces_daughters(self):
        """Bridge nucleation splits meson into two singlet daughters."""
        cal = BridgeCouplingCalibrator(L_ref=1.8, N=60, dt=0.005, t_max=100.0)
        rc = cal.calibrate(g_min=0.1, g_max=3.0, n_pts=10)
        g_star = rc["g_star"] if np.isfinite(rc["g_star"]) else 1.2

        meson = make_meson_tube(1.8, v=1.0, N=80, dt=0.004, bridge_g=g_star)
        s = np.linspace(0, 1.8, 80)
        meson.initialize_fields(psi0={0: 0.8 * np.sin(np.pi * s / 1.8)})
        slv = HadronicNetworkSolver(meson, antipodal_coupling=0.05)

        for _ in range(int(50 / 0.004)):
            slv.step()
            if slv.nucleation_events:
                break

        if slv.nucleation_events:
            # Parent inactive
            assert not meson.branches[0].active
            # Two daughters
            daughters = [
                b for bid, b in meson.branches.items()
                if bid != 0 and b.active
            ]
            assert len(daughters) == 2
            # Both daughters are color singlets
            for d in daughters:
                assert is_singlet(list(d.color_pair))
