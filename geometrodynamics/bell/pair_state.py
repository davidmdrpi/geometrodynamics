"""
Bell pair state — the singlet-like entangled state from throat topology.

A Bell pair is a ConjugatePair (one non-orientable throat) equipped with
cavity memory. The key property: the pair is ONE connected topological
object, not two particles with hidden values.

The Bell layer now uses the real antipodal cavity as a detector-conditioned
history engine. Detector settings feed into local boundary responses, these
responses drive the cavity modes through a time-symmetric emit/confirm
sequence, and the relative 0-branch / π-branch weights are read off from the
resulting packet history plus the residual mode amplitudes.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass

import numpy as np

from geometrodynamics.constants import (
    CAVITY_ALPHA_ADV,
    CAVITY_ALPHA_EMIT,
)
from geometrodynamics.embedding.topology import ConjugatePair, make_singlet_pair
from geometrodynamics.transaction.cavity import AntipodalCavity, make_cavity


@dataclass(frozen=True)
class BranchWeightResult:
    """Detector-conditioned cavity branch weights.

    Parameters
    ----------
    w0, wpi : float
        Normalised 0-branch / π-branch weights.
    packet_weight_0, packet_weight_pi : float
        Weight deposited by discrete packet firings on each branch.
    residual_weight_0, residual_weight_pi : float
        Remaining closure-compatible weight carried by the still-ringing modes.
    n_packets_0, n_packets_pi : int
        Number of fired packets on each branch.
    energy_final : float
        Total cavity energy after the detector-conditioned history evolution.
    phase_spin : float
        Detector-induced SU(2)/Hopf phase offset used in closure checks.
    """

    w0: float
    wpi: float
    packet_weight_0: float
    packet_weight_pi: float
    residual_weight_0: float
    residual_weight_pi: float
    n_packets_0: int
    n_packets_pi: int
    energy_final: float
    phase_spin: float


@dataclass
class BellPair:
    """A Bell-experiment-ready conjugate pair with cavity memory.

    Combines the topological pair state with the persistent cavity that
    determines closure branch weights.
    """

    pair: ConjugatePair
    cavity: AntipodalCavity

    @property
    def tau_semi(self) -> float:
        return self.pair.tau_semi

    def _detector_mode_response(self, theta: float) -> dict[int, float]:
        """Local detector boundary response by cavity channel.

        Derived from the Hopf fibration:
          - Mode 0 (l=1): driven by cos²(θ/2) = (holonomy/π)²
            where holonomy = πcos(χ) at hyper-latitude χ = θ/2.
            This is the fiber-aligned component.
          - Mode 1 (l=3): driven by sin²(θ/2) = (curvature × 2)²
            where curvature = ½sin(χ).
            This is the fiber-transverse component.

        The cos²/sin² decomposition is the standard SU(2) measurement
        basis — it follows directly from the Hopf bundle structure.
        """
        c2 = float(np.cos(theta / 2.0) ** 2)
        s2 = float(np.sin(theta / 2.0) ** 2)
        return {
            0: c2,
            1: s2,
        }

    def evolve_history(
        self,
        theta_a: float,
        theta_b: float,
        *,
        n_steps: int = 192,
        persist: bool = False,
    ) -> BranchWeightResult:
        """Evolve the detector-conditioned cavity history and return branch weights.

        The phase_spin is now derived from Hopf holonomy (not tuned):

            phase_spin = α_spin × θ_AB + π[cos(θ_a) − cos(θ_b)] / 2

        where the first term is the geodesic spin transport phase
        and the second is the relative detector holonomy from the
        Hopf connection A = ½cos(χ)dφ.

        Parameters
        ----------
        theta_a, theta_b : float
            Detector settings for the two mouths.
        n_steps : int
            Number of history steps over one half-round-trip time.
        persist : bool
            If True, mutate the stored cavity state. By default a copy is used
            so Bell predictions are deterministic across repeated calls.
        """
        from geometrodynamics.bell.hopf_phases import derived_phase_spin

        cavity = self.cavity if persist else deepcopy(self.cavity)
        tau = max(self.tau_semi, 1e-9)
        dt = tau / max(n_steps, 1)
        sigma_t = 0.18 * tau
        t_a = 0.25 * tau
        t_b = 0.75 * tau
        anti_quality = self.pair.antipodal_quality

        # Derived from Hopf holonomy — no tuned constants
        theta_transport = self.pair.geodesic_separation
        phase_spin = derived_phase_spin(theta_a, theta_b, theta_transport)

        resp_a = self._detector_mode_response(theta_a)
        resp_b = self._detector_mode_response(theta_b)

        packet0 = 0.0
        packet_pi = 0.0
        n_packets_0 = 0
        n_packets_pi = 0

        def gauss(t: float, centre: float) -> float:
            x = (t - centre) / max(sigma_t, 1e-12)
            return float(np.exp(-0.5 * x * x))

        for i in range(n_steps):
            t = (i + 0.5) * dt
            env_a = gauss(t, t_a)
            env_b = gauss(t, t_b)
            env_a_adv = gauss(t, tau - t_a)
            env_b_adv = gauss(t, tau - t_b)

            emit_srcs: dict[int, float] = {}
            adv_srcs: dict[int, float] = {}
            for mode_n in cavity.modes:
                # Local detector responses feed the retarded/advanced source terms.
                emit_srcs[mode_n] = CAVITY_ALPHA_EMIT * (
                    resp_a.get(mode_n, 0.0) * env_a +
                    resp_b.get(mode_n, 0.0) * env_b
                )
                adv_srcs[mode_n] = CAVITY_ALPHA_ADV * anti_quality * (
                    resp_b.get(mode_n, 0.0) * env_a_adv +
                    resp_a.get(mode_n, 0.0) * env_b_adv
                )

            cavity.step(emit_srcs, adv_srcs, dt)
            # The cavity's instantaneous phase is the throat phase entering the
            # branch-closure rule.
            _, phi_throat = cavity.confirm_payload()
            pkts = cavity.check_and_fire_packets(
                tau_semi=tau,
                phi_spin=phase_spin,
                phi_throat=phi_throat,
                t_now=t,
            )
            for pkt in pkts:
                if pkt.branch == 0:
                    packet0 += pkt.delta_b
                    n_packets_0 += 1
                else:
                    packet_pi += pkt.delta_b
                    n_packets_pi += 1

        residual0 = 0.0
        residual_pi = 0.0
        _, phi_throat = cavity.confirm_payload()
        sigma_phase = 0.6
        for mode in cavity.modes.values():
            raw = (mode.omega * tau + phase_spin + phi_throat + np.pi) % (2.0 * np.pi) - np.pi
            mm0 = abs(raw)
            mmpi = abs(((raw - np.pi) + np.pi) % (2.0 * np.pi) - np.pi)
            amp = abs(mode.b)
            residual0 += amp * np.exp(-(mm0 ** 2) / (2.0 * sigma_phase ** 2))
            residual_pi += amp * np.exp(-(mmpi ** 2) / (2.0 * sigma_phase ** 2))

        # Small branch prior from the actual mode content avoids numerical collapse
        # to exactly zero while still letting the history dominate.
        energy = cavity.energy()
        prior = 1e-6 * max(energy, 1.0)
        score0 = packet0 + residual0 + prior
        score_pi = packet_pi + residual_pi + prior
        total = score0 + score_pi
        if total <= 1e-30:
            w0 = wpi = 0.5
        else:
            w0 = score0 / total
            wpi = score_pi / total

        if persist:
            self.pair.shared_phase = phi_throat

        return BranchWeightResult(
            w0=float(w0),
            wpi=float(wpi),
            packet_weight_0=float(packet0),
            packet_weight_pi=float(packet_pi),
            residual_weight_0=float(residual0),
            residual_weight_pi=float(residual_pi),
            n_packets_0=n_packets_0,
            n_packets_pi=n_packets_pi,
            energy_final=float(energy),
            phase_spin=phase_spin,
        )

    @property
    def branch_weights(self) -> tuple[float, float]:
        """Compatibility property using the default equal-setting history."""
        res = self.evolve_history(0.0, 0.0)
        return (res.w0, res.wpi)


def make_bell_pair(
    chi_a: float = 0.5,
    theta_a: float = 0.3,
    phi_a: float = 0.0,
    omega_0: float = 1.055,
    omega_1: float = 1.219,
) -> BellPair:
    """Create a Bell pair at the specified S³ location.

    The cavity frequencies are set to the Tangherlini l=1 and l=3
    eigenvalues, matching the physical throat spectrum.
    """
    pair = make_singlet_pair(chi_a, theta_a, phi_a)
    cavity = make_cavity(pair_key=(0, 1), omega_0=omega_0, omega_1=omega_1)
    return BellPair(pair=pair, cavity=cavity)
