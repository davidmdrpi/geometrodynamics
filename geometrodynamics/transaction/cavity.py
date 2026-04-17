"""
Antipodal S³ cavity — persistent resonant state between wormhole pairs.

Ported from v39.  The S³ spatial topology is treated as a resonant cavity.
Each antipodal particle-antiparticle pair has a set of CavityMode
oscillators b_n(t) that accumulate energy from local throat excitation
(retarded, S_emit) and advanced absorber response (S_adv).

A discrete momentum packet fires when the Bohr-like resonance condition

    ω_n τ_semi + φ_spin + φ_throat  ≡  0 or π  (mod 2π)

is satisfied and the mode amplitude exceeds CAVITY_BMIN.

This is the bridge between one-shot transaction logic and globally closed
Bell-capable histories.  The cavity supplies shared history, persistent
memory, and the global closure condition that Bell correlations depend on.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from geometrodynamics.constants import (
    CAVITY_GAMMA,
    CAVITY_BMIN,
    CAVITY_LAMBDA,
    CAVITY_SOFT_COUP,
    CAVITY_ALPHA_EMIT,
    CAVITY_ALPHA_ADV,
    CAVITY_PACKET_FRAC,
    PHASE_MATCH_MAX,
    C_GW,
)
from geometrodynamics.transaction.s3_geometry import geo4, antipode4, s3_tangent_direction
from geometrodynamics.transaction.particles import MouthState


def _wrap_phase(phi: float) -> float:
    """Wrap angle to [−π, π)."""
    return float((phi + np.pi) % (2.0 * np.pi) - np.pi)


@dataclass
class CavityMode:
    """One resonant mode of the S³ antipodal cavity for a specific pair.

    ODE:  b̈_n + 2γ_n ḃ_n + ω_n² b_n  =  S_n^emit(t) + S_n^adv(t)

    S^emit is sourced by local throat mode velocity (retarded, post GW hit).
    S^adv  is sourced by the absorber's throat velocity (causal-gate-gated).

    The mode is a harmonic oscillator with the same eigenfrequency as the
    corresponding Tangherlini throat mode — same quantum numbers, expressed
    as a global cavity standing wave rather than a local wormhole oscillation.
    """

    n: int                     # cavity index: 0 → l=1 throat, 1 → l=3 throat
    omega: float               # eigenfrequency (matched to Tangherlini)
    gamma: float = CAVITY_GAMMA
    b: float = 0.0             # real envelope amplitude
    bdot: float = 0.0
    n_packets: int = 0         # discrete packets fired so far

    def step(self, src_emit: float, src_adv: float, dt: float) -> None:
        """Integrate the cavity mode ODE one timestep."""
        acc = (
            src_emit + src_adv
            - 2.0 * self.gamma * self.bdot
            - self.omega ** 2 * self.b
        )
        self.bdot += dt * acc
        self.b += dt * self.bdot

    def energy(self) -> float:
        """Mode energy ½(ḃ² + ω²b²)."""
        return 0.5 * (self.bdot ** 2 + self.omega ** 2 * self.b ** 2)

    def instantaneous_phase(self) -> float:
        """Phase of the cavity oscillator in the (b, ḃ/ω) plane."""
        return float(
            np.angle(self.b + 1j * self.bdot / max(self.omega, 1e-9))
        )

    def closure_check(
        self,
        tau_semi: float,
        phi_spin: float,
        phi_throat: float,
        branch_tol: float = PHASE_MATCH_MAX,
    ) -> tuple[float, int, bool]:
        """Bohr-like resonance condition for the S³ cavity half-round-trip.

            ω_n τ_semi + φ_spin + φ_throat  ≡  0   (mod 2π)   0-branch
            ω_n τ_semi + φ_spin + φ_throat  ≡  π   (mod 2π)   π-branch

        Returns (mismatch_rad, branch_int, is_closed).

        The π-branch naturally handles the SU(2) holonomy offset for
        antipodal transport: the spinor picks up a half-angle phase of π/2
        at θ=π, making the sum ≈ π not ≈ 0.
        """
        raw = _wrap_phase(self.omega * tau_semi + phi_spin + phi_throat)
        mm0 = abs(raw)
        mmpi = abs(_wrap_phase(raw - np.pi))
        if mm0 <= mmpi:
            return mm0, 0, (mm0 <= branch_tol and abs(self.b) > CAVITY_BMIN)
        return mmpi, 1, (mmpi <= branch_tol and abs(self.b) > CAVITY_BMIN)


@dataclass
class CavityPacket:
    """Immutable record of one discrete momentum packet transfer."""

    t: float
    pair_key: tuple
    mode_n: int
    branch: int        # 0 = direct (0-phase), 1 = π-shifted
    delta_b: float     # cavity amplitude consumed
    delta_p: float     # momentum kick applied to each end
    mismatch: float    # resonance mismatch at fire time
    phi_spin: float    # SU(2) spin holonomy at fire time


@dataclass
class AntipodalCavity:
    """Persistent resonant S³ cavity between one antipodal pair.

    Contains one CavityMode per Tangherlini quantum number.  Driven by
    local throat excitation (S^emit) and the advanced absorber response
    (S^adv, gated by the causal constraint).  Transfers discrete momentum
    packets when the resonance condition closes; logs each transfer.

    Between packet transfers the cavity provides a continuous soft force
    via CAVITY_SOFT_COUP — the mean-field analogue of Coulomb attraction
    while the quantum exchange accumulates.
    """

    pair_key: tuple
    modes: Dict[int, CavityMode] = field(default_factory=dict)
    packet_log: List[CavityPacket] = field(default_factory=list)

    def energy(self) -> float:
        """Total energy across all cavity modes."""
        return sum(m.energy() for m in self.modes.values())

    def step(
        self,
        emit_srcs: Dict[int, float],
        adv_srcs: Dict[int, float],
        dt: float,
    ) -> None:
        """Integrate all cavity mode ODEs for one timestep."""
        for n, mode in self.modes.items():
            mode.step(emit_srcs.get(n, 0.0), adv_srcs.get(n, 0.0), dt)

    def complex_amplitude(self) -> complex:
        """Combined phasor of all modes."""
        z = 0.0 + 0.0j
        for mode in self.modes.values():
            z += mode.b + 1j * mode.bdot / max(mode.omega, 1e-9)
        return z

    def confirm_payload(self) -> tuple[float, float]:
        """Cavity-mediated confirmation factor for transactions.

        Returns (envelope, phase) where the envelope rises with cavity
        amplitude and the phase is the instantaneous cavity phase.
        """
        z = self.complex_amplitude()
        mag = abs(z)
        if mag <= 1e-12:
            return 0.0, 0.0
        scale = 0.25 * CAVITY_BMIN
        env = float(mag / (mag + scale))
        return env, float(np.angle(z))

    def check_and_fire_packets(
        self,
        tau_semi: float,
        phi_spin: float,
        phi_throat: float,
        t_now: float,
    ) -> List[CavityPacket]:
        """Check all modes for resonance closure; fire packets where closed.

        Returns list of fired packets (may be empty).
        """
        fired = []
        for n, mode in self.modes.items():
            mismatch, branch, closed = mode.closure_check(
                tau_semi, phi_spin, phi_throat,
            )
            if not closed:
                continue

            branch_sign = 1.0 if branch == 0 else -1.0
            delta_b = CAVITY_PACKET_FRAC * abs(mode.b)
            delta_p = CAVITY_LAMBDA * delta_b * branch_sign

            # Reduce cavity amplitude
            mode.b *= (1.0 - CAVITY_PACKET_FRAC)
            mode.bdot *= (1.0 - CAVITY_PACKET_FRAC)
            mode.n_packets += 1

            pkt = CavityPacket(
                t=t_now,
                pair_key=self.pair_key,
                mode_n=n,
                branch=branch,
                delta_b=delta_b,
                delta_p=delta_p,
                mismatch=mismatch,
                phi_spin=phi_spin,
            )
            self.packet_log.append(pkt)
            fired.append(pkt)
        return fired


# ── Cavity factory ───────────────────────────────────────────────────────────

def make_cavity(
    pair_key: tuple,
    omega_0: float = 1.055,
    omega_1: float = 1.219,
) -> AntipodalCavity:
    """Create a two-mode antipodal cavity for a particle-antiparticle pair.

    Default frequencies are the l=1 and l=3 Tangherlini mode eigenvalues.
    """
    return AntipodalCavity(
        pair_key=pair_key,
        modes={
            0: CavityMode(n=0, omega=omega_0),
            1: CavityMode(n=1, omega=omega_1),
        },
    )
