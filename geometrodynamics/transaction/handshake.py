"""
Wheeler–Feynman absorber-theory handshake on S³.

Retarded offer → advanced confirm → phase-closure transaction.
Implements the full amplitude algebra: hit envelope, antipodal match,
mode overlap, SU(2) spin phase, and time-symmetric phase closure.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

from geometrodynamics.constants import (
    EPS_HIT,
    SIGMA_ANTI,
    OFFER_TTL,
    CONFIRM_TTL,
    PHASE_MATCH_SIGMA,
    R_MID,
    S3_GREEN_EPS,
)
from geometrodynamics.transaction.s3_geometry import geo4, antipode4
from geometrodynamics.transaction.particles import MouthState, Particle4


# ── Helper functions ─────────────────────────────────────────────────────────

def gw_hit_envelope(psi_diff: float, eps_hit: float = EPS_HIT) -> float:
    """Gaussian hit envelope for GW shell / particle overlap."""
    return float(np.exp(-(psi_diff ** 2) / (2.0 * eps_hit ** 2)))


def antipodal_match_weight(c_p4, d_p4, sigma: float = SIGMA_ANTI) -> float:
    """How close candidate and destination are to being antipodal."""
    anti_err = geo4(c_p4, antipode4(d_p4))
    return float(np.exp(-(anti_err ** 2) / (2.0 * sigma ** 2)))


def complex_mode_overlap(mouth_c: MouthState, mouth_d: MouthState) -> complex:
    """Normalised complex inner product of two mouths' mode spectra."""
    keys = set(mouth_c.modes).intersection(mouth_d.modes)
    if not keys:
        return 0.0 + 0.0j
    num = 0.0 + 0.0j
    den_c = den_d = 0.0
    for k in keys:
        ac = mouth_c.modes[k].complex_amplitude()
        ad = mouth_d.modes[k].complex_amplitude()
        num += np.conj(ac) * ad
        den_c += abs(ac) ** 2
        den_d += abs(ad) ** 2
    if den_c <= 1e-18 or den_d <= 1e-18:
        return 0.0 + 0.0j
    return num / np.sqrt(den_c * den_d)


def mode_overlap(mouth_c: MouthState, mouth_d: MouthState) -> float:
    return float(abs(complex_mode_overlap(mouth_c, mouth_d)))


def su2_spin_phase(
    theta_transport: float,
    mouth_c: Optional[MouthState] = None,
    mouth_d: Optional[MouthState] = None,
    alpha_spin: float = 0.5,
) -> complex:
    """SU(2) spin phase from geodesic transport + relative mouth phase."""
    rel_phase = 0.0
    if mouth_c is not None and mouth_d is not None:
        coh = complex_mode_overlap(mouth_c, mouth_d)
        if abs(coh) > 1e-18:
            rel_phase = float(np.angle(coh))
        else:
            rel_phase = mouth_d.mean_phase() - mouth_c.mean_phase()
    return complex(np.exp(1j * (alpha_spin * theta_transport + 0.5 * rel_phase)))


def wrap_phase(phi: float) -> float:
    return float((phi + np.pi) % (2.0 * np.pi) - np.pi)


def mouth_activity(mouth: MouthState) -> float:
    """RMS amplitude of all active throat modes."""
    if not mouth.modes:
        return 0.0
    total = 0.0
    for m in mouth.modes.values():
        scale = max(m.omega, 1e-9)
        total += m.a ** 2 + (m.adot / scale) ** 2
    return float(np.sqrt(max(total, 0.0)))


# ── Signal dataclasses ───────────────────────────────────────────────────────

@dataclass
class OfferSignal:
    """Retarded offer from source, via candidate, toward absorber.

    Stores the source particle's geometric charge at offer-creation time
    so that ``complete_transaction`` can populate ``Transaction.q_src``
    from the true emitter, not the relay candidate.
    """

    grav_id: int
    src_pid: int
    cand_pid: int
    t_birth: float
    hit_weight: float
    theta_src_cand: float
    response: float
    amp: complex
    q_src: float = 0.0        # source charge, set by make_offer
    geom_phase: float = 0.0
    ttl: float = OFFER_TTL

    def is_alive(self, t_now: float) -> bool:
        return (t_now - self.t_birth) <= self.ttl and abs(self.amp) > 1e-12


@dataclass
class Transaction:
    """Confirmed handshake: offer + confirm amplitudes pass phase closure."""

    grav_id: int
    src_pid: int
    cand_pid: int
    dst_pid: int
    t_birth: float
    amp: complex
    offer_amp: complex
    confirm_amp: complex
    q_src: float
    q_dst: float
    field_at_dst: float
    overlap: float
    anti_weight: float
    spin_phase: float
    hit_weight: float
    pair_kernel: float
    phase_mismatch: float
    phase_match_weight: float
    is_confirmed: bool = True
    ttl: float = CONFIRM_TTL

    def is_alive(self, t_now: float) -> bool:
        return (
            self.is_confirmed
            and (t_now - self.t_birth) <= self.ttl
            and abs(self.amp) > 1e-12
        )


# ── Amplitude constructors ───────────────────────────────────────────────────

def retarded_offer_amplitude(
    src: Particle4,
    cand: Particle4,
    hit_weight: float,
    gw_radius: float,
) -> tuple[complex, float, float, float]:
    """Compute retarded offer amplitude from source through candidate."""
    theta_src_cand = geo4(src.p4, cand.p4)
    response = mouth_activity(cand.mouth)
    response = max(response, 0.10 * hit_weight)
    geom_phase = 0.5 * (gw_radius + theta_src_cand)
    phase = geom_phase + cand.mouth.mean_phase()
    amp = hit_weight * response * np.exp(1j * phase)
    return complex(amp), theta_src_cand, float(response), float(geom_phase)


def advanced_confirm_amplitude(
    cand: Particle4,
    dst: Particle4,
    dst_field_sol: Optional[dict],
    sample_radius: float = R_MID + S3_GREEN_EPS,
) -> tuple[complex, float, float, float, float, float]:
    """Advanced (time-reversed) confirmation from absorber dst.

    Parameters
    ----------
    dst_field_sol : dict or None
        Radial field solution, as returned by
        ``tangherlini.maxwell.solve_maxwell_from_eigenmode()``.
        Must contain keys ``"r"`` (radial grid) and ``"E_r"``
        (electric field).  If None, falls back to mouth_activity.
    sample_radius : float
        Radius at which to sample the absorber's self-field.
        Default is R_MID + S3_GREEN_EPS (just outside the throat,
        where the regularised Green function is well-defined).
    """
    anti_w = antipodal_match_weight(cand.p4, dst.p4)
    coh = complex_mode_overlap(cand.mouth, dst.mouth)
    overlap = float(abs(coh))
    theta_cand_dst = geo4(cand.p4, dst.p4)
    spin_phase = su2_spin_phase(theta_cand_dst, cand.mouth, dst.mouth)

    if dst_field_sol is not None:
        field_at_dst = float(
            np.interp(sample_radius, dst_field_sol["r"], dst_field_sol["E_r"])
        )
    else:
        field_at_dst = float(mouth_activity(dst.mouth)) + 1e-6

    adv_phase = np.conj(spin_phase) * np.exp(
        -1j * (0.5 * theta_cand_dst + dst.mouth.mean_phase())
    )
    amp = anti_w * overlap * field_at_dst * adv_phase
    return (
        complex(amp),
        anti_w,
        overlap,
        field_at_dst,
        float(np.angle(spin_phase)),
        theta_cand_dst,
    )


def retro_phase_match(
    offer_amp: complex,
    confirm_amp: complex,
    sigma: float = PHASE_MATCH_SIGMA,
) -> tuple[float, float]:
    """Time-symmetric phase closure scoring both 0 and π branches."""
    total = wrap_phase(np.angle(offer_amp) + np.angle(confirm_amp))
    mismatch_0 = abs(total)
    mismatch_pi = abs(wrap_phase(total - np.pi))
    mismatch = min(mismatch_0, mismatch_pi)
    weight = float(np.exp(-(mismatch ** 2) / (2.0 * sigma ** 2)))
    return mismatch, weight


# ── Transaction lifecycle helpers ────────────────────────────────────────────

def make_offer(
    src: Particle4,
    cand: Particle4,
    *,
    grav_id: int,
    hit_weight: float,
    gw_radius: float,
    t_now: float,
) -> OfferSignal:
    """Create an OfferSignal from a retarded offer amplitude computation.

    Captures the source particle's geometric charge at creation time so
    that downstream ``complete_transaction`` records the true emitter's
    charge, not the relay candidate's.

    Parameters
    ----------
    src, cand : Particle4
        Source and candidate particles.
    grav_id : int
        ID of the gravitational wave that triggered this offer.
    hit_weight : float
        GW shell / candidate overlap weight.
    gw_radius : float
        Current angular radius of the expanding GW shell.
    t_now : float
        Simulation time at which the offer is created.

    Returns
    -------
    OfferSignal
    """
    amp, theta_sc, response, geom_phase = retarded_offer_amplitude(
        src, cand, hit_weight, gw_radius,
    )
    return OfferSignal(
        grav_id=grav_id,
        src_pid=src.pid,
        cand_pid=cand.pid,
        t_birth=t_now,
        hit_weight=hit_weight,
        theta_src_cand=theta_sc,
        response=response,
        amp=amp,
        q_src=src.mouth.q_geom(),
        geom_phase=geom_phase,
    )


def complete_transaction(
    offer: OfferSignal,
    cand: Particle4,
    dst: Particle4,
    dst_field_sol: Optional[dict],
    *,
    t_now: float,
    min_phase_weight: float = 0.01,
) -> Optional[Transaction]:
    """Attempt to close a transaction from a live offer and an absorber.

    Runs the advanced confirm amplitude, checks phase closure, and if the
    closure weight exceeds ``min_phase_weight``, instantiates and returns
    a confirmed ``Transaction``.  Returns ``None`` if closure fails.

    This is the missing orchestration step between the raw amplitude
    functions and a completed handshake.  A full simulation engine would
    call this for every (offer, candidate absorber) pair.

    Parameters
    ----------
    offer : OfferSignal
        A live (not expired) retarded offer.
    cand : Particle4
        The candidate particle that relayed the offer.
    dst : Particle4
        The proposed absorber particle.
    dst_field_sol : dict or None
        Radial field solution from ``solve_maxwell_from_eigenmode``,
        or None to fall back to mouth activity.
    t_now : float
        Current simulation time.
    min_phase_weight : float
        Minimum phase-closure weight to accept the transaction.

    Returns
    -------
    Transaction or None
    """
    if not offer.is_alive(t_now):
        return None

    confirm_amp, anti_w, overlap, field_at_dst, spin_ang, theta_cd = (
        advanced_confirm_amplitude(cand, dst, dst_field_sol)
    )

    # A zero-amplitude confirm has no physical content (e.g. anti_w ≈ 0
    # for non-antipodal pairs).  Reject before phase-closure check, which
    # would produce a meaningless phase angle from a zero complex number.
    if abs(confirm_amp) < 1e-15:
        return None

    mismatch, weight = retro_phase_match(offer.amp, confirm_amp)

    if weight < min_phase_weight:
        return None

    # Combined amplitude: product of offer and confirm, weighted by closure
    total_amp = offer.amp * confirm_amp * weight

    # Pair kernel: product of hit weight and antipodal match
    pair_kernel = offer.hit_weight * anti_w

    return Transaction(
        grav_id=offer.grav_id,
        src_pid=offer.src_pid,
        cand_pid=offer.cand_pid,
        dst_pid=dst.pid,
        t_birth=t_now,
        amp=total_amp,
        offer_amp=offer.amp,
        confirm_amp=confirm_amp,
        q_src=offer.q_src,
        q_dst=dst.mouth.q_geom(),
        field_at_dst=field_at_dst,
        overlap=overlap,
        anti_weight=anti_w,
        spin_phase=spin_ang,
        hit_weight=offer.hit_weight,
        pair_kernel=pair_kernel,
        phase_mismatch=mismatch,
        phase_match_weight=weight,
        is_confirmed=True,
    )
