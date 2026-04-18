"""
Bell phase derivation from Hopf holonomy.

Derives the phase_spin entering the cavity closure condition from
the existing Hopf connection and SU(2) transport infrastructure,
replacing the tuned 0.125 factor with a geometric formula.

The derived phase has two contributions:

1. Geodesic transport phase:  α_spin × θ_AB
   For spin-½ (α_spin = 0.5) and an antipodal pair (θ_AB = π):
   baseline = π/2.

2. Detector holonomy:  Δφ_det = [holonomy(θ_a) − holonomy(θ_b)] / 2
   From the Hopf connection A = ½cos(χ)dφ, the holonomy at
   hyper-latitude χ is πcos(χ).  The detector setting θ maps to
   χ = θ (the measurement axis selects a Hopf fiber), giving:
   Δφ_det = π[cos(θ_a) − cos(θ_b)] / 2

The total:
   phase_spin = π/2 + π[cos(θ_a) − cos(θ_b)] / 2

This produces setting-dependent branch weights that are derived
from the geometry rather than tuned.

Note: The derived formula may produce different CHSH values than
the tuned version.  Both are tested and compared.
"""

from __future__ import annotations

import numpy as np

from geometrodynamics.hopf.connection import hopf_holonomy
from geometrodynamics.transaction.s3_geometry import geo4


def geodesic_spin_phase(theta_transport: float, alpha_spin: float = 0.5) -> float:
    """SU(2) spin phase from geodesic parallel transport.

    For a geodesic of angle θ on S³, a spin-½ spinor accumulates
    a geometric phase of α_spin × θ.

    For α_spin = 0.5 (spin-½) and θ = π (antipodal):
        phase = π/2
    """
    return alpha_spin * theta_transport


def detector_holonomy_phase(theta_a: float, theta_b: float) -> float:
    """Relative Hopf holonomy between two detector settings.

    Each detector at angle θ projects onto the Hopf fiber at
    hyper-latitude χ = θ.  The holonomy at that fiber is πcos(θ).

    The relative phase between detectors A and B:
        Δφ = [πcos(θ_a) − πcos(θ_b)] / 2

    The factor of 1/2 comes from the spin-½ representation:
    the detector phase enters through the SU(2) representation
    of the Hopf fiber rotation.
    """
    hol_a = float(hopf_holonomy(theta_a))  # πcos(θ_a)
    hol_b = float(hopf_holonomy(theta_b))  # πcos(θ_b)
    return (hol_a - hol_b) / 2.0


def derived_phase_spin(
    theta_a: float,
    theta_b: float,
    theta_transport: float = np.pi,
    alpha_spin: float = 0.5,
) -> float:
    """Full Hopf-derived phase_spin for the Bell closure condition.

    phase_spin = geodesic_spin_phase + detector_holonomy_phase
               = α_spin × θ_AB + π[cos(θ_a) − cos(θ_b)] / 2

    For an antipodal pair (θ_AB = π) with equal settings (θ_a = θ_b):
        phase_spin = π/2

    Parameters
    ----------
    theta_a, theta_b : float
        Detector settings in radians.
    theta_transport : float
        Geodesic distance between the pair mouths on S³.
        Default π (antipodal).
    alpha_spin : float
        Spin coupling constant.  0.5 for spin-½.
    """
    geo_phase = geodesic_spin_phase(theta_transport, alpha_spin)
    det_phase = detector_holonomy_phase(theta_a, theta_b)
    raw = geo_phase + det_phase
    # Wrap to [−π, π)
    return float(((raw + np.pi) % (2.0 * np.pi)) - np.pi)


def compare_phase_formulas(
    theta_a: float,
    theta_b: float,
) -> dict:
    """Compare derived vs tuned phase_spin values.

    Returns both values and their difference for diagnostic purposes.
    """
    derived = derived_phase_spin(theta_a, theta_b)
    tuned = float(
        ((0.125 * (theta_a - theta_b) + np.pi) % (2.0 * np.pi)) - np.pi
    )
    return {
        "theta_a": theta_a,
        "theta_b": theta_b,
        "derived": derived,
        "tuned": tuned,
        "difference": derived - tuned,
    }
