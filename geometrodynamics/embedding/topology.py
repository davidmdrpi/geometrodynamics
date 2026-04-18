"""
Embedding topology — non-orientable throat structure.

Encodes the topological defect class that underlies both particles and
Bell pairs.  Every particle is a non-orientable wormhole throat in the
5D embedding.  The key structures:

  ThroatDefect   — one throat with orientation, wrap parity, conjugacy
  ConjugatePair  — a particle-antiparticle pair sharing one throat
  transport operators — one-pass (orientation flip) and two-pass (identity)

The non-orientable wrap means:
  • one pass through the throat flips orientation (particle ↔ antiparticle)
  • two passes return to the original sector (identity)
  • the pair is ONE connected topological object, not two separate particles

For Bell: the shared topology means measurement on one mouth constrains
the global closure condition for the whole pair, producing correlations
without hidden variables or nonlocal signaling.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Tuple

import numpy as np

from geometrodynamics.transaction.s3_geometry import (
    nrm4, antipode4, geo4, hsp,
)


@dataclass
class ThroatDefect:
    """A single non-orientable wormhole throat in the embedding.

    Attributes
    ----------
    orientation : int
        +1 or −1.  The Möbius-like non-orientable topology assigns each
        mouth an orientation sign.  The two mouths of one throat carry
        opposite signs.
    wrap_parity : int
        +1 for zero-wrap sector (even passes), −1 for one-wrap sector
        (odd passes).  Determines the sign of the spinor transport.
    species_id : str
        Defect class identifier.  All electrons share the same species_id
        because the wrap distance is universal — "all electrons are
        identical" is a topological statement, not a coincidence.
    """

    orientation: int = +1
    wrap_parity: int = +1
    species_id: str = "electron"

    @property
    def is_particle(self) -> bool:
        """Convention: orientation +1 = particle."""
        return self.orientation == +1

    @property
    def conjugate_orientation(self) -> int:
        """Orientation of the conjugate (other mouth)."""
        return -self.orientation

    def one_pass_transport(self) -> ThroatDefect:
        """Transport through the throat once: flips orientation + wrap parity.

        This is the non-orientable map.  Going through the throat once
        turns a particle into an antiparticle (orientation flip) and
        changes the wrap parity (one more pass through the Möbius strip).
        """
        return ThroatDefect(
            orientation=-self.orientation,
            wrap_parity=-self.wrap_parity,
            species_id=self.species_id,
        )

    def two_pass_transport(self) -> ThroatDefect:
        """Transport through the throat twice: identity.

        Two passes through a Möbius strip return to the original sector.
        This is the topological origin of the 4π periodicity of spin-½:
        one geometric rotation (2π) = one throat pass = sign flip;
        two rotations (4π) = two passes = identity.
        """
        return ThroatDefect(
            orientation=self.orientation,
            wrap_parity=self.wrap_parity,
            species_id=self.species_id,
        )

    def spinor_sign(self) -> int:
        """Sign acquired by a spinor after transport through this throat.

        One-pass: −1 (the SU(2) sign flip from Hopf holonomy at χ=0).
        Two-pass: +1 (identity, 4π periodicity).
        """
        return self.wrap_parity


@dataclass
class ConjugatePair:
    """A particle-antiparticle pair = one connected non-orientable throat.

    The two mouths (A and B) sit at antipodal points on S³ and share
    the same throat topology.  They are NOT two independent particles
    with hidden variables — they are two views of one topological object.

    For Bell: the shared cavity memory of this pair determines which
    global closure branches are admissible under given detector settings.

    Attributes
    ----------
    mouth_a, mouth_b : ThroatDefect
        The two mouths.  Always have opposite orientation.
    p4_a, p4_b : ndarray
        Positions on S³ (unit 4-vectors).  Should be approximately antipodal.
    shared_phase : float
        Global phase of the pair state (from cavity memory).
    """

    mouth_a: ThroatDefect
    mouth_b: ThroatDefect
    p4_a: np.ndarray = field(default_factory=lambda: np.array([0, 0, 0, 1.0]))
    p4_b: np.ndarray = field(default_factory=lambda: np.array([0, 0, 0, -1.0]))
    shared_phase: float = 0.0

    def __post_init__(self):
        assert self.mouth_a.orientation == -self.mouth_b.orientation, (
            "Conjugate pair mouths must have opposite orientation"
        )
        assert self.mouth_a.species_id == self.mouth_b.species_id, (
            "Conjugate pair mouths must be the same species"
        )

    @property
    def antipodal_quality(self) -> float:
        """How close to perfect antipodal alignment (1.0 = exact)."""
        return float(
            np.exp(-(geo4(self.p4_a, antipode4(self.p4_b))) ** 2 / 0.02)
        )

    @property
    def geodesic_separation(self) -> float:
        """Geodesic angle between the two mouths on S³."""
        return geo4(self.p4_a, self.p4_b)

    @property
    def tau_semi(self) -> float:
        """Half-round-trip time for a signal on S³ between the mouths.

        τ_semi = θ_AB / C_GW where C_GW is the GW propagation speed.
        """
        from geometrodynamics.constants import C_GW
        return self.geodesic_separation / C_GW


def make_singlet_pair(
    chi_a: float = 0.5,
    theta_a: float = 0.3,
    phi_a: float = 0.0,
    species: str = "electron",
) -> ConjugatePair:
    """Create an antipodal singlet-like pair at specified S³ coordinates.

    Mouth A is placed at (χ, θ, φ); mouth B at the exact antipode.
    Orientations are +1/−1 (particle/antiparticle).
    """
    p4_a = hsp(chi_a, theta_a, phi_a)
    p4_b = antipode4(p4_a)

    return ConjugatePair(
        mouth_a=ThroatDefect(orientation=+1, wrap_parity=+1, species_id=species),
        mouth_b=ThroatDefect(orientation=-1, wrap_parity=-1, species_id=species),
        p4_a=p4_a,
        p4_b=p4_b,
        shared_phase=0.0,
    )
