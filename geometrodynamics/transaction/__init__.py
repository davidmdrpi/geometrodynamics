"""
Retrocausal Wheeler–Feynman transaction protocol on S³.

Implements the absorber-theory handshake: retarded offer → advanced
confirm → phase-closure transaction, with S³ Green function propagation
and SU(2) spin-phase matching.
"""

from geometrodynamics.transaction.particles import (
    ThroatMode,
    MouthState,
    Particle4,
    GravWave,
)
from geometrodynamics.transaction.s3_geometry import (
    nrm4,
    antipode4,
    geo4,
    hsp,
    s3_green_potential,
    s3_green_field_kernel,
)
from geometrodynamics.transaction.handshake import (
    OfferSignal,
    Transaction,
    retarded_offer_amplitude,
    advanced_confirm_amplitude,
    retro_phase_match,
    make_offer,
    complete_transaction,
)
