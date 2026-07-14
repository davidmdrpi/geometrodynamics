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
from geometrodynamics.transaction.network import (
    NetworkMouth,
    NetworkThroat,
    NetworkTraversal,
    TraversalLeg,
    closure_offset,
    emergent_frequency,
    network_confirmation,
    projected_kernel,
    s3_leg,
    traverse_throat,
)
from geometrodynamics.transaction.network import (
    MouthPort,
    emergence_train,
    transparent_port,
)
from geometrodynamics.transaction.network import (
    effective_green,
    loop_eigenvalue,
)
