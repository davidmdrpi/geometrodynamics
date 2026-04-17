"""
Embedding subpackage — non-orientable throat topology.

Encodes the extra-dimensional wormhole structure that makes particles
identical, gives them spin-½, and connects particle-antiparticle pairs
through a shared non-orientable throat.
"""

from geometrodynamics.embedding.topology import (
    ThroatDefect,
    ConjugatePair,
    make_singlet_pair,
)

from geometrodynamics.embedding.transport import (
    derive_throat_transport,
    verify_transport_properties,
    orientation_reversal_on_s3,
    verify_hopf_preservation,
    derive_singlet_from_transport,
)

__all__ = [
    "ThroatDefect",
    "ConjugatePair",
    "make_singlet_pair",
    "derive_throat_transport",
    "verify_transport_properties",
    "orientation_reversal_on_s3",
    "verify_hopf_preservation",
    "derive_singlet_from_transport",
]
