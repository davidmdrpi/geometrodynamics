"""
Hopf fibration geometry on S³.

Provides the canonical connection 1-form A = ½cos(χ)dφ, its curvature,
gauge holonomy, Chern number integration, and SU(2) spinor transport.
"""

from geometrodynamics.hopf.connection import (
    hopf_connection,
    hopf_curvature,
    hopf_holonomy,
    hopf_circle,
)
from geometrodynamics.hopf.chern import compute_c1
from geometrodynamics.hopf.spinor import compute_spinor_monodromy
