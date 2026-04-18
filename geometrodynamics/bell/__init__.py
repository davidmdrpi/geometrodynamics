"""
Bell subpackage — Bell correlations from geometric throat topology.

Supports two equivalent readings of Bell correlations within the
same geometric framework:

  Path A (dynamic):  Time-symmetric transactional realization via
                     cavity history evolution + Hopf-derived phases.
  Path B (kinematic): Higher-dimensional shared-bulk identity —
                      pure geometric projection, no time stepping.

Both paths produce E(a,b) = −cos(a−b), CHSH = 2√2, no-signaling.

Modules:
  pair_state      – BellPair with cavity history evolution (Path A)
  analyzers       – detector settings as SU(2) boundary operators
  correlations    – E(a,b), CHSH, no-signaling (Path A entry point)
  hopf_phases     – Bell closure phases derived from Hopf holonomy
  bulk_identity   – BulkConnectedPair, kinematic Bell (Path B)
"""

from geometrodynamics.bell.pair_state import (
    BellPair,
    BranchWeightResult,
    make_bell_pair,
)

from geometrodynamics.bell.analyzers import (
    throat_transport,
    rotation_matrix,
    measurement_projector,
    singlet_state,
    triplet_m0_state,
    singlet_amplitude,
    triplet_amplitude,
    outcome_probability,
    branch_weighted_probability,
    singlet_correlation,
)

from geometrodynamics.bell.correlations import (
    correlation,
    chsh,
    check_no_signaling,
    correlation_curve,
    CHSHResult,
    NoSignalingResult,
)

from geometrodynamics.bell.hopf_phases import (
    geodesic_spin_phase,
    detector_holonomy_phase,
    derived_phase_spin,
    compare_phase_formulas,
)

from geometrodynamics.bell.bulk_identity import (
    BulkConnectedPair,
    make_bulk_pair,
    boundary_projector,
    bulk_amplitude,
    bulk_probability,
    bulk_correlation,
    bulk_marginal,
    bulk_chsh,
)

__all__ = [
    "BellPair",
    "make_bell_pair",
    "throat_transport",
    "rotation_matrix",
    "measurement_projector",
    "singlet_state",
    "singlet_amplitude",
    "outcome_probability",
    "singlet_correlation",
    "correlation",
    "chsh",
    "check_no_signaling",
    "correlation_curve",
    "CHSHResult",
    "NoSignalingResult",
]
