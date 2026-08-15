"""Wave-collision kinematics: where pair creation can and cannot happen."""

from geometrodynamics.waves.pair_creation import (  # noqa: F401
    WavePair,
    breit_wheeler_cross_section,
    crossing_locus,
    crossing_window,
    invariant_along_the_crossing,
    mandelstam_s,
    opening_angle,
    outgoing_momentum,
    threshold_windows,
)
from geometrodynamics.waves.pair_history import (  # noqa: F401
    PairHistorySystem,
    Throat,
    closure_residual,
    feasible_delay_band,
    geodesic_distance,
)
from geometrodynamics.waves.field_solve import (  # noqa: F401
    branch_arrivals,
    esu_frequency,
    field_peaks,
    image_field,
    spectral_field,
    through_throat_arrivals,
)
from geometrodynamics.waves.branch_coupling import (  # noqa: F401
    CoupledThroat,
    branch_pair_matrix,
    coupled_arrivals,
    coupled_propagator,
    coupled_waveform,
    series_radius,
    stability_threshold,
    resonance_poles,
    branch_labels,
    free_branch_propagator,
    leg_branches,
    mouth_transfer,
)
from geometrodynamics.waves.throat_operator import (  # noqa: F401
    BoundaryCondition,
    DirectionalThroat,
    MouthPair,
    free_green,
    is_stable,
    mouth_active_spectrum,
    regularized_green,
    spectrum_by_channel,
    stability_thresholds,
)
from geometrodynamics.waves.throat_positivity import (  # noqa: F401
    apex,
    cone_coordinates,
    inertia_below,
    is_non_negative,
    threshold_matrix,
)
