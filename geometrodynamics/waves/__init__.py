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
from geometrodynamics.waves.two_source import (  # noqa: F401
    defect_of_pair,
    disconnection_defect,
    energy_functional,
    interaction_energy,
    invisible_partner,
    is_real_field_compatible,
    mouth_channel_invariant,
    recover_boundary,
    recover_complex_response,
    recover_response,
    response_matrix,
    static_response,
)
from geometrodynamics.waves.two_wave import (  # noqa: F401
    GaussianPulse,
    RetardedGrid,
    TwoWaveSetup,
    arrival_directions,
    contract_stress,
    gamma_omega,
    green_omega,
    normalized_invariant,
    solve_field,
    stress_tensor,
    wkb_invariant,
)
from geometrodynamics.waves.finite_throat import (  # noqa: F401
    FiniteThroat,
    bounce_delays,
    causal_onset,
    dtn_matrix,
    green_identity_residual,
    impulse_response,
    interior_profile,
    response_spectrum,
)
from geometrodynamics.waves.finite_mouth import (  # noqa: F401
    FiniteMouthThroat,
    mouth_green,
    regular_radial,
    screened_products,
    shell_average_cross,
    shell_average_self,
)
from geometrodynamics.waves.neck import (  # noqa: F401
    WORKING_NECK,
    NeckThroat,
    exterior_dtn,
    exterior_dtn_monopole,
    exterior_log_derivative,
    rayleigh_quotient,
)
