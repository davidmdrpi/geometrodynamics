"""
Black-hole subpackage — coherent wormhole-throat condensate model.

A black hole is modelled as a macroscopic coherent condensate of the
same non-orientable wormhole throats that appear microscopically as
particle-antiparticle pairs.  This subpackage provides:

  condensate   – CoherentCondensate, ThroatState, constructors
  interior     – Hayward regular metric, geodesics, horizon finder
  entropy      – Bekenstein-Hawking from throat counting, first law
  derivation   – Bridging condensate → metric via Einstein equations
"""

from geometrodynamics.blackhole.condensate import (
    CoherentCondensate,
    ThroatState,
    build_schwarzschild_condensate,
    build_charged_condensate,
    L_PLANCK_SQ,
    A_THROAT_MIN,
)

from geometrodynamics.blackhole.interior import (
    f_schwarzschild,
    f_hayward,
    df_hayward_dr,
    find_horizons,
    critical_core_scale,
    surface_gravity,
    hawking_temperature,
    kretschner_hayward,
    integrate_radial_geodesic,
    tortoise_hayward,
    RadialGeodesic,
)

from geometrodynamics.blackhole.entropy import (
    compute_entropy_balance,
    check_first_law,
    information_capacity_bits,
    page_time,
    evaporation_step,
    EntropyBalance,
    FirstLawCheck,
)

from geometrodynamics.blackhole.derivation import (
    throat_density,
    verify_throat_normalisation,
    stress_energy_from_throat_density,
    mass_function_from_density,
    mass_function_analytic,
    derive_metric_from_condensate,
    derive_core_scale,
    derive_temperature_from_modes,
    check_energy_conditions,
    full_derivation,
    StressEnergy,
    MetricDerivation,
    TemperatureDerivation,
    EnergyConditions,
    CoreScaleConstraints,
    FullDerivation,
)

__all__ = [
    # condensate
    "CoherentCondensate",
    "ThroatState",
    "build_schwarzschild_condensate",
    "build_charged_condensate",
    "L_PLANCK_SQ",
    "A_THROAT_MIN",
    # interior
    "f_schwarzschild",
    "f_hayward",
    "df_hayward_dr",
    "find_horizons",
    "critical_core_scale",
    "surface_gravity",
    "hawking_temperature",
    "kretschner_hayward",
    "integrate_radial_geodesic",
    "tortoise_hayward",
    "RadialGeodesic",
    # entropy
    "compute_entropy_balance",
    "check_first_law",
    "information_capacity_bits",
    "page_time",
    "evaporation_step",
    "EntropyBalance",
    "FirstLawCheck",
]
