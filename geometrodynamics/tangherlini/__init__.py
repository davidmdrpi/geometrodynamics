"""
5D Tangherlini wormhole eigenmodes and Maxwell solver.

Provides Chebyshev spectral solvers for the radial eigenvalue problem
on the 5D Schwarzschild-Tangherlini background, the sourced Maxwell BVP
that validates the Coulomb law from eigenfunction throat flux, and the
derived α_q coupling ratios.
"""

from geometrodynamics.tangherlini.radial import (
    solve_radial_modes,
    V_tangherlini,
    r_to_rstar,
    rstar_to_r,
)
from geometrodynamics.tangherlini.maxwell import solve_maxwell_from_eigenmode
from geometrodynamics.tangherlini.alpha_q import derive_alpha_q, throat_du_dr
from geometrodynamics.tangherlini.lepton_spectrum import (
    Crossing,
    LadderFit,
    LEPTON_BASELINE_DEPTHS,
    LEPTON_BASELINE_PHASE,
    LEPTON_BASELINE_PINHOLE,
    LEPTON_BASELINE_RESISTANCE,
    LEPTON_BASELINE_TRANSPORT,
    S3_ACTION_BASE,
    TAU_BETA_50PI,
    calibrate_electron_predict_heavier,
    compute_knotted_lepton_spectrum,
    compute_tunneling_envelope,
    derive_geometric_beta,
    solved_lepton_masses_mev,
    tau_uplift_2pi_quanta,
    tune_transport_and_resistance,
)
