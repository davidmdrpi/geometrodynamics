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
    DEFAULT_ASSIGNMENT_ANGULAR,
    DEFAULT_ASSIGNMENT_RADIAL,
    HBAR_C_MEV_FM,
    PDG_LEPTON_MASSES_MEV,
    LeptonAssignment,
    LeptonMode,
    LeptonSpectrumReport,
    compute_lepton_spectrum,
    format_report,
)
from geometrodynamics.tangherlini.multipass import (
    AMPLIFIERS,
    KNOT_DEPTHS,
    MultiPassFit,
    amp_exponential,
    amp_linear,
    amp_power,
    evaluate_multipass,
    format_multipass_fit,
    format_multipass_scan,
    scan_multipass,
)
