"""
Derived source law: α_q(l, n) from eigenfunction throat flux.

The electromagnetic 4-current d*F = J is localised at the wormhole throat
shell r = R_MID.  The coupling for every mode is the dimensionless ratio

    α_q(l, n) = (du_{l,n}/dr)|_throat / |du_{1,0}/dr|_throat

so α_q(1, 0) = ±1 by construction.
"""

import numpy as np

from geometrodynamics.constants import R_MID

_N_THROAT_FIT = 8  # near-throat points for forced-origin fit


def throat_du_dr(fn: dict, r_mid: float = R_MID) -> float:
    """Throat radial derivative du/dr|_{r=R_MID} via forced-origin slope.

    Because u(R_MID) = 0 exactly (Dirichlet BC), u ≈ A·(r − R_MID) near
    the throat.  The optimal A is

        A = Σ u_k · Δr_k / Σ Δr_k²,  Δr_k = r_phys[k] − R_MID

    using the first _N_THROAT_FIT non-zero grid points.
    """
    r = fn["r_phys"]
    u = fn["u_half"]
    dr = r[1:_N_THROAT_FIT] - r_mid
    return float(np.dot(dr, u[1:_N_THROAT_FIT]) / np.dot(dr, dr))


def derive_alpha_q(modes: dict[int, dict], r_mid: float = R_MID) -> dict:
    """Build the α_q table for all computed (l, n=0) ground-state modes.

    Parameters
    ----------
    modes : dict mapping l → {'omega': ..., 'funcs': [...], ...}
        As returned by computing solve_radial_modes for each l.

    Returns
    -------
    dict mapping (l, 0) → float
        Signed ratio A_{l,0} / |A_{1,0}|.
    """
    A_ref = abs(throat_du_dr(modes[1]["funcs"][0], r_mid))
    table = {}
    for l, data in modes.items():
        fn = data["funcs"][0]  # n=0 ground state
        A = throat_du_dr(fn, r_mid)
        table[(l, 0)] = float(A / A_ref)
    return table
