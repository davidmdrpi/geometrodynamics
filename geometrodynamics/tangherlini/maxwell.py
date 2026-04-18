"""
Sourced Maxwell solver on the wormhole background.

Extracts charge Q from the l=1 eigenmode throat derivative, then validates
d*F = J by solving the Coulomb BVP and comparing to Q/r.

Steps:
  1. Q = R_MID² |du_{1,0}/dr|_{r=R_MID}  (no free parameter)
  2. Solve (1/r²)d/dr[r² dA/dr] = 0 with Neumann BC and Dirichlet A=0
  3. Compare A_solved to exact A = Q/r − Q/R_OUTER
"""

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

from geometrodynamics.constants import R_MID, R_OUTER
from geometrodynamics.tangherlini.alpha_q import throat_du_dr


def solve_maxwell_from_eigenmode(
    modes: dict,
    N_sol: int = 8000,
    r_mid: float = R_MID,
    r_outer: float = R_OUTER,
) -> dict:
    """Validate the Coulomb law from eigenfunction throat flux.

    Parameters
    ----------
    modes : dict
        Must contain key 1 with 'funcs' list (l=1 eigenmodes).
    N_sol : int
        Number of radial grid points for the BVP.

    Returns
    -------
    dict with Q, error metrics, solved and exact potential arrays.
    """
    fn1 = modes[1]["funcs"][0]
    du_dr = throat_du_dr(fn1, r_mid)
    Q = r_mid ** 2 * abs(du_dr)

    # Sparse BVP solve
    r = np.linspace(r_mid, r_outer, N_sol)
    dr = r[1] - r[0]
    r2 = r ** 2
    r2hp = (r + 0.5 * dr) ** 2
    r2hm = (r - 0.5 * dr) ** 2

    L = lil_matrix((N_sol, N_sol))
    for i in range(1, N_sol - 1):
        L[i, i] = -(r2hp[i] + r2hm[i]) / (r2[i] * dr ** 2)
        L[i, i - 1] = r2hm[i] / (r2[i] * dr ** 2)
        L[i, i + 1] = r2hp[i] / (r2[i] * dr ** 2)

    # Neumann at i=0
    L[0, 0] = -3 / (2 * dr)
    L[0, 1] = 4 / (2 * dr)
    L[0, 2] = -1 / (2 * dr)
    # Dirichlet at i=-1
    L[-1, -1] = 1.0

    rhs = np.zeros(N_sol)
    rhs[0] = -Q / r[0] ** 2
    rhs[-1] = 0.0

    L_csr = L.tocsr()
    A_sol = spsolve(L_csr, rhs)
    A_exact = Q / r - Q / r_outer

    diff = A_sol - A_exact
    err_max = float(np.abs(diff).max())
    err_rms = float(np.sqrt(np.mean(diff ** 2)))
    rel_err = err_max / float(np.abs(A_exact).max())

    # Linear-system residual: ||L·A_sol − rhs||_∞ / ||rhs||_∞
    # This is the correct acceptance metric for the sparse BVP solve.
    # The previous PDE residual (raw FD Laplacian of A_sol) was misleading:
    # it grows with N near the throat due to stencil cancellation, not
    # solver inaccuracy.  The linear-system residual is mesh-independent.
    linsys_res_vec = L_csr @ A_sol - rhs
    linsys_res_abs = float(np.abs(linsys_res_vec).max())
    rhs_norm = float(np.abs(rhs).max()) + 1e-30
    linsys_res = linsys_res_abs / rhs_norm

    # Compute radial electric field E_r = -dA/dr for downstream consumers
    # (e.g. transaction.handshake.advanced_confirm_amplitude).
    E_r = np.empty_like(A_sol)
    E_r[1:-1] = -(A_sol[2:] - A_sol[:-2]) / (2.0 * dr)
    E_r[0] = -(A_sol[1] - A_sol[0]) / dr
    E_r[-1] = -(A_sol[-1] - A_sol[-2]) / dr

    return dict(
        Q=Q,
        du_dr=du_dr,
        err_max=err_max,
        err_rms=err_rms,
        rel_err=rel_err,
        linsys_res=linsys_res,
        linsys_res_abs=linsys_res_abs,
        r=r,            # canonical key used by transaction subsystem
        r_grid=r,       # backward-compatible alias
        E_r=E_r,        # radial electric field
        A_solved=A_sol,
        A_exact=A_exact,
    )
