"""
5D Tangherlini radial eigenmodes via Chebyshev spectral method.

The Chebyshev grid is laid in tortoise-coordinate (r*) space on the
interval [r*(R_MID + ε), r*(R_OUTER − ε)].  The effective potential
V(r, l) uses physical r obtained by inverting the tortoise map.
"""

import numpy as np
from scipy.optimize import brentq
from scipy.linalg import eig as scipy_eig

from geometrodynamics.constants import R_MID, R_OUTER


# ── Coordinate maps ─────────────────────────────────────────────────────────

def r_to_rstar(r: float, rs: float = R_MID) -> float:
    """5D tortoise coordinate r*(r) = r + (rs/2) ln|(r−rs)/(r+rs)|."""
    return r + rs / 2.0 * np.log(np.abs((r - rs) / (r + rs) + 1e-15))


def rstar_to_r(rstar: float, rs: float = R_MID, tol: float = 1e-12) -> float:
    """Invert the 5D tortoise coordinate via bisection."""
    def f(r):
        if r <= rs:
            return -1e30
        return r + rs / 2.0 * np.log(np.abs((r - rs) / (r + rs) + 1e-30)) - rstar

    r_lo = rs + 1e-10
    r_hi = max(abs(rstar) + rs + 10.0, 2.0 * rs)
    try:
        return brentq(f, r_lo, r_hi, xtol=tol)
    except Exception:
        return rs + 1e-8


# ── Effective potential ──────────────────────────────────────────────────────
#
# TWO OPERATORS LIVE HERE, AND THE DIFFERENCE IS A BUG, NOT A CONVENTION.
#
# `V_tangherlini_legacy` is what this module carried from its first commit until
# PR #271.  It is NOT the master potential of a minimally coupled massless
# scalar, and its generic old name implied that it was.
#
# `V_scalar_tangherlini` is that master potential, derived rather than recalled
# (PR #270, `tangherlini.dynamics.derive_the_field_equations`) and agreeing
# bit-for-bit with `dynamics.master_potential`.  For
#
#     ds² = −A dt² + A⁻¹dr² + r²dΩ_n² ,     A = 1 − (r_h/r)^{n−1}
#
# a minimally coupled massless scalar Φ = e^{−iωt}Y_ℓ(Ω)R(r) obeys
#
#     (1/r^n) ∂_r(r^n A R') + (ω²/A − ℓ(ℓ+n−1)/r²) R = 0 ,
#
# and the unique first-derivative-free Schrödinger form, with dr* = dr/A and
# ψ = r^{n/2}R, carries
#
#     V = A[ ℓ(ℓ+n−1)/r² + n(n−2)A/(4r²) + n A'/(2r) ] .
#
# At n = 3 that is A[(ℓ(ℓ+2) + 3/4)/r² + (9/4)r_h²/r⁴], so
#
#     V_scalar − V_legacy = 3A²/(4r²) ,
#
# an ℓ-INDEPENDENT shift.  Two consequences worth keeping in view:
#   * differences V_{ℓ₂} − V_{ℓ₁} are unchanged EXACTLY, so any result built on
#     the cross-ℓ operator survives algebraically;
#   * barrier heights and their sums are not, so γ-type quantities move.
#
# The flat limit settles which is which with no appeal to authority: at
# r_h → 0 the regular solution is φ = J_{ℓ+1}(ωr)/r, so ψ = r^{1/2}J_{ℓ+1}(ωr),
# and Bessel's equation gives V → ((ℓ+1)² − ¼)/r² = (ℓ(ℓ+2) + 3/4)/r².
#
# `V_tangherlini` now delegates to the corrected operator.  Archived runs are
# reproduced by calling `V_tangherlini_legacy` explicitly.

def V_tangherlini_legacy(r: float | np.ndarray, l: int, rs: float = R_MID):
    """**Frozen.** The pre-#271 potential, kept only to reproduce archived runs.

    ``V(r, l) = A(r)·[l(l+2)/r² + 3·rs²/r⁴]``, ``A = 1 − (rs/r)²``.

    This is *not* the master potential of a minimally coupled massless scalar
    with ``ψ = r^{3/2}R``; it is short of that by ``3A²/(4r²)``.  Do not use it
    for new work — see `V_scalar_tangherlini`.
    """
    f = 1.0 - (rs / r) ** 2
    return f * (l * (l + 2) / r ** 2 + 3.0 * rs ** 2 / r ** 4)


def V_scalar_tangherlini(r: float | np.ndarray, l: int, rs: float = R_MID,
                         n: int = 3):
    """The master potential of a minimally coupled massless scalar, general ``n``.

    ``V = A[ l(l+n−1)/r² + n(n−2)A/(4r²) + n A'/(2r) ]`` with
    ``A = 1 − (rs/r)^{n−1}`` and ``ψ = r^{n/2}R``.

    Verified symbolically at ``n = 2 … 6`` and against the flat-limit Bessel
    form; identical bit-for-bit to `tangherlini.dynamics.master_potential`.
    """
    r = np.asarray(r, dtype=float) if not np.isscalar(r) else float(r)
    A = 1.0 - (rs / r) ** (n - 1)
    A_prime = (n - 1) * rs ** (n - 1) / r ** n
    return A * (l * (l + n - 1) / r ** 2
                + (n * (n - 2) / 4.0) * A / r ** 2
                + (n / 2.0) * A_prime / r)


def V_tangherlini(r: float | np.ndarray, l: int, rs: float = R_MID):
    """The canonical radial potential — the minimally coupled scalar one.

    Delegates to `V_scalar_tangherlini`.  Before PR #271 this name carried the
    operator now called `V_tangherlini_legacy`, which differs by ``3A²/(4r²)``;
    `experiments/closure_ledger/scalar_operator_audit_probe.py` records what
    that changed and what it left exactly invariant.
    """
    return V_scalar_tangherlini(r, l, rs)


# ── Chebyshev differentiation ───────────────────────────────────────────────

def _cheb_diff(N: int):
    """Chebyshev differentiation matrix on N+1 Gauss–Lobatto nodes."""
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1)
    c[0] = 2
    c[N] = 2
    c *= (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1))
    dX = X - X.T
    D = (c[:, None] / c[None, :]) / (dX + np.eye(N + 1))
    D -= np.diag(np.sum(D, axis=1))
    return x, D


# ── Eigenvalue solver ────────────────────────────────────────────────────────

def solve_radial_modes(
    l: int,
    N: int = 80,
    n_modes: int = 4,
    rs: float = R_MID,
    r_outer: float = R_OUTER,
) -> tuple[np.ndarray, list[dict], np.ndarray]:
    """Chebyshev eigenvalue solver for 5D Tangherlini radial modes.

    Parameters
    ----------
    l : int
        Angular momentum quantum number (λ = l(l+2) on S³).
    N : int
        Number of Chebyshev interior nodes.
    n_modes : int
        Number of lowest eigenfrequencies to return.

    Returns
    -------
    omegas : ndarray of shape (n_modes,)
        Eigenfrequencies ω_n.
    funcs : list of dict
        Each dict contains 'u_half', 'r_half', 'u_full', 'r_full', 'r_phys'.
    r_grid : ndarray
        Physical-r Chebyshev grid (outer half).
    """
    rs_min = r_to_rstar(rs + 5e-4, rs)
    rs_max = r_to_rstar(r_outer - 5e-4, rs)
    x, D = _cheb_diff(N)
    D2 = D @ D
    L = (rs_max - rs_min) / 2.0
    rsg = rs_min + L * (1.0 - x)
    rg = np.array([rstar_to_r(s, rs) for s in rsg])
    Vg = V_tangherlini(rg, l, rs)

    H = -(1.0 / L ** 2) * D2 + np.diag(Vg)
    H_int = H[1:N, 1:N]
    ev, evec = scipy_eig(H_int)
    ev = np.real(ev)
    evec = np.real(evec)

    pos = np.where(ev > 0)[0]
    ev = ev[pos]
    evec = evec[:, pos]
    idx = np.argsort(ev)
    n_ret = min(n_modes, len(idx))
    ev = ev[idx[:n_ret]]
    evec = evec[:, idx[:n_ret]]

    oms = np.sqrt(ev)
    funcs = []
    for k in range(n_ret):
        u = np.zeros(N + 1)
        u[1:N] = evec[:, k]
        if abs(u.min()) > u.max():
            u = -u
        u /= abs(u).max() + 1e-12
        r_full = np.concatenate([2 * rs - rg[::-1], rg])
        u_full = np.concatenate([-u[::-1], u])
        funcs.append(
            dict(u_half=u, r_half=rg, u_full=u_full, r_full=r_full, r_phys=rg)
        )

    return oms, funcs, rg
