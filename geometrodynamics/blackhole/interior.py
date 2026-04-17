"""
Nonsingular black-hole interior from a wormhole-network core.

The standard Schwarzschild metric has f(r) = 1 - 2M/r → singular at r=0.
The wormhole-network core replaces the singularity with a regular de
Sitter-like interior parameterised by the core length scale l, which
is determined by the collective throat structure.

The effective metric function is the Hayward form:

    f(r) = 1 - 2M r² / (r³ + 2M l²)

Properties:
  - f(r) → 1 - 2M/r  for r >> l   (Schwarzschild exterior)
  - f(r) → 1 - r²/l²  for r << l   (de Sitter core, regular)
  - f(0) = 1  (no singularity)
  - Horizon exists at r_H ≈ 2M for l << M

The physical interpretation: the de Sitter core is the effective
geometry of the coherent wormhole-network interior, where throat
repulsion prevents collapse to a point.

For a Kerr-like extension, we use the Ghosh-Kumar rotating regular
metric that reduces to Hayward in the non-rotating limit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp


# ── Metric functions ─────────────────────────────────────────────────────────

def f_schwarzschild(r: float | np.ndarray, M: float) -> float | np.ndarray:
    """Standard Schwarzschild metric function f(r) = 1 - 2M/r."""
    r = np.asarray(r, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(r > 0, 1.0 - 2.0 * M / r, -np.inf)


def f_hayward(
    r: float | np.ndarray,
    M: float,
    l: float,
) -> float | np.ndarray:
    """Hayward regular metric function.

    f(r) = 1 - 2M r² / (r³ + 2M l²)

    Parameters
    ----------
    r : float or array
        Radial coordinate.
    M : float
        Mass parameter.
    l : float
        Core length scale from the wormhole-network interior.
        l → 0 recovers Schwarzschild.
    """
    r = np.asarray(r, dtype=float)
    denom = r ** 3 + 2.0 * M * l ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(denom > 0, 1.0 - 2.0 * M * r ** 2 / denom, 1.0)


def df_hayward_dr(
    r: float | np.ndarray,
    M: float,
    l: float,
) -> float | np.ndarray:
    """Radial derivative df/dr of the Hayward metric function.

    Used for surface gravity κ = ½|f'(r_H)| and geodesic equations.
    """
    r = np.asarray(r, dtype=float)
    r3 = r ** 3
    ml2 = 2.0 * M * l ** 2
    denom = r3 + ml2
    with np.errstate(divide="ignore", invalid="ignore"):
        num = -2.0 * M * r * (2.0 * denom - 3.0 * r ** 3) / denom ** 2
        return np.where(denom > 0, num, 0.0)


# ── Horizon finder ───────────────────────────────────────────────────────────

def find_horizons(
    M: float,
    l: float,
    r_min: float = 1e-6,
    r_max: float | None = None,
    n_scan: int = 2000,
) -> list[float]:
    """Find radii where f(r) = 0 (horizons) for the Hayward metric.

    Returns all horizon radii in ascending order.  For l << M there are
    two horizons (outer ≈ 2M, inner ≈ l²/M^{1/3}).  For l above a critical
    value l_crit, there are no horizons (no trapped surface).

    Parameters
    ----------
    r_max : float or None
        Upper bound of the scan window.  Default ``None`` auto-scales
        to ``max(10 * M, 50)`` so that the outer horizon at ≈ 2M is
        always inside the window regardless of mass.
    """
    if r_max is None:
        r_max = max(10.0 * M, 50.0)
    rs = np.linspace(r_min, r_max, n_scan)
    fs = f_hayward(rs, M, l)
    horizons = []
    for i in range(len(fs) - 1):
        if fs[i] * fs[i + 1] < 0:
            try:
                rh = brentq(lambda r: float(f_hayward(r, M, l)), rs[i], rs[i + 1])
                horizons.append(float(rh))
            except ValueError:
                pass
    return sorted(horizons)


def critical_core_scale(M: float) -> float:
    """Critical core scale l_crit above which no horizon forms.

    For the Hayward metric the degenerate-root condition gives
    l_crit ≈ 0.77 · M  (linear in M, found numerically).

    We compute it by binary search for the largest l at which
    ``find_horizons`` still returns at least one root.
    """
    def has_horizon(l_test):
        return len(find_horizons(M, l_test)) > 0

    l_lo = 0.0
    l_hi = 5.0 * M
    for _ in range(80):
        l_mid = 0.5 * (l_lo + l_hi)
        if has_horizon(l_mid):
            l_lo = l_mid
        else:
            l_hi = l_mid
    return 0.5 * (l_lo + l_hi)


# ── Surface gravity and temperature ─────────────────────────────────────────

def surface_gravity(M: float, l: float) -> float:
    """Surface gravity κ = ½ |f'(r_H)| at the outer horizon.

    Returns 0 if no horizon exists (l > l_crit).
    """
    horizons = find_horizons(M, l)
    if not horizons:
        return 0.0
    r_h = horizons[-1]  # outer horizon
    return 0.5 * abs(float(df_hayward_dr(r_h, M, l)))


def hawking_temperature(M: float, l: float) -> float:
    """Hawking temperature T = κ/(2π) in geometric units (ℏ = 1).

    For l → 0: T → 1/(8πM) (standard Schwarzschild).
    For l → l_crit: T → 0 (extremal regular BH, zero temperature).
    """
    return surface_gravity(M, l) / (2.0 * np.pi)


# ── Kretschner scalar (curvature invariant) ──────────────────────────────────

def kretschner_hayward(r: float | np.ndarray, M: float, l: float) -> float | np.ndarray:
    """Kretschner scalar K = R_{abcd} R^{abcd} for the Hayward metric.

    Key property: K is finite everywhere, including r = 0.
    For Schwarzschild (l=0): K = 48 M²/r⁶ → ∞ at r = 0.
    For Hayward (l>0): K(0) = 24/(l⁴) — finite.
    """
    r = np.asarray(r, dtype=float)
    # Numerical computation via finite differences of the metric
    eps = max(float(np.min(np.abs(r[r > 0]))) * 1e-4, 1e-10) if np.any(r > 0) else 1e-10

    def _f(x):
        return f_hayward(x, M, l)

    def _fp(x):
        return float(df_hayward_dr(x, M, l))

    def _fpp(x):
        return float((_fp(x + eps) - _fp(x - eps)) / (2 * eps))

    result = np.zeros_like(r, dtype=float)
    for i, ri in enumerate(np.atleast_1d(r)):
        if ri < eps:
            # Analytic limit: K(0) = 24/l⁴ for Hayward
            if l > 0:
                result[i] = 24.0 / l ** 4
            else:
                result[i] = np.inf
        else:
            fi = float(_f(ri))
            fpi = _fp(ri)
            fppi = _fpp(ri)
            # Kretschner for spherically symmetric metric:
            # K = (f'')² + 2(f'/r)² + 2((1-f)/r²)²   (simplified form)
            # Full expression for ds² = -f dt² + dr²/f + r² dΩ²
            term1 = fppi ** 2
            term2 = 2.0 * (fpi / ri) ** 2
            term3 = 2.0 * ((1.0 - fi) / ri ** 2) ** 2
            result[i] = term1 + term2 + term3

    if result.ndim == 0:
        return float(result)
    return result


# ── Radial geodesics ─────────────────────────────────────────────────────────

@dataclass
class RadialGeodesic:
    """Result of integrating a radial geodesic through the interior."""
    tau: np.ndarray        # proper time
    r: np.ndarray          # radial coordinate
    t: np.ndarray          # coordinate time (may diverge at horizon)
    is_complete: bool      # True if geodesic reaches r > 0 everywhere
    r_min: float           # minimum radius reached


def integrate_radial_geodesic(
    M: float,
    l: float,
    r_start: float,
    dr_start: float = -1.0,
    tau_max: float = 50.0,
    E: float = 1.0,
    L: float = 0.0,
) -> RadialGeodesic:
    """Integrate a radial timelike geodesic in the Hayward spacetime.

    Uses proper time τ as the affine parameter.  For l = 0 (Schwarzschild),
    infalling geodesics hit r = 0 in finite proper time.  For l > 0
    (Hayward/wormhole core), geodesics bounce at r ~ l and the spacetime
    is geodesically complete.

    Completeness is detected by checking whether the radial velocity
    changes sign (bounce) before r reaches a singularity floor.

    Parameters
    ----------
    r_start : float
        Initial radial coordinate (should be > r_H for infall).
    dr_start : float
        Initial dr/dτ (negative for infall).
    E : float
        Conserved energy per unit mass.
    L : float
        Conserved angular momentum (0 for radial).
    """
    # Clamp floor to prevent division by zero while preserving dynamics
    R_FLOOR = 1e-10

    def rhs(tau, y):
        r, rdot = y
        r_eff = max(r, R_FLOOR)
        fp = float(df_hayward_dr(r_eff, M, l))
        rddot = -0.5 * fp
        if L != 0.0:
            f = float(f_hayward(r_eff, M, l))
            rddot += L ** 2 * f / r_eff ** 3
        # If we hit the floor, stop infall
        if r <= R_FLOOR and rdot < 0:
            return [0.0, 0.0]
        return [rdot, rddot]

    # Set initial dr/dτ from energy conservation
    f0 = float(f_hayward(r_start, M, l))
    if dr_start == -1.0 and E ** 2 > f0:
        dr_start = -np.sqrt(E ** 2 - f0)

    y0 = [r_start, dr_start]
    sol = solve_ivp(
        rhs, [0, tau_max], y0,
        method="RK45",
        max_step=0.005,
        dense_output=True,
        rtol=1e-10, atol=1e-13,
    )

    r_arr = sol.y[0]
    rdot_arr = sol.y[1]
    r_min = float(np.min(r_arr))

    # Completeness criterion based on velocity at minimum radius.
    #
    # Schwarzschild (l=0): near r=0, f→-∞ so (dr/dτ)²→+∞.
    #   The geodesic crashes into the singularity with large |rdot|.
    #
    # Hayward (l>0): near r=0, f→1=E² so (dr/dτ)²→0.
    #   The geodesic asymptotically approaches r=0 with rdot→0,
    #   never actually reaching it (takes infinite proper time).
    #
    # Detection: at the point of smallest r, check if rdot is
    # still large (singular crash) or has decayed to ~0 (bounce/asymptote).
    i_min = int(np.argmin(r_arr))
    rdot_at_min = abs(rdot_arr[i_min])

    # For E=1 infall, the Schwarzschild singularity velocity scales
    # as |rdot| ~ √(2M/r) → large.  The Hayward asymptote gives
    # |rdot| ~ r/l → tiny.  Threshold: if |rdot| < 0.1 at r_min,
    # the geodesic is decelerating to a stop (complete).
    is_complete = rdot_at_min < 0.1

    # Reconstruct coordinate time by post-integration quadrature
    t_arr = np.zeros_like(sol.t)
    for i in range(1, len(sol.t)):
        dtau = sol.t[i] - sol.t[i - 1]
        r_mid = max(0.5 * (r_arr[i] + r_arr[i - 1]), R_FLOOR)
        f_mid = float(f_hayward(r_mid, M, l))
        if abs(f_mid) > 1e-10:
            t_arr[i] = t_arr[i - 1] + E / f_mid * dtau
        else:
            t_arr[i] = t_arr[i - 1]

    return RadialGeodesic(
        tau=sol.t,
        r=r_arr,
        t=t_arr,
        is_complete=is_complete,
        r_min=r_min,
    )


# ── Penrose diagram coordinates ──────────────────────────────────────────────

def tortoise_hayward(
    r: float | np.ndarray,
    M: float,
    l: float,
    dr: float = 0.001,
) -> float | np.ndarray:
    """Tortoise coordinate r* = ∫ dr/f(r) for the Hayward metric.

    Computed by numerical quadrature from r to a reference point.
    """
    r = np.atleast_1d(np.asarray(r, dtype=float))
    result = np.zeros_like(r)
    r_ref = 10.0 * M  # reference point where r* ≈ r

    for i, ri in enumerate(r):
        if ri <= 0:
            result[i] = -np.inf
            continue
        n_pts = max(int(abs(ri - r_ref) / dr), 100)
        rs = np.linspace(ri, r_ref, n_pts)
        fs = f_hayward(rs, M, l)
        # Avoid division by zero at horizons
        fs = np.where(np.abs(fs) < 1e-12, np.sign(fs) * 1e-12, fs)
        integrand = 1.0 / fs
        result[i] = -float(np.trapz(integrand, rs))

    if len(result) == 1:
        return float(result[0])
    return result
