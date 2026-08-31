"""The finite-mouth scalar-flat handle in a round `S⁴_R` spatial bulk.

WHAT IS ASSUMED, AND WHAT IS FORCED
───────────────────────────────────
**One** geometric assumption: the observed closed `S³` universe is the totally
geodesic equator of a round four-sphere spatial bulk ``Sigma_bulk = S^4_R``. It
adds no scale and no dimensionless shape parameter.

Everything else here is forced. Excise two geodesic four-balls of angular
radius ``a < pi/2`` and join their ``S^3`` boundaries by ``[-S,S] x S^3`` with

    ds^2 = -N(s)^2 dt^2 + ds^2 + f(s)^2 dOmega_3^2

Demanding ``{}^{(4)}R = 0`` on the time-symmetric slice gives ``f f'' = 1-f'^2``,
one integration gives ``f'^2 = 1 - b^2/f^2``, and hence ``f = sqrt(s^2+b^2)``.
Darmois matching to the exterior (``f_m = R sin a``, ``f'_+ = cos a``) then fixes
**both** remaining constants at once:

    S = R sin a cos a ,   b = R sin^2 a ,   L = 2S = R sin 2a

There is no independent neck radius, throat length, or tube area to choose.
That is the whole point: PR #263-#265 spent three rounds discovering that a
chosen tube area was carrying the answer, and here there is nothing to choose.

THE SPATIAL GEOMETRY IS SHARED; ONLY THE LAPSE DIFFERS
──────────────────────────────────────────────────────
``ds^2 + (s^2+b^2)dOmega_3^2`` is *also* a time-symmetric slice of 5D
Schwarzschild-Tangherlini with ``F(r) = 1 - b^2/r^2`` and ``r^2 = s^2+b^2``. So
the repository's two throat pictures are not different spatial geometries — they
are one spatial geometry with two lapses:

    N_vac(s)  = |s|/sqrt(s^2+b^2)   horizon at the neck, vacuum, non-traversable
    N_trav(s) = 1                   ultrastatic, traversable, NEC-violating

``lapse_vacuum`` and ``lapse_ultrastatic`` below are exactly those, evaluated on
the *same* ``throat_radius``. The physical fork is the single number ``N(0)``.

WHY A DtN OPERATOR AND NOT AN S-MATRIX
──────────────────────────────────────
There is no infinity here, so a transmission amplitude defined against Jost
phases at ``s = +-infinity`` is not the physical observable — and its constant
phase is exactly what dissolved PR #276's closure verdict. The physical object
is the Dirichlet-to-Neumann map at the two *actual* ``S^3`` mouth surfaces,

    Y_l(omega) : (phi_A, phi_B) -> (q_A, q_B) ,  q = 2 pi^2 f_m^3 n^s d_s phi

which gives the field an absolute reference surface. Closure then becomes a
finite-boundary determinant rather than an asymptotic phase test, and a change
of mode basis conjugates the operator instead of shifting a free constant.

WHAT THIS MODULE DELIBERATELY DOES NOT DO
─────────────────────────────────────────
It does not touch the discrete BAM identification. ``Phi_spatial``, ``(-1)^l``,
``eta_orientation``, ``eta_wrap`` and ``U_spin`` are five separate operations
and none of them is folded into ``f(s)``, a sign convention, or a harmonic
phase here. They may be multiplied into a boundary condition only after the
boundary lift is fixed — see ``docs/finite_mouth_prereg.md``.

Predictions P1-P6 are frozen in that document, which was committed *before*
this module existed.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

__all__ = [
    "BULK_RADIUS",
    "MOUTH_ANGLE",
    "ThroatGeometry",
    "geometry",
    "throat_radius",
    "throat_radius_derivative",
    "spatial_ricci_scalar",
    "misner_sharp",
    "lapse_ultrastatic",
    "lapse_vacuum",
    "stress_tensor",
    "extrinsic_curvature_angular",
    "junction_jumps",
    "null_energy_at_neck",
    "null_energy_integral",
    "static_admittance",
    "monopole_conductance",
    "static_resistance",
    "solve_admittance",
]

#: Bulk four-sphere radius. Sets the only scale.
BULK_RADIUS = 1.0

#: Mouth angular radius on the bulk `S^4`. The only shape parameter.
MOUTH_ANGLE = 0.30


class ThroatGeometry:
    """The forced geometry of one finite handle, from ``R`` and ``a`` alone.

    Every attribute is derived. Nothing here is settable, because nothing in
    the construction is free once the ``S^4_R`` completion is assumed.
    """

    __slots__ = ("bulk_radius", "mouth_angle", "neck_radius", "half_length",
                 "mouth_radius", "proper_length", "rapidity")

    def __init__(self, bulk_radius: float = BULK_RADIUS,
                 mouth_angle: float = MOUTH_ANGLE) -> None:
        if not 0.0 < mouth_angle < 0.5 * math.pi:
            raise ValueError("the mouth angular radius must satisfy "
                             f"0 < a < pi/2; got {mouth_angle}")
        if bulk_radius <= 0.0:
            raise ValueError(f"the bulk radius must be positive; got {bulk_radius}")
        self.bulk_radius = float(bulk_radius)
        self.mouth_angle = float(mouth_angle)
        # P1 -- forced by Darmois matching, not chosen
        self.neck_radius = bulk_radius * math.sin(mouth_angle) ** 2   # b
        self.half_length = bulk_radius * math.sin(mouth_angle) * math.cos(mouth_angle)
        self.mouth_radius = bulk_radius * math.sin(mouth_angle)       # f_m
        self.proper_length = 2.0 * self.half_length                   # L = R sin 2a
        # X = arcosh(1/sin a); the natural throat rapidity, with tanh X = cos a
        self.rapidity = math.acosh(1.0 / math.sin(mouth_angle))

    def __repr__(self) -> str:                                    # pragma: no cover
        return (f"ThroatGeometry(R={self.bulk_radius:g}, a={self.mouth_angle:g}, "
                f"b={self.neck_radius:.6g}, L={self.proper_length:.6g})")


@lru_cache(maxsize=32)
def geometry(bulk_radius: float = BULK_RADIUS,
             mouth_angle: float = MOUTH_ANGLE) -> ThroatGeometry:
    """Cached :class:`ThroatGeometry`."""
    return ThroatGeometry(bulk_radius, mouth_angle)


# ── the spatial profile ─────────────────────────────────────────────────────

def throat_radius(s, neck_radius: float):
    """``f(s) = sqrt(s^2 + b^2)`` — the unique scalar-flat profile."""
    s = np.asarray(s, dtype=float)
    return np.sqrt(s * s + neck_radius * neck_radius)


def throat_radius_derivative(s, neck_radius: float):
    """``f'(s) = s/f``, so ``f'(0) = 0``: the neck is where ``f'`` vanishes."""
    s = np.asarray(s, dtype=float)
    return s / throat_radius(s, neck_radius)


def spatial_ricci_scalar(s, neck_radius: float):
    """``{}^{(4)}R = -6 f''/f + 6(1-f'^2)/f^2``, which must vanish identically.

    Computed from the closed forms rather than from the profile it defines, so
    that it is a genuine check and not a tautology: ``f'' = b^2/f^3``.
    """
    s = np.asarray(s, dtype=float)
    f = throat_radius(s, neck_radius)
    fp = throat_radius_derivative(s, neck_radius)
    fpp = neck_radius * neck_radius / f ** 3
    return -6.0 * fpp / f + 6.0 * (1.0 - fp * fp) / (f * f)


def misner_sharp(s, neck_radius: float):
    """``mu = f^2(1 - f'^2)``, the Tangherlini/Misner-Sharp mass parameter.

    Constant ``= b^2`` throughout the throat, which is P2's inside half.
    """
    s = np.asarray(s, dtype=float)
    f = throat_radius(s, neck_radius)
    fp = throat_radius_derivative(s, neck_radius)
    return f * f * (1.0 - fp * fp)


def misner_sharp_exterior(chi, bulk_radius: float):
    """``mu_ext(chi) = R^2 sin^4 chi`` on the round ``S^4_R``."""
    chi = np.asarray(chi, dtype=float)
    return bulk_radius ** 2 * np.sin(chi) ** 4


# ── the two lapses: the whole physical fork ─────────────────────────────────

def lapse_ultrastatic(s, neck_radius: float):
    """``N = 1``. Traversable, and joins the ultrastatic exterior shell-free."""
    return np.ones_like(np.asarray(s, dtype=float))


def lapse_vacuum(s, neck_radius: float):
    """``N = |s|/sqrt(s^2+b^2)``, the Tangherlini lapse. ``N(0) = 0``.

    Vacuum and non-traversable: a Killing horizon sits at the neck. Note this
    is applied to the **same** ``throat_radius`` as the ultrastatic branch —
    the two pictures share their spatial geometry exactly.
    """
    s = np.asarray(s, dtype=float)
    return np.abs(s) / throat_radius(s, neck_radius)


# ── stress, from the general-lapse Einstein equations ───────────────────────

def stress_tensor(s, neck_radius: float,
                  lapse: Optional[Callable] = None,
                  lapse_derivative: Optional[Callable] = None,
                  lapse_second: Optional[Callable] = None,
                  step: float = 1e-5) -> Dict[str, np.ndarray]:
    """Orthonormal-frame stress for an arbitrary lapse on the forced profile.

        8 pi G_5 rho   = 0                                (any lapse)
        8 pi G_5 p_s   = 3 f'N'/(fN) - 3 b^2/f^4
        8 pi G_5 p_Om  = N''/N + 2 f'N'/(fN) + b^2/f^4

    ``rho = 0`` is not an assumption but Gauss-Codazzi on a time-symmetric
    slice: ``16 pi rho = {}^{(4)}R`` with ``K = 0``, and the profile was chosen
    to make ``{}^{(4)}R`` vanish. It therefore holds for **every** lapse, which
    is what makes the neck NEC result below unavoidable.

    Derivatives of the lapse default to centred differences, so a caller may
    supply any smooth ``N`` without also supplying its derivatives.
    """
    s = np.asarray(s, dtype=float)
    if lapse is None:
        lapse = lapse_ultrastatic
    f = throat_radius(s, neck_radius)
    fp = throat_radius_derivative(s, neck_radius)
    n = np.asarray(lapse(s, neck_radius), dtype=float)
    if lapse_derivative is not None:
        dn = np.asarray(lapse_derivative(s, neck_radius), dtype=float)
    else:
        dn = (np.asarray(lapse(s + step, neck_radius), dtype=float)
              - np.asarray(lapse(s - step, neck_radius), dtype=float)) / (2.0 * step)
    if lapse_second is not None:
        ddn = np.asarray(lapse_second(s, neck_radius), dtype=float)
    else:
        ddn = (np.asarray(lapse(s + step, neck_radius), dtype=float)
               - 2.0 * n
               + np.asarray(lapse(s - step, neck_radius), dtype=float)) / (step * step)
    with np.errstate(divide="ignore", invalid="ignore"):
        shear = np.where(n != 0.0, fp * dn / (f * n), 0.0)
        curve = np.where(n != 0.0, ddn / n, 0.0)
    quartic = f ** 4
    return {
        "density": np.zeros_like(f),
        "radial_pressure": 3.0 * shear - 3.0 * neck_radius ** 2 / quartic,
        "angular_pressure": curve + 2.0 * shear + neck_radius ** 2 / quartic,
        "radial_nec": 3.0 * shear - 3.0 * neck_radius ** 2 / quartic,
    }


def null_energy_at_neck(neck_radius: float) -> float:
    """``8 pi G_5 (rho + p_s)|_0 = -3/b^2``, for **every** smooth ``N(0) > 0``.

    The lapse cannot help: its only appearance in ``p_s`` is through
    ``3 f'N'/(fN)``, and ``f'(0) = 0`` is precisely what makes ``s = 0`` a neck.
    So the term vanishes there whatever ``N'`` does — the result needs no
    reflection symmetry, which is stronger than requiring ``N'(0) = 0``.

    Smooth **and** traversable therefore implies radial NEC violation at the
    neck. The only escape is ``N(0) = 0``: the Tangherlini horizon branch.
    """
    return -3.0 / (neck_radius * neck_radius)


def null_energy_integral(bulk_radius: float, mouth_angle: float) -> float:
    """``8 pi G_5 int T_ab k^a k^b dlambda`` across the finite ultrastatic throat.

        = -(3/R)[ cot a + (pi/2 - a) csc^2 a ]

    Finite, unlike PR #276's two-infinite-end version. For small mouths it
    diverges as ``-3 pi/(2 R a^2)``, i.e. the point-mouth limit is singularly
    expensive — which is consistent with the repository's independent finding
    that the point-mouth constraint problem is singular and that finite mouths
    are not merely numerical regulators.
    """
    a = mouth_angle
    return -3.0 * (1.0 / math.tan(a)
                   + (0.5 * math.pi - a) / math.sin(a) ** 2) / bulk_radius


# ── the junction ────────────────────────────────────────────────────────────

def extrinsic_curvature_angular(f_value: float, f_prime: float) -> float:
    """``K^A_B = (f'/f) delta^A_B`` — the angular block, as a scalar."""
    return f_prime / f_value


def junction_jumps(bulk_radius: float = BULK_RADIUS,
                   mouth_angle: float = MOUTH_ANGLE) -> Dict[str, object]:
    """P3 — the seam carries no Israel layer.

    Compares the induced metric and extrinsic curvature computed *separately*
    from the throat side and the exterior side. The throat side is evaluated
    from ``f, f'`` at ``s = S``; the exterior side from ``R sin chi, cos chi``
    at ``chi = a``. Nothing is shared between the two evaluations, so agreement
    is a check rather than an identity.
    """
    g = geometry(bulk_radius, mouth_angle)
    s_edge = g.half_length
    f_in = float(throat_radius(s_edge, g.neck_radius))
    fp_in = float(throat_radius_derivative(s_edge, g.neck_radius))
    f_out = bulk_radius * math.sin(mouth_angle)
    fp_out = math.cos(mouth_angle)
    k_in = extrinsic_curvature_angular(f_in, fp_in)
    k_out = extrinsic_curvature_angular(f_out, fp_out)
    # normal pressure on each side, ultrastatic
    p_in = -3.0 * g.neck_radius ** 2 / f_in ** 4
    p_out = -3.0 / bulk_radius ** 2
    return {
        "areal_radius_inside": f_in,
        "areal_radius_outside": f_out,
        "induced_metric_jump": abs(f_in - f_out),
        "extrinsic_curvature_inside": k_in,
        "extrinsic_curvature_outside": k_out,
        "extrinsic_curvature_jump": abs(k_in - k_out),
        "lapse_jump": 0.0,                       # ultrastatic on both sides
        "surface_stress_vanishes": bool(abs(f_in - f_out) < 1e-12
                                        and abs(k_in - k_out) < 1e-12),
        "normal_pressure_inside": p_in,
        "normal_pressure_outside": p_out,
        "normal_pressure_jump": abs(p_in - p_out),
        "exterior_density": 6.0 / bulk_radius ** 2,
        "second_derivative_jumps": (
            "f'' = b^2/f^3 inside and -R sin chi outside do NOT agree, so the "
            "geometry is C^1 and not C^2. That is a finite STEP in bulk stress, "
            "not a delta function: the Israel layer depends on [K_ab], which "
            "vanishes. Tangential pressure and density jump, which is allowed "
            "because they involve second normal derivatives."),
    }


# ── the static Dirichlet-to-Neumann map ─────────────────────────────────────

def static_admittance(ell: int, bulk_radius: float = BULK_RADIUS,
                      mouth_angle: float = MOUTH_ANGLE) -> np.ndarray:
    """P5 — the closed-form finite-mouth admittance ``Y_l(0)``.

    With ``s = b sinh x`` and ``u = y/cosh x`` the static radial equation
    collapses to ``y'' = k^2 y`` with ``k = l+1``, giving

        Y_l(0) = 2 pi^2 F^2 [[ k coth 2kX - cos a ,  -k csch 2kX ],
                             [ -k csch 2kX       ,   k coth 2kX - cos a ]]

    for ``X = arcosh(1/sin a)`` and ``F = R sin a``. The flux convention is
    ``q = 2 pi^2 f_m^3 n^s d_s phi`` with ``n^s`` the outward normal, so ``q``
    is the flux **out of the throat** at each mouth.
    """
    g = geometry(bulk_radius, mouth_angle)
    k = float(ell + 1)
    x = g.rapidity
    prefactor = 2.0 * math.pi ** 2 * g.mouth_radius ** 2
    diagonal = k / math.tanh(2.0 * k * x) - math.cos(mouth_angle)
    off = -k / math.sinh(2.0 * k * x)
    return prefactor * np.array([[diagonal, off], [off, diagonal]])


def monopole_conductance(bulk_radius: float = BULK_RADIUS,
                         mouth_angle: float = MOUTH_ANGLE) -> float:
    """``G = pi^2 R^2 sin^4 a / cos a`` — the ``l = 0`` admittance scale."""
    return (math.pi ** 2 * bulk_radius ** 2 * math.sin(mouth_angle) ** 4
            / math.cos(mouth_angle))


def static_resistance(bulk_radius: float = BULK_RADIUS,
                      mouth_angle: float = MOUTH_ANGLE) -> float:
    """``I_3 = int ds/f^3 = 2 cos a/(R^2 sin^4 a)``, with ``G = 2 pi^2 / I_3``."""
    return (2.0 * math.cos(mouth_angle)
            / (bulk_radius ** 2 * math.sin(mouth_angle) ** 4))


def solve_admittance(ell: int, bulk_radius: float = BULK_RADIUS,
                     mouth_angle: float = MOUTH_ANGLE,
                     steps: int = 4000) -> np.ndarray:
    """The same ``Y_l(0)``, built numerically and sharing no algebra with it.

    Solves the static radial equation ``(f^3 u')' = l(l+2) f u`` as a two-point
    Dirichlet problem for the two independent boundary data ``(1,0)`` and
    ``(0,1)``, then reads off the outward fluxes. It never uses the
    ``sinh``/``cosh`` reduction that produced the closed form, so agreement is
    a regression rather than an identity.

    **Why a direct BVP solve and not shooting.** The two homogeneous solutions
    behave as ``e^{+-kx}`` with ``k = l+1`` across a rapidity span ``2kX``; at
    ``l = 5, a = 0.3`` that is ``e^23 ~ 1e10``. A shooting basis marched from
    one end is swamped by the growing mode, and the cancellation needed to
    impose the far boundary condition loses about ten digits — an error that
    gets *worse* with refinement rather than better, which is how it was
    caught. The conservative tridiagonal discretisation below has no growing
    mode to amplify and converges cleanly at second order.
    """
    g = geometry(bulk_radius, mouth_angle)
    lam = float(ell * (ell + 2))
    s_grid = np.linspace(-g.half_length, g.half_length, steps + 1)
    h = float(s_grid[1] - s_grid[0])
    f_node = throat_radius(s_grid, g.neck_radius)
    f_half = throat_radius(0.5 * (s_grid[:-1] + s_grid[1:]), g.neck_radius)
    w_half = f_half ** 3                      # w_{i+1/2}

    # interior unknowns u_1..u_{M-1}; conservative stencil
    #   [w_{i+1/2}(u_{i+1}-u_i) - w_{i-1/2}(u_i-u_{i-1})]/h^2 = lam f_i u_i
    m = steps - 1
    lower = w_half[:-1][1:] / h ** 2                       # coeff of u_{i-1}
    upper = w_half[1:][:-1] / h ** 2                       # coeff of u_{i+1}
    diag = -(w_half[:-1] + w_half[1:]) / h ** 2 - lam * f_node[1:-1]

    banded = np.zeros((3, m))
    banded[0, 1:] = upper
    banded[1, :] = diag
    banded[2, :-1] = lower

    from scipy.linalg import solve_banded
    columns = []
    for phi_a, phi_b in ((1.0, 0.0), (0.0, 1.0)):
        rhs = np.zeros(m)
        rhs[0] -= w_half[0] / h ** 2 * phi_a
        rhs[-1] -= w_half[-1] / h ** 2 * phi_b
        interior = solve_banded((1, 1), banded, rhs)
        u = np.concatenate(([phi_a], interior, [phi_b]))
        # v = f^3 u' at the half points, then half a cell of v' = lam f u to
        # reach the boundary itself -- second order, and conservative.
        v_left = w_half[0] * (u[1] - u[0]) / h - 0.5 * h * lam * f_node[0] * u[0]
        v_right = w_half[-1] * (u[-1] - u[-2]) / h + 0.5 * h * lam * f_node[-1] * u[-1]
        # q = 2 pi^2 f^3 u' n^s: outward normal is -d_s at s=-S, +d_s at s=+S
        columns.append([-2.0 * math.pi ** 2 * v_left,
                        2.0 * math.pi ** 2 * v_right])
    return np.array(columns).T
