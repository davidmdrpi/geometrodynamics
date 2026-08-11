"""
The continuous ℝ³ realisation of a spin-2 wave on S².
=====================================================

`spin2_tidal.py` samples the tensor as discrete tangent-space ellipses.  That
is faithful but flat: it never touches the embedding, so it gives no intuition
for how tidal shear meets a bulk.  This module projects `h_ab` into a
**continuous surface deformation** `r(θ, φ, t)` instead.

The obstruction, and the way through it
───────────────────────────────────────
A purely radial deformation cannot do it.  For `X = r n̂` the induced metric is
`g_ab = r² ĝ_ab + ∂_a r ∂_b r`, and the gradient term is *second* order, so at
first order

    δg_ab = 2ρ ĝ_ab      —  purely **conformal**.

Shape alone carries only the trace.  A trace-free `h_ab` therefore cannot be a
height field, which is exactly why the ellipse sampling existed.

The full displacement has a tangential part as well, and then

    X = (1 + ερ) n̂ + ε ξ,     δg_ab = 2ρ ĝ_ab + 2∇₍ₐξ_b₎ .

Setting `ξ = ∇W` and demanding tracelessness fixes the radial part completely:

    ρ = −½ ΔW ,        h_ab = [2∇₍ₐ∇_b₎W]^TF ,

so **one potential carries both**: its trace-free Hessian is the shear you
measure, and minus half its Laplacian is the radial displacement you can see.
The shear lives in how the material slides; the shape carries what is left.

Solving for the potential
─────────────────────────
For an axisymmetric field `h(d)` the Hessian condition `W'' − cot d W' = h` is
first order in `ψ ≡ W'`, with integrating factor `sin d`:

    ψ(d) = sin d · [ C + ∫₀^d h(x)/sin x dx ] ,
    ρ(d) = −½h(d) − cos d · [ C + ∫₀^d h(x)/sin x dx ] .

No derivative of the solved field is ever taken — only one integral, and its
integrand is regular because `h ~ sin²d` at both poles.

What the free constant is
─────────────────────────
`C` is the whole kernel of the construction: it is the `ℓ = 1` mode, and an
`ℓ = 1` displacement is a **rigid translation** of the sphere.  Fixing `C` by
removing the dipole makes the deformation unique up to nothing.  The `ℓ = 0`
mode cannot appear at all — `∫ΔW dA = 0` — so a spin-2 wave can never breathe
the sphere's area.

What comes out
──────────────
`ℓ = 2` gives `ρ = P₂(cos d)` exactly: the quadrupole tidal wave *is* the
prolate–oblate deformation of the textbook picture, derived rather than drawn.
And the deformed surface has the same area as the round one to first order,
because that is what trace-free means.

Scope
─────
Still a spin-2 field on a fixed background — the deformation is a faithful
*representation* of `h_ab` in the embedding, not backreaction.  What it buys is
that the wave now has an extrinsic amplitude, so it can be watched reaching
into a bulk.
"""

from __future__ import annotations

import math
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from geometrodynamics.viz.spin2_tidal import (
    ANTIPODAL_TIME,
    RETURN_TIME,
    Spin2WaveSim,
)
from geometrodynamics.viz.warped_sphere import NestedShells

__all__ = [
    "EmbeddedTidalSurface",
    "MaterialPatch",
    "measure_patch_area_invariance",
    "measure_patch_shape_history",
    "measure_induced_metric",
    "measure_quadrupole_shape",
    "measure_area_and_multipoles",
    "measure_bulk_reach",
]


class EmbeddedTidalSurface:
    """A spin-2 wave drawn as a continuous deformation of the embedded sphere.

    Parameters
    ----------
    shells
        The Russian-doll bulk the deformation reaches into.  Defaults to the
        ``radial_caustic`` vacuole, as in ``warped_sphere``.
    gain
        Display amplitude ``ε``.  ``None`` calibrates it from the run's own
        peak radial profile so the surface uses a fixed fraction of the gap.
    """

    def __init__(self, sim: Optional[Spin2WaveSim] = None,
                 shells: Optional[NestedShells] = None,
                 source_theta: float = 0.5 * math.pi,
                 source_phi: float = 0.0,
                 n_theta: int = 121, n_phi: int = 181,
                 gain: Optional[float] = None,
                 fill: float = 0.75, **sim_kwargs) -> None:
        self.sim = sim if sim is not None else Spin2WaveSim(**sim_kwargs)
        self.shells = shells or NestedShells()
        self.source_theta = float(source_theta)
        self.source_phi = float(source_phi)
        self.n_theta = int(n_theta)
        self.n_phi = int(n_phi)

        s = np.array([math.sin(self.source_theta) * math.cos(self.source_phi),
                      math.sin(self.source_theta) * math.sin(self.source_phi),
                      math.cos(self.source_theta)])
        self.source_direction = s
        tmp = np.array([0.0, 0.0, 1.0])
        if abs(float(np.dot(tmp, s))) > 0.9:
            tmp = np.array([1.0, 0.0, 0.0])
        e1 = tmp - float(np.dot(tmp, s)) * s
        self._e1 = e1 / np.linalg.norm(e1)
        self._e2 = np.cross(s, self._e1)

        # render lattice in geodesic-polar coordinates about the source
        self.d_grid = np.linspace(0.0, math.pi, self.n_theta)
        self.psi_grid = np.linspace(0.0, 2.0 * math.pi, self.n_phi)
        self._D, self._A = np.meshgrid(self.d_grid, self.psi_grid, indexing="ij")

        self.fill = float(fill)
        self.gain = float(gain) if gain is not None else self._calibrate()

    # ── clock ───────────────────────────────────────────────────────────────
    @property
    def t(self) -> float:
        return self.sim.t

    def reset(self) -> None:
        self.sim.reset()

    def advance_to(self, t_target: float) -> None:
        self.sim.advance_to(t_target)

    def _calibrate(self, samples: int = 200,
                   t_end: float = RETURN_TIME) -> float:
        """Fix ε so the peak radial excursion fills ``fill`` of the half-gap."""
        self.sim.reset()
        peak = 0.0
        for i in range(samples):
            self.sim.advance_to((i + 1) * t_end / samples)
            peak = max(peak, float(np.max(np.abs(self.radial_profile()))))
        self.sim.reset()
        self.peak_radial_profile = peak
        return (self.fill * self.shells.delta / peak) if peak > 0.0 else 1.0

    # ── the potential and its two parts ─────────────────────────────────────
    def _integral(self) -> np.ndarray:
        """``∫₀^d h/sin`` on the solver grid — regular, since ``h ~ sin²d``."""
        d = self.sim.d
        f = self.sim.h / np.sin(d)
        step = self.sim.dd
        out = np.empty_like(f)
        out[0] = 0.5 * f[0] * step          # from the pole to the first centre
        out[1:] = out[0] + np.cumsum(0.5 * (f[1:] + f[:-1]) * step)
        return out

    def _dipole_constant(self, integral: np.ndarray) -> float:
        """The ``C`` that removes the rigid translation from the deformation."""
        d = self.sim.d
        w = np.sin(d) * self.sim.dd                      # measure on the sphere
        cos = np.cos(d)
        # ⟨ρ, cos⟩ = 0 with ρ = −h/2 − cos·(C + I)
        num = -0.5 * float(np.sum(self.sim.h * cos * w)) \
            - float(np.sum(cos ** 2 * integral * w))
        den = float(np.sum(cos ** 2 * w))
        return num / den if den > 0.0 else 0.0

    def profiles(self) -> Dict[str, np.ndarray]:
        """``h`` (shear), ``ψ`` (tangential slide) and ``ρ`` (radial), on one grid."""
        d = self.sim.d
        h = self.sim.h
        integral = self._integral()
        c = self._dipole_constant(integral)
        bracket = c + integral
        return {
            "d": d,
            "shear": h,
            "slide": np.sin(d) * bracket,
            "radial": -0.5 * h - np.cos(d) * bracket,
            "dipole_constant": c,
        }

    def radial_profile(self) -> np.ndarray:
        return self.profiles()["radial"]

    def slide_profile(self) -> np.ndarray:
        return self.profiles()["slide"]

    # ── the embedding ───────────────────────────────────────────────────────
    def _frames(self, D: np.ndarray, A: np.ndarray
                ) -> Tuple[np.ndarray, np.ndarray]:
        s, c = np.sin(D), np.cos(D)
        ca, sa = np.cos(A), np.sin(A)
        u = ca[..., None] * self._e1[None, None, :] \
            + sa[..., None] * self._e2[None, None, :]
        n_hat = c[..., None] * self.source_direction[None, None, :] + s[..., None] * u
        e_d = -s[..., None] * self.source_direction[None, None, :] + c[..., None] * u
        return n_hat, e_d

    def point_at(self, d: float, psi: float) -> np.ndarray:
        """The undeformed unit point at geodesic distance ``d``, azimuth ``psi``."""
        return (math.cos(d) * self.source_direction
                + math.sin(d) * (math.cos(psi) * self._e1
                                 + math.sin(psi) * self._e2))

    def positions(self, D: np.ndarray, A: np.ndarray,
                  gain: Optional[float] = None) -> np.ndarray:
        """``X = (R + ερ) n̂ + ε ψ ê_d`` — the continuous deformed surface."""
        eps = self.gain if gain is None else float(gain)
        p = self.profiles()
        rho = np.interp(D, p["d"], p["radial"])
        slide = np.interp(D, p["d"], p["slide"])
        n_hat, e_d = self._frames(D, A)
        r = self.shells.r_mid + eps * rho
        return r[..., None] * n_hat + eps * slide[..., None] * e_d

    def mesh(self, gain: Optional[float] = None
             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        X = self.positions(self._D, self._A, gain=gain)
        return X[..., 0], X[..., 1], X[..., 2]

    def shear_on_mesh(self) -> np.ndarray:
        """The trace-free amplitude at each mesh point — colour, not shape."""
        p = self.profiles()
        return np.interp(self._D, p["d"], p["shear"])

    def material_latitudes(self, n_rings: int = 15, n_points: int = 240,
                           gain: Optional[float] = None) -> list:
        """Rings of material points, evenly spaced *before* the wave arrives.

        Drawn deformed, their spacing is the tidal stretch and squeeze — the
        continuous version of what the ellipses sampled.
        """
        out = []
        for d0 in np.linspace(0.0, math.pi, n_rings + 2)[1:-1]:
            a = np.linspace(0.0, 2.0 * math.pi, n_points)
            D = np.full_like(a, d0)
            out.append(self.positions(D[None, :], a[None, :], gain=gain)[0])
        return out

    def material_meridians(self, n_lines: int = 12, n_points: int = 300,
                           gain: Optional[float] = None) -> list:
        out = []
        for a0 in np.linspace(0.0, 2.0 * math.pi, n_lines, endpoint=False):
            d = np.linspace(0.0, math.pi, n_points)
            A = np.full_like(d, a0)
            out.append(self.positions(d[None, :], A[None, :], gain=gain)[0])
        return out

    # ── arbitrary points on the surface ─────────────────────────────────────
    def coordinates_of(self, unit_points: np.ndarray
                       ) -> Tuple[np.ndarray, np.ndarray]:
        """``(d, ψ)`` of unit vectors, in the source's geodesic-polar frame."""
        p = np.asarray(unit_points, dtype=float)
        d = np.arccos(np.clip(p @ self.source_direction, -1.0, 1.0))
        psi = np.arctan2(p @ self._e2, p @ self._e1)
        return d, psi

    def map_points(self, unit_points: np.ndarray,
                   gain: Optional[float] = None) -> np.ndarray:
        """Carry reference points of the round sphere to the deformed surface."""
        d, psi = self.coordinates_of(unit_points)
        return self.positions(d[None, :], psi[None, :], gain=gain)[0]

    # ── principal axes ──────────────────────────────────────────────────────
    def principal_axes(self, d, h_cross: float = 0.0) -> Dict[str, np.ndarray]:
        """Eigenvalues and eigen-directions of ``h_ab`` in the ``(ê_d, ê_ψ)`` dyad.

        For the trace-free matrix ``[[h₊, h_ˣ], [h_ˣ, −h₊]]`` the eigenvalues
        are ``±√(h₊² + h_ˣ²)`` — **equal in magnitude, always**, which is what
        trace-free means and why the two bars of a principal cross are the
        same length.  The stretch axis is at ``β = ½ atan2(h_ˣ, h₊)``, so it
        swaps between ``ê_d`` and ``ê_ψ`` when ``h₊`` changes sign rather than
        staying put — the axes are computed, never assumed.
        """
        p = self.profiles()
        hp = np.interp(np.asarray(d, dtype=float), p["d"], p["shear"])
        hx = np.full_like(hp, float(h_cross))
        lam = np.hypot(hp, hx)
        beta = 0.5 * np.arctan2(hx, hp)
        return {"lambda_plus": lam, "lambda_minus": -lam, "angle": beta,
                "h_plus": hp, "h_cross": hx}

    def tangent_basis(self, d: float, psi: float,
                      gain: Optional[float] = None,
                      step: float = 1e-4) -> Tuple[np.ndarray, np.ndarray]:
        """An orthonormal tangent basis **of the deformed surface** at a point.

        To first order in ``ε`` this is the ``(ê_d, ê_ψ)`` dyad the tensor is
        written in; taking it from the drawn surface keeps the bars lying in
        the surface at any gain.
        """
        eps = self.gain if gain is None else float(gain)

        def X(dd, aa):
            return self.positions(np.array([[dd]]), np.array([[aa]]),
                                  gain=eps)[0, 0]

        x_d = (X(d + step, psi) - X(d - step, psi)) / (2.0 * step)
        x_a = (X(d, psi + step) - X(d, psi - step)) / (2.0 * step)
        e1 = x_d / max(float(np.linalg.norm(x_d)), 1e-30)
        e2 = x_a - float(x_a @ e1) * e1
        e2 = e2 / max(float(np.linalg.norm(e2)), 1e-30)
        return e1, e2

    def principal_cross(self, d: float, psi: float, size: float = 0.08,
                        reference: Optional[float] = None,
                        gain: Optional[float] = None,
                        lift: float = 0.0) -> Dict[str, object]:
        """Two short tangent bars along the eigen-directions of ``h_ab``.

        The stretch bar (positive eigenvalue) and the squeeze bar (negative
        one), each centred on the surface point, each with half-length
        proportional to ``|λ|``.  ``reference`` is the ``|λ|`` that maps to
        ``size``; it defaults to the current peak on the sphere.
        """
        axes = self.principal_axes(np.array([d]))
        lam = float(axes["lambda_plus"][0])
        beta = float(axes["angle"][0])
        ref = reference if reference is not None else float(
            np.max(np.abs(self.profiles()["shear"])))
        scale = size * (lam / ref) if ref > 0.0 else 0.0
        e1, e2 = self.tangent_basis(d, psi, gain=gain)
        v_plus = math.cos(beta) * e1 + math.sin(beta) * e2
        v_minus = -math.sin(beta) * e1 + math.cos(beta) * e2
        centre = self.positions(np.array([[d]]), np.array([[psi]]),
                                gain=gain)[0, 0]
        if lift:
            normal = np.cross(e1, e2)
            n = float(np.linalg.norm(normal))
            if n > 0.0:
                normal = normal / n
                if float(normal @ centre) < 0.0:
                    normal = -normal          # outward
                centre = centre + lift * normal
        return {
            "centre": centre,
            "stretch": np.stack([centre - scale * v_plus,
                                 centre + scale * v_plus]),
            "squeeze": np.stack([centre - scale * v_minus,
                                 centre + scale * v_minus]),
            "lambda_plus": lam,
            "lambda_minus": -lam,
            "angle": beta,
            "half_length": scale,
            "stretch_direction": v_plus,
            "squeeze_direction": v_minus,
        }

    # ── what the deformation actually is ────────────────────────────────────
    def induced_metric_perturbation(self, d0: float, a0: float = 0.7,
                                    gain: Optional[float] = None,
                                    step: float = 1e-5) -> Dict[str, float]:
        """Measure ``δg_ab`` of the drawn surface, by finite differences.

        The point of the module in one number: the trace-free part comes back
        as the solved ``h``, and the trace comes back as zero.
        """
        eps = self.gain if gain is None else float(gain)
        R = self.shells.r_mid

        def X(d, a):
            return self.positions(np.array([[d]]), np.array([[a]]),
                                  gain=eps)[0, 0]

        X_d = (X(d0 + step, a0) - X(d0 - step, a0)) / (2.0 * step)
        X_a = (X(d0, a0 + step) - X(d0, a0 - step)) / (2.0 * step)
        g_dd = float(X_d @ X_d)
        g_aa = float(X_a @ X_a)
        g_da = float(X_d @ X_a)
        sin2 = math.sin(d0) ** 2
        dg_dd = (g_dd - R ** 2) / (eps * R ** 2)
        dg_aa = (g_aa - R ** 2 * sin2) / (eps * R ** 2 * sin2)
        p = self.profiles()
        return {
            "distance": d0,
            "trace": dg_dd + dg_aa,
            "h_plus_measured": 0.5 * (dg_dd - dg_aa),
            "h_plus_solved": float(np.interp(d0, p["d"], p["shear"])),
            "off_diagonal": g_da / (eps * R ** 2),
        }

    def surface_area(self, n_d: int = 400, n_a: int = 8,
                     gain: Optional[float] = None) -> float:
        """Area of the drawn surface — unchanged at first order, by tracelessness."""
        eps = self.gain if gain is None else float(gain)
        d = np.linspace(1e-4, math.pi - 1e-4, n_d)
        a = np.linspace(0.0, 2.0 * math.pi, n_a, endpoint=False)
        D, A = np.meshgrid(d, a, indexing="ij")
        step = 1e-5
        Xd = (self.positions(D + step, A, gain=eps)
              - self.positions(D - step, A, gain=eps)) / (2.0 * step)
        Xa = (self.positions(D, A + step, gain=eps)
              - self.positions(D, A - step, gain=eps)) / (2.0 * step)
        gdd = np.einsum("...i,...i->...", Xd, Xd)
        gaa = np.einsum("...i,...i->...", Xa, Xa)
        gda = np.einsum("...i,...i->...", Xd, Xa)
        root = np.sqrt(np.clip(gdd * gaa - gda ** 2, 0.0, None))
        return float(np.mean(root, axis=1).sum() * (d[1] - d[0]) * 2.0 * math.pi)

    def excursion(self) -> Dict[str, float]:
        """How far the surface reaches toward each doll, as a fraction."""
        disp = self.gain * self.radial_profile()
        return {
            "outward_fraction": float(np.max(disp)) / self.shells.delta,
            "inward_fraction": float(-np.min(disp)) / self.shells.delta,
            "peak_distance": float(self.sim.d[int(np.argmax(np.abs(disp)))]),
        }


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_induced_metric(surface: Optional[EmbeddedTidalSurface] = None,
                           t: float = 1.2,
                           distances=(0.4, 0.9, 1.5, 2.2, 2.8),
                           gain: float = 1e-4) -> Dict[str, object]:
    """The theorem: the drawn surface's own metric perturbation **is** ``h_ab``.

    Measured at a small gain, where the first-order statement is the whole
    statement; the residuals are the finite-difference error and ``O(ε)``.
    """
    s = surface or EmbeddedTidalSurface()
    s.reset()
    s.advance_to(t)
    raw = [s.induced_metric_perturbation(d0, gain=gain) for d0 in distances]
    # normalise against the largest amplitude on the sphere, not the local one:
    # where the wave has not arrived yet both numbers are ~0 and a local ratio
    # would divide zero by zero
    scale = max(max(abs(r["h_plus_solved"]) for r in raw),
                float(np.max(np.abs(s.profiles()["shear"]))), 1e-12)
    rows = [{**r,
             "relative_error": abs(r["h_plus_measured"]
                                   - r["h_plus_solved"]) / scale,
             "relative_trace": abs(r["trace"]) / scale} for r in raw]
    return {
        "time": s.t,
        "gain": gain,
        "rows": rows,
        "worst_relative_error": max(r["relative_error"] for r in rows),
        "worst_relative_trace": max(r["relative_trace"] for r in rows),
        "worst_off_diagonal": max(abs(r["off_diagonal"]) for r in rows),
        "reproduces_h": bool(max(r["relative_error"] for r in rows) < 5e-3),
    }


def measure_quadrupole_shape(n: int = 3000) -> Dict[str, object]:
    """``ℓ = 2`` gives ``ρ = P₂(cos d)`` — the textbook prolate–oblate shape.

    Derived, not drawn: the quadrupole mode is fed in as the exact ``h`` and
    the radial profile comes out as the Legendre polynomial.
    """
    sim = Spin2WaveSim(n=n)
    sim.start_from(np.ones_like(sim.d))          # q = 1 is exactly ℓ = 2
    s = EmbeddedTidalSurface(sim=sim, gain=1.0)
    p = s.profiles()
    d = p["d"]
    rho = p["radial"]
    legendre = 0.5 * (3.0 * np.cos(d) ** 2 - 1.0)
    scale = float(np.max(np.abs(rho)))
    ratio = scale / float(np.max(np.abs(legendre)))
    err = float(np.max(np.abs(rho - ratio * legendre))) / max(scale, 1e-30)
    w = np.sin(d) * sim.dd
    return {
        "amplitude_ratio": ratio,
        "shape_error": err,
        "is_legendre_p2": bool(err < 1e-3),
        "dipole_constant": float(p["dipole_constant"]),
        "residual_dipole": float(np.sum(rho * np.cos(d) * w)),
        "residual_monopole": float(np.sum(rho * w)),
    }


def measure_area_and_multipoles(surface: Optional[EmbeddedTidalSurface] = None,
                                t: float = 1.2, gain: float = 1e-3
                                ) -> Dict[str, object]:
    """No monopole, no dipole, and the area is unchanged at first order."""
    s = surface or EmbeddedTidalSurface()
    s.reset()
    s.advance_to(t)
    p = s.profiles()
    d, rho = p["d"], p["radial"]
    w = np.sin(d) * s.sim.dd
    area = s.surface_area(gain=gain)
    round_area = 4.0 * math.pi * s.shells.r_mid ** 2
    scale = float(np.max(np.abs(rho))) or 1.0
    return {
        "time": s.t,
        "gain": gain,
        "monopole": float(np.sum(rho * w)) / scale,
        "dipole": float(np.sum(rho * np.cos(d) * w)) / scale,
        "area": area,
        "round_area": round_area,
        "relative_area_change": abs(area - round_area) / round_area,
        "area_is_first_order_unchanged": bool(
            abs(area - round_area) / round_area < gain ** 2 * 50.0),
    }


def measure_bulk_reach(surface: Optional[EmbeddedTidalSurface] = None,
                       frames: int = 200) -> Dict[str, object]:
    """How far into the bulk the tensor wave actually reaches, and where."""
    s = surface or EmbeddedTidalSurface()
    s.reset()
    best = {"out": 0.0, "in": 0.0, "t_out": 0.0, "t_in": 0.0, "where": 0.0}
    for i in range(frames):
        s.advance_to((i + 1) * RETURN_TIME / frames)
        ex = s.excursion()
        if ex["outward_fraction"] > best["out"]:
            best["out"], best["t_out"] = ex["outward_fraction"], s.t
        if ex["inward_fraction"] > best["in"]:
            best["in"], best["t_in"] = ex["inward_fraction"], s.t
            best["where"] = ex["peak_distance"]
    return {
        "max_outward_fraction": best["out"],
        "at_time_outward": best["t_out"],
        "max_inward_fraction": best["in"],
        "at_time_inward": best["t_in"],
        "peak_distance": best["where"],
        "r_inner": s.shells.r_inner,
        "r_outer": s.shells.r_outer,
        "stays_between_the_dolls": bool(best["out"] < 1.0 and best["in"] < 1.0),
    }


# ════════════════════════════════════════════════════════════════════════════
# A MATERIAL PATCH
# ════════════════════════════════════════════════════════════════════════════
class MaterialPatch:
    """A disc of material points, carried by the displacement.

    The patch is a geodesic disc of the *reference* sphere — a set of labelled
    particles — drawn wherever the deformation puts them.  Because the
    perturbation is trace-free, the patch **changes shape without changing
    size**: that is the spin-2 identity made visible rather than numerical, and
    it is checked here by measuring the drawn patch's own area.

    Its shape is summarised by the eigenvalues of the area-weighted second
    moment about its centroid, so the aspect ratio and the long axis are
    measured from the drawn patch rather than read off the tensor.
    """

    def __init__(self, surface: EmbeddedTidalSurface, centre_distance: float,
                 centre_azimuth: float = 0.0, radius: float = 0.20,
                 n_rings: int = 6, n_spokes: int = 48) -> None:
        if not 0.0 < radius < 0.5 * math.pi:
            raise ValueError("patch radius must lie in (0, π/2)")
        self.surface = surface
        self.radius = float(radius)
        self.centre_distance = float(centre_distance)
        self.centre_azimuth = float(centre_azimuth)

        c = surface.point_at(self.centre_distance, self.centre_azimuth)
        tmp = np.array([0.0, 0.0, 1.0])
        if abs(float(tmp @ c)) > 0.9:
            tmp = np.array([1.0, 0.0, 0.0])
        a = tmp - float(tmp @ c) * c
        a = a / np.linalg.norm(a)
        b = np.cross(c, a)
        self.centre_direction = c

        alpha = np.linspace(0.0, 2.0 * math.pi, n_spokes, endpoint=False)
        rings = np.linspace(0.0, self.radius, n_rings + 1)[1:]
        pts = [c]
        for r in rings:
            ring = (math.cos(r) * c[None, :]
                    + math.sin(r) * (np.cos(alpha)[:, None] * a[None, :]
                                     + np.sin(alpha)[:, None] * b[None, :]))
            pts.append(ring)
        self.reference = np.vstack([np.atleast_2d(p) for p in pts])
        self.n_rings = int(n_rings)
        self.n_spokes = int(n_spokes)
        # the boundary is the last ring, closed
        self.boundary_reference = np.vstack([self.reference[-n_spokes:],
                                             self.reference[-n_spokes:][:1]])
        self.reference_area = 2.0 * math.pi * (1.0 - math.cos(self.radius)) \
            * surface.shells.r_mid ** 2

    # ── drawn ───────────────────────────────────────────────────────────────
    def points(self, gain: Optional[float] = None) -> np.ndarray:
        return self.surface.map_points(self.reference, gain=gain)

    def boundary(self, gain: Optional[float] = None) -> np.ndarray:
        return self.surface.map_points(self.boundary_reference, gain=gain)

    def triangles(self, gain: Optional[float] = None) -> np.ndarray:
        """Triangles tiling the drawn patch, for area and for filling it in."""
        p = self.points(gain=gain)
        ns = self.n_spokes
        tris = []
        centre = p[0]
        first = p[1:1 + ns]
        for k in range(ns):
            tris.append([centre, first[k], first[(k + 1) % ns]])
        for r in range(self.n_rings - 1):
            inner = p[1 + r * ns: 1 + (r + 1) * ns]
            outer = p[1 + (r + 1) * ns: 1 + (r + 2) * ns]
            for k in range(ns):
                k2 = (k + 1) % ns
                tris.append([inner[k], outer[k], outer[k2]])
                tris.append([inner[k], outer[k2], inner[k2]])
        return np.asarray(tris)

    def area(self, gain: Optional[float] = None) -> float:
        """Area of the drawn patch — the number that should not move."""
        t = self.triangles(gain=gain)
        cross = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
        return 0.5 * float(np.sum(np.linalg.norm(cross, axis=1)))

    def shape(self, gain: Optional[float] = None) -> Dict[str, object]:
        """Aspect ratio and long axis, from the drawn patch's second moment.

        The shape is read in the patch's **own** plane, which is the least
        variance eigenvector of its full 3×3 second moment.  That matters: at
        a display gain the surface tilts away from radial, so projecting
        against the radial direction instead would leave part of the tilt in
        the answer and rotate the measured long axis out of the surface.
        """
        p = self.points(gain=gain)
        centroid = p.mean(axis=0)
        rel = p - centroid
        cov = rel.T @ rel / len(rel)
        vals, vecs = np.linalg.eigh(cov)                # ascending
        normal = vecs[:, 0]                             # thinnest direction
        long_axis, short_axis = vecs[:, 2], vecs[:, 1]
        ratio = math.sqrt(max(vals[2], 0.0) / max(vals[1], 1e-30))
        return {"centroid": centroid, "aspect_ratio": ratio,
                "long_axis": long_axis, "short_axis": short_axis,
                "normal": normal, "flatness": float(vals[0] / max(vals[1], 1e-30)),
                "moments": vals[:0:-1]}


def measure_patch_area_invariance(surface: Optional[EmbeddedTidalSurface] = None,
                                  centre_distance: float = 1.2,
                                  gain: float = 1e-2, frames: int = 120,
                                  t_end: float = RETURN_TIME
                                  ) -> Dict[str, object]:
    """The patch changes shape without changing size, all the way through.

    Measured on the drawn patch itself.  The residual is second order in the
    gain, which is exactly what trace-free buys.
    """
    s = surface or EmbeddedTidalSurface()
    s.reset()
    patch = MaterialPatch(s, centre_distance=centre_distance)
    areas, ratios, shown = [], [], []
    for i in range(frames):
        s.advance_to((i + 1) * t_end / frames)
        areas.append(patch.area(gain=gain))
        ratios.append(patch.shape(gain=gain)["aspect_ratio"])
        # the same patch at the *display* gain, where the distortion is visible
        shown.append(patch.shape(gain=None)["aspect_ratio"])
    a0 = patch.reference_area
    swing = (max(areas) - min(areas)) / a0
    return {
        "reference_area": a0,
        "min_area": min(areas),
        "max_area": max(areas),
        "relative_area_swing": swing,
        "max_aspect_ratio": max(ratios),
        "max_aspect_ratio_at_display_gain": max(shown),
        "display_gain": s.gain,
        "gain": gain,
        "area_is_invariant": bool(swing < 50.0 * gain ** 2 + 1e-6),
        "shape_does_change": bool(max(shown) > 1.05),
    }


def measure_patch_shape_history(surface: Optional[EmbeddedTidalSurface] = None,
                                centre_distance: float = 1.2,
                                gain: Optional[float] = None, frames: int = 200,
                                radius: float = 0.12,
                                t_end: float = RETURN_TIME) -> Dict[str, object]:
    """When the patch is most distorted, and along which axis.

    The patch's long axis is read from its own second moment, then compared
    with the tensor's stretch eigenvector at the same place — so the picture
    and the algebra are checked against each other rather than assumed equal.

    Two things are reported honestly rather than hidden.  ``flatness`` is the
    patch's out-of-plane moment as a fraction of its short in-plane one: the
    drawn patch is a curved cap, not a disc, and its shape is only meaningful
    in its own plane.  And the alignment degrades with patch *size* near the
    focus, where the shear scale is the pulse width — a patch that straddles
    the focal ring averages over a sign change, so it cannot report a single
    eigenvector.  That is a statement about finite patches, not about the
    field.
    """
    s = surface or EmbeddedTidalSurface()
    s.reset()
    patch = MaterialPatch(s, centre_distance=centre_distance, radius=radius)
    best = {"ratio": 0.0, "time": 0.0, "alignment": 0.0, "flatness": 0.0,
            "lambda": 0.0}
    for i in range(frames):
        s.advance_to((i + 1) * t_end / frames)
        sh = patch.shape(gain=gain)
        if sh["aspect_ratio"] > best["ratio"]:
            axes = s.principal_axes(np.array([centre_distance]))
            e1, e2 = s.tangent_basis(centre_distance, patch.centre_azimuth,
                                     gain=gain)
            beta = float(axes["angle"][0])
            stretch = math.cos(beta) * e1 + math.sin(beta) * e2
            best = {"ratio": float(sh["aspect_ratio"]), "time": s.t,
                    "alignment": abs(float(sh["long_axis"] @ stretch)),
                    "flatness": float(sh["flatness"]),
                    "lambda": float(axes["lambda_plus"][0])}
    return {
        "centre_distance": centre_distance,
        "patch_radius": float(radius),
        "max_aspect_ratio": best["ratio"],
        "at_time": best["time"],
        "eigenvalue_there": best["lambda"],
        "out_of_plane_fraction": best["flatness"],
        "long_axis_alignment": best["alignment"],
        "aligns_with_the_stretch_axis": bool(best["alignment"] > 0.98),
    }


def measure_patch_size_convergence(centre_distance: float = 2.97,
                                   radii: Sequence[float] = (0.24, 0.18, 0.12,
                                                             0.08, 0.05),
                                   frames: int = 200) -> Dict[str, object]:
    """The patch reports the point eigenvector only in the limit of a point.

    Near the antipodal focus the shear turns over on the scale of the pulse,
    so a patch wide enough to straddle the focal ring averages a sign change
    and its long axis tilts away from the local stretch axis.  Shrinking it
    recovers the eigenvector — which is the check that the disagreement is
    the patch's size and not the construction.
    """
    surface = EmbeddedTidalSurface()
    rows = []
    for r in radii:
        m = measure_patch_shape_history(surface, centre_distance=centre_distance,
                                        radius=float(r), frames=frames)
        rows.append({"radius": float(r),
                     "aspect_ratio": m["max_aspect_ratio"],
                     "alignment": m["long_axis_alignment"]})
    align = [row["alignment"] for row in rows]
    return {
        "centre_distance": centre_distance,
        "rows": rows,
        "worst_alignment": min(align),
        "best_alignment": max(align),
        "smallest_patch_alignment": align[-1],
        "converges_to_the_eigenvector": bool(align[-1] > 0.99),
        "improves_as_the_patch_shrinks": bool(align[-1] > align[0]),
    }


def measure_where_the_area_law_fails(surface: Optional[EmbeddedTidalSurface] = None,
                                     smooth_distance: float = 1.20,
                                     focal_distance: Optional[float] = None,
                                     radius: float = 0.12, frames: int = 200,
                                     gains: Sequence[float] = (0.5, 1.0, 2.0, 4.0),
                                     t_end: float = RETURN_TIME
                                     ) -> Dict[str, object]:
    """Trace-free preserves area at first order — and the focus is where that runs out.

    Two identical patches are carried through a full return: one at mid-latitude
    where the field is smooth, one on the focal ring.  Each is asked, at its own
    worst moment, how much its drawn area moved.

    The contrast is the point.  At the *display* gain the smooth patch moves its
    area by about 2%, while the focal patch — same size, same gain — moves by
    about 26%.  Nothing is wrong with either: the area law is first order in
    ``ε``, and its residual is second order *times the local gradient of the
    field*.  Away from the focus that gradient is the wavelength; on the focal
    ring it is the pulse width, so the same term is an order of magnitude
    larger there.  At a gain small enough for the linear statement to be the
    whole statement both collapse — see ``measure_patch_area_invariance``,
    which finds ``1.9e-07`` at ``ε = 1e-2``.

    The residual is confirmed to scale as ``ε²`` (the fitted exponent is
    returned), which is what makes this the boundary of the linear description
    rather than a numerical defect.  Physically it is the honest statement of
    where a refocusing wave stops being describable by ``h_ab`` alone.
    """
    s = surface or EmbeddedTidalSurface()
    focal = (math.pi - 0.94 * s.sim.pulse_width if focal_distance is None
             else float(focal_distance))
    smooth = MaterialPatch(s, centre_distance=smooth_distance, radius=radius,
                           n_rings=8, n_spokes=64)
    focused = MaterialPatch(s, centre_distance=focal, radius=radius,
                            n_rings=8, n_spokes=64)

    s.reset()
    worst = {"smooth": [0.0, 0.0], "focal": [0.0, 0.0]}   # [|dA|, aspect]
    for i in range(frames):
        s.advance_to((i + 1) * t_end / frames)
        for key, patch in (("smooth", smooth), ("focal", focused)):
            flat = patch.area(gain=0.0)
            da = abs(patch.area() / flat - 1.0)
            worst[key][0] = max(worst[key][0], da)
            worst[key][1] = max(worst[key][1],
                                float(patch.shape()["aspect_ratio"]))

    # ...and the residual is second order in the gain, on the focal ring
    s.reset()
    s.advance_to(focal_time(s))
    resid = [abs(focused.area(gain=g) / focused.area(gain=0.0) - 1.0)
             for g in gains]
    lg, lr = np.log(np.asarray(gains, dtype=float)), np.log(np.asarray(resid))
    exponent = float(np.polyfit(lg, lr, 1)[0])

    return {
        "smooth_distance": float(smooth_distance),
        "focal_distance": focal,
        "display_gain": s.gain,
        "smooth_worst_area_change": worst["smooth"][0],
        "smooth_max_aspect_ratio": worst["smooth"][1],
        "focal_worst_area_change": worst["focal"][0],
        "focal_max_aspect_ratio": worst["focal"][1],
        "gains": [float(g) for g in gains],
        "focal_area_residuals": [float(r) for r in resid],
        "residual_exponent": exponent,
        "focal_over_smooth_area_change": worst["focal"][0] / max(worst["smooth"][0],
                                                                1e-30),
        "focal_over_smooth_distortion": (worst["focal"][1] - 1.0) / max(
            worst["smooth"][1] - 1.0, 1e-30),
        "the_focus_distorts_harder": bool(
            worst["focal"][1] - 1.0 > 3.0 * (worst["smooth"][1] - 1.0)),
        "the_area_law_fails_first_at_the_focus": bool(
            worst["focal"][0] > 3.0 * worst["smooth"][0]),
        "residual_is_second_order": bool(abs(exponent - 2.0) < 0.35),
    }


def focal_time(surface: EmbeddedTidalSurface) -> float:
    """When the pulse reaches the focal ring: one antipodal transit, less its width."""
    return float(ANTIPODAL_TIME - 0.94 * surface.sim.pulse_width)
