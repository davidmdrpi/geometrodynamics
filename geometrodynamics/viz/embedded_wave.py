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
from typing import Dict, Optional, Tuple

import numpy as np

from geometrodynamics.viz.spin2_tidal import (
    ANTIPODAL_TIME,
    RETURN_TIME,
    Spin2WaveSim,
)
from geometrodynamics.viz.warped_sphere import NestedShells

__all__ = [
    "EmbeddedTidalSurface",
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
