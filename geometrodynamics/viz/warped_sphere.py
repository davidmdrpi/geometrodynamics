"""
One continuous S², warped by the wave it is carrying.
=====================================================

This restores the geometry the archive rendered before the programme drifted
into algebraic probes: a **single closed surface** whose *radius* carries the
field, nested between two fixed shells like Russian dolls.

::

    r(θ, φ, t) = R_mid + Δ · tanh( g · u(θ, φ, t) / u_ref )

    R_inner = R_mid − Δ  ·······  the inner doll
    R_mid   ···················  the surface the wave lives on
    R_outer = R_mid + Δ  ·······  the outer doll

Three things are worth being exact about.

**The surface is closed.**  It carries its own poles and its φ seam matches to
the last bit, so it is one continuous manifold — not a strip, not a capped
sphere with mouths cut out of it.  That is the whole point: a pulse on it
sweeps every point once and *fills its own void*, which is why a wormhole
cannot be made this way and why a ring is needed (``radial_caustic``).

**The warp is the solved field, not a drawn mound.**  ``archive/
geometrodynamics_v39.py`` displaced the radius with a prescribed ``sech²``
envelope plus a hand-grown mound.  Here the displacement is
:class:`~geometrodynamics.viz.throat_wavefront.BareSphereSim` — a real wave
solve — so the mound at the antipode appears because the wave focuses there,
not because a growth function was keyed to the clock.

**The display gain is honest about being a display gain.**  ``tanh`` is
strictly increasing, so it preserves the sign and the ordering of every
amplitude and never invents a feature; what it does not preserve is the ratio
of amplitudes, and it guarantees the surface stays strictly between the two
shells.  Nothing here is backreaction: the field does not feel the geometry it
is being drawn on.  That step is still open.

The default radii are ``R_inner = 0.74`` and ``R_outer = 1.26`` — deliberately
the same vacuole as :mod:`geometrodynamics.viz.radial_caustic`, so the shell
the ring caustic lands on is the same inner doll drawn here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np

from geometrodynamics.viz.throat_wavefront import BareSphereSim

__all__ = [
    "NestedShells",
    "WarpedSphere",
    "measure_containment",
    "measure_focus",
]

TWO_PI = 2.0 * math.pi
ANTIPODAL_TIME = math.pi          # a point's front reaches the antipode at t = π
RETURN_TIME = TWO_PI              # ...and is home again at t = 2π


# ════════════════════════════════════════════════════════════════════════════
# THE DOLLS
# ════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class NestedShells:
    """Three concentric spheres: the warped one, and the two it lives between.

    ``delta`` is the half-gap, so the middle surface can travel exactly as far
    outward as inward before touching a doll.
    """

    r_mid: float = 1.0
    delta: float = 0.26

    def __post_init__(self) -> None:
        if self.r_mid <= 0.0:
            raise ValueError("r_mid must be positive")
        if not 0.0 < self.delta < self.r_mid:
            raise ValueError("delta must lie in (0, r_mid) — the inner doll "
                             "has to have a positive radius")

    @property
    def r_inner(self) -> float:
        return self.r_mid - self.delta

    @property
    def r_outer(self) -> float:
        return self.r_mid + self.delta

    def contains(self, r) -> bool:
        """Is every radius strictly between the two dolls?"""
        a = np.asarray(r, dtype=float)
        return bool(np.all(a > self.r_inner) and np.all(a < self.r_outer))

    def unit_sphere(self, n_theta: int = 60, n_phi: int = 90
                    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """A closed unit mesh, poles included and seam repeated."""
        th = np.linspace(0.0, math.pi, n_theta)
        ph = np.linspace(0.0, TWO_PI, n_phi)
        TH, PH = np.meshgrid(th, ph, indexing="ij")
        return (np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH))


# ════════════════════════════════════════════════════════════════════════════
# THE WARPED SURFACE
# ════════════════════════════════════════════════════════════════════════════
class WarpedSphere:
    """A closed S² whose radius is displaced by a solved wave.

    Parameters
    ----------
    shells
        The three radii.  Defaults to the ``radial_caustic`` vacuole.
    n_theta, n_phi
        Render resolution of the *surface mesh*.  It is independent of the
        solver: the field is a function of geodesic distance alone, so any
        mesh — including one carrying the poles — is sampled exactly.
    gain
        Display gain inside the ``tanh``.  Larger makes weak late structure
        visible; it cannot change a sign or an ordering.
    """

    def __init__(
        self,
        shells: Optional[NestedShells] = None,
        n_theta: int = 121,
        n_phi: int = 181,
        source_theta: float = 0.5 * math.pi,
        source_phi: float = 0.0,
        pulse_width: float = 0.18,
        gain: float = 1.6,
        n_radial: int = 900,
    ) -> None:
        if n_theta < 8 or n_phi < 8:
            raise ValueError("n_theta and n_phi must be at least 8")
        if gain <= 0.0:
            raise ValueError("gain must be positive")
        self.shells = shells or NestedShells()
        self.gain = float(gain)
        self.n_theta = int(n_theta)
        self.n_phi = int(n_phi)

        self.sim = BareSphereSim(
            n_theta=self.n_theta, n_phi=self.n_phi,
            source_theta=source_theta, source_phi=source_phi,
            pulse_width=pulse_width, n_radial=n_radial)

        # render mesh: poles included, seam repeated so the surface closes
        self.theta = np.linspace(0.0, math.pi, self.n_theta)
        self.phi = np.linspace(0.0, TWO_PI, self.n_phi)
        self._TH, self._PH = np.meshgrid(self.theta, self.phi, indexing="ij")
        self._dist = self.sim.geodesic_distance(self._TH, self._PH)

        # source frame, for drawing circles of constant geodesic distance
        s = np.array([math.sin(source_theta) * math.cos(source_phi),
                      math.sin(source_theta) * math.sin(source_phi),
                      math.cos(source_theta)])
        self.source_direction = s
        tmp = np.array([0.0, 0.0, 1.0])
        if abs(float(np.dot(tmp, s))) > 0.9:
            tmp = np.array([1.0, 0.0, 0.0])
        e1 = tmp - float(np.dot(tmp, s)) * s
        self._e1 = e1 / np.linalg.norm(e1)
        self._e2 = np.cross(s, self._e1)

        self.reference_amplitude = 1.0
        self.calibrate()

    # ── the clock ───────────────────────────────────────────────────────────
    @property
    def t(self) -> float:
        return self.sim.t

    def reset(self) -> None:
        self.sim.reset()

    def advance_to(self, t_target: float) -> None:
        self.sim.advance_to(t_target)

    def calibrate(self, t_end: float = RETURN_TIME, samples: int = 240) -> float:
        """Fix ``u_ref`` from the run's own peak amplitude, then rewind.

        Done on the 1-D solve, so it costs one extra pass over a vector and
        makes the display scale a property of the physics rather than a knob.
        """
        self.sim.reset()
        peak = 0.0
        for i in range(samples):
            self.sim.advance_to((i + 1) * t_end / samples)
            peak = max(peak, float(np.max(np.abs(self.sim._sim.u))))
        self.sim.reset()
        self.reference_amplitude = peak if peak > 0.0 else 1.0
        return self.reference_amplitude

    # ── the field ───────────────────────────────────────────────────────────
    def field(self) -> np.ndarray:
        """The solved field on the render mesh."""
        return self.sim.field_at_distance(self._dist)

    def normalised_field(self) -> np.ndarray:
        """The field divided by the run's peak — in ``[-1, 1]``, sign intact."""
        return self.field() / self.reference_amplitude

    def displacement(self) -> np.ndarray:
        """Radial displacement, strictly inside ``(−Δ, Δ)``."""
        return self.shells.delta * np.tanh(self.gain * self.normalised_field())

    def radius(self) -> np.ndarray:
        """The warped radius — strictly between the two dolls, always."""
        return self.shells.r_mid + self.displacement()

    def radius_at_distance(self, d) -> np.ndarray:
        """The warped radius at a given geodesic distance from the source."""
        u = self.sim.field_at_distance(d) / self.reference_amplitude
        return self.shells.r_mid + self.shells.delta * np.tanh(self.gain * u)

    def mesh(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(X, Y, Z)`` of the warped surface — closed, poles included."""
        r = self.radius()
        return (r * np.sin(self._TH) * np.cos(self._PH),
                r * np.sin(self._TH) * np.sin(self._PH),
                r * np.cos(self._TH))

    # ── the fronts ──────────────────────────────────────────────────────────
    def front_distances(self) -> Dict[str, float]:
        """Where the leading crest and the trailing trough currently are.

        Measured from the 1-D solve — the distance of its extrema — rather
        than assumed to be ``t``, so a front that has passed the antipode and
        is contracting reports the contracting distance.
        """
        s = self.sim._sim
        u = s.u
        i_hi = int(np.argmax(u))
        i_lo = int(np.argmin(u))
        return {
            "crest_distance": float(s.theta[i_hi]),
            "crest_amplitude": float(u[i_hi]),
            "trough_distance": float(s.theta[i_lo]),
            "trough_amplitude": float(u[i_lo]),
        }

    def geodesic_circle(self, d: float, n: int = 240) -> np.ndarray:
        """The circle at geodesic distance ``d``, drawn *on* the warped surface.

        The field depends on ``d`` alone, so the whole circle sits at one
        radius — the ring is a genuine curve of the surface, not a curve
        floating near it.
        """
        psi = np.linspace(0.0, TWO_PI, n)
        dirs = (math.cos(d) * self.source_direction[None, :]
                + math.sin(d) * (np.cos(psi)[:, None] * self._e1[None, :]
                                 + np.sin(psi)[:, None] * self._e2[None, :]))
        return float(self.radius_at_distance(d)) * dirs

    def marker(self, d: float) -> np.ndarray:
        """A point of the warped surface at geodesic distance ``d`` (on the seam)."""
        p = (math.cos(d) * self.source_direction + math.sin(d) * self._e1)
        return float(self.radius_at_distance(d)) * p

    # ── diagnostics ─────────────────────────────────────────────────────────
    def excursion(self) -> Dict[str, float]:
        """How far the surface has travelled toward each doll, as a fraction."""
        disp = self.displacement()
        return {
            "outward_fraction": float(np.max(disp)) / self.shells.delta,
            "inward_fraction": float(-np.min(disp)) / self.shells.delta,
            "peak_distance": float(self._dist.flat[int(np.argmax(np.abs(disp)))]),
        }

    def is_closed(self, tol: float = 1e-12) -> bool:
        """Does the mesh actually close — seam matched, poles single-valued?"""
        X, Y, Z = self.mesh()
        seam = (abs(float(np.max(np.abs(X[:, 0] - X[:, -1])))) < tol
                and abs(float(np.max(np.abs(Y[:, 0] - Y[:, -1])))) < tol
                and abs(float(np.max(np.abs(Z[:, 0] - Z[:, -1])))) < tol)
        poles = (float(np.ptp(Z[0, :])) < tol and float(np.ptp(Z[-1, :])) < tol)
        return bool(seam and poles)


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_containment(surface: WarpedSphere, t_end: float = RETURN_TIME,
                        frames: int = 120) -> Dict[str, object]:
    """The surface never touches either doll, over a whole run."""
    surface.reset()
    r_lo, r_hi = math.inf, -math.inf
    for i in range(frames):
        surface.advance_to((i + 1) * t_end / frames)
        r = surface.radius()
        r_lo = min(r_lo, float(np.min(r)))
        r_hi = max(r_hi, float(np.max(r)))
    sh = surface.shells
    return {
        "r_min": r_lo,
        "r_max": r_hi,
        "r_inner": sh.r_inner,
        "r_outer": sh.r_outer,
        "contained": bool(r_lo > sh.r_inner and r_hi < sh.r_outer),
        "closest_approach_inner": r_lo - sh.r_inner,
        "closest_approach_outer": sh.r_outer - r_hi,
    }


def measure_focus(surface: WarpedSphere, t_end: float = 1.15 * ANTIPODAL_TIME,
                  frames: int = 160) -> Dict[str, object]:
    """When and where the surface deforms most — the antipodal focus.

    The deformation is largest when the front arrives at the antipode, and it
    is *there* that it arrives, both measured rather than imposed.
    """
    surface.reset()
    best = {"amplitude": -math.inf, "time": 0.0, "distance": 0.0}
    for i in range(frames):
        t = (i + 1) * t_end / frames
        surface.advance_to(t)
        disp = np.abs(surface.displacement())
        a = float(np.max(disp))
        if a > best["amplitude"]:
            best = {
                "amplitude": a,
                "time": surface.t,
                "distance": float(surface._dist.flat[int(np.argmax(disp))]),
            }
    return {
        "peak_time": best["time"],
        "peak_distance": best["distance"],
        "peak_displacement": best["amplitude"],
        "antipodal_time": ANTIPODAL_TIME,
        "time_error": abs(best["time"] - ANTIPODAL_TIME),
        "distance_error": abs(best["distance"] - math.pi),
    }
