"""
Ring wavefronts on a closed surface with a throat — the watchable model.

This module restores the *visual* thread of the programme: instead of
probing the algebra of the throat, it **runs the geometry and watches what
the wave does**.  Everything here is a live classical wave solve on a real
Riemannian surface; nothing is fitted, and no quantum ingredient is
inserted anywhere.

The surface
───────────
A unit 2-sphere with both polar caps ``θ < a`` and ``θ > π − a`` removed,
and the two mouth circles joined by a **catenoidal neck** — a surface of
revolution ``r(s) = b·cosh((s − L/2)/b)`` in arclength gauge.  Requiring a
``C¹`` join (the circumference *and* its slope continuous at both mouths)
fixes the neck completely from the single parameter ``a``:

```
r(0) = sin a,  r'(0) = −cos a   ⟹   b = sin a / √(1 + cos²a),
                                    L = 2b·asinh(cos a)
```

The neck has constant Gauss curvature ``K = −1/b²`` (in arclength gauge
``r'' = r/b²``), and the two pieces cancel exactly in Gauss–Bonnet:

```
∫_sphere K dA = 4π cos a,   ∫_neck K dA = −4π cos a   ⟹   χ = 0
```

so the glued surface is a **torus** (azimuth-preserving gluing) or a
**Klein bottle** (azimuth-reversing gluing) — the orientable and
non-orientable throat, selected by :attr:`ThroatWaveSim.twist_parity`.

Why the picture is not the sphere's picture
───────────────────────────────────────────
On the bare sphere a pulse launched from a point expands as a single ring,
reaches the great circle, contracts, and refocuses at the antipode: it
never crosses itself.  Open the throat and the geometry acquires a second,
*shorter* closed path through the neck,

```
ℓ_throat = 2(π/2 − a) + L        vs.       ℓ_sphere = 2π
```

so a piece of the ring dives into the north mouth, crosses the bulk as an
embedded shell, and re-emerges from the south mouth **while the exterior
ring is still running**.  The two fronts then cross: the wavefront
self-intersects on the surface, which is impossible without the handle.
The re-emergent front also arrives back at its *own source* before the
antipodal refocus — a precursor that exists only because the handle exists.

The asymmetries are real
────────────────────────
The sphere is maximally symmetric intrinsically, but the glued surface is
not: the inner (neck) route and the outer (sphere) route between the same
two mouth circles have different lengths, different curvature sign, and —
when ``twist_parity = −1`` — opposite induced orientations.  Gluing the
tube frame ``(∂_s, ∂_ψ)`` to the sphere frame ``(∂_θ, ∂_φ)`` gives
``ds∧dψ = −dθ∧dφ`` at the north mouth and ``−ε·dθ∧dφ`` at the south, so
``ε = +1`` is the orientable handle and ``ε = −1`` reverses orientation on
the throat loop.  A twist offset ``τ = twist_steps · dφ`` rotates where the
re-emergent front lands, steering the refocus without changing any energy.

What is in here
───────────────
* :class:`ThroatGeometry` — the closed-form neck, its curvature and area.
* :class:`ThroatWaveSim` — the coupled solve.  ``mode="throat"`` opens the
  neck; ``mode="plugged"`` seals both mouths with a perfect mirror and is
  the *controlled* baseline: identical grid, identical domain, only the
  handle removed.
* diagnostics: :meth:`~ThroatWaveSim.energy`, :func:`wavefront_components`
  (does the front cross itself?), :func:`measure_transmission`,
  :func:`measure_arrivals`.
* renderers: :func:`draw_hemispheres` (near/far orthographic discs — the
  antipode is *visible*), :func:`draw_neck_strip`, :func:`draw_map`,
  :func:`plot_wavefront_panel`, :func:`run_animation`.
* :func:`export_frames` — quantised frame stacks for external players.

Honest scope
────────────
A linear scalar wave on a *fixed* classical surface: geometry → field, with
no backreaction, so a focus can be sharp but cannot nucleate a new throat
here.  The join is ``C¹`` and not ``C²``, so the mouth carries a curvature
ring which scatters a little on its own; :func:`measure_transmission`
reports the mouth reflection so that effect is visible rather than hidden.
The 2-surface is the *section* of the ``S³`` picture (point → ring → great
circle → antipode is the section of point → shell → maximal shell →
antipode), not the ``S³`` itself.
"""

from __future__ import annotations

import base64
import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure


# ── Palette (matches antipodal_focusing / geometry_panels) ──────────────────
_PAL: Dict[str, str] = dict(
    bg="#0a0e17",
    panel="#111827",
    border="#1e293b",
    text="#e2e8f0",
    dim="#64748b",
    sphere="#fbbf24",     # the exterior surface
    neck="#f472b6",       # the throat / bulk route
    caustic="#fb7185",    # the refocus flash
    plug="#60a5fa",       # the sealed control
    hole="#050810",       # the removed caps
)

_FIELD_CMAP = LinearSegmentedColormap.from_list(
    "throat_field",
    ["#2b6cb0", "#1a3050", "#0a0e17", "#3a2c12", "#d98a2b", "#fbbf24"],
)

C: float = 1.0                  # wave speed (units R = c = 1)
SPHERE_LOOP: float = 2.0 * math.pi   # the bare sphere's closed geodesic
ANTIPODAL_TIME: float = math.pi      # bare-sphere antipodal refocus, t = πR

_CFL: float = 0.40
_MODES: Tuple[str, ...] = ("throat", "plugged")


# ════════════════════════════════════════════════════════════════════════════
# GEOMETRY — the catenoidal neck fixed by the mouth angle alone
# ════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class ThroatGeometry:
    """The neck joining the two mouth circles of a capped unit sphere.

    The mouths sit at polar angle ``a`` and ``π − a``; each is a circle of
    circumference ``2π sin a``.  Continuity of the circumference and of its
    arclength derivative across both mouths determines the neck uniquely.

    Parameters
    ----------
    mouth_angle
        Geodesic radius ``a`` of each removed polar cap, in radians.  Must
        satisfy ``0 < a < π/2``.
    """

    mouth_angle: float

    def __post_init__(self) -> None:
        a = self.mouth_angle
        if not (0.0 < a < 0.5 * math.pi):
            raise ValueError("mouth_angle must lie in (0, π/2)")

    # ── the two closed-form constants ───────────────────────────────────────
    @property
    def neck_radius(self) -> float:
        """Minimum radius ``b`` of the neck (the waist)."""
        a = self.mouth_angle
        return math.sin(a) / math.sqrt(1.0 + math.cos(a) ** 2)

    @property
    def length(self) -> float:
        """Arclength ``L`` of the neck from mouth to mouth."""
        return 2.0 * self.neck_radius * math.asinh(math.cos(self.mouth_angle))

    @property
    def mouth_radius(self) -> float:
        """Circumference radius at either mouth, ``sin a``."""
        return math.sin(self.mouth_angle)

    # ── profile ─────────────────────────────────────────────────────────────
    def radius(self, s: np.ndarray | float) -> np.ndarray | float:
        """Circumference radius ``r(s) = b cosh((s − L/2)/b)``."""
        b = self.neck_radius
        return b * np.cosh((np.asarray(s, dtype=float) - 0.5 * self.length) / b)

    def radius_slope(self, s: np.ndarray | float) -> np.ndarray | float:
        """``r'(s) = sinh((s − L/2)/b)`` — the flare of the neck."""
        b = self.neck_radius
        return np.sinh((np.asarray(s, dtype=float) - 0.5 * self.length) / b)

    def height(self, s: np.ndarray) -> np.ndarray:
        """Axial coordinate ``z(s)`` of the isometric surface of revolution.

        ``z' = √(1 − r'²)``, real because ``|r'| ≤ cos a ≤ 1`` — the neck is
        always embeddable as a surface of revolution in ``ℝ³``.
        """
        s = np.asarray(s, dtype=float)
        dz = np.sqrt(np.clip(1.0 - self.radius_slope(s) ** 2, 0.0, None))
        return np.concatenate([[0.0], np.cumsum(0.5 * (dz[1:] + dz[:-1]) * np.diff(s))])

    # ── curvature / area (closed form) ──────────────────────────────────────
    @property
    def gauss_curvature(self) -> float:
        """Constant Gauss curvature ``K = −1/b²`` of the neck."""
        return -1.0 / self.neck_radius ** 2

    @property
    def area(self) -> float:
        """Neck area ``4π b² cos a``."""
        return 4.0 * math.pi * self.neck_radius ** 2 * math.cos(self.mouth_angle)

    @property
    def sphere_area(self) -> float:
        """Area of the capped sphere, ``4π cos a``."""
        return 4.0 * math.pi * math.cos(self.mouth_angle)

    # ── the two routes between the mouths ───────────────────────────────────
    @property
    def outer_route(self) -> float:
        """Mouth-to-mouth distance *over* the sphere, ``π − 2a``."""
        return math.pi - 2.0 * self.mouth_angle

    @property
    def inner_route(self) -> float:
        """Mouth-to-mouth distance *through* the neck, ``L``."""
        return self.length

    @property
    def shortcut_ratio(self) -> float:
        """How much shorter the bulk route is: ``(π − 2a) / L``."""
        return self.outer_route / self.inner_route

    def free_flight(self, source_theta: float = 0.5 * math.pi) -> float:
        """Time before the expanding ring first reaches a mouth.

        Until then the front is a single unbroken circle on a closed
        surface — it cannot cross itself, whatever the topology.
        """
        return min(abs(source_theta - self.mouth_angle),
                   abs(math.pi - self.mouth_angle - source_theta))

    def mirror_echo(self, source_theta: float = 0.5 * math.pi) -> float:
        """Return time of the *reflected* pulse when the mouth is plugged.

        Straight up the meridian to the mouth and straight back:
        ``2(θ₀ − a)``.  No handle involved.
        """
        return 2.0 * abs(source_theta - self.mouth_angle)

    def throat_loop(self, source_theta: float = 0.5 * math.pi) -> float:
        """Length of the closed geodesic through the neck from a source.

        From a source at polar angle ``θ₀`` on the sphere: up the meridian
        to the north mouth, through the neck, back down to ``θ₀``.  This
        loop generates ``π₁`` and does not exist without the handle.  It
        exceeds :meth:`mirror_echo` by **exactly** the neck length ``L``,
        which is what makes the delay a direct read-out of the bulk.
        """
        return self.mirror_echo(source_theta) + self.length

    def total_curvature(self) -> float:
        """``∫K dA`` over the whole glued surface — exactly ``0`` (χ = 0)."""
        return self.sphere_area * 1.0 + self.gauss_curvature * self.area

    def euler_characteristic(self) -> float:
        """``χ = ∫K dA / 2π`` — ``0`` for both the torus and the Klein bottle."""
        return self.total_curvature() / (2.0 * math.pi)


def orientation_parity(twist_parity: int) -> int:
    """Orientation holonomy around the throat loop: ``+1`` keeps, ``−1`` flips.

    Gluing the neck frame ``(∂_s, ∂_ψ)`` to the sphere frame
    ``(∂_θ, ∂_φ)`` gives ``ds∧dψ = −dθ∧dφ`` at the north mouth (there
    ``∂_s = −∂_θ``, ``∂_ψ = ∂_φ``) and ``ds∧dψ = −ε·dθ∧dφ`` at the south
    (there ``∂_s = −∂_θ`` again, ``∂_ψ = ε ∂_φ``).  A globally consistent
    orientation exists iff the two signs agree, i.e. iff ``ε = +1``.
    """
    if twist_parity not in (1, -1):
        raise ValueError("twist_parity must be +1 or -1")
    return int(twist_parity)


def surface_name(twist_parity: int) -> str:
    """``"torus"`` for the orientable handle, ``"klein_bottle"`` otherwise."""
    return "torus" if orientation_parity(twist_parity) == 1 else "klein_bottle"


def mirror_symmetry_broken(twist_steps: int, n_phi: int) -> bool:
    """Whether the twist offset breaks the meridian mirror of a point source.

    A point source on the sphere is symmetric under the reflection
    ``R: φ → −φ`` through its own meridian.  ``R`` survives the gluing
    ``φ ↦ ε φ + τ`` iff ``R∘g = g∘R``, i.e. ``−(εφ + τ) = ε(−φ) + τ``,
    i.e. ``τ ≡ −τ``, i.e. ``τ ∈ {0, π}``.  For those two offsets the
    ``ε = +1`` and ``ε = −1`` gluings are carried into one another by the
    mirror, so **a mirror-symmetric source cannot see the orientation of
    the throat at all**.  Everywhere else it can.

    This is the precise sense in which the throat's inner/outer asymmetry
    is real but hidden: it takes a twist to expose it.
    """
    tau2 = (2 * int(twist_steps)) % int(n_phi)
    return tau2 != 0


# ════════════════════════════════════════════════════════════════════════════
# THE SOLVE — sphere and neck stepped together
# ════════════════════════════════════════════════════════════════════════════
class ThroatWaveSim:
    """Scalar wave ``u_tt = Δ_g u`` on the capped sphere, with or without a neck.

    Both pieces are cell-centred (no grid point ever lands on a mouth or a
    pole), so the Laplace–Beltrami operator is used in conservative form

    ```
    Δu = (1/r) ∂_x ( r ∂_x u ) + (1/r²) ∂_ψψ u
    ```

    with ``x = θ, r = sin θ`` on the sphere and ``x = s, r = r(s)`` on the
    neck.  The two domains exchange one ghost ring each, linearly
    interpolated onto the other domain's first cell centre, with the mouth
    azimuths identified by ``ψ = ε·φ + τ``.

    Parameters
    ----------
    mode
        ``"throat"`` — the neck is open.  ``"plugged"`` — both mouths are
        perfect mirrors (Neumann); the *controlled* baseline, identical in
        every other respect.
    mouth_angle
        Cap radius ``a`` in radians.
    n_theta, n_phi
        Sphere grid.  ``n_phi`` must be even.
    twist_parity
        ``+1`` azimuth-preserving (orientable handle → torus), ``−1``
        azimuth-reversing (non-orientable → Klein bottle).
    twist_steps
        Azimuthal offset of the gluing in whole ``dφ`` steps, ``τ =
        twist_steps · dφ``.  Rotates where the bulk route re-emerges.
    source_theta, source_phi
        Launch point of the Gaussian cap pulse.
    pulse_width
        Geodesic width of the launch pulse, in radians.
    """

    def __init__(
        self,
        mode: str = "throat",
        mouth_angle: float = 0.5,
        n_theta: int = 96,
        n_phi: int = 128,
        twist_parity: int = 1,
        twist_steps: int = 0,
        source_theta: float = 0.5 * math.pi,
        source_phi: float = 0.0,
        pulse_width: float = 0.14,
        cfl: float = _CFL,
    ) -> None:
        if mode not in _MODES:
            raise ValueError(f"mode must be one of {_MODES}")
        if n_phi % 2:
            raise ValueError("n_phi must be even")
        self.mode = mode
        self.geom = ThroatGeometry(mouth_angle)
        self.twist_parity = orientation_parity(twist_parity)
        self.n_theta = int(n_theta)
        self.n_phi = int(n_phi)
        self.pulse_width = float(pulse_width)
        self.source_theta = float(source_theta)
        self.source_phi = float(source_phi)

        a = mouth_angle
        self.dth = (math.pi - 2.0 * a) / self.n_theta
        self.dphi = 2.0 * math.pi / self.n_phi
        self.twist_steps = int(twist_steps) % self.n_phi
        self.theta = a + (np.arange(self.n_theta) + 0.5) * self.dth
        self.phi = np.arange(self.n_phi) * self.dphi

        # sphere metric: cell radii and face radii
        self._r = np.sin(self.theta)[:, None]
        self._rp = np.sin(self.theta + 0.5 * self.dth)[:, None]   # upper face
        self._rm = np.sin(self.theta - 0.5 * self.dth)[:, None]   # lower face

        # neck grid (arclength), only used when the throat is open
        self.n_s = max(4, int(round(self.geom.length / self.dth)))
        self.ds = self.geom.length / self.n_s
        self.s = (np.arange(self.n_s) + 0.5) * self.ds
        self._q = np.asarray(self.geom.radius(self.s), dtype=float)[:, None]
        self._qp = np.asarray(
            self.geom.radius(self.s + 0.5 * self.ds), dtype=float)[:, None]
        self._qm = np.asarray(
            self.geom.radius(self.s - 0.5 * self.ds), dtype=float)[:, None]

        # ghost interpolation weights (each domain's ghost sits at the other
        # domain's first cell-centre offset, mirrored through the mouth)
        self._w_sphere = (0.5 * self.ds - 0.5 * self.dth) / self.dth
        self._w_neck = (0.5 * self.dth - 0.5 * self.ds) / self.ds

        r_min = min(float(np.min(self._r)), self.geom.neck_radius)
        self.dt = cfl * min(self.dth, self.ds, self.dphi * r_min)

        self.t = 0.0
        self.reset()

    # ── setup ───────────────────────────────────────────────────────────────
    def reset(self) -> None:
        """Launch the Gaussian cap pulse; zero initial velocity."""
        d = self._geodesic_distance_from_source()
        self.u = np.exp(-((d / self.pulse_width) ** 2))
        self.u_prev = self.u.copy()
        self.v = np.zeros((self.n_s, self.n_phi))     # neck field
        self.v_prev = self.v.copy()
        self.t = 0.0
        self._e0 = self.energy()

    def _geodesic_distance_from_source(self) -> np.ndarray:
        th = self.theta[:, None]
        ph = self.phi[None, :]
        cos_d = (
            math.cos(self.source_theta) * np.cos(th)
            + math.sin(self.source_theta) * np.sin(th) * np.cos(ph - self.source_phi)
        )
        return np.arccos(np.clip(cos_d, -1.0, 1.0))

    # ── coupling ────────────────────────────────────────────────────────────
    def _south_index(self) -> np.ndarray:
        """Sphere azimuth index → neck azimuth index at the south mouth."""
        k = np.arange(self.n_phi)
        return (self.twist_parity * k + self.twist_steps) % self.n_phi

    def _south_index_inverse(self) -> np.ndarray:
        """Neck azimuth index → sphere azimuth index at the south mouth."""
        m = np.arange(self.n_phi)
        return (self.twist_parity * (m - self.twist_steps)) % self.n_phi

    def _ghosts(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Ghost rings ``(sphere_north, sphere_south, neck_north, neck_south)``."""
        if self.mode == "plugged":
            return self.u[0], self.u[-1], self.v[0], self.v[-1]
        w_s, w_n = self._w_sphere, self._w_neck
        # sphere's ghost beyond the north mouth ← neck, azimuth identity
        g_sn = (1.0 - w_s) * self.v[0] + w_s * self.v[1]
        # sphere's ghost beyond the south mouth ← neck's far end, twisted
        idx = self._south_index()
        g_ss = ((1.0 - w_s) * self.v[-1] + w_s * self.v[-2])[idx]
        # neck's ghost before s=0 ← sphere just inside the north mouth
        g_nn = (1.0 - w_n) * self.u[0] + w_n * self.u[1]
        # neck's ghost past s=L ← sphere just inside the south mouth, untwisted
        inv = self._south_index_inverse()
        g_ns = ((1.0 - w_n) * self.u[-1] + w_n * self.u[-2])[inv]
        return g_sn, g_ss, g_nn, g_ns

    # ── operators ───────────────────────────────────────────────────────────
    @staticmethod
    def _lap(
        f: np.ndarray, lo: np.ndarray, hi: np.ndarray,
        r: np.ndarray, rp: np.ndarray, rm: np.ndarray,
        dx: float, dphi: float,
    ) -> np.ndarray:
        up = np.vstack([f[1:], hi[None, :]])
        dn = np.vstack([lo[None, :], f[:-1]])
        radial = (rp * (up - f) - rm * (f - dn)) / (r * dx * dx)
        ang = (np.roll(f, -1, axis=1) - 2.0 * f + np.roll(f, 1, axis=1)) / (
            r * r * dphi * dphi)
        return radial + ang

    def laplacian_sphere(self) -> np.ndarray:
        g_sn, g_ss, _, _ = self._ghosts()
        return self._lap(self.u, g_sn, g_ss, self._r, self._rp, self._rm,
                         self.dth, self.dphi)

    def laplacian_neck(self) -> np.ndarray:
        _, _, g_nn, g_ns = self._ghosts()
        return self._lap(self.v, g_nn, g_ns, self._q, self._qp, self._qm,
                         self.ds, self.dphi)

    # ── stepping ────────────────────────────────────────────────────────────
    def step(self) -> None:
        """Advance one leapfrog step of both pieces together."""
        k = (C * self.dt) ** 2
        u_new = 2.0 * self.u - self.u_prev + k * self.laplacian_sphere()
        if self.mode == "throat":
            v_new = 2.0 * self.v - self.v_prev + k * self.laplacian_neck()
            self.v_prev, self.v = self.v, v_new
        self.u_prev, self.u = self.u, u_new
        self.t += self.dt

    def run(self, n_steps: int) -> None:
        for _ in range(int(n_steps)):
            self.step()

    def advance_to(self, t_target: float) -> None:
        while self.t < t_target:
            self.step()

    # ── diagnostics ─────────────────────────────────────────────────────────
    def energy(self) -> float:
        """Total energy ``½∫(u_t² + |∇u|²) dA`` over sphere **and** neck."""
        g_sn, g_ss, g_nn, g_ns = self._ghosts()
        e = self._piece_energy(self.u, self.u_prev, g_sn, g_ss, self._r,
                               self._rp, self.dth)
        if self.mode == "throat":
            e += self._piece_energy(self.v, self.v_prev, g_nn, g_ns, self._q,
                                    self._qp, self.ds)
        return float(e)

    def _piece_energy(
        self, f: np.ndarray, f_prev: np.ndarray, lo: np.ndarray, hi: np.ndarray,
        r: np.ndarray, rp: np.ndarray, dx: float,
    ) -> float:
        f_t = (f - f_prev) / self.dt
        up = np.vstack([f[1:], hi[None, :]])
        d_x = (up - f) / dx                       # on the upper faces
        d_p = (np.roll(f, -1, axis=1) - f) / self.dphi
        kinetic = 0.5 * np.sum(f_t ** 2 * r) * dx * self.dphi
        grad = 0.5 * np.sum((d_x ** 2) * rp) * dx * self.dphi
        grad += 0.5 * np.sum((d_p ** 2) / r) * dx * self.dphi
        return float(kinetic + grad)

    def energy_drift(self) -> float:
        """Relative departure of the total energy from its launch value."""
        return abs(self.energy() - self._e0) / max(self._e0, 1e-30)

    def neck_energy_fraction(self) -> float:
        """Fraction of the total energy currently inside the neck."""
        if self.mode != "throat":
            return 0.0
        _, _, g_nn, g_ns = self._ghosts()
        e_neck = self._piece_energy(self.v, self.v_prev, g_nn, g_ns, self._q,
                                    self._qp, self.ds)
        return float(e_neck / max(self.energy(), 1e-30))

    def energy_density(self) -> np.ndarray:
        """Local energy density ``u_t² + |∇u|²`` on the sphere piece.

        Positive across a whole wavefront (unlike ``|u|``, which vanishes
        between the crest and the trough of one pulse), so it is what the
        front-counting in :func:`wavefront_components` uses.
        """
        _, g_ss, _, _ = self._ghosts()
        u_t = (self.u - self.u_prev) / self.dt
        up = np.vstack([self.u[1:], g_ss[None, :]])
        d_th = (up - self.u) / self.dth
        d_ph = (np.roll(self.u, -1, axis=1) - self.u) / (self.dphi * self._r)
        return u_t ** 2 + d_th ** 2 + d_ph ** 2

    def sample(self, theta: float, phi: float) -> float:
        """Field on the sphere at ``(θ, φ)`` (nearest cell)."""
        j = int(np.clip(round((theta - self.geom.mouth_angle) / self.dth - 0.5),
                        0, self.n_theta - 1))
        k = int(round((phi % (2.0 * math.pi)) / self.dphi)) % self.n_phi
        return float(self.u[j, k])

    def peak(self) -> Tuple[float, float, float]:
        """``(|u|max, θ, φ)`` of the sphere field."""
        j, k = np.unravel_index(int(np.argmax(np.abs(self.u))), self.u.shape)
        return float(abs(self.u[j, k])), float(self.theta[j]), float(self.phi[k])

    def is_finite(self) -> bool:
        return bool(np.all(np.isfinite(self.u)) and np.all(np.isfinite(self.v)))


# ════════════════════════════════════════════════════════════════════════════
# WAVEFRONT TOPOLOGY — does the front cross itself?
# ════════════════════════════════════════════════════════════════════════════
def _label_periodic_phi(mask: np.ndarray) -> int:
    """Count connected components of ``mask`` with ``φ`` (axis 1) periodic."""
    from scipy import ndimage

    lab, n = ndimage.label(mask)
    if n == 0:
        return 0
    parent = list(range(n + 1))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[max(rx, ry)] = min(rx, ry)

    left, right = lab[:, 0], lab[:, -1]
    for lv, rv in zip(left, right):
        if lv and rv:
            union(int(lv), int(rv))
    return len({find(i) for i in range(1, n + 1)})


def wavefront_components(
    sim: ThroatWaveSim, frac: float = 0.15, smooth: float = 2.0,
) -> int:
    """Number of separate wavefront pieces on the sphere.

    The front is located by the **energy density**
    ``u_t² + |∇u|²`` rather than by ``|u|``: the pulse has a crest and a
    trough separated by a zero, so a ``|u|`` level set splits one physical
    ring into concentric shards, while the energy density is positive
    across the whole front.  A light Gaussian smoothing (``smooth`` cells,
    periodic in ``φ``) removes grid speckle.

    A single expanding ring on a closed surface returns ``1`` for its whole
    free flight: it grows to the great circle and shrinks to the antipode
    without ever meeting itself.  Counts above one mean the surface is
    carrying more than one front, and on a closed surface those must cross.
    """
    from scipy import ndimage

    e = sim.energy_density()
    if smooth > 0.0:
        e = ndimage.gaussian_filter(e, sigma=(smooth, smooth),
                                    mode=("nearest", "wrap"))
    peak = float(e.max())
    if peak <= 0.0:
        return 0
    return _label_periodic_phi(e > frac * peak)


def track_wavefront(
    sim: ThroatWaveSim, t_end: float, n_samples: int = 120, frac: float = 0.15,
    t_min: Optional[float] = None,
) -> Dict[str, object]:
    """Component count of the wavefront sampled over a run.

    Returns the sampled series plus the time at which a second front first
    appears — the moment the surface stops carrying a single ring.

    Samples before ``t_min`` (default three launch-pulse widths) are
    recorded but excluded from the split time.  A cap released from rest is
    d'Alembert data: it splits into an outgoing *and* an ingoing ring, and
    the ingoing one collapses through the source and re-expands within
    about one pulse width.  That is a launch artefact of releasing from
    rest, not a wavefront crossing, so the free-flight claim is scored
    after it has resolved.
    """
    t0 = 3.0 * sim.pulse_width if t_min is None else float(t_min)
    sim.reset()
    times: List[float] = []
    counts: List[int] = []
    t_split: Optional[float] = None
    for i in range(n_samples):
        sim.advance_to((i + 1) * t_end / n_samples)
        n = wavefront_components(sim, frac=frac)
        times.append(sim.t)
        counts.append(n)
        if t_split is None and n > 1 and sim.t >= t0:
            t_split = sim.t
    arr = np.asarray(counts)
    scored = arr[np.asarray(times) >= t0]
    return {
        "t_min": t0,
        "max_components": int(scored.max()) if scored.size else 0,
        "t_split": t_split,
        "single_ring_until": float(t_split) if t_split is not None else float(t_end),
        "times": times,
        "counts": counts,
    }


# ════════════════════════════════════════════════════════════════════════════
# ARRIVALS AND TRANSMISSION
# ════════════════════════════════════════════════════════════════════════════
@dataclass
class Arrival:
    """A local maximum of ``|u|`` at a watched point."""

    time: float
    amplitude: float


def watch_point(
    sim: ThroatWaveSim, theta: float, phi: float, t_end: float,
    n_samples: int = 900,
) -> Tuple[np.ndarray, np.ndarray]:
    """Run the sim and return ``(times, u(θ,φ,t))`` at a watched point."""
    sim.reset()
    ts = np.empty(n_samples)
    a = np.empty(n_samples)
    for i in range(n_samples):
        sim.advance_to((i + 1) * t_end / n_samples)
        ts[i] = sim.t
        a[i] = sim.sample(theta, phi)
    return ts, a


def peak_in_window(
    times: np.ndarray, amps: np.ndarray, t_lo: float, t_hi: float,
) -> Arrival:
    """Largest ``|amps|`` inside ``(t_lo, t_hi)``, as an :class:`Arrival`.

    All arrival times measured this way carry the same leading-edge bias of
    roughly the launch pulse's half width — the peak of a finite pulse is
    not its geodesic front.  The bias is common to every watched point, so
    it cancels in *differences* of arrival times, which is how the
    load-bearing measurements below are taken.

    The peak is refined to sub-sample resolution by fitting a parabola to
    the three samples straddling it, so the answer is not quantised by the
    sampling interval.
    """
    m = (times > t_lo) & (times < t_hi)
    if not np.any(m):
        return Arrival(time=float("nan"), amplitude=0.0)
    t = times[m]
    w = np.abs(amps)[m]
    i = int(np.argmax(w))
    if 0 < i < len(w) - 1:
        y0, y1, y2 = w[i - 1], w[i], w[i + 1]
        den = y0 - 2.0 * y1 + y2
        if den < 0.0:                     # a genuine interior maximum
            shift = 0.5 * (y0 - y2) / den
            if abs(shift) <= 1.0:
                dt = float(t[i + 1] - t[i - 1]) * 0.5
                return Arrival(time=float(t[i] + shift * dt),
                               amplitude=float(y1 - 0.25 * (y0 - y2) * shift))
    return Arrival(time=float(t[i]), amplitude=float(w[i]))


def measure_arrivals(
    sim: ThroatWaveSim,
    theta: float,
    phi: float,
    t_end: float,
    n_samples: int = 900,
    t_min: float = 0.0,
    min_amp_frac: float = 0.12,
) -> List[Arrival]:
    """Times at which a wavefront reaches the watched point.

    Peaks are local maxima of ``|u(θ,φ,t)|`` after ``t_min`` and above
    ``min_amp_frac`` of the largest excursion in that window.  Set
    ``t_min`` past the launch transient when watching the source itself,
    or the initial cap swamps the threshold.
    """
    ts, a = watch_point(sim, theta, phi, t_end, n_samples)
    w = np.abs(a)
    keep = ts >= t_min
    if not np.any(keep) or w[keep].max() <= 0.0:
        return []
    thresh = min_amp_frac * w[keep].max()
    out: List[Arrival] = []
    for i in range(1, len(w) - 1):
        if ts[i] < t_min:
            continue
        if w[i] >= w[i - 1] and w[i] > w[i + 1] and w[i] > thresh:
            out.append(Arrival(time=float(ts[i]), amplitude=float(w[i])))
    return out


def measure_echo_delay(
    mouth_angle: float = 0.5,
    n_theta: int = 144,
    n_phi: int = 192,
    pulse_width: float = 0.18,
    t_end: float = 3.6,
    n_samples: int = 900,
    window: float = 0.45,
) -> Dict[str, float]:
    """Read the neck length off the return echo at the source.

    Two runs on the *same* grid from the *same* source, differing only in
    whether the throat is open:

    * plugged — the pulse runs to the mouth, mirrors, and comes back at
      ``2(θ₀ − a)``;
    * throat — the pulse enters the mouth, crosses the bulk, leaves the
      other mouth and comes back at ``2(θ₀ − a) + L``.

    The two paths share every segment except the crossing, so the delay
    between the echoes is the neck arclength ``L`` and the common pulse
    bias cancels.  Also reports how far the open throat suppresses the
    mirror echo — energy that goes *through* the hole instead of bouncing.
    """
    geom = ThroatGeometry(mouth_angle)
    th0 = 0.5 * math.pi
    kw = dict(mouth_angle=mouth_angle, n_theta=n_theta, n_phi=n_phi,
              pulse_width=pulse_width, source_theta=th0, source_phi=0.0)
    t_mir = geom.mirror_echo(th0)
    t_loop = geom.throat_loop(th0)

    ts_p, a_p = watch_point(ThroatWaveSim(mode="plugged", **kw), th0, 0.0,
                            t_end, n_samples)
    ts_t, a_t = watch_point(ThroatWaveSim(mode="throat", **kw), th0, 0.0,
                            t_end, n_samples)

    mirror = peak_in_window(ts_p, a_p, t_mir - window, t_mir + window)
    bulk = peak_in_window(ts_t, a_t, t_loop - window, t_loop + window)
    # the leak window must stop short of the bulk return, or its rising
    # flank is scored as a reflection that never happened
    half = 0.4 * geom.length
    leak = peak_in_window(ts_t, a_t, t_mir - half, t_mir + half)

    delay = bulk.time - mirror.time
    return {
        "neck_length": float(geom.length),
        "t_mirror_predicted": float(t_mir),
        "t_mirror_measured": float(mirror.time),
        "t_bulk_predicted": float(t_loop),
        "t_bulk_measured": float(bulk.time),
        "delay_measured": float(delay),
        "delay_error": float(delay - geom.length),
        "delay_rel_error": float(abs(delay - geom.length) / geom.length),
        "mirror_amplitude": float(mirror.amplitude),
        "bulk_amplitude": float(bulk.amplitude),
        "mirror_suppression": float(1.0 - leak.amplitude / max(mirror.amplitude, 1e-30)),
    }


def measure_bulk_precursor(
    mouth_angle: float = 0.5,
    n_theta: int = 144,
    n_phi: int = 192,
    pulse_width: float = 0.18,
    t_end: float = 3.6,
    n_samples: int = 900,
) -> Dict[str, float]:
    """Does the throat beat the geodesic to the antipode, and can a twist aim it?

    With the gluing offset ``τ = π`` the bulk route runs source → north
    mouth → neck → south mouth → **antipode**, of length
    ``2(θ₀ − a) + L``, which for a small mouth is shorter than the
    geodesic ``π``.  With ``τ = 0`` the same route returns to the source
    instead and the antipode sees nothing extra.  The plugged run is the
    common baseline, so the reported precursor is the throat's own signal.
    """
    geom = ThroatGeometry(mouth_angle)
    th0 = 0.5 * math.pi
    kw = dict(mouth_angle=mouth_angle, n_theta=n_theta, n_phi=n_phi,
              pulse_width=pulse_width, source_theta=th0, source_phi=0.0)
    t_route = geom.throat_loop(th0)
    ts, base = watch_point(ThroatWaveSim(mode="plugged", **kw), th0, math.pi,
                           t_end, n_samples)
    _, aimed = watch_point(
        ThroatWaveSim(mode="throat", twist_steps=n_phi // 2, **kw),
        th0, math.pi, t_end, n_samples)
    _, unaimed = watch_point(
        ThroatWaveSim(mode="throat", twist_steps=0, **kw),
        th0, math.pi, t_end, n_samples)

    lo, hi = t_route - 0.45, t_route + 0.25
    p_aimed = peak_in_window(ts, aimed - base, lo, hi)
    p_unaimed = peak_in_window(ts, unaimed - base, lo, hi)
    focus = peak_in_window(ts, base, t_route + 0.25, t_end)
    return {
        "t_route_predicted": float(t_route),
        "t_geodesic_focus": float(ANTIPODAL_TIME),
        "t_precursor_measured": float(p_aimed.time),
        "t_focus_measured": float(focus.time),
        "lead_time": float(focus.time - p_aimed.time),
        "precursor_aimed": float(p_aimed.amplitude),
        "precursor_unaimed": float(p_unaimed.amplitude),
        "aiming_contrast": float(p_aimed.amplitude / max(p_unaimed.amplitude, 1e-30)),
    }


def measure_orientation_visibility(
    mouth_angle: float = 0.5,
    n_theta: int = 96,
    n_phi: int = 128,
    pulse_width: float = 0.18,
    t_probe: float = 3.4,
) -> Dict[str, object]:
    """Compare the torus and Klein-bottle gluings at several twist offsets.

    Returns, per offset, the largest pointwise difference between the two
    fields relative to the field scale, alongside the mirror-symmetry
    prediction of :func:`mirror_symmetry_broken`.  The two must agree: the
    orientation is invisible exactly at ``τ ∈ {0, π}``.
    """
    rows: List[Dict[str, object]] = []
    for steps in (0, n_phi // 8, n_phi // 4, 3 * n_phi // 8, n_phi // 2):
        fields = []
        for parity in (1, -1):
            s = ThroatWaveSim(mode="throat", mouth_angle=mouth_angle,
                              n_theta=n_theta, n_phi=n_phi,
                              twist_parity=parity, twist_steps=steps,
                              pulse_width=pulse_width)
            s.advance_to(t_probe)
            fields.append(s.u.copy())
        scale = max(float(np.max(np.abs(fields[0]))), 1e-30)
        rows.append({
            "twist_steps": int(steps),
            "tau_over_pi": float(2.0 * steps / n_phi),
            "relative_difference": float(np.max(np.abs(fields[0] - fields[1])) / scale),
            "mirror_broken": bool(mirror_symmetry_broken(steps, n_phi)),
        })
    return {"t_probe": float(t_probe), "rows": rows}


def measure_transmission(
    sim: ThroatWaveSim, t_end: float, n_samples: int = 400,
) -> Dict[str, float]:
    """Energy budget of the first encounter with the mouth.

    Returns the peak fraction of the total energy that is simultaneously
    inside the neck (``transmitted``), the complement left outside
    (``reflected``), the time of that peak, and the energy drift over the
    run (the numerical honesty check).
    """
    sim.reset()
    best, t_best = 0.0, 0.0
    for i in range(n_samples):
        sim.advance_to((i + 1) * t_end / n_samples)
        f = sim.neck_energy_fraction()
        if f > best:
            best, t_best = f, sim.t
    return {
        "transmitted": float(best),
        "reflected": float(1.0 - best),
        "t_peak": float(t_best),
        "energy_drift": float(sim.energy_drift()),
    }


def measure_focus(
    sim: ThroatWaveSim,
    theta: float,
    phi: float,
    t_end: float,
    n_samples: int = 600,
) -> Dict[str, float]:
    """Largest excursion at a watched point and when it happens."""
    sim.reset()
    best, t_best = 0.0, 0.0
    for i in range(n_samples):
        sim.advance_to((i + 1) * t_end / n_samples)
        a = abs(sim.sample(theta, phi))
        if a > best:
            best, t_best = a, sim.t
    return {"peak": float(best), "t_peak": float(t_best)}


# ════════════════════════════════════════════════════════════════════════════
# RENDERING
# ════════════════════════════════════════════════════════════════════════════
def _style(ax: Axes, title: str = "") -> None:
    ax.set_facecolor(_PAL["panel"])
    for sp in ax.spines.values():
        sp.set_color(_PAL["border"])
    ax.tick_params(colors=_PAL["dim"], labelsize=7)
    if title:
        ax.set_title(title, color=_PAL["text"], fontsize=9, pad=6)


def _orthographic_grid(n_px: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    lin = np.linspace(-1.0, 1.0, n_px)
    X, Y = np.meshgrid(lin, lin)
    rr = X * X + Y * Y
    inside = rr <= 1.0
    Z = np.sqrt(np.clip(1.0 - rr, 0.0, None))
    return X, Y, Z, inside


def sphere_image(
    sim: ThroatWaveSim, n_px: int = 220, far_side: bool = False,
    view_theta: Optional[float] = None, view_phi: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Orthographic image of the sphere field, and the in-disc mask.

    The view axis defaults to the source direction, so ``far_side=False``
    shows the launch hemisphere and ``far_side=True`` shows the hemisphere
    containing the antipode — the refocus is watchable without rotating.
    Cells inside the removed caps are masked out (they are holes).
    """
    vt = sim.source_theta if view_theta is None else view_theta
    vp = sim.source_phi if view_phi is None else view_phi
    n = np.array([math.sin(vt) * math.cos(vp), math.sin(vt) * math.sin(vp),
                  math.cos(vt)])
    up = np.array([0.0, 0.0, 1.0])
    e1 = np.cross(up, n)
    if np.linalg.norm(e1) < 1e-9:
        e1 = np.array([1.0, 0.0, 0.0])
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(n, e1)

    X, Y, Z, inside = _orthographic_grid(n_px)
    sign = -1.0 if far_side else 1.0
    xs = (sign * X)[..., None] * e1 + Y[..., None] * e2 + (sign * Z)[..., None] * n
    th = np.arccos(np.clip(xs[..., 2], -1.0, 1.0))
    ph = np.arctan2(xs[..., 1], xs[..., 0]) % (2.0 * math.pi)

    j = np.clip(np.round((th - sim.geom.mouth_angle) / sim.dth - 0.5).astype(int),
                0, sim.n_theta - 1)
    k = (np.round(ph / sim.dphi).astype(int)) % sim.n_phi
    img = sim.u[j, k]
    a = sim.geom.mouth_angle
    hole = (th < a) | (th > math.pi - a)
    return np.where(inside & ~hole, img, np.nan), inside & ~hole


def draw_hemispheres(
    ax_near: Axes, ax_far: Axes, sim: ThroatWaveSim, n_px: int = 220,
    vmax: Optional[float] = None,
) -> None:
    """Paint the launch hemisphere and the antipodal hemisphere side by side."""
    v = vmax if vmax is not None else max(float(np.max(np.abs(sim.u))), 1e-9)
    for ax, far, label in ((ax_near, False, "launch side"),
                           (ax_far, True, "antipodal side")):
        img, _ = sphere_image(sim, n_px=n_px, far_side=far)
        ax.clear()
        ax.imshow(img, cmap=_FIELD_CMAP, vmin=-v, vmax=v, origin="lower",
                  extent=(-1, 1, -1, 1), interpolation="bilinear")
        ax.set_xlim(-1.05, 1.05)
        ax.set_ylim(-1.05, 1.05)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        _style(ax, label)


def draw_neck_strip(ax: Axes, sim: ThroatWaveSim, vmax: Optional[float] = None) -> None:
    """The neck field unrolled onto ``(s, ψ)`` — the shell inside the throat."""
    ax.clear()
    if sim.mode != "throat":
        ax.text(0.5, 0.5, "neck sealed (control)", color=_PAL["plug"],
                ha="center", va="center", fontsize=9, transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        _style(ax, "throat interior")
        return
    v = vmax if vmax is not None else max(float(np.max(np.abs(sim.u))), 1e-9)
    ax.imshow(sim.v.T, cmap=_FIELD_CMAP, vmin=-v, vmax=v, origin="lower",
              aspect="auto", interpolation="bilinear",
              extent=(0.0, sim.geom.length, 0.0, 2.0 * math.pi))
    ax.set_xlabel("neck arclength s", color=_PAL["dim"], fontsize=8)
    ax.set_ylabel("azimuth ψ", color=_PAL["dim"], fontsize=8)
    _style(ax, f"throat interior  (L = {sim.geom.length:.3f}, b = "
               f"{sim.geom.neck_radius:.3f})")


def draw_map(ax: Axes, sim: ThroatWaveSim, vmax: Optional[float] = None) -> None:
    """Equirectangular map of the sphere field — the whole surface at once."""
    ax.clear()
    v = vmax if vmax is not None else max(float(np.max(np.abs(sim.u))), 1e-9)
    ax.imshow(sim.u, cmap=_FIELD_CMAP, vmin=-v, vmax=v, origin="lower",
              aspect="auto", interpolation="bilinear",
              extent=(0.0, 2.0 * math.pi, sim.geom.mouth_angle,
                      math.pi - sim.geom.mouth_angle))
    ax.set_xlabel("φ", color=_PAL["dim"], fontsize=8)
    ax.set_ylabel("θ", color=_PAL["dim"], fontsize=8)
    n_front = wavefront_components(sim)
    _style(ax, f"surface map — fronts: {n_front}")


def draw_geometry(ax: Axes, geom: ThroatGeometry) -> None:
    """Profile of the capped sphere and its neck, drawn to scale."""
    ax.clear()
    a = geom.mouth_angle
    th = np.linspace(a, math.pi - a, 200)
    ax.plot(np.sin(th), np.cos(th), color=_PAL["sphere"], lw=2.0,
            label="capped sphere")
    ax.plot(-np.sin(th), np.cos(th), color=_PAL["sphere"], lw=2.0)
    s = np.linspace(0.0, geom.length, 200)
    r = np.asarray(geom.radius(s))
    z = geom.height(s)
    z = z - 0.5 * z[-1]
    off = 1.55
    ax.plot(off + r, z, color=_PAL["neck"], lw=2.0, label="catenoidal neck")
    ax.plot(off - r, z, color=_PAL["neck"], lw=2.0)
    ax.axhline(0.0, color=_PAL["border"], lw=0.6)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    leg = ax.legend(loc="lower center", fontsize=7, frameon=False)
    for t in leg.get_texts():
        t.set_color(_PAL["dim"])
    _style(ax, f"geometry  (a = {a:.2f}, χ = {geom.euler_characteristic():.2g})")


def plot_wavefront_panel(sim: ThroatWaveSim, figsize: Tuple[float, float] = (11.0, 7.0)) -> Figure:
    """A full still: both hemispheres, the neck interior, the map, the geometry."""
    fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
    gs = fig.add_gridspec(2, 3, hspace=0.32, wspace=0.25)
    draw_hemispheres(fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), sim)
    draw_geometry(fig.add_subplot(gs[0, 2]), sim.geom)
    draw_map(fig.add_subplot(gs[1, :2]), sim)
    draw_neck_strip(fig.add_subplot(gs[1, 2]), sim)
    fig.suptitle(
        f"{sim.mode}  ·  t = {sim.t:.3f}  ·  {surface_name(sim.twist_parity)}",
        color=_PAL["text"], fontsize=11)
    return fig


def run_animation(
    sim: Optional[ThroatWaveSim] = None, frames: int = 300, t_end: float = 6.0,
    interval: int = 40,
):
    """Live plugged-vs-throat animation (requires an interactive backend)."""
    from matplotlib import animation

    throat = sim if sim is not None else ThroatWaveSim(mode="throat")
    plug = ThroatWaveSim(
        mode="plugged", mouth_angle=throat.geom.mouth_angle,
        n_theta=throat.n_theta, n_phi=throat.n_phi,
        source_theta=throat.source_theta, source_phi=throat.source_phi,
        pulse_width=throat.pulse_width)
    throat.reset()
    plug.reset()

    fig = plt.figure(figsize=(12.0, 6.4), facecolor=_PAL["bg"])
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.22)
    axes = [fig.add_subplot(gs[r, c]) for r in (0, 1) for c in (0, 1, 2)]
    dt_frame = t_end / frames

    def update(i: int):
        throat.advance_to((i + 1) * dt_frame)
        plug.advance_to((i + 1) * dt_frame)
        draw_hemispheres(axes[0], axes[1], plug)
        draw_neck_strip(axes[2], plug)
        draw_hemispheres(axes[3], axes[4], throat)
        draw_neck_strip(axes[5], throat)
        axes[0].set_ylabel("plugged", color=_PAL["plug"], fontsize=9)
        axes[3].set_ylabel("throat", color=_PAL["neck"], fontsize=9)
        fig.suptitle(f"t = {throat.t:.3f}", color=_PAL["text"], fontsize=11)
        return axes

    return animation.FuncAnimation(fig, update, frames=frames,
                                   interval=interval, blit=False)


# ════════════════════════════════════════════════════════════════════════════
# EXPORT — quantised frame stacks for external players
# ════════════════════════════════════════════════════════════════════════════
def _downsample(a: np.ndarray, rows: int, cols: int) -> np.ndarray:
    ri = np.linspace(0, a.shape[0] - 1, rows).round().astype(int)
    ci = np.linspace(0, a.shape[1] - 1, cols).round().astype(int)
    return a[np.ix_(ri, ci)]


def export_frames(
    sim: ThroatWaveSim, t_end: float = 6.0, frames: int = 160,
    rows: int = 72, cols: int = 120, neck_rows: int = 20,
    compand: str = "linear",
) -> Dict[str, object]:
    """Run the sim and return a compact, base64 uint8 frame stack.

    The field is symmetrically quantised about zero against a single global
    scale so that every frame shares one colour mapping; ``scale`` and
    ``compand`` are returned so a player can recover physical amplitudes.

    Parameters
    ----------
    compand
        ``"linear"`` — bytes are proportional to the amplitude.  The launch
        pulse then sets the scale and everything the wave does afterwards
        lands in the bottom few percent of the range, which is faithful but
        unwatchable.  ``"signed_sqrt"`` stores
        ``sign(u)·√(|u|/scale)`` instead, which keeps the sign and the
        ordering of amplitudes exactly while lifting the late, weak
        structure into view — the display convention, not a change to the
        physics.
    """
    if compand not in ("linear", "signed_sqrt"):
        raise ValueError("compand must be 'linear' or 'signed_sqrt'")
    sim.reset()
    sph: List[np.ndarray] = []
    nek: List[np.ndarray] = []
    times: List[float] = []
    neck_frac: List[float] = []
    for i in range(frames):
        sim.advance_to((i + 1) * t_end / frames)
        sph.append(_downsample(sim.u, rows, cols))
        nek.append(_downsample(sim.v, neck_rows, cols))
        times.append(sim.t)
        neck_frac.append(sim.neck_energy_fraction())

    stack = np.stack(sph)
    neck = np.stack(nek)
    if compand == "linear":
        scale = float(np.percentile(np.abs(stack), 99.5)) or 1.0
        shape = lambda a: a / scale
    else:
        scale = float(np.max(np.abs(stack))) or 1.0
        shape = lambda a: np.sign(a) * np.sqrt(np.abs(a) / scale)
    q = np.clip(np.round(shape(stack) * 127.0) + 128, 0, 255).astype(np.uint8)
    qn = np.clip(np.round(shape(neck) * 127.0) + 128, 0, 255).astype(np.uint8)
    return {
        "compand": compand,
        "mode": sim.mode,
        "surface": surface_name(sim.twist_parity),
        "mouth_angle": sim.geom.mouth_angle,
        "neck_length": sim.geom.length,
        "neck_radius": sim.geom.neck_radius,
        "twist_parity": sim.twist_parity,
        "twist_steps": sim.twist_steps,
        "rows": rows, "cols": cols, "neck_rows": neck_rows,
        "frames": frames, "t_end": t_end, "scale": scale,
        "times": times,
        "neck_energy_fraction": neck_frac,
        "sphere_b64": base64.b64encode(q.tobytes()).decode("ascii"),
        "neck_b64": base64.b64encode(qn.tobytes()).decode("ascii"),
    }
