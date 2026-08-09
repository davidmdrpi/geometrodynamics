"""
Ring wavefronts on a closed surface with a throat — the watchable model.

This module restores the *visual* thread of the programme: instead of
probing the algebra of the throat, it **runs the geometry and watches what
the wave does**.  Everything here is a live classical wave solve on a real
Riemannian surface; nothing is fitted, and no quantum ingredient is
inserted anywhere.

Three surfaces, one clock
─────────────────────────
The point of the module is the comparison, and it needs all three terms:

1. :class:`BareSphereSim` — the **uncut unit 2-sphere**.  The canonical
   picture: point → expanding ring → great circle → contracting ring →
   antipodal focus, with nothing to reflect from and nothing to fall
   through.
2. ``ThroatWaveSim(mode="plugged")`` — the same sphere with both polar caps
   ``θ < a``, ``θ > π − a`` cut out and the mouths **sealed by a mirror**.
   This isolates what merely cutting holes does.
3. ``ThroatWaveSim(mode="throat")`` — the same cut sphere with the mouths
   **joined by a neck**.  This isolates what opening a second geometric
   route does.

The neck is a genuine **catenoid** — the minimal surface of revolution,
``H = 0``, ``r = b cosh(z/b)``.  Requiring the circumference *and* its
arclength slope to match the sphere at the mouth fixes it completely from
the single parameter ``a``:

```
r = sin a,  dr/ds = −cos a   ⟹   b = sin²a,   L = sin 2a
```

Its Gauss curvature ``K = −b²/r⁴`` **varies**: exactly ``−1`` at each
mouth — the sphere's ``+1`` with its sign flipped, no jump in magnitude —
deepening to ``−1/sin⁴a`` at the waist.  So the wave crosses a real
minimal-surface neck of continuously changing negative curvature, not a
manufactured constant-curvature strip.

On Gauss–Bonnet, honestly
─────────────────────────
For *any* surface of revolution ``K = −r''/r`` and ``dA = 2πr ds``, so
``∫K dA = −2π[r']`` depends **only on the boundary slopes**.  The ``C¹``
join pins those, so ``∫K dA = 4π cos a − 4π cos a = 0`` and ``χ = 0`` for
the catenoid and for any other ``C¹``-matched profile alike.  The closure
is a check on the *join*, never evidence for a particular neck.

The glued surface is a **torus** (azimuth-preserving gluing) or a **Klein
bottle** (azimuth-reversing).  Gluing the neck frame ``(∂_s, ∂_ψ)`` to the
sphere frame ``(∂_θ, ∂_φ)`` gives ``ds∧dψ = −dθ∧dφ`` at the north mouth and
``−ε·dθ∧dφ`` at the south, so ``ε = +1`` is the orientable handle.  A twist
offset ``τ = twist_steps · dφ`` rotates where the bulk route re-emerges,
steering the refocus at no energy cost.

How the two pieces are coupled
──────────────────────────────
Each mouth is **one finite-volume face shared** by a sphere cell and a neck
cell.  Its flux is evaluated once and handed to both with opposite signs,
so the discrete divergence theorem holds across the mouth and the discrete
energy is conserved to round-off (``energy_drift ~ 1e-15``).  The launch is
a purely **outgoing** ring: a cap released from rest is d'Alembert data and
splits into an outgoing *and* an ingoing front, which would put two fronts
on the surface for reasons that have nothing to do with the geometry.  Its
zero mode is then pinned — a closed surface has no boundary, so
``d²/dt² ∫u dA = 0`` and the mean field ramps linearly unless
``∫u_t dA = 0``, and a one-way launch is a monopole as well as a ring.

What is in here
───────────────
* :class:`ThroatGeometry` — the catenoid in closed form, its varying
  curvature, its area, and the two routes between the mouths.
* :class:`ThroatWaveSim`, :class:`BareSphereSim` — the three solves.
* :func:`measure_arrival_multiplicity` — how many times a front crosses
  each point; the honest self-intersection diagnostic.
* :func:`measure_mouth_budget` — transmission by **integrated flux**, not
  by a snapshot of what happens to be stored in the neck.
* :func:`measure_echo_delay`, :func:`measure_bulk_precursor`,
  :func:`measure_orientation_visibility`.
* renderers: :func:`draw_hemispheres` (near/far orthographic discs — the
  antipode is *visible*), :func:`draw_neck_strip`, :func:`draw_map`,
  :func:`plot_wavefront_panel`, :func:`run_animation`.
* :func:`export_frames` — quantised frame stacks for external players.

Honest scope
────────────
A linear scalar wave on a *fixed* classical surface: geometry → field, with
no backreaction, so a focus can be sharp but cannot nucleate a new throat
here.  The join is ``C¹`` and not ``C²``, so each mouth carries a curvature
ring which scatters on its own; that scattering is inside the reported
mouth budget rather than removed from it.  The 2-surface is the *section*
of the ``S³`` picture (point → ring → great circle → antipode is the
section of point → shell → maximal shell → antipode), not the ``S³``
itself.
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
# GEOMETRY — a true catenoid, fixed by the mouth angle alone
# ════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class ThroatGeometry:
    """The catenoid joining the two mouth circles of a capped unit sphere.

    The mouths sit at polar angle ``a`` and ``π − a``; each is a circle of
    circumference ``2π sin a``.  The neck is a genuine **catenoid** — the
    minimal surface of revolution, ``H = 0`` — written ``r = b cosh(z/b)``
    in its axial coordinate ``z``, so ``dr/ds = tanh(z/b)`` in arclength.

    Requiring the circumference *and* its arclength slope to match the
    sphere at the mouth,

    ```
    r = sin a,   dr/ds = −cos a   ⟹   tanh(z₀/b) = cos a,  b cosh(z₀/b) = sin a
    ```

    has the closed-form solution

    ```
    b = sin²a,        L = sin 2a,        z₀ = b·artanh(cos a)
    ```

    with ``L`` the arclength from mouth to mouth.  The Gauss curvature
    ``K = −b²/r⁴`` is **not** constant: it is exactly ``−1`` at each mouth
    — the sphere's ``+1`` flipped in sign, with no jump in magnitude — and
    deepens to ``−1/sin⁴a`` at the waist.

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

    # ── the closed-form constants ───────────────────────────────────────────
    @property
    def neck_radius(self) -> float:
        """Waist radius ``b = sin²a`` — the narrowest circle of the catenoid."""
        return math.sin(self.mouth_angle) ** 2

    @property
    def length(self) -> float:
        """Arclength ``L = sin 2a`` from mouth to mouth."""
        return math.sin(2.0 * self.mouth_angle)

    @property
    def axial_half_height(self) -> float:
        """``z₀ = b·artanh(cos a)`` — half the catenoid's axial extent."""
        return self.neck_radius * math.atanh(math.cos(self.mouth_angle))

    @property
    def mouth_radius(self) -> float:
        """Circumference radius at either mouth, ``sin a``."""
        return math.sin(self.mouth_angle)

    # ── profile, in arclength ───────────────────────────────────────────────
    def _sigma(self, s: np.ndarray | float) -> np.ndarray:
        """``σ = sinh(z/b)``, the natural variable: ``σ = s/b − cot a``."""
        s = np.asarray(s, dtype=float)
        return s / self.neck_radius - 1.0 / math.tan(self.mouth_angle)

    def radius(self, s: np.ndarray | float) -> np.ndarray | float:
        """``r(s) = b√(1+σ²)`` — the catenoid in arclength gauge."""
        return self.neck_radius * np.sqrt(1.0 + self._sigma(s) ** 2)

    def radius_slope(self, s: np.ndarray | float) -> np.ndarray | float:
        """``dr/ds = σ/√(1+σ²) = tanh(z/b)`` — always in ``(−1, 1)``."""
        sig = self._sigma(s)
        return sig / np.sqrt(1.0 + sig ** 2)

    def height(self, s: np.ndarray | float) -> np.ndarray | float:
        """Axial coordinate ``z(s) = b·asinh(σ)``, exact (no quadrature)."""
        return self.neck_radius * np.arcsinh(self._sigma(s))

    # ── curvature / area ────────────────────────────────────────────────────
    def curvature(self, s: np.ndarray | float) -> np.ndarray | float:
        """Gauss curvature ``K(s) = −b²/r⁴`` — varying along the neck."""
        r = np.asarray(self.radius(s), dtype=float)
        return -self.neck_radius ** 2 / r ** 4

    @property
    def curvature_at_mouth(self) -> float:
        """``K = −1`` exactly: the sphere's ``+1`` with its sign flipped."""
        return float(self.curvature(0.0))

    @property
    def curvature_at_waist(self) -> float:
        """``K = −1/b² = −1/sin⁴a`` — the deepest curvature of the neck."""
        return -1.0 / self.neck_radius ** 2

    @property
    def area(self) -> float:
        """Neck area ``π b² (2w + sinh 2w)``, with ``w = z₀/b = artanh(cos a)``."""
        b, w = self.neck_radius, self.axial_half_height / self.neck_radius
        return math.pi * b * b * (2.0 * w + math.sinh(2.0 * w))

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
        """Mouth-to-mouth distance *through* the neck, ``L = sin 2a``."""
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

    # ── Gauss–Bonnet ────────────────────────────────────────────────────────
    def neck_total_curvature(self) -> float:
        """``∫K dA`` over the neck — ``−4π cos a``, whatever the profile.

        For *any* surface of revolution ``K = −r''/r`` and ``dA = 2πr ds``,
        so ``∫K dA = −2π[r']`` depends only on the boundary slopes.  The
        ``C¹`` join pins those to ``∓cos a``, so the neck's total curvature
        is fixed by the join alone.  **The χ = 0 closure is therefore a
        check on the join, not evidence for any particular neck** — the
        catenoid and any other ``C¹``-matched profile close it equally.
        """
        return -4.0 * math.pi * math.cos(self.mouth_angle)

    def total_curvature(self) -> float:
        """``∫K dA`` over the whole glued surface — exactly ``0`` (χ = 0)."""
        return self.sphere_area + self.neck_total_curvature()

    def euler_characteristic(self) -> float:
        """``χ = ∫K dA / 2π`` — ``0`` for both the torus and the Klein bottle."""
        return self.total_curvature() / (2.0 * math.pi)


def _outgoing_velocity(
    f: np.ndarray, d: np.ndarray, w: float, weight: np.ndarray,
) -> np.ndarray:
    """Initial velocity for a purely outgoing ring with no monopole.

    Two conditions have to hold at once, and the obvious launches each
    violate one of them.

    * A cap released from rest (``u_t = 0``) is d'Alembert data: it splits
      into an outgoing *and* an ingoing ring, so the surface carries two
      fronts for reasons that have nothing to do with the geometry.  The
      one-way condition is ``u_t = −∂f/∂d``.
    * A closed surface has no boundary, so ``d²/dt² ∫u dA = ∫Δu dA = 0``
      and the mean field ramps **linearly** unless ``∫u_t dA = 0``.  The
      one-way launch is a monopole as well as a ring and fails this.

    Subtracting the *mean* of ``u_t`` fixes the second at the cost of
    something worse: it gives every point of the surface — including
    everywhere the front has not reached — a non-zero initial velocity, so
    the far side starts moving before anything could have arrived.  The
    correction has to stay inside the pulse.  Adding a multiple of the
    profile itself does that:

    ```
    u_t = −∂f/∂d + c·f,      c = −∫(−∂f/∂d) dA / ∫f dA
    ```

    which is one-way, zero-monopole, and identically zero wherever ``f``
    is, so the rest of the surface is left alone.
    """
    grad = 2.0 * d / (w * w) * f            # = −∂f/∂d for a Gaussian in d
    c = -float(np.sum(grad * weight)) / float(np.sum(f * weight))
    return grad + c * f


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

        # sphere metric: cell radii, and the n+1 face radii (both mouths at sin a)
        self._r = np.sin(self.theta)[:, None]
        self._rf = np.sin(a + np.arange(self.n_theta + 1) * self.dth)[:, None]

        # neck grid (arclength), only used when the throat is open
        self.n_s = max(4, int(round(self.geom.length / self.dth)))
        self.ds = self.geom.length / self.n_s
        self.s = (np.arange(self.n_s) + 0.5) * self.ds
        self._q = np.asarray(self.geom.radius(self.s), dtype=float)[:, None]
        self._qf = np.asarray(
            self.geom.radius(np.arange(self.n_s + 1) * self.ds),
            dtype=float)[:, None]

        r_min = min(float(np.min(self._r)), self.geom.neck_radius)
        self.dt = cfl * min(self.dth, self.ds, self.dphi * r_min)

        self.t = 0.0
        self.reset()

    # ── setup ───────────────────────────────────────────────────────────────
    def reset(self) -> None:
        """Launch a purely **outgoing** Gaussian ring pulse.

        Releasing a cap from rest is d'Alembert data: it splits into an
        outgoing *and* an ingoing ring, and the ingoing one collapses
        through the source and re-expands, so the surface carries two
        fronts forever after.  That is a property of the launch, not of the
        geometry, and it wrecks any honest count of how many times a front
        crosses a point.  Setting the previous time level to the same
        profile displaced *outward* by one step, ``u⁻ = f(d + c·dt)``, makes
        the pulse one-way to leading order and leaves a single ring.
        """
        d = self._geodesic_distance_from_source()
        w = self.pulse_width
        f = np.exp(-((d / w) ** 2))
        self.u = f
        self.u_prev = f - self.dt * _outgoing_velocity(f, d, w, self._r)
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

    def interface_flux(self) -> Tuple[np.ndarray, np.ndarray]:
        """Shared mouth fluxes ``(north, south)``, as the **sphere** sees them.

        Each mouth is one finite-volume face shared by a sphere cell and a
        neck cell.  Its flux is evaluated once,

        ```
        F = r_mouth · (f_outer − f_inner) / h,     h = ½(dθ + ds)
        ```

        and handed to both cells with opposite signs, so whatever leaves one
        domain enters the other exactly — the discrete divergence theorem
        holds across the mouth and energy is conserved to round-off.  (The
        earlier scheme gave each domain its own interpolated ghost ring;
        those two ghosts disagreed at ``O(dθ)`` and leaked.)

        A sealed mouth is the zero-flux face, which is the same statement as
        a perfect mirror and is conservative for the same reason.
        """
        if self.mode != "throat":
            z = np.zeros(self.n_phi)
            return z, z
        R = self.geom.mouth_radius
        h = 0.5 * (self.dth + self.ds)
        north = R * (self.u[0] - self.v[0]) / h
        south = R * (self.v[-1][self._south_index()] - self.u[-1]) / h
        return north, south

    # ── operators ───────────────────────────────────────────────────────────
    @staticmethod
    def _radial(
        f: np.ndarray, f_lo: np.ndarray, f_up: np.ndarray,
        r: np.ndarray, rf: np.ndarray, dx: float,
    ) -> np.ndarray:
        """``(1/r)∂_x(r ∂_x f)`` with the two boundary face fluxes supplied."""
        inner = rf[1:-1] * (f[1:] - f[:-1]) / dx
        flux = np.concatenate([f_lo[None, :], inner, f_up[None, :]], axis=0)
        return (flux[1:] - flux[:-1]) / (r * dx)

    def _angular(self, f: np.ndarray, r: np.ndarray) -> np.ndarray:
        return (np.roll(f, -1, axis=1) - 2.0 * f + np.roll(f, 1, axis=1)) / (
            r * r * self.dphi * self.dphi)

    def laplacian_sphere(self) -> np.ndarray:
        north, south = self.interface_flux()
        return (self._radial(self.u, north, south, self._r, self._rf, self.dth)
                + self._angular(self.u, self._r))

    def laplacian_neck(self) -> np.ndarray:
        north, south = self.interface_flux()
        inv = self._south_index_inverse()
        return (self._radial(self.v, -north, -south[inv], self._q, self._qf,
                             self.ds)
                + self._angular(self.v, self._q))

    # ── stepping ────────────────────────────────────────────────────────────
    def step(self) -> None:
        """Advance one leapfrog step of both pieces together."""
        k = (C * self.dt) ** 2
        lap_u = self.laplacian_sphere()
        lap_v = self.laplacian_neck() if self.mode == "throat" else None
        u_new = 2.0 * self.u - self.u_prev + k * lap_u
        if lap_v is not None:
            v_new = 2.0 * self.v - self.v_prev + k * lap_v
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
    def _piece_energy(
        self, f: np.ndarray, f_prev: np.ndarray, r: np.ndarray,
        rf: np.ndarray, dx: float,
    ) -> float:
        """Cell kinetic + interior-face gradient energy of one domain.

        The gradient term is the **cross** product between the two stored
        time levels, ``⟨∇fⁿ, ∇fⁿ⁻¹⟩``, not ``|∇fⁿ|²``.  That is the exact
        invariant of the leapfrog step: the velocity ``(fⁿ − fⁿ⁻¹)/dt``
        lives at the half step, so pairing it with a same-time gradient
        leaves an ``O(dt²)`` wobble that peaks whenever a sharp front is
        crossing something.  With the cross term the discrete energy is
        conserved to round-off, which makes ``energy_drift`` a real check
        on the scheme instead of a measure of its own staggering.
        """
        f_t = (f - f_prev) / self.dt
        kinetic = 0.5 * np.sum(f_t ** 2 * r) * dx * self.dphi
        d_x = (f[1:] - f[:-1]) / dx
        d_x0 = (f_prev[1:] - f_prev[:-1]) / dx
        grad = 0.5 * np.sum(d_x * d_x0 * rf[1:-1]) * dx * self.dphi
        d_p = (np.roll(f, -1, axis=1) - f) / self.dphi
        d_p0 = (np.roll(f_prev, -1, axis=1) - f_prev) / self.dphi
        grad += 0.5 * np.sum(d_p * d_p0 / r) * dx * self.dphi
        return float(kinetic + grad)

    def _interface_energy(self) -> float:
        """Gradient energy stored on the two shared mouth faces."""
        if self.mode != "throat":
            return 0.0
        R = self.geom.mouth_radius
        h = 0.5 * (self.dth + self.ds)
        idx = self._south_index()
        dn = (self.u[0] - self.v[0]) / h
        dn0 = (self.u_prev[0] - self.v_prev[0]) / h
        ds_ = (self.v[-1][idx] - self.u[-1]) / h
        ds0 = (self.v_prev[-1][idx] - self.u_prev[-1]) / h
        return float(0.5 * R * (np.sum(dn * dn0) + np.sum(ds_ * ds0))
                     * h * self.dphi)

    def energy(self) -> float:
        """Total energy ``½∫(u_t² + |∇u|²) dA`` over sphere **and** neck.

        Includes the gradient energy held on the two shared mouth faces, so
        the sum is the exact discrete invariant of the coupled scheme.
        """
        e = self._piece_energy(self.u, self.u_prev, self._r, self._rf, self.dth)
        if self.mode == "throat":
            e += self._piece_energy(self.v, self.v_prev, self._q, self._qf,
                                    self.ds)
            e += self._interface_energy()
        return float(e)

    def neck_energy(self) -> float:
        """Energy currently stored inside the neck."""
        if self.mode != "throat":
            return 0.0
        return self._piece_energy(self.v, self.v_prev, self._q, self._qf, self.ds)

    def neck_energy_fraction(self) -> float:
        """Fraction of the total energy currently stored inside the neck.

        A *storage* fraction, not a transmission coefficient — the wave is
        still going in and out.  For the throughput see
        :func:`measure_mouth_budget`.
        """
        if self.mode != "throat":
            return 0.0
        return float(self.neck_energy() / max(self.energy(), 1e-30))

    def mouth_power(self) -> Tuple[float, float]:
        """Instantaneous power ``(north, south)`` flowing sphere → neck.

        Read off the same shared face flux the stepper uses.  Abel summation
        on the discrete operator gives, for the neck alone,

        ```
        dE_neck/dt = dφ ( north·v_t[0] − south[inv]·v_t[−1] )
        ```

        so these two terms are exactly the neck's energy budget and they
        integrate to its stored energy — the check that the mouth coupling
        is conservative rather than merely plausible.
        """
        if self.mode != "throat":
            return 0.0, 0.0
        north, south = self.interface_flux()
        inv = self._south_index_inverse()
        v_t = (self.v - self.v_prev) / self.dt
        p_n = float(np.sum(north * v_t[0]) * self.dphi)
        p_s = float(-np.sum(south[inv] * v_t[-1]) * self.dphi)
        return p_n, p_s

    def inward_power_at(self, theta_face: float) -> float:
        """Power crossing the sphere circle ``θ_face`` toward the north mouth.

        Used to normalise the mouth budget: it is the energy actually
        *offered* to the mouth, which on a closed surface is only a part of
        the wave.
        """
        j = int(np.clip(round((theta_face - self.geom.mouth_angle) / self.dth),
                        1, self.n_theta - 1))
        f = self._rf[j, 0] * (self.u[j] - self.u[j - 1]) / self.dth
        u_t = 0.5 * ((self.u[j] - self.u_prev[j])
                     + (self.u[j - 1] - self.u_prev[j - 1])) / self.dt
        return float(np.sum(f * u_t) * self.dphi)

    def energy_density(self) -> np.ndarray:
        """Local energy density ``u_t² + |∇u|²`` on the sphere piece."""
        u_t = (self.u - self.u_prev) / self.dt
        d_th = np.zeros_like(self.u)
        d_th[:-1] = (self.u[1:] - self.u[:-1]) / self.dth
        d_th[-1] = d_th[-2]
        d_ph = (np.roll(self.u, -1, axis=1) - self.u) / (self.dphi * self._r)
        return u_t ** 2 + d_th ** 2 + d_ph ** 2

    def energy_drift(self) -> float:
        """Relative departure of the total energy from its launch value."""
        return abs(self.energy() - self._e0) / max(self._e0, 1e-30)

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
# THE BARE SPHERE — the canonical picture, no holes at all
# ════════════════════════════════════════════════════════════════════════════
class BareSphereSim:
    """The uncut unit 2-sphere, presented on the same ``(θ, φ)`` grid.

    This is the control the other two are *departures from*: point →
    expanding ring → great circle → contracting ring → antipodal focus, with
    nothing cut out and nothing to reflect from.

    A point source on a sphere makes the field a function of the geodesic
    distance from that source alone, so the solve is **exactly** the
    axisymmetric one already in
    :mod:`geometrodynamics.viz.antipodal_focusing` — no new physics and no
    polar-grid pathology.  This class runs that 1-D solver and maps it onto
    the 2-D grid through ``d(θ, φ)``, which is exact, not an interpolation
    of a coarser 2-D run.

    The public surface matches :class:`ThroatWaveSim` closely enough for the
    renderers, the exporter and the diagnostics to take either.
    """

    mode = "bare"
    has_neck = False

    def __init__(
        self,
        n_theta: int = 96,
        n_phi: int = 128,
        source_theta: float = 0.5 * math.pi,
        source_phi: float = 0.0,
        pulse_width: float = 0.18,
        n_radial: int = 900,
        **_ignored,
    ) -> None:
        from geometrodynamics.viz.antipodal_focusing import SphereWaveSim

        self.n_theta = int(n_theta)
        self.n_phi = int(n_phi)
        self.pulse_width = float(pulse_width)
        self.source_theta = float(source_theta)
        self.source_phi = float(source_phi)
        self.geom = None
        self.n_s = 0

        self.dth = math.pi / self.n_theta
        self.dphi = 2.0 * math.pi / self.n_phi
        self.theta = (np.arange(self.n_theta) + 0.5) * self.dth
        self.phi = np.arange(self.n_phi) * self.dphi
        self.theta_min, self.theta_max = 0.0, math.pi

        self._sim = SphereWaveSim(n=n_radial)
        self.dt = self._sim.dt
        th = self.theta[:, None]
        ph = self.phi[None, :]
        cos_d = (math.cos(self.source_theta) * np.cos(th)
                 + math.sin(self.source_theta) * np.sin(th)
                 * np.cos(ph - self.source_phi))
        self._dist = np.arccos(np.clip(cos_d, -1.0, 1.0))
        self.reset()

    # ── setup / stepping ────────────────────────────────────────────────────
    def reset(self) -> None:
        """Launch a purely outgoing ring at the source, at this width."""
        s = self._sim
        w = self.pulse_width
        f = np.exp(-((s.theta / w) ** 2))
        s.u = f
        s.u_prev = f - s.dt * _outgoing_velocity(f, s.theta, w, np.sin(s.theta))
        s.t = 0.0
        self._e0 = self.energy()

    @property
    def t(self) -> float:
        return self._sim.t

    def step(self) -> None:
        self._sim.step()

    def run(self, n_steps: int) -> None:
        self._sim.run(int(n_steps))

    def advance_to(self, t_target: float) -> None:
        self._sim.advance_to(t_target)

    # ── fields ──────────────────────────────────────────────────────────────
    @property
    def u(self) -> np.ndarray:
        return self._sim.sample(self._dist)

    @property
    def u_prev(self) -> np.ndarray:
        return np.interp(np.clip(self._dist, 0.0, math.pi),
                         self._sim.theta, self._sim.u_prev)

    @property
    def v(self) -> np.ndarray:
        """Empty neck field — there is no neck."""
        return np.zeros((0, self.n_phi))

    def sample(self, theta: float, phi: float) -> float:
        cos_d = (math.cos(self.source_theta) * math.cos(theta)
                 + math.sin(self.source_theta) * math.sin(theta)
                 * math.cos(phi - self.source_phi))
        d = math.acos(max(-1.0, min(1.0, cos_d)))
        return float(self._sim.sample(np.array([d]))[0])

    def energy_density(self) -> np.ndarray:
        """``u_t² + |∇u|²``, evaluated in 1-D and mapped — no polar blow-up.

        The field depends only on ``d``, so ``|∇u|² = (∂u/∂d)²`` exactly and
        the ``1/sin θ`` factor that would wreck a polar grid never appears.
        """
        s = self._sim
        u_t = (s.u - s.u_prev) / s.dt
        u_d = np.gradient(s.u, s.dth)
        dens = u_t ** 2 + u_d ** 2
        return np.interp(np.clip(self._dist, 0.0, math.pi), s.theta, dens)

    # ── diagnostics (same names, trivially satisfied) ───────────────────────
    def energy(self) -> float:
        return float(self._sim.energy())

    def energy_drift(self) -> float:
        return abs(self.energy() - self._e0) / max(self._e0, 1e-30)

    def neck_energy_fraction(self) -> float:
        return 0.0

    def is_finite(self) -> bool:
        return bool(np.all(np.isfinite(self._sim.u)))


# ════════════════════════════════════════════════════════════════════════════
# WAVEFRONT TOPOLOGY — does the front cross itself?
# ════════════════════════════════════════════════════════════════════════════
@dataclass
class Multiplicity:
    """Outcome of :func:`measure_arrival_multiplicity`."""

    counts: np.ndarray            # arrivals per grid cell over the window
    max_arrivals: int
    area_fraction_multi: float    # area-weighted fraction with ≥ 2 arrivals
    source_side_fraction: float   # ...of that, within 60° of the source meridian
    first_multi_time: Optional[float]
    first_multi_theta: Optional[float]
    t_window: float


def measure_arrival_multiplicity(
    sim, t_window: float, hi: float = 0.50, lo: float = 0.15,
    t_min: Optional[float] = None,
) -> Multiplicity:
    """How many times a wavefront crosses each point of the surface.

    This replaces counting connected components of a level set, which could
    not tell the two interesting cases apart: a hole *cutting* one ring into
    arcs raises the component count without any front ever meeting another
    front, and that is exactly what both the plugged and the open surface do
    the moment the ring reaches a mouth.

    Arrival multiplicity asks the question directly.  Each cell runs a
    **hysteresis trigger** on its own energy density ``u_t² + |∇u|²``: an
    arrival is counted when the density rises past ``hi × (that cell's peak
    over the window)``, and the cell cannot count another until it has
    fallen back below ``lo ×`` the same peak.  Plain local-maximum counting
    does not survive here — a wave in 2+1 dimensions violates Huygens'
    principle, so every front drags a slowly decaying wake whose ripples are
    local maxima too, and they get counted as arrivals that never happened.

    The thresholds are calibrated on the case whose answer is known: the
    bare sphere, where the front provably passes once.  Any ``hi ≥ 0.5``
    returns *exactly* zero second arrivals there, and the sealed-vs-open
    contrast is stable across that whole range, so the conclusion is not an
    artefact of the pair chosen.  Below it the wake starts tripping the
    trigger and the bare sphere reports fronts that do not exist.

    On a closed surface with no boundary the front is the geodesic circle of
    radius ``t``: it sweeps each point exactly once, so within a half period
    the count is ``1`` everywhere and no two fronts ever meet.  A count of
    ``2`` means two fronts crossed the same point — and on a closed surface
    two fronts that coexist must cross.

    Runs the sim twice: once to learn each cell's peak, once to trigger
    against it.  ``t_min`` (default two pulse widths) skips the launch.
    """
    t0 = 2.0 * sim.pulse_width if t_min is None else float(t_min)
    shape = (sim.n_theta, sim.n_phi)

    sim.reset()
    peaks = np.zeros(shape)
    while sim.t < t_window:
        peaks = np.maximum(peaks, sim.energy_density())
        sim.step()
    peaks = np.maximum(peaks, peaks.max() * 1e-9)

    sim.reset()
    counts = np.zeros(shape, dtype=np.int32)
    armed = np.ones(shape, dtype=bool)          # ready to register an arrival
    first_t: Optional[float] = None
    first_row: Optional[int] = None
    while sim.t < t_window:
        e = sim.energy_density()
        if sim.t >= t0:
            fire = armed & (e > hi * peaks)
            if np.any(fire):
                was = counts.copy()
                counts[fire] += 1
                armed[fire] = False
                new_multi = (counts >= 2) & (was < 2)
                if first_t is None and np.any(new_multi):
                    first_t = sim.t
                    first_row = int(np.argmax(new_multi.any(axis=1)))
            armed |= e < lo * peaks
        else:
            armed &= e < lo * peaks             # let the launch settle
        sim.step()

    w = np.sin(sim.theta)[:, None] * np.ones((1, sim.n_phi))
    multi = counts >= 2
    frac_multi = float(np.sum(w * multi) / np.sum(w))
    # the discriminating half: a reflection sends a second front back toward
    # the source, a bulk crossing sends it out of the far mouth instead
    dphi_src = np.abs((sim.phi - sim.source_phi + math.pi) % (2 * math.pi) - math.pi)
    near = (dphi_src < math.pi / 3.0)[None, :] & np.ones((sim.n_theta, 1), bool)
    frac_src = float(np.sum(w * multi * near) / np.sum(w * near))
    return Multiplicity(
        counts=counts,
        max_arrivals=int(counts.max()),
        area_fraction_multi=frac_multi,
        source_side_fraction=frac_src,
        first_multi_time=first_t,
        first_multi_theta=(float(sim.theta[first_row])
                           if first_row is not None else None),
        t_window=float(t_window),
    )


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


def measure_mouth_budget(
    sim: ThroatWaveSim, t_window: Optional[float] = None, n_steps: int = 4000,
    reference_cells: int = 6,
) -> Dict[str, float]:
    """Energy budget of the first encounter with the mouth, by **flux**.

    The earlier version reported the peak *stored* energy fraction inside
    the neck and called it "transmitted".  That is not a transmission
    coefficient: it is a snapshot of what happens to be in the neck at one
    instant, it depends on how long the neck is, and it ignores everything
    that has already passed through.  This measures the throughput instead,
    by integrating the actual power crossing surfaces:

    * ``offered`` — the energy that crosses a reference circle a few cells
      inside the mouth, heading for it.  On a closed surface only part of
      the wave ever reaches the mouth at all, so this, not the total energy,
      is the right denominator.
    * ``through`` — the energy that crosses the mouth face itself into the
      neck, from the same shared flux the stepper uses.
    * ``transmission = through / offered``, and ``reflection = 1 −
      transmission``.

    Both are running maxima over the first encounter, before the wave can
    come back and cross the same surfaces again.
    """
    if sim.mode != "throat":
        raise ValueError("mouth budget needs mode='throat'")
    geom = sim.geom
    if t_window is None:
        t_window = geom.free_flight(sim.source_theta) + 1.2 * geom.length
    theta_ref = geom.mouth_angle + reference_cells * sim.dth

    sim.reset()
    q_in = q_ref = 0.0
    best_in = best_ref = 0.0
    stored = 0.0
    t_best = 0.0
    while sim.t < t_window:
        p_n, _ = sim.mouth_power()
        q_in += p_n * sim.dt
        q_ref += sim.inward_power_at(theta_ref) * sim.dt
        best_ref = max(best_ref, q_ref)
        if q_in > best_in:
            best_in, t_best = q_in, sim.t
        stored = max(stored, sim.neck_energy_fraction())
        sim.step()

    trans = best_in / best_ref if best_ref > 0 else 0.0
    return {
        "offered": float(best_ref),
        "through": float(best_in),
        "transmission": float(trans),
        "reflection": float(1.0 - trans),
        "peak_stored_fraction": float(stored),
        "t_through_peak": float(t_best),
        "t_window": float(t_window),
        "theta_reference": float(theta_ref),
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
    _style(ax, "surface map")


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
    has_neck = getattr(sim, "n_s", 0) > 0
    sph: List[np.ndarray] = []
    nek: List[np.ndarray] = []
    times: List[float] = []
    neck_frac: List[float] = []
    for i in range(frames):
        sim.advance_to((i + 1) * t_end / frames)
        sph.append(_downsample(sim.u, rows, cols))
        nek.append(_downsample(sim.v, neck_rows, cols) if has_neck
                   else np.zeros((neck_rows, cols)))
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
    g = sim.geom
    return {
        "compand": compand,
        "mode": sim.mode,
        "surface": surface_name(getattr(sim, "twist_parity", 1)),
        "mouth_angle": g.mouth_angle if g else 0.0,
        "neck_length": g.length if g else 0.0,
        "neck_radius": g.neck_radius if g else 0.0,
        "theta_min": float(sim.theta[0] - 0.5 * sim.dth),
        "theta_max": float(sim.theta[-1] + 0.5 * sim.dth),
        "has_neck": bool(has_neck),
        "twist_parity": getattr(sim, "twist_parity", 1),
        "twist_steps": getattr(sim, "twist_steps", 0),
        "rows": rows, "cols": cols, "neck_rows": neck_rows,
        "frames": frames, "t_end": t_end, "scale": scale,
        "times": times,
        "neck_energy_fraction": neck_frac,
        "sphere_b64": base64.b64encode(q.tobytes()).decode("ascii"),
        "neck_b64": base64.b64encode(qn.tobytes()).decode("ascii"),
    }
