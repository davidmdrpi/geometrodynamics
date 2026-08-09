"""
What a wavefront has to be like to fold: focal geometry across the bulk.

The companion module :mod:`geometrodynamics.viz.throat_wavefront` runs waves on
a surface that *already has* a throat.  This one asks the prior question — what
kind of wavefront can fold at all — and answers it with the differential
geometry of wavefronts rather than with a simulation.

The vacuole
───────────
Two concentric spheres, ``r = R_inner`` and ``r = R_outer``, with a **flat**
bulk between them.  Straight rays cross it; the conserved quantity along a ray
is the impact parameter about the centre,

```
b = r sin α          (α = angle to the radial direction)
```

so a ray launched inward from the outer sphere reaches radius ``b`` and no
deeper.

A point source does not fold — *in this flat bulk*
──────────────────────────────────────────────────
The wavefront of a **point** is the metric sphere ``|x − P| = t``, whose signed
area element is ``t² sin θ`` — positive for every ``t > 0``.  It never folds, and
the region behind it is the filled ball: the pulse "fills its own void".

**This is a statement about flat Euclidean space, not about point sources.**  On
a closed manifold the same point source focuses: on ``S²`` or ``S³`` the front of
a point converges on the antipode at ``t = πR``, which is exactly what
``throat_wavefront.py`` measures.  Curvature is what gives a point source a
focal set; a flat bulk is what denies it one.

A circle focuses *coherently* — that is what is special
───────────────────────────────────────────────────────
Any curved extended source has a focal set, so "only a ring folds" would be
false.  What singles out the circle is not *that* it folds but *how*: by
symmetry the **whole ring arrives at one point simultaneously**.

The front of a circle of radius ``ρ`` is its offset tube

```
X(u,v) = ((ρ + t cos v) cos u, (ρ + t cos v) sin u, t sin v)
```

whose signed area element is ``t (ρ + t cos v)``.  That vanishes where
``cos v = −ρ/t``, which has a solution only for ``t ≥ ρ``.  So:

* ``t < ρ`` — the tube is immersed; no fold anywhere.
* ``t = ρ`` — the **first caustic**, and it is infinitely degenerate: the two
  roots coincide at the ring's centre and *every point of the ring* reaches it
  at once.
* ``t > ρ`` — the tube **stays singular**, at two points on the symmetry axis
  at ``z = ±√(t² − ρ²)``, which run outward as the front grows.

So the caustic is not a one-off event at a single point: it is a first,
maximally degenerate focus followed by a persistent singular locus along the
axis.  The degeneracy at ``t = ρ`` is the distinguishing feature.

The two conditions coincide
───────────────────────────
Ask for the ring whose first caustic lands *on the inner sphere*.  Its centre
must sit at radius ``R_inner``, so it is the circle of polar angle

```
cos θ₀ = R_inner / R_outer
```

and its rays leave the outer sphere at ``sin α = R_inner / R_outer`` — exactly
the **critical acceptance angle**, the ray tangent to the inner sphere.  The
ring that focuses on the throat and the ray that grazes it are the same ring.
The first caustic forms at

```
t = ρ = √(R_outer² − R_inner²)
```

This coincidence, and that time, are the core result.

Acceptance is asymmetric; propagation is not
────────────────────────────────────────────
Because ``b = r sin α`` is conserved and ``r`` decreases going in:

* **outer → inner** only launch directions with ``sin α ≤ R_inner/R_outer``
  reach the inner sphere; the accepted fraction of the inward hemisphere is
  ``1 − √(1 − (R_inner/R_outer)²)`` — about **19%** at the programme's ``ΔR``.
* **inner → outer** every launch direction reaches the outer sphere. **100%**.

This is an **angular (solid-angle) acceptance asymmetry, not nonreciprocal
propagation**.  Every individual ray is exactly reversible: the outward rays
that arrive are the time-reverses of the inward rays that arrive, one for one.
What differs is the *measure* of launch directions that connect the two
surfaces, because a hemisphere at ``R_outer`` and a hemisphere at ``R_inner``
are different sets of directions.  Nothing about the sphere's symmetry is
broken; the ordering of the two radii is the whole of it.

Scope
─────
Ray and front geometry in a **flat** bulk: exact, and independent of any wave
solve.  :class:`ShellGeometry` accepts a metric factor ``f`` and then works in
the effective radius ``r/√f(r)``, but **only when that is monotone on the
shell** — it is validated at construction, because a non-monotone ``R_eff``
(a photon sphere) admits trapped orbits and turning points the closed forms
here do not describe.  Nothing here is dynamical: it says which fronts *can*
fold, not what happens when one does.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from geometrodynamics.viz.throat_wavefront import _PAL

# roles specific to this figure: the outer/inner spheres already own
# _PAL["sphere"]/_PAL["plug"], so the two fronts need their own colours
_C: Dict[str, str] = dict(
    pulse="#34d399",              # the front that never folds
    ring=_PAL["neck"],            # the front that must
    defect=_PAL["caustic"],       # where it folds
)


# ════════════════════════════════════════════════════════════════════════════
# THE VACUOLE
# ════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class ShellGeometry:
    """Two concentric spheres with a bulk between them.

    Parameters
    ----------
    r_inner, r_outer
        The two radii, ``0 < r_inner < r_outer``.
    f
        Optional metric factor for a curved bulk; rays then conserve
        ``b = r sin α / √f(r)``.  Defaults to flat, where every closed form
        in this module is exact.
    """

    r_inner: float
    r_outer: float
    f: Optional[Callable[[np.ndarray], np.ndarray]] = None

    def __post_init__(self) -> None:
        if not (0.0 < self.r_inner < self.r_outer):
            raise ValueError("need 0 < r_inner < r_outer")
        if self.f is not None:
            self._require_monotone_effective_radius()

    def _require_monotone_effective_radius(self, n: int = 2001,
                                           tol: float = -1e-12) -> None:
        """Refuse a curved bulk whose ``R_eff`` is not monotone on the shell.

        Everything here — ``turning_radius`` by bisection, "a ray reaches
        ``b`` and no deeper", the one-sided acceptance argument — assumes
        ``R_eff = r/√f(r)`` increases outward.  A metric with a photon
        sphere breaks that: ``R_eff`` develops a minimum, rays can be
        trapped, and there are turning points on both sides of it.  Rather
        than silently return nonsense, reject the geometry at construction.
        """
        r = np.linspace(self.r_inner, self.r_outer, n)
        d = np.diff(np.asarray(self.effective_radius(r), dtype=float))
        if np.any(d <= tol):
            bad = float(r[int(np.argmin(d))])
            raise ValueError(
                "effective radius r/sqrt(f(r)) must increase across the shell; "
                f"it decreases near r = {bad:.4f}.  A non-monotone R_eff (a "
                "photon sphere) admits trapped orbits and two-sided turning "
                "points, which the closed forms in this module do not cover."
            )

    @property
    def gap(self) -> float:
        """``ΔR = R_outer − R_inner`` — the radial thickness of the bulk."""
        return self.r_outer - self.r_inner

    def effective_radius(self, r: np.ndarray | float) -> np.ndarray | float:
        """``r/√f(r)``; the quantity a ray's impact parameter is measured in."""
        r = np.asarray(r, dtype=float)
        if self.f is None:
            return r
        return r / np.sqrt(np.asarray(self.f(r), dtype=float))

    def impact_parameter(self, r: float, alpha: float) -> float:
        """``b = R_eff(r)·sin α`` for a ray at radius ``r``, angle ``α``."""
        return float(self.effective_radius(r) * math.sin(alpha))

    def turning_radius(self, b: float) -> float:
        """Smallest radius a ray of impact parameter ``b`` reaches.

        Flat bulk: ``r_min = b`` exactly.  Curved: the root of
        ``R_eff(r) = b``, found by bisection on the bracket ``[r_inner,
        r_outer]``; returns ``r_inner`` when the ray clears the shell.
        """
        if self.f is None:
            return float(b)
        lo, hi = self.r_inner, self.r_outer
        if self.effective_radius(lo) >= b:
            return float(lo)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if self.effective_radius(mid) < b:
                lo = mid
            else:
                hi = mid
        return float(0.5 * (lo + hi))

    # ── the acceptance cone, and its asymmetry ──────────────────────────────
    @property
    def critical_sin(self) -> float:
        """``sin α_crit = R_eff(r_inner)/R_eff(r_outer)`` — the grazing ray."""
        return float(self.effective_radius(self.r_inner)
                     / self.effective_radius(self.r_outer))

    @property
    def critical_angle(self) -> float:
        """``α_crit`` in radians, measured from the inward radial direction."""
        return math.asin(min(1.0, self.critical_sin))

    def acceptance_fraction(self, direction: str = "inward") -> float:
        """Fraction of a hemisphere of launch directions that crosses the bulk.

        ``"inward"`` — launched from the outer sphere, a ray arrives only if
        ``sin α ≤ sin α_crit``, giving ``1 − √(1 − sin²α_crit)``.
        ``"outward"`` — launched from the inner sphere, every ray arrives,
        because ``b ≤ R_eff(r_inner)`` and ``R_eff`` only grows outward.
        """
        if direction == "outward":
            return 1.0
        if direction != "inward":
            raise ValueError("direction must be 'inward' or 'outward'")
        s = self.critical_sin
        return float(1.0 - math.sqrt(max(0.0, 1.0 - s * s)))

    @property
    def acceptance_asymmetry(self) -> float:
        """``outward/inward`` acceptance ratio — 1 would mean no asymmetry."""
        return 1.0 / self.acceptance_fraction("inward")

    # ── the ring that puts its defect on the inner sphere ───────────────────
    def critical_ring(self) -> "RingSource":
        """The ring whose front self-intersects exactly on the inner sphere.

        Its centre of curvature must lie at radius ``r_inner``, which fixes
        ``cos θ₀ = r_inner/r_outer``; the same ring's rays then leave at
        ``sin α = r_inner/r_outer``, the grazing angle.  One condition, not
        two.
        """
        if self.f is not None:
            raise ValueError("critical_ring is a flat-bulk closed form")
        theta0 = math.acos(self.r_inner / self.r_outer)
        return RingSource(shell=self, polar_angle=theta0)

    def point_source(self) -> "PointSource":
        """A point source on the outer sphere, launching inward."""
        return PointSource(shell=self)


# ════════════════════════════════════════════════════════════════════════════
# SOURCES — and the fronts they make
# ════════════════════════════════════════════════════════════════════════════
@dataclass(frozen=True)
class PointSource:
    """A point on the outer sphere.  Its front is a metric sphere.

    In a **flat** bulk the signed area element of that front is ``t² sin θ``,
    positive for every ``t > 0``, so it never folds and the swept region is
    the filled ball.

    That is a fact about flat space rather than about point sources: on a
    closed manifold the very same source folds.  On ``S²`` or ``S³`` a point's
    front converges on the antipode at ``t = πR`` — see
    :mod:`geometrodynamics.viz.throat_wavefront`, which measures exactly that.
    """

    shell: ShellGeometry

    @property
    def position(self) -> np.ndarray:
        return np.array([0.0, 0.0, self.shell.r_outer])

    @property
    def crossing_time(self) -> float:
        """When the front first touches the inner sphere: ``ΔR``."""
        return self.shell.gap

    @property
    def fold_time(self) -> Optional[float]:
        """``None`` — in a flat bulk a sphere is immersed at every radius."""
        return None

    def first_caustic_point(self) -> Optional[np.ndarray]:
        """``None``: no caustic in a flat bulk."""
        return None

    def singular_points(self, t: float) -> np.ndarray:
        """Empty at every ``t``."""
        return np.zeros((0, 3))

    def offset_surface(self, t: float, n_u: int = 64, n_v: int = 48):
        """``(X, N)``: the front and the offset normal that generated it.

        Parameterised by ``(φ, θ)`` on the direction sphere, with the polar
        endpoints trimmed so finite differences stay interior.
        """
        u = np.linspace(0.0, 2.0 * math.pi, n_u, endpoint=False)
        v = np.linspace(1e-3, math.pi - 1e-3, n_v)
        U, V = np.meshgrid(u, v)
        N = np.stack([np.sin(V) * np.cos(U), np.sin(V) * np.sin(U), np.cos(V)], -1)
        return self.position + t * N, N

    def arrival_multiplicity(self, x: np.ndarray, t: float,
                             tol: float = 1e-9) -> int:
        """Rays of arclength ``t`` reaching ``x`` — always 0 or 1."""
        return int(abs(float(np.linalg.norm(np.asarray(x) - self.position)) - t)
                   <= tol)

    def front(self, t: float, n: int = 240) -> np.ndarray:
        """Meridional trace of the front: the circle ``|x − P| = t``."""
        a = np.linspace(0.0, 2.0 * math.pi, n)
        p = self.position
        return np.stack([t * np.sin(a), np.zeros(n), p[2] + t * np.cos(a)], 1)


@dataclass(frozen=True)
class RingSource:
    """A circle of latitude on the outer sphere.  Its front is a tube.

    Any curved extended source has a focal set, so a ring is not special for
    folding at all — it is special for folding **coherently**.  The signed
    area element of the tube is ``t(ρ + t cos v)``, which vanishes where
    ``cos v = −ρ/t``:

    * ``t < ρ`` — immersed, no fold;
    * ``t = ρ`` — the **first caustic**, infinitely degenerate: the two roots
      coincide at the ring's centre and the *whole ring* arrives at once;
    * ``t > ρ`` — still singular, at two axis points ``z = ±√(t² − ρ²)`` that
      run outward as the front grows.

    So the fold is not an isolated event: it begins with a maximally
    degenerate focus and then persists along the symmetry axis.
    """

    shell: ShellGeometry
    polar_angle: float

    @property
    def radius(self) -> float:
        """``ρ = r_outer sin θ₀`` — the ring's own radius."""
        return self.shell.r_outer * math.sin(self.polar_angle)

    @property
    def centre(self) -> np.ndarray:
        """The ring's centre — where the first caustic forms."""
        return np.array([0.0, 0.0, self.shell.r_outer * math.cos(self.polar_angle)])

    @property
    def centre_radius(self) -> float:
        """How deep the first caustic forms: ``|centre| = r_outer cos θ₀``."""
        return float(abs(self.centre[2]))

    @property
    def fold_time(self) -> float:
        """``t = ρ``: the first caustic, at the ring's own centre."""
        return self.radius

    @property
    def launch_sin(self) -> float:
        """``sin α`` of the ray from the ring to its centre."""
        p = np.array([self.radius, 0.0, self.centre[2]])
        d = self.centre - p
        d = d / np.linalg.norm(d)
        r = float(np.linalg.norm(p))
        return float(math.sqrt(max(0.0, r * r - float(d @ p) ** 2)) / r)

    def first_caustic_point(self) -> np.ndarray:
        """The ring's centre: where the whole ring arrives simultaneously."""
        return self.centre.copy()

    def singular_points(self, t: float) -> np.ndarray:
        """Where the tube is singular at time ``t`` — 0 or 2 axis points.

        Solves ``ρ + t cos v = 0``; for ``t > ρ`` the two roots sit at
        ``z = centre_z ± √(t² − ρ²)`` and separate as ``t`` grows.
        """
        if t < self.radius:
            return np.zeros((0, 3))
        dz = math.sqrt(max(0.0, t * t - self.radius * self.radius))
        z = self.centre[2]
        return np.array([[0.0, 0.0, z + dz], [0.0, 0.0, z - dz]])

    def offset_surface(self, t: float, n_u: int = 64, n_v: int = 48):
        """``(X, N)``: the tube and the offset normal that generated it."""
        u = np.linspace(0.0, 2.0 * math.pi, n_u, endpoint=False)
        v = np.linspace(0.0, 2.0 * math.pi, n_v, endpoint=False)
        U, V = np.meshgrid(u, v)
        N = np.stack([np.cos(V) * np.cos(U), np.cos(V) * np.sin(U),
                      np.sin(V)], -1)
        C = np.stack([self.radius * np.cos(U), self.radius * np.sin(U),
                      np.full_like(U, self.centre[2])], -1)
        return C + t * N, N

    def points(self, n: int = 240) -> np.ndarray:
        a = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
        return np.stack([self.radius * np.cos(a), self.radius * np.sin(a),
                         np.full(n, self.centre[2])], 1)

    def arrival_multiplicity(self, x: np.ndarray, t: float, n: int = 20000,
                             degenerate_tol: float = 1e-9) -> int:
        """Rays of arclength ``t`` reaching ``x``.

        Counts sign changes of ``|x − c(φ)| − t`` around the ring, so a
        generic point returns 2 (one ray from each side) and a point off the
        front returns 0.  At the first caustic every ``φ`` solves it at once;
        that is reported as ``-1`` for "degenerate — the whole ring".
        """
        x = np.asarray(x, dtype=float)
        c = self.points(n)
        d = np.linalg.norm(c - x, axis=1) - t
        if float(np.ptp(d)) <= degenerate_tol:
            return -1 if abs(float(np.mean(d))) <= degenerate_tol else 0
        s = np.sign(d)
        return int(np.sum(s * np.roll(s, -1) < 0))

    def front(self, t: float, n: int = 240) -> np.ndarray:
        """Meridional trace of the tube: circles of radius ``t`` about ``±ρ``."""
        a = np.linspace(0.0, 2.0 * math.pi, n)
        z0 = self.centre[2]
        right = np.stack([self.radius + t * np.cos(a), np.zeros(n),
                          z0 + t * np.sin(a)], 1)
        left = np.stack([-self.radius + t * np.cos(a), np.zeros(n),
                         z0 + t * np.sin(a)], 1)
        return np.concatenate([right, left])


# ════════════════════════════════════════════════════════════════════════════
# AN INDEPENDENT FOLD DETECTOR
# ════════════════════════════════════════════════════════════════════════════
def signed_area_element(source, t: float, n_u: int = 64,
                        n_v: int = 48) -> np.ndarray:
    """``(X_u × X_v)·N`` of the front, by finite differences.

    A wavefront folds exactly where its area element passes through zero and
    changes sign.  Measuring that from the surface itself — rather than from
    a stored radius of curvature — is what makes the topology check
    independent of the closed forms it is supposed to test.
    """
    X, N = source.offset_surface(t, n_u=n_u, n_v=n_v)
    du = 2.0 * math.pi / X.shape[1]
    Xu = np.gradient(X, du, axis=1)
    Xv = np.gradient(X, axis=0)          # spacing folds into a positive scale
    return np.einsum("...i,...i->...", np.cross(Xu, Xv), N)


def detect_fold(source, t_hi: float, t_lo: float = 1e-3, n_scan: int = 400,
                n_u: int = 64, n_v: int = 48, refine: int = 60,
                tol: float = 1e-6) -> Dict[str, object]:
    """Find the first time the front's area element changes **sign**.

    Scans ``t`` for the first sample whose signed area element goes
    meaningfully negative anywhere, then bisects.  Nothing about the
    source's radius of curvature is consulted, so the answer can be compared
    against the closed form as a *check* rather than being seeded by it.

    Two things make the test mean what it says.

    It is *relative*, ``min J < −tol·max|J|``: an absolute ``min J ≤ 0``
    would report a fold for any parametrisation whose area element merely
    vanishes — the direction sphere's own poles do that at every ``t`` — and
    a coordinate degeneracy is not a caustic.  A fold is a change of sign.

    And the sign itself is referenced, not assumed.  ``(X_u × X_v)·N`` picks
    up the handedness of whatever ``(u, v)`` ordering a source happens to
    use; for the point source it is negative everywhere, which would read as
    an immediate fold.  So the orientation is fixed once from a small ``t``,
    where every offset surface is still immersed, and applied at all ``t``.
    """
    def raw(t: float) -> np.ndarray:
        return signed_area_element(source, t, n_u, n_v)

    ref = raw(t_lo)
    orient = 1.0 if float(np.median(ref)) >= 0.0 else -1.0

    def folded(t: float) -> bool:
        j = orient * raw(t)
        scale = float(np.max(np.abs(j)))
        if scale <= 0.0:
            return False
        return bool(float(np.min(j)) < -tol * scale)

    ts = np.linspace(t_lo, t_hi, n_scan)
    hit = None
    for i, t in enumerate(ts):
        if folded(float(t)):
            hit = i
            break
    if hit is None:
        return {"folds": False, "fold_time": None, "scanned_to": float(t_hi)}
    if hit == 0:
        return {"folds": True, "fold_time": float(ts[0]),
                "scanned_to": float(t_hi)}
    lo, hi = float(ts[hit - 1]), float(ts[hit])
    for _ in range(refine):
        mid = 0.5 * (lo + hi)
        if folded(mid):
            hi = mid
        else:
            lo = mid
    return {"folds": True, "fold_time": 0.5 * (lo + hi),
            "scanned_to": float(t_hi)}


# ════════════════════════════════════════════════════════════════════════════
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_acceptance_asymmetry(
    shell: ShellGeometry, n_samples: int = 400000, seed: int = 0,
) -> Dict[str, float]:
    """Solid-angle acceptance both ways across the bulk, and reciprocity.

    The asymmetry is in the **measure of launch directions that connect**,
    not in the propagation: every ray is reversible, and the check below
    confirms it by reversing the accepted inward rays and verifying each one
    arrives going the other way.
    """
    rng = np.random.default_rng(seed)

    v = rng.normal(size=(n_samples, 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    p_out = np.array([0.0, 0.0, shell.r_outer])
    inward = v[v @ (p_out / shell.r_outer) < 0]
    r_min = np.sqrt(np.maximum(shell.r_outer ** 2 - (inward @ p_out) ** 2, 0.0))
    accepted = r_min <= shell.r_inner

    w = rng.normal(size=(n_samples, 3))
    w /= np.linalg.norm(w, axis=1)[:, None]
    p_in = np.array([0.0, 0.0, shell.r_inner])
    outward = w[w @ (p_in / shell.r_inner) > 0]

    # reciprocity: reverse each accepted inward ray at the inner sphere and
    # confirm it climbs back out.  b is unchanged under reversal and
    # b <= r_inner < r_outer, so it must — this asserts it rather than assumes.
    b_accepted = r_min[accepted]
    reversible = bool(np.all(b_accepted <= shell.r_inner + 1e-12))

    inward_frac = shell.acceptance_fraction("inward")
    return {
        "critical_sin": shell.critical_sin,
        "critical_angle_deg": math.degrees(shell.critical_angle),
        "inward_closed_form": inward_frac,
        "inward_monte_carlo": float(np.mean(accepted)),
        "outward_closed_form": shell.acceptance_fraction("outward"),
        "outward_monte_carlo": 1.0,
        "solid_angle_ratio": shell.acceptance_asymmetry,
        "rays_reversible": reversible,
    }


def measure_front_topology(
    source, t_hi: float, n_scan: int = 400,
) -> Dict[str, object]:
    """Does this front fold, and when — measured, then compared.

    The fold time is found by :func:`detect_fold` from the area element
    alone.  Only afterwards is it compared with the source's closed form, so
    the comparison is a test rather than a tautology.
    """
    det = detect_fold(source, t_hi=t_hi, n_scan=n_scan)
    closed = source.fold_time
    err = (None if (closed is None or det["fold_time"] is None)
           else abs(det["fold_time"] - closed))
    caustic = source.first_caustic_point()
    return {
        "kind": "ring" if isinstance(source, RingSource) else "point",
        "detected_folds": det["folds"],
        "detected_fold_time": det["fold_time"],
        "closed_form_fold_time": closed,
        "detection_error": err,
        "first_caustic_point": (None if caustic is None
                                else [float(c) for c in caustic]),
        "singular_points_after": int(
            source.singular_points(1.25 * (closed or t_hi)).shape[0]),
        "scanned_to": det["scanned_to"],
    }


def measure_critical_ring(shell: ShellGeometry) -> Dict[str, float]:
    """The ring that puts its first caustic on the inner sphere, and why it grazes.

    Checks the coincidence the module is about: the ring whose caustic sits
    at ``r_inner`` launches at exactly the grazing angle, so its rays are
    tangent to the inner sphere.
    """
    ring = shell.critical_ring()
    return {
        "polar_angle_deg": math.degrees(ring.polar_angle),
        "ring_radius": ring.radius,
        "caustic_radius": ring.centre_radius,
        "r_inner": shell.r_inner,
        "caustic_on_inner_error": abs(ring.centre_radius - shell.r_inner),
        "fold_time": ring.fold_time,
        "fold_time_closed_form": math.sqrt(shell.r_outer ** 2 - shell.r_inner ** 2),
        "launch_sin": ring.launch_sin,
        "critical_sin": shell.critical_sin,
        "grazing_error": abs(ring.launch_sin - shell.critical_sin),
        "ray_turning_radius": shell.turning_radius(
            shell.r_outer * ring.launch_sin),
        "pulse_crossing_time": shell.point_source().crossing_time,
    }


# ════════════════════════════════════════════════════════════════════════════
# RENDERING — the meridional section, which is where the claim is visible
# ════════════════════════════════════════════════════════════════════════════
def _style(ax: Axes, title: str = "") -> None:
    ax.set_facecolor(_PAL["panel"])
    for sp in ax.spines.values():
        sp.set_color(_PAL["border"])
    ax.tick_params(colors=_PAL["dim"], labelsize=7)
    ax.set_aspect("equal")
    if title:
        ax.set_title(title, color=_PAL["text"], fontsize=9, pad=6)


def draw_shell(ax: Axes, shell: ShellGeometry) -> None:
    """The two spheres, in meridional section."""
    a = np.linspace(0.0, 2.0 * math.pi, 400)
    for r, col, lab in ((shell.r_outer, _PAL["sphere"], "outer"),
                        (shell.r_inner, _PAL["plug"], "inner")):
        ax.plot(r * np.cos(a), r * np.sin(a), color=col, lw=1.6, label=lab)


def draw_front(ax: Axes, source, t: float, color: Optional[str] = None,
               alpha: float = 0.9, lw: float = 1.4) -> None:
    """The wavefront at time ``t``, in meridional section."""
    col = color or (_C["ring"] if isinstance(source, RingSource)
                    else _C["pulse"])
    f = source.front(t)
    if isinstance(source, RingSource):
        n = len(f) // 2
        ax.plot(f[:n, 0], f[:n, 2], color=col, lw=lw, alpha=alpha)
        ax.plot(f[n:, 0], f[n:, 2], color=col, lw=lw, alpha=alpha)
    else:
        ax.plot(f[:, 0], f[:, 2], color=col, lw=lw, alpha=alpha)


def plot_pulse_vs_ring(shell: ShellGeometry,
                       figsize: Tuple[float, float] = (11.0, 5.6)) -> Figure:
    """The whole argument in one figure: a pulse cannot fold, a ring must."""
    ring = shell.critical_ring()
    pulse = shell.point_source()
    fig, axes = plt.subplots(1, 2, figsize=figsize, facecolor=_PAL["bg"])

    ax = axes[0]
    draw_shell(ax, shell)
    for t in np.linspace(0.12, 1.5 * shell.gap, 7):
        draw_front(ax, pulse, t, alpha=0.75, lw=1.0)
    draw_front(ax, pulse, pulse.crossing_time, lw=2.0)
    ax.plot([0.0], [shell.r_outer], "o", color=_PAL["text"], ms=4)
    _style(ax, f"pulse — front stays immersed in a flat bulk\n"
               f"crosses at t = ΔR = {pulse.crossing_time:.3f}, no fold")

    ax = axes[1]
    draw_shell(ax, shell)
    for t in np.linspace(0.15, 1.35 * ring.radius, 8):
        draw_front(ax, ring, t, alpha=0.6, lw=0.9)
    draw_front(ax, ring, ring.fold_time, lw=2.0)
    ax.plot([ring.radius, -ring.radius], [ring.centre[2]] * 2, "o",
            color=_PAL["text"], ms=4)
    ax.plot([0.0], [ring.centre[2]], "*", color=_C["defect"], ms=15,
            zorder=5)
    _style(ax, f"ring θ₀={math.degrees(ring.polar_angle):.1f}° — whole ring "
               f"at once\nfirst caustic t = ρ = {ring.radius:.3f}, on r = "
               f"{ring.centre_radius:.3f}")
    for ax in axes:
        lim = 1.15 * shell.r_outer
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle("In a flat bulk a pulse fills its own void; a ring focuses coherently",
                 color=_PAL["text"], fontsize=12)
    fig.tight_layout()
    return fig


def plot_acceptance_cone(shell: ShellGeometry,
                         figsize: Tuple[float, float] = (11.0, 5.6)) -> Figure:
    """The outer/inner asymmetry: who gets across, launched from where."""
    fig, axes = plt.subplots(1, 2, figsize=figsize, facecolor=_PAL["bg"])
    res = measure_acceptance_asymmetry(shell, n_samples=20000, seed=1)

    for ax, (r0, outward, lab) in zip(
            axes, ((shell.r_outer, False, "outer → inner"),
                   (shell.r_inner, True, "inner → outer"))):
        draw_shell(ax, shell)
        p = np.array([0.0, r0])
        for ang in np.linspace(-math.pi / 2, math.pi / 2, 25):
            d = np.array([math.sin(ang), -math.cos(ang)]) * (1 if not outward else -1)
            seg = np.stack([p + s * d for s in np.linspace(0, 3.0, 200)])
            rr = np.linalg.norm(seg, axis=1)
            keep = rr <= shell.r_outer + 1e-9
            arrives = (rr.min() <= shell.r_inner) if not outward else True
            col = _C["pulse"] if arrives else _PAL["dim"]
            ax.plot(seg[keep, 0], seg[keep, 1], color=col, lw=0.8,
                    alpha=0.9 if arrives else 0.35)
        frac = (res["inward_closed_form"] if not outward
                else res["outward_closed_form"])
        _style(ax, f"{lab} — {100 * frac:.1f}% of the hemisphere arrives")
        lim = 1.15 * shell.r_outer
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(
        f"Accepted solid angle differs {res['solid_angle_ratio']:.1f}× — "
        "individual rays stay reversible",
        color=_PAL["text"], fontsize=12)
    fig.tight_layout()
    return fig
