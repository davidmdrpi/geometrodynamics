"""
Why a throat needs a *ring*: front topology across the bulk.

The companion module :mod:`geometrodynamics.viz.throat_wavefront` runs waves on
a surface that *already has* a throat.  This one asks the prior question — what
kind of wavefront can make one — and answers it with the differential geometry
of wavefronts rather than with a simulation.

The vacuole
───────────
Two concentric spheres, ``r = R_inner`` and ``r = R_outer``, with the bulk
between them (the static bound vacuole of the programme; ``ΔR = R_outer −
R_inner``).  Straight rays cross the bulk; the conserved quantity along a ray
is the impact parameter about the centre,

```
b = r sin α          (α = angle to the radial direction)
```

so a ray launched inward from the outer sphere reaches radius ``b`` and no
deeper.

A pulse cannot make a defect
────────────────────────────
The wavefront of a **point** source is the metric sphere ``|x − P| = t``.  A
sphere is embedded for every ``t``: the front never touches itself, and the
region behind it is the filled ball.  In the programme's own words the pulse
"just fills its own void".  The reason is exact and is a statement about focal
sets: *the focal set of a point is empty*, so the front has nothing to fold on.
It crosses the bulk at ``t = ΔR`` and that is all it does.

A ring must
───────────
The wavefront of a **circle** of radius ``ρ`` is the offset (tube) surface of
that circle, and a curve *does* have a focal set: the locus of its centres of
curvature.  For a circle every point shares the *same* centre of curvature — the
centre — so the focal set collapses to a single point and the whole ring arrives
there at once, at

```
t = ρ
```

That is a degenerate caustic of infinite multiplicity: not a smooth focus but a
point where the front ceases to be embedded.  A codimension-2 defect of the
wavefront, made by geometry alone.

The two conditions coincide
───────────────────────────
Ask for the ring whose defect lands *on the inner sphere*.  Its centre must sit
at radius ``R_inner``, so it is the circle of polar angle

```
cos θ₀ = R_inner / R_outer
```

and its rays leave the outer sphere at ``sin α = R_inner / R_outer`` — which is
exactly the **critical acceptance angle**, the ray tangent to the inner sphere.
The ring that focuses on the throat and the ray that grazes it are the same
condition.  The defect forms at

```
t = ρ = √(R_outer² − R_inner²)
```

The asymmetry between the two surfaces
──────────────────────────────────────
Because ``b = r sin α`` is conserved and ``r`` decreases going in:

* **outer → inner** only rays with ``sin α ≤ R_inner/R_outer`` arrive at all;
  the rest turn around in the bulk.  The accepted fraction of the inward
  hemisphere is ``1 − √(1 − (R_inner/R_outer)²)`` — about **19%** at the
  programme's ``ΔR``.
* **inner → outer** every ray arrives, because ``b ≤ R_inner < R_outer`` always.
  **100%**.

Same bulk, same rays, opposite directions, and the two surfaces do not see each
other the same way.  That is the inner/outer asymmetry in its plainest form: it
is not a broken symmetry of the sphere, it is the ordering of the two radii.

Scope
─────
Pure ray/front geometry in a flat bulk: exact, and independent of any wave
solve.  A curved bulk replaces ``r`` by the effective radius ``r/√f(r)``
everywhere below — :class:`ShellGeometry` accepts an ``f`` and the structure
carries over — but the closed forms quoted here are the flat ones.  Nothing
here is dynamical: it says which fronts *can* fold, not what happens when one
does.
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

    The focal set of a point is empty, so the front never folds: it is an
    embedded sphere for every ``t``, and the swept region is the filled ball.
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
    def self_intersection_time(self) -> Optional[float]:
        """``None`` — a sphere is embedded at every radius."""
        return None

    def focal_set(self) -> np.ndarray:
        """Empty: a point has no centre of curvature."""
        return np.zeros((0, 3))

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

    Every point of a circle shares one centre of curvature, so the focal set
    is that single point and the *entire* ring arrives there simultaneously
    at ``t = ρ``.  The front stops being embedded — a degenerate caustic of
    infinite multiplicity, which is the defect.
    """

    shell: ShellGeometry
    polar_angle: float

    @property
    def radius(self) -> float:
        """``ρ = r_outer sin θ₀`` — the ring's own radius."""
        return self.shell.r_outer * math.sin(self.polar_angle)

    @property
    def centre(self) -> np.ndarray:
        """The ring's centre — its whole focal set."""
        return np.array([0.0, 0.0, self.shell.r_outer * math.cos(self.polar_angle)])

    @property
    def centre_radius(self) -> float:
        """How deep the defect forms: ``|centre| = r_outer cos θ₀``."""
        return float(abs(self.centre[2]))

    @property
    def self_intersection_time(self) -> float:
        """``t = ρ``: when the front folds onto its own centre."""
        return self.radius

    @property
    def launch_sin(self) -> float:
        """``sin α`` of the ray from the ring to its centre."""
        p = np.array([self.radius, 0.0, self.centre[2]])
        d = self.centre - p
        d = d / np.linalg.norm(d)
        r = float(np.linalg.norm(p))
        return float(math.sqrt(max(0.0, r * r - float(d @ p) ** 2)) / r)

    def focal_set(self) -> np.ndarray:
        """A single point: the ring's centre."""
        return self.centre.reshape(1, 3)

    def points(self, n: int = 240) -> np.ndarray:
        a = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
        return np.stack([self.radius * np.cos(a), self.radius * np.sin(a),
                         np.full(n, self.centre[2])], 1)

    def arrival_multiplicity(self, x: np.ndarray, t: float, n: int = 20000,
                             degenerate_tol: float = 1e-9) -> int:
        """Rays of arclength ``t`` reaching ``x``.

        Counts sign changes of ``|x − c(φ)| − t`` around the ring, so a
        generic point returns 2 (one ray from each side) and a point off the
        front returns 0.  On the focal point every ``φ`` is a solution at
        once; that is reported as ``-1`` for "degenerate — the whole ring".
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
# MEASUREMENTS
# ════════════════════════════════════════════════════════════════════════════
def measure_acceptance_asymmetry(
    shell: ShellGeometry, n_samples: int = 400000, seed: int = 0,
) -> Dict[str, float]:
    """Closed form against Monte-Carlo, both directions across the bulk."""
    rng = np.random.default_rng(seed)

    def sample(r_launch: float, outward: bool) -> float:
        v = rng.normal(size=(n_samples, 3))
        v /= np.linalg.norm(v, axis=1)[:, None]
        p = np.array([0.0, 0.0, r_launch])
        nrm = p / r_launch
        keep = (v @ nrm > 0) if outward else (v @ nrm < 0)
        v = v[keep]
        r_min = np.sqrt(np.maximum(r_launch ** 2 - (v @ p) ** 2, 0.0))
        if outward:
            return 1.0                      # R_eff grows outward: all escape
        return float(np.mean(r_min <= shell.r_inner))

    inward = shell.acceptance_fraction("inward")
    return {
        "critical_sin": shell.critical_sin,
        "critical_angle_deg": math.degrees(shell.critical_angle),
        "inward_closed_form": inward,
        "inward_monte_carlo": sample(shell.r_outer, outward=False),
        "outward_closed_form": shell.acceptance_fraction("outward"),
        "outward_monte_carlo": sample(shell.r_inner, outward=True),
        "asymmetry_ratio": shell.acceptance_asymmetry,
    }


def measure_front_topology(
    source, t_values: np.ndarray, off_axis: float = 0.02,
) -> Dict[str, object]:
    """Does this front ever stop being embedded, and when?

    Samples the arrival multiplicity just off the axis and on it.  A point
    source answers "never"; a ring answers with its own radius.
    """
    ring = isinstance(source, RingSource)
    ts = np.asarray(t_values, dtype=float)
    # the fold is a measure-zero event in t, so sample it rather than hope
    if source.self_intersection_time is not None:
        ts = np.unique(np.append(ts, source.self_intersection_time))
    rows: List[dict] = []
    for t in ts:
        if ring:
            axis = source.centre.copy()
            near = axis + np.array([off_axis, 0.0, 0.0])
            rows.append({"t": float(t),
                         "on_axis": source.arrival_multiplicity(axis, t),
                         "off_axis": source.arrival_multiplicity(near, t)})
        else:
            p = source.position + np.array([0.0, 0.0, -t])
            rows.append({"t": float(t), "on_axis": 1,
                         "off_axis": source.arrival_multiplicity(p, t)})
    degenerate = [r["t"] for r in rows if r["on_axis"] == -1]
    return {
        "kind": "ring" if ring else "point",
        "self_intersection_time": source.self_intersection_time,
        "focal_set_size": int(source.focal_set().shape[0]),
        "degenerate_times": degenerate,
        "ever_self_intersects": bool(source.self_intersection_time is not None),
        "rows": rows,
    }


def measure_critical_ring(shell: ShellGeometry) -> Dict[str, float]:
    """The ring that lands its defect on the inner sphere, and why it grazes.

    Checks the coincidence the module is about: the ring whose focal point
    sits at ``r_inner`` launches at exactly the grazing angle, so its rays
    are tangent to the inner sphere.
    """
    ring = shell.critical_ring()
    return {
        "polar_angle_deg": math.degrees(ring.polar_angle),
        "ring_radius": ring.radius,
        "defect_radius": ring.centre_radius,
        "r_inner": shell.r_inner,
        "defect_on_inner_error": abs(ring.centre_radius - shell.r_inner),
        "defect_time": ring.self_intersection_time,
        "defect_time_closed_form": math.sqrt(shell.r_outer ** 2 - shell.r_inner ** 2),
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
    _style(ax, f"pulse — front is a sphere, focal set empty\n"
               f"crosses at t = ΔR = {pulse.crossing_time:.3f}, never folds")

    ax = axes[1]
    draw_shell(ax, shell)
    for t in np.linspace(0.15, 1.35 * ring.radius, 8):
        draw_front(ax, ring, t, alpha=0.6, lw=0.9)
    draw_front(ax, ring, ring.self_intersection_time, lw=2.0)
    ax.plot([ring.radius, -ring.radius], [ring.centre[2]] * 2, "o",
            color=_PAL["text"], ms=4)
    ax.plot([0.0], [ring.centre[2]], "*", color=_C["defect"], ms=15,
            zorder=5)
    _style(ax, f"ring θ₀={math.degrees(ring.polar_angle):.1f}° — focal set is "
               f"one point\ndefect at t = ρ = {ring.radius:.3f}, on r = "
               f"{ring.centre_radius:.3f}")
    for ax in axes:
        lim = 1.15 * shell.r_outer
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle("A pulse fills its own void; a ring folds onto its centre",
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
        f"Same bulk, opposite directions: {res['asymmetry_ratio']:.1f}× asymmetry",
        color=_PAL["text"], fontsize=12)
    fig.tight_layout()
    return fig
