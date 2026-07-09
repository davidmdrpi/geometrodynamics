"""
Antipodal focusing — open vs closed wave refocusing (a visual study).

A wave launched on an *open* membrane runs to the edge and disperses into
its own vacated space; nothing forms.  A wave on a *closed* surface cannot
escape — it is forced to reconverge at the antipode, and the geometry
itself decides whether that focus disperses or nucleates an object.  This
module makes the contrast watchable and testable with two live classical
wave solves:

* :class:`PlaneWaveSim` — the 2-D radial wave ``u_tt = u_rr + (1/r) u_r``
  on an open disk (a sponge layer drains the rim).  The ring expands, hits
  the edge and leaves; the peak decays and never recovers.
* :class:`SphereWaveSim` — the axisymmetric wave on a 2-sphere,
  ``u_tt = u_θθ + cot θ · u_θ`` (the Laplace–Beltrami operator).  A cap
  pulse at the north pole grows to the equator, then *splays inward* to the
  south pole and refocuses at ``t = πR`` (the great-circle half-period) —
  the analogue of the zonal ``S³`` result ``ψ = f / sin χ`` refocusing at
  the antipode ``π−χ₀`` (PR #166).

The physical amplitude carries the geometric focusing factor
:func:`focusing_factor` ``= 1/sin θ``, which diverges at the poles — the
caustic.  Above a critical mass the focus persists into a non-orientable
throat (the ``C`` inner/outer swap; ``Σc₁ = 0`` conjugate pairs, PR #58);
below it the field revives.

The sim classes run headlessly and are directly testable.  The
``draw_*`` functions render into an existing :class:`~matplotlib.figure.Figure`;
``plot_*_panel`` wrappers create one; :func:`run_animation` drives a live
plane-vs-sphere view.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure


# ── Palette (matches geometry_panels) ───────────────────────────────────────
_PAL: Dict[str, str] = dict(
    bg="#0a0e17",
    panel="#111827",
    border="#1e293b",
    text="#e2e8f0",
    dim="#64748b",
    plane="#60a5fa",     # open geometry, cool
    sphere="#fbbf24",    # closed geometry, focus
    caustic="#fb7185",   # the refocus flash
    object="#f472b6",    # the throat / C-swap, non-orientable
    persist="#34d399",   # above threshold
)

# Signed-field diverging colormap: trough (blue) → rest (near-bg) → crest (amber)
_FIELD_CMAP = LinearSegmentedColormap.from_list(
    "antipodal_field",
    ["#2b6cb0", "#1a3050", "#0a0e17", "#3a2c12", "#d98a2b", "#fbbf24"],
)

# ── Physical / numerical constants ──────────────────────────────────────────
C: float = 1.0                 # wave speed (units R = c = 1)
PERIOD: float = 2.0 * math.pi  # great-circle period on the unit sphere
HALF_PERIOD: float = math.pi   # antipodal refocus time  t = πR

_N_PLANE: int = 260            # radial grid points (open disk)
_N_SPHERE: int = 300           # polar grid points (closed shell)
_CFL: float = 0.45             # Courant factor
_PULSE_R: float = 0.055        # launch pulse width (plane, in units of R_max)
_PULSE_TH: float = 0.10        # launch pulse width (sphere, radians)
_SPONGE_R0: float = 0.82       # plane sponge inner edge
_SPONGE_K: float = 0.10        # plane sponge strength


def focusing_factor(theta: np.ndarray | float) -> np.ndarray | float:
    """Geometric focusing factor ``1/sin θ`` on the sphere.

    A latitude wavefront at polar angle ``θ`` has circumference
    ``2π sin θ``; the same energy squeezed onto that shrinking circle makes
    the physical amplitude grow as ``1/sin θ`` — divergent (a caustic) at
    both poles.  On the open plane the analogue is the wavefront
    circumference ``2π r``, which only grows: no forced focus.
    """
    return 1.0 / np.sin(theta)


# ════════════════════════════════════════════════════════════════════════════
# PLANE — radial wave on an open disk (leapfrog, sponge rim)
# ════════════════════════════════════════════════════════════════════════════
class PlaneWaveSim:
    """Radial wave ``u_tt = u_rr + (1/r) u_r`` on the open disk ``r ∈ [0,1]``.

    A Gaussian cap launched at the centre expands outward; a sponge layer
    near the rim drains the energy (an open boundary).  The peak amplitude
    decays monotonically — there is no return path.

    Parameters
    ----------
    n
        Number of radial grid points.
    """

    def __init__(self, n: int = _N_PLANE) -> None:
        self.n = n
        self.dr = 1.0 / n
        self.dt = _CFL * self.dr
        self.r = np.linspace(0.0, 1.0, n + 1)
        self._sponge = np.where(
            self.r > _SPONGE_R0,
            1.0 - _SPONGE_K * (self.r - _SPONGE_R0) / (1.0 - _SPONGE_R0),
            1.0,
        )
        self.t = 0.0
        self.reset()

    def reset(self) -> None:
        """Launch the cap pulse at the centre; zero initial velocity."""
        self.u = np.exp(-((self.r / _PULSE_R) ** 2))
        self.u_prev = self.u.copy()
        self.t = 0.0

    def _laplacian(self, u: np.ndarray) -> np.ndarray:
        lap = np.zeros_like(u)
        r = self.r
        dr = self.dr
        lap[1:-1] = (
            (u[2:] - 2 * u[1:-1] + u[:-2]) / dr ** 2
            + (1.0 / r[1:-1]) * (u[2:] - u[:-2]) / (2 * dr)
        )
        lap[0] = 4.0 * (u[1] - u[0]) / dr ** 2          # regular pole limit
        lap[-1] = lap[-2]
        return lap

    def step(self) -> None:
        """Advance one leapfrog step (drained at the rim)."""
        u_new = (
            2 * self.u - self.u_prev + (C * self.dt) ** 2 * self._laplacian(self.u)
        ) * self._sponge
        self.u_prev = self.u
        self.u = u_new
        self.t += self.dt

    def run(self, n_steps: int) -> None:
        for _ in range(n_steps):
            self.step()

    def advance_to(self, t_target: float) -> None:
        while self.t < t_target:
            self.step()

    def peak(self) -> Tuple[float, float]:
        """Return ``(|u|_max, r_at_max)``."""
        idx = int(np.argmax(np.abs(self.u)))
        return float(np.abs(self.u[idx])), float(self.r[idx])

    def sample(self, x: np.ndarray) -> np.ndarray:
        """Interpolate ``u`` at radii ``x`` (clamped to the grid)."""
        return np.interp(np.clip(x, 0.0, 1.0), self.r, self.u)


# ════════════════════════════════════════════════════════════════════════════
# SPHERE — axisymmetric wave on a 2-sphere (leapfrog, closed)
# ════════════════════════════════════════════════════════════════════════════
class SphereWaveSim:
    """Axisymmetric wave ``u_tt = u_θθ + cot θ · u_θ`` on the unit 2-sphere.

    A Gaussian cap launched at the north pole ``θ=0`` grows to the equator,
    then contracts (splays outer → inner) to the south pole ``θ=π`` and
    refocuses there — the caustic — at ``t = πR``, before reviving.  The
    surface is closed: energy is conserved and the wave is handed back to
    itself.

    Parameters
    ----------
    n
        Number of polar-angle grid points.
    """

    def __init__(self, n: int = _N_SPHERE) -> None:
        self.n = n
        self.dth = math.pi / n
        self.dt = _CFL * self.dth
        self.theta = np.linspace(0.0, math.pi, n + 1)
        self._sin = np.sin(self.theta)
        with np.errstate(divide="ignore", invalid="ignore"):
            self._cot = np.cos(self.theta) / self._sin
        self.t = 0.0
        self.reset()

    def reset(self) -> None:
        """Launch the cap pulse at the north pole; zero initial velocity."""
        self.u = np.exp(-((self.theta / _PULSE_TH) ** 2))
        self.u_prev = self.u.copy()
        self.t = 0.0

    def _laplacian(self, u: np.ndarray) -> np.ndarray:
        lap = np.zeros_like(u)
        dth = self.dth
        lap[1:-1] = (
            (u[2:] - 2 * u[1:-1] + u[:-2]) / dth ** 2
            + self._cot[1:-1] * (u[2:] - u[:-2]) / (2 * dth)
        )
        lap[0] = 4.0 * (u[1] - u[0]) / dth ** 2         # north-pole limit
        lap[-1] = 4.0 * (u[-2] - u[-1]) / dth ** 2      # south-pole limit
        return lap

    def step(self) -> None:
        """Advance one leapfrog step (closed — no boundary loss)."""
        u_new = 2 * self.u - self.u_prev + (C * self.dt) ** 2 * self._laplacian(self.u)
        self.u_prev = self.u
        self.u = u_new
        self.t += self.dt

    def run(self, n_steps: int) -> None:
        for _ in range(n_steps):
            self.step()

    def advance_to(self, t_target: float) -> None:
        while self.t < t_target:
            self.step()

    def peak(self) -> Tuple[float, float]:
        """Return ``(|u|_max, θ_at_max)``."""
        idx = int(np.argmax(np.abs(self.u)))
        return float(np.abs(self.u[idx])), float(self.theta[idx])

    def sample(self, th: np.ndarray) -> np.ndarray:
        """Interpolate ``u`` at polar angles ``th`` (clamped to the grid)."""
        return np.interp(np.clip(th, 0.0, math.pi), self.theta, self.u)

    def focus_intensity(self, eps: float = 0.11) -> float:
        """Physical caustic intensity ``|u|/√(sin θ)`` just off the far pole.

        Measured at ``θ = π − eps`` to avoid the exact-pole grid
        singularity; this is the amplitude that must exceed the critical
        mass for a throat to nucleate.
        """
        th = math.pi - eps
        return float(abs(self.sample(np.array([th]))[0]) / math.sqrt(math.sin(th)))

    def energy(self) -> float:
        """Total wave energy ``∫(u_t² + u_θ²) sin θ dθ`` (conserved: closed)."""
        u_t = (self.u - self.u_prev) / self.dt
        u_th = np.gradient(self.u, self.dth)
        dens = (u_t ** 2 + u_th ** 2) * self._sin
        return float(np.trapezoid(dens, self.theta))


# ── Diagnostics (pure-numeric; no rendering) ────────────────────────────────
@dataclass
class RefocusResult:
    """Outcome of :func:`measure_refocus`."""

    refocus_time: float
    refocus_theta: float
    peak_intensity: float
    plane_peak_ratio: float     # plane |u|_max(t_refocus) / |u|_max(0)


def measure_refocus(
    n_periods: float = 0.6, sphere_n: int = _N_SPHERE, plane_n: int = _N_PLANE
) -> RefocusResult:
    """Run both sims and locate the sphere's antipodal refocus.

    Returns the time and location of the first far-hemisphere peak (the
    caustic), its physical intensity, and — for contrast — how far the
    plane's peak has decayed by then.  The refocus time is expected near
    ``t = πR`` and the location near ``θ = π`` (PR #166's ``ψ = f/sinχ``
    refocus at the antipode).
    """
    sph = SphereWaveSim(n=sphere_n)
    pln = PlaneWaveSim(n=plane_n)
    p0, _ = pln.peak()
    t_stop = n_periods * PERIOD
    best_t, best_th, best_pk, best_I = 0.0, 0.0, 0.0, 0.0
    while sph.t < t_stop:
        sph.step()
        pk, th = sph.peak()
        # a genuine refocus: peak in the far hemisphere, after leaving the launch pole
        if th > 0.5 * math.pi and sph.t > 1.0 and pk > best_pk:
            best_pk, best_t, best_th = pk, sph.t, th
            best_I = sph.focus_intensity()
    pln.advance_to(best_t if best_t > 0 else t_stop)
    p_now, _ = pln.peak()
    return RefocusResult(
        refocus_time=best_t,
        refocus_theta=best_th,
        peak_intensity=best_I,
        plane_peak_ratio=p_now / p0,
    )


# ════════════════════════════════════════════════════════════════════════════
# Rendering helpers
# ════════════════════════════════════════════════════════════════════════════
def _style_ax(ax: Axes) -> Axes:
    ax.set_facecolor(_PAL["bg"])
    for spine in ax.spines.values():
        spine.set_color(_PAL["border"])
    ax.tick_params(colors=_PAL["dim"], labelsize=8)
    ax.title.set_color(_PAL["text"])
    return ax


def _style_fig(fig: Figure) -> Figure:
    fig.patch.set_facecolor(_PAL["bg"])
    return fig


def _disk_field(sim: PlaneWaveSim, px: int = 220) -> np.ndarray:
    """Revolve the radial profile into a square field image (NaN outside)."""
    ax = np.linspace(-1.0, 1.0, px)
    xx, yy = np.meshgrid(ax, ax)
    rr = np.sqrt(xx ** 2 + yy ** 2)
    img = sim.sample(rr)
    img[rr > 1.0] = np.nan
    return img


def _sphere_field(sim: SphereWaveSim, px: int = 220, tilt: float = 0.5) -> np.ndarray:
    """Orthographic projection, oriented so the *antipode* (the refocus
    point) faces the viewer: the launch pole sits on the far side, so the
    wave appears at the limb, converges to the front-centre caustic at
    ``t = πR``, and re-expands.  Lambert-shaded; NaN off-disk."""
    ax = np.linspace(-1.0, 1.0, px)
    xx, yy = np.meshgrid(ax, ax)
    rr2 = xx ** 2 + yy ** 2
    inside = rr2 <= 1.0
    zz = np.sqrt(np.clip(1.0 - rr2, 0.0, None))
    px_, py_, pz_ = xx, -yy, zz
    # pole axis points AWAY from the viewer (−z): antipode θ=π faces us
    pax = np.array([0.0, math.sin(tilt), -math.cos(tilt)])
    dot = np.clip(px_ * pax[0] + py_ * pax[1] + pz_ * pax[2], -1.0, 1.0)
    theta = np.arccos(dot)
    img = sim.sample(theta)
    shade = 0.5 + 0.5 * np.clip(px_ * (-0.35) + py_ * 0.45 + pz_ * 0.8, 0.0, None)
    img = img * shade
    img[~inside] = np.nan
    return img


def _peak_strip(ax: Axes, times: np.ndarray, peaks: np.ndarray, color: str,
                mark_pi: bool) -> None:
    _style_ax(ax)
    ax.plot(times, peaks, color=color, lw=1.8)
    ax.set_xlim(0, PERIOD)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("t / R", color=_PAL["dim"], fontsize=8)
    if mark_pi:
        ax.axvline(HALF_PERIOD, color=_PAL["caustic"], lw=1.0, ls="--", alpha=0.8)
        ax.text(HALF_PERIOD + 0.06, 0.9, "πR", color=_PAL["caustic"],
                fontsize=8, family="monospace")


def _vlim(sim) -> float:
    # track the current field with a low floor so the (small-amplitude but
    # sharply localized) antipodal caustic still reads brightly
    return max(0.10, float(np.max(np.abs(sim.u))))


# ── Panel: plane vs sphere contrast ─────────────────────────────────────────
def draw_contrast(
    fig: Figure,
    plane: Optional[PlaneWaveSim] = None,
    sphere: Optional[SphereWaveSim] = None,
    *,
    history: Optional[Dict[str, list]] = None,
) -> Dict[str, Axes]:
    """Render the plane-vs-sphere field comparison into ``fig``.

    ``plane`` / ``sphere`` default to fresh sims advanced to the antipodal
    refocus (``t ≈ πR``) so a still frame already tells the story.  Pass
    ``history`` = ``{"t": [...], "plane": [...], "sphere": [...]}`` (peak
    amplitude tracks) to draw the strip charts; otherwise a single point is
    shown.
    """
    fig.clear()
    _style_fig(fig)
    if plane is None:
        plane = PlaneWaveSim()
        plane.advance_to(HALF_PERIOD)
    if sphere is None:
        sphere = SphereWaveSim()
        sphere.advance_to(HALF_PERIOD)

    gs = fig.add_gridspec(2, 2, height_ratios=(3.0, 1.0), hspace=0.32, wspace=0.16)
    ax_p = fig.add_subplot(gs[0, 0])
    ax_s = fig.add_subplot(gs[0, 1])
    ax_ps = fig.add_subplot(gs[1, 0])
    ax_ss = fig.add_subplot(gs[1, 1])

    for ax in (ax_p, ax_s):
        ax.set_facecolor(_PAL["bg"])
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color(_PAL["border"])

    ax_p.imshow(_disk_field(plane), cmap=_FIELD_CMAP, vmin=-_vlim(plane),
                vmax=_vlim(plane), origin="lower", extent=(-1, 1, -1, 1))
    ax_p.add_patch(plt.Circle((0, 0), 1.0, fill=False, ec=_PAL["plane"], lw=1.2, alpha=0.4))
    pk_p, loc_p = plane.peak()
    ax_p.set_title(f"plane · open disk   |u|={pk_p:.2f}", color=_PAL["plane"])

    ax_s.imshow(_sphere_field(sphere), cmap=_FIELD_CMAP, vmin=-_vlim(sphere),
                vmax=_vlim(sphere), origin="lower", extent=(-1, 1, -1, 1))
    ax_s.add_patch(plt.Circle((0, 0), 1.0, fill=False, ec=_PAL["sphere"], lw=1.2, alpha=0.5))
    pk_s, loc_s = sphere.peak()
    ax_s.set_title(f"sphere · closed shell   |u|={pk_s:.2f}", color=_PAL["sphere"])

    if history and history.get("t"):
        t = np.asarray(history["t"])
        _peak_strip(ax_ps, t, np.asarray(history["plane"]), _PAL["plane"], False)
        _peak_strip(ax_ss, t, np.asarray(history["sphere"]), _PAL["sphere"], True)
    else:
        for ax, sim, col, mk in ((ax_ps, plane, _PAL["plane"], False),
                                  (ax_ss, sphere, _PAL["sphere"], True)):
            pk, _ = sim.peak()
            _peak_strip(ax, np.array([sim.t]), np.array([pk]), col, mk)
    ax_ps.set_title("plane peak — escapes", color=_PAL["dim"], fontsize=9)
    ax_ss.set_title("sphere peak — refocus at πR", color=_PAL["dim"], fontsize=9)

    fig.text(0.02, 0.015,
             "open geometry disperses · closed geometry refocuses at the antipode",
             color=_PAL["text"], family="monospace", fontsize=9)
    return dict(ax_plane=ax_p, ax_sphere=ax_s, ax_plane_strip=ax_ps, ax_sphere_strip=ax_ss)


def plot_contrast_panel(*, t: float = HALF_PERIOD, fig: Optional[Figure] = None) -> Figure:
    """Standalone wrapper: plane vs sphere at time ``t`` (default ``πR``)."""
    fig = fig if fig is not None else plt.figure(figsize=(11, 6.2))
    plane = PlaneWaveSim()
    plane.advance_to(t)
    sphere = SphereWaveSim()
    sphere.advance_to(t)
    draw_contrast(fig, plane, sphere)
    return fig


# ── Panel: the geometric focusing factor ────────────────────────────────────
def draw_focusing_factor(fig: Figure) -> Dict[str, Axes]:
    """Render the wavefront-circumference / ``1/sin θ`` comparison."""
    fig.clear()
    _style_fig(fig)
    ax = fig.add_subplot(111)
    _style_ax(ax)
    ax.grid(True, color=_PAL["border"], lw=0.4, alpha=0.5)
    u = np.linspace(0, 1, 240)
    ax.plot(u, u, color=_PAL["plane"], lw=2.2,
            label="plane  2πr  → disperses off the edge")
    th = np.linspace(1e-3, math.pi - 1e-3, 400)
    ax.plot(th / math.pi, np.sin(th), color=_PAL["sphere"], lw=2.2,
            label="sphere  2π sinθ  → collapses to the antipode")
    ff = np.clip(focusing_factor(th), 0, 3.0) / 3.0
    ax.plot(th / math.pi, ff, color=_PAL["caustic"], lw=1.8, ls="--",
            label="1/sinθ → ∞  (caustic at both poles)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["θ=0 (pole)", "π/2 (equator)", "π (antipode)"])
    ax.set_title("Focusing factor — wavefront circumference & 1/sinθ", color=_PAL["text"])
    leg = ax.legend(loc="upper center", facecolor=_PAL["panel"],
                    edgecolor=_PAL["border"], fontsize=8.5)
    for txt in leg.get_texts():
        txt.set_color(_PAL["text"])
    return dict(ax=ax)


def plot_focusing_panel(*, fig: Optional[Figure] = None) -> Figure:
    """Standalone wrapper for :func:`draw_focusing_factor`."""
    fig = fig if fig is not None else plt.figure(figsize=(9, 4.2))
    draw_focusing_factor(fig)
    return fig


# ── Panel: the object at the focus (disperse vs throat) ─────────────────────
def draw_focus_object(fig: Figure, nucleated: bool = True) -> Dict[str, Axes]:
    """Cross-section schematic at the antipode: disperse or nucleate a throat.

    ``nucleated`` chooses the state — above the critical mass the inner and
    outer sheets pinch into a non-orientable throat (the ``C`` inner/outer
    swap, ``Σc₁ = 0`` conjugate pair); below it the focus passes through and
    the field revives.
    """
    fig.clear()
    _style_fig(fig)
    ax = fig.add_subplot(111)
    ax.set_facecolor(_PAL["bg"])
    ax.set_axis_off()
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.75, 1.15)
    ax.set_aspect("equal")

    p = 1.0 if nucleated else 0.0
    R = 1.0
    ax_top = 2.15                       # opening half-angle of the shell arc
    neck = R * (1.0 - 0.86 * p)
    drop = R * 0.6 * p

    # outer + inner sheets: an arc that pinches into a throat at the bottom
    for radius, color, ang in ((R, _PAL["sphere"], ax_top),
                               (R - 0.10, _PAL["object"], ax_top - 0.15)):
        a = np.linspace(-math.pi / 2 - ang, -math.pi / 2 + ang, 200)
        ax.plot(radius * np.cos(a), radius * np.sin(a), color=color, lw=2.4,
                solid_capstyle="round")
        # descending lips → the neck
        for sgn in (-1, 1):
            x0, y0 = sgn * radius * math.sin(ang), -radius * math.cos(ang)
            mouth_y = -R - drop
            t = np.linspace(0, 1, 60)
            bx = (1 - t) ** 2 * x0 + 2 * (1 - t) * t * (sgn * R * 0.45) + t ** 2 * (sgn * neck * 0.4)
            by = (1 - t) ** 2 * y0 + 2 * (1 - t) * t * (-R + 0.1) + t ** 2 * mouth_y
            ax.plot(bx, by, color=color, lw=2.2)

    if p > 0.15:                        # the throat mouth ring
        from matplotlib.patches import Ellipse
        mouth_y = -R - drop
        ax.add_patch(Ellipse((0, mouth_y), neck * 0.8, neck * 0.26,
                             fill=True, fc=_PAL["object"], ec=_PAL["object"],
                             alpha=0.35, lw=2.0))
        ax.text(0, mouth_y - 0.28, "throat · RP² cross-cap", ha="center",
                color=_PAL["object"], fontsize=9, family="monospace")
    else:                               # incoming fronts converging (disperse)
        for sgn in (-1, 1):
            ax.annotate("", xy=(sgn * 0.14, -R + 0.12), xytext=(sgn * 0.6, -R + 0.75),
                        arrowprops=dict(arrowstyle="->", color=_PAL["caustic"], lw=2.0))
        ax.text(-0.62, -R + 0.85, "outer", color=_PAL["sphere"], fontsize=9,
                family="monospace", ha="center")
        ax.text(0.62, -R + 0.85, "inner", color=_PAL["object"], fontsize=9,
                family="monospace", ha="center")

    ax.text(0, -R - drop - (0.5 if p > 0.15 else 0.12), "antipodal focus  θ = π",
            ha="center", color=_PAL["dim"], fontsize=8.5, family="monospace")
    state = "above threshold — throat nucleates" if nucleated else \
        "below threshold — focus disperses, field revives"
    ax.set_title(state, color=(_PAL["persist"] if nucleated else _PAL["sphere"]))
    return dict(ax=ax)


def plot_object_panel(nucleated: bool = True, *, fig: Optional[Figure] = None) -> Figure:
    """Standalone wrapper for :func:`draw_focus_object`."""
    fig = fig if fig is not None else plt.figure(figsize=(5.2, 5.2))
    draw_focus_object(fig, nucleated=nucleated)
    return fig


# ════════════════════════════════════════════════════════════════════════════
# Live animation — plane vs sphere, one clock
# ════════════════════════════════════════════════════════════════════════════
def run_animation(
    n_steps: int = 1200,
    interval_ms: int = 20,
    show: bool = True,
):
    """Animate the plane and sphere sims side by side on one clock.

    Both start from an identical cap pulse and run at the same speed; the
    peak-amplitude strip charts accumulate over one great-circle period.

    Parameters
    ----------
    n_steps
        Number of animation frames to schedule.
    interval_ms
        Frame interval in milliseconds.
    show
        If ``True``, call :func:`plt.show`.

    Returns
    -------
    (matplotlib.animation.FuncAnimation, (PlaneWaveSim, SphereWaveSim))
    """
    from matplotlib.animation import FuncAnimation

    plane = PlaneWaveSim()
    sphere = SphereWaveSim()
    hist: Dict[str, list] = {"t": [], "plane": [], "sphere": []}
    fig = plt.figure(figsize=(11, 6.2))
    dt_frame = 0.02

    def _frame(_i: int):
        plane.advance_to(min(plane.t + dt_frame, PERIOD))
        sphere.advance_to(sphere.t + dt_frame)
        if sphere.t > PERIOD:
            plane.reset(); sphere.reset(); hist["t"].clear()
            hist["plane"].clear(); hist["sphere"].clear()
        hist["t"].append(sphere.t)
        hist["plane"].append(plane.peak()[0])
        hist["sphere"].append(sphere.peak()[0])
        draw_contrast(fig, plane, sphere, history=hist)
        return []

    anim = FuncAnimation(fig, _frame, frames=n_steps, interval=interval_ms, blit=False)
    if show:
        plt.show()
    return anim, (plane, sphere)
