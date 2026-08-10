#!/usr/bin/env python3
"""
Geometrodynamic QED — v43: spin 0 against spin 2, on one clock
==============================================================

v41 restored the archive's picture: one continuous S² whose *radius* carries a
solved field.  But that field is a **scalar**, displayed extrinsically as a
height.  A metric perturbation is not that kind of object — `h_ab` is symmetric
and trace-free, spin 2, and it does not push a surface outward at all.  It
shears it.

Both panels run the same pulse width on the same clock, from the same point.

left · spin 0
    the scalar, as a radial height.  The surface breathes: local area changes,
    the amplitude is free to peak anywhere, and at `t = π` it piles up **on**
    the antipode.

right · spin 2
    the tensor, as tidal ellipses — what a small ring of test particles
    becomes.  Red is stretched along the geodesic radial direction and squeezed
    across it; blue is the other way round.  The sphere itself is **not**
    deformed, because a trace-free perturbation does not change area: to first
    order every ellipse has exactly the area of the circle it came from.

What to watch for
─────────────────
* the tensor field **vanishes at both poles** — spin weight 2 forbids a
  well-defined amplitude where the frame degenerates.  So at the focus it is a
  *ring* around the antipode, never a peak on it, and the smallest source it
  admits is a ring too;
* the ellipse pattern repeats every **180°**, not 360°.  That is the spin
  weight, visible directly;
* through the caustic the polarisation does not simply invert.  One passage is
  a **quarter turn** — the Gouy shift, Maslov index 1, measured here as a
  correlation of 0.82 with the Hilbert transform of the inbound waveform
  against −0.35 with its inverse.  The stretch and compression axes really do
  swap, but it takes **two** focal passages: at `t = 2π`, after the antipode
  and back home, the field is minus what it started as (correlation 0.997).

Honest scope
────────────
A spin-2 field on a **fixed** S², not a solution of the 4-D linearised Einstein
equations — general relativity in 2+1 dimensions has no propagating tensor
modes at all.  What is faithful here is the polarisation structure: two
components, spin weight 2, `ℓ ≥ 2` only, area-preserving shear, and the
behaviour at a focus.  The ellipse sizes are a display gain; their *shapes*
are the solved field.

Usage
─────
    python scripts/geometrodynamics_v43_tidal_sphere.py             # animate
    python scripts/geometrodynamics_v43_tidal_sphere.py --still out.png
    python scripts/geometrodynamics_v43_tidal_sphere.py --save out.gif
"""

from __future__ import annotations

import argparse
import math
from typing import List, Optional, Tuple

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from geometrodynamics.viz.spin2_tidal import (
    ANTIPODAL_TIME,
    RETURN_TIME,
    Spin2WaveSim,
    TidalField,
)
from geometrodynamics.viz.warped_sphere import NestedShells, WarpedSphere

FPS = 14
PULSE_WIDTH = 0.18

_SHELL_CMAP = LinearSegmentedColormap.from_list(
    "shell", ["#2a55c8", "#6a86d2", "#9fa3ab", "#b4b4b4", "#d2907c",
              "#e0553c", "#cc2412"])

_PAL = {
    "bg": "#000000",
    "text": "#e8ecf4",
    "dim": "#8892a4",
    "radial": "#e0483c",      # stretched along ê_d
    "transverse": "#3f6ff0",  # stretched across it
    "ghost": "#3a3f4a",
}


def _view_direction(elev: float, azim: float) -> np.ndarray:
    e, a = math.radians(elev), math.radians(azim)
    return np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a),
                     math.sin(e)])


# ════════════════════════════════════════════════════════════════════════════
class TidalFigure:
    """Two panels, one clock: the scalar as height, the tensor as shear."""

    def __init__(self, figsize=(13.4, 7.2), t_end: float = RETURN_TIME,
                 n_rings: int = 11, n_azimuth: int = 14,
                 ellipse_size: float = 0.135) -> None:
        self.t_end = float(t_end)
        self.scalar = WarpedSphere(n_theta=101, n_phi=151,
                                   pulse_width=PULSE_WIDTH)
        self.tensor = TidalField(sim=Spin2WaveSim(n=1200,
                                                  pulse_width=PULSE_WIDTH))
        self.size = float(ellipse_size)
        self.samples: List[Tuple[float, float]] = [
            (d, psi)
            for d in np.linspace(0.16, math.pi - 0.16, n_rings)
            for psi in np.linspace(0.0, 2.0 * math.pi, n_azimuth, endpoint=False)
        ]
        self.gain = self._calibrate()

        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        self.ax_s = self.fig.add_subplot(121, projection="3d",
                                         facecolor=_PAL["bg"])
        self.ax_t = self.fig.add_subplot(122, projection="3d",
                                         facecolor=_PAL["bg"])
        self.fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0,
                                 wspace=0.0)
        self.unit = NestedShells().unit_sphere(61, 91)

    def _calibrate(self, samples: int = 240) -> float:
        """Display gain from the run's own peak strain — a measured number."""
        self.tensor.reset()
        peak = 0.0
        for i in range(samples):
            self.tensor.advance_to((i + 1) * self.t_end / samples)
            peak = max(peak, abs(self.tensor.sim.peak()[1]))
        self.tensor.reset()
        self.peak_strain = peak
        return 0.85 / peak if peak > 0.0 else 1.0

    # ── panels ──────────────────────────────────────────────────────────────
    def _draw_scalar(self, ax) -> None:
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        X, Y, Z = self.scalar.mesh()
        c = np.tanh(self.scalar.gain * self.scalar.normalised_field())
        ax.plot_surface(X, Y, Z, facecolors=_SHELL_CMAP(0.5 * (1.0 + c)),
                        rstride=1, cstride=1, antialiased=True, shade=True,
                        alpha=0.95, linewidth=0)
        ax.plot_wireframe(X, Y, Z, color="white", alpha=0.05,
                          rstride=10, cstride=10, linewidth=0.2)

    def _draw_tensor(self, ax, view: np.ndarray) -> None:
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        X, Y, Z = self.unit
        ax.plot_surface(X, Y, Z, color="#15171c", rstride=2, cstride=2,
                        antialiased=True, shade=True, alpha=0.98, linewidth=0)
        ax.plot_wireframe(X, Y, Z, color=_PAL["ghost"], alpha=0.30,
                          rstride=10, cstride=10, linewidth=0.25)
        scale = self.peak_strain or 1.0
        for d, psi in self.samples:
            h = self.tensor.ring_amplitude(d)
            a = abs(h) / scale
            if a < 0.012:
                continue
            p = self.tensor.point(d, psi)
            if float(np.dot(p, view)) < 0.06:      # far side — cull
                continue
            ring = self.tensor.ellipse(
                d, psi, eps=self.gain, n=40,
                size=self.size * max(math.sin(d), 0.22), lift=1.012)
            ax.plot(ring[:, 0], ring[:, 1], ring[:, 2],
                    color=_PAL["radial"] if h > 0 else _PAL["transverse"],
                    alpha=float(np.clip(0.30 + 0.65 * a, 0.0, 0.95)),
                    linewidth=0.8 + 1.5 * a, zorder=40)

    def _frame_axes(self, ax, elev: float, azim: float, zoom: float) -> None:
        lim = 1.02 * (self.scalar.shells.r_outer if ax is self.ax_s else 1.13)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        ax.set_box_aspect((1, 1, 1), zoom=zoom)
        ax.set_axis_off()
        ax.view_init(elev=elev, azim=azim)

    def _text(self) -> None:
        f = self.fig
        for artist in list(f.texts):
            artist.remove()
        peak_d, peak_h = self.tensor.sim.peak()
        ex = self.scalar.excursion()
        f.suptitle("Geometrodynamic QED  v43  —  the same wave as a scalar and "
                   "as a tensor, on one clock",
                   color=_PAL["text"], fontsize=12.5, fontweight="bold",
                   y=0.972)
        f.text(0.5, 0.930,
               f"t = {self.tensor.t:5.3f}      antipodal focus at t = π = "
               f"{ANTIPODAL_TIME:.3f}      home at t = 2π = {RETURN_TIME:.3f}",
               color=_PAL["dim"], fontsize=9.5, family="monospace",
               ha="center")
        f.text(0.25, 0.893, "spin 0  ·  scalar u, as radial height",
               color=_PAL["text"], fontsize=10.5, ha="center")
        f.text(0.75, 0.893, "spin 2  ·  tensor h_ab, as tidal ellipses",
               color=_PAL["text"], fontsize=10.5, ha="center")
        f.text(0.25, 0.093,
               f"area changes: the surface breathes\n"
               f"out {100 * ex['outward_fraction']:5.1f}%   "
               f"in {100 * ex['inward_fraction']:5.1f}%  of the gap\n"
               f"free to peak at the antipode",
               color=_PAL["dim"], fontsize=8.4, family="monospace",
               ha="center", va="top")
        ring = ANTIPODAL_TIME - peak_d
        where = (f"ring of radius {ring:.3f} about the antipode"
                 if peak_d > 0.5 * ANTIPODAL_TIME
                 else f"ring of radius {peak_d:.3f} about the source")
        f.text(0.75, 0.093,
               f"area preserved: pure shear, no breathing\n"
               f"peak |h| at d = {peak_d:.3f} — a {where}\n"
               f"h = 0 at both poles: spin weight 2 forbids it",
               color=_PAL["dim"], fontsize=8.4, family="monospace",
               ha="center", va="top")
        f.text(0.75, 0.022,
               "red: stretched along ê_d   ·   blue: stretched across it",
               color=_PAL["dim"], fontsize=8.0, family="monospace",
               ha="center")

    # ── one frame ───────────────────────────────────────────────────────────
    def draw(self, t: float, azim: Optional[float] = None) -> None:
        if t < self.tensor.t:
            self.tensor.reset()
        if t < self.scalar.t:
            self.scalar.reset()
        self.tensor.advance_to(t)
        self.scalar.advance_to(t)
        az = (40.0 + 360.0 * t / self.t_end) if azim is None else azim
        elev = 13.0
        self._draw_scalar(self.ax_s)
        self._draw_tensor(self.ax_t, _view_direction(elev, az))
        self._frame_axes(self.ax_s, elev, az, 1.42)
        self._frame_axes(self.ax_t, elev, az, 1.42)
        self._text()


# ════════════════════════════════════════════════════════════════════════════
def still(path: str) -> None:
    """Three moments: outbound, at the focus, and after it."""
    holder = TidalFigure(figsize=(13.4, 7.2))
    times = [0.9, ANTIPODAL_TIME, 1.35 * ANTIPODAL_TIME]
    frames = []
    for t in times:
        holder.draw(t)
        holder.fig.canvas.draw()
        frames.append(np.asarray(holder.fig.canvas.buffer_rgba()).copy())
    plt.close(holder.fig)

    fig, axes = plt.subplots(len(frames), 1,
                             figsize=(11.0, 5.9 * len(frames)),
                             facecolor=_PAL["bg"])
    for ax, img, t in zip(np.atleast_1d(axes), frames, times):
        ax.imshow(img)
        ax.set_axis_off()
        ax.set_title(f"t = {t:.3f}" + ("   antipodal focus"
                                       if abs(t - ANTIPODAL_TIME) < 1e-9 else ""),
                     color=_PAL["text"], fontsize=11)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.01, hspace=0.04)
    fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def animate(save: Optional[str] = None, frames: int = 120):
    from matplotlib import animation

    holder = TidalFigure()

    def update(i: int):
        holder.draw((i + 1) * holder.t_end / frames)
        return []

    anim = animation.FuncAnimation(holder.fig, update, frames=frames,
                                   interval=1000 // FPS, blit=False)
    if save:
        anim.save(save, fps=FPS, dpi=100,
                  savefig_kwargs={"facecolor": _PAL["bg"]})
        print(f"wrote {save}")
    else:
        plt.show()
    return anim


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG",
                    help="write a three-moment sheet and exit")
    ap.add_argument("--save", metavar="FILE",
                    help="render the animation to a file")
    ap.add_argument("--frames", type=int, default=120)
    a = ap.parse_args(argv)

    if a.still:
        matplotlib.use("Agg")
        still(a.still)
        return 0
    if a.save:
        matplotlib.use("Agg")
    animate(save=a.save, frames=a.frames)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
