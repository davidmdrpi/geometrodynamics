#!/usr/bin/env python3
"""
Geometrodynamic QED — v41: the warped sphere, restored
======================================================

The archive's original picture, put back — and this time the warp is solved.

``archive/geometrodynamics_v39.py`` rendered **one continuous S²** whose radius
carried the field, nested between two fixed shells like Russian dolls.  That
picture is what made the BAM intuition watchable, and #242/#243 replaced it
with maps, strips and sections: correct, measurable, and no longer the thing
you could see.  This restores the geometry first; projections and slices come
after, not instead.

What is on screen
─────────────────
* the **middle surface** at ``R_mid = 1.00``, closed — its own poles, its seam
  matched to machine precision.  Nothing is cut out of it, so a pulse on it
  sweeps every point exactly once and *fills its own void*;
* the **two dolls** at ``R_inner = 0.74`` and ``R_outer = 1.26``, the same
  vacuole ``radial_caustic`` puts the ring caustic on.  They brighten as the
  surface comes toward them;
* the radius itself, displaced by the **solved** field — red pushed out toward
  the outer doll, blue pulled in toward the inner one;
* the leading **crest** (red) and trailing **trough** (blue) as rings *of* the
  surface, at their measured geodesic distances rather than at ``t``.

What v39 did differently
────────────────────────
v39 displaced the radius with a prescribed ``sech²`` envelope plus a mound
grown on a clock.  Here the displacement is
:class:`~geometrodynamics.viz.throat_wavefront.BareSphereSim`, an actual wave
solve, so the deformation at the antipode appears because the wave focuses
there.  And it appears **inward**: the focus arrives inverted and pulls the
surface toward the inner doll rather than pushing a spike out of it.

Honest labelling
────────────────
``r = R_mid + Δ·tanh(g·u/u_ref)`` is a *display* of the field as a radial
displacement, not backreaction — the wave does not feel the surface it is
deforming.  ``tanh`` is strictly increasing, so every sign and every ordering
of amplitudes survives it; ratios do not, and the bound is what keeps the
surface strictly between the two dolls.

Usage
─────
    python scripts/geometrodynamics_v41_warped_sphere.py             # animate
    python scripts/geometrodynamics_v41_warped_sphere.py --still out.png
    python scripts/geometrodynamics_v41_warped_sphere.py --save out.gif
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from geometrodynamics.viz.warped_sphere import (
    ANTIPODAL_TIME,
    RETURN_TIME,
    NestedShells,
    WarpedSphere,
)

FPS = 14
N_THETA, N_PHI = 121, 181

_SHELL_CMAP = LinearSegmentedColormap.from_list(
    "shell", ["#2a55c8", "#6a86d2", "#9fa3ab", "#b4b4b4", "#d2907c",
              "#e0553c", "#cc2412"])

_PAL = {
    "bg": "#000000",
    "text": "#e8ecf4",
    "dim": "#8892a4",
    "crest": "#e0483c",
    "trough": "#3f6ff0",
    "outer": "#7a2020",
    "inner": "#1b3168",
    "outer_hot": "#ff4444",
    "inner_hot": "#4488ff",
}


# ════════════════════════════════════════════════════════════════════════════
class WarpedSphereFigure:
    """The single 3-D panel: two dolls, one warped surface, two rings."""

    def __init__(self, surface: Optional[WarpedSphere] = None,
                 figsize=(8.6, 8.6), t_end: float = RETURN_TIME) -> None:
        self.s = surface or WarpedSphere(n_theta=N_THETA, n_phi=N_PHI)
        self.t_end = float(t_end)
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        self.ax = self.fig.add_subplot(111, projection="3d",
                                       facecolor=_PAL["bg"])
        self.fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
        self.doll = self.s.shells.unit_sphere(31, 49)

    # ── pieces ──────────────────────────────────────────────────────────────
    def _draw_doll(self, radius: float, colour: str, hot: str,
                   heat: float, stride: int = 5) -> None:
        h = float(np.clip(heat, 0.0, 1.0))
        c = ((1.0 - h) * np.array(matplotlib.colors.to_rgb(colour))
             + h * np.array(matplotlib.colors.to_rgb(hot)))
        X, Y, Z = self.doll
        self.ax.plot_wireframe(radius * X, radius * Y, radius * Z,
                               color=tuple(c), alpha=0.085 + 0.16 * h,
                               rstride=stride, cstride=stride,
                               linewidth=0.22 + 0.7 * h)

    def _draw_surface(self) -> None:
        X, Y, Z = self.s.mesh()
        c = np.tanh(self.s.gain * self.s.normalised_field())
        self.ax.plot_surface(X, Y, Z,
                             facecolors=_SHELL_CMAP(0.5 * (1.0 + c)),
                             rstride=1, cstride=1, antialiased=True,
                             shade=True, alpha=0.90, linewidth=0)
        self.ax.plot_wireframe(X, Y, Z, color="white", alpha=0.05,
                               rstride=10, cstride=10, linewidth=0.2)

    def _draw_rings(self) -> None:
        f = self.s.front_distances()
        ref = self.s.reference_amplitude
        for key, colour in (("crest", _PAL["crest"]), ("trough", _PAL["trough"])):
            d = f[f"{key}_distance"]
            a = abs(f[f"{key}_amplitude"]) / ref
            if a < 0.02 or d < 1e-3 or d > math.pi - 1e-3:
                continue
            ring = self.s.geodesic_circle(d)
            self.ax.plot(ring[:, 0], ring[:, 1], ring[:, 2], color=colour,
                         alpha=float(np.clip(0.35 + 0.65 * a, 0.0, 0.95)),
                         linewidth=1.0 + 2.6 * a, zorder=60)

    def _draw_poles(self) -> None:
        for d, colour, label in ((0.0, "#ffffff", "source"),
                                 (math.pi, "#ffd66b", "antipode")):
            p = self.s.marker(d)
            self.ax.plot([p[0]], [p[1]], [p[2]], "o", color=colour, ms=4.0,
                         alpha=0.9, zorder=70)
            self.ax.text(p[0] * 1.12, p[1] * 1.12, p[2] * 1.12, label,
                         color=colour, fontsize=7.5, alpha=0.75)

    def _draw_text(self) -> None:
        ex = self.s.excursion()
        f = self.s.front_distances()
        ax = self.ax
        ax.text2D(0.035, 0.965,
                  "Geometrodynamic QED  v41  —  one continuous S², warped by "
                  "the wave it carries",
                  transform=ax.transAxes, color=_PAL["text"], fontsize=10.5,
                  fontweight="bold")
        ax.text2D(0.035, 0.932,
                  f"t = {self.s.t:5.3f}     antipodal focus at t = π = "
                  f"{ANTIPODAL_TIME:.3f}     home at t = 2π = {RETURN_TIME:.3f}",
                  transform=ax.transAxes, color=_PAL["dim"], fontsize=9,
                  family="monospace")
        lines = [
            (f"R_outer  {self.s.shells.r_outer:.2f}   surface is "
             f"{100 * ex['outward_fraction']:5.1f}% of the way out",
             _PAL["crest"]),
            (f"R_inner  {self.s.shells.r_inner:.2f}   surface is "
             f"{100 * ex['inward_fraction']:5.1f}% of the way in",
             _PAL["trough"]),
            (f"crest   at d = {f['crest_distance']:.3f}", _PAL["dim"]),
            (f"trough  at d = {f['trough_distance']:.3f}", _PAL["dim"]),
        ]
        y = 0.115
        for txt, col in lines:
            ax.text2D(0.035, y, txt, transform=ax.transAxes, color=col,
                      fontsize=8.5, family="monospace")
            y -= 0.026
        ax.text2D(0.62, 0.055,
                  "radius = R_mid + Δ·tanh(g·u/u_ref)\n"
                  "the field is solved; the displacement\n"
                  "is a display of it, not backreaction",
                  transform=ax.transAxes, color=_PAL["dim"], fontsize=7.2,
                  family="monospace", va="bottom")

    def azimuth(self, t: float) -> float:
        """Follow the wave: the camera faces the source at ``t = 0`` and the
        antipode at ``t = π``, kept 40° off so the warp reads on the limb."""
        return 40.0 + 360.0 * float(t) / self.t_end

    # ── one frame ───────────────────────────────────────────────────────────
    def draw(self, t: float, azim: Optional[float] = None) -> None:
        if t < self.s.t:
            self.s.reset()
        self.s.advance_to(t)
        ax = self.ax
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        ex = self.s.excursion()
        self._draw_doll(self.s.shells.r_outer, _PAL["outer"], _PAL["outer_hot"],
                        ex["outward_fraction"])
        self._draw_doll(self.s.shells.r_inner, _PAL["inner"], _PAL["inner_hot"],
                        ex["inward_fraction"])
        self._draw_surface()
        self._draw_rings()
        self._draw_poles()
        self._draw_text()

        lim = 1.02 * self.s.shells.r_outer
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        ax.set_box_aspect((1, 1, 1), zoom=1.48)
        ax.set_axis_off()
        ax.view_init(elev=13.0, azim=self.azimuth(self.s.t)
                     if azim is None else azim)


# ════════════════════════════════════════════════════════════════════════════
def still(path: str, n: int = 4) -> None:
    """A contact sheet through the run: launch, ring, focus, return."""
    surface = WarpedSphere(n_theta=N_THETA, n_phi=N_PHI)
    times = [0.45, 0.5 * ANTIPODAL_TIME, ANTIPODAL_TIME, 1.55 * ANTIPODAL_TIME]
    times = times[:n]
    fig = plt.figure(figsize=(3.9 * len(times), 4.3), facecolor=_PAL["bg"])
    fig.subplots_adjust(left=0.005, right=0.995, top=0.90, bottom=0.005,
                        wspace=0.0)
    for k, t in enumerate(times):
        ax = fig.add_subplot(1, len(times), k + 1, projection="3d",
                             facecolor=_PAL["bg"])
        holder = WarpedSphereFigure.__new__(WarpedSphereFigure)
        holder.s = surface
        holder.fig = fig
        holder.ax = ax
        holder.doll = surface.shells.unit_sphere(31, 49)
        holder.t_end = RETURN_TIME
        surface.reset()
        surface.advance_to(t)
        ex = surface.excursion()
        holder._draw_doll(surface.shells.r_outer, _PAL["outer"],
                          _PAL["outer_hot"], ex["outward_fraction"])
        holder._draw_doll(surface.shells.r_inner, _PAL["inner"],
                          _PAL["inner_hot"], ex["inward_fraction"])
        holder._draw_surface()
        holder._draw_rings()
        lim = 1.02 * surface.shells.r_outer
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        ax.set_box_aspect((1, 1, 1), zoom=1.48)
        ax.set_axis_off()
        ax.view_init(elev=13.0, azim=40.0 + 360.0 * t / RETURN_TIME)
        ax.set_title(f"t = {t:.3f}" + ("   antipodal focus"
                                       if abs(t - ANTIPODAL_TIME) < 1e-9 else ""),
                     color=_PAL["text"], fontsize=10, y=0.99)
    fig.suptitle("one continuous surface, warped by a solved wave — the focus "
                 "pulls it in toward the inner doll",
                 color=_PAL["text"], fontsize=12, y=0.985)
    fig.savefig(path, dpi=120, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def animate(save: Optional[str] = None, frames: int = 120):
    from matplotlib import animation

    holder = WarpedSphereFigure()
    holder.s.reset()

    def update(i: int):
        t = (i + 1) * holder.t_end / frames
        holder.draw(t)
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
                    help="write a contact sheet through the run and exit")
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
