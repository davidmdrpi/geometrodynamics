#!/usr/bin/env python3
"""
Geometrodynamic QED — v44: the tensor wave as a continuous embedded surface
===========================================================================

v43 drew the spin-2 field as discrete tangent-space ellipses.  Faithful, but
flat: it never touches the embedding, so it says nothing about how tidal shear
meets a bulk.  This projects `h_ab` into a **continuous** surface deformation
`r(θ, φ, t)` in ℝ³ instead, and the projection is forced rather than chosen.

Why a height field alone cannot do it
─────────────────────────────────────
For `X = r n̂` the induced metric is `g_ab = r² ĝ_ab + ∂_a r ∂_b r`, whose
gradient term is *second* order.  So at first order a radial deformation gives
`δg_ab = 2ρ ĝ_ab` — purely conformal.  Shape carries the trace and nothing
else, which is exactly why the ellipses existed.

The displacement that does
──────────────────────────
Add the tangential part.  With `ξ = ∇W`,

    X = (R + ερ) n̂ + ε ξ ,      δg_ab = 2ρ ĝ_ab + 2∇₍ₐξ_b₎ ,

and demanding tracelessness fixes the radial part completely: `ρ = −½ΔW`,
while `h_ab = [2∇₍ₐ∇_b₎W]^TF`.  **One potential carries both** — the shear you
measure is its trace-free Hessian, the shape you see is minus half its
Laplacian.  The shear lives in how the material *slides*; the shape carries
what is left over.

What is on screen
─────────────────
* the **deformed surface**, continuous, its radius carrying `ρ`;
* its **material lattice** — rings and meridians of test particles, evenly
  spaced before the wave arrives.  Where they bunch and spread is the tidal
  stretch and squeeze, continuously, with nothing sampled;
* colour is the **shear** `h` itself: red stretched along `ê_d`, blue across it;
* **principal-axis bars** at sparse points — two short tangent vectors along the
  eigenvectors of `h_ab`, red for the positive eigenvalue and blue for the
  negative, each of length proportional to `|λ|`.  They are always the same
  length as each other, because trace-free in two dimensions means
  `λ± = ±|h|`, and the stretch axis *swaps* between `ê_d` and `ê_ψ` wherever
  `h` changes sign;
* the two **shells** of the vacuole, so the reach into the bulk is legible;
* the profiles panel: one potential, its two parts.

What comes out of it
────────────────────
`ℓ = 2` gives `ρ = P₂(cos d)` to `1.0e-07` — the quadrupole tidal wave *is* the
prolate–oblate deformation of the textbook picture, derived here rather than
drawn.  And the surface keeps its area to first order, because that is what
trace-free means.

Honest scope
────────────
Still a spin-2 field on a fixed background: this is a faithful *representation*
of `h_ab` in the embedding, not backreaction.  What it buys is an extrinsic
amplitude, so the wave can be watched reaching into a bulk.  The gain `ε` is a
display choice; the shape, at any gain, is the solved field.

Usage
─────
    python scripts/geometrodynamics_v44_embedded_wave.py             # animate
    python scripts/geometrodynamics_v44_embedded_wave.py --still out.png
    python scripts/geometrodynamics_v44_embedded_wave.py --save out.gif
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from geometrodynamics.viz.embedded_wave import EmbeddedTidalSurface
from geometrodynamics.viz.spin2_tidal import (
    ANTIPODAL_TIME,
    RETURN_TIME,
    Spin2WaveSim,
)

FPS = 14
PULSE_WIDTH = 0.18

_SHEAR_CMAP = LinearSegmentedColormap.from_list(
    "shear", ["#2a55c8", "#6a86d2", "#9fa3ab", "#b4b4b4", "#d2907c",
              "#e0553c", "#cc2412"])

_PAL = {
    "bg": "#000000",
    "panel": "#07090f",
    "text": "#e8ecf4",
    "dim": "#8892a4",
    "radial": "#e0483c",
    "transverse": "#3f6ff0",
    "lattice": "#dfe6f2",
    "outer": "#7a2020",
    "inner": "#1b3168",
    "outer_hot": "#ff4444",
    "inner_hot": "#4488ff",
    "border": "#232a38",
}


class EmbeddedFigure:
    """One continuous surface, its material lattice, and the two profiles."""

    def __init__(self, figsize=(12.6, 8.0), t_end: float = RETURN_TIME) -> None:
        self.t_end = float(t_end)
        self.s = EmbeddedTidalSurface(
            sim=Spin2WaveSim(n=1200, pulse_width=PULSE_WIDTH),
            n_theta=121, n_phi=181)
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(2, 2, width_ratios=[2.5, 1.0],
                                   height_ratios=[1.35, 1.0],
                                   left=0.005, right=0.975, top=0.905,
                                   bottom=0.035, wspace=0.06, hspace=0.30)
        self.ax3d = self.fig.add_subplot(gs[:, 0], projection="3d",
                                         facecolor=_PAL["bg"])
        self.ax_p = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_t = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["bg"])
        self.doll = self.s.shells.unit_sphere(31, 49)
        self.scale = self.s.peak_radial_profile or 1.0
        self.cross_samples = [
            (d, psi)
            for d in np.linspace(0.18, math.pi - 0.18, 14)
            for psi in np.linspace(0.0, 2.0 * math.pi, 13, endpoint=False)
        ]
        self.scale_shear = self._calibrate_shear()

    def _calibrate_shear(self, samples: int = 160) -> float:
        """The run's own peak |λ|, so bar lengths mean the same thing all through."""
        self.s.reset()
        peak = 0.0
        for i in range(samples):
            self.s.advance_to((i + 1) * self.t_end / samples)
            peak = max(peak, float(np.max(np.abs(self.s.profiles()["shear"]))))
        self.s.reset()
        return peak

    # ── pieces ──────────────────────────────────────────────────────────────
    def _doll(self, radius: float, colour: str, hot: str, heat: float) -> None:
        h = float(np.clip(heat, 0.0, 1.0))
        c = ((1.0 - h) * np.array(matplotlib.colors.to_rgb(colour))
             + h * np.array(matplotlib.colors.to_rgb(hot)))
        X, Y, Z = self.doll
        self.ax3d.plot_wireframe(radius * X, radius * Y, radius * Z,
                                 color=tuple(c), alpha=0.075 + 0.16 * h,
                                 rstride=5, cstride=5,
                                 linewidth=0.22 + 0.7 * h)

    def _surface(self) -> None:
        X, Y, Z = self.s.mesh()
        shear = self.s.shear_on_mesh()
        peak = float(np.max(np.abs(self.s.profiles()["shear"]))) or 1.0
        c = np.tanh(2.2 * shear / peak)
        self.ax3d.plot_surface(X, Y, Z, facecolors=_SHEAR_CMAP(0.5 * (1.0 + c)),
                               rstride=1, cstride=1, antialiased=True,
                               shade=True, alpha=0.93, linewidth=0)

    def _lattice(self) -> None:
        for curve in self.s.material_latitudes(n_rings=13, n_points=220):
            self.ax3d.plot(curve[:, 0], curve[:, 1], curve[:, 2],
                           color=_PAL["lattice"], alpha=0.34, linewidth=0.6,
                           zorder=30)
        for curve in self.s.material_meridians(n_lines=10, n_points=240):
            self.ax3d.plot(curve[:, 0], curve[:, 1], curve[:, 2],
                           color=_PAL["lattice"], alpha=0.13, linewidth=0.4,
                           zorder=30)

    def _crosses(self, view: np.ndarray) -> None:
        """Principal-axis bars at sparse points: ↔ stretch, ↕ squeeze.

        Aligned with the eigenvectors of ``h_ab`` and scaled by ``|λ|``.  The
        two bars are always the same length — trace-free means ``λ± = ±|h|``,
        so an asymmetric cross would be a bug rather than a feature.
        """
        ref = float(np.max(np.abs(self.s.profiles()["shear"]))) or 1.0
        run_peak = self.scale_shear or ref
        for d0, psi in self.cross_samples:
            lam = abs(float(self.s.principal_axes(np.array([d0]))["lambda_plus"][0]))
            if lam < 0.02 * run_peak:
                continue
            centre = self.s.positions(np.array([[d0]]), np.array([[psi]]))[0, 0]
            if float(centre @ view) < 0.06 * float(np.linalg.norm(centre)):
                continue                       # far side
            c = self.s.principal_cross(d0, psi, size=0.125, reference=run_peak,
                                       lift=0.014)
            a = float(np.clip(lam / run_peak, 0.0, 1.0))
            for key, colour in (("stretch", _PAL["radial"]),
                                ("squeeze", _PAL["transverse"])):
                bar = c[key]
                self.ax3d.plot(bar[:, 0], bar[:, 1], bar[:, 2], color=colour,
                               alpha=float(np.clip(0.45 + 0.5 * a, 0.0, 1.0)),
                               linewidth=1.0 + 1.9 * a,
                               solid_capstyle="round", zorder=50)

    def _profiles(self) -> None:
        ax = self.ax_p
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        p = self.s.profiles()
        d = p["d"]
        ax.axhline(0.0, color=_PAL["border"], lw=0.7)
        ax.plot(d, p["shear"] / self.scale, color=_PAL["radial"], lw=1.4,
                label="shear  h  (trace-free)")
        ax.plot(d, p["radial"] / self.scale, color="#e8ecf4", lw=1.4,
                label="radial  ρ = −½ΔW")
        ax.plot(d, p["slide"] / self.scale, color=_PAL["transverse"], lw=1.1,
                alpha=0.85, label="slide  ξ = ∇W")
        ax.set_xlim(0.0, math.pi)
        ax.set_ylim(-1.55, 1.55)
        ax.set_xticks([0.0, math.pi / 2, math.pi])
        ax.set_xticklabels(["source", "π/2", "antipode"], fontsize=7.5)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["border"])
        leg = ax.legend(loc="upper right", fontsize=6.6, frameon=False)
        for txt in leg.get_texts():
            txt.set_color(_PAL["dim"])
        ax.set_title("one potential, two parts", color=_PAL["dim"],
                     fontsize=8.5, pad=4)

    def _readout(self) -> None:
        ax = self.ax_t
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        ax.set_axis_off()
        ex = self.s.excursion()
        p = self.s.profiles()
        m = self.s.induced_metric_perturbation(
            float(np.clip(ex["peak_distance"], 0.15, math.pi - 0.15)),
            gain=1e-4)
        peak = float(np.max(np.abs(p["shear"])))
        lines = [
            (f"t = {self.s.t:5.3f}", _PAL["text"]),
            ("", None),
            ("THE DRAWN SURFACE'S OWN METRIC", _PAL["text"]),
            (f"  h+ measured   {m['h_plus_measured']:+.6f}", _PAL["dim"]),
            (f"  h+ solved     {m['h_plus_solved']:+.6f}", _PAL["dim"]),
            (f"  trace         {m['trace']:+.2e}", _PAL["dim"]),
            ("", None),
            ("REACH INTO THE BULK", _PAL["text"]),
            (f"  toward R_outer {100 * ex['outward_fraction']:5.1f}%",
             _PAL["radial"]),
            (f"  toward R_inner {100 * ex['inward_fraction']:5.1f}%",
             _PAL["transverse"]),
            ("", None),
            (f"peak shear |h| = {peak:.4f}", _PAL["dim"]),
            (f"display gain ε = {self.s.gain:.2f}", _PAL["dim"]),
        ]
        y = 0.97
        for txt, col in lines:
            if txt:
                ax.text(0.02, y, txt, color=col or _PAL["dim"],
                        family="monospace", fontsize=8.2,
                        transform=ax.transAxes, va="top")
            y -= 0.063
        ax.text(0.02, 0.005,
                "the shear is the trace-free Hessian of W,\n"
                "the shape is −½ΔW — a representation of\n"
                "h_ab in the embedding, not backreaction",
                color=_PAL["dim"], family="monospace", fontsize=7.0,
                transform=ax.transAxes, va="bottom")

    # ── one frame ───────────────────────────────────────────────────────────
    def draw(self, t: float, azim: Optional[float] = None) -> None:
        if t < self.s.t:
            self.s.reset()
        self.s.advance_to(t)
        ax = self.ax3d
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        ex = self.s.excursion()
        self._doll(self.s.shells.r_outer, _PAL["outer"], _PAL["outer_hot"],
                   ex["outward_fraction"])
        self._doll(self.s.shells.r_inner, _PAL["inner"], _PAL["inner_hot"],
                   ex["inward_fraction"])
        self._surface()
        self._lattice()
        az = ((75.0 + 360.0 * t / self.t_end) if azim is None else azim)
        e, a = math.radians(13.0), math.radians(az)
        self._crosses(np.array([math.cos(e) * math.cos(a),
                                math.cos(e) * math.sin(a), math.sin(e)]))
        lim = 1.02 * self.s.shells.r_outer
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        ax.set_box_aspect((1, 1, 1), zoom=1.5)
        ax.set_axis_off()
        ax.view_init(elev=13.0, azim=az)
        self._profiles()
        self._readout()
        for artist in list(self.fig.texts):
            artist.remove()
        self.fig.suptitle(
            "Geometrodynamic QED  v44  —  a spin-2 wave as a continuous "
            "deformation of the embedded sphere",
            color=_PAL["text"], fontsize=12.0, fontweight="bold", y=0.968)
        self.fig.text(0.5, 0.936,
                      "shape = −½ΔW    shear = trace-free Hessian of W    "
                      "bars are the eigenvectors of h_ab, length ∝ |λ|",
                      color=_PAL["dim"], fontsize=9.0, family="monospace",
                      ha="center")


# ════════════════════════════════════════════════════════════════════════════
def still(path: str) -> None:
    holder = EmbeddedFigure()
    times = [0.9, ANTIPODAL_TIME, 1.4 * ANTIPODAL_TIME]
    frames = []
    for t in times:
        holder.draw(t)
        holder.fig.canvas.draw()
        frames.append(np.asarray(holder.fig.canvas.buffer_rgba()).copy())
    plt.close(holder.fig)
    fig, axes = plt.subplots(len(frames), 1, figsize=(11.0, 7.0 * len(frames)),
                             facecolor=_PAL["bg"])
    for ax, img, t in zip(np.atleast_1d(axes), frames, times):
        ax.imshow(img)
        ax.set_axis_off()
        ax.set_title(f"t = {t:.3f}" + ("   antipodal focus"
                                       if abs(t - ANTIPODAL_TIME) < 1e-9 else ""),
                     color=_PAL["text"], fontsize=11)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.01, hspace=0.05)
    fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def animate(save: Optional[str] = None, frames: int = 110):
    from matplotlib import animation

    holder = EmbeddedFigure()

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
    ap.add_argument("--still", metavar="PNG")
    ap.add_argument("--save", metavar="FILE")
    ap.add_argument("--frames", type=int, default=110)
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
