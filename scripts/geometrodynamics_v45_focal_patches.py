#!/usr/bin/env python3
"""
Geometrodynamic QED — v45: what refocusing does to a trace-free deformation
===========================================================================

v44 established the projection: one potential `W` carries both the shape
`ρ = −½ΔW` and the shear `h_ab = [2∇₍ₐ∇_b₎W]^TF`, so the drawn surface *has*
the solved metric perturbation as its own geometry.  But the spin-2 identity
was still mostly numerical — it lived in the trace-free Hessian, not on screen.

This asks the next question instead.  Every point of the sphere runs its own
principal-strain history; those histories are all launched from one source and
they all arrive back at the antipode together.  **What happens to a continuous
trace-free deformation where they refocus?**

Three things on screen, and nothing else
────────────────────────────────────────
1. the **continuous deformed embedding surface**, coloured by the signed shear
   `h` — red where the stretch axis is `ê_d`, blue where it has swapped to
   `ê_ψ`, so the polarisation structure is the colour;
2. **sparse principal-strain vectors** — at scattered points, two short tangent
   bars along the eigenvectors of `h_ab`, red for `λ₊` and blue for `λ₋`, each
   of length proportional to `|λ|`;
3. **two advected material patches** — discs of labelled test particles, one at
   mid-latitude and one sitting on the focal ring, carried by the displacement
   and drawn wherever it puts them.

No algebraic panel.  The patches are the argument: they are the only objects
here that can change shape *and* be checked for size at the same time.

What the refocus does
─────────────────────
The strains arrive on a **ring, not a point**.  A spin-2 field must vanish at
both poles of its own geodesic-polar frame — it is a weight-2 object, so it has
no `ℓ = 0` or `ℓ = 1` piece and no value at the focus itself — and the measured
amplitude on the antipode is `~1e-08` of the peak while the ring around it is
the brightest thing on the sphere.  Its radius tracks the pulse width at
`≈ 0.94 × width`.

On that ring the patch distorts about **eight times harder** than the identical
patch at mid-latitude — `1.397` against `1.047` in aspect ratio at the same
display gain.  Its **area is held by the trace-free condition**, but only at
first order in `ε`: across a full return the relative swing is `1.9e-07` at
`ε = 1e-2`, where the linear statement is the whole statement.

And that is where the interesting thing happens.  Push to the display gain and
the residual — second order in `ε`, times the *local gradient* of the field —
grows to `2%` for the smooth patch and `26%` for the focal one.  Away from the
focus the gradient scale is the wavelength; on the focal ring it is the pulse
width, so the same term is an order of magnitude larger there.  The measured
exponent is `2.004`, so this is the linear description running out, not a
numerical defect: **the area law fails first, and fails hardest, exactly on the
refocus.**  Both numbers are on screen, and both vanish with `ε`.

Honest scope
────────────
* The peak effective energy density `∝ ḣ_ab ḣ^ab` amplifies by a finite `~2.3×`
  at the focus in this linear model — **and the scalar does the same** (`~2.0×`).
  The factor is not a spin-2 protection mechanism.  What belongs to the spin is
  the *node* and the *ring*, not the number.
* No singularity forms here, and none can: this is a linear field on a fixed
  round background with no backreaction and no bulk crossing rule.  Representing
  the focusing is the prerequisite for those, not a substitute.
* The gain `ε` is a display choice.  The shape at any gain is the solved field.

Usage
─────
    python scripts/geometrodynamics_v45_focal_patches.py             # animate
    python scripts/geometrodynamics_v45_focal_patches.py --still out.png
    python scripts/geometrodynamics_v45_focal_patches.py --save out.gif
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from geometrodynamics.viz.embedded_wave import EmbeddedTidalSurface, MaterialPatch
from geometrodynamics.viz.spin2_tidal import (
    ANTIPODAL_TIME,
    RETURN_TIME,
    Spin2WaveSim,
)

FPS = 14
PULSE_WIDTH = 0.18
PATCH_RADIUS = 0.12

_SHEAR_CMAP = LinearSegmentedColormap.from_list(
    "shear", ["#2a55c8", "#6a86d2", "#9fa3ab", "#b4b4b4", "#d2907c",
              "#e0553c", "#cc2412"])

_PAL = {
    "bg": "#000000",
    "text": "#e8ecf4",
    "dim": "#8892a4",
    "stretch": "#e0483c",
    "squeeze": "#3f6ff0",
    "patch_far": "#f2e3a4",
    "patch_focus": "#7bf0c8",
}


class FocalPatchFigure:
    """One surface, sparse principal-strain bars, two constant-area patches."""

    def __init__(self, figsize=(11.0, 8.4), t_end: float = RETURN_TIME) -> None:
        self.t_end = float(t_end)
        self.s = EmbeddedTidalSurface(
            sim=Spin2WaveSim(n=1200, pulse_width=PULSE_WIDTH),
            n_theta=121, n_phi=181)
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        self.ax3d = self.fig.add_axes((0.0, 0.0, 1.0, 0.92), projection="3d",
                                      facecolor=_PAL["bg"])

        # The camera sweeps through the x–y plane from the source to the
        # antipode, so the patches are put on the meridian it faces.
        self.patch_azimuth = 1.5 * math.pi
        self.focus_distance = math.pi - 0.94 * PULSE_WIDTH
        self.discs = [
            (name, colour,
             MaterialPatch(self.s, centre_distance=d,
                           centre_azimuth=self.patch_azimuth,
                           radius=PATCH_RADIUS, n_rings=8, n_spokes=64))
            for name, d, colour in (
                ("mid-latitude", 1.20, _PAL["patch_far"]),
                ("on the focal ring", self.focus_distance, _PAL["patch_focus"]))
        ]

        # Roughly area-uniform sampling: the azimuthal count follows sin d, so
        # the bars do not pile into a rosette at the source or at the focus.
        self.cross_samples = []
        for d in np.linspace(0.16, math.pi - 0.16, 15):
            n_psi = max(4, int(round(14.0 * math.sin(d))))
            for psi in np.linspace(0.0, 2.0 * math.pi, n_psi, endpoint=False):
                self.cross_samples.append((float(d), float(psi + 0.4 * d)))
        self.scale_shear = self._calibrate_shear()

    def _calibrate_shear(self, samples: int = 160) -> float:
        """The run's own peak ``|λ|``, so bar lengths mean one thing throughout."""
        self.s.reset()
        peak = 0.0
        for i in range(samples):
            self.s.advance_to((i + 1) * self.t_end / samples)
            peak = max(peak, float(np.max(np.abs(self.s.profiles()["shear"]))))
        self.s.reset()
        return peak

    # ── 1. the continuous surface ───────────────────────────────────────────
    def _surface(self) -> None:
        X, Y, Z = self.s.mesh()
        shear = self.s.shear_on_mesh()
        peak = float(np.max(np.abs(self.s.profiles()["shear"]))) or 1.0
        c = np.tanh(2.2 * shear / peak)
        self.ax3d.plot_surface(X, Y, Z, facecolors=_SHEAR_CMAP(0.5 * (1.0 + c)),
                               rstride=1, cstride=1, antialiased=True,
                               shade=True, alpha=1.0, linewidth=0)

    # ── 2. sparse principal-strain vectors ──────────────────────────────────
    def _crosses(self, view: np.ndarray) -> None:
        """Two bars per sample, along the eigenvectors of ``h_ab``, length ∝ |λ|.

        They are always equal in length to each other: trace-free in two
        dimensions means ``λ± = ±|h|``, so an asymmetric cross would be a bug.
        """
        run_peak = self.scale_shear or 1.0
        for d0, psi in self.cross_samples:
            lam = abs(float(self.s.principal_axes(np.array([d0]))["lambda_plus"][0]))
            if lam < 0.03 * run_peak:
                continue
            centre = self.s.positions(np.array([[d0]]), np.array([[psi]]))[0, 0]
            if float(centre @ view) < 0.10 * float(np.linalg.norm(centre)):
                continue                                       # far side
            # Near either pole the ring is shorter than a full-length bar would
            # be, so a tangential bar would overshoot its own circle; taper it.
            taper = min(1.0, math.sin(d0) / 0.42)
            c = self.s.principal_cross(d0, psi, size=0.175 * taper,
                                       reference=run_peak, lift=0.022)
            a = float(np.clip(lam / run_peak, 0.0, 1.0))
            for key, colour in (("stretch", _PAL["stretch"]),
                                ("squeeze", _PAL["squeeze"])):
                bar = c[key]
                self.ax3d.plot(bar[:, 0], bar[:, 1], bar[:, 2], color=colour,
                               alpha=float(np.clip(0.55 + 0.45 * a, 0.0, 1.0)),
                               linewidth=1.5 + 2.8 * a,
                               solid_capstyle="round", zorder=50)

    # ── 3. the advected constant-area patches ───────────────────────────────
    def _patches(self, view: np.ndarray) -> list:
        rows = []
        for name, colour, disc in self.discs:
            lift = 1.0 + 0.020
            facing = float(disc.points().mean(axis=0) @ view)
            near = facing > 0.0
            poly = Poly3DCollection(disc.triangles() * lift, facecolors=colour,
                                    edgecolors="none",
                                    alpha=0.62 if near else 0.06,
                                    zorder=60)
            self.ax3d.add_collection3d(poly)
            b = disc.boundary() * (lift + 0.004)
            self.ax3d.plot(b[:, 0], b[:, 1], b[:, 2], color=colour,
                           linewidth=2.4 if near else 0.7,
                           alpha=1.0 if near else 0.22, zorder=61)
            if near:                       # the patch's own long axis
                sh = disc.shape()
                c0, u = sh["centroid"] * lift, sh["long_axis"]
                seg = np.vstack([c0 - 0.11 * u, c0 + 0.11 * u])
                self.ax3d.plot(seg[:, 0], seg[:, 1], seg[:, 2], color="#101010",
                               linewidth=1.4, alpha=0.75, zorder=62)
            sh = disc.shape()
            rows.append((name, colour, float(sh["aspect_ratio"]),
                         disc.area() / disc.area(gain=0.0) - 1.0))
        return rows

    # ── one frame ───────────────────────────────────────────────────────────
    def draw(self, t: float, azim: Optional[float] = None) -> None:
        if t < self.s.t:
            self.s.reset()
        self.s.advance_to(t)
        ax = self.ax3d
        ax.clear()
        ax.set_facecolor(_PAL["bg"])

        # The camera follows the wave: azimuth 0 looks at the source, 180 at
        # the antipode, so the refocus happens facing us rather than behind.
        az = (180.0 * t / math.pi) if azim is None else azim
        e, a = math.radians(15.0), math.radians(az)
        view = np.array([math.cos(e) * math.cos(a),
                         math.cos(e) * math.sin(a), math.sin(e)])

        self._surface()
        self._crosses(view)
        rows = self._patches(view)

        lim = 1.16
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        ax.set_box_aspect((1, 1, 1), zoom=1.62)
        ax.set_axis_off()
        ax.view_init(elev=15.0, azim=az)
        self._caption(t, rows)

    def _caption(self, t: float, rows: list) -> None:
        for artist in list(self.fig.texts):
            artist.remove()
        f = self.fig
        f.suptitle("Geometrodynamic QED  v45  —  a trace-free deformation "
                   "through its own antipodal refocus",
                   color=_PAL["text"], fontsize=12.4, fontweight="bold", y=0.977)
        f.text(0.5, 0.944,
               "the strains refocus on a ring, never on the point — and the "
               "area law fails there first",
               color=_PAL["dim"], fontsize=9.2, family="monospace", ha="center")

        phase = ("approaching the focus" if t < ANTIPODAL_TIME - 0.25
                 else "at the antipodal focus" if t < ANTIPODAL_TIME + 0.25
                 else "past the focus, returning")
        f.text(0.018, 0.908, f"t = {t:5.3f}   {phase}", color=_PAL["text"],
               family="monospace", fontsize=9.6, va="top")
        f.text(0.982, 0.908,
               "bars ∝ |λ| along the eigenvectors of h_ab\n"
               "red λ₊ (stretch)    blue λ₋ (squeeze)",
               color=_PAL["dim"], family="monospace", fontsize=8.4,
               va="top", ha="right")
        # ...and the patch readout, low enough to stay clear of the sphere
        y = 0.088
        for name, colour, ratio, darea in rows:
            f.text(0.018, y,
                   f"patch · {name:<18s} aspect {ratio:5.3f}   "
                   f"δarea {100 * darea:+6.2f}%",
                   color=colour, family="monospace", fontsize=9.0, va="bottom")
            y -= 0.030
        f.text(0.018, y - 0.004,
               "δarea is second order in the display gain ε, and vanishes with it",
               color=_PAL["dim"], family="monospace", fontsize=7.8, va="bottom")


# ════════════════════════════════════════════════════════════════════════════
def still(path: str) -> None:
    holder = FocalPatchFigure()
    # in flight; the refocus; and the mid-latitude patch's own worst moment,
    # so the two patches can be compared at the peak of each
    times = [1.35, 3.079, 4.964]
    # The movie's camera tracks the wave; the sheet instead faces each frame's
    # subject, so the last panel looks at the mid-latitude patch rather than
    # at wherever the front happens to be.
    azims = [None, None, 69.0]
    frames = []
    for t, az in zip(times, azims):
        holder.draw(t, azim=az)
        holder.fig.canvas.draw()
        frames.append(np.asarray(holder.fig.canvas.buffer_rgba()).copy())
    plt.close(holder.fig)
    fig, axes = plt.subplots(len(frames), 1, figsize=(10.0, 7.7 * len(frames)),
                             facecolor=_PAL["bg"])
    labels = ["in flight", "the refocus", "the mid-latitude patch's own peak"]
    for ax, img, t, lab in zip(np.atleast_1d(axes), frames, times, labels):
        ax.imshow(img)
        ax.set_axis_off()
        ax.set_title(f"t = {t:.3f}   {lab}", color=_PAL["text"], fontsize=11)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.01, hspace=0.05)
    fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def animate(save: Optional[str] = None, frames: int = 120):
    from matplotlib import animation

    holder = FocalPatchFigure()

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
