#!/usr/bin/env python3
"""
Geometrodynamic QED — v40: the ring defect
==========================================

Successor to ``archive/geometrodynamics_v39.py``.  v39 rendered the vacuole,
the Hopf structure and an exchange cycle, and drew the wormhole as a *given*
object.  This one renders the step before that: **where a throat can come
from**, and why a ring is needed to make one.

It is driven entirely by :mod:`geometrodynamics.viz.radial_caustic`, whose
claims are closed-form and covered by
``experiments/closure_ledger/ring_caustic_defect_probe.py`` (8/8) — so every
line drawn here is a computed quantity, not an artist's impression.

What you are looking at
───────────────────────
Two fronts are launched into the same bulk at the same speed:

* a **pulse** from a point on the outer sphere.  Its front is the metric
  sphere ``|x − P| = t``.  The focal set of a point is empty, so the front is
  embedded at every ``t``: it never touches itself, and behind it is just the
  filled ball.  It crosses the bulk at ``t = ΔR`` and does nothing else.  It
  fills its own void.

* a **ring** — the circle of latitude ``cos θ₀ = R_in/R_out``.  Its front is
  the offset tube of a circle, and a circle *does* have a focal set: every one
  of its points shares the same centre of curvature.  So the focal set is a
  single point, the entire ring arrives there at once at ``t = ρ``, and the
  front stops being embedded.  That is a codimension-2 defect of the
  wavefront, made by geometry alone.

The ring is chosen so its defect lands on the inner sphere — and that same
ring turns out to launch at exactly the grazing angle ``sin α = R_in/R_out``,
tangent to the inner sphere.  The ring that focuses on the throat and the ray
that grazes it are the same ring.

The panel also shows the asymmetry the two surfaces have across the same
bulk: outward from the inner sphere every ray escapes, inward from the outer
sphere only ``1 − √(1 − (R_in/R_out)²)`` of the hemisphere arrives.

Honest labelling
────────────────
The **wavefront defect is computed**.  The throat that would nucleate there is
**not**: this programme has no backreaction, so the geometry cannot respond to
the focus.  Where the renderer marks a would-be throat it says ``schematic``,
and the shell does not actually deform in any solved sense — the deformation
drawn is a cue keyed to the computed defect time, nothing more.

Usage
─────
    python scripts/geometrodynamics_v40_ring_defect.py            # animate
    python scripts/geometrodynamics_v40_ring_defect.py --still out.png
    python scripts/geometrodynamics_v40_ring_defect.py --save out.mp4
"""

from __future__ import annotations

import argparse
import math
import sys
from typing import Optional, Tuple

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)

from geometrodynamics.viz.radial_caustic import (
    RingSource,
    ShellGeometry,
    _C,
    draw_front,
    draw_shell,
    measure_acceptance_asymmetry,
    measure_critical_ring,
)
from geometrodynamics.viz.throat_wavefront import _PAL

# ── the programme's vacuole (v39 constants: R_MID = 1, ΔR/2 = 0.26) ─────────
R_MID, DELTA = 1.00, 0.26
R_INNER, R_OUTER = R_MID - DELTA, R_MID + DELTA

FPS = 30
T_MAX_FACTOR = 1.30          # run to 1.3 ρ so the fold is well past
FLASH_WIDTH = 0.05           # how long the defect marker flares, in t


# ════════════════════════════════════════════════════════════════════════════
# geometry helpers
# ════════════════════════════════════════════════════════════════════════════
def sphere_wire(r: float, n_u: int = 24, n_v: int = 16):
    u = np.linspace(0, 2 * np.pi, n_u)
    v = np.linspace(0, np.pi, n_v)
    x = r * np.outer(np.cos(u), np.sin(v))
    y = r * np.outer(np.sin(u), np.sin(v))
    z = r * np.outer(np.ones_like(u), np.cos(v))
    return x, y, z


def _clip_to_bulk(x, y, z, shell: ShellGeometry):
    """Blank whatever leaves the bulk — the front only exists in the shell."""
    r = np.sqrt(x * x + y * y + z * z)
    out = (r > shell.r_outer) | (r < shell.r_inner - 1e-9)
    x, y, z = x.copy(), y.copy(), z.copy()
    x[out] = np.nan
    y[out] = np.nan
    z[out] = np.nan
    return x, y, z


def tube(ring: RingSource, t: float, shell: ShellGeometry,
         n_u: int = 72, n_v: int = 34):
    """Offset tube of the ring — the ring's wavefront at time ``t``.

    Self-intersects once ``t > ρ``; that is the point of the picture, so the
    surface is drawn as it comes out rather than trimmed, then clipped to the
    bulk because the front does not exist outside it.
    """
    u = np.linspace(0, 2 * np.pi, n_u)
    v = np.linspace(0, 2 * np.pi, n_v)
    U, V = np.meshgrid(u, v)
    rad = ring.radius + t * np.cos(V)
    return _clip_to_bulk(rad * np.cos(U), rad * np.sin(U),
                         ring.centre[2] + t * np.sin(V), shell)


def pulse_sphere(centre: np.ndarray, t: float, shell: ShellGeometry,
                 n_u: int = 48, n_v: int = 26):
    x, y, z = sphere_wire(t, n_u, n_v)
    return _clip_to_bulk(x + centre[0], y + centre[1], z + centre[2], shell)


def shell_deformation(ring: RingSource, t: float, amp: float = 0.055) -> float:
    """A *cue*, not a solution: how far to dimple the inner sphere.

    Zero before the computed defect time, then saturating.  There is no
    backreaction in this programme, so this is presentation only and is
    labelled as such wherever it is drawn.
    """
    tf = ring.self_intersection_time
    if t <= tf:
        return 0.0
    return float(amp * (1.0 - math.exp(-(t - tf) / (0.25 * tf))))


# ════════════════════════════════════════════════════════════════════════════
# figure
# ════════════════════════════════════════════════════════════════════════════
class RingDefectFigure:
    """Three panels: the bulk in 3-D, the meridional section, the numbers."""

    def __init__(self, shell: ShellGeometry, figsize=(15.0, 7.4)) -> None:
        self.shell = shell
        self.ring = shell.critical_ring()
        self.pulse = shell.point_source()
        self.crit = measure_critical_ring(shell)
        self.acc = measure_acceptance_asymmetry(shell, n_samples=60000, seed=7)
        self.t_max = T_MAX_FACTOR * self.ring.self_intersection_time

        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 0.95],
                                   wspace=0.06)
        self.ax3d = self.fig.add_subplot(gs[0, 0], projection="3d")
        self.ax2d = self.fig.add_subplot(gs[0, 1])
        self.axtx = self.fig.add_subplot(gs[0, 2])
        self.fig.suptitle(
            "Geometrodynamic QED  v40  —  the ring defect",
            color=_PAL["text"], fontsize=13, y=0.97)

    # ── static furniture ────────────────────────────────────────────────────
    def _style3d(self) -> None:
        ax = self.ax3d
        ax.set_facecolor(_PAL["bg"])
        lim = 1.15 * self.shell.r_outer
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        ax.set_box_aspect((1, 1, 1))
        ax.set_axis_off()

    def _style2d(self, t: float) -> None:
        ax = self.ax2d
        ax.set_facecolor(_PAL["panel"])
        for sp in ax.spines.values():
            sp.set_color(_PAL["border"])
        lim = 1.15 * self.shell.r_outer
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("meridional section  ·  both fronts, one clock",
                     color=_PAL["text"], fontsize=9, pad=6)

    # ── per-frame ───────────────────────────────────────────────────────────
    def draw(self, t: float) -> None:
        sh, ring, pulse = self.shell, self.ring, self.pulse
        folded = t >= ring.self_intersection_time
        dent = shell_deformation(ring, t)

        # ---------------- 3-D ----------------
        ax = self.ax3d
        ax.clear()
        self._style3d()
        x, y, z = sphere_wire(sh.r_outer)
        ax.plot_wireframe(x, y, z, color=_PAL["sphere"], lw=0.35, alpha=0.30)
        x, y, z = sphere_wire(sh.r_inner - dent, 30, 20)
        ax.plot_surface(x, y, z, color=_PAL["plug"], alpha=0.16,
                        linewidth=0, shade=False)

        rp = ring.points(140)
        ax.plot(rp[:, 0], rp[:, 1], rp[:, 2], color=_C["ring"], lw=2.0)

        tx, ty, tz = tube(ring, t, sh)
        ax.plot_wireframe(tx, ty, tz, color=_C["ring"], lw=0.4,
                          alpha=0.55 if not folded else 0.75)

        px, py, pz = pulse_sphere(pulse.position, t, sh)
        ax.plot_wireframe(px, py, pz, color=_C["pulse"], lw=0.4, alpha=0.45)

        c = ring.centre
        if folded:
            flare = max(0.0, 1.0 - (t - ring.self_intersection_time) / FLASH_WIDTH)
            ax.scatter([c[0]], [c[1]], [c[2]], s=90 + 320 * flare,
                       color=_C["defect"], depthshade=False, zorder=10)
        ax.view_init(elev=18, azim=-58)

        # ---------------- meridional ----------------
        ax = self.ax2d
        ax.clear()
        self._style2d(t)
        a = np.linspace(0, 2 * np.pi, 400)
        ax.plot(sh.r_outer * np.cos(a), sh.r_outer * np.sin(a),
                color=_PAL["sphere"], lw=1.5)
        ax.plot((sh.r_inner - dent) * np.cos(a), (sh.r_inner - dent) * np.sin(a),
                color=_PAL["plug"], lw=1.5)

        # the grazing rays that define the critical ring
        for s in (+1, -1):
            p = np.array([s * ring.radius, ring.centre[2]])
            ax.plot([p[0], 0.0], [p[1], ring.centre[2]],
                    color=_PAL["dim"], lw=0.8, ls=":")
        draw_front(ax, pulse, t, lw=1.8)
        draw_front(ax, ring, t, lw=1.8)
        ax.plot([0.0], [sh.r_outer], "o", color=_C["pulse"], ms=5)
        ax.plot([ring.radius, -ring.radius], [ring.centre[2]] * 2, "o",
                color=_C["ring"], ms=5)
        if folded:
            ax.plot([0.0], [ring.centre[2]], "*", color=_C["defect"], ms=17,
                    zorder=6)

        # ---------------- numbers ----------------
        ax = self.axtx
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        for sp in ax.spines.values():
            sp.set_color(_PAL["border"])
        ax.set_xticks([])
        ax.set_yticks([])
        mono = dict(family="monospace", fontsize=8.3,
                    transform=ax.transAxes, va="top")
        y = 0.955
        ax.text(0.5, y, "THE RING DEFECT", color=_PAL["text"], ha="center",
                family="monospace", fontsize=10, transform=ax.transAxes,
                va="top")
        y -= 0.050
        lines = [
            ("", ""),
            (f"t = {t:6.3f}", _PAL["text"]),
            ("", ""),
            ("PULSE  point source, outer sphere", _C["pulse"]),
            (f"  front  = sphere |x-P| = t", _PAL["dim"]),
            (f"  focal set  = empty", _PAL["dim"]),
            (f"  self-intersects  = never", _C["pulse"]),
            (f"  crosses bulk at  dR = {sh.gap:.4f}", _PAL["dim"]),
            ("", ""),
            ("RING   circle of latitude", _C["ring"]),
            (f"  theta0 = {self.crit['polar_angle_deg']:.3f} deg", _PAL["dim"]),
            (f"  rho    = {ring.radius:.4f}", _PAL["dim"]),
            (f"  focal set  = 1 point (its centre)", _PAL["dim"]),
            (f"  folds at   t = rho = {ring.radius:.4f}", _C["ring"]),
            ("", ""),
            ("THE COINCIDENCE", _PAL["text"]),
            (f"  defect radius  = {self.crit['defect_radius']:.6f}", _PAL["dim"]),
            (f"  R_inner        = {self.crit['r_inner']:.6f}", _PAL["dim"]),
            (f"  error          = {self.crit['defect_on_inner_error']:.1e}",
             _PAL["dim"]),
            (f"  launch sin a   = {self.crit['launch_sin']:.6f}", _PAL["dim"]),
            (f"  grazing sin a  = {self.crit['critical_sin']:.6f}", _PAL["dim"]),
            (f"  error          = {self.crit['grazing_error']:.1e}", _PAL["dim"]),
            ("  -> the focusing ring IS the grazing ring", _PAL["text"]),
            ("", ""),
            ("ASYMMETRY across the same bulk", _PAL["text"]),
            (f"  outer->inner  {100*self.acc['inward_closed_form']:5.1f}% "
             f"(MC {100*self.acc['inward_monte_carlo']:.1f}%)", _PAL["dim"]),
            (f"  inner->outer  {100*self.acc['outward_closed_form']:5.1f}%",
             _PAL["dim"]),
            (f"  ratio         {self.acc['asymmetry_ratio']:.2f}x", _PAL["text"]),
        ]
        for txt, col in lines:
            if txt:
                ax.text(0.045, y, txt, color=col or _PAL["dim"], **mono)
            y -= 0.0258

        state = ("FRONT EMBEDDED  —  no defect" if not folded
                 else "FRONT FOLDED  —  defect on the inner sphere")
        ax.text(0.045, y - 0.022, state,
                color=_C["defect"] if folded else _PAL["dim"],
                family="monospace", fontsize=9, transform=ax.transAxes,
                va="top")
        ax.text(0.045, 0.030,
                "the wavefront defect is computed; the throat\n"
                "that would form there is not — no backreaction\n"
                "in this model.  the shell dimple is a schematic\n"
                "cue keyed to the computed defect time.",
                color=_PAL["dim"], family="monospace", fontsize=7.0,
                transform=ax.transAxes, va="bottom")


# ════════════════════════════════════════════════════════════════════════════
def build(shell: Optional[ShellGeometry] = None) -> RingDefectFigure:
    return RingDefectFigure(shell or ShellGeometry(R_INNER, R_OUTER))


def still(path: str, n: int = 4) -> None:
    """A contact sheet: the meridional panel at ``n`` times through the fold."""
    fig_holder = build()
    sh, ring = fig_holder.shell, fig_holder.ring
    plt.close(fig_holder.fig)

    times = np.linspace(0.35 * ring.radius, 1.25 * ring.radius, n)
    fig, axes = plt.subplots(1, n, figsize=(3.5 * n, 3.9), facecolor=_PAL["bg"])
    for ax, t in zip(np.atleast_1d(axes), times):
        folded = t >= ring.self_intersection_time
        ax.set_facecolor(_PAL["panel"])
        draw_shell(ax, sh)
        draw_front(ax, sh.point_source(), t, lw=1.6)
        draw_front(ax, ring, t, lw=1.6)
        if folded:
            ax.plot([0.0], [ring.centre[2]], "*", color=_C["defect"], ms=15)
        lim = 1.15 * sh.r_outer
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color(_PAL["border"])
        ax.set_title(f"t = {t:.3f}" + ("   defect" if folded else ""),
                     color=_C["defect"] if folded else _PAL["text"], fontsize=9)
    fig.suptitle(
        f"pulse never folds  ·  ring folds onto its centre at t = ρ = "
        f"{ring.radius:.4f}, on r = {ring.centre_radius:.4f} = R_inner",
        color=_PAL["text"], fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=120, facecolor=fig.get_facecolor())
    print(f"wrote {path}")


def animate(save: Optional[str] = None, frames: int = 150):
    from matplotlib import animation

    holder = build()

    def update(i: int):
        holder.draw((i + 1) * holder.t_max / frames)
        return []

    anim = animation.FuncAnimation(holder.fig, update, frames=frames,
                                   interval=1000 // FPS, blit=False)
    if save:
        anim.save(save, fps=FPS, dpi=110,
                  savefig_kwargs={"facecolor": _PAL["bg"]})
        print(f"wrote {save}")
    else:
        plt.show()
    return anim


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG",
                    help="write a contact sheet through the fold and exit")
    ap.add_argument("--save", metavar="FILE", help="render the animation to a file")
    ap.add_argument("--frames", type=int, default=150)
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
    sys.exit(main())
