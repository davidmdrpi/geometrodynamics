#!/usr/bin/env python3
"""
Geometrodynamic QED — v49: draw the vectors, not their tips, and they intersect
===============================================================================

Four rounds of this series established that a height field cannot wind and
cannot self-intersect at any gap or gain. All four were about the same object:
**the graph of the displacement's tips**, `r = f(σ)`, which is embedded by
construction.

Draw the **vectors themselves** and the obstruction is gone, for a reason that is
entirely classical: *neighbouring normals to a curve meet at its centre of
curvature*. The normal family has an envelope — the evolute — and a normal of
length `L` crosses its neighbours as soon as `L` exceeds the local radius of
curvature `ρ = 1/κ`. Nothing is added by hand.

The same wave, the same instant, two objects: the graph gives `0`
self-intersections, the normal field gives `520`.

What is on screen
─────────────────
* **the graph of the tips**, with its self-intersection count sitting at zero
  however hard the wave is driven;
* **the normal field**, the same wave drawn as vectors, with every crossing
  marked — those dots are the caustic;
* **the reset**, where a vector leaving through `R_outer` comes back at
  `R_inner` at the angle it left and crosses things it could never have reached;
* **the curvature history**, `ρ_min` falling as the ring converges, against the
  length at which the drawn bundle first crosses.

The threshold is the ring concentration, at last visible
────────────────────────────────────────────────────────
`ρ_min` falls `0.1408 → 0.0540` between mid-flight and the focus — the
converging ring **sharpens its own surface** by `2.61×` — and the first drawn
crossing falls with it, `0.492 → 0.189`. This is where the ring's concentration
finally shows up. Not as height, which barely moves and never beats the launch,
but as *curvature*, which is what decides whether the normals meet.

And the gap matters again
─────────────────────────
The vector length is what spans the gap, so `L` and `δ` are one knob. At `L = δ`:
`δ = 0.40` gives `382` crossings, `δ = 0.09` gives `0` from the normals alone and
`472` once they re-enter. Reducing the shell separation now *produces*
intersections rather than being unable to.

Honest scope
────────────
The vector length `L` is a display choice, like every gain in this series. The
directions and the curvature are the surface's own, and the crossings are
counted with the same segment predicate used for the earlier negative results.

The drawn bundle is deliberately **sparse** (`stride = 26`) so the vectors read
as vectors, so its crossing counts are far below the dense counts the module
measures — the same geometry, sampled for legibility. The dense bundle is what
`normal_field_probe` reports.

Usage
─────
    python scripts/geometrodynamics_v49_normal_field.py             # animate
    python scripts/geometrodynamics_v49_normal_field.py --still out.png
    python scripts/geometrodynamics_v49_normal_field.py --gaps out.png
    python scripts/geometrodynamics_v49_normal_field.py --save out.gif
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.viz.circle_slice import ANTIPODAL_TIME, TWO_PI
from geometrodynamics.viz.normal_field import (
    NormalField,
    segment_crossing_points,
)
from geometrodynamics.viz.slice_folding import self_intersections

FPS = 14
PULSE_WIDTH = 0.18
DELTA = 0.26
GAIN = 0.30
LENGTH = 0.35
STRIDE = 26                    # sparse enough to read as vectors

_PAL = {
    "bg": "#05070c",
    "panel": "#080b13",
    "text": "#e8ecf4",
    "dim": "#7d8798",
    "faint": "#2a3244",
    "shell": "#39445c",
    "mid": "#4a5675",
    "curve": "#ffd166",
    "vec": "#7cc4ff",
    "stub": "#ff9f43",
    "cross": "#ff3b6b",
    "rho": "#c08cff",
    "source": "#8ef0c0",
    "antipode": "#ff8ad0",
}


class NormalFieldFigure:
    """The graph that cannot cross, the vectors that do, and why."""

    def __init__(self, figsize=(13.4, 8.6),
                 t_end: float = ANTIPODAL_TIME * 1.10) -> None:
        self.t_end = float(t_end)
        self.nf = NormalField(delta=DELTA, pulse_width=PULSE_WIDTH,
                              n_sigma=4001, gain=GAIN, stride=STRIDE)
        self.dense = NormalField(delta=DELTA, pulse_width=PULSE_WIDTH,
                                 n_sigma=4001, gain=GAIN, stride=8)
        self.shells = self.nf.shells
        self.history = []

        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 3, height_ratios=[1.0, 0.72],
            left=0.04, right=0.975, top=0.858, bottom=0.115,
            wspace=0.16, hspace=0.30)
        self.ax_graph = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["bg"])
        self.ax_norm = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["bg"])
        self.ax_reset = self.fig.add_subplot(gs[0, 2], facecolor=_PAL["bg"])
        self.ax_rho = self.fig.add_subplot(gs[1, :2], facecolor=_PAL["panel"])
        self.ax_count = self.fig.add_subplot(gs[1, 2], facecolor=_PAL["panel"])

    def advance_to(self, t: float) -> None:
        for f in (self.nf, self.dense):
            if t < f.t:
                f.reset()
            f.advance_to(t)

    # ── furniture ───────────────────────────────────────────────────────────
    def _annulus(self, ax) -> None:
        a = np.linspace(0.0, TWO_PI, 361)
        for r, colour, ls in ((self.shells.r_outer, _PAL["shell"], "-"),
                              (self.shells.r_inner, _PAL["shell"], "-"),
                              (self.shells.r_mid, _PAL["mid"], (0, (4, 4)))):
            ax.plot(r * np.cos(a), r * np.sin(a), color=colour, lw=1.1, ls=ls,
                    alpha=0.85, zorder=4)

    def _slice_curve(self, ax, lw: float = 1.6, alpha: float = 1.0) -> None:
        X, _, _ = self.nf.frame()
        ax.plot(X[:, 0], X[:, 1], color=_PAL["curve"], lw=lw, alpha=alpha,
                zorder=10)

    def _square(self, ax, lim: float) -> None:
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_axis_off()

    # ── panels ──────────────────────────────────────────────────────────────
    def _graph(self) -> None:
        """The object four rounds of negative results were about."""
        ax = self.ax_graph
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        self._annulus(ax)
        self._slice_curve(ax, lw=2.0)
        pts = self.nf.slice.points(gain=GAIN)
        hits = self_intersections(pts)
        self._square(ax, self.shells.r_outer * 1.22)
        ax.set_title("the graph of the tips", color=_PAL["text"],
                     fontsize=10.2, pad=6)
        ax.text(0.0, -self.shells.r_outer * 1.12,
                f"self-intersections {hits}", color=_PAL["curve"],
                family="monospace", fontsize=9.4, ha="center",
                fontweight="bold")

    def _normals(self) -> None:
        """The same wave, drawn as vectors.  The dots are the caustic."""
        ax = self.ax_norm
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        self._annulus(ax)
        self._slice_curve(ax, lw=1.2, alpha=0.55)
        a, b = self.nf.vectors(LENGTH)
        segs = np.stack([a, b], axis=1)
        for s in segs:
            ax.plot(s[:, 0], s[:, 1], color=_PAL["vec"], lw=1.0, alpha=0.85,
                    zorder=12)
        pts = segment_crossing_points(a, b)
        if len(pts):
            ax.plot(pts[:, 0], pts[:, 1], ".", ms=4.2, color=_PAL["cross"],
                    zorder=25)
        self._square(ax, self.shells.r_outer * 1.22)
        ax.set_title(f"the same wave as vectors   L = {LENGTH:.2f}",
                     color=_PAL["text"], fontsize=10.2, pad=6)
        ax.text(0.0, -self.shells.r_outer * 1.12,
                f"crossings {len(pts)}", color=_PAL["cross"],
                family="monospace", fontsize=9.4, ha="center",
                fontweight="bold")

    def _reset(self) -> None:
        """Vectors clipped at R_outer, and the stubs that come back inside."""
        ax = self.ax_reset
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        self._annulus(ax)
        self._slice_curve(ax, lw=1.2, alpha=0.55)
        base, clipped, stub_base, stub_tip = self.nf.vectors_with_reset(LENGTH)
        for p, q in zip(base, clipped):
            ax.plot([p[0], q[0]], [p[1], q[1]], color=_PAL["vec"], lw=1.0,
                    alpha=0.75, zorder=12)
        for p, q in zip(stub_base, stub_tip):
            ax.plot([p[0], q[0]], [p[1], q[1]], color=_PAL["stub"], lw=1.2,
                    alpha=0.95, zorder=14)
        a = np.vstack([base, stub_base])
        b = np.vstack([clipped, stub_tip])
        pts = segment_crossing_points(a, b)
        if len(pts):
            ax.plot(pts[:, 0], pts[:, 1], ".", ms=4.2, color=_PAL["cross"],
                    zorder=25)
        self._square(ax, self.shells.r_outer * 1.22)
        ax.set_title("...with the reset at R_inner", color=_PAL["text"],
                     fontsize=10.2, pad=6)
        ax.text(0.0, -self.shells.r_outer * 1.12,
                f"crossings {len(pts)}    wrapped {len(stub_base)}",
                color=_PAL["stub"], family="monospace", fontsize=9.4,
                ha="center", fontweight="bold")

    def _rho_history(self, t: float) -> None:
        """The converging ring sharpening its own surface."""
        ax = self.ax_rho
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        if not self.history:
            probe = NormalField(delta=DELTA, pulse_width=PULSE_WIDTH,
                                n_sigma=2001, gain=GAIN, stride=8)
            probe.reset()
            for k in range(140):
                probe.advance_to((k + 1) * self.t_end / 140)
                self.history.append((probe.t, probe.envelope_distance()))
        shown = [h for h in self.history if h[0] <= t + 1e-9]
        if shown:
            ax.plot([h[0] for h in shown], [h[1] for h in shown],
                    color=_PAL["rho"], lw=1.7)
            ax.plot([shown[-1][0]], [shown[-1][1]], "o", ms=5.5,
                    color=_PAL["cross"])
        ax.axhline(LENGTH, color=_PAL["vec"], lw=1.1, ls=(0, (4, 3)))
        ax.text(0.06, LENGTH * 1.06, f"drawn vector length L = {LENGTH:.2f}",
                color=_PAL["vec"], family="monospace", fontsize=7.8)
        ax.axvline(ANTIPODAL_TIME, color=_PAL["faint"], lw=1.0)
        ax.text(ANTIPODAL_TIME - 0.06, 0.62, "focus", color=_PAL["dim"],
                family="monospace", fontsize=7.8, ha="right")
        ax.set_xlim(0.0, self.t_end)
        ax.set_ylim(0.0, 0.70)
        ax.set_xlabel("t", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["faint"])
        ax.set_title("ρ_min — the ring sharpens its own surface, and the "
                     "normals cross below the dashed line",
                     color=_PAL["text"], fontsize=9.4, pad=5)

    def _count(self) -> None:
        """Crossings against vector length, with and without the reset."""
        ax = self.ax_count
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        Ls = np.linspace(0.04, 0.60, 15)
        plain = [self.nf.crossings(float(L)) for L in Ls]
        both = [self.nf.crossings(float(L), with_reset=True) for L in Ls]
        ax.plot(Ls, plain, color=_PAL["vec"], lw=1.6, label="normals")
        ax.plot(Ls, both, color=_PAL["stub"], lw=1.6, label="with reset")
        rho = self.nf.envelope_distance()
        ax.axvline(rho, color=_PAL["rho"], lw=1.1, ls=(0, (3, 3)))
        ax.text(rho + 0.012, max(max(both), 1) * 0.90, "ρ_min",
                color=_PAL["rho"], family="monospace", fontsize=7.8)
        ax.axvline(LENGTH, color=_PAL["faint"], lw=1.0)
        ax.set_xlim(0.0, 0.60)
        ax.set_ylim(0.0, max(max(both), 1) * 1.15)
        ax.set_xlabel("vector length L", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["faint"])
        leg = ax.legend(loc="upper left", fontsize=7.4, frameon=False)
        for txt in leg.get_texts():
            txt.set_color(_PAL["dim"])
        ax.set_title("crossings against length", color=_PAL["text"],
                     fontsize=9.4, pad=5)

    # ── one frame ───────────────────────────────────────────────────────────
    def draw(self, t: float) -> None:
        self.advance_to(t)
        self._graph()
        self._normals()
        self._reset()
        self._rho_history(t)
        self._count()
        self._caption(t)

    def _caption(self, t: float) -> None:
        for artist in list(self.fig.texts):
            artist.remove()
        f = self.fig
        f.suptitle("Geometrodynamic QED  v49  —  draw the vectors, not their "
                   "tips, and they intersect",
                   color=_PAL["text"], fontsize=12.4, fontweight="bold", y=0.972)
        f.text(0.5, 0.938,
               "neighbouring normals meet at the centre of curvature — so the "
               "threshold is ρ, and the focusing drives it down",
               color=_PAL["dim"], fontsize=9.2, family="monospace", ha="center")
        phase = ("outbound" if t < 0.5 * math.pi
                 else "converging" if t < ANTIPODAL_TIME - 0.12
                 else "at the focus" if t < ANTIPODAL_TIME + 0.12
                 else "past it")
        f.text(0.04, 0.908, f"t = {t:5.3f}   {phase}", color=_PAL["text"],
               family="monospace", fontsize=9.4, va="top")
        f.text(0.975, 0.908,
               f"ρ_min {self.nf.envelope_distance():.4f}    "
               f"δ = {DELTA:.2f}    gain {GAIN:.2f}",
               color=_PAL["dim"], family="monospace", fontsize=8.8,
               va="top", ha="right")
        f.text(0.04, 0.028,
               "the graph and the vectors are the same field at the same "
               "instant — one is embedded by construction, the other has an "
               f"envelope · the bundle is sampled every {STRIDE} points so the "
               "vectors read as vectors, so these counts sit below the dense ones",
               color=_PAL["dim"], family="monospace", fontsize=8.2, va="bottom")


# ════════════════════════════════════════════════════════════════════════════
def still(path: str) -> None:
    holder = NormalFieldFigure()
    times = [1.60, 2.70, 3.06]
    labels = ["outbound", "converging", "at the focus"]
    frames = []
    for t in times:
        holder.draw(t)
        holder.fig.canvas.draw()
        frames.append(np.asarray(holder.fig.canvas.buffer_rgba()).copy())
    plt.close(holder.fig)
    fig, axes = plt.subplots(len(frames), 1, figsize=(11.8, 7.6 * len(frames)),
                             facecolor=_PAL["bg"])
    for ax, img, t, lab in zip(np.atleast_1d(axes), frames, times, labels):
        ax.imshow(img)
        ax.set_axis_off()
        ax.set_title(f"t = {t:.3f}   {lab}", color=_PAL["text"], fontsize=11)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.01, hspace=0.05)
    fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def gaps(path: str, deltas=(0.40, 0.26, 0.16, 0.09), t: float = 3.0) -> None:
    """The gap sweep: at ``L = δ`` the reset takes over as the gap closes."""
    fig, axes = plt.subplots(1, len(deltas), figsize=(4.0 * len(deltas), 5.0),
                             facecolor=_PAL["bg"])
    for ax, delta in zip(np.atleast_1d(axes), deltas):
        nf = NormalField(delta=float(delta), pulse_width=PULSE_WIDTH,
                         n_sigma=4001, gain=GAIN, stride=STRIDE)
        nf.reset()
        nf.advance_to(t)
        ax.set_facecolor(_PAL["bg"])
        a = np.linspace(0.0, TWO_PI, 361)
        for r, colour, ls in ((nf.shells.r_outer, _PAL["shell"], "-"),
                              (nf.shells.r_inner, _PAL["shell"], "-"),
                              (nf.shells.r_mid, _PAL["mid"], (0, (4, 4)))):
            ax.plot(r * np.cos(a), r * np.sin(a), color=colour, lw=1.0, ls=ls,
                    alpha=0.85)
        X, _, _ = nf.frame()
        ax.plot(X[:, 0], X[:, 1], color=_PAL["curve"], lw=1.1, alpha=0.55)
        L = float(delta)
        base, clipped, sb, st = nf.vectors_with_reset(L)
        for p, q in zip(base, clipped):
            ax.plot([p[0], q[0]], [p[1], q[1]], color=_PAL["vec"], lw=0.9,
                    alpha=0.75)
        for p, q in zip(sb, st):
            ax.plot([p[0], q[0]], [p[1], q[1]], color=_PAL["stub"], lw=1.1,
                    alpha=0.95)
        pts = segment_crossing_points(np.vstack([base, sb]),
                                      np.vstack([clipped, st]))
        if len(pts):
            ax.plot(pts[:, 0], pts[:, 1], ".", ms=3.4, color=_PAL["cross"])
        lim = nf.shells.r_outer * 1.16
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.set_title(f"δ = {delta:.2f}   L = δ\nnormals alone "
                     f"{nf.crossings(L)}   with reset {len(pts)}",
                     color=_PAL["text"], fontsize=9.4, pad=8)
    fig.suptitle("The gap is a knob on intersections again — as it closes, the "
                 "reset takes over",
                 color=_PAL["text"], fontsize=12.0, fontweight="bold", y=0.985)
    fig.text(0.5, 0.045,
             "blue: the normals, clipped at R_outer · orange: the stub that "
             "re-enters at R_inner · red: where they cross",
             color=_PAL["dim"], family="monospace", fontsize=8.6, ha="center")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.84, bottom=0.10, wspace=0.05)
    fig.savefig(path, dpi=120, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def animate(save: Optional[str] = None, frames: int = 100):
    from matplotlib import animation

    holder = NormalFieldFigure()

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
    ap.add_argument("--gaps", metavar="PNG")
    ap.add_argument("--save", metavar="FILE")
    ap.add_argument("--frames", type=int, default=100)
    a = ap.parse_args(argv)
    if a.still or a.gaps:
        matplotlib.use("Agg")
        if a.still:
            still(a.still)
        if a.gaps:
            gaps(a.gaps)
        return 0
    if a.save:
        matplotlib.use("Agg")
    animate(save=a.save, frames=a.frames)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
