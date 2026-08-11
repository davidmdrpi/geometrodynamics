#!/usr/bin/env python3
"""
Geometrodynamic QED — v48: the ring reaches across, and what it takes to fold
=============================================================================

The scalar wave is not only its focal pulse. A ring leaves the source, thins to
a minimum at the equator, then **grows** as it converges — `0.156 → 0.933`, a
factor of `5.97`, following `1/√(sin d)` to within `3.4%`. All of that happens
*before* the focal pulse. So: can that ring reach across the vacuole, and can it
ever intersect?

Those are two questions, and they have different answers.

Reaching across — yes, easily
─────────────────────────────
The threshold is exactly `δ / max|u|`, so either knob buys it. Shrinking the gap
buys something extra: **lead**, the head start the converging ring gets on the
focal pulse.

    δ = 0.26, ε = 0.40   spans from d = 3.14   lead 0.21
    δ = 0.16, ε = 0.80   spans from d = 2.44   lead 0.81
    δ = 0.09, ε = 0.80   spans from d = 1.83   lead 1.41

At the bottom of that table the ring is spanning the gap from just past the
equator, for the whole converging leg — not an instant at the focus but a
sustained state. The ring gets there long before the pulse.

Intersecting — no, and not for want of trying
─────────────────────────────────────────────
Swept over gap and gain with a real segment-intersection test: **zero**, at up to
`346` seam crossings. A curve drawn as `r = f(σ)` with `f` single-valued puts
exactly one radius at each direction, so it is *embedded* by construction. Two of
its points cannot occupy the same place. This is the winding obstruction seen
from the side — a graph can no more cross itself than it can wind — and no gap
and no energy changes it.

What does fold
──────────────
Let the material move **sideways**. If each point carries a tangential
displacement as well as a radial one,

    σ(σ₀) = σ₀ + λ ε ∂_σ u ,     r(σ₀) = R_mid·exp(ε u) ,

the map `σ₀ ↦ σ` folds where `∂σ/∂σ₀ = 1 + λ ε ∂²_σ u < 0`. Predicted threshold
`λε = 1/max(−∂²_σ u)`, found by bisection to a relative `1.8e-12`. Past it the
curve is multivalued in `σ`, stops being a graph, and self-intersects; below it,
it does not.

**And that threshold does not know the gap exists** — spread `0.0` across
`δ = 0.26, 0.12, 0.05`, while the spanning threshold scales with `δ` directly.
The two knobs are orthogonal: the gap sets when the wave *reaches across*, the
front's curvature sets when it *folds*. Shrinking the vacuole can never buy an
intersection.

What it does scale with is the pulse: `λε ≈ 0.385 w²` across a threefold range of
widths, spread `3.7%`. Narrow fronts fold sooner, because folding is about how
sharply the front turns, not how tall it is.

The ring was the right place to look
────────────────────────────────────
`∂²_σ u` peaks at the steep converging front, so the fold appears there first —
measured at `0.0157` from the antipode, on the ring, at the moment of tightest
focus. The intuition was right about *where*; what was missing is that the
freedom it needs is tangential, and the gap is the wrong knob for buying it.

Honest scope
────────────
The mixing `λ` is a modelling choice, not derived from the scalar equation: it
says how much of the displacement is along the surface rather than across it.
`λ = 0` is exactly the height field of v46. What is derived is the *threshold*
given `λ`, and its independence from the gap.

Usage
─────
    python scripts/geometrodynamics_v48_ring_and_fold.py             # animate
    python scripts/geometrodynamics_v48_ring_and_fold.py --still out.png
    python scripts/geometrodynamics_v48_ring_and_fold.py --save out.gif
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.viz.circle_slice import (
    ANTIPODAL_TIME,
    RETURN_TIME,
    TWO_PI,
    BulkAnnulus,
    CircleSlice,
)
from geometrodynamics.viz.slice_folding import LagrangianSlice, self_intersections
from geometrodynamics.viz.warped_sphere import NestedShells

FPS = 14
PULSE_WIDTH = 0.18
NARROW_DELTA = 0.09           # the reduced gap, where the ring spans early
SPAN_GAIN = 0.80
FOLD_MIXING = 1.0

_PAL = {
    "bg": "#05070c",
    "panel": "#080b13",
    "text": "#e8ecf4",
    "dim": "#7d8798",
    "faint": "#2a3244",
    "shell": "#39445c",
    "mid": "#4a5675",
    "wave": "#ffd166",
    "span": "#ff9f43",
    "hot": "#ff5c5c",
    "cool": "#5cc8ff",
    "fold": "#c08cff",
    "cross": "#ff3b6b",
    "source": "#8ef0c0",
    "antipode": "#ff8ad0",
}


class RingFoldFigure:
    """The ring reaching across, and the fold that finally intersects."""

    def __init__(self, figsize=(13.2, 8.4), t_end: float = ANTIPODAL_TIME * 1.12
                 ) -> None:
        self.t_end = float(t_end)
        self.shells = NestedShells(r_mid=1.0, delta=NARROW_DELTA)

        self.span = CircleSlice(
            bulk=BulkAnnulus(self.shells, mode="conformal"),
            radial_law="multiplicative", pulse_width=PULSE_WIDTH,
            n_sigma=1441, gain=SPAN_GAIN)
        self.graph = CircleSlice(
            bulk=BulkAnnulus(self.shells, mode="conformal"),
            radial_law="multiplicative", pulse_width=PULSE_WIDTH,
            n_sigma=1441, gain=6.0)
        self.lag = LagrangianSlice(mixing=FOLD_MIXING, pulse_width=PULSE_WIDTH,
                                   n_sigma=2001, n_radial=2400,
                                   shells=self.shells)

        self.fold_gain = self._fold_gain()
        self.history = []

        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 3, width_ratios=[1.0, 1.0, 1.12], height_ratios=[1.0, 0.80],
            left=0.045, right=0.975, top=0.862, bottom=0.115,
            wspace=0.22, hspace=0.30)
        self.ax_span = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["bg"])
        self.ax_graph = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["bg"])
        self.ax_fold = self.fig.add_subplot(gs[0, 2], facecolor=_PAL["bg"])
        self.ax_ring = self.fig.add_subplot(gs[1, :2], facecolor=_PAL["panel"])
        self.ax_jac = self.fig.add_subplot(gs[1, 2], facecolor=_PAL["panel"])

    def _fold_gain(self, frames: int = 120) -> float:
        self.lag.reset()
        best = math.inf
        for i in range(frames):
            self.lag.advance_to((i + 1) * ANTIPODAL_TIME / frames)
            best = min(best, self.lag.predicted_fold_product())
        self.lag.reset()
        return 3.0 * best / FOLD_MIXING           # comfortably past threshold

    def advance_to(self, t: float) -> None:
        for s in (self.span, self.graph, self.lag):
            if t < s.t:
                s.reset()
            s.advance_to(t)

    # ── furniture ───────────────────────────────────────────────────────────
    def _annulus(self, ax, lw: float = 1.2) -> None:
        a = np.linspace(0.0, TWO_PI, 361)
        for r, colour, ls in ((self.shells.r_outer, _PAL["shell"], "-"),
                              (self.shells.r_inner, _PAL["shell"], "-"),
                              (self.shells.r_mid, _PAL["mid"], (0, (4, 4)))):
            ax.plot(r * np.cos(a), r * np.sin(a), color=colour, lw=lw, ls=ls,
                    alpha=0.9, zorder=5)

    def _square(self, ax, lim: float) -> None:
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_axis_off()

    def _folded_curve(self, ax, xy: np.ndarray, colour: str) -> None:
        ax.plot(xy[:, 0], xy[:, 1], color=colour, lw=1.9,
                solid_capstyle="round", zorder=15)

    # ── panels ──────────────────────────────────────────────────────────────
    def _spanning(self) -> None:
        """The converging ring, reaching across a reduced gap."""
        ax = self.ax_span
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        self._annulus(ax)
        s = self.span
        r_raw = self.shells.r_mid * np.exp(
            SPAN_GAIN * s.field())                # unwrapped, to judge spanning
        drawn, sheet = s.bulk.wrap(r_raw)
        xy = np.stack([drawn * np.cos(s.sigma), drawn * np.sin(s.sigma)], -1)
        spans = (r_raw >= self.shells.r_outer) | (r_raw <= self.shells.r_inner)
        start = 0
        for c in list(np.nonzero(np.diff(sheet) != 0)[0]) + [len(sheet) - 1]:
            seg = xy[start:c + 1]
            if len(seg) > 1:
                hot = bool(np.any(spans[start:c + 1]))
                ax.plot(seg[:, 0], seg[:, 1],
                        color=_PAL["span"] if hot else _PAL["wave"],
                        lw=2.4 if hot else 1.8, solid_capstyle="round",
                        zorder=15)
            start = c + 1
        ax.plot([self.shells.r_mid], [0.0], "o", ms=4, color=_PAL["source"])
        ax.plot([-self.shells.r_mid], [0.0], "o", ms=4, color=_PAL["antipode"])
        self._square(ax, self.shells.r_outer * 1.20)
        ax.set_title(f"the ring reaching across   δ = {NARROW_DELTA:.2f}",
                     color=_PAL["text"], fontsize=10.0, pad=6)
        if bool(np.any(spans)):
            ax.text(0.0, -self.shells.r_outer * 1.10, "SPANNING",
                    color=_PAL["span"], family="monospace", fontsize=9.0,
                    ha="center", fontweight="bold")

    def _graph_never(self) -> None:
        """The same graph driven absurdly hard — still embedded."""
        ax = self.ax_graph
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        self._annulus(ax)
        s = self.graph
        drawn, sheet = s.bulk.wrap(s.radius(gain=6.0))
        xy = np.stack([drawn * np.cos(s.sigma), drawn * np.sin(s.sigma)], -1)
        start = 0
        for c in list(np.nonzero(np.diff(sheet) != 0)[0]) + [len(sheet) - 1]:
            seg = xy[start:c + 1]
            if len(seg) > 1:
                ax.plot(seg[:, 0], seg[:, 1],
                        color=_PAL["hot"] if sheet[start] != 0 else _PAL["wave"],
                        lw=1.5, solid_capstyle="round", zorder=15)
            start = c + 1
        hits = self_intersections(xy)
        self._square(ax, self.shells.r_outer * 1.20)
        ax.set_title("the same graph at gain 6.0", color=_PAL["text"],
                     fontsize=10.0, pad=6)
        ax.text(0.0, -self.shells.r_outer * 1.10,
                f"crossings {s.seam_crossings(gain=6.0)['unsigned']}    "
                f"self-intersections {hits}",
                color=_PAL["cool"], family="monospace", fontsize=8.4,
                ha="center")

    def _fold(self) -> None:
        """Tangential freedom: the map folds and the curve crosses itself."""
        ax = self.ax_fold
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        self._annulus(ax)
        xy = self.lag.points(gain=self.fold_gain)
        self._folded_curve(ax, xy, _PAL["fold"])
        hits = self_intersections(xy)
        folded = self.lag.is_folded(gain=self.fold_gain)
        if folded:
            j = self.lag.jacobian(gain=self.fold_gain) < 0.0
            ax.plot(xy[j, 0], xy[j, 1], ".", ms=2.6, color=_PAL["cross"],
                    zorder=20)
        ax.plot([-self.shells.r_mid], [0.0], "o", ms=4, color=_PAL["antipode"])
        # zoomed on the antipode: the fold is ~0.016 rad wide, invisible whole
        w = 0.30
        ax.set_xlim(-self.shells.r_outer - 0.03, -self.shells.r_outer + 2 * w)
        ax.set_ylim(-w, w)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.set_title(f"tangential freedom   λ = {FOLD_MIXING:.1f}   "
                     "(zoom on the antipode)",
                     color=_PAL["text"], fontsize=9.6, pad=6)
        ax.text(-self.shells.r_outer + w, -w * 0.86,
                f"{'FOLDED' if folded else 'not folded'}    "
                f"self-intersections {hits}",
                color=_PAL["cross"] if hits else _PAL["dim"],
                family="monospace", fontsize=8.8, ha="center",
                fontweight="bold" if hits else "normal")

    def _ring_history(self, t: float) -> None:
        """The ring's height as it converges — the thing that grows."""
        ax = self.ax_ring
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        d = np.linspace(0.0, math.pi, 900)
        if not self.history:
            probe = CircleSlice(pulse_width=PULSE_WIDTH, n_sigma=17,
                                n_radial=1200, gain=1.0)
            probe.reset()
            for k in range(220):
                probe.advance_to((k + 1) * self.t_end / 220)
                v = np.abs(probe.sim.field_at_distance(d))
                m = int(np.argmax(v))
                self.history.append((probe.t, float(d[m]), float(v[m])))
        shown = [h for h in self.history if h[0] <= t + 1e-9]
        u = np.abs(self.span.sim.field_at_distance(d))
        j = int(np.argmax(u))
        hd = [h[1] for h in shown]
        hh = [h[2] for h in shown]
        ax.plot(hd, hh, color=_PAL["wave"], lw=1.6)
        ax.plot([d[j]], [u[j]], "o", ms=5.5, color=_PAL["span"])
        ax.axhline(NARROW_DELTA / SPAN_GAIN, color=_PAL["span"], lw=1.1,
                   ls=(0, (4, 3)))
        ax.text(0.06, NARROW_DELTA / SPAN_GAIN + 0.03,
                "spans the gap above here", color=_PAL["span"],
                family="monospace", fontsize=7.6)
        ax.axvline(0.5 * math.pi, color=_PAL["faint"], lw=1.0)
        ax.text(0.5 * math.pi + 0.04, 0.90, "equator", color=_PAL["dim"],
                family="monospace", fontsize=7.6)
        ax.set_xlim(0.0, math.pi)
        ax.set_ylim(0.0, 1.05)
        ax.set_xticks([0.0, 0.5 * math.pi, math.pi])
        ax.set_xticklabels(["source", "π/2", "antipode"], fontsize=7.6)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["faint"])
        ax.set_title("the ring's height as it converges — it thins, then grows",
                     color=_PAL["text"], fontsize=9.6, pad=5)

    def _jacobian(self) -> None:
        """Where the map folds: ∂σ/∂σ₀ dipping below zero."""
        ax = self.ax_jac
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        jac = self.lag.jacobian(gain=self.fold_gain)
        ax.plot(self.lag.sigma0, jac, color=_PAL["fold"], lw=1.5)
        ax.axhline(0.0, color=_PAL["cross"], lw=1.2, ls=(0, (3, 3)))
        bad = jac < 0.0
        if np.any(bad):
            ax.plot(self.lag.sigma0[bad], jac[bad], ".", ms=2.4,
                    color=_PAL["cross"])
        ax.set_xlim(-math.pi, math.pi)
        ax.set_ylim(-1.6, 2.6)
        ax.set_xticks([-math.pi, 0.0, math.pi])
        ax.set_xticklabels(["antipode", "source", "antipode"], fontsize=7.6)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["faint"])
        ax.set_title("∂σ/∂σ₀ — below zero is a fold", color=_PAL["text"],
                     fontsize=9.6, pad=5)

    # ── one frame ───────────────────────────────────────────────────────────
    def draw(self, t: float) -> None:
        self.advance_to(t)
        self._spanning()
        self._graph_never()
        self._fold()
        self._ring_history(t)
        self._jacobian()
        self._caption(t)

    def _caption(self, t: float) -> None:
        for artist in list(self.fig.texts):
            artist.remove()
        f = self.fig
        f.suptitle("Geometrodynamic QED  v48  —  the ring reaches across, and "
                   "what it takes to fold",
                   color=_PAL["text"], fontsize=12.4, fontweight="bold", y=0.972)
        f.text(0.5, 0.938,
               "the gap decides reaching · the front's curvature decides "
               "folding · they are not the same knob",
               color=_PAL["dim"], fontsize=9.2, family="monospace", ha="center")
        phase = ("outbound" if t < 0.5 * math.pi
                 else "converging" if t < ANTIPODAL_TIME - 0.12
                 else "at the focus" if t < ANTIPODAL_TIME + 0.12
                 else "past it")
        f.text(0.045, 0.908, f"t = {t:5.3f}   {phase}", color=_PAL["text"],
               family="monospace", fontsize=9.4, va="top")
        f.text(0.975, 0.908,
               f"fold threshold λε = {self.fold_gain * FOLD_MIXING / 1.45:.5f}"
               f"   driven at {self.fold_gain * FOLD_MIXING:.5f}",
               color=_PAL["dim"], family="monospace", fontsize=8.6,
               va="top", ha="right")
        f.text(0.045, 0.028,
               "a height field is embedded at any gap and any gain — it is the "
               "tangential freedom that folds it, and that threshold does not "
               "know the gap exists",
               color=_PAL["dim"], family="monospace", fontsize=8.2, va="bottom")


# ════════════════════════════════════════════════════════════════════════════
def still(path: str) -> None:
    holder = RingFoldFigure()
    times = [2.30, ANTIPODAL_TIME - 0.10, ANTIPODAL_TIME + 0.30]
    labels = ["converging", "at the focus", "past it"]
    frames = []
    for t in times:
        holder.draw(t)
        holder.fig.canvas.draw()
        frames.append(np.asarray(holder.fig.canvas.buffer_rgba()).copy())
    plt.close(holder.fig)
    fig, axes = plt.subplots(len(frames), 1, figsize=(11.6, 7.4 * len(frames)),
                             facecolor=_PAL["bg"])
    for ax, img, t, lab in zip(np.atleast_1d(axes), frames, times, labels):
        ax.imshow(img)
        ax.set_axis_off()
        ax.set_title(f"t = {t:.3f}   {lab}", color=_PAL["text"], fontsize=11)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.01, hspace=0.05)
    fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def animate(save: Optional[str] = None, frames: int = 100):
    from matplotlib import animation

    holder = RingFoldFigure()

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
    ap.add_argument("--frames", type=int, default=100)
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
