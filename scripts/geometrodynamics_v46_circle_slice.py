#!/usr/bin/env python3
"""
Geometrodynamic QED — v46: a circle slice, and a bulk the wave can wrap through
===============================================================================

Everything so far has been the whole sphere. This cuts it with the great circle
through the source and its antipode, and watches the cross-section instead.

The slice
─────────
Parametrised by `σ ∈ [−π, π)`, the geodesic distance from the source along that
circle is just `d = |σ|` — so one circle carries **both** halves of the wave:
two lobes running in opposite directions, meeting head-on at `σ = ±π`, which is
the antipode. Nothing is re-solved. The field is the 2-D solve sampled at
`d(σ)` (verified to `1.4e-14`), so the slice inherits the real caustic rather
than the `2×` superposition a 1-D wave on a circle would give.

The bulk
────────
The slice lives in the vacuole — the annulus between `R_inner` and `R_outer`.
The crossing rule is the obvious one: a radius that would pass `R_outer`
re-enters at `R_inner`. So the wave that reaches up into the bulk **comes back
inside the circle**. That glues the two boundaries, makes the radial coordinate
periodic with period `gap = R_outer − R_inner`, and turns the space the curve
lives on into a torus `S¹_σ × S¹_ρ`.

What the bulk gives, and what it does not
─────────────────────────────────────────
On a torus a closed curve has integer winding numbers, and integers are exactly
the stable objects a crossing rule is meant to produce.

This rule does not produce one, and the reason is worth having on screen. The
curve is a **graph** `r = f(σ)` with `f` single-valued, so its radial winding
number is identically zero — at every amplitude, at every time. Every outward
crossing of the seam is paid for by an inward one: the signed count is `0` while
the unsigned count climbs `0 → 2 → 4 → 6` with the gain.

So the seam can be crossed as often as you like without accumulating any
topological charge. **A height field cannot wind.** That is a real constraint on
a crossing rule of this kind, and it says where a stable topological object
would have to come from instead — not from the amplitude of a scalar height,
but from a curve free to stop being a graph.

What separates different waves
──────────────────────────────
Driven at one common gain, pulses from `0.36` down to `0.08` all cross the seam
the same number of times — the launch amplitude is `1` for all of them and their
peaks barely differ. What varies is **how much of the circle rides the far
sheet**: `0.155 → 0.033` of the circumference, a `4.7×` spread. Divided by the
pulse width that is `2.61` for all of them — the far-sheet arc simply *is* the
pulse. None of them winds.

What is on screen
─────────────────
* the **slice** in its annulus, at a gain below threshold — two lobes leaving
  the source and running to the antipode, stretching the circle as they go;
* a **zoom** on the antipode, where the two arms arrive;
* the **same wave driven past threshold**, wrapping through the seam and
  reappearing inside the circle, coloured by which sheet it is on;
* the **unrolled torus** `(σ, ρ)`, where the seam is a line and the crossings
  are countable, with the running ledger.

Honest scope
────────────
The crossing rule is a *representation* choice, not a derived boundary
condition: nothing here makes the wave dynamically aware of the seam, and the
field is a linear scalar on a fixed round background. The zero winding number is
a fact about graphs, so it constrains this construction and any other that draws
the bulk excursion as a height.

Usage
─────
    python scripts/geometrodynamics_v46_circle_slice.py             # animate
    python scripts/geometrodynamics_v46_circle_slice.py --still out.png
    python scripts/geometrodynamics_v46_circle_slice.py --waves out.png
    python scripts/geometrodynamics_v46_circle_slice.py --save out.gif
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
    CircleSlice,
)

FPS = 14
PULSE_WIDTH = 0.18
OVER_THRESHOLD = 3.6          # how hard the wrapped panel is driven

_PAL = {
    "bg": "#05070c",
    "panel": "#080b13",
    "text": "#e8ecf4",
    "dim": "#7d8798",
    "faint": "#2a3244",
    "shell": "#39445c",
}
_PAL.update({
    "mid": "#4a5675",
    "sheet0": "#ffd166",       # the wave on its own sheet
    "sheet_out": "#ff5c5c",    # ...having wrapped outward through R_outer
    "sheet_in": "#5cc8ff",     # ...having wrapped inward through R_inner
    "source": "#8ef0c0",
    "antipode": "#ff8ad0",
    "seam": "#ff5c5c",
})


def _sheet_colour(k: int) -> str:
    if k == 0:
        return _PAL["sheet0"]
    return _PAL["sheet_out"] if k > 0 else _PAL["sheet_in"]


class SliceFigure:
    """The slice, a zoom on the antipode, the wrapped version, and the torus."""

    def __init__(self, figsize=(13.0, 8.2), t_end: float = RETURN_TIME) -> None:
        self.t_end = float(t_end)
        self.s = CircleSlice(pulse_width=PULSE_WIDTH, n_sigma=1441)
        self.b = self.s.bulk

        peak = self._run_peak()
        self.threshold = (self.b.r_outer - self.b.r_mid) / max(peak, 1e-30)
        self.hot_gain = OVER_THRESHOLD * self.threshold

        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 3, width_ratios=[1.35, 1.0, 1.0], height_ratios=[1.0, 0.82],
            left=0.035, right=0.975, top=0.885, bottom=0.135,
            wspace=0.20, hspace=0.24)
        self.ax_main = self.fig.add_subplot(gs[:, 0], facecolor=_PAL["bg"])
        self.ax_zoom = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["bg"])
        self.ax_wrap = self.fig.add_subplot(gs[0, 2], facecolor=_PAL["bg"])
        self.ax_tor = self.fig.add_subplot(gs[1, 1:], facecolor=_PAL["panel"])

    def _run_peak(self, samples: int = 200) -> float:
        self.s.reset()
        peak = 0.0
        for i in range(samples):
            self.s.advance_to((i + 1) * self.t_end / samples)
            peak = max(peak, float(np.max(np.abs(self.s.field()))))
        self.s.reset()
        return peak

    # ── shared furniture ────────────────────────────────────────────────────
    def _annulus(self, ax, caption: str = "") -> None:
        a = np.linspace(0.0, TWO_PI, 361)
        ca, sa = np.cos(a), np.sin(a)
        for r, colour, lw, al in ((self.b.r_outer, _PAL["shell"], 1.3, 0.95),
                                  (self.b.r_inner, _PAL["shell"], 1.3, 0.95),
                                  (self.b.r_mid, _PAL["mid"], 0.8, 0.55)):
            ax.plot(r * ca, r * sa, color=colour, lw=lw, alpha=al,
                    ls="-" if r != self.b.r_mid else (0, (4, 4)), zorder=5)
        if caption:
            # in data coordinates, so the label sits on the circle rather than
            # floating at the top of a tall axes box
            ax.text(0.0, self.b.r_outer + 0.07, caption, color=_PAL["text"],
                    fontsize=10.0, ha="center", va="bottom")

    def _poles(self, ax) -> None:
        ax.plot([self.b.r_mid], [0.0], marker="o", ms=4.5,
                color=_PAL["source"], zorder=20)
        ax.plot([-self.b.r_mid], [0.0], marker="o", ms=4.5,
                color=_PAL["antipode"], zorder=20)

    def _curve(self, ax, gain: float, lw: float = 2.0) -> None:
        """The slice, split at the seam and coloured by which sheet it is on."""
        sheet = self.s.sheet(gain=gain)
        drawn, _ = self.b.wrap(self.s.radius(gain=gain))
        xy = np.stack([drawn * np.cos(self.s.sigma),
                       drawn * np.sin(self.s.sigma)], axis=-1)
        cuts = list(np.nonzero(np.diff(sheet) != 0)[0])
        start = 0
        for c in cuts + [len(sheet) - 1]:
            seg = xy[start:c + 1]
            if len(seg) > 1:
                ax.plot(seg[:, 0], seg[:, 1], color=_sheet_colour(sheet[start]),
                        lw=lw, solid_capstyle="round", zorder=15)
            start = c + 1

    # ── the four panels ─────────────────────────────────────────────────────
    def _main(self) -> None:
        ax = self.ax_main
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        self._annulus(ax, "the slice — one circle, both halves of the wave")
        self._curve(ax, self.s.gain, lw=2.1)
        self._poles(ax)
        lim = self.b.r_outer * 1.22
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.text(self.b.r_mid + 0.10, 0.055, "source", color=_PAL["source"],
                family="monospace", fontsize=8)
        ax.text(-self.b.r_mid - 0.10, 0.055, "antipode", color=_PAL["antipode"],
                family="monospace", fontsize=8, ha="right")

    def _zoom(self) -> None:
        ax = self.ax_zoom
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        self._annulus(ax)
        self._curve(ax, self.s.gain, lw=2.4)
        self._poles(ax)
        x0, x1 = -self.b.r_outer - 0.05, -self.b.r_inner + 0.28
        ax.set_xlim(x0, x1)
        ax.set_ylim(-0.5 * (x1 - x0), 0.5 * (x1 - x0))
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.set_title("the antipode, close up", color=_PAL["text"],
                     fontsize=10.0, pad=6)

    def _wrapped(self) -> None:
        ax = self.ax_wrap
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        self._annulus(ax)
        self._curve(ax, self.hot_gain, lw=2.0)
        self._poles(ax)
        lim = self.b.r_outer * 1.24
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.set_title(f"the same wave at {OVER_THRESHOLD:.1f}× threshold",
                     color=_PAL["text"], fontsize=10.0, pad=6)
        ax.text(0.0, -lim * 0.90,
                "red = wrapped out through R_outer",
                color=_PAL["sheet_out"], family="monospace", fontsize=7.6,
                ha="center")

    def _torus(self) -> None:
        ax = self.ax_tor
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        sigma, rho = self.s.unrolled(gain=self.hot_gain)
        sheet = self.s.sheet(gain=self.hot_gain)
        for y in (0.0, 1.0):
            ax.axhline(y, color=_PAL["seam"], lw=1.1, alpha=0.75)
        ax.axhline((self.b.r_mid - self.b.r_inner) / self.b.gap,
                   color=_PAL["mid"], lw=0.8, ls=(0, (4, 4)), alpha=0.6)
        cuts = list(np.nonzero(np.diff(sheet) != 0)[0])
        start = 0
        for c in cuts + [len(sheet) - 1]:
            sl = slice(start, c + 1)
            if c + 1 - start > 1:
                ax.plot(sigma[sl], rho[sl], color=_sheet_colour(sheet[start]),
                        lw=1.9, solid_capstyle="round")
            start = c + 1
        ax.set_xlim(-math.pi, math.pi)
        ax.set_ylim(-0.06, 1.06)
        ax.set_xticks([-math.pi, 0.0, math.pi])
        ax.set_xticklabels(["antipode", "source", "antipode"], fontsize=7.5)
        ax.set_yticks([0.0, 1.0])
        ax.set_yticklabels(["R_inner", "R_outer"], fontsize=7.5)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["faint"])
        ax.set_title("unrolled: the seam is a line, the crossings are countable",
                     color=_PAL["text"], fontsize=9.6, pad=5)

    # ── one frame ───────────────────────────────────────────────────────────
    def draw(self, t: float) -> None:
        if t < self.s.t:
            self.s.reset()
        self.s.advance_to(t)
        self._main()
        self._zoom()
        self._wrapped()
        self._torus()
        self._caption(t)

    def _caption(self, t: float) -> None:
        for artist in list(self.fig.texts):
            artist.remove()
        f = self.fig
        c = self.s.seam_crossings(gain=self.hot_gain)
        w = self.s.winding_number(gain=self.hot_gain)
        f.suptitle("Geometrodynamic QED  v46  —  a circle slice, and a bulk "
                   "the wave can wrap through",
                   color=_PAL["text"], fontsize=12.4, fontweight="bold", y=0.972)
        f.text(0.5, 0.938,
               "the seam is crossed in pairs — a height field cannot wind",
               color=_PAL["dim"], fontsize=9.2, family="monospace", ha="center")

        phase = ("running to the antipode" if t < ANTIPODAL_TIME - 0.25
                 else "arriving at the antipode" if t < ANTIPODAL_TIME + 0.25
                 else "past it, on the way home")
        f.text(0.035, 0.914, f"t = {t:5.3f}   {phase}", color=_PAL["text"],
               family="monospace", fontsize=9.4, va="top")
        f.text(0.975, 0.038,
               f"crossings  out {c['outward']}   in {c['inward']}   "
               f"signed {c['signed']:+d}   winding {w:+d}",
               color=_PAL["text"] if c["signed"] == 0 else _PAL["seam"],
               family="monospace", fontsize=9.6, va="bottom", ha="right")
        f.text(0.975, 0.014,
               f"display gain {self.s.gain:.3f}   wrap threshold "
               f"{self.threshold:.3f}   driven at {self.hot_gain:.3f}",
               color=_PAL["dim"], family="monospace", fontsize=8.0,
               va="bottom", ha="right")
        f.text(0.035, 0.038,
               "the crossing rule glues R_outer to R_inner, so the radial\n"
               "direction is a circle and (σ, ρ) is a torus — the wave that\n"
               "reaches into the bulk comes back inside the slice",
               color=_PAL["dim"], family="monospace", fontsize=8.0, va="bottom")


# ════════════════════════════════════════════════════════════════════════════
def still(path: str) -> None:
    holder = SliceFigure()
    times = [0.85, ANTIPODAL_TIME - 0.16, ANTIPODAL_TIME + 0.75]
    labels = ["in flight", "at the antipode", "past it"]
    frames = []
    for t in times:
        holder.draw(t)
        holder.fig.canvas.draw()
        frames.append(np.asarray(holder.fig.canvas.buffer_rgba()).copy())
    plt.close(holder.fig)
    fig, axes = plt.subplots(len(frames), 1, figsize=(11.5, 7.3 * len(frames)),
                             facecolor=_PAL["bg"])
    for ax, img, t, lab in zip(np.atleast_1d(axes), frames, times, labels):
        ax.imshow(img)
        ax.set_axis_off()
        ax.set_title(f"t = {t:.3f}   {lab}", color=_PAL["text"], fontsize=11)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.01, hspace=0.05)
    fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def waves(path: str, widths=(0.36, 0.24, 0.14, 0.08), gain: float = 0.80,
          frames: int = 200) -> None:
    """Four different scalar waves meeting the same bulk at the same gain."""
    fig, axes = plt.subplots(1, len(widths), figsize=(4.0 * len(widths), 4.9),
                             facecolor=_PAL["bg"])
    for ax, w in zip(np.atleast_1d(axes), widths):
        s = CircleSlice(pulse_width=float(w), gain=gain, n_sigma=1441)
        b = s.bulk
        s.reset()
        best_t, best_arc = 0.0, -1.0
        for i in range(frames):
            s.advance_to((i + 1) * RETURN_TIME / frames)
            arc = float(np.mean(s.sheet(gain=gain) != 0))
            if arc > best_arc:
                best_arc, best_t = arc, s.t
        s.reset()
        s.advance_to(best_t)

        ax.set_facecolor(_PAL["bg"])
        a = np.linspace(0.0, TWO_PI, 361)
        for r, colour, ls in ((b.r_outer, _PAL["shell"], "-"),
                              (b.r_inner, _PAL["shell"], "-"),
                              (b.r_mid, _PAL["mid"], (0, (4, 4)))):
            ax.plot(r * np.cos(a), r * np.sin(a), color=colour, lw=1.1, ls=ls,
                    alpha=0.85)
        sheet = s.sheet(gain=gain)
        drawn, _ = b.wrap(s.radius(gain=gain))
        xy = np.stack([drawn * np.cos(s.sigma), drawn * np.sin(s.sigma)], -1)
        start = 0
        for c in list(np.nonzero(np.diff(sheet) != 0)[0]) + [len(sheet) - 1]:
            seg = xy[start:c + 1]
            if len(seg) > 1:
                ax.plot(seg[:, 0], seg[:, 1], color=_sheet_colour(sheet[start]),
                        lw=1.8, solid_capstyle="round")
            start = c + 1
        ax.plot([b.r_mid], [0.0], "o", ms=4, color=_PAL["source"])
        ax.plot([-b.r_mid], [0.0], "o", ms=4, color=_PAL["antipode"])
        lim = b.r_outer * 1.12
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_axis_off()
        cr = s.seam_crossings(gain=gain)
        ax.set_title(f"pulse {w:.2f}\narc on the far sheet {best_arc:.3f}   "
                     f"crossings {cr['unsigned']}   winding "
                     f"{s.winding_number(gain=gain):+d}",
                     color=_PAL["text"], fontsize=9.2, pad=8)
    fig.suptitle("Four waves, one bulk, one gain — the arc differs, the "
                 "winding does not", color=_PAL["text"], fontsize=12.0,
                 fontweight="bold", y=0.985)
    fig.text(0.5, 0.045,
             "each shown at its own moment of widest far-sheet arc — that arc is "
             "2.61 × the pulse width for all of them,\nwhile the crossing count "
             "and the winding number are the same for all of them",
             color=_PAL["dim"], family="monospace", fontsize=8.6, ha="center")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.86, bottom=0.10, wspace=0.05)
    fig.savefig(path, dpi=120, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def animate(save: Optional[str] = None, frames: int = 120):
    from matplotlib import animation

    holder = SliceFigure()

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
    ap.add_argument("--waves", metavar="PNG")
    ap.add_argument("--save", metavar="FILE")
    ap.add_argument("--frames", type=int, default=120)
    a = ap.parse_args(argv)
    if a.still or a.waves:
        matplotlib.use("Agg")
        if a.still:
            still(a.still)
        if a.waves:
            waves(a.waves)
        return 0
    if a.save:
        matplotlib.use("Agg")
    animate(save=a.save, frames=a.frames)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
