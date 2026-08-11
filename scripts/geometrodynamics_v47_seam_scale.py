#!/usr/bin/env python3
"""
Geometrodynamic QED — v47: the scaling at the seam is a choice, and it shows
============================================================================

v46 folded the bulk by **translation**: a radius past `R_outer` came back as
`r − gap`. That carries a radial *offset* across unchanged — and the two
boundary circles do not have the same circumference. A feature emerging at
`R_inner` keeps its full radial height while sitting on an arc that is shorter
by `R_outer/R_inner = 1.703`, so it comes back **squashed**. The emerging wave
is not the same wave.

Two rules
─────────
`translate`   `r → r − gap`, a translation in `r`.
`conformal`   `r → r · (R_inner/R_outer)`, a translation in `ln r`.

They agree to first order in the excursion and part company as it grows.

What the choice decides
───────────────────────
* **Shape.** Crossing outward, the translate rule multiplies a feature's aspect
  ratio by exactly `1.7027`; the conformal rule leaves it at `1.0000`, because
  height and arc length shrink by the same factor. Conformal returns a faithful
  scaled copy — translate returns a caricature.

* **Whether the radius survives.** Translate sheets are arithmetically spaced,
  so going inward they march `0.74 → 0.22 → −0.30 → −0.82`: straight through the
  origin into negative radius. Subtracting a fixed `gap` has nothing to stop it.
  Conformal sheets are geometric — `0.74 → 0.435 → 0.255 → 0.150` — accumulating
  at the origin and never arriving. They pair with a multiplicative radial law
  `r = R_mid·exp(εu)`, which is positive by construction.

* **What a winding number would look like.** Take a curve that genuinely winds —
  a ramp climbing one full period as `σ` goes around, which is a logarithmic
  spiral on the conformal seam. It returns to the same point of the quotient
  **magnified by `1.7027`**, the same factor from any starting radius (spread
  `2e-16`). On the translate seam the same curve returns *displaced* by `gap`,
  and the ratio drifts with where you started (spread `0.22`) — it is not a
  scale at all. So the conformal gluing turns topological charge into an
  observable magnification, and the translate gluing hides it.

What the choice does not decide
───────────────────────────────
The winding number. Rebuilt on a conformal seam with a multiplicative radial
law — a different identification, a different sheet structure, a different
notion of size — it is **still identically zero**, at gains driving up to `274`
unsigned crossings. `ρ(σ)` comes from a single-valued function on the circle
whichever coordinate the seam translates in, so its degree is zero either way.

That is worth having: the v46 result was not an artefact of an arbitrary
scaling. A height field cannot wind, and no choice of gluing rescues it. What
the conformal rule adds is that the winding you cannot have would have been
*visible* — as a factor of `1.703` per turn.

What is on screen
─────────────────
* the same wave at the same gain, folded both ways, side by side;
* a curve that does wind, on the conformal annulus, coming back magnified;
* the emerging feature close up under each rule — squashed against faithful;
* the sheet ladder, arithmetic against geometric, with `r = 0` marked.

Honest scope
────────────
Both rules are representation choices, not derived boundary conditions: nothing
here makes the wave dynamically aware of the seam. The conformal rule is the one
that respects the geometry it is folding, which is an argument from consistency
rather than from dynamics.

Usage
─────
    python scripts/geometrodynamics_v47_seam_scale.py             # animate
    python scripts/geometrodynamics_v47_seam_scale.py --still out.png
    python scripts/geometrodynamics_v47_seam_scale.py --save out.gif
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
    winding_curve,
)

FPS = 14
PULSE_WIDTH = 0.18
GAIN = 0.62                    # ~3× the multiplicative wrap threshold

_PAL = {
    "bg": "#05070c",
    "panel": "#080b13",
    "text": "#e8ecf4",
    "dim": "#7d8798",
    "faint": "#2a3244",
    "shell": "#39445c",
    "mid": "#4a5675",
    "sheet0": "#ffd166",
    "sheet_out": "#ff5c5c",
    "sheet_in": "#5cc8ff",
    "source": "#8ef0c0",
    "antipode": "#ff8ad0",
    "spiral": "#c08cff",
    "zero": "#ff5c5c",
}


def _sheet_colour(k: int) -> str:
    if k == 0:
        return _PAL["sheet0"]
    return _PAL["sheet_out"] if k > 0 else _PAL["sheet_in"]


class SeamScaleFigure:
    """The same wave folded two ways, and what the difference costs."""

    def __init__(self, figsize=(13.4, 8.6), t_end: float = RETURN_TIME) -> None:
        self.t_end = float(t_end)
        # One wave, one radial law; only the seam rule differs between panels.
        self.slices = {
            mode: CircleSlice(bulk=BulkAnnulus(mode=mode),
                              radial_law="multiplicative",
                              pulse_width=PULSE_WIDTH, n_sigma=1441, gain=GAIN)
            for mode in ("translate", "conformal")
        }
        self.b = self.slices["conformal"].bulk

        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 3, height_ratios=[1.0, 0.72],
            left=0.035, right=0.972, top=0.858, bottom=0.135,
            wspace=0.16, hspace=0.28)
        self.ax = {
            "translate": self.fig.add_subplot(gs[0, 0], facecolor=_PAL["bg"]),
            "conformal": self.fig.add_subplot(gs[0, 1], facecolor=_PAL["bg"]),
            "spiral": self.fig.add_subplot(gs[0, 2], facecolor=_PAL["bg"]),
            "zoom_t": self.fig.add_subplot(gs[1, 0], facecolor=_PAL["bg"]),
            "zoom_c": self.fig.add_subplot(gs[1, 1], facecolor=_PAL["bg"]),
            "ladder": self.fig.add_subplot(gs[1, 2], facecolor=_PAL["panel"]),
        }

    def advance_to(self, t: float) -> None:
        for s in self.slices.values():
            if t < s.t:
                s.reset()
            s.advance_to(t)

    # ── furniture ───────────────────────────────────────────────────────────
    def _annulus(self, ax) -> None:
        a = np.linspace(0.0, TWO_PI, 361)
        for r, colour, ls in ((self.b.r_outer, _PAL["shell"], "-"),
                              (self.b.r_inner, _PAL["shell"], "-"),
                              (self.b.r_mid, _PAL["mid"], (0, (4, 4)))):
            ax.plot(r * np.cos(a), r * np.sin(a), color=colour, lw=1.2, ls=ls,
                    alpha=0.9, zorder=5)

    def _curve(self, ax, s: CircleSlice, lw: float = 2.0) -> None:
        sheet = s.sheet(gain=GAIN)
        drawn, _ = s.bulk.wrap(s.radius(gain=GAIN))
        xy = np.stack([drawn * np.cos(s.sigma), drawn * np.sin(s.sigma)], -1)
        start = 0
        for c in list(np.nonzero(np.diff(sheet) != 0)[0]) + [len(sheet) - 1]:
            seg = xy[start:c + 1]
            if len(seg) > 1:
                ax.plot(seg[:, 0], seg[:, 1], color=_sheet_colour(sheet[start]),
                        lw=lw, solid_capstyle="round", zorder=15)
            start = c + 1
        ax.plot([self.b.r_mid], [0.0], "o", ms=4, color=_PAL["source"], zorder=20)
        ax.plot([-self.b.r_mid], [0.0], "o", ms=4, color=_PAL["antipode"],
                zorder=20)

    def _square(self, ax, lim: float) -> None:
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_axis_off()

    # ── panels ──────────────────────────────────────────────────────────────
    def _fold(self, mode: str, title: str) -> None:
        ax = self.ax[mode]
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        self._annulus(ax)
        self._curve(ax, self.slices[mode])
        self._square(ax, self.b.r_outer * 1.18)
        ax.set_title(title, color=_PAL["text"], fontsize=10.2, pad=6)

    def _spiral(self) -> None:
        ax = self.ax["spiral"]
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        self._annulus(ax)
        b = BulkAnnulus(mode="conformal")
        sigma, r = winding_curve(b, turns=1, n=1441, r_start=b.r_inner * 1.001)
        drawn, sheet = b.wrap(r)
        xy = np.stack([drawn * np.cos(sigma), drawn * np.sin(sigma)], -1)
        ax.plot(xy[:, 0], xy[:, 1], color=_PAL["spiral"], lw=2.2,
                solid_capstyle="round", zorder=15)
        ax.plot([xy[0, 0]], [xy[0, 1]], "o", ms=5, color=_PAL["spiral"])
        ax.annotate("", xy=xy[-1], xytext=xy[-40],
                    arrowprops=dict(arrowstyle="-|>", color=_PAL["spiral"],
                                    lw=1.6))
        self._square(ax, self.b.r_outer * 1.18)
        ax.set_title("what winding would look like  (w = 1)",
                     color=_PAL["text"], fontsize=10.2, pad=6)
        ax.text(0.0, -self.b.r_outer * 1.10,
                f"returns magnified ×{b.r_outer / b.r_inner:.3f}\n"
                "— and it is not a graph",
                color=_PAL["spiral"], family="monospace", fontsize=8.0,
                ha="center", va="top")

    def _zoom(self, key: str, mode: str, title: str) -> None:
        """The emerging feature at R_inner: squashed, or a faithful copy."""
        ax = self.ax[key]
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        s = self.slices[mode]
        a = np.linspace(0.0, TWO_PI, 361)
        for r, colour, ls in ((self.b.r_inner, _PAL["shell"], "-"),
                              (self.b.r_mid, _PAL["mid"], (0, (4, 4)))):
            ax.plot(r * np.cos(a), r * np.sin(a), color=colour, lw=1.1, ls=ls,
                    alpha=0.9)
        sheet = s.sheet(gain=GAIN)
        drawn, _ = s.bulk.wrap(s.radius(gain=GAIN))
        xy = np.stack([drawn * np.cos(s.sigma), drawn * np.sin(s.sigma)], -1)
        start = 0
        for c in list(np.nonzero(np.diff(sheet) != 0)[0]) + [len(sheet) - 1]:
            seg = xy[start:c + 1]
            if len(seg) > 1:
                ax.plot(seg[:, 0], seg[:, 1], color=_sheet_colour(sheet[start]),
                        lw=2.6, solid_capstyle="round")
            start = c + 1
        ax.plot([-self.b.r_mid], [0.0], "o", ms=4, color=_PAL["antipode"])
        w = 0.46
        ax.set_xlim(-self.b.r_mid - 0.06, -self.b.r_mid + 2 * w - 0.06)
        ax.set_ylim(-w, w)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.set_title(title, color=_PAL["text"], fontsize=9.6, pad=5)

    def _ladder(self) -> None:
        """Where the sheets sit going inward: arithmetic vs geometric."""
        ax = self.ax["ladder"]
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        for i, (mode, colour) in enumerate((("translate", _PAL["sheet_out"]),
                                            ("conformal", _PAL["sheet_in"]))):
            edges = BulkAnnulus(mode=mode).sheet_edges(n=4)["inward"]
            y = 1.0 - i
            ax.plot(edges, [y] * len(edges), "o", ms=6, color=colour)
            ax.plot([self.b.r_inner] + edges, [y] * (len(edges) + 1),
                    color=colour, lw=1.1, alpha=0.55)
            for k, e in enumerate(edges):
                dy = (9, 20) if i == 0 else (-15, -26)
                ax.annotate(f"{e:.2f}", (e, y), textcoords="offset points",
                            xytext=(0, dy[k % 2]), ha="center",
                            color=colour, family="monospace", fontsize=7.4)
            ax.text(1.30, y, mode, color=colour, family="monospace",
                    fontsize=8.4, va="center")
        ax.axvline(0.0, color=_PAL["zero"], lw=1.4, ls=(0, (3, 3)))
        ax.text(0.0, 1.62, "r = 0", color=_PAL["zero"], family="monospace",
                fontsize=8.2, ha="center")
        ax.axvline(self.b.r_inner, color=_PAL["shell"], lw=1.0)
        ax.set_xlim(-1.5, 1.75)
        ax.set_ylim(-0.95, 2.05)
        ax.set_yticks([])
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["faint"])
        ax.set_title("inward sheets — one rule walks through the origin",
                     color=_PAL["text"], fontsize=9.4, pad=5)

    # ── one frame ───────────────────────────────────────────────────────────
    def draw(self, t: float) -> None:
        self.advance_to(t)
        self._fold("translate", "translate:  r → r − gap")
        self._fold("conformal", "conformal:  r → r · (R_inner/R_outer)")
        self._spiral()
        self._zoom("zoom_t", "translate", "emerging feature — squashed ×1.703")
        self._zoom("zoom_c", "conformal", "emerging feature — faithful ×1.000")
        self._ladder()
        self._caption(t)

    def _caption(self, t: float) -> None:
        for artist in list(self.fig.texts):
            artist.remove()
        f = self.fig
        ct = self.slices["translate"].seam_crossings(gain=GAIN)
        cc = self.slices["conformal"].seam_crossings(gain=GAIN)
        wt = self.slices["translate"].winding_number(gain=GAIN)
        wc = self.slices["conformal"].winding_number(gain=GAIN)
        f.suptitle("Geometrodynamic QED  v47  —  the scaling at the seam is a "
                   "choice, and it shows",
                   color=_PAL["text"], fontsize=12.4, fontweight="bold", y=0.972)
        f.text(0.5, 0.938,
               "one wave, one gain, two gluings — the shape changes, the "
               "winding number does not",
               color=_PAL["dim"], fontsize=9.2, family="monospace", ha="center")
        f.text(0.035, 0.912, f"t = {t:5.3f}     gain {GAIN:.2f}     "
               "radial law  r = R_mid·exp(εu)",
               color=_PAL["text"], family="monospace", fontsize=9.2, va="top")
        f.text(0.972, 0.062,
               f"translate  crossings {ct['unsigned']}  winding {wt:+d}       "
               f"conformal  crossings {cc['unsigned']}  winding {wc:+d}",
               color=_PAL["text"], family="monospace", fontsize=9.4,
               va="bottom", ha="right")
        f.text(0.035, 0.030,
               "the conformal seam returns a scaled copy and can never reach a "
               "negative radius; it also makes a winding number visible as a "
               "magnification.\nwhat it does not do is give the wave one — a "
               "single-valued height has degree zero whichever coordinate the "
               "seam translates in.",
               color=_PAL["dim"], family="monospace", fontsize=8.2, va="bottom")


# ════════════════════════════════════════════════════════════════════════════
def still(path: str) -> None:
    holder = SeamScaleFigure()
    times = [ANTIPODAL_TIME - 0.16, ANTIPODAL_TIME + 0.45]
    labels = ["at the antipode, wrapping", "past it"]
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


def animate(save: Optional[str] = None, frames: int = 110):
    from matplotlib import animation

    holder = SeamScaleFigure()

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
