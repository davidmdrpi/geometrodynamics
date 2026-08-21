#!/usr/bin/env python3
"""
Geometrodynamic QED — v66: two waves, and where they connect inner to outer
===========================================================================

v46 put one scalar wave on the great circle through a source and its antipode,
drew it as a radial height in the vacuole, and glued R_outer to R_inner so the
radial direction is a circle. Its result was negative and sharp: the curve is a
GRAPH r = f(sigma), so its radial winding is identically zero. Every outward
crossing of the seam is paid for by an inward one. A HEIGHT FIELD CANNOT WIND,
and one wave running to its antipode never meets itself.

That is a statement about ONE wave. This is the two-wave version:

    one wave pulsing OUTWARD and one pulsing INWARD, both refocusing at the
    antipode -- do they connect, inner to outer, and where?

YES. At the antipode, on the seam, at the refocus.

And the threshold is not a new number. A single wave crosses the seam when
eps*u = gap/2; the pair spans |delta| = 2 eps |u| of the radial circle and
touches through the seam when THAT reaches gap -- the same inequality. So v46's
"the wave comes back inside the circle" and this round's "the two pulses connect
inner to outer" are one event described twice, and the two thresholds agree to
zero.

What the panels show
--------------------
**Top left - the slice, below threshold.** Both curves in the annulus: the
outward-driven wave and the inward-driven one, mirror images about R_mid. They
coincide wherever the field is quiet and open into a lens around each pulse.
Nothing connects.

**Top centre - the antipode, close up.** Where the two arms arrive, and the
lens between them at its widest.

**Top right - at threshold.** The lens has closed THROUGH THE SEAM. One curve
sits exactly at R_outer and the other exactly at R_inner, which after gluing is
the same point. Note which is which: the antipodal refocus is a RAREFACTION, so
it is the INWARD-driven wave that bulges out to R_outer and the outward-driven
one that dips to R_inner.

**Bottom right - unrolled, with the band.** The torus cut open: the seam is a
line, both curves are graphs, and the shaded band between them is the thing that
matters. Below threshold it misses the seam; at threshold it touches; above it,
the band covers the WHOLE radial circle on an arc around the antipode -- no
radius is left outside the pair. That is what two graphs can do and one cannot.

**Bottom left - the threshold, over a return period.** The gain needed to
connect, against time, for both source configurations. The cheapest connection
is always at a REFOCUS. Meeting mid-flight -- the two travelling pulses crossing
at the quarter points -- is the WORST place, because they partially cancel.

Honest scope
------------
Everything v46 put in is still put in. The crossing rule is a REPRESENTATION
choice, not a derived boundary condition; nothing makes either wave dynamically
aware of the seam or of the other wave; the field is a LINEAR scalar on a FIXED
round background, so the two waves do not interact at all -- they are drawn on
the same torus and the question is only whether their images meet. The gain is a
DISPLAY amplitude, as in v46, so every threshold quoted is in those units.

Usage
-----
    python scripts/geometrodynamics_v66_two_wave_slice.py --still out.png
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.viz.circle_slice import ANTIPODAL_TIME, RETURN_TIME
from geometrodynamics.viz.two_wave_slice import (ANTIPODAL_SOURCES, CO_LOCATED,
                                                 TwoWaveSlice)

_PAL = {
    "bg": "#010106", "panel": "#02020a", "grid": "#0e1a2a", "rule": "#1a2838",
    "text": "#e8ecf4", "dim": "#6a8aad",
    "out": "#ffc857", "inn": "#5cc8ff", "band": "#b388ff",
    "hot": "#ff6b8a", "ok": "#7cff9e",
}


def _style(ax, title, xlabel="", ylabel=""):
    ax.set_title(title, color=_PAL["text"], fontsize=9.4, pad=7,
                 family="monospace")
    ax.set_xlabel(xlabel, color=_PAL["dim"], fontsize=7.6, family="monospace")
    ax.set_ylabel(ylabel, color=_PAL["dim"], fontsize=7.6, family="monospace")
    ax.tick_params(colors=_PAL["dim"], labelsize=7.0)
    for s in ax.spines.values():
        s.set_color(_PAL["rule"])


class TwoWaveFigure:
    """The lens, the antipode, the connection, the band, and the threshold."""

    def __init__(self, figsize=(14.4, 9.0)) -> None:
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 3, left=0.055, right=0.975, top=0.815, bottom=0.095,
            wspace=0.28, hspace=0.40, height_ratios=[1.0, 0.85])
        self.ax_slice = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_zoom = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_conn = self.fig.add_subplot(gs[0, 2], facecolor=_PAL["panel"])
        self.ax_thr = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_roll = self.fig.add_subplot(gs[1, 1:], facecolor=_PAL["panel"])

        self.pair = TwoWaveSlice(offset=CO_LOCATED)
        self.pair.slice_.reset()
        self.pair.slice_.advance_to(ANTIPODAL_TIME)
        self.threshold = self.pair.contact_gain()
        self.below = 0.86 * self.threshold

    # ── helpers ─────────────────────────────────────────────────────────────
    def _rings(self, ax, zoom=False):
        b = self.pair.slice_.bulk
        th = np.linspace(0, 2 * math.pi, 721)
        for r, lab in ((b.r_inner, "R_inner"), (b.r_outer, "R_outer")):
            ax.plot(r * np.cos(th), r * np.sin(th), lw=1.0,
                    color=_PAL["dim"], alpha=0.55)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)

    def _curves(self, ax, gain, lw=1.9):
        """Both curves, the inward one dashed.

        They coincide wherever the field is quiet -- which is most of the
        circle -- so a solid pair would draw one on top of the other and the
        picture would show a single wave.
        """
        seg_a, seg_b = self.pair.segments(gain=gain)
        styles = ((seg_a, _PAL["out"], "-", lw * 1.5, "outward-driven"),
                  (seg_b, _PAL["inn"], (0, (4, 3)), lw, "inward-driven"))
        for segs, col, ls, width, lab in styles:
            for i, p in enumerate(segs):
                ax.plot(p[:, 0], p[:, 1], lw=width, color=col, ls=ls,
                        solid_capstyle="round", label=lab if i == 0 else None)

    # ── panels ──────────────────────────────────────────────────────────────
    def _slice(self):
        ax = self.ax_slice
        self._rings(ax)
        self._curves(ax, self.below)
        b = self.pair.slice_.bulk
        ax.plot([b.r_outer * 1.0], [0.0], "o", ms=5, color=_PAL["ok"])
        ax.text(b.r_outer * 1.06, 0.02, "source", color=_PAL["ok"],
                fontsize=7.0, family="monospace")
        ax.plot([-b.r_outer], [0.0], "o", ms=5, color=_PAL["hot"])
        ax.text(-b.r_outer * 1.52, 0.02, "antipode", color=_PAL["hot"],
                fontsize=7.0, family="monospace")
        _style(ax, "below threshold — a lens, and no connection")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.6, loc="lower center",
                  ncol=2, framealpha=0.9)
        ax.text(0.5, 0.985, f"gain {self.below:.3f}  =  0.86 x threshold",
                transform=ax.transAxes, ha="center", va="top",
                color=_PAL["dim"], fontsize=6.8, family="monospace")

    def _zoom(self):
        ax = self.ax_zoom
        self._rings(ax)
        self._curves(ax, self.below, lw=2.4)
        b = self.pair.slice_.bulk
        ax.set_xlim(-b.r_outer * 1.12, -b.r_inner * 0.72)
        ax.set_ylim(-0.34, 0.34)
        ax.plot([-b.r_outer], [0.0], "o", ms=5, color=_PAL["hot"])
        _style(ax, "the antipode, close up")
        ax.text(0.5, 0.055, "the two arms arrive; the lens is widest here",
                transform=ax.transAxes, ha="center", color=_PAL["dim"],
                fontsize=6.8, family="monospace",
                bbox=dict(facecolor=_PAL["panel"], edgecolor="none",
                          alpha=0.85, pad=2.0))

    def _connection(self):
        ax = self.ax_conn
        self._rings(ax)
        self._curves(ax, self.threshold, lw=2.2)
        b = self.pair.slice_.bulk
        got = self.pair.contact(gain=self.threshold)
        sg = got["sigma_of_closest_approach"]
        for r, col in ((got["radius_a"], _PAL["out"]),
                       (got["radius_b"], _PAL["inn"])):
            ax.plot([r * math.cos(sg)], [r * math.sin(sg)], "o", ms=9,
                    color=col, mec=_PAL["text"], mew=0.8, zorder=6)
        ax.set_xlim(-b.r_outer * 1.30, b.r_outer * 1.30)
        ax.set_ylim(-b.r_outer * 1.30, b.r_outer * 1.30)
        _style(ax, "at threshold — connected, through the seam")
        ax.text(0.5, 0.985,
                f"gain {self.threshold:.3f}   sigma = {got['sigma_over_pi']:+.0f} pi",
                transform=ax.transAxes, ha="center", va="top",
                color=_PAL["text"], fontsize=6.8, family="monospace")
        ax.text(0.5, 0.045,
                f"outward-driven -> {got['radius_a']:.3f} = R_inner\n"
                f"inward-driven  -> {got['radius_b']:.3f} = R_outer\n"
                "glued, these are ONE point",
                transform=ax.transAxes, ha="center", color=_PAL["hot"],
                fontsize=6.4, family="monospace",
                bbox=dict(facecolor=_PAL["panel"], edgecolor="none",
                          alpha=0.85, pad=2.0))

    def _unrolled(self):
        """The torus cut open, recentred on the antipode.

        Centred there rather than on the source because that is where the two
        pulses are and where the band closes; on a source-centred axis the
        whole event sits split across the two edges.
        """
        ax = self.ax_roll
        p = self.pair
        b = p.slice_.bulk
        # sigma' = sigma - pi, wrapped: the antipode at 0, the source at +-pi
        raw = p.sigma
        shifted = np.mod(raw - math.pi + math.pi, 2 * math.pi) - math.pi
        order = np.argsort(shifted)
        x = shifted[order]
        gain = 1.22 * self.threshold
        r_a, r_b = (v[order] for v in p.radii(gain=gain))
        lo, hi = np.minimum(r_a, r_b), np.maximum(r_a, r_b)
        ax.fill_between(x, lo, hi, color=_PAL["band"], alpha=0.32, lw=0,
                        label="the band between them")
        ax.plot(x, r_a, lw=2.0, color=_PAL["out"], label="outward-driven")
        ax.plot(x, r_b, lw=1.7, color=_PAL["inn"], ls=(0, (4, 3)),
                label="inward-driven")
        for r, lab in ((b.r_inner, "R_inner"), (b.r_outer, "R_outer")):
            ax.axhline(r, color=_PAL["hot"], lw=1.1)
            ax.text(math.pi * 0.995, r, f"{lab} ", color=_PAL["hot"],
                    fontsize=6.6, family="monospace", va="center", ha="right")
        covered = np.abs(r_a - r_b) >= b.gap
        if covered.any():
            # drawn last, and edged, or it disappears under the band it is
            # a subset of -- the whole point is that this arc is different
            ax.fill_between(x, b.r_inner, b.r_outer, where=covered,
                            color=_PAL["ok"], alpha=0.55, lw=0, zorder=4,
                            label="band covers EVERY radius")
            lo_x = float(x[covered].min())
            hi_x = float(x[covered].max())
            for xv in (lo_x, hi_x):
                ax.axvline(xv, color=_PAL["ok"], lw=1.0, ls=":", zorder=5)
            ax.annotate("", xy=(lo_x, b.r_outer + 0.03),
                        xytext=(hi_x, b.r_outer + 0.03),
                        arrowprops=dict(arrowstyle="<->", color=_PAL["ok"],
                                        lw=1.0))
            ax.text(0.5 * (lo_x + hi_x), b.r_outer + 0.05,
                    f"the connected arc, {hi_x - lo_x:.3f} rad",
                    color=_PAL["ok"], fontsize=6.4, family="monospace",
                    ha="center", va="bottom")
        ax.set_xlim(-math.pi, math.pi)
        ax.set_ylim(min(lo.min(), b.r_inner) - 0.05,
                    max(hi.max(), b.r_outer) + 0.05)
        ax.set_xticks([-math.pi, -math.pi / 2, 0, math.pi / 2, math.pi])
        ax.set_xticklabels(["source", "", "ANTIPODE", "", "source"],
                           fontsize=7.0, family="monospace", color=_PAL["dim"])
        _style(ax, "unrolled: above threshold the band leaves no radius out",
               "", "radius")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.4, loc="upper left",
                  framealpha=0.92)
        ax.text(0.5, 0.045,
                "gain 1.22 x threshold  —  on the green arc every radius lies "
                "between the two waves;\nthe pair has spanned the whole radial "
                "circle, which one graph can never do",
                transform=ax.transAxes, ha="center", color=_PAL["text"],
                fontsize=6.6, family="monospace",
                bbox=dict(facecolor=_PAL["panel"], edgecolor="none",
                          alpha=0.85, pad=2.0))

    def _threshold(self):
        ax = self.ax_thr
        for offset, col, lab in ((CO_LOCATED, _PAL["ok"], "co-located sources"),
                                 (ANTIPODAL_SOURCES, _PAL["hot"],
                                  "antipodal sources")):
            pr = TwoWaveSlice(offset=offset)
            rows = pr.scan(samples=160)
            ts = np.array([r["t"] for r in rows]) / math.pi
            gs = np.array([r["contact_gain"] for r in rows])
            ax.plot(ts, gs, lw=1.6, color=col, label=lab)
            k = int(np.argmin(gs))
            ax.plot([ts[k]], [gs[k]], "*", ms=13, color=col,
                    mec=_PAL["text"], mew=0.6, zorder=5)
        for x in (1.0, 2.0):
            ax.axvline(x, color=_PAL["dim"], lw=0.8, ls=":")
        ax.text(1.02, 0.86, "refocus", transform=ax.get_xaxis_transform(),
                color=_PAL["dim"], fontsize=6.4, family="monospace")
        ax.axvline(0.5, color=_PAL["dim"], lw=0.8, ls="--")
        ax.text(0.52, 0.30, "mid-flight\ncrossing",
                transform=ax.get_xaxis_transform(), color=_PAL["dim"],
                fontsize=6.4, family="monospace")
        ax.set_yscale("log")
        ax.set_ylim(0.15, 12.0)
        ax.grid(True, color=_PAL["grid"], lw=0.5, alpha=0.7)
        _style(ax, "the gain needed to connect", "t / pi", "contact gain")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.6, loc="upper center",
                  framealpha=0.9)
        ax.text(0.5, 0.045,
                "stars = cheapest.  both are at a REFOCUS;\n"
                "meeting mid-flight costs 7-9x more",
                transform=ax.transAxes, ha="center", color=_PAL["text"],
                fontsize=6.6, family="monospace")

    # ── assembly ────────────────────────────────────────────────────────────
    def draw(self):
        self._slice()
        self._zoom()
        self._connection()
        self._threshold()
        self._unrolled()
        self.fig.text(0.5, 0.955,
                      "v66  -  two waves, and where they connect inner to outer",
                      color=_PAL["text"], fontsize=15.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.917,
                      "v46: one wave is a GRAPH, so it cannot wind and never "
                      "meets itself.  two graphs bound a BAND, and a band can "
                      "cover the radial circle.",
                      color=_PAL["dim"], fontsize=8.2, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.886,
                      "they connect AT THE ANTIPODE, ON THE SEAM, AT THE "
                      "REFOCUS  --  and at exactly the amplitude where one "
                      "wave would have wrapped: both are eps u = gap/2",
                      color=_PAL["ok"], fontsize=7.6, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.858,
                      "the antipodal refocus is a RAREFACTION, so it is the "
                      "inward-driven wave that reaches R_outer",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.021,
                      "the crossing rule is a REPRESENTATION choice, not a "
                      "derived boundary condition   -   LINEAR scalar on a "
                      "FIXED round background: the two waves do not interact   "
                      "-   the gain is a DISPLAY amplitude",
                      color="#3d5570", fontsize=6.8, ha="center",
                      family="monospace")


def still(path: str) -> None:
    fig = TwoWaveFigure()
    fig.draw()
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v66.png")
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
