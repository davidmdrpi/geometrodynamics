#!/usr/bin/env python3
"""
Geometrodynamic QED — v68: how ONE surface answers to two focusing fronts
=========================================================================

v67 established the object: one scalar deformation of one surface,

    u = s_A u_A + s_B u_B          r = R_mid + eps u

and there is only ever one closed curve. This is the offset question asked
*inside* that object: as the two antipodal axes move apart, where do wave A and
wave B individually sit on that same surface, what are their signs, where do
they interfere, and how does the surface answer?

The rule the whole figure obeys
-------------------------------
Only ONE thing in this figure is ever drawn as a closed curve in the annulus:
the surface. The two contributions c_A = s_A u_A and c_B = s_B u_B are
components of that one deformation, so they appear only on GRAPHS of field
value against sigma, where nothing invites reading them as separate objects.
Drawing them as closed curves is exactly the v66 error, and it is not repeated
here. What the annulus panels do instead is COLOUR the single surface by which
contribution is live along it.

What the panels show
--------------------
**Row 1 - the surface, at five offsets.** One closed curve each, coloured by
which front owns that arc: A only, B only, or both present. The two source axes
and the bisector are marked. At alpha = 0 the curve is a perfect circle,
because the two contributions cancel identically.

**Row 2 - the same five, decomposed.** c_A and c_B thin and signed, their sum
thick, over sigma. The shaded band is where both are present -- the only place
interference happens at all.

**Row 3 left - the interference map.** sigma against alpha: where both fronts
are live, and whether they share a sign there. The two source axes and the
bisector run through it, so the geometry of the overlap is visible directly.

**Row 3 centre - what the offset buys.** The overlap arc collapses to zero once
the offset exceeds the pulse width, and the amplification max|u| / max|c_A|
sits at 1.01-1.02 throughout. The offset does not turn interference ON; it
turns the CANCELLATION off, and what is left is two nearly independent dents.

**Row 3 right - the two fronts arriving.** sigma against time at alpha = pi/2:
each source launches two arms, they cross, and the surface answers at each
crossing. The refocus is where the deformation is largest.

Honest scope
------------
Unchanged. The crossing rule that glues R_outer to R_inner is a REPRESENTATION
choice, not a derived boundary condition. The field is a LINEAR scalar on a
FIXED round background, so the two contributions do not act on each other --
they add, and that is all "interference" means here. The gain is a DISPLAY
amplitude.

Usage
-----
    python scripts/geometrodynamics_v68_two_fronts.py --still v68.png
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from geometrodynamics.viz.circle_slice import RETURN_TIME, TWO_PI
from geometrodynamics.viz.one_surface import (
    OPPOSED, OneSurfaceSlice, decompose,
    measure_how_one_surface_answers_two_fronts)

_PAL = {
    "bg": "#010106", "panel": "#02020a", "grid": "#0e1a2a", "rule": "#1a2838",
    "text": "#e8ecf4", "dim": "#6a8aad",
    "a": "#ffc857", "b": "#5cc8ff", "both": "#b388ff",
    "hot": "#ff6b8a", "ok": "#7cff9e", "quiet": "#2b3a4e",
}

OFFSETS = (0.0, 0.15, 0.35, 0.65, 1.0)


def _style(ax, title, xlabel="", ylabel="", size=9.0):
    ax.set_title(title, color=_PAL["text"], fontsize=size, pad=6,
                 family="monospace")
    ax.set_xlabel(xlabel, color=_PAL["dim"], fontsize=7.2, family="monospace")
    ax.set_ylabel(ylabel, color=_PAL["dim"], fontsize=7.2, family="monospace")
    ax.tick_params(colors=_PAL["dim"], labelsize=6.6)
    for s in ax.spines.values():
        s.set_color(_PAL["rule"])


class TwoFrontsFigure:
    """One surface, two contributions, five offsets, and the overlap."""

    def __init__(self, figsize=(15.2, 12.4)) -> None:
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            3, 15, left=0.045, right=0.982, top=0.855, bottom=0.058,
            wspace=1.5, hspace=0.34, height_ratios=[1.0, 0.72, 0.86])
        self.ax_ring = [self.fig.add_subplot(gs[0, 3 * i:3 * i + 3],
                                             facecolor=_PAL["panel"])
                        for i in range(5)]
        self.ax_dec = [self.fig.add_subplot(gs[1, 3 * i:3 * i + 3],
                                            facecolor=_PAL["panel"])
                       for i in range(5)]
        self.ax_map = self.fig.add_subplot(gs[2, 0:5], facecolor=_PAL["panel"])
        self.ax_buy = self.fig.add_subplot(gs[2, 5:10], facecolor=_PAL["panel"])
        self.ax_time = self.fig.add_subplot(gs[2, 10:15],
                                            facecolor=_PAL["panel"])

        # one clock for everything, at the return focus where the fronts are
        # back on their own axes and the deformation is largest
        self.host = OneSurfaceSlice(offset=0.0, signs=OPPOSED)
        self.t_focus = 1.99 * math.pi
        self.gain = 0.34

    def _surf(self, alpha):
        q = OneSurfaceSlice(offset=alpha, signs=OPPOSED)
        q.slice_ = self.host.slice_
        return q

    # ── row 1: the one surface, coloured by which front owns each arc ───────
    def _ring(self, ax, f):
        alpha = f * math.pi
        q = self._surf(alpha)
        sl = q.slice_
        sl.reset()
        sl.advance_to(self.t_focus)
        d = decompose(q)
        b = sl.bulk
        s = q.sigma
        r = q.radius(gain=self.gain)

        th = np.linspace(0, TWO_PI, 721)
        for rr in (b.r_inner, b.r_outer):
            ax.plot(rr * np.cos(th), rr * np.sin(th), lw=1.0,
                    color=_PAL["hot"], alpha=0.55, ls=(0, (6, 5)))
        ax.plot(b.r_mid * np.cos(th), b.r_mid * np.sin(th), lw=0.7,
                color=_PAL["rule"])

        # colour each segment of the ONE curve by which contribution is live
        pa = np.abs(d["contribution_a"])
        pb = np.abs(d["contribution_b"])
        cut = 0.02 * max(float(pa.max()), float(pb.max()), 1e-30)
        live_a, live_b = pa > cut, pb > cut
        cols = np.where(live_a & live_b, _PAL["both"],
                        np.where(live_a, _PAL["a"],
                                 np.where(live_b, _PAL["b"], _PAL["quiet"])))
        xy = np.stack([r * np.cos(s), r * np.sin(s)], axis=-1)
        xy = np.concatenate([xy, xy[:1]])
        segs = np.stack([xy[:-1], xy[1:]], axis=1)
        ax.add_collection(LineCollection(segs, colors=list(cols), linewidths=2.6,
                                         capstyle="round", zorder=5))

        # the two source axes and the bisector
        coincident = abs(alpha) < 1e-9
        marks = ([(0.0, _PAL["both"], "A=B")] if coincident
                 else [(0.0, _PAL["a"], "A"), (alpha, _PAL["b"], "B")])
        for ang, col, lab in marks:
            ax.plot([0, b.r_outer * 1.30 * math.cos(ang)],
                    [0, b.r_outer * 1.30 * math.sin(ang)], lw=0.8,
                    color=col, alpha=0.5, ls=":")
            ax.text(b.r_outer * 1.42 * math.cos(ang),
                    b.r_outer * 1.42 * math.sin(ang), lab, color=col,
                    fontsize=8.4, family="monospace", ha="center",
                    va="center", weight="bold")
        bis = 0.5 * alpha
        ax.plot([0, b.r_outer * 1.22 * math.cos(bis)],
                [0, b.r_outer * 1.22 * math.sin(bis)], lw=1.1,
                color=_PAL["both"], ls=(0, (4, 4)), alpha=0.85)

        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_xlim(-b.r_outer * 1.52, b.r_outer * 1.52)
        ax.set_ylim(-b.r_outer * 1.52, b.r_outer * 1.52)
        tag = ("alpha = 0   (they cancel)" if f == 0
               else f"alpha = {f:.2f} pi")
        _style(ax, tag, size=8.6)
        ax.text(0.5, 0.02, f"max|u| {d['peak_field']:.3f}   "
                           f"overlap {d['overlap_arc']:.3f} rad",
                transform=ax.transAxes, ha="center", color=_PAL["dim"],
                fontsize=6.2, family="monospace")
        sl.reset()

    # ── row 2: the same five, as a decomposition ────────────────────────────
    def _dec(self, ax, f, first=False):
        alpha = f * math.pi
        q = self._surf(alpha)
        sl = q.slice_
        sl.reset()
        sl.advance_to(self.t_focus)
        d = decompose(q)
        x = q.sigma / math.pi
        order = np.argsort(x)
        x = x[order]
        ca = d["contribution_a"][order]
        cb = d["contribution_b"][order]
        u = d["field"][order]
        live = d["both_present"][order]

        ax.axhline(0, color=_PAL["rule"], lw=0.8)
        if live.any():
            ax.fill_between(x, -1.4, 1.4, where=live, color=_PAL["both"],
                            alpha=0.22, lw=0)
        ax.plot(x, ca, lw=1.2, color=_PAL["a"], label="c_A")
        ax.plot(x, cb, lw=1.2, color=_PAL["b"], label="c_B")
        ax.plot(x, u, lw=2.2, color=_PAL["ok"], label="u = c_A + c_B")
        for ang, col in ((0.0, _PAL["a"]), (f, _PAL["b"])):
            ax.axvline(ang, color=col, lw=0.7, ls=":", alpha=0.6)
        ax.axvline(0.5 * f, color=_PAL["both"], lw=1.0, ls=(0, (4, 4)),
                   alpha=0.8)
        ax.set_ylim(-1.4, 1.4)
        ax.set_xlim(-1, 1)
        ax.set_xticks([-1, 0, 1])
        _style(ax, "", "sigma / pi", "field" if first else "", size=8.0)
        if not first:
            ax.set_yticklabels([])
        if first:
            ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                      labelcolor=_PAL["text"], fontsize=5.8, loc="upper left",
                      framealpha=0.95)
        ax.text(0.97, 0.94,
                f"x{d['amplification']:.3f}" if f else "u = 0 EXACTLY",
                transform=ax.transAxes, ha="right", va="top",
                color=_PAL["ok"] if f else _PAL["hot"], fontsize=7.0,
                family="monospace")
        sl.reset()

    # ── row 3 ───────────────────────────────────────────────────────────────
    def _map(self):
        ax = self.ax_map
        n_a = 260
        afs = np.linspace(0.0, 1.0, n_a)
        n_s = len(self.host.sigma)   # the circle, endpoint dropped
        img = np.zeros((n_a, n_s))
        sl = self.host.slice_
        for i, f in enumerate(afs):
            q = self._surf(f * math.pi)
            sl.reset()
            sl.advance_to(self.t_focus)
            d = decompose(q)
            pa, pb = np.abs(d["contribution_a"]), np.abs(d["contribution_b"])
            cut = 0.02 * max(float(pa.max()), float(pb.max()), 1e-30)
            la, lb = pa > cut, pb > cut
            row = np.where(la & lb, np.where(d["reinforcing"], 2.0, -2.0),
                           np.where(la, 1.0, np.where(lb, -1.0, 0.0)))
            img[i] = row[np.argsort(q.sigma)]
        sl.reset()
        from matplotlib.colors import ListedColormap, BoundaryNorm
        cmap = ListedColormap([_PAL["both"], _PAL["b"], _PAL["quiet"],
                               _PAL["a"], _PAL["ok"]])
        norm = BoundaryNorm([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5], cmap.N)
        ax.imshow(img, aspect="auto", origin="lower", cmap=cmap, norm=norm,
                  extent=[-1, 1, 0, 1], interpolation="nearest")
        ax.plot([0, 0], [0, 1], color=_PAL["a"], lw=1.0, ls=":")
        ax.plot(afs, afs, color=_PAL["b"], lw=1.0, ls=":")
        ax.plot(0.5 * afs, afs, color=_PAL["text"], lw=1.2, ls=(0, (4, 4)))
        ax.text(-0.06, 0.90, "A axis", color=_PAL["a"], fontsize=6.4,
                family="monospace", ha="right")
        ax.text(0.80, 0.90, "B axis", color=_PAL["b"], fontsize=6.4,
                family="monospace", ha="right")
        ax.text(0.44, 0.72, "bisector", color=_PAL["text"], fontsize=6.4,
                family="monospace", ha="left")
        _style(ax, "where each front is live, against offset",
               "sigma / pi", "offset alpha / pi")
        ax.text(0.5, 0.045,
                "yellow = A only   blue = B only   green = both, reinforcing\n"
                "the overlap closes as soon as the offset clears the pulse",
                transform=ax.transAxes, ha="center", color=_PAL["text"],
                fontsize=6.4, family="monospace",
                bbox=dict(facecolor=_PAL["panel"], edgecolor="none",
                          alpha=0.85, pad=2.0))

    def _buy(self):
        ax = self.ax_buy
        got = measure_how_one_surface_answers_two_fronts(
            offsets=tuple(np.linspace(0.0, 1.0, 21)), samples=120)
        a = [r["offset_over_pi"] for r in got["rows"]]
        amp = [r["amplification"] for r in got["rows"]]
        ov = [r["overlap_arc"] for r in got["rows"]]
        ax.plot(a, amp, "o-", lw=1.8, ms=4, color=_PAL["ok"],
                mec=_PAL["text"], mew=0.4, label="max|u| / max|c_A|")
        ax.axhline(1.0, color=_PAL["rule"], lw=1.0, ls=":")
        ax2 = ax.twinx()
        ax2.plot(a, ov, "s--", lw=1.5, ms=4, color=_PAL["both"],
                 mec=_PAL["text"], mew=0.4, label="overlap arc (rad)")
        ax2.set_ylabel("overlap arc (rad)", color=_PAL["both"], fontsize=6.8,
                       family="monospace", labelpad=1)
        ax2.tick_params(colors=_PAL["both"], labelsize=6.6)
        for sp in ax2.spines.values():
            sp.set_color(_PAL["rule"])
        w = self.host.slice_.sim.pulse_width / math.pi
        ax.axvspan(0, w, color=_PAL["hot"], alpha=0.15)
        ax.text(w * 1.2, 0.5, "one pulse width",
                transform=ax.get_xaxis_transform(), color=_PAL["hot"],
                fontsize=6.2, family="monospace")
        ax.set_ylim(-0.05, 1.25)
        ax.grid(True, color=_PAL["grid"], lw=0.5, alpha=0.7)
        _style(ax, "what the offset actually buys", "offset alpha / pi",
               "amplification")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, facecolor=_PAL["panel"],
                  edgecolor=_PAL["rule"], labelcolor=_PAL["text"],
                  fontsize=6.2, loc="center right", framealpha=0.95)
        lo, hi = got["amplification_when_apart"]
        ax.text(0.52, 0.30,
                f"apart, the total is {lo:.3f}-{hi:.3f} x ONE contribution.\n"
                "the offset does not turn interference on -- it turns the\n"
                "CANCELLATION off, leaving two near-independent dents",
                transform=ax.transAxes, ha="center", color=_PAL["text"],
                fontsize=6.4, family="monospace",
                bbox=dict(facecolor=_PAL["panel"], edgecolor="none",
                          alpha=0.9, pad=2.0))

    def _time(self):
        ax = self.ax_time
        alpha = 0.5 * math.pi
        q = self._surf(alpha)
        sl = q.slice_
        n_t = 260
        img = np.zeros((n_t, len(q.sigma)))
        sl.reset()
        order = np.argsort(q.sigma)
        for i in range(n_t):
            sl.advance_to((i + 1) * RETURN_TIME / n_t)
            img[i] = q.field()[order]
        sl.reset()
        lim = float(np.max(np.abs(img)))
        ax.imshow(img, aspect="auto", origin="lower", cmap="RdBu_r",
                  vmin=-lim, vmax=lim, extent=[-1, 1, 0, 2],
                  interpolation="nearest")
        ax.set_ylim(0, 2.14)
        for ang, col, lab in ((0.0, _PAL["a"], "A"), (0.5, _PAL["b"], "B")):
            ax.axvline(ang, color=col, lw=1.0, ls=":")
            ax.text(ang + 0.03, 1.84, lab, color=col, fontsize=7.4,
                    family="monospace", weight="bold")
        ax.axvline(0.25, color=_PAL["text"], lw=1.2, ls=(0, (4, 4)))
        ax.text(0.28, 1.62, "bisector", color=_PAL["text"], fontsize=6.4,
                family="monospace")
        for tt in (1.0, 2.0):
            ax.axhline(tt, color=_PAL["text"], lw=0.8, ls=":", alpha=0.6)
        ax.text(-0.96, 1.04, "antipodal refocus", color=_PAL["text"],
                fontsize=6.4, family="monospace")
        _style(ax, "the two fronts arriving  (alpha = pi/2)",
               "sigma / pi", "t / pi")
        ax.text(0.5, 0.045,
                "each source launches two arms; the surface answers\n"
                "wherever a front is, and most at a refocus",
                transform=ax.transAxes, ha="center", color=_PAL["text"],
                fontsize=6.4, family="monospace",
                bbox=dict(facecolor=_PAL["panel"], edgecolor="none",
                          alpha=0.85, pad=2.0))

    # ── assembly ────────────────────────────────────────────────────────────
    def draw(self):
        for i, f in enumerate(OFFSETS):
            self._ring(self.ax_ring[i], f)
            self._dec(self.ax_dec[i], f, first=(i == 0))
        self._map()
        self._buy()
        self._time()
        self.fig.text(0.5, 0.972,
                      "v68  -  how ONE surface answers to two focusing fronts",
                      color=_PAL["text"], fontsize=15.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.948,
                      "only ONE thing here is ever a closed curve in the "
                      "annulus: the surface.  c_A and c_B are COMPONENTS of "
                      "its deformation, so they appear only on graphs of field "
                      "against sigma.",
                      color=_PAL["dim"], fontsize=8.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.925,
                      "ROW 1  the one surface, coloured by which front owns "
                      "each arc      ROW 2  the same five, decomposed      "
                      "ROW 3  the overlap, its price, and the fronts in time",
                      color=_PAL["ok"], fontsize=7.6, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.900,
                      "at alpha = 0 the surface is a PERFECT CIRCLE: the two "
                      "contributions cancel identically, so there is nothing "
                      "for the bulk to connect",
                      color=_PAL["hot"], fontsize=7.6, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.594,
                      "ROW 1: an INWARD dent is a negative contribution and an "
                      "OUTWARD one positive -- so the sign of each front is "
                      "read straight off the surface",
                      color=_PAL["dim"], fontsize=7.2, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.877,
                      "past one pulse width the contributions stop "
                      "overlapping entirely and the total is 1.01-1.02x ONE "
                      "of them  --  two near-independent dents, not a "
                      "reinforced one",
                      color=_PAL["dim"], fontsize=7.2, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.012,
                      "the crossing rule is a REPRESENTATION choice, not a "
                      "derived boundary condition   -   LINEAR scalar on a "
                      "FIXED round background: the contributions ADD, and that "
                      "is all 'interference' means here   -   the gain is a "
                      "DISPLAY amplitude",
                      color="#3d5570", fontsize=6.6, ha="center",
                      family="monospace")


def still(path: str) -> None:
    fig = TwoFrontsFigure()
    fig.draw()
    fig.fig.savefig(path, dpi=108, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v68.png")
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
