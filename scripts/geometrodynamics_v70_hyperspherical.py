#!/usr/bin/env python3
"""
Geometrodynamic QED - v70: what higher dimension does to the bulk picture
=========================================================================

3-D intuition models an extra dimension as "the same sphere, with another
direction available". That is wrong in ways that matter here.

What the panels show
--------------------
**Top left - the shell measure collapses onto the equator.** sin^(n-1)(chi),
normalised, for a range of n. Broad on S^2, violently peaked by n = 50. Writing
chi = pi/2 + delta gives exp(-(n-1) delta^2/2), so the band width is ~1/sqrt(n)
-- and the measured std(chi) x sqrt(n) runs 0.9669 -> 1.000000 over n = 2..1000.
ALMOST EVERY POINT IS NEARLY 90 DEGREES FROM ANY CHOSEN POINT.

**Top right - so the antipode is non-generic.** The distribution of the angle
between two random directions, by ambient dimension. It piles up at pi/2 and the
near-antipodal tail vanishes: 3.2e-04 of pairs have alpha > 0.99 pi on S^2, and
none at all in 2e+05 samples by n = 10. Selecting x <-> -x is NOT "pairing a
point with the far one" -- it picks a vanishing-measure relation out of an
enormous nearly-orthogonal majority, and gets MORE non-generic with dimension.
(That does not make the identification correct. It removes the bland reading.)

**Bottom left - THE COLLAPSE IS f^n, NOT f.** For ds^2 + f^2 dOmega_n^2 the
transverse measure of an angular patch is f^n dOmega_n, so squeezing from F to
f0 costs (f0/F)^n. The 2-D drawing's l ~ f understates the S^3 case by a FACTOR
OF A MILLION at f0/F = 1e-03. THE ANGULAR OVERLAP CAN STAY FINITE WHILE THE
PHYSICAL OVERLAP COLLAPSES AS f0^n -- so the finite bearing is much less like
two ribbons squeezing together and much more like a vast angular configuration
space packed into an extremely small proper region. The yes/no overlap
criterion is untouched: it was always angular.

**Bottom right - orientability is a PARITY effect, and this repo uses two
quotients that are always opposite.** RP^n is orientable iff n is odd, since the
antipodal map is the restriction of -I on R^(n+1) and det(-I) = (-1)^(n+1). The
repo carries the SPATIAL quotient S^d/+- = RP^d and the TWO-BODY EXCHANGE space
(R^d\\0)/+- ~ RP^(d-1) -- one apart, hence ALWAYS of opposite parity. At d = 3
the spatial RP^3 is orientable while the exchange RP^2 is not, and it is that
RP^2 whose Pin- structure makes the throat a spinor. Raising the spatial
dimension by one SWAPS WHICH IS WHICH.

What is put in
--------------
MEASURE, ORIENTATION AND FRAMES ON ROUND SPHERES. No field equation is solved,
nothing evolves, and no throat profile appears except through the f^n weighting
the bottom-left panel shares with v69. The bottom-right panel is a statement
about what WOULD have to be re-derived if the dimension moved -- not an argument
that it should. And S^3 is not a generic sphere (SU(2), parallelizable, Hopf),
which is a standing warning against extrapolating any of this smoothly.

Usage
-----
    python scripts/geometrodynamics_v70_hyperspherical.py --still out.png
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.viz.hyperspherical import (patch_collapse,
                                                 projective_space_is_orientable,
                                                 shell_measure)

_PAL = {
    "bg": "#010106", "panel": "#02020a", "grid": "#0e1a2a", "rule": "#1a2838",
    "text": "#e8ecf4", "dim": "#6a8aad",
    "cool": "#5cc8ff", "warm": "#ffb347", "good": "#7cff9e",
    "bad": "#ff6b8a", "kern": "#b388ff", "hot": "#ff5ec7",
}
_DIMS = (2, 3, 5, 10, 50)


def _style(ax, title, xlabel, ylabel):
    ax.set_title(title, color=_PAL["text"], fontsize=9.6, pad=8,
                 family="monospace")
    ax.set_xlabel(xlabel, color=_PAL["dim"], fontsize=8, family="monospace")
    ax.set_ylabel(ylabel, color=_PAL["dim"], fontsize=8, family="monospace")
    ax.tick_params(colors=_PAL["dim"], labelsize=7.2)
    for s in ax.spines.values():
        s.set_color(_PAL["rule"])
    ax.grid(True, color=_PAL["grid"], lw=0.5, alpha=0.7)


class HypersphereFigure:
    """Concentration, non-genericity, the f^n collapse, and parity."""

    def __init__(self, figsize=(13.8, 8.9)) -> None:
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, left=0.070, right=0.975, top=0.805, bottom=0.075,
            wspace=0.24, hspace=0.44)
        self.ax_shell = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_ang = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_coll = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_par = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])

    # ── the shell measure ───────────────────────────────────────────────────
    def _shell(self) -> None:
        ax = self.ax_shell
        chi = np.linspace(1e-6, math.pi - 1e-6, 1400)
        cols = [_PAL["cool"], _PAL["good"], _PAL["warm"], _PAL["kern"],
                _PAL["hot"]]
        for n, col in zip(_DIMS, cols):
            w = np.array([shell_measure(float(c), n) for c in chi])
            w /= np.trapezoid(w, chi) if hasattr(np, "trapezoid") \
                else np.trapz(w, chi)
            ax.plot(chi / math.pi, w, lw=1.9, color=col,
                    label=f"S^{n}   width ~ {1/math.sqrt(n):.2f}")
        ax.axvline(0.5, color=_PAL["rule"], lw=0.9, ls=":")
        ax.set_xlim(0, 1)
        _style(ax, "the shell measure collapses onto the equator",
               "geodesic angle from a pole   chi / pi",
               "normalised measure  sin^(n-1)(chi)")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.8, loc="upper right",
                  framealpha=0.93)
        ax.text(0.03, 0.72,
                "chi = pi/2 + delta  ->  exp(-(n-1) delta^2 / 2)\n"
                "so the band width is ~ 1/sqrt(n)\n"
                "measured std(chi) x sqrt(n): 0.9669 -> 1.000000",
                transform=ax.transAxes, ha="left", va="top",
                color=_PAL["dim"], fontsize=6.8, family="monospace")

    # ── the angle between random directions ─────────────────────────────────
    def _angles(self) -> None:
        ax = self.ax_ang
        rng = np.random.default_rng(11)
        cols = [_PAL["cool"], _PAL["good"], _PAL["warm"], _PAL["kern"],
                _PAL["hot"]]
        for n, col in zip((3, 4, 10, 50, 500), cols):
            a = rng.normal(size=(120000, n))
            a /= np.linalg.norm(a, axis=1, keepdims=True)
            b = rng.normal(size=(120000, n))
            b /= np.linalg.norm(b, axis=1, keepdims=True)
            ang = np.arccos(np.clip(np.einsum("ij,ij->i", a, b), -1, 1))
            h, e = np.histogram(ang / math.pi, bins=220, range=(0, 1),
                                density=True)
            ax.plot(0.5 * (e[1:] + e[:-1]), h, lw=1.8, color=col,
                    label=f"ambient n = {n}")
        ax.axvline(0.5, color=_PAL["rule"], lw=0.9, ls=":")
        ax.axvspan(0.99, 1.0, color=_PAL["bad"], alpha=0.22, lw=0)
        ax.text(0.975, 0.30, "near-antipodal\n3.2e-04 at n=3,\n0 by n=10",
                transform=ax.transAxes, ha="right", va="center",
                color=_PAL["bad"], fontsize=6.6, family="monospace")
        ax.set_xlim(0, 1)
        _style(ax, "so the antipode is a vanishing-measure relation",
               "angle between two random directions   alpha / pi",
               "density")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.8, loc="upper left",
                  framealpha=0.93)

    # ── the f^n collapse ────────────────────────────────────────────────────
    def _collapse(self) -> None:
        ax = self.ax_coll
        squeeze = np.geomspace(1e-4, 1.0, 200)
        cols = {1: _PAL["cool"], 2: _PAL["good"], 3: _PAL["warm"],
                4: _PAL["hot"]}
        for n, col in cols.items():
            ax.plot(squeeze, squeeze ** n, lw=2.2 if n == 3 else 1.6,
                    color=col,
                    label=f"S^{n}  cross-section:  (f0/F)^{n}")
        ax.axvline(1e-3, color=_PAL["rule"], lw=0.9, ls=":")
        ax.plot([1e-3], [1e-3], "o", ms=6, color=_PAL["cool"], zorder=6)
        ax.plot([1e-3], [1e-9], "o", ms=6, color=_PAL["warm"], zorder=6)
        ax.annotate("", xy=(1e-3, 1e-9), xytext=(1e-3, 1e-3),
                    arrowprops=dict(arrowstyle="<->", color=_PAL["bad"],
                                    lw=1.3))
        ax.text(1.6e-3, 2.4e-6,
                "a MILLION-fold\nunderstatement\nby the 2-D drawing",
                color=_PAL["bad"], fontsize=6.9, family="monospace",
                va="center", ha="left")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(1e-14, 2.0)
        _style(ax, "the collapse is f^n, not f",
               "squeeze   f0 / F", "physical measure of a FIXED angular patch")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.8, loc="lower right",
                  framealpha=0.93)
        ax.text(0.03, 0.06,
                "the ANGULAR footprint is unchanged along every curve;\n"
                "only the PHYSICAL measure carries the exponent.\n"
                "the yes/no overlap criterion was always angular.",
                transform=ax.transAxes, ha="left", va="bottom",
                color=_PAL["dim"], fontsize=6.8, family="monospace")

    # ── parity ──────────────────────────────────────────────────────────────
    def _parity(self) -> None:
        ax = self.ax_par
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_color(_PAL["rule"])
        ax.set_title("orientability flips with PARITY, and the repo uses two",
                     color=_PAL["text"], fontsize=9.6, pad=8,
                     family="monospace")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        dims = (2, 3, 4, 5)
        xs = np.linspace(0.17, 0.83, len(dims))
        for x, d in zip(xs, dims):
            here = (d == 3)
            for y, n, lab in ((0.63, d, "spatial   RP^%d" % d),
                              (0.34, d - 1, "exchange  RP^%d" % (d - 1))):
                ok = projective_space_is_orientable(n)
                col = _PAL["good"] if ok else _PAL["bad"]
                ax.add_patch(plt.Rectangle((x - 0.085, y - 0.075), 0.17, 0.15,
                                           facecolor=col, alpha=0.20,
                                           edgecolor=col,
                                           lw=2.0 if here else 1.0))
                ax.text(x, y + 0.028, lab, ha="center", color=_PAL["text"],
                        fontsize=6.9, family="monospace")
                ax.text(x, y - 0.036,
                        "orientable" if ok else "NON-orientable",
                        ha="center", color=col, fontsize=6.9,
                        family="monospace")
            ax.text(x, 0.84, "spatial d = %d%s" % (d, "  <- here" if here else ""),
                    ha="center", color=_PAL["text"] if here else _PAL["dim"],
                    fontsize=7.4 if here else 7.0, family="monospace")
            ax.text(x, 0.235, "pi_1 = %s" % ("Z (braid)" if d == 2 else "Z_2"),
                    ha="center", color=_PAL["dim"], fontsize=6.5,
                    family="monospace")
        ax.text(0.5, 0.035,
                "RP^n orientable iff n ODD:  det(-I_(n+1)) = (-1)^(n+1).  the "
                "two quotients are ONE APART,\nhence ALWAYS OPPOSITE -- "
                "raising the spatial dimension SWAPS which is non-orientable.\n"
                "at d = 3 the Pin- structure lives on the exchange RP^2.",
                ha="center", va="bottom", color=_PAL["dim"], fontsize=6.7,
                family="monospace")

    # ── assembly ────────────────────────────────────────────────────────────
    def draw(self) -> "HypersphereFigure":
        self._shell()
        self._angles()
        self._collapse()
        self._parity()
        self.fig.text(0.5, 0.962,
                      "v70  -  what higher dimension does to the bulk picture",
                      color=_PAL["text"], fontsize=15.5, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.925,
                      "3-D intuition models an extra dimension as 'the same "
                      "sphere with another direction available'.  It is not.",
                      color=_PAL["dim"], fontsize=8.2, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.894,
                      "THE ANGULAR OVERLAP CAN STAY FINITE WHILE THE PHYSICAL "
                      "OVERLAP COLLAPSES AS f0^n  -  a fixed angular footprint, "
                      "a millionfold smaller proper region at n = 3",
                      color=_PAL["good"], fontsize=7.6, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.866,
                      "and the peak everyone quotes is the unit BALL's (d = 5), "
                      "not the sphere's (d = 7, S^6) -- and it moves to d = 100 "
                      "at R = 4, so it is not a fact about dimension alone",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.840,
                      "MEASURE, ORIENTATION AND FRAMES ONLY -- no field "
                      "equation, nothing evolving; and S^3 is exceptional "
                      "(SU(2), parallelizable, Hopf), so none of this "
                      "extrapolates smoothly",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.018,
                      "the bottom-right panel says what would have to be "
                      "RE-DERIVED if the dimension moved  -  it is not an "
                      "argument that it should",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")
        return self


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--still", default="v70.png")
    ap.add_argument("--dpi", type=int, default=110)
    args = ap.parse_args(argv)
    matplotlib.use("Agg")
    fig = HypersphereFigure().draw()
    fig.fig.savefig(args.still, dpi=args.dpi, facecolor=_PAL["bg"])
    print(f"wrote {args.still}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
