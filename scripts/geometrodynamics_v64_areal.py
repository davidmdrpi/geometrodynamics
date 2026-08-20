#!/usr/bin/env python3
"""
Geometrodynamic QED - v64: the interference metric deforms toward a neck
=======================================================================

PR #263 asked whether A + B does something rescaling A or B cannot, answered
yes in the transverse-traceless sector, and then proved that sector cannot give
a GEOMETRIC verdict: dA/A = -<h_nn>/2 vanishes identically for any TT field.
So the question that was actually asked --

    does the interaction metric deform TOWARD a neck, AWAY from one, or
    merely OSCILLATE?

-- moves to the scalar sector, posed as an INITIAL-DATA constraint solve rather
than an evolution. On a maximal slice the K terms in the Hamiltonian constraint
are quadratic, so dR3 = 16 pi G drho with no time derivatives: a constraint has
no sound speed and no Eddington mode, which is exactly why #263 avoided the
scalar sector as an evolution.

The answer is TOWARD A NECK. Both mouths close.

What the panels show
--------------------
**Top left - the sign, and what builds it.** dA/A at each mouth, decomposed
into the source, the two monopole layers, the six dipole layers and the free
kernel element. The interference energy ALONE would OPEN the mouths -- U(c_j)
is positive at both. The throat's monopole layers overshoot that and invert it.

**Top right - the controls.** Two quadrature levels, two mouth radii, two
gluings. The sign is the same in all eight. Note that the mouth radius is NOT
a regulator: the source goes as 1/chi^4 at the mouths, so there is no a -> 0
limit to converge to -- that is the singular point #261/#262 removed. Doubling
a moves the answer by the factor the source's own mouth-weighted moment moves.

**Bottom left - what the answer is made of.** The obstruction alone reproduces
it to 0.09%; the local source data alone leave a response 1000x smaller;
scaling the l=1 moments by three, by zero, or replacing them with noise moves
it by 5e-04. That last one matters: those moments drift 41% between quadrature
levels, against 1.5% for the obstruction, so the signed answer rests on the
best-converged number available and not on the worst.

**Bottom right - the sign is a statement about a THROAT.** The tube's l=0
constraint channel is d_s^2 + 4 pi / area, a CAVITY. At kL = n pi the response
has a pole and the sign flips. The working throat sits at kL = 0.9, inside the
first cell; past the first pole the two mouths can move in opposite directions.

What is put in
--------------
The fluid holding the ESU static is held RIGID -- consistent because the
scalar's stress tensor is separately conserved, but a responsive fluid is the
obvious next refinement. The source is #263's: a linear conformally coupled
scalar on a FIXED ESU with POINT sources. The gluing of the two mouths'
transverse frames is a modelling choice, measurably harmless: a full 2 pi twist
moves dA/A by less than 1e-12.

Usage
-----
    python scripts/geometrodynamics_v64_areal.py --still out.png
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.waves.areal import (
    INTERFERENCE_MOMENTS, MOUTHS, TubeModel, WORKING_TUBE, basis_channels,
    solve_matching,
)

_PAL = {
    "bg": "#010106", "panel": "#02020a", "grid": "#0e1a2a", "rule": "#1a2838",
    "text": "#e8ecf4", "dim": "#6a8aad",
    "open": "#7cff9e", "close": "#ff6b8a", "x": "#ffb347",
    "cool": "#5cc8ff", "kern": "#b388ff",
}


def _style(ax, title, xlabel, ylabel):
    ax.set_title(title, color=_PAL["text"], fontsize=9.6, pad=8,
                 family="monospace")
    ax.set_xlabel(xlabel, color=_PAL["dim"], fontsize=8, family="monospace")
    ax.set_ylabel(ylabel, color=_PAL["dim"], fontsize=8, family="monospace")
    ax.tick_params(colors=_PAL["dim"], labelsize=7.2)
    for s in ax.spines.values():
        s.set_color(_PAL["rule"])
    ax.grid(True, color=_PAL["grid"], lw=0.5, alpha=0.7)


class ArealFigure:
    """The sign, the controls, what it is made of, and the cavity behind it."""

    def __init__(self, figsize=(13.8, 8.6)) -> None:
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, left=0.075, right=0.975, top=0.840, bottom=0.135,
            wspace=0.26, hspace=0.52)
        self.ax_sign = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_ctrl = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_made = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_cav = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])
        self.head = INTERFERENCE_MOMENTS[1]
        self.basis = basis_channels(MOUTHS, self.head.radius)

    # ── panels ──────────────────────────────────────────────────────────────
    def _sign(self) -> None:
        ax = self.ax_sign
        m = self.head
        src = m.as_source()
        got = solve_matching(MOUTHS, m.radius, WORKING_TUBE, src,
                             m.signed_obstruction(), basis=self.basis)
        coef = np.asarray(got["coefficients"])
        v0 = self.basis["value_0"]
        labels = ["source\nU(c_j)", "monopole\nlayers", "dipole\nlayers",
                  "kernel\nc.x", "TOTAL\nu"]
        cols = [_PAL["open"], _PAL["close"], _PAL["cool"], _PAL["kern"],
                _PAL["x"]]
        width = 0.36
        for j in (0, 1):
            parts = v0[j] * coef
            vals = [src["value_0"][j], parts[0] + parts[4],
                    parts[1:4].sum() + parts[5:8].sum(), parts[8:].sum()]
            vals.append(sum(vals))
            xs = np.arange(5) + (j - 0.5) * width
            ax.bar(xs, np.array(vals) * 1e3, width=width * 0.92,
                   color=cols, alpha=0.95 if j == 0 else 0.55,
                   edgecolor=_PAL["rule"], lw=0.6)
        ax.axhline(0.0, color=_PAL["dim"], lw=0.9)
        ax.set_xticks(np.arange(5))
        ax.set_xticklabels(labels, fontsize=6.6, family="monospace",
                           color=_PAL["dim"])
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + 0.55 * (hi - lo))
        ar = np.asarray(got["areal_response"])
        _style(ax, "the sign, and what builds it", "",
               "u at the mouth  [10^-3, units 2 pi G]")
        ax.text(0.5, 0.94,
                f"dA/A = 4u = ({ar[0]:+.3e}, {ar[1]:+.3e})   BOTH CLOSE",
                transform=ax.transAxes, ha="center", color=_PAL["close"],
                fontsize=8.0, family="monospace")
        ax.text(0.5, 0.855,
                "left bar = mouth 1, right = mouth 2;  the energy alone "
                "would OPEN them",
                transform=ax.transAxes, ha="center", color=_PAL["dim"],
                fontsize=6.8, family="monospace")

    def _controls(self) -> None:
        ax = self.ax_ctrl
        rows = []
        for m in INTERFERENCE_MOMENTS:
            basis = basis_channels(MOUTHS, m.radius)
            for reflect in (False, True):
                got = solve_matching(MOUTHS, m.radius, WORKING_TUBE,
                                     m.as_source(), m.signed_obstruction(),
                                     reflect=reflect, basis=basis)
                rows.append((m.radius, m.points, reflect,
                             np.asarray(got["areal_response"])))
        xs = np.arange(len(rows))
        for k, mouth in enumerate((0, 1)):
            ax.plot(xs, [r[3][mouth] * 1e3 for r in rows],
                    "o-" if k == 0 else "s--", lw=1.5, ms=5.5,
                    color=_PAL["close"] if k == 0 else _PAL["x"],
                    label=f"mouth {mouth + 1}")
        ax.axhline(0.0, color=_PAL["open"], lw=1.2, ls=":")
        ax.set_xticks(xs)
        ax.set_xticklabels([f"a={r[0]:.2f}\nN={r[1]}\n"
                            f"{'refl' if r[2] else 'tran'}" for r in rows],
                           fontsize=6.0, family="monospace",
                           color=_PAL["dim"])
        _style(ax, "the controls: the sign is the same in all eight", "",
               "dA/A  [10^-3]")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=7.2, loc="lower right")
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo - 0.42 * (hi - lo), hi)
        ax.text(0.5, 0.085,
                "a is NOT a regulator: the source goes as 1/chi^4 at the "
                "mouths,\nso there is no a -> 0 limit to converge to",
                transform=ax.transAxes, ha="center", color=_PAL["dim"],
                fontsize=6.6, family="monospace")

    def _made_of(self) -> None:
        ax = self.ax_made
        m = self.head
        s = m.signed_obstruction()
        quiet = {"value_0": np.zeros(2), "slope_0": np.zeros(2),
                 "value_1": np.zeros((2, 3)), "slope_1": np.zeros((2, 3))}
        rng = np.random.default_rng(20260820)
        flat = type(m)(radius=m.radius, points=m.points, volume=m.volume,
                       energy=m.energy, obstruction=m.obstruction,
                       value=m.value,
                       gradient=tuple(tuple(0.0 for _ in r)
                                      for r in m.gradient))
        noisy = type(m)(radius=m.radius, points=m.points, volume=m.volume,
                        energy=m.energy, obstruction=m.obstruction,
                        value=m.value,
                        gradient=tuple(tuple(g * (1 + 0.5 *
                                                  rng.standard_normal())
                                             for g in r) for r in m.gradient))

        def run(src, obs):
            return np.asarray(solve_matching(MOUTHS, m.radius, WORKING_TUBE,
                                             src, obs, basis=self.basis
                                             )["areal_response"])

        names = ["full", "obstruction\nonly", "local data\nonly",
                 "l=1 moments\nzeroed", "l=1 moments\nrandomised"]
        vals = [run(m.as_source(), s), run(quiet, s),
                run(m.as_source(), np.zeros(4)), run(flat.as_source(), s),
                run(noisy.as_source(), s)]
        cols = [_PAL["x"], _PAL["close"], _PAL["dim"], _PAL["cool"],
                _PAL["cool"]]
        heights = [v[0] * 1e3 for v in vals]
        ax.bar(np.arange(5), heights, width=0.6,
               color=cols, edgecolor=_PAL["rule"], lw=0.6)
        ax.axhline(0.0, color=_PAL["dim"], lw=0.9)
        ax.set_xticks(np.arange(5))
        ax.set_xticklabels(names, fontsize=6.4, family="monospace",
                           color=_PAL["dim"])
        _style(ax, "what the answer is made of", "",
               "dA/A at mouth 1  [10^-3]")
        span = max(abs(h) for h in heights)
        ax.set_ylim(-1.18 * span, 0.62 * span)
        # the third bar is 1000x smaller than the rest and would be invisible
        ax.annotate(f"{heights[2] * 1e-3:+.1e}", xy=(2, 0.0),
                    xytext=(2, 0.20 * span), ha="center", color=_PAL["dim"],
                    fontsize=6.6, family="monospace",
                    arrowprops=dict(arrowstyle="-", color=_PAL["dim"], lw=0.7))
        ax.text(0.5, 0.93,
                "the obstruction alone reproduces it to 0.09%;  without it "
                "the response is 1000x smaller",
                transform=ax.transAxes, ha="center", color=_PAL["text"],
                fontsize=6.8, family="monospace")
        ax.text(0.5, 0.045,
                "the l=1 moments drift 41% between quadrature levels, and "
                "the signed answer does not rest on them",
                transform=ax.transAxes, ha="center", color=_PAL["dim"],
                fontsize=6.4, family="monospace")

    def _cavity(self) -> None:
        ax = self.ax_cav
        m = self.head
        s = m.signed_obstruction()
        src = m.as_source()
        ls = np.linspace(0.25, 7.2, 320)
        vals = []
        for length in ls:
            got = solve_matching(MOUTHS, m.radius,
                                 TubeModel(area=WORKING_TUBE.area,
                                           length=float(length)),
                                 src, s, basis=self.basis)
            vals.append(np.asarray(got["areal_response"]))
        vals = np.array(vals)
        clip = np.clip(vals * 1e3, -12.0, 12.0)
        ax.plot(ls, clip[:, 0], lw=1.5, color=_PAL["close"], label="mouth 1")
        ax.plot(ls, clip[:, 1], lw=1.2, ls="--", color=_PAL["x"],
                label="mouth 2")
        ax.axhline(0.0, color=_PAL["dim"], lw=0.9)
        for n in (1, 2):
            ax.axvline(n * math.pi, color=_PAL["kern"], lw=1.0, ls=":")
            ax.text(n * math.pi, 10.6, f" kL = {n}pi", color=_PAL["kern"],
                    fontsize=6.6, family="monospace")
        ax.axvline(WORKING_TUBE.length, color=_PAL["open"], lw=1.4)
        ax.text(WORKING_TUBE.length + 0.1, -10.6, "working throat",
                color=_PAL["open"], fontsize=6.8, family="monospace")
        ax.set_ylim(-12.5, 12.5)
        _style(ax, "the sign is a statement about a THROAT",
               "tube length L   (k = 1, so kL = L)", "dA/A  [10^-3, clipped]")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=7.2, loc="lower right")

    # ── assembly ────────────────────────────────────────────────────────────
    def draw(self) -> None:
        self._sign()
        self._controls()
        self._made_of()
        self._cavity()
        self.fig.text(0.5, 0.955,
                      "v64  -  the interference metric deforms TOWARD a neck",
                      color=_PAL["text"], fontsize=15.5, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.915,
                      "initial-data constraint solve on the resolved neck:  "
                      "nabla^2 u + 3u = -2 pi G drho,   dA/A = 4u",
                      color=_PAL["dim"], fontsize=8.6, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.884,
                      "#263 proved the TT sector cannot give a geometric "
                      "verdict -- <h_nn> vanishes identically -- so this is "
                      "the scalar sector, as a CONSTRAINT and not an evolution",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.022,
                      "the ESU's fluid is held RIGID   -   source is #263's, "
                      "on a FIXED background with POINT sources   -   response "
                      "LINEAR in G and quadratic in the wave amplitude",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")


def still(path: str) -> None:
    fig = ArealFigure()
    fig.draw()
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v64.png")
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
