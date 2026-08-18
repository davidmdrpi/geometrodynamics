#!/usr/bin/env python3
"""
Geometrodynamic QED — v60: the throat has an interior, and the interior is
the delay
==========================================================================

PRs #253-#259 all carried the same disclaimer: the throat is point-supported --
no interior, no proper length, no delay -- and what stood in for one was a
rank-one mouth-transfer model, lossy for kappa < 1. This round replaces it with
a tube whose Dirichlet-to-Neumann map is exact, and the boundary condition
becomes FREQUENCY-DEPENDENT. That dependence is the interior.

What the panels show
--------------------
**Top left - the delay.** The onset of the throat's cross-mouth response
against tube length. Slope 1 below the ambient path and flat above it: the tube
takes L to traverse, the ambient takes d, and `min(L, d)` decides which arrives
first. The point throat sits at zero, because that is what a point throat is.

**Top right - the impulse response itself.** The throat operator's own response
to a pulse at a mouth, source and observer removed from the problem. Same mouth
in and out starts instantly - a mouth reflects - and opposite mouths starts at
the traversal time.

**Bottom left - there is no point limit, there is a band.** Freezing A at one
frequency is exact there and 121% wrong at 3 lam_0. The two channels fail
differently: the antisymmetric one converges to -L/2A like L^2, the symmetric
one diverges like 2/(A lam L), because a massless tube holds a zero mode and a
point cannot.

**Bottom right - and that zero mode breaks PR #258's tomography.** det S goes
to zero linearly in lam, so the static response is rank one and the
disconnection defect diverges. Give the tube an interior mass and the rank
returns -- with W = -beta(lam) exact to 1e-13.

What is put in
--------------
Two lines: the tube's DtN map and the value/flux matching at the mouths. The
interior is one-dimensional and the mouths are still points; the growing mode
at sigma* = 2 sqrt(pi/A) is what that costs, and it is a single-mouth object.
**No backreaction.**

Usage
-----
    python scripts/geometrodynamics_v60_finite_throat.py --still out.png
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.waves.finite_throat import (
    FOUR_PI,
    FiniteThroat,
    causal_onset,
    impulse_response,
    measure_the_growing_mode_belongs_to_the_mouth,
    measure_the_point_throat_is_a_single_frequency_match,
    measure_the_static_limit_is_rank_one_and_the_defect_diverges,
    measure_the_throat_transmits_at_the_traversal_time,
)
from geometrodynamics.waves.two_wave import RetardedGrid

_PAL = {
    "bg": "#010106",
    "panel": "#02020a",
    "grid": "#0e1a2a",
    "rule": "#1a2838",
    "text": "#e8ecf4",
    "dim": "#6a8aad",
    "free": "#63798f",
    "sym": "#7cff9e",          # the tube
    "anti": "#ffb347",         # the ambient path
    "hot": "#ff6b8a",          # the point throat / what is lost
    "cool": "#5cc8ff",
}


class FiniteThroatFigure:
    """Delay, impulse response, the band, and the rank collapse."""

    def __init__(self, figsize=(13.8, 8.6)) -> None:
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, width_ratios=[1.0, 1.05], height_ratios=[1.0, 1.0],
            left=0.068, right=0.975, top=0.845, bottom=0.135,
            wspace=0.26, hspace=0.62)
        self.ax_del = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_imp = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_bnd = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_rnk = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])
        self._delay = measure_the_throat_transmits_at_the_traversal_time()
        self._band = measure_the_point_throat_is_a_single_frequency_match()
        self._rank = measure_the_static_limit_is_rank_one_and_the_defect_diverges()
        self._mode = measure_the_growing_mode_belongs_to_the_mouth()

    # -- the delay -----------------------------------------------------------
    def _draw_delay(self) -> None:
        ax = self.ax_del
        rows = self._delay["rows"]
        d = 1.3
        ls = np.array([r["length"] for r in rows])
        on = np.array([r["onset_opposite"] for r in rows])
        off = float(np.mean([min(r["length"], d) - r["onset_opposite"]
                             for r in rows if r["length"] < d]))
        grid = np.linspace(0.0, 3.4, 200)
        ax.plot(grid, np.minimum(grid, d) - off, lw=1.4, color=_PAL["dim"],
                ls=(0, (4, 3)), zorder=3, label="min(L, d)  −  probe width")
        ax.plot(ls, on, "o", ms=10.0, color=_PAL["sym"], mec=_PAL["bg"],
                mew=0.8, zorder=6, label="the tube's own arrival")
        ax.plot([0.9], [self._delay["point_throat_onset_opposite"]], "X",
                ms=13.0, color=_PAL["hot"], mec=_PAL["bg"], mew=0.8, zorder=7,
                label="the SAME throat, A frozen: instantaneous")
        ax.annotate("", xy=(0.9, self._delay["point_throat_onset_opposite"]),
                    xytext=(0.9, 0.74), arrowprops=dict(
                        arrowstyle="<->", color=_PAL["hot"], lw=0.9,
                        shrinkA=3, shrinkB=3))
        ax.axvline(d, color=_PAL["anti"], lw=1.0, ls=(0, (2, 4)), zorder=2)
        ax.annotate("the ambient path, d = 1.3\n(there whether or not\nthe "
                    "mouths are joined)", xy=(d + 0.08, 0.62),
                    color=_PAL["anti"], fontsize=6.5, family="monospace",
                    linespacing=1.7)
        ax.annotate(f"slope below d:  "
                    f"{self._delay['slope_below_the_ambient_path']:.4f}\n"
                    f"predicted:      1\n"
                    f"spread above d: "
                    f"{self._delay['onset_spread_above_it']:.1e}",
                    xy=(0.52, 0.30), xycoords="axes fraction", va="top",
                    color=_PAL["text"], fontsize=6.8, family="monospace",
                    linespacing=1.8)
        ax.set_xlabel("tube proper length  L", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("onset of the cross-mouth response",
                      color=_PAL["dim"], fontsize=8)
        ax.set_xlim(-0.15, 3.35)
        ax.set_ylim(-0.08, 1.45)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper left", fontsize=6.2, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title("the throat transmits at the traversal time",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # -- the impulse response ------------------------------------------------
    def _draw_impulse(self) -> None:
        ax = self.ax_imp
        t = FiniteThroat(separation=1.3, length=0.6, area=FOUR_PI)
        sig = float(t.growing_modes()["symmetric"])
        grid = RetardedGrid(n=1 << 17, span=300.0, eps=sig + 0.8)
        imp = impulse_response(t, grid, width=0.03)
        ts = imp["times"]
        sl = ts < 3.4
        for key, col, lab in (("same_mouth", _PAL["anti"],
                               "same mouth in and out"),
                              ("opposite_mouths", _PAL["sym"],
                               "opposite mouths")):
            y = imp[key] / np.abs(imp[key][sl]).max()
            ax.plot(ts[sl], y[sl], lw=1.5, color=col, zorder=5, label=lab)
        frozen = t.boundary_matrix(complex(1.0)).real.astype(complex)
        pt = impulse_response(t, grid, width=0.03, constant=frozen)
        y = pt["opposite_mouths"] / np.abs(pt["opposite_mouths"][sl]).max()
        ax.plot(ts[sl], y[sl], lw=1.1, color=_PAL["hot"], alpha=0.85, zorder=4,
                label="point throat, opposite mouths")
        for x, c, lab in ((0.0, _PAL["anti"], "0"),
                          (t.length, _PAL["sym"], "L"),
                          (2.0 * t.length, _PAL["anti"], "2L"),
                          (t.separation, _PAL["cool"], "d")):
            ax.axvline(x, color=c, lw=0.8, ls=(0, (2, 4)), zorder=2)
            ax.annotate(lab, xy=(x + 0.03, 1.02), color=c, fontsize=6.6,
                        family="monospace")
        ax.annotate("the operator's own response — no source,\n"
                    "no observer, so this is the throat's\n"
                    "causal structure and nothing else",
                    xy=(0.28, 0.44), xycoords="axes fraction", va="top",
                    color=_PAL["dim"], fontsize=6.2, family="monospace",
                    linespacing=1.75)
        ax.set_xlabel("t", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("r_ij(t),  normalized", color=_PAL["dim"], fontsize=8)
        ax.set_ylim(-1.15, 1.15)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="lower right", fontsize=6.4, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title("a mouth reflects instantly; the tube takes L",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # -- the band ------------------------------------------------------------
    def _draw_band(self) -> None:
        ax = self.ax_bnd
        rows = self._band["band"]
        xs = np.array([r["lambda"] for r in rows])
        err = np.array([max(r["relative_error"], 1e-17) for r in rows])
        ax.semilogy(xs, err, "o-", ms=5.0, lw=1.5, color=_PAL["hot"],
                    zorder=5, label="error of a FROZEN A(λ₀)")
        ax.axhline(1.0, color=_PAL["dim"], lw=0.9, ls=(0, (3, 3)), zorder=2)
        ax.annotate("100%", xy=(2.55, 1.25), color=_PAL["dim"], fontsize=6.4,
                    family="monospace")
        short = self._band["short_tubes"]
        ls = np.array([r["length"] for r in short])
        ax2 = ax.twiny()
        ax2.semilogy(ls, [abs(r["anti_error"]) for r in short], "s--", ms=4.5,
                     lw=1.2, color=_PAL["sym"], zorder=5,
                     label="antisymmetric: |A − (−L/2A)|,  O(L²)")
        ax2.semilogy(ls, [abs(r["sym"]) for r in short], "^--", ms=4.5, lw=1.2,
                     color=_PAL["anti"], zorder=5,
                     label="symmetric: A_sym,  diverges as 2/(AλL)")
        ax2.set_xlabel("tube length L   (green / amber)", color=_PAL["dim"],
                       fontsize=7)
        ax2.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax2.spines.values():
            sp.set_color(_PAL["rule"])
        ax.annotate("a point throat is a finite throat read at ONE λ.\n"
                    "the antisymmetric channel has a limit; the\n"
                    "symmetric one cannot — a massless tube holds a\n"
                    "zero mode and a point has nowhere to put it.",
                    xy=(0.30, 0.70), xycoords="axes fraction", va="top",
                    color=_PAL["dim"], fontsize=6.3, family="monospace",
                    linespacing=1.75)
        ax.set_xlabel("λ / λ₀   (pink)", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("relative error / channel value", color=_PAL["dim"],
                      fontsize=8)
        ax.set_ylim(1e-6, 3e1)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        leg = ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=6.2,
                        framealpha=0.0, labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], which="both")
        ax.set_title("there is no point limit — there is a band",
                     color=_PAL["text"], fontsize=8.6, pad=30)

    # -- the rank collapse ---------------------------------------------------
    def _draw_rank(self) -> None:
        ax = self.ax_rnk
        rows = self._rank["rows"]
        lam = np.array([r["lambda"] for r in rows])
        det = np.array([abs(r["det_S"]) for r in rows])
        ax.loglog(lam, det, "o-", ms=5.5, lw=1.5, color=_PAL["sym"], zorder=5,
                  label="|det S|,  massless tube")
        ax.loglog(lam, self._rank["linear_coefficient"] * lam, lw=1.0,
                  color=_PAL["dim"], ls=(0, (3, 3)), zorder=3,
                  label=f"{self._rank['linear_coefficient']:.2f} · λ")
        massive = self._rank["massive"]
        ms = np.array([r["mass"] for r in massive])
        ax.loglog(ms ** 2, [abs(r["det_S"]) for r in massive], "s", ms=6.0,
                  color=_PAL["anti"], mec=_PAL["bg"], mew=0.8, zorder=6,
                  label="|det S| at λ = 0,  interior mass m  (x = m²)")
        ax.annotate("rank 1 — the symmetric channel is EMPTY,\n"
                    "so PR #258's defect W = S₁₂/det S − G₀\n"
                    f"diverges like 1/λ:  {rows[0]['defect']:.2e} at λ = 1e-08\n"
                    "\n"
                    "a POINT throat is statically rank 2.\n"
                    "that is a falsifiable difference, from a\n"
                    "STATIC measurement.",
                    xy=(0.035, 0.97), xycoords="axes fraction", va="top",
                    color=_PAL["text"], fontsize=6.4, family="monospace",
                    linespacing=1.8)
        ax.annotate("and off the collapse, W = −β(λ) EXACTLY\n"
                    f"(worst error {self._rank['worst_defect_error']:.1e}) — "
                    "PR #258's theorem\n"
                    "survives, returning the interior's amplitude.",
                    xy=(0.40, 0.46), xycoords="axes fraction", va="top",
                    color=_PAL["dim"], fontsize=6.4, family="monospace",
                    linespacing=1.8)
        ax.set_xlabel("λ   (or m² for the squares)", color=_PAL["dim"],
                      fontsize=8)
        ax.set_ylabel("|det S|", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="lower right", fontsize=6.3, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], which="both")
        ax.set_title("the static response collapses to rank one",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # -- frame ---------------------------------------------------------------
    def draw(self) -> None:
        self._draw_delay()
        self._draw_impulse()
        self._draw_band()
        self._draw_rank()
        self.fig.suptitle("v60 — THE THROAT HAS AN INTERIOR, AND THE INTERIOR "
                          "IS THE DELAY",
                          color=_PAL["text"], fontsize=13.2, y=0.962,
                          family="monospace")
        self.fig.text(0.5, 0.908,
                      "two lines put in — the tube's Dirichlet-to-Neumann map "
                      "and value/flux matching at the mouths   ·   everything "
                      "else follows",
                      color=_PAL["dim"], fontsize=8.4, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.070,
                      "the boundary condition is now FREQUENCY-DEPENDENT, and "
                      "that dependence is the interior:  A(λ) = −N(λ)⁻¹,  "
                      "self-adjoint at every λ (defect 0.0) against 0.30 for "
                      "the rank-one transfer model it replaces",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        mode = self._mode
        self.fig.text(0.5, 0.047,
                      f"the price of point mouths: one growing mode per "
                      f"channel at σ* = 2√(π/A), matched to "
                      f"{mode['worst_closed_form_error']:.1e} with NO L in it "
                      f"and blind to the mouth separation to "
                      f"{mode['separation_spread_far']:.0e} — a SINGLE-MOUTH "
                      f"object, at the scale √(A/4π)",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.022,
                      "conformally coupled scalar on a fixed ESU   ·   the "
                      "interior is one-dimensional and L, A, m are chosen, not "
                      "derived   ·   NO backreaction",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")


def still(path: str) -> None:
    fig = FiniteThroatFigure()
    fig.draw()
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v60.png")
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
