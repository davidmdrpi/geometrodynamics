#!/usr/bin/env python3
"""
Geometrodynamic QED - v69: the centre as a finite bearing, not a point
======================================================================

Every picture in this arc has put a POINT in the middle. A point is where the
clock-hands story works -- two arms meet and change direction for free, because
at a point there is no angular direction left to change -- and it is bought
with f = 0, where the geometry stops existing.

Blow the point up. Keep dl^2 = ds^2 + f(s)^2 dOmega^2 and set f_min = f0 > 0.
The middle becomes a small circle (in the 2-D cross-section) or the space of
radial directions S^(d-1) / RP^(d-1) in the bulk. Nothing is singular, the two
arms are free to differ, and the hinge costs something finite.

What the panels show
--------------------
**Top left - the profile, with two unequal arms.** f(s) against proper distance,
outer arm to the left and inner to the right, pinched to f0 in the middle. The
two arms are DIFFERENT LENGTHS and different widths, which the concentric
inner/outer picture could not represent. L(F) is closed-form and reduces to
physical_throat's own length() when the two are made equal.

**Top right - scale transport.** The same route drawn as a strip whose width is
the PHYSICAL width w(s) = f(s) dtheta of a feature of fixed ANGULAR width. The
angular width is what is carried and is constant; the physical width is
squeezed into the bearing by f_o/f0 and let out again at f_i/f0. The end-to-end
ratio is f_i/f_o and does not involve f0.

**Bottom left - the bearing, and what intersection becomes.** The central circle
with two fronts landing on it at their angular positions. They meet iff their
angular extents overlap -- a question about ANGLES, with f0 nowhere in it. What
f0 sets is how BIG the meeting is: f0 x (overlap angle). Two configurations are
drawn, one that meets and one that misses.

**Bottom right - the hinge cost, and the correction.** T(alpha) against alpha
on log-log, so the claim is a SLOPE: the arc f0*alpha has slope 1, the geodesic
has slope 2. Four curves -- the integrated geodesic, the leading law
alpha^2/(2 I2), the same law with its I4 correction, and the linear guess the
corner route would charge. The geodesic is QUADRATIC and far below the arc:
1.25% of it at alpha = 0.1 and 36% at pi.

THE MOMENT HIERARCHY. At angular dimension q = 2 -- the physical case, an S^2
cross-section -- I2 = int ds/f^2 is ALSO physical_throat's resistance, so the
geometric hinge and the monopole channel are one integral there. They are one
DIRICHLET FORM int w phi'^2 ds read twice, but the WEIGHTS differ in general:
the azimuth's is the metric coefficient f^2 in any dimension, the monopole's is
the volume element f^q. They match at q = 2 and nowhere else.
What is universal is the LEADING FUNCTIONAL FORM alpha^2/(2 I2); I2 itself is
not (4/f0 here, pi/f0 on a hyperbolic neck). I4 = int ds/f^4 is the first
ADDITIONAL INDEPENDENT moment, and where the neck's SHAPE is first felt:
shape = 1 - alpha^2 I4/(4 I2^3), i.e. 1 - alpha^2/120 here against
1 - alpha^2/(8 pi^2) hyperbolic. I6 and beyond enter at O(alpha^6).

What is put in
--------------
GEOMETRY only: a metric, its geodesics, and an angular width transported along
them. No field equation is solved. The scalar-flat profile f'^2 = 1 - f0/f is a
worked example chosen because every closed form in it is checkable against
physical_throat -- not a claim that the bearing must be that profile.

Usage
-----
    python scripts/geometrodynamics_v69_regularized_center.py --still out.png
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Circle, Wedge

from geometrodynamics.viz.regularized_center import (RegularizedCenter,
                                                     hyperbolic_neck)

_PAL = {
    "bg": "#010106", "panel": "#02020a", "grid": "#0e1a2a", "rule": "#1a2838",
    "text": "#e8ecf4", "dim": "#6a8aad",
    "outer": "#5cc8ff", "inner": "#ffb347", "neck": "#7cff9e",
    "warn": "#ff6b8a", "kern": "#b388ff", "hot": "#ff5ec7",
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


class BearingFigure:
    """Two unequal arms, a finite hinge, and what it costs to turn."""

    def __init__(self, centre: Optional[RegularizedCenter] = None,
                 figsize=(13.8, 8.9)) -> None:
        self.c = centre or RegularizedCenter(neck=0.12, outer=1.0, inner=0.42)
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, left=0.072, right=0.975, top=0.800, bottom=0.075,
            wspace=0.24, hspace=0.42)
        self.ax_prof = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_strip = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_ring = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_cost = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])

    # ── the profile ─────────────────────────────────────────────────────────
    def _profile(self) -> None:
        ax, c = self.ax_prof, self.c
        s, f = c.profile(n=600)
        ax.plot(s, f, lw=2.1, color=_PAL["neck"])
        ax.plot(s, -f, lw=2.1, color=_PAL["neck"])
        ax.fill_between(s, -f, f, color=_PAL["neck"], alpha=0.09)
        ax.axhline(0.0, color=_PAL["rule"], lw=0.7, ls=":")
        ax.plot([0.0], [c.neck], "o", ms=6, color=_PAL["hot"], zorder=6)
        ax.annotate(f"f0 = {c.neck:.3g}", xy=(0.0, c.neck),
                    xytext=(0.10 * c.outer_length(), 0.42 * c.outer),
                    color=_PAL["hot"], fontsize=7.4, family="monospace",
                    arrowprops=dict(arrowstyle="-", color=_PAL["hot"], lw=0.8))
        for x, lab, col in ((-c.outer_length(), "outer", _PAL["outer"]),
                            (c.inner_length(), "inner", _PAL["inner"])):
            ax.axvline(x, color=col, lw=1.3, ls="--")
        ax.text(-c.outer_length() * 0.98, -c.outer * 0.86,
                f"L_o = {c.outer_length():.4f}\nf_o = {c.outer:.3g}",
                color=_PAL["outer"], fontsize=7.2, family="monospace",
                ha="left", va="bottom")
        ax.text(c.inner_length() * 0.96, -c.outer * 0.86,
                f"L_i = {c.inner_length():.4f}\nf_i = {c.inner:.3g}",
                color=_PAL["inner"], fontsize=7.2, family="monospace",
                ha="right", va="bottom")
        ax.text(0.5, 0.95,
                f"L_o / L_i = {c.outer_length()/c.inner_length():.3f}"
                "   -- the arms need NOT match",
                transform=ax.transAxes, ha="center", color=_PAL["text"],
                fontsize=7.6, family="monospace")
        ax.set_ylim(-c.outer * 1.22, c.outer * 1.22)
        _style(ax, "the profile: a finite neck, two unequal arms",
               "proper distance from the bearing", "transverse scale  f(s)")

    # ── scale transport ─────────────────────────────────────────────────────
    def _strip(self) -> None:
        ax, c = self.ax_strip, self.c
        dtheta = 0.55
        s, f = c.profile(n=600)
        w = f * dtheta
        ax.fill_between(s, -0.5 * w, 0.5 * w, color=_PAL["kern"], alpha=0.30,
                        lw=0.0)
        ax.plot(s, 0.5 * w, lw=1.7, color=_PAL["kern"])
        ax.plot(s, -0.5 * w, lw=1.7, color=_PAL["kern"])
        for x, val, col, lab in (
                (-c.outer_length(), c.outer * dtheta, _PAL["outer"], "w_o"),
                (0.0, c.neck * dtheta, _PAL["hot"], "w_min"),
                (c.inner_length(), c.inner * dtheta, _PAL["inner"], "w_i")):
            ax.plot([x, x], [-0.5 * val, 0.5 * val], lw=3.0, color=col,
                    solid_capstyle="butt", zorder=5)
            ax.text(x, 0.5 * val + 0.035 * c.outer * dtheta,
                    f"{lab} = {val:.3g}", color=col, fontsize=7.0,
                    family="monospace", ha="center")
        ax.text(0.5, 0.035,
                f"fixed ANGULAR width dtheta = {dtheta:.2f}\n"
                f"squeezed {c.outer/c.neck:.1f}x into the bearing, then let out\n"
                f"w_i/w_o = {c.inner/c.outer:.4f} = f_i/f_o   (no f0 in it)",
                transform=ax.transAxes, ha="center", va="bottom",
                color=_PAL["dim"], fontsize=6.8, family="monospace")
        ax.set_ylim(-1.28 * c.outer * dtheta, 0.78 * c.outer * dtheta)
        _style(ax, "scale transport: the width follows f(s)",
               "proper distance from the bearing", "physical width  w = f dtheta")

    # ── the bearing ─────────────────────────────────────────────────────────
    def _ring(self) -> None:
        ax, c = self.ax_ring, self.c
        ax.set_aspect("equal")
        ax.set_xlim(-2.05, 2.05)
        ax.set_ylim(-1.62, 1.62)
        for spine in ax.spines.values():
            spine.set_color(_PAL["rule"])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("the bearing: intersection becomes overlap",
                     color=_PAL["text"], fontsize=9.6, pad=8,
                     family="monospace")

        for cx, (sep, wa, wb), tag in ((-1.00, (0.30, 0.55, 0.50), "they MEET"),
                                       (1.00, (1.35, 0.40, 0.35), "they MISS")):
            r = 0.70
            ax.add_patch(Circle((cx, 0.22), r, facecolor="none",
                                edgecolor=_PAL["dim"], lw=1.1, ls=":"))
            a_c, b_c = math.degrees(-0.5 * sep), math.degrees(0.5 * sep)
            for centre_deg, width, col in ((a_c, wa, _PAL["outer"]),
                                           (b_c, wb, _PAL["inner"])):
                half = math.degrees(0.5 * width)
                ax.add_patch(Wedge((cx, 0.22), r, centre_deg - half,
                                   centre_deg + half, width=0.17,
                                   facecolor=col, alpha=0.75, lw=0.0))
            reach = 0.5 * (wa + wb)
            overlap = max(reach - sep, 0.0)
            if overlap > 0.0:
                ax.add_patch(Wedge((cx, 0.22), r + 0.10,
                                   math.degrees(-0.5 * overlap),
                                   math.degrees(0.5 * overlap), width=0.09,
                                   facecolor=_PAL["hot"], lw=0.0))
            ax.plot([cx], [0.22], "o", ms=3.0, color=_PAL["rule"])
            ax.text(cx, -0.68, tag, color=_PAL["text"], fontsize=8.2,
                    family="monospace", ha="center")
            detail = (f"separation {sep:.2f} rad\n"
                      f"reach (w_a+w_b)/2 = {reach:.3f}\n")
            if overlap > 0.0:
                detail += (f"overlap {overlap:.3f} rad\n"
                           f"= {c.neck*overlap:.2e} across")
            else:
                detail += (f"gap {sep-reach:.3f} rad\n"
                           f"= {c.neck*(sep-reach):.2e} across")
            ax.text(cx, -0.84, detail, color=_PAL["dim"], fontsize=6.6,
                    family="monospace", ha="center", va="top")
        ax.text(0.0, 1.52,
                "WHETHER they meet is a question about ANGLES -- f0 is not in "
                "it.\nHOW BIG the meeting is, is f0 x (overlap angle).  "
                "f0 -> 0 shrinks\nthe overlap AND the gap together: the "
                "verdict survives, the length does not.",
                color=_PAL["dim"], fontsize=6.8, family="monospace",
                ha="center", va="top")

    # ── the hinge cost ──────────────────────────────────────────────────────
    def _cost(self) -> None:
        ax, c = self.ax_cost, self.c
        alphas = np.geomspace(0.02, math.pi, 42)
        geo = np.array([c.turn_cost(float(a)) for a in alphas])
        arc = c.neck * alphas
        # log-log, so the whole claim is a SLOPE: the arc is 1, the hinge is 2
        ax.plot(alphas, arc, lw=1.9, ls="--", color=_PAL["warn"],
                label="the arc  f0 alpha   (corner route, slope 1)")
        ax.plot(alphas, geo, lw=2.1, color=_PAL["neck"],
                label="the GEODESIC turn cost   (slope 2)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        _style(ax, "the hinge cost: quadratic, not linear",
               "turn angle  alpha  [rad]", "T(alpha)")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.9, loc="upper left",
                  framealpha=0.93)
        f_small = c.turn_cost(0.1) / (c.neck * 0.1)
        f_pi = c.turn_cost(math.pi) / (c.neck * math.pi)
        ax.text(0.028, 0.545,
                f"I2 = {c.resistance():.4e} = physical_throat's resistance\n"
                f"the geodesic spends {f_small*100:.2f}% of the arc at "
                f"alpha = 0.1,\n{f_pi*100:.0f}% at pi;  T(pi)/(L_o+L_i) = "
                f"{c.turn_cost(math.pi)/c.arm_length_sum():.2e}\n"
                f"and pi is the LARGEST separation there is",
                transform=ax.transAxes, ha="left", va="top",
                color=_PAL["dim"], fontsize=6.6, family="monospace")

        # The I4 correction is ~1% and invisible against the leading term on
        # a log axis, so it gets its own linear SHAPE axis, where it is the
        # whole signal.  Both profiles are evaluated at the SAME small neck --
        # 1/120 and 1/(8 pi^2) are LONG-ARM limits, and at the drawing's own
        # fat f0 the scalar-flat coefficient is 1/85, which would make the two
        # look alike for the wrong reason.
        ins = ax.inset_axes([0.455, 0.075, 0.520, 0.375])
        ins.set_facecolor("#04040e")
        thin = RegularizedCenter(neck=1e-4, outer=1.0, inner=1.0)
        lin = np.linspace(0.02, math.pi, 22)
        shape = np.array([thin.turn_cost(float(a))
                          / thin.turn_cost_small_angle(float(a)) for a in lin])
        ins.plot(lin, shape, lw=1.9, color=_PAL["neck"],
                 label="scalar-flat:  1 - a^2/120")
        ins.plot(lin, 1.0 - lin ** 2 / 120.0, lw=0.9, ls=":",
                 color=_PAL["kern"])
        coarse = np.linspace(0.02, math.pi, 9)
        hyp = np.array([hyperbolic_neck(1e-4, 1.0, 1.0, float(a))["shape"]
                        for a in coarse])
        ins.plot(coarse, hyp, lw=1.6, ls="--", color=_PAL["inner"],
                 label="hyperbolic:   1 - a^2/79")
        ins.plot(lin, 1.0 - lin ** 2 / (8.0 * math.pi ** 2), lw=0.9, ls=":",
                 color=_PAL["kern"])
        ins.axhline(1.0, color=_PAL["rule"], lw=0.8, ls=":")
        ins.set_ylim(0.815, 1.025)
        ins.tick_params(colors=_PAL["dim"], labelsize=5.6)
        for sp in ins.spines.values():
            sp.set_color(_PAL["rule"])
        ins.text(0.5, 0.055, "shape T/(a^2/2 I2), f0 = 1e-04:  I4 is the "
                             "whole signal",
                 transform=ins.transAxes, ha="center", va="bottom",
                 color=_PAL["text"], fontsize=5.9, family="monospace")
        ins.legend(facecolor="#04040e", edgecolor=_PAL["rule"],
                   labelcolor=_PAL["text"], fontsize=5.5, loc="upper right",
                   framealpha=0.9)

    # ── assembly ────────────────────────────────────────────────────────────
    def draw(self) -> "BearingFigure":
        self._profile()
        self._strip()
        self._ring()
        self._cost()
        c = self.c
        self.fig.text(0.5, 0.960,
                      "v69  -  the centre as a finite bearing, not a point",
                      color=_PAL["text"], fontsize=15.5, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.921,
                      "dl^2 = ds^2 + f(s)^2 dOmega^2 with f_min = f0 > 0 -- the "
                      "point is blown up into the space of radial directions, "
                      "S^(d-1) or RP^(d-1)",
                      color=_PAL["dim"], fontsize=8.2, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.890,
                      "THE HINGE SURVIVES REGULARISATION, AND IS CHEAPER THAN "
                      "PROPOSED: the geodesic turn cost is QUADRATIC, "
                      "T(alpha) = alpha^2/(2 I2) + O(alpha^4), not the arc "
                      "f0 alpha -- the saving is PYTHAGORAS, not leverage",
                      color=_PAL["neck"], fontsize=7.6, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.862,
                      "monopole flux and infinitesimal rotation are ONE "
                      "Dirichlet form int w phi'^2 ds -- but the weights match "
                      "only at q = 2 (metric f^2 vs volume f^q), which is the "
                      "physical case and where I2 is also the resistance",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.838,
                      "the LEADING FORM alpha^2/(2 I2) is universal (I2 is "
                      "not); I4 is the first extra moment, where the shape "
                      "shows.  GEOMETRY ONLY -- no field equation is solved.",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.018,
                      f"f0 = {c.neck:g}   f_o = {c.outer:g}   f_i = {c.inner:g}"
                      f"   -   drawn at a deliberately LARGE f0 so the bearing "
                      f"is visible; the working point is f0 = 1e-03",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")
        return self


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--still", default="v69.png")
    ap.add_argument("--neck", type=float, default=0.12)
    ap.add_argument("--outer", type=float, default=1.0)
    ap.add_argument("--inner", type=float, default=0.42)
    ap.add_argument("--dpi", type=int, default=110)
    args = ap.parse_args(argv)
    matplotlib.use("Agg")
    fig = BearingFigure(RegularizedCenter(neck=args.neck, outer=args.outer,
                                          inner=args.inner)).draw()
    fig.fig.savefig(args.still, dpi=args.dpi, facecolor=_PAL["bg"])
    print(f"wrote {args.still}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
