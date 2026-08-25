#!/usr/bin/env python3
"""
Geometrodynamic QED - v71: the first evolved Einstein equations
===============================================================

Every gravity result before this one is stationary, weak-field or linearized.
This is D = 5, spherical symmetry, one massless scalar, evolved in ingoing
Eddington-Finkelstein coordinates -- horizon-penetrating by construction.

What the panels show
--------------------
**Top left - the constraint the code never solves.** The hierarchy SOLVES the rr
and vr Einstein equations on every slice, so their residuals are identically
zero and testing them would be circular. vv is the one independent component
left over, it carries d_v A, and the code never forms d_v A for any other
purpose. Its residual converges at SECOND ORDER: 1.989 -> 1.997 -> 1.999 over a
four-fold refinement. This is the characteristic-scheme analogue of a
Hamiltonian/momentum constraint test and is labelled as an analogue, not as one.

**Top right - two exact solutions, before anything is claimed.** Tangherlini
comes back at machine precision (1.6e-15) because with phi = 0 the A quadrature
is int 2s ds = r^2 exactly. And phi = cos(w(v-r)) J_1(wr)/r solves the flat D = 5
wave equation in these coordinates, so d_v phi is known in closed form -- the psi
quadrature converges at rate 2.003 against it. That is what pins the scheme's
order; nothing else does.

**Bottom left - A REGULAR CENTRE FORBIDS A TRAPPED SURFACE.** The vr equation is
an exact quadrature, r^(n-1) e^delta A = (n-1) int_0^r s^(n-2) e^delta ds: a
positive integrand over a positive interval. So A > 0 strictly, IDENTICALLY, and
no trapped surface can sit on a regular-centred ingoing null slice. Four profile
families driven to min A = 5.6e-03 confirm the code obeys the proof. HORIZON
FORMATION IS NOT OBSERVABLE IN THIS GAUGE -- a statement about the chart, not
the physics. Collapse still happens; the trapped region is reached once the
centre stops being regular, which is why production characteristic codes use
OUTGOING null cones or excise the centre.

**Bottom right - what this round did NOT earn.** Two horizon-penetrating
time-domain constructions for a test scalar on a fixed background. Both stable,
both CONVERGED, and they disagree: real parts within 0.3% at l = 1, damping
rates apart by 37%. So no quasinormal frequency is reported and the transfer
function is not built. A converged number is not a correct number.

What is put in
--------------
CLASSICAL, SPHERICALLY SYMMETRIC, ONE MASSLESS SCALAR, second order and stated
as such. No matter model from the rest of the arc; no charge, no winding, no S^3
harmonics above the monopole. Horizon PERSISTENCE is shown only on a seeded
background, where it is exact. The perturbation spectrum and the retarded
outer-to-inner transfer function are NOT delivered.

Usage
-----
    python scripts/geometrodynamics_v71_tangherlini_dynamics.py --still out.png
"""

from __future__ import annotations

import argparse
import warnings
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.tangherlini.dynamics import (evolve, hierarchy,
                                                   tangherlini_A)

_PAL = {
    "bg": "#010106", "panel": "#02020a", "grid": "#0e1a2a", "rule": "#1a2838",
    "text": "#e8ecf4", "dim": "#6a8aad",
    "cool": "#5cc8ff", "warm": "#ffb347", "good": "#7cff9e",
    "bad": "#ff6b8a", "kern": "#b388ff", "hot": "#ff5ec7",
}

# measured, and carried here rather than recomputed at render time
_CONSTRAINT = [(400, 0.05013, 1.55108e-04), (800, 0.02503, 3.90695e-05),
               (1600, 0.01251, 9.78623e-06), (3200, 0.00625, 2.44780e-06)]
_FLATMODE = [(400, 1.470e-03), (800, 3.649e-04), (1600, 9.078e-05),
             (3200, 2.262e-05), (6400, 5.644e-06)]
_TRAPPED = [("centred gaussian", 12.0, 7.256e-02), ("thin shell r=2", 5.0, 1.513e-02),
            ("thin shell r=2", 12.0, 5.627e-03), ("r^2 e^{-r^2/2}", 12.0, 2.703e-02),
            ("oscillatory", 12.0, 5.743e-03)]


def _style(ax, title, xlabel, ylabel):
    ax.set_title(title, color=_PAL["text"], fontsize=9.6, pad=8,
                 family="monospace")
    ax.set_xlabel(xlabel, color=_PAL["dim"], fontsize=8, family="monospace")
    ax.set_ylabel(ylabel, color=_PAL["dim"], fontsize=8, family="monospace")
    ax.tick_params(colors=_PAL["dim"], labelsize=7.2)
    for s in ax.spines.values():
        s.set_color(_PAL["rule"])
    ax.grid(True, color=_PAL["grid"], lw=0.5, alpha=0.7)


class DynamicsFigure:
    """Constraint convergence, exactness, the positivity bound, and the miss."""

    def __init__(self, figsize=(13.8, 8.9)) -> None:
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, left=0.075, right=0.975, top=0.800, bottom=0.075,
            wspace=0.26, hspace=0.46)
        self.ax_con = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_ex = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_pos = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_miss = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])

    # ── the monitored constraint ────────────────────────────────────────────
    def _constraint(self) -> None:
        ax = self.ax_con
        h = np.array([x[1] for x in _CONSTRAINT])
        e = np.array([x[2] for x in _CONSTRAINT])
        ax.loglog(h, e, "o-", color=_PAL["good"], lw=2.0, ms=7,
                  label="max |vv residual|")
        ref = e[0]*(h/h[0])**2
        ax.loglog(h, ref, ":", color=_PAL["dim"], lw=1.4, label="exact h^2")
        for (N, hh, ee), i in zip(_CONSTRAINT, range(4)):
            if i:
                rate = np.log2(_CONSTRAINT[i-1][2]/ee)
                ax.annotate(f"{rate:.3f}", xy=(hh, ee), xytext=(-46, 11),
                            textcoords="offset points", color=_PAL["good"],
                            fontsize=7.4, family="monospace")
        ax.text(0.97, 0.05,
                "rr and vr are SOLVED on every slice -- their residuals\n"
                "are identically zero and testing them is circular.\n"
                "vv is the one component left over, and it carries d_v A.",
                transform=ax.transAxes, color=_PAL["dim"], fontsize=7.0,
                family="monospace", va="bottom", ha="right")
        _style(ax, "the Einstein equation the code never solves",
               "grid spacing   h", "max | vv residual |")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=7.2, loc="upper left",
                  framealpha=0.93)

    # ── the two exact solutions ─────────────────────────────────────────────
    def _exact(self) -> None:
        ax = self.ax_ex
        r = np.linspace(0.5, 9.0, 700)
        hh = float(r[1]-r[0])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, A, _, _ = hierarchy(r, hh, np.zeros(r.size),
                                   interior_mass=r[0]**2 - 1.0)
        ax.plot(r, tangherlini_A(r, 1.0), lw=3.2, color=_PAL["dim"],
                label="exact  A = 1 - r_h^2/r^2")
        ax.plot(r[::12], A[::12], "o", ms=4.5, color=_PAL["cool"],
                label="the quadrature  (err 1.6e-15)")
        ax.axhline(0.0, color=_PAL["rule"], lw=0.9)
        ax.axvline(1.0, color=_PAL["bad"], lw=1.1, ls=":")
        ax.text(1.06, -0.55, "horizon\nr = r_h", color=_PAL["bad"], fontsize=7.0,
                family="monospace", va="center")
        ax.set_ylim(-1.6, 1.15)
        axi = ax.inset_axes([0.47, 0.27, 0.48, 0.29])
        hs = np.array([12.0/(N-1) for N, _ in _FLATMODE])
        es = np.array([x[1] for x in _FLATMODE])
        axi.loglog(hs, es, "o-", color=_PAL["warm"], lw=1.6, ms=4)
        axi.loglog(hs, es[0]*(hs/hs[0])**2, ":", color=_PAL["dim"], lw=1.2)
        axi.set_title("flat mode: rate 2.003", color=_PAL["warm"], fontsize=6.8,
                      family="monospace", pad=3)
        axi.tick_params(colors=_PAL["dim"], labelsize=5.6)
        for s in axi.spines.values():
            s.set_color(_PAL["rule"])
        axi.set_facecolor(_PAL["panel"])
        _style(ax, "Tangherlini is an exact fixed point",
               "areal radius   r", "metric function   A")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=7.2, loc="upper left",
                  framealpha=0.93)

    # ── the positivity bound ────────────────────────────────────────────────
    def _positivity(self) -> None:
        ax = self.ax_pos
        cols = [_PAL["cool"], _PAL["good"], _PAL["warm"], _PAL["kern"],
                _PAL["hot"]]
        profiles = {
            "thin shell r=2, amp 2": lambda r: 2.0*(np.exp(-((r-2.0)/0.5)**2)
                                                    + np.exp(-((r+2.0)/0.5)**2)),
            "thin shell r=2, amp 5": lambda r: 5.0*(np.exp(-((r-2.0)/0.5)**2)
                                                    + np.exp(-((r+2.0)/0.5)**2)),
            "thin shell r=2, amp 12": lambda r: 12.0*(np.exp(-((r-2.0)/0.5)**2)
                                                      + np.exp(-((r+2.0)/0.5)**2)),
            "oscillatory, amp 12": lambda r: 12.0*np.exp(-r**2/8.0)*np.cos(3.0*r),
        }
        r = np.linspace(0.0, 8.0, 1600)
        hh = float(r[1]-r[0])
        for (name, prof), col in zip(profiles.items(), cols):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                _, A, _, _ = hierarchy(r, hh, prof(r))
            ax.semilogy(r[1:], A[1:], lw=1.8, color=col, label=name)
        ax.axhline(1.0, color=_PAL["rule"], lw=0.9, ls="--")
        ax.set_ylim(1e-3, 2.5e2)
        ax.text(0.5, 0.975,
                "A = 0 IS UNREACHABLE.  the vr quadrature reads\n"
                "  r^(n-1) e^d A = (n-1) int_0^r s^(n-2) e^d ds\n"
                "a positive integrand over a positive interval,\n"
                "so A > 0 strictly for r > 0 -- identically.",
                transform=ax.transAxes, color=_PAL["bad"], fontsize=7.0,
                family="monospace", va="top", ha="center")
        ax.text(0.025, 0.03,
                "min A reached across four profile families: 5.6e-03.\n"
                "so horizon FORMATION is not observable in this chart;\n"
                "the transition is the loss of central regularity.",
                transform=ax.transAxes, color=_PAL["dim"], fontsize=6.9,
                family="monospace", va="bottom")
        _style(ax, "a regular centre forbids a trapped surface",
               "areal radius   r", "metric function   A   (log scale)")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.8, loc="lower right",
                  framealpha=0.93)

    # ── the shortfall ───────────────────────────────────────────────────────
    def _miss(self) -> None:
        ax = self.ax_miss
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        for s in ax.spines.values():
            s.set_color(_PAL["rule"])
        ax.set_title("what this round did NOT earn", color=_PAL["text"],
                     fontsize=9.6, pad=8, family="monospace")
        rows = [("constraint convergence", True, "second order, 1.999"),
                ("horizon formation", None, "resolved as a CHART obstruction"),
                ("horizon persistence", True, "exact, on a seeded background"),
                ("perturbation spectrum", False, "two converged codes disagree"),
                ("retarded transfer function", False, "not built -- same signals")]
        for i, (name, ok, note) in enumerate(rows):
            y = 0.90 - i*0.147
            col = _PAL["good"] if ok else (_PAL["bad"] if ok is False
                                           else _PAL["warm"])
            mark = "YES" if ok else ("NO" if ok is False else "~")
            ax.add_patch(plt.Rectangle((0.045, y-0.052), 0.91, 0.104,
                                       facecolor=col, alpha=0.13,
                                       edgecolor=col, lw=1.1))
            ax.text(0.075, y+0.012, name, color=_PAL["text"], fontsize=8.0,
                    family="monospace", va="center")
            ax.text(0.905, y+0.012, mark, color=col, fontsize=8.0, ha="right",
                    family="monospace", va="center", weight="bold")
            ax.text(0.075, y-0.030, note, color=_PAL["dim"], fontsize=6.8,
                    family="monospace", va="center")
        ax.text(0.5, 0.022,
                "Kerr-Schild  1.01622 - 0.36231i     tortoise  1.01876 - 0.26404i\n"
                "real parts within 0.3%;  DAMPING RATES APART BY 37%\n"
                "A CONVERGED NUMBER IS NOT A CORRECT NUMBER",
                transform=ax.transAxes, ha="center", va="bottom",
                color=_PAL["bad"], fontsize=7.0, family="monospace")

    # ── assembly ────────────────────────────────────────────────────────────
    def draw(self) -> "DynamicsFigure":
        self._constraint()
        self._exact()
        self._positivity()
        self._miss()
        self.fig.text(0.5, 0.962,
                      "v71  -  the first evolved Einstein equations",
                      color=_PAL["text"], fontsize=15.5, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.925,
                      "D = 5, spherical symmetry, one massless scalar, in "
                      "ingoing Eddington-Finkelstein coordinates  -  "
                      "horizon-penetrating by construction",
                      color=_PAL["dim"], fontsize=8.2, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.894,
                      "THE UNUSED EINSTEIN EQUATION CONVERGES AT SECOND ORDER "
                      "(1.989 -> 1.999)  -  AND A REGULAR CENTRE FORBIDS A "
                      "TRAPPED SURFACE, IDENTICALLY",
                      color=_PAL["good"], fontsize=7.6, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.866,
                      "the vr equation is an exact QUADRATURE, not an ODE -- "
                      "which is what makes the geometry machine-precise on "
                      "Tangherlini and what makes A = 0 unreachable",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.840,
                      "CLASSICAL, SPHERICALLY SYMMETRIC, ONE MASSLESS SCALAR, "
                      "SECOND ORDER  -  the perturbation spectrum and the "
                      "retarded transfer function are NOT delivered",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.018,
                      "the bottom-left result is about the CHART, not the "
                      "physics  -  collapse still happens, and the trapped "
                      "region is reached once the centre stops being regular",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")
        return self


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--still", default="v71.png")
    ap.add_argument("--dpi", type=int, default=110)
    args = ap.parse_args(argv)
    matplotlib.use("Agg")
    fig = DynamicsFigure().draw()
    fig.fig.savefig(args.still, dpi=args.dpi, facecolor=_PAL["bg"])
    print(f"wrote {args.still}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
