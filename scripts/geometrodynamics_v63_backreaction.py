#!/usr/bin/env python3
"""
Geometrodynamic QED - v63: A + B does what rescaling A or B cannot
=================================================================

PR #260 gated backreaction on a growing mode; #261 and #262 closed that gate.
So this round asks the roadmap's first GR question -- and deliberately not
"does spacetime pinch off?":

    does A + B produce a metric response that rescaling A or B alone cannot?

The structure makes it a linear-algebra question. The field equation is linear
so the fields add; T is quadratic so T[A+B] = T[A] + T[B] + dT with dT
bilinear; linearized Einstein is linear so the responses add. Rescaling A -> cA
sends beta_A -> c^2 beta_A, so everything reachable is the two-parameter family
{c^2 beta_A + d^2 beta_B}, and the question is whether beta[dT] lies in it.

What the panels show
--------------------
**Top left - the channel is never on resonance.** The response is a driven
oscillator, so its transfer function is exactly 1/(w3^2 - w^2), a filter at
w3 = 2 sqrt(2). That number is DERIVED here, not quoted: Cartan about the ESU
in the homogeneous anisotropy gives dG^TT = beta'' + (8/a^2) beta. The measured
source spectrum is drawn against it, with the integers marked -- the conformal
scalar on the ESU has spectrum w_n = n+1 and, on a compact static space, rings
on it forever; T is quadratic and integers are closed under sums and
differences, so the source rings on integers too. But 2 sqrt(2) is IRRATIONAL,
0.172 from the nearest integer, so this channel is off resonance BY
CONSTRUCTION. (A first draft of this round said instead that T being quadratic
puts the power at 2*w0 and chose the carrier to match. That is wrong: the
dominant peak sits at w = 6 whatever the carrier is, because the ringing is the
background's, not the pulse's.)

**Top right - the control, which is the point.** A first attempt at this round
reported 0.982 unreachable and it was PURE QUADRATURE NOISE: independent rules
for the same quantity correlated at -0.04. Two refinement levels now agree in
correlation AND magnitude. Nothing here would be quotable without this panel.

**Bottom left - the three responses.** beta_A, beta_B and the interference
beta_x, in one shear component. The interference term is comparable in size to
the single-wave ones, not a small correction to them.

**Bottom right - the residual.** beta_x against the best fit that any rescaling
of A and B could produce. What is left over is the answer.

What is put in
--------------
The n = 3 harmonic only -- the homogeneous shear, not the full TT tower. A FIXED
ESU, POINT sources and PR #257's POINT throat rather than the resolved mouths of
#261/#262. The response is LINEAR: first-order response, not a solved coupled
system.

Usage
-----
    python scripts/geometrodynamics_v63_backreaction.py --still out.png
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.waves.backreaction import (
    TENSOR_MODE_FREQUENCY, ShearQuadrature, WORKING_BACKREACTION,
    shear_response, unreachable_fraction,
)

_PAL = {
    "bg": "#010106", "panel": "#02020a", "grid": "#0e1a2a", "rule": "#1a2838",
    "text": "#e8ecf4", "dim": "#6a8aad",
    "a": "#7cff9e", "b": "#5cc8ff", "x": "#ffb347",
    "hot": "#ff6b8a", "cool": "#b388ff",
}

_COARSE = ShearQuadrature(bulk=(16, 10, 20), ball=(8, 6, 12))
_FINE = ShearQuadrature(bulk=(20, 12, 24), ball=(10, 8, 16))
_WINDOW = (4.0, 30.0)


class BackreactionFigure:
    """The channel, the control, the responses, and what is left over."""

    def __init__(self, figsize=(13.8, 8.6)) -> None:
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, left=0.068, right=0.975, top=0.845, bottom=0.135,
            wspace=0.26, hspace=0.52)
        self.ax_ch = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_cv = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_bt = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_rs = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])

        s = WORKING_BACKREACTION
        self.setup = s
        self.dt = s.grid.dt
        self.times = s.grid.times
        self.coarse = s.shear_sources(_COARSE)
        self.fine = s.shear_sources(_FINE)
        self.b_coarse = {k: shear_response(self.coarse[k], self.dt)
                         for k in ("A", "B", "cross")}
        self.b_fine = {k: shear_response(self.fine[k], self.dt)
                       for k in ("A", "B", "cross")}
        self.res = unreachable_fraction(
            self.b_fine["A"], self.b_fine["B"], self.b_fine["cross"],
            self.times, _WINDOW)

    # -- the channel ---------------------------------------------------------
    def _draw_channel(self) -> None:
        ax = self.ax_ch
        n = self.setup.grid.n
        freqs = np.fft.rfftfreq(n, d=self.dt) * 2.0 * math.pi
        spec = np.abs(np.fft.rfft(self.fine["A"].reshape(n, 9),
                                  axis=0)).sum(axis=1)
        k = (freqs > 0.15) & (freqs < 40.0)
        ax.loglog(freqs[k], spec[k] / spec[k].max(), lw=1.5, color=_PAL["a"],
                  zorder=6, label="measured source spectrum  |S(w)|")
        w = freqs[k]
        tr = 1.0 / np.abs(TENSOR_MODE_FREQUENCY ** 2 - w ** 2 + 1e-9)
        ax.loglog(w, tr / tr.max(), lw=1.3, color=_PAL["cool"],
                  ls=(0, (4, 2)), zorder=5,
                  label="transfer function  1/(w3^2 - w^2)")
        for m in range(1, 14):
            ax.axvline(m, color=_PAL["dim"], lw=0.5, alpha=0.35, zorder=2)
        ax.axvline(TENSOR_MODE_FREQUENCY, color=_PAL["hot"], lw=1.2,
                   ls=(0, (2, 3)), zorder=4)
        ax.annotate(f"w3 = 2*sqrt(2)\n= {TENSOR_MODE_FREQUENCY:.4f}",
                    xy=(TENSOR_MODE_FREQUENCY * 1.06, 1.5e-5),
                    color=_PAL["hot"], fontsize=6.2, family="monospace",
                    linespacing=1.6)
        ax.annotate("faint lines: the integers.  the conformal scalar on the\n"
                    "ESU has spectrum w_n = n+1 and rings on it forever, and T\n"
                    "is quadratic, so the source rings on integers too.  but\n"
                    "w3 = 2*sqrt(2) is IRRATIONAL -- 0.172 from the nearest --\n"
                    "so this channel is off resonance BY CONSTRUCTION.",
                    xy=(0.035, 0.235), xycoords="axes fraction", va="top",
                    color=_PAL["text"], fontsize=6.0, family="monospace",
                    linespacing=1.75)
        ax.set_ylim(3e-6, 8.0)
        ax.set_xlabel("w", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("normalised", color=_PAL["dim"], fontsize=8)
        self._style(ax, "the channel is never on resonance", "upper right")

    # -- the control ---------------------------------------------------------
    def _draw_control(self) -> None:
        ax = self.ax_cv
        sl = (self.times > _WINDOW[0]) & (self.times < _WINDOW[1])
        names = ("A", "B", "cross")
        corr, ratio = [], []
        for key in names:
            v1 = self.b_coarse[key][sl].reshape(-1)
            v2 = self.b_fine[key][sl].reshape(-1)
            corr.append(float(v1 @ v2 / np.linalg.norm(v1)
                              / np.linalg.norm(v2)))
            ratio.append(float(np.linalg.norm(v2) / np.linalg.norm(v1)))
        xs = np.arange(3)
        ax.bar(xs - 0.18, corr, width=0.34, color=_PAL["a"], zorder=6,
               label="correlation between refinements")
        ax.bar(xs + 0.18, ratio, width=0.34, color=_PAL["b"], zorder=6,
               label="magnitude ratio")
        ax.axhline(1.0, color=_PAL["dim"], lw=0.9, zorder=4)
        ax.axhline(-0.04, color=_PAL["hot"], lw=1.2, ls=(0, (4, 2)), zorder=5,
                   label="the first draft's correlation:  -0.04")
        ax.set_xticks(xs)
        ax.set_xticklabels([f"beta_{n}" for n in names])
        ax.set_ylim(-0.25, 2.05)
        ax.annotate("a first attempt reported 0.982 unreachable and it was\n"
                    "PURE QUADRATURE NOISE.  nothing about the run looked\n"
                    "wrong -- the tell only appeared on recomputing the same\n"
                    "quantity with an independent rule.",
                    xy=(0.035, 0.985), xycoords="axes fraction", va="top",
                    color=_PAL["text"], fontsize=6.0, family="monospace",
                    linespacing=1.75)
        ax.set_ylabel("agreement", color=_PAL["dim"], fontsize=8)
        self._style(ax, "the control - without which none of this is quotable",
                    "center right")

    # -- the responses -------------------------------------------------------
    def _draw_betas(self) -> None:
        ax = self.ax_bt
        k = (self.times > 0.0) & (self.times < 32.0)
        t = self.times[k]
        for key, col, lab in (("A", _PAL["a"], "beta_A"),
                              ("B", _PAL["b"], "beta_B"),
                              ("cross", _PAL["x"], "beta_x  (interference)")):
            ax.plot(t, self.b_fine[key][k, 0, 1], lw=1.4, color=col, zorder=6,
                    label=lab)
        ax.axhline(0.0, color=_PAL["dim"], lw=0.8, zorder=3)
        lim = 1.9 * float(np.max(np.abs(self.b_fine["cross"][k, 0, 1])))
        ax.set_ylim(-lim, lim)
        ax.annotate("the interference response is COMPARABLE to the\n"
                    f"single-wave ones -- |beta_x|/|beta_A| = "
                    f"{self.res['cross_over_single']:.2f} -- not a\n"
                    "small correction to them.",
                    xy=(0.035, 0.985), xycoords="axes fraction", va="top",
                    color=_PAL["text"], fontsize=6.0, family="monospace",
                    linespacing=1.75)
        ax.set_xlabel("t", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("shear response  beta_12", color=_PAL["dim"], fontsize=8)
        self._style(ax, "the three responses", "lower left")

    # -- the residual --------------------------------------------------------
    def _draw_residual(self) -> None:
        ax = self.ax_rs
        sl = (self.times > _WINDOW[0]) & (self.times < _WINDOW[1])
        t = self.times[sl]
        va = self.b_fine["A"][sl].reshape(len(t), 9)
        vb = self.b_fine["B"][sl].reshape(len(t), 9)
        vx = self.b_fine["cross"][sl].reshape(len(t), 9)
        m = np.stack([va.reshape(-1), vb.reshape(-1)], axis=1)
        coef, *_ = np.linalg.lstsq(m, vx.reshape(-1), rcond=None)
        fit = (coef[0] * va + coef[1] * vb)[:, 1]
        ax.plot(t, vx[:, 1], lw=1.5, color=_PAL["x"], zorder=6,
                label="beta_x  (what A+B actually does)")
        ax.plot(t, fit, lw=1.2, color=_PAL["cool"], ls=(0, (4, 2)), zorder=5,
                label="best fit from rescaling A and B")
        ax.plot(t, vx[:, 1] - fit, lw=1.1, color=_PAL["hot"], zorder=7,
                label="what no rescaling can reach")
        ax.axhline(0.0, color=_PAL["dim"], lw=0.8, zorder=3)
        lim = 2.0 * float(np.max(np.abs(vx[:, 1])))
        ax.set_ylim(-lim, lim)
        ax.annotate(f"UNREACHABLE FRACTION = "
                    f"{self.res['unreachable_off_the_span']:.3f}\n"
                    "measured off the full linear SPAN, which strictly\n"
                    "contains the physical cone -- so this is the\n"
                    "conservative figure.  it is not a constant: 0.88-1.00\n"
                    "over windows, 0.56-0.99 over carriers.",
                    xy=(0.035, 0.985), xycoords="axes fraction", va="top",
                    color=_PAL["text"], fontsize=6.0, family="monospace",
                    linespacing=1.75)
        ax.set_xlabel("t", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("shear response  beta_12", color=_PAL["dim"], fontsize=8)
        self._style(ax, "what is left over is the answer", "lower left")

    # -- shared --------------------------------------------------------------
    def _style(self, ax, title: str, legend: str) -> None:
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc=legend, fontsize=6.3, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], which="both")
        ax.set_title(title, color=_PAL["text"], fontsize=8.6, pad=6)

    def draw(self) -> None:
        self._draw_channel()
        self._draw_control()
        self._draw_betas()
        self._draw_residual()
        self.fig.suptitle("v63 - A + B DOES WHAT RESCALING A OR B CANNOT",
                          color=_PAL["text"], fontsize=13.2, y=0.962,
                          family="monospace")
        self.fig.text(0.5, 0.908,
                      "the gate #260 set is closed, so this is the roadmap's "
                      "first GR question   -   and it has a number",
                      color=_PAL["dim"], fontsize=8.4, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.070,
                      "rescaling A -> cA sends beta_A -> c^2 beta_A, so "
                      "everything reachable is {c^2 beta_A + d^2 beta_B}   -   "
                      "the question is whether beta_x lies in it",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.047,
                      "the TT channel is the only one whose answer does not "
                      "depend on the fluid holding the ESU static: a perfect "
                      "fluid carries no anisotropic stress",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.022,
                      "n = 3 harmonic only   -   FIXED ESU, POINT sources and "
                      "#257's POINT throat, not #261/#262's resolved mouths   "
                      "-   LINEAR response, not a solved coupled system",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")


def still(path: str) -> None:
    fig = BackreactionFigure()
    fig.draw()
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v63.png")
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
