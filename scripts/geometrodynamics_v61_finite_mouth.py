#!/usr/bin/env python3
"""
Geometrodynamic QED — v61: the negative mode does not survive a finite mouth
============================================================================

PR #260 gated the roadmap on one question: does its growing mode survive a
finite-radius mouth? It does not, and the statement is STRUCTURAL. At lam < 0
the tube's channel functions are strictly negative and the resolved mouth's are
strictly positive, so their difference has no zero for ANY parameter choice.

What produced PR #260's mode was the LINEARIZATION of the mouth's self-energy:
freezing a screened quantity at its unscreened leading term. The mode it
manufactured sits at kappa*a = 1, exactly the edge of that approximation.

What the panels show
--------------------
**Top left - the two sides never meet.** The tube's symmetric channel against
the ambient's, down the imaginary axis. One is negative everywhere, the other
positive everywhere, so det(A - G) = 0 has no root there. The linearized mouth's
ambient curve is drawn too: it crosses, and it crosses at kappa*a = 1.

**Top right - the screening is the whole difference.** The exact G(a,-kappa^2)
decays like e^{-kappa a}; PR #260's constant 1/(4 pi a) does not. The two agree
to 0.8% for kappa*a < 0.1 and disagree in SIGN beyond it.

**Bottom left - where the mode went.** Not away: below the gap and above zero,
at lam_0 -> 8 pi a/(A L), two mouth capacitances restoring a tube of volume A L.
The point limit drives it to zero FROM ABOVE; the linearized mouth put it on the
other side.

**Bottom right - the delay survives.** Slope 1.0010 in L, saturating at the
ambient path d exactly, with the mouth adding only a sub-leading O(a) shift.

What is put in
--------------
Mouths are SPHERES in a fixed ambient coupled through ONE channel each, so only
l = 0 talks to the tube; the dropped multipoles are screened as (a/d)^l, PR
#250's law. Not a solved neck geometry. **No backreaction.**

Usage
-----
    python scripts/geometrodynamics_v61_finite_mouth.py --still out.png
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.waves.finite_mouth import (
    FOUR_PI,
    FiniteMouthThroat,
    measure_the_delay_survives_with_a_radius_correction,
    measure_the_instability_was_the_linearization,
    measure_the_mode_became_soft_and_positive,
    measure_the_negative_mode_does_not_survive,
    mouth_green,
)

_PAL = {
    "bg": "#010106", "panel": "#02020a", "grid": "#0e1a2a", "rule": "#1a2838",
    "text": "#e8ecf4", "dim": "#6a8aad",
    "sym": "#7cff9e",      # the ambient, resolved
    "anti": "#ffb347",     # the tube
    "hot": "#ff6b8a",      # the linearization / what was wrong
    "cool": "#5cc8ff",
}


class FiniteMouthFigure:
    """Signs, screening, the soft mode, and the surviving delay."""

    def __init__(self, figsize=(13.8, 8.6)) -> None:
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, width_ratios=[1.0, 1.0], height_ratios=[1.0, 1.0],
            left=0.068, right=0.975, top=0.845, bottom=0.135,
            wspace=0.26, hspace=0.52)
        self.ax_sign = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_scr = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_soft = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_del = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])
        self._sweep = measure_the_negative_mode_does_not_survive()
        self._lin = measure_the_instability_was_the_linearization()
        self._soft = measure_the_mode_became_soft_and_positive()
        self._delay = measure_the_delay_survives_with_a_radius_correction()

    # -- the signs -----------------------------------------------------------
    def _draw_signs(self) -> None:
        ax = self.ax_sign
        a = 0.05
        t = FiniteMouthThroat(separation=1.3, length=0.9, area=FOUR_PI,
                              radius=a)
        sig = np.geomspace(0.05, 300.0, 500)
        tube = np.array([t.signed_parts(float(s))[0][0] for s in sig])
        amb = np.array([t.signed_parts(float(s))[1][0] for s in sig])
        lin = np.array([t.linearized_channels(float(s))[0] for s in sig])
        ax.semilogx(sig, tube, lw=1.8, color=_PAL["anti"], zorder=6,
                    label="the TUBE   (always < 0)")
        ax.semilogx(sig, amb, lw=1.8, color=_PAL["sym"], zorder=6,
                    label="the resolved MOUTH   (always > 0)")
        ax.semilogx(sig, lin, lw=1.3, color=_PAL["hot"], ls=(0, (4, 2)),
                    zorder=5, label="PR #260:  A − (1/4πa + g)   (crosses)")
        ax.axhline(0.0, color=_PAL["dim"], lw=0.9, zorder=3)
        ax.axvline(1.0 / a, color=_PAL["hot"], lw=0.9, ls=(0, (2, 4)), zorder=2)
        ax.annotate("κ = 1/a", xy=(1.05 / a, -0.75), color=_PAL["hot"],
                    fontsize=6.6, family="monospace")
        ax.annotate("a difference of a negative and a positive\n"
                    "number has no zero — so no parameter\n"
                    f"choice can grow.  {self._sweep['samples']} samples, "
                    f"{self._sweep['positives']} roots.",
                    xy=(0.04, 0.30), xycoords="axes fraction", va="top",
                    color=_PAL["text"], fontsize=6.6, family="monospace",
                    linespacing=1.8)
        ax.set_xlabel("σ,   λ = −σ²", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("symmetric channel", color=_PAL["dim"], fontsize=8)
        ax.set_ylim(-1.1, 1.1)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="lower right", fontsize=6.3, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], which="both")
        ax.set_title("the two sides never meet", color=_PAL["text"],
                     fontsize=8.6, pad=6)

    # -- the screening -------------------------------------------------------
    def _draw_screening(self) -> None:
        ax = self.ax_scr
        a = 0.05
        kap = np.geomspace(0.05, 300.0, 500)
        exact = np.array([mouth_green(a, -float(k) ** 2) for k in kap])
        const = np.full_like(kap, 1.0 / (FOUR_PI * a))
        ax.loglog(kap * a, np.maximum(exact, 1e-8), lw=1.8, color=_PAL["sym"],
                  zorder=6, label="exact  G(a, −κ²)   ~ e^{−κa}/4πa")
        ax.loglog(kap * a, const, lw=1.3, color=_PAL["hot"], ls=(0, (4, 2)),
                  zorder=5, label="PR #260's constant  1/(4πa)")
        ax.axvline(1.0, color=_PAL["dim"], lw=0.9, ls=(0, (2, 4)), zorder=2)
        ax.annotate("κa = 1", xy=(1.1, 2e-6), color=_PAL["dim"], fontsize=6.6,
                    family="monospace")
        ax.annotate("the mouth's self-energy is SCREENED.\n"
                    "freezing it at the unscreened value is\n"
                    "what manufactured the growing mode —\n"
                    "and the mode sits at κa ≈ 1, the edge\n"
                    "of the approximation that produced it.",
                    xy=(0.04, 0.40), xycoords="axes fraction", va="top",
                    color=_PAL["dim"], fontsize=6.4, family="monospace",
                    linespacing=1.8)
        roots = self._lin["linearized_roots_times_radius"]
        ax.plot(roots, [1.0 / (FOUR_PI * a)] * len(roots), "X", ms=9.0,
                color=_PAL["hot"], mec=_PAL["bg"], mew=0.8, zorder=7,
                label="where the linearized model roots")
        ax.set_xlabel("κ·a", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("mouth self-energy", color=_PAL["dim"], fontsize=8)
        ax.set_ylim(1e-7, 1e2)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper right", fontsize=6.3, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], which="both")
        ax.set_title("the screening is the whole difference",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # -- the soft mode -------------------------------------------------------
    def _draw_soft(self) -> None:
        ax = self.ax_soft
        rows = self._soft["rows"]
        radii = np.array([r["radius"] for r in rows])
        lam = np.array([r["lambda_0"] for r in rows])
        closed = np.array([r["closed_form"] for r in rows])
        ax.loglog(radii, lam, "o", ms=8.0, color=_PAL["sym"], mec=_PAL["bg"],
                  mew=0.8, zorder=6, label="λ₀,  resolved mouth")
        ax.loglog(radii, closed, lw=1.3, color=_PAL["dim"], ls=(0, (4, 3)),
                  zorder=4, label="8πa/(AL)  — two capacitances")
        lin = 1.0 / radii ** 2
        ax.loglog(radii, lin, lw=1.3, color=_PAL["hot"], ls=(0, (2, 3)),
                  zorder=4, label="|λ| of PR #260's mode,  ~1/a²   (λ < 0)")
        ax.axhline(1.0, color=_PAL["cool"], lw=0.9, ls=(0, (3, 3)), zorder=3)
        ax.annotate("the free ESU gap, λ = 1", xy=(0.006, 1.2),
                    color=_PAL["cool"], fontsize=6.4, family="monospace")
        ax.annotate("the mode did not vanish — it went SOFT.\n"
                    "the point limit drives it to zero FROM\n"
                    "ABOVE; PR #260 put it on the other side.",
                    xy=(0.04, 0.55), xycoords="axes fraction", va="top",
                    color=_PAL["text"], fontsize=6.5, family="monospace",
                    linespacing=1.8)
        ax.set_xlabel("mouth radius  a", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("|λ| of the softest mode", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="lower right", fontsize=6.2, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], which="both")
        ax.set_title("where the mode went: soft, and positive",
                     color=_PAL["text"], fontsize=8.6, pad=6)

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
        ax.plot(ls, on, "o", ms=9.0, color=_PAL["sym"], mec=_PAL["bg"],
                mew=0.8, zorder=6, label="resolved mouth")
        ax.axvline(d, color=_PAL["anti"], lw=1.0, ls=(0, (2, 4)), zorder=2)
        ax.annotate("the ambient path, d", xy=(d + 0.07, 0.12),
                    color=_PAL["anti"], fontsize=6.5, family="monospace")
        ax.annotate(f"slope in L:   "
                    f"{self._delay['slope_in_length']:.4f}   (predicted 1)\n"
                    f"spread above d: "
                    f"{self._delay['onset_spread_above_d']:.1e}\n"
                    f"slope in a:   "
                    f"{self._delay['slope_in_radius']:.2f}   (sub-leading)\n"
                    f"contour ε = {self._delay['contour']:.1f}, where PR #260\n"
                    f"needed ε > σ* ≈ 2 to clear the mode",
                    xy=(0.05, 0.62), xycoords="axes fraction", va="top",
                    color=_PAL["text"], fontsize=6.5, family="monospace",
                    linespacing=1.8)
        ax.set_xlabel("tube proper length  L", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("onset of the cross-mouth response", color=_PAL["dim"],
                      fontsize=8)
        ax.set_xlim(-0.15, 3.35)
        ax.set_ylim(-0.05, 1.45)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="lower right", fontsize=6.3, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title("and the traversal delay survives", color=_PAL["text"],
                     fontsize=8.6, pad=6)

    # -- frame ---------------------------------------------------------------
    def draw(self) -> None:
        self._draw_signs()
        self._draw_screening()
        self._draw_soft()
        self._draw_delay()
        self.fig.suptitle("v61 — THE NEGATIVE MODE DOES NOT SURVIVE A FINITE "
                          "MOUTH", color=_PAL["text"], fontsize=13.2, y=0.962,
                          family="monospace")
        self.fig.text(0.5, 0.908,
                      "PR #260 gated the roadmap on this question   ·   the "
                      "answer is structural, not parametric",
                      color=_PAL["dim"], fontsize=8.4, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.070,
                      "a sphere needs no 1/(4πχ) subtraction, so its "
                      "self-energy is the UNSUBTRACTED G(a,λ) — positive and "
                      "screened — where a point keeps the renormalized g(λ) < 0",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.047,
                      "the mean-value identities it rests on are checked "
                      "against quadrature on S³, not assumed:  ⟨G(·,c₂)⟩ = "
                      "f(a)G(d) to 1.0e-10,  ⟨⟨G⟩⟩ = f(a)G(a) to 4.1e-04 "
                      "(grid-limited)",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.022,
                      "conformally coupled scalar on a fixed ESU   ·   spheres "
                      "in a fixed ambient with MONOPOLE coupling, not a solved "
                      "neck   ·   NO backreaction",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")


def still(path: str) -> None:
    fig = FiniteMouthFigure()
    fig.draw()
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v61.png")
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
