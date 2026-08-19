#!/usr/bin/env python3
"""
Geometrodynamic QED — v62: the balls removed, and the answer made a theorem
===========================================================================

PR #261 answered PR #260's gate — the growing mode does not survive a
finite-radius mouth — and named its own weakest point in doing so: those mouths
were SPHERES IN A FIXED AMBIENT. The balls were never removed, and only l = 0
coupled to the tube.

This round removes them. The ambient is S^3 minus two balls, the tube is glued
along the boundary spheres, and the binary question is re-asked.

The answer is still no, and it is now a THEOREM. With the balls removed there is
no subtraction anywhere, so E = Int (|grad phi|^2 + phi^2) + A Int (|u'|^2 + m^2
|u|^2) is a sum of non-negative terms; E = 0 forces phi = 0, matching forces
u(0) = u(L) = 0, and Poincare finishes it. So lam > 0 for EVERY configuration —
all multipoles, no truncation, no sweep.

What the panels show
--------------------
**Top left — the exterior map has a sign, in every channel.** N_l(lam) computed
by shooting the radial equation from the far pole, positive for every l and
every lam < 0, and INCREASING in l: the monopole is the softest channel, so the
higher ones cannot be the first to go soft. That is why PR #261's monopole
truncation was never a stability limitation.

**Top right — the theorem, and the object it is about.** Rayleigh quotients of
explicit trial configurations, all positive and all above the computed lowest
mode. The degenerate purely-interior case lands on the Poincare floor pi^2/L^2
exactly, which is the case the theorem turns on.

**Bottom left — what the fixed ambient cost.** PR #261's f(a)G(a) against this
round's 1/N_0. They agree to 1.3e-04 at a = 0.02 and part company at large
radius and deep lam, reaching 11% at a = 0.35, lam = -4 — the fraction of the
sphere wrongly left in.

**Bottom right — the soft mode is FORCED, and there can be more than one.**
F_sym runs from +inf at lam = 0 to -inf at lam = 1, so a state exists without
scanning for it. Above L = pi the tube's own harmonics enter the gap and each
brings another — a pole is a sign change with NO zero, and PR #261's "exactly
one state" was a statement about L < pi.

What is put in
--------------
The tube has ONE transverse channel, so A is a coupling and the neck is a
quantum-graph edge rather than a solved cross-section. The ambient metric is
FIXED. **No backreaction.**

Usage
-----
    python scripts/geometrodynamics_v62_neck.py --still out.png
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.waves.finite_mouth import (
    FOUR_PI, mouth_green, regular_radial,
)
from geometrodynamics.waves.neck import (
    NeckThroat,
    exterior_dtn,
    exterior_dtn_monopole,
    measure_the_negative_mode_does_not_survive_the_neck,
    measure_the_quadratic_form_is_positive,
    measure_the_soft_mode_is_forced_by_the_two_ends,
    measure_what_the_fixed_ambient_cost,
)

_PAL = {
    "bg": "#010106", "panel": "#02020a", "grid": "#0e1a2a", "rule": "#1a2838",
    "text": "#e8ecf4", "dim": "#6a8aad",
    "sym": "#7cff9e",      # positive / the exterior
    "anti": "#ffb347",     # the tube
    "hot": "#ff6b8a",      # what was wrong, or what a pole is
    "cool": "#5cc8ff",
}
_ELL_COLORS = ("#7cff9e", "#5cc8ff", "#b388ff", "#ffb347", "#ff6b8a")


class NeckFigure:
    """The exterior map, the theorem, the price of the old model, the count."""

    def __init__(self, figsize=(13.8, 8.6)) -> None:
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, width_ratios=[1.0, 1.0], height_ratios=[1.0, 1.0],
            left=0.068, right=0.975, top=0.845, bottom=0.135,
            wspace=0.26, hspace=0.52)
        self.ax_dtn = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_form = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_cost = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_count = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])
        self._sweep = measure_the_negative_mode_does_not_survive_the_neck()
        self._form = measure_the_quadratic_form_is_positive(trials=40)
        self._cost = measure_what_the_fixed_ambient_cost()
        self._count = measure_the_soft_mode_is_forced_by_the_two_ends()

    # -- the exterior map ----------------------------------------------------
    def _draw_dtn(self) -> None:
        ax = self.ax_dtn
        a = 0.05
        sig = np.geomspace(0.05, 200.0, 220)
        for i, ell in enumerate((0, 1, 2, 3, 5)):
            vals = [exterior_dtn(a, -float(s) ** 2, ell) for s in sig]
            ax.loglog(sig, vals, lw=1.7 if ell == 0 else 1.2,
                      color=_ELL_COLORS[i], zorder=6 - i, label=f"l = {ell}")
        ax.annotate("l = 0 — the only channel a one-mode tube drives",
                    xy=(0.9, 0.30), color=_PAL["sym"], fontsize=6.4,
                    family="monospace")
        ax.annotate("N_l(λ) > 0 in every channel, and INCREASING in l\n"
                    "— the monopole is the softest, so the higher\n"
                    "channels cannot be the first to go soft.\n"
                    "shooting vs the l=0 closed form:  1.7e-14",
                    xy=(0.045, 0.94), xycoords="axes fraction", va="top",
                    color=_PAL["text"], fontsize=6.6, family="monospace",
                    linespacing=1.8)
        ax.axhline(FOUR_PI * a, color=_PAL["dim"], lw=0.9, ls=(0, (2, 4)),
                   zorder=2)
        ax.annotate("4πa — the capacitance", xy=(0.06, FOUR_PI * a * 1.06),
                    color=_PAL["dim"], fontsize=6.4, family="monospace")
        ax.set_ylim(0.15, 90.0)
        ax.set_xlabel("σ,   λ = −σ²", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("exterior DtN  N_l(λ)", color=_PAL["dim"], fontsize=8)
        self._style(ax, "the exterior map has a sign — with the ball gone",
                    legend="lower right")

    # -- the quadratic form --------------------------------------------------
    def _draw_form(self) -> None:
        ax = self.ax_form
        m = self._form
        # the measurement reports only the first eight rows, so the scatter is
        # regenerated here from the same seed and the same recipe
        rng = np.random.default_rng(20260819)
        t = NeckThroat(separation=1.3, length=0.9, area=FOUR_PI, radius=0.05)
        from geometrodynamics.waves.neck import rayleigh_quotient
        vals, fracs = [], []
        for _ in range(40):
            c = rng.normal(size=4)
            k = float(abs(rng.normal()) + 0.4)

            def profile(chi, c=c, k=k, a0=0.05):
                x = (chi - a0) / (math.pi - a0)
                return float(math.exp(-k * x)
                             * (c[0] + c[1] * x + c[2] * x ** 2
                                + c[3] * x ** 3))

            end = profile(0.05)
            b = float(rng.normal())
            r = rayleigh_quotient(
                t, profile,
                lambda s, e=end, b=b: float(e + b * math.sin(math.pi * s / 0.9)),
                n=1201)
            vals.append(r["quotient"])
            fracs.append(r["ambient_norm"] / r["norm"])
        ax.scatter(fracs, vals, s=13, color=_PAL["sym"], alpha=0.85, zorder=6,
                   edgecolors="none", label="trial configurations")
        ax.axhline(m["lowest_computed_mode"], color=_PAL["cool"], lw=1.3,
                   zorder=5,
                   label=f"lowest computed mode  {m['lowest_computed_mode']:.4f}")
        ax.axhline(m["poincare_floor"], color=_PAL["anti"], lw=1.1,
                   ls=(0, (4, 2)), zorder=4,
                   label=f"Poincaré floor π²/L²  {m['poincare_floor']:.3f}")
        ax.axhline(0.0, color=_PAL["hot"], lw=1.0, zorder=3)
        ax.set_yscale("log")
        ax.set_ylim(0.012, 900.0)
        ax.annotate("E = ∫(|∇φ|²+φ²) + A∫(|u'|²+m²|u|²)\n"
                    "every term ≥ 0, nothing is subtracted.\n"
                    "E = 0 => φ ≡ 0 => u(0)=u(L)=0 => Poincaré.\n"
                    "so λ > 0 — all multipoles, no truncation.",
                    xy=(0.045, 0.955), xycoords="axes fraction", va="top",
                    color=_PAL["text"], fontsize=6.6, family="monospace",
                    linespacing=1.8)
        ax.set_xlabel("fraction of the norm outside the neck",
                      color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("Rayleigh quotient  E/‖·‖²", color=_PAL["dim"],
                      fontsize=8)
        self._style(ax, "the theorem — and the object it is about",
                    legend="lower left")

    # -- what the fixed ambient cost -----------------------------------------
    def _draw_cost(self) -> None:
        ax = self.ax_cost
        radii = np.geomspace(0.01, 0.45, 90)
        for lam, colour, tag in ((0.0, _PAL["sym"], "λ = 0"),
                                 (-4.0, _PAL["hot"], "λ = −4")):
            err = [abs(regular_radial(float(a), lam) * mouth_green(float(a), lam)
                       - 1.0 / exterior_dtn_monopole(float(a), lam))
                   / abs(1.0 / exterior_dtn_monopole(float(a), lam))
                   for a in radii]
            ax.loglog(radii, err, lw=1.7, color=colour, zorder=6,
                      label=f"PR #261 vs the neck,   {tag}")
        ax.axvline(0.05, color=_PAL["dim"], lw=0.9, ls=(0, (2, 4)), zorder=2)
        ax.annotate("the working radius", xy=(0.052, 2e-5), color=_PAL["dim"],
                    fontsize=6.4, family="monospace", rotation=90)
        ax.annotate("PR #261 used f(a)G(a) — a shell average taken\n"
                    "with the balls STILL IN.  this round uses 1/N₀,\n"
                    "the inverse exterior map with them gone.\n"
                    "the error is the fraction wrongly left in.",
                    xy=(0.045, 0.94), xycoords="axes fraction", va="top",
                    color=_PAL["text"], fontsize=6.6, family="monospace",
                    linespacing=1.8)
        ax.set_xlabel("mouth radius  a", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("relative disagreement", color=_PAL["dim"], fontsize=8)
        self._style(ax, "what the fixed ambient cost, priced",
                    legend="lower right")

    # -- the count -----------------------------------------------------------
    def _draw_count(self) -> None:
        ax = self.ax_count
        lams = np.linspace(1e-4, 1.0 - 1e-4, 1400)
        for length, colour in ((0.9, _PAL["sym"]), (8.0, _PAL["cool"])):
            t = NeckThroat(separation=1.3, length=length, area=FOUR_PI,
                           radius=0.05)
            vals = np.array([t.channel_functions(float(x))[0] for x in lams])
            ax.plot(lams, np.arcsinh(vals), lw=1.5, color=colour, zorder=6,
                    label=f"F_sym,  L = {length}")
            b = t.bound_states(n=4000)
            for root in b["symmetric_roots"]:
                ax.scatter([root], [0.0], s=34, color=colour, zorder=8,
                           edgecolors=_PAL["bg"], linewidths=0.6)
            for pole in b["symmetric_poles"]:
                ax.axvline(pole, color=_PAL["hot"], lw=1.0, ls=(0, (2, 3)),
                           zorder=4)
        ax.axhline(0.0, color=_PAL["dim"], lw=0.9, zorder=3)
        ax.annotate("+∞ at λ→0⁺  ·  −∞ at λ→1⁻  =>  a state is FORCED\n"
                    "(N₀ → 2π(π−a+sin a cos a)(1−λ) at the gap edge)\n"
                    "dashed red = a POLE: a sign change with NO zero.\n"
                    "above L = π the tube's harmonics enter the gap,\n"
                    "so PR #261's 'exactly one state' meant L < π.",
                    xy=(0.045, 0.33), xycoords="axes fraction", va="top",
                    color=_PAL["text"], fontsize=6.4, family="monospace",
                    linespacing=1.8)
        ax.set_ylim(-11.0, 8.0)
        ax.set_xlabel("λ   (0 to the free ESU gap at 1)", color=_PAL["dim"],
                      fontsize=8)
        ax.set_ylabel("asinh F_sym(λ)", color=_PAL["dim"], fontsize=8)
        self._style(ax, "the soft mode is forced — and there can be more",
                    legend="upper right")

    # -- shared styling ------------------------------------------------------
    def _style(self, ax, title: str, legend: str) -> None:
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc=legend, fontsize=6.3, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], which="both")
        ax.set_title(title, color=_PAL["text"], fontsize=8.6, pad=6)

    # -- frame ---------------------------------------------------------------
    def draw(self) -> None:
        self._draw_dtn()
        self._draw_form()
        self._draw_cost()
        self._draw_count()
        self.fig.suptitle("v62 — THE BALLS REMOVED, AND THE ANSWER MADE A "
                          "THEOREM", color=_PAL["text"], fontsize=13.2,
                          y=0.962, family="monospace")
        self.fig.text(0.5, 0.908,
                      "PR #261 answered the gate with spheres in a FIXED "
                      "ambient   ·   this round removes them, and the answer "
                      "does not move",
                      color=_PAL["dim"], fontsize=8.4, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.070,
                      "with the balls removed there is NO SUBTRACTION "
                      "anywhere, so the energy is a sum of non-negative terms "
                      "— positivity of the form, not a sign on a 2×2",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.047,
                      f"{self._sweep['samples']} sweep samples agree with the "
                      f"theorem ({self._sweep['positives']} roots, worst "
                      f"approach {self._sweep['worst_approach_to_zero']:.1e})  "
                      "·  and one correction: 'exactly one state below the "
                      "gap' was a statement about L < π",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.022,
                      "conformally coupled scalar on a FIXED S³ with two balls "
                      "cut out   ·   ONE transverse tube channel, so A is a "
                      "coupling, not a solved cross-section   ·   NO "
                      "backreaction",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")


def still(path: str) -> None:
    fig = NeckFigure()
    fig.draw()
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v62.png")
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
