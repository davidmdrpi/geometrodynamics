#!/usr/bin/env python3
"""
Geometrodynamic QED — v59: the two-wave invariant is branch-resolved
====================================================================

PR #258 built a static kernel and said plainly it was not this: no local null
momenta, so it could not tell equal-energy collinear from counterpropagating
waves. This round solves the time-dependent field on the throated ESU, builds
the improved conformal stress tensors, and applies exactly that control.

The known limit is the control, not the result. WKB says

```
N  =  (T_A:T_B)/(T_A^00 T_B^00)  =  (1 - cos th)^2       0 collinear, 4 head-on
```

and the exact field returns 3.99995 and 1.8e-10. What the round is about is the
difference: **the same two sources at the same event give three different
answers, selected by which branch has arrived.**

What the panels show
--------------------
**Top left - the branch-resolved invariant.** Sources and observer fixed on one
great circle. Direct branches: collinear, N at its numerical floor. A on its
long-way winding image: that branch runs the other way round the sphere, its
arrival direction is reversed, and the same pair reads head-on, N = 4. Then the
explicit `(i, j)` audit: all four two-leg paths, `j` the mouth the source drives
and `i` the mouth the signal leaves from. The prediction depends on `i` alone,
so the four paths carry two values, and the field picks the right one at each
delay to 0.08% - not fitted.

The hollow rings are the control that scopes the whole panel: the same channels
with **`beta = 0`**, the mouths *disconnected*. They land on the connected
values to a part in 1e6. They have to: N is amplitude-normalized, and a single
channel is a single arrival direction, so the connection cannot enter. It is the
*mouths* that create the branch. What sees the connection is still `W = -beta`,
from the low-frequency limit of the same solve.

**Top right - the solved waveform.** The field at the observer, free and with
the throat. The free field lives strictly on the light cones - S3 x R is
conformally flat, so the conformal scalar obeys Huygens exactly - and the throat
adds two-leg arrivals plus a ring-down between them. Every tail in this model is
the throat's.

**Bottom left - how the WKB limit is approached.** Head-on to 4 and collinear to
0, against carrier frequency. The collinear null is far stronger than leading
order because the two wavefronts share their normal exactly; that is what makes
the multipath effect above a large one rather than a competition between small
ones. The under-resolved contour is drawn too, because it is wrong by four
orders and looks plausible.

The same panel carries the interference tensor `dT = T[phi_A + phi_B] - T[phi_A]
- T[phi_B]`, and it says something the invariant does not: `dT^00` normalized by
`sqrt(T_A^00 T_B^00)` is 2.000 collinear - the coherent maximum - against 1.044
head-on. The interference energy is *largest* exactly where `T_A:T_B` is null,
so the invariant is not a proxy for it. A backreaction estimate driven by
`T_A:T_B` would look at the collinear case, see nothing, and be wrong about its
own source by the size of the whole effect.

**Bottom right - the caustic, and where geometric optics stops.** Exact
amplitude against the WKB `1/(4 pi sin chi)` near the antipode, for three
carriers. The ratio is a function of `omega*e` alone - the curves collapse to
`|sin(omega e)|` - so the caustic is cut off at `e* ~ 1/omega` and the exact
amplitude saturates at `omega/4pi` where WKB diverges.

What is put in
--------------
A linear conformally coupled field on a fixed background, the mouth positions,
and the boundary data - four numbers chosen, not derived, strictly inside PR
#257's cone. **No backreaction**: the stress tensor is computed from the field
and never fed back.

Usage
-----
    python scripts/geometrodynamics_v59_two_wave.py --still out.png
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.waves.two_source import geodesic, mouth_positions
from geometrodynamics.waves.two_wave import (
    OBSERVER_REACH,
    SOURCE_GAP,
    TWO_PI,
    GaussianPulse,
    RetardedGrid,
    TwoWaveSetup,
    arrival_directions,
    circle_point,
    green_omega,
    measure_multipath_destroys_the_collinear_null,
    measure_the_cross_mouth_channels_are_labelled_by_the_exit_mouth,
    measure_the_interference_tensor_is_largest_where_the_invariant_is_null,
    measure_the_wkb_collinear_head_on_result_is_recovered,
    solve_field,
    working_pair,
)

_PAL = {
    "bg": "#010106",
    "panel": "#02020a",
    "grid": "#0e1a2a",
    "rule": "#1a2838",
    "text": "#e8ecf4",
    "dim": "#6a8aad",
    "free": "#63798f",
    "sym": "#7cff9e",          # the WKB target / free field
    "anti": "#ffb347",         # the throat
    "hot": "#ff6b8a",          # the winding image / what is wrong
    "cool": "#5cc8ff",
}


class TwoWaveFigure:
    """Branches, waveform, the limit, and the caustic."""

    def __init__(self, figsize=(13.8, 8.6)) -> None:
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, width_ratios=[1.0, 1.05], height_ratios=[1.0, 1.0],
            left=0.065, right=0.975, top=0.845, bottom=0.135,
            wspace=0.24, hspace=0.46)
        self.ax_br = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_wf = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_lim = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_ca = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])
        self._multipath = measure_multipath_destroys_the_collinear_null()
        self._limit = measure_the_wkb_collinear_head_on_result_is_recovered()
        self._audit = (
            measure_the_cross_mouth_channels_are_labelled_by_the_exit_mouth())
        self._cross = (
            measure_the_interference_tensor_is_largest_where_the_invariant_is_null())

    # -- the branch-resolved invariant ---------------------------------------
    def _draw_branches(self) -> None:
        ax = self.ax_br
        rows = self._multipath["rows"]
        labels = ["A direct\n+ B direct", "A LONG-WAY\nimage + B direct",
                  "A direct\n+ B via a mouth"]
        vals = [max(abs(r["invariant"]), 1e-11) for r in rows]
        preds = [max(r["geometric_prediction"], 1e-11) for r in rows]
        cols = [_PAL["sym"], _PAL["hot"], _PAL["anti"]]
        xs = np.arange(3.0)
        for x, v, pr, c in zip(xs, vals, preds, cols):
            ax.vlines(x, 1e-11, v, color=c, lw=1.2, alpha=0.55, zorder=3)
            ax.plot([x], [v], "o", ms=11.0, color=c, mec=_PAL["bg"], mew=0.8,
                    zorder=6)
            ax.plot([x], [pr], "_", ms=22.0, mew=2.4, color=_PAL["text"],
                    alpha=0.85, zorder=7)
            ax.annotate(f"{rows[int(x)]['invariant']:.4g}", xy=(x, v),
                        xytext=(13, 2), textcoords="offset points",
                        ha="left", va="center", color=c, fontsize=7.2,
                        family="monospace")
        # the (i, j) audit: two resolvable channels, each with its β = 0 control
        audit = self._audit["rows"]
        for k, r in enumerate(audit):
            x = 3.0 + k
            v = max(abs(r["measured_invariant"]), 1e-11)
            ax.vlines(x, 1e-11, v, color=_PAL["cool"], lw=1.2, alpha=0.55,
                      zorder=3)
            ax.plot([x], [v], "o", ms=11.0, color=_PAL["cool"],
                    mec=_PAL["bg"], mew=0.8, zorder=6)
            ax.plot([x], [max(r["predicted_invariant"], 1e-11)], "_", ms=22.0,
                    mew=2.4, color=_PAL["text"], alpha=0.85, zorder=7)
            ax.plot([x], [max(abs(r["control_beta_zero"]), 1e-11)], "o",
                    ms=17.0, mfc="none", mec=_PAL["dim"], mew=1.3, zorder=8)
            ax.annotate(f"{r['measured_invariant']:.4g}", xy=(x, v),
                        xytext=(13, 2), textcoords="offset points",
                        ha="left", va="center", color=_PAL["cool"],
                        fontsize=7.2, family="monospace")
        labels += [f"(i={int(r['exit_mouth'])},\nj={int(r['entry_mouth'])})"
                   for r in audit]
        xs = np.arange(3.0 + len(audit))
        ax.plot([], [], "o", ms=8, color=_PAL["dim"],
                label="exact solved field")
        ax.plot([], [], "_", ms=14, mew=2.2, color=_PAL["text"],
                label="geometry:  (1 − n̂_A·n̂_B)²")
        ax.plot([], [], "o", ms=9, mfc="none", mec=_PAL["dim"], mew=1.3,
                label="β = 0 control  (mouths disconnected)")
        ax.set_yscale("log")
        ax.set_ylim(1e-11, 1e5)
        ax.set_xlim(-0.55, 4.75)
        ax.set_xticks(xs)
        ax.set_xticklabels(labels, color=_PAL["dim"], fontsize=6.8)
        ax.axhline(4.0, color=_PAL["dim"], lw=0.9, ls=(0, (3, 3)), zorder=2)
        ax.annotate("WKB head-on = 4", xy=(-0.45, 6.5), color=_PAL["dim"],
                    fontsize=6.4, family="monospace", ha="left")
        ax.annotate("same two sources, same observation point.\n"
                    "only the BRANCH changes.  each (i,j) channel is\n"
                    "predicted by its EXIT mouth i alone — and the\n"
                    "β = 0 control lands on it, to a part in 1e6.",
                    xy=(0.30, 0.03), xycoords="axes fraction", va="bottom",
                    color=_PAL["text"], fontsize=6.5, family="monospace",
                    linespacing=1.7)
        ax.set_ylabel("N = (T_A:T_B)/(T_A⁰⁰ T_B⁰⁰)", color=_PAL["dim"],
                      fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper left", fontsize=6.4, framealpha=0.0,
                        labelcolor=_PAL["dim"], ncol=2)
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], axis="y", which="both")
        ax.set_title("the invariant is branch-resolved",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # -- the solved waveform -------------------------------------------------
    def _draw_waveform(self) -> None:
        ax = self.ax_wf
        grid = RetardedGrid(n=1 << 17, span=600.0, eps=0.05)
        obs = circle_point(SOURCE_GAP + OBSERVER_REACH)
        src = circle_point(0.0)
        width = 0.06
        pulse = GaussianPulse(amplitude=1.0 / (width * math.sqrt(TWO_PI)),
                              carrier=0.0, width=width)
        sols = {}
        for throat in (False, True):
            setup = TwoWaveSetup(pair=working_pair(), source_a=src,
                                 source_b=src, observer=obs, pulse_a=pulse,
                                 pulse_b=pulse, grid=grid,
                                 with_throat=throat)
            sols[throat] = solve_field(setup, src, pulse)
        ts = grid.times
        sl = (ts > 0.2) & (ts < 8.2)
        gain = 5.0
        only = sols[True]["phi"] - sols[False]["phi"]
        ax.plot(ts[sl], sols[False]["phi"][sl], lw=1.5, color=_PAL["sym"],
                zorder=5, label="free field")
        ax.plot(ts[sl], gain * only[sl], lw=1.2, color=_PAL["anti"],
                alpha=0.95, zorder=6,
                label=f"the throat's own field  (×{gain:.0f})")
        chi = geodesic(obs, src)
        for t_arr in (chi, TWO_PI - chi, TWO_PI + chi):
            ax.axvline(t_arr, color=_PAL["free"], lw=0.7, ls=(0, (2, 4)),
                       zorder=2)
        c1, c2 = mouth_positions(working_pair().separation)
        two_leg = sorted({round(geodesic(src, cj) + geodesic(ci, obs), 4)
                          for ci in (c1, c2) for cj in (c1, c2)})
        for t_arr in two_leg:
            if t_arr < 8.2:
                ax.axvline(t_arr, color=_PAL["anti"], lw=0.7,
                           ls=(0, (1, 3)), zorder=2)
        ax.annotate("dotted grey: free branch ledger  χ, 2π−χ, 2π+χ\n"
                    "dotted amber: the throat's two-leg arrivals\n"
                    "between them the free field is 1.4e-08 of its peak\n"
                    "(Huygens is exact here) and the throat's is 8.1e-02",
                    xy=(0.50, 0.34), xycoords="axes fraction", va="top",
                    color=_PAL["dim"], fontsize=6.3, family="monospace",
                    linespacing=1.75)
        ax.set_xlabel("t", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("φ at the observer", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper right", fontsize=6.6, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title("the solved field, and where it has support",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # -- the WKB limit -------------------------------------------------------
    def _draw_limit(self) -> None:
        ax = self.ax_lim
        rows = self._limit["rows"]
        cs = np.array([r["carrier"] for r in rows])
        head = np.array([abs(r["head_on"] - 4.0) for r in rows])
        col = np.array([abs(r["collinear"]) for r in rows])
        ax.loglog(cs, head, "o-", ms=4.0, lw=1.5, color=_PAL["hot"],
                  zorder=5, label="|N − 4|,  head-on")
        ax.loglog(cs, col, "s-", ms=4.0, lw=1.5, color=_PAL["sym"],
                  zorder=5, label="|N|,  collinear  (WKB says 0)")
        bad = abs(self._limit["under_resolved_contour_value"])
        good = abs(self._limit["converged_value_there"])
        ax.loglog([32.0], [bad], "x", ms=9.0, mew=2.0, color=_PAL["anti"],
                  zorder=6, label="collinear, contour under-resolved")
        ax.annotate(f"ε ≈ 2π/span:  {bad:.1e}\n"
                    f"converged:    {good:.1e}\n"
                    "→ four orders, and it looks plausible",
                    xy=(32.0, bad), xytext=(-10, -16),
                    textcoords="offset points", ha="right", va="top",
                    color=_PAL["anti"], fontsize=6.3, family="monospace",
                    linespacing=1.7)
        x = self._cross
        ax.annotate("the collinear null is stronger than leading order:\n"
                    "the two wavefronts share their normal EXACTLY,\n"
                    "so amplitude gradients cannot tilt either k.\n"
                    "but it is NOT an absence of interference —\n"
                    f"ΔT⁰⁰/√(T_A⁰⁰T_B⁰⁰) = "
                    f"{x['collinear_interference']:.3f} collinear, "
                    f"{x['head_on_interference']:.3f} head-on.\n"
                    "T_A:T_B is smallest where the interference energy\n"
                    "is LARGEST, so it is the wrong proxy for it.",
                    xy=(0.03, 0.44), xycoords="axes fraction", va="top",
                    color=_PAL["dim"], fontsize=6.3, family="monospace",
                    linespacing=1.75)
        ax.set_xlabel("carrier ω₀", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("departure from the WKB value", color=_PAL["dim"],
                      fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper right", fontsize=6.4, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], which="both")
        ax.set_title("how the known limit is approached",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # -- the caustic ---------------------------------------------------------
    def _draw_caustic(self) -> None:
        ax = self.ax_ca
        scaled = np.linspace(0.05, 8.0, 400)
        for k, (c, col, mk) in enumerate(((8.5, _PAL["sym"], "o"),
                                          (16.5, _PAL["anti"], "s"),
                                          (32.5, _PAL["cool"], "^"))):
            ratios = np.array(
                [abs(green_omega(math.pi - u / c, complex(c, 0.0)))
                 * 4.0 * math.pi * math.sin(u / c) for u in scaled])
            # staggered markers: the three carriers land on the *same* curve,
            # so drawing them as lines would show one and hide two
            ax.plot(scaled[4 * k::14], ratios[4 * k::14], mk, ms=5.0,
                    color=col, alpha=0.95, zorder=5, label=f"ω = {c}")
        ax.plot(scaled, np.abs(np.sin(scaled)), lw=1.2,
                color=_PAL["text"], alpha=0.45, zorder=3, label="|sin(ωe)|")
        ax.axvline(1.0, color=_PAL["dim"], lw=0.8, ls=(0, (2, 4)))
        ax.annotate("e* ~ 1/ω", xy=(1.0, 1.02), color=_PAL["dim"],
                    fontsize=6.6, family="monospace")
        ax.annotate("exact / WKB depends on ωe ALONE — the three\n"
                    "carriers collapse to 6.6e-15.  WKB diverges at the\n"
                    "antipode; the exact amplitude saturates at ω/4π.",
                    xy=(0.035, 0.235), xycoords="axes fraction", va="top",
                    color=_PAL["dim"], fontsize=6.3, family="monospace",
                    linespacing=1.75)
        ax.set_xlabel("ω·e,   e = π − χ", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("exact amplitude / WKB", color=_PAL["dim"], fontsize=8)
        ax.set_ylim(0.0, 1.25)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper right", fontsize=6.4, framealpha=0.0,
                        labelcolor=_PAL["dim"], ncol=2)
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title("the caustic, and where geometric optics stops",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # -- frame ---------------------------------------------------------------
    def draw(self) -> None:
        self._draw_branches()
        self._draw_waveform()
        self._draw_limit()
        self._draw_caustic()
        m = self._multipath
        self.fig.suptitle("v59 — THE TWO-WAVE INVARIANT IS BRANCH-RESOLVED",
                          color=_PAL["text"], fontsize=13.2, y=0.962,
                          family="monospace")
        self.fig.text(0.5, 0.908,
                      "the WKB collinear/head-on result is the CONTROL, not "
                      "the result   ·   the exact multipath field departs from "
                      "it at O(1)",
                      color=_PAL["dim"], fontsize=8.4, ha="center",
                      family="monospace")
        a = self._audit
        self.fig.text(0.5, 0.070,
                      f"all four (i,j) channels carry just TWO predicted "
                      f"invariants, one per EXIT mouth "
                      f"({a['distinct_predictions'][0]:.6f}, "
                      f"{a['distinct_predictions'][1]:.6f}) — the field picks "
                      f"the right one to "
                      f"{100 * a['worst_relative_error']:.2f}%, not fitted",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.047,
                      f"the MOUTHS create the branch; the connection between "
                      f"them does not — N moves {a['beta_sweep_spread']:.1e} "
                      f"over β ∈ [0, 0.26], while its weight moves "
                      f"{100 * a['the_weight_moves_instead']:.1f}%.  "
                      f"W = −β still sees it",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.022,
                      f"conformally coupled scalar on a fixed ESU, throat "
                      f"strictly inside PR #257's cone (Löwner margin "
                      f"{m['loewner_margin']:.3f})   ·   NO backreaction: T is "
                      "computed from the field, never fed back",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")


def still(path: str) -> None:
    fig = TwoWaveFigure()
    fig.draw()
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v59.png")
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
