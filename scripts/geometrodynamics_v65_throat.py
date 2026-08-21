#!/usr/bin/env python3
"""
Geometrodynamic QED - v65: which throat is physical, and the sign reverses
=========================================================================

PR #264 measured a signed dA/A on a product tube of area 4 pi and length 0.9 --
both free parameters -- and found the mouths CLOSE. Then, asked to match the
tube's area to its own mouths', it found the sign REVERSES. So:

    A and L were free parameters. WHICH VALUES ARE PHYSICAL?

They were never free. On a maximal slice the background constraint is
R_hat = 16 pi G rho, so a profile does not choose its matter -- the matter is
whatever the profile implies. Ask instead for a throat needing NO matter
(R_hat = 0) that glues on with NO surface layer, and both conditions are used
up: f0 = sin^3 a is forced, and L and I follow in closed form. Nothing is left
to choose.

On that throat the interference OPENS both mouths.

What the panels show
--------------------
**Top left - the profile, and what each throat holds.** The forced vacuum
profile f'^2 = 1 - f0/f against the round sphere it glues to, with the neck
marked. That profile IS the time-symmetric Schwarzschild slice, so f0 = 2M and
the mass is derived: M = sin^3(a)/2. The bar inset is the matter each candidate tube would have to contain:
rho/3 for #261-#264's, 133 rho for the matched one, and ZERO for this one.

**Top right - there is no cavity.** The l=0 tube operator is nabla^2 + R_hat/2.
With matter in the tube that is a wave equation and the throat is a resonant
cavity, with poles at kL = n pi -- the structure that made #264's sign
conditional. With R_hat = 0 it is the plain Laplacian, (f^2 u')' = 0, and the
solutions are monotone. Nothing to flip a sign across.

**Bottom left - the mechanism, and it is one number.** Split any symmetric
two-port as Y = G [[-1,1],[1,-1]] + shunt I. The conductance is scanned over
eight orders and never changes the sign. The shunt passes through a pole near
2e-03 and flips it. A vacuum tube has shunt ZERO BY IDENTITY: (f^2 u')' = 0
leaves nowhere to put monopole flux.

**Bottom right - the answer, and what it costs.** dA/A on the physical throat
against #264's, at both mouth radii. The sign reverses. The magnitude is 3000x
larger and grows as a^-3 -- not noise, but the physics of a throat that barely
lifts the constraint's exact degeneracy.

What is put in
--------------
The source is #263's, on a FIXED background with POINT sources; the ESU's fluid
is held RIGID; the exterior is the round S^3 with two balls removed, so the
ambient's own response to carrying a handle is not modelled. A vacuum tube is
the SIMPLEST acceptable matter, not the only one -- any rho(s) >= 0 that glues
smoothly is admissible.

Usage
-----
    python scripts/geometrodynamics_v65_throat.py --still out.png
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.waves.areal import (INTERFERENCE_MOMENTS, MOUTHS,
                                          TubeModel, basis_channels,
                                          solve_matching)
from geometrodynamics.waves.physical_throat import (VacuumThroat,
                                                    product_tube_density_ratio)

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


class _Fixed:
    def __init__(self, y0, y1):
        self._y0, self._y1 = y0, y1

    def admittance(self, ell):
        return self._y0 if int(ell) == 0 else self._y1


class ThroatFigure:
    """The profile, the missing cavity, the mechanism, and the reversal."""

    def __init__(self, figsize=(13.8, 8.6)) -> None:
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, left=0.075, right=0.975, top=0.820, bottom=0.135,
            wspace=0.26, hspace=0.52)
        self.ax_prof = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_cav = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_mech = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_ans = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])
        self.m = INTERFERENCE_MOMENTS[1]
        self.a = self.m.radius
        self.basis = basis_channels(MOUTHS, self.a)
        self.throat = VacuumThroat(mouth_radius=self.a)

    # ── panels ──────────────────────────────────────────────────────────────
    def _profile(self) -> None:
        ax = self.ax_prof
        a = self.a
        f0, mouth = self.throat.neck_radius(), math.sin(a)
        f = np.linspace(f0, mouth, 900)
        s = (f0 * np.arccosh(np.sqrt(f / f0))
             + np.sqrt(f * (f - f0)))                 # proper length from neck
        ax.plot(np.concatenate([-s[::-1], s]),
                np.concatenate([f[::-1], f]) / mouth,
                lw=2.0, color=_PAL["open"], label="vacuum throat (forced)")
        half = self.throat.length() / 2.0
        for side in (-1, 1):
            chi = np.linspace(a, a + 0.06, 200)
            ax.plot(side * (half + (chi - a)), np.sin(chi) / mouth,
                    lw=1.6, ls="--", color=_PAL["cool"],
                    label="ambient S^3" if side == 1 else None)
        ax.plot([0.0], [f0 / mouth], "o", ms=7, color=_PAL["x"], zorder=5)
        ax.annotate(f"neck   f0 = sin^3 a = {f0:.2e}\n"
                    f"       = 2M,  so M = {f0/2:.2e}", xy=(0.0, f0 / mouth),
                    xytext=(0.016, 0.40), color=_PAL["x"], fontsize=6.8,
                    family="monospace",
                    arrowprops=dict(arrowstyle="-", color=_PAL["x"], lw=0.7))
        ax.axhline(1.0, color=_PAL["dim"], lw=0.8, ls=":")
        _style(ax, "the profile the constraint forces",
               "proper length from the neck", "f(s) / sin a")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=7.0, loc="center right",
                  framealpha=0.92)
        ax.text(0.5, 0.93,
                f"L = {self.throat.length():.5f} = {self.throat.length()/a:.2f} a"
                "   -- forced, not chosen",
                transform=ax.transAxes, ha="center", color=_PAL["text"],
                fontsize=7.2, family="monospace")
        inset = ax.inset_axes([0.09, 0.50, 0.36, 0.34])
        names = ["#264\n4pi", "matched\n4pi sin^2a", "vacuum\n(this)"]
        vals = [product_tube_density_ratio(4 * math.pi),
                product_tube_density_ratio(4 * math.pi * math.sin(a) ** 2),
                0.0]
        inset.bar(range(3), np.maximum(vals, 1e-3),
                  color=[_PAL["close"], _PAL["close"], _PAL["open"]],
                  edgecolor=_PAL["rule"], lw=0.5)
        inset.set_yscale("log")
        inset.axhline(1.0, color=_PAL["dim"], lw=0.9, ls="--")
        inset.set_xticks(range(3))
        inset.set_xticklabels(names, fontsize=5.2, family="monospace",
                              color=_PAL["dim"])
        inset.set_title("matter the tube must hold, / ambient",
                        color=_PAL["dim"], fontsize=6.0, family="monospace")
        inset.tick_params(colors=_PAL["dim"], labelsize=5.4)
        inset.set_facecolor(_PAL["panel"])
        for sp in inset.spines.values():
            sp.set_color(_PAL["rule"])

    def _cavity(self) -> None:
        ax = self.ax_cav
        s = np.linspace(0.0, 1.0, 400)
        for kl, col, lab in ((0.9, _PAL["dim"], "matter tube, kL = 0.9"),
                             (2.6, _PAL["close"], "matter tube, kL = 2.6"),
                             (5.9, _PAL["kern"], "matter tube, kL = 5.9")):
            ax.plot(s, np.cos(kl * s), lw=1.3, color=col, label=lab)
        # the real vacuum mode: u = A + B int ds/f^2, built on the profile
        f0, mouth = self.throat.neck_radius(), math.sin(self.a)
        t = np.linspace(0.0, math.sqrt(mouth - f0), 4000)
        ff = f0 + t ** 2
        ds = 2.0 * np.sqrt(ff)                      # ds/dt on the vacuum profile
        arc = np.concatenate([[0.0], np.cumsum(0.5 * (ds[1:] + ds[:-1])
                                               * np.diff(t))])
        pot = np.concatenate([[0.0], np.cumsum(0.5 * (ds[1:] / ff[1:] ** 2
                                                      + ds[:-1] / ff[:-1] ** 2)
                                               * np.diff(t))])
        half = np.concatenate([-arc[::-1], arc[1:]]) / (2.0 * arc[-1]) + 0.5
        mode = np.concatenate([-pot[::-1], pot[1:]]) / pot[-1]
        ax.plot(half, mode, lw=2.4, color=_PAL["open"],
                label="vacuum throat: A + B int ds/f^2")
        ax.axhline(0.0, color=_PAL["dim"], lw=0.8)
        ax.set_ylim(-1.35, 2.15)
        _style(ax, "there is no cavity", "position along the tube  s / L",
               "the l=0 tube mode")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.6, loc="lower right",
                  framealpha=0.92)
        ax.text(0.5, 0.975,
                "nabla^2 + R/2 :  with matter a WAVE equation (poles at "
                "kL = n pi)\nwith R = 0 the plain Laplacian, (f^2 u')' = 0 -- "
                "MONOTONE\n(its whole rise sits at the neck, where f is "
                "smallest)",
                transform=ax.transAxes, ha="center", va="top",
                color=_PAL["text"], fontsize=6.6, family="monospace")

    def _mechanism(self) -> None:
        ax = self.ax_mech
        j = np.array([[-1.0, 1.0], [1.0, -1.0]])
        y1 = self.throat.admittance(1)
        g = self.throat.conductance()

        def run(cond, shunt):
            return float(np.asarray(solve_matching(
                MOUTHS, self.a, _Fixed(cond * j + shunt * np.eye(2), y1),
                self.m.as_source(), self.m.signed_obstruction(),
                basis=self.basis)["areal_response"], float)[0])

        mults = np.geomspace(1e-3, 1e5, 40)
        ax.plot(mults, [run(g * x, 0.0) for x in mults], lw=1.8,
                color=_PAL["open"], label="scan the CONDUCTANCE (shunt = 0)")
        shunts = np.geomspace(1e-6, 3e1, 60)
        ax.plot(shunts, [run(g, x) for x in shunts], lw=1.8,
                color=_PAL["close"], label="scan the SHUNT (conductance fixed)")
        ax.axhline(0.0, color=_PAL["dim"], lw=1.0)
        ax.axvline(float(TubeModel().shunt()), color=_PAL["kern"], lw=1.2,
                   ls=":")
        ax.text(float(TubeModel().shunt()) * 1.15, -40.0,
                "#264's\nshunt", color=_PAL["kern"], fontsize=6.4,
                family="monospace")
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=1e-2, linscale=0.5)
        ax.set_ylim(-3e2, 3e2)
        _style(ax, "the mechanism, and it is one number",
               "multiplier (conductance)  /  shunt", "dA/A mouth 1  [symlog]")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.6, loc="lower left")
        ax.text(0.97, 0.93,
                "a vacuum tube has shunt ZERO\nby identity: (f^2 u')' = 0",
                transform=ax.transAxes, ha="right", va="top",
                color=_PAL["text"], fontsize=6.8, family="monospace",
                bbox=dict(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                          alpha=0.92, pad=3.0))

    def _answer(self) -> None:
        ax = self.ax_ans
        wide = TubeModel()
        rows = []
        for mom in INTERFERENCE_MOMENTS:
            if mom.points != max(x.points for x in INTERFERENCE_MOMENTS
                                 if x.radius == mom.radius):
                continue
            basis = basis_channels(MOUTHS, mom.radius)
            vac = VacuumThroat(mouth_radius=mom.radius)
            for lab, throat in (("#264 product tube", wide),
                                ("vacuum throat", vac)):
                v = np.asarray(solve_matching(
                    MOUTHS, mom.radius, throat, mom.as_source(),
                    mom.signed_obstruction(), basis=basis)["areal_response"],
                    float)
                rows.append((mom.radius, lab, v))
        xs = np.arange(len(rows))
        cols = [_PAL["close"] if "264" in r[1] else _PAL["open"] for r in rows]
        for k, off in ((0, -0.17), (1, 0.17)):
            ax.bar(xs + off, [r[2][k] for r in rows], width=0.32,
                   color=cols, alpha=0.95 if k == 0 else 0.55,
                   edgecolor=_PAL["rule"], lw=0.6)
        ax.axhline(0.0, color=_PAL["dim"], lw=1.0)
        ax.set_yscale("symlog", linthresh=1e-3, linscale=0.5)
        ax.set_xticks(xs)
        ax.set_xticklabels([f"a={r[0]:.2f}\n{r[1]}" for r in rows],
                           fontsize=6.0, family="monospace", color=_PAL["dim"])
        _style(ax, "the answer, and what it costs", "",
               "dA/A  [symlog, units 2 pi G]")
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi * 12.0)
        ax.text(0.5, 0.95,
                "the sign REVERSES on the throat that is forced",
                transform=ax.transAxes, ha="center", va="top",
                color=_PAL["open"], fontsize=8.0, family="monospace")
        ax.text(0.5, 0.055,
                "3000x larger and growing as a^-3: a zero-shunt throat barely\n"
                "lifts the constraint's exact degeneracy.  the sign is robust;\n"
                "the amplitude at which linearising holds is the open question.",
                transform=ax.transAxes, ha="center", color=_PAL["dim"],
                fontsize=6.4, family="monospace")

    # ── assembly ────────────────────────────────────────────────────────────
    def draw(self) -> None:
        self._profile()
        self._cavity()
        self._mechanism()
        self._answer()
        self.fig.text(0.5, 0.955,
                      "v65  -  which throat is physical, and the sign reverses",
                      color=_PAL["text"], fontsize=15.5, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.915,
                      "R_hat = 16 pi G rho on a maximal slice, so a profile "
                      "does not choose its matter -- the matter is whatever "
                      "the profile implies",
                      color=_PAL["dim"], fontsize=8.4, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.884,
                      "R_hat = 0 + K = 0 + a spherical neck IS the "
                      "time-symmetric Schwarzschild slice, so f0 = 2M and "
                      "M = sin^3(a)/2 -- the mass from the mouth, nothing "
                      "left to choose",
                      color=_PAL["open"], fontsize=7.4, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.856,
                      "quasi-local (no asymptotic region, so no ADM mass), "
                      "dimensionless M/R, a handle in ONE S^3 -- and the neck "
                      "is a minimal surface, hence an APPARENT HORIZON",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.022,
                      "source is #263's, FIXED background, POINT sources   -   "
                      "the ESU's fluid held RIGID   -   a vacuum tube is the "
                      "SIMPLEST acceptable matter, not the only one",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")


def still(path: str) -> None:
    fig = ThroatFigure()
    fig.draw()
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v65.png")
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
