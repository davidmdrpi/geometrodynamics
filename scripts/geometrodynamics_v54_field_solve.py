#!/usr/bin/env python3
"""
Geometrodynamic QED — v54: the field reproduces the ledger, and signs it
=======================================================================

v53 built a **ray** ledger — short way, long way, winding — and ended by
conceding that rank counting had reached its limit. This round solves the field
instead and asks whether the ledger survives contact with it.

It does, and the branches are **exact support** rather than stationary-phase
contributions. On the Einstein static universe the scalar Laplacian has
eigenvalues `−n(n+2)` and `R = 6`, so the **conformally** coupled massless field
has `ω² = n(n+2) + 1 = (n+1)²` — integer frequencies. The retarded Green
function is exactly periodic, and in closed form it is a sum of images:

```
G(χ,t) = 1/(4π sin χ) [ Σ_k δ(t − χ − 2πk) − Σ_k δ(t + χ − 2πk) ]
```

What the panels show
────────────────────
**Left — the solved field over `(χ, t)`.** Colour is *signed*: amber positive,
cyan negative. The dashed curves are v53's branch times, drawn from the ray
ledger and never fitted. Two things are visible at once — the field lives only
on those curves, and **the sign flips at every focal crossing**. That striping
is the Maslov index, and it is the first quantity in this arc that the ray
picture could not in principle have carried.

**Top right — one time slice.** The peaks of the *solved* field against the ray
ledger's arrival times, with the predicted signs.

**Bottom right — which field the ledger belongs to.** The **minimally** coupled
scalar has `ω = √(n(n+2))`, irrational, so no images and no sharp branches: 63%
of the peak amplitude sits *between* the arrivals against `4e-08` for the
conformal field. v53 never said which field its ledger described, because rays
cannot tell the two apart.

What is put in
──────────────
A linear scalar field on a fixed background, with the throat still an
**identification map** rather than a solution — `shells.junction` priced it and
the bill is inherited, unpaid. **Not done:** no backreaction, no topology
change, no rate, and no two-source invariant yet.

Usage
─────
    python scripts/geometrodynamics_v54_field_solve.py --still out.png
"""

from __future__ import annotations

import argparse
import math
from typing import List, Optional, Tuple

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

from geometrodynamics.waves.field_solve import (
    TWO_PI,
    branch_arrivals,
    field_peaks,
    spectral_field,
)

_PAL = {
    "bg": "#010106",
    "panel": "#02020a",
    "grid": "#0e1a2a",
    "rule": "#1a2838",
    "text": "#e8ecf4",
    "dim": "#6a8aad",
    "faint": "#1b2d42",
    "pos": "#ffb347",          # positive field
    "neg": "#5cc8ff",          # negative field — a sign flip, not a colour whim
    "branch": "#7cff9e",
    "minimal": "#ff6b8a",
}

# A diverging map in the series' own two hues, so the sign reads as a sign.
_CMAP = LinearSegmentedColormap.from_list(
    "signed", [_PAL["neg"], "#0a1420", "#02020a", "#160f04", _PAL["pos"]])


# ════════════════════════════════════════════════════════════════════════════
def _field_grid(chis: np.ndarray, ts: np.ndarray, width: float,
                conformal: bool = True) -> np.ndarray:
    """Vectorised mode sum over the whole ``(χ, t)`` plane.

    The Gaussian source width sets the truncation: ``e^{−(mw)²/2}`` is already
    below ``1e-13`` by ``m = 8/w``, so the sum is complete rather than cut off.
    """
    n_max = max(64, int(8.0 / width))
    m = np.arange(1, n_max + 1)
    om = m.astype(float) if conformal else np.sqrt((m - 1.0) * (m + 1.0))
    damp = np.exp(-(om * width) ** 2 / 2.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        amp = np.where(om > 1e-12, damp / np.where(om > 1e-12, om, 1.0), 0.0)
    # zonal weight m·sin(mχ) / (2π² sinχ), then Σ_m amp·sin(ωt)
    sc = (m[None, :] * np.sin(np.outer(chis, m))
          / (2.0 * math.pi ** 2 * np.sin(chis))[:, None])
    st = np.sin(np.outer(ts, om))
    return (sc * amp[None, :]) @ st.T


# ════════════════════════════════════════════════════════════════════════════
class FieldFigure:
    """The solved field, its branch curves, and the coupling that owns them."""

    def __init__(self, width: float = 0.045, chi_slice: float = 1.1,
                 t_max: float = 4.0 * math.pi, figsize=(13.6, 8.2)) -> None:
        self.width = width
        self.chi_slice = chi_slice
        self.t_max = t_max
        self.chis = np.linspace(0.12, math.pi - 0.12, 420)
        self.ts = np.linspace(0.02, t_max, 900)
        self.grid = _field_grid(self.chis, self.ts, width)

        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, width_ratios=[1.25, 1.0], height_ratios=[1.0, 0.85],
            left=0.055, right=0.975, top=0.845, bottom=0.135,
            wspace=0.20, hspace=0.34)
        self.ax_map = self.fig.add_subplot(gs[:, 0], facecolor=_PAL["panel"])
        self.ax_slice = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_cmp = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])

    # ── the (χ, t) map ──────────────────────────────────────────────────────
    def _draw_map(self) -> None:
        ax = self.ax_map
        v = float(np.percentile(np.abs(self.grid), 99.6))
        ax.imshow(self.grid.T, origin="lower", aspect="auto",
                  extent=[self.chis[0], self.chis[-1], self.ts[0],
                          self.ts[-1]],
                  cmap=_CMAP, norm=TwoSlopeNorm(vmin=-v, vcenter=0.0, vmax=v),
                  interpolation="nearest")

        # v53's ray ledger, drawn on top and never fitted
        cc = np.linspace(self.chis[0], self.chis[-1], 200)
        k = 0
        while True:
            drew = False
            for curve, lab in ((cc + TWO_PI * k, 0),
                               (TWO_PI * (k + 1) - cc, 1)):
                if curve.min() <= self.t_max:
                    ax.plot(cc, curve, ls=(0, (4, 3)), lw=0.9,
                            color=_PAL["branch"], alpha=0.55, zorder=3)
                    drew = True
            if not drew:
                break
            k += 1
            if k > 6:
                break
        ax.axvline(self.chi_slice, color=_PAL["text"], lw=0.8, alpha=0.45,
                   zorder=4)
        for j in range(1, 5):
            ax.axhline(math.pi * j, color=_PAL["faint"], lw=0.6, ls=":",
                       alpha=0.8, zorder=2)
            ax.annotate(f"focus {j}", xy=(self.chis[0], math.pi * j),
                        xytext=(4, 3), textcoords="offset points",
                        color=_PAL["faint"], fontsize=6.0, ha="left",
                        family="monospace", zorder=5)
        ax.set_xlabel("χ — geodesic distance from the source",
                      color=_PAL["dim"], fontsize=8.5)
        ax.set_ylabel("t", color=_PAL["dim"], fontsize=8.5)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        ax.set_title("the solved field — colour is the SIGN, dashes are v53's "
                     "branch times",
                     color=_PAL["text"], fontsize=9.4, pad=6)

    # ── one slice ───────────────────────────────────────────────────────────
    def _draw_slice(self) -> None:
        ax = self.ax_slice
        ts = np.linspace(0.02, self.t_max, 3000)
        f = _field_grid(np.array([self.chi_slice]), ts, self.width)[0]
        ax.axhline(0.0, color=_PAL["faint"], lw=0.7)
        ax.plot(ts, f, lw=1.3, color=_PAL["text"], alpha=0.9, zorder=4)
        ax.fill_between(ts, 0, f, where=f > 0, color=_PAL["pos"], alpha=0.35,
                        zorder=3)
        ax.fill_between(ts, 0, f, where=f < 0, color=_PAL["neg"], alpha=0.35,
                        zorder=3)
        for r in branch_arrivals(self.chi_slice, self.t_max):
            ax.axvline(r["t"], color=_PAL["branch"], lw=0.8, ls=(0, (3, 3)),
                       alpha=0.7, zorder=2)
            ax.annotate(f"{'+' if r['sign'] > 0 else '−'}",
                        xy=(r["t"], 0.86), xycoords=("data", "axes fraction"),
                        color=_PAL["pos"] if r["sign"] > 0 else _PAL["neg"],
                        fontsize=10, ha="center", family="monospace")
        ax.set_xlim(0, self.t_max)
        ax.set_xlabel("t", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel(f"φ at χ = {self.chi_slice}", color=_PAL["dim"],
                      fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title("peaks of the SOLVED field against the ray ledger — "
                     "and the sign it predicts",
                     color=_PAL["text"], fontsize=8.3, pad=6)

    # ── which coupling owns the ledger ──────────────────────────────────────
    def _draw_cmp(self) -> None:
        ax = self.ax_cmp
        ts = np.linspace(0.02, TWO_PI, 2200)
        fc = _field_grid(np.array([self.chi_slice]), ts, 0.06, True)[0]
        fm = _field_grid(np.array([self.chi_slice]), ts, 0.06, False)[0]
        ax.axhline(0.0, color=_PAL["faint"], lw=0.7)
        ax.plot(ts, fc, lw=1.4, color=_PAL["branch"], alpha=0.95,
                label="conformal  ω = n+1", zorder=4)
        ax.plot(ts, fm, lw=1.2, color=_PAL["minimal"], alpha=0.9,
                label="minimal  ω = √(n(n+2))", zorder=3)
        for r in branch_arrivals(self.chi_slice, TWO_PI):
            ax.axvline(r["t"], color=_PAL["faint"], lw=0.8, ls=(0, (3, 3)),
                       zorder=1)
        ax.set_xlim(0, TWO_PI)
        ax.set_xlabel("t", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("φ", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper right", fontsize=6.8, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title("only the conformal field has branches — the minimal one "
                     "fills the gaps",
                     color=_PAL["text"], fontsize=8.3, pad=6)

    # ── frame ───────────────────────────────────────────────────────────────
    def draw(self) -> None:
        self._draw_map()
        self._draw_slice()
        self._draw_cmp()
        self.fig.suptitle("v54 — THE FIELD REPRODUCES THE LEDGER, AND SIGNS IT",
                          color=_PAL["text"], fontsize=13.2, y=0.962,
                          family="monospace")
        self.fig.text(0.5, 0.908,
                      "ω = n+1 on the Einstein static universe, so the Green "
                      "function is a sum of images   ·   the branches are EXACT "
                      "support, not stationary phase",
                      color=_PAL["dim"], fontsize=8.4, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.062,
                      "the sign is (−1) per focal crossing — the Maslov index, "
                      "and the first quantity in this arc a path-length ledger "
                      "could not have carried",
                      color=_PAL["dim"], fontsize=7.8, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.030,
                      "put in: a linear field on a fixed background, the throat "
                      "still an identification map   ·   not done: backreaction, "
                      "topology change, and the two-source invariant",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")


# ════════════════════════════════════════════════════════════════════════════
def still(path: str) -> None:
    fig = FieldFigure()
    fig.draw()
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v54.png")
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
