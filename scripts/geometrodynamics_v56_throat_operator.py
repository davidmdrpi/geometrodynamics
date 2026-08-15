#!/usr/bin/env python3
"""
Geometrodynamic QED — v56: a conserving throat cannot ring up
=============================================================

v55 solved a mouth relation self-consistently and said plainly what it was not:
a relation between field **values** carried by the free Green function, with no
normal-derivative matching, no reflected channel, a `1×1` mouth object where a
conserving junction needs `2×2` unitary, and `κ²` power throughput. It also
found its poles **off** the real axis and had to separate three thresholds to
say what that meant.

This round replaces the relation with the real object. A point-supported throat
is a **self-adjoint extension** of the Laplacian on `S³ ∖ {M⁺, M⁻}`, and von
Neumann's theorem parametrizes those by a unitary between the deficiency spaces
— `U(2)` — equivalently, by Krein's formula, a Hermitian `2×2` boundary matrix:

```
M(ω) q = φ_in|_mouths ,   M = A − Γ(ω) ,   Γ = [[g, G_d], [G_d, g]]

G(χ,ω) = sin(ω(π−χ)) / (4π sin χ sin(πω)) ,     g(ω) = −(ω/4π) cot(πω)
```

What the panels show
────────────────────
**Top left — the secular function.** `det M(ω)`, with the free spectrum as its
poles (dotted, `ω = n+1`) and the coupled spectrum as its zeros. In the
exchange-symmetric case it factorizes into `g ± G_d = α ± β`, both monotone from
`−∞` to `+∞` across every gap, so there are **exactly two** coupled frequencies
strictly between consecutive free ones. Interlacing, not merely shifting.

**Top right — the headline.** Newton from a grid of *complex* seeds, the same
method v55 used to find its poles. For a Hermitian `A` every converged root
lands on the real axis to `10⁻¹⁸`. For v55's directional relation, **every**
root is off it, several in the lower half plane where the mode grows — and that
is still true at `κ = 1`, so it is the **directionality**, not the loss.

**Bottom left — the channels v55 did not have.** `|r|²` and `|t|²` from the
Cayley transform `S = (A−ic)(A+ic)⁻¹`, summing to `1` at every reference scale.
v55's model sits off the unit circle at `(0, κ²)` — outside `U(2)` unless
`κ = 1`.

**Bottom right — switching the throat off.** The shift to the nearest free
frequency against `1/‖A‖`, slope `1`. Off is `‖A‖ → ∞`, not `A → 0`: the
diagonal of `A` is an *inverse* scattering length, so `α = 0` is a resonant
throat rather than no throat.

What is put in
──────────────
A linear scalar field on a fixed background, and the boundary matrix itself —
four real numbers chosen, not derived. `shells.junction` is what would fix them
from a matter model, and nothing here computes the exotic-matter bill. The
throat is **point-supported**: no interior, no proper length, and therefore no
delay — the `Δ` of v51–v55 is not a parameter of a point extension and does not
survive into one, which is a real loss of structure relative to those rounds.

**Not done:** no backreaction, no stress tensor, no topology change, no rate,
and no two-source invariant.

Usage
─────
    python scripts/geometrodynamics_v56_throat_operator.py --still out.png
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.waves.throat_operator import (
    DirectionalThroat,
    MouthPair,
    _krein_complex,
    complex_root_search,
    coupled_spectrum,
    spectrum_by_channel,
)

_PAL = {
    "bg": "#010106",
    "panel": "#02020a",
    "grid": "#0e1a2a",
    "rule": "#1a2838",
    "text": "#e8ecf4",
    "dim": "#6a8aad",
    "faint": "#1b2d42",
    "free": "#63798f",         # the uncoupled ESU spectrum
    "sym": "#7cff9e",          # symmetric channel
    "anti": "#ffb347",         # antisymmetric channel
    "hot": "#ff6b8a",          # the non-conserving control
    "cool": "#5cc8ff",
}

SEP = 1.3
ALPHA, BETA = 0.05, 0.03
N_GAPS = 6


# ════════════════════════════════════════════════════════════════════════════
class OperatorFigure:
    """The spectrum, the roots, the channels, and the decoupling limit."""

    def __init__(self, figsize=(13.8, 8.4)) -> None:
        self.pair = MouthPair(SEP, ALPHA, ALPHA, BETA)
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, width_ratios=[1.25, 1.0], height_ratios=[1.0, 0.9],
            left=0.058, right=0.975, top=0.845, bottom=0.135,
            wspace=0.22, hspace=0.40)
        self.ax_sec = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_cpx = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_uni = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_dec = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])

    # ── the secular function ────────────────────────────────────────────────
    def _draw_secular(self) -> None:
        ax = self.ax_sec
        a, b = ALPHA, BETA
        for m in range(1, N_GAPS + 1):
            t = np.linspace(0.0, 1.0, 900)[1:-1]
            ws = m + 0.5 * (1.0 - np.cos(math.pi * t))
            sym = np.array([self.pair.channel_functions(float(w))[0]
                            for w in ws])
            anti = np.array([self.pair.channel_functions(float(w))[1]
                             for w in ws])
            ax.plot(ws, sym, lw=1.4, color=_PAL["sym"], alpha=0.95, zorder=4,
                    label="g + G_d   (symmetric)" if m == 1 else None)
            ax.plot(ws, anti, lw=1.4, color=_PAL["anti"], alpha=0.95, zorder=4,
                    label="g − G_d   (antisymmetric)" if m == 1 else None)
        ax.axhline(a + b, color=_PAL["sym"], lw=0.9, ls=(0, (4, 3)),
                   alpha=0.7, zorder=3)
        ax.axhline(a - b, color=_PAL["anti"], lw=0.9, ls=(0, (4, 3)),
                   alpha=0.7, zorder=3)
        for n in range(1, N_GAPS + 2):
            ax.axvline(float(n), color=_PAL["free"], lw=0.9, ls=":", zorder=2)
        for r in coupled_spectrum(self.pair, N_GAPS):
            ax.plot([r["omega"]], [a + b if r["omega"] else 0], ".", ms=0)
        rows = spectrum_by_channel(self.pair, N_GAPS)["rows"]
        for r in rows:
            for k, col in (("symmetric", _PAL["sym"]),
                           ("antisymmetric", _PAL["anti"])):
                if r[k] is None:
                    continue
                ax.plot([r[k]], [a + b if k == "symmetric" else a - b], "o",
                        ms=4.6, color=col, mec=_PAL["bg"], mew=0.6, zorder=6)
        ax.annotate("dotted verticals: free spectrum ω = n+1 (the poles)   ·   "
                    "dots: coupled spectrum (the zeros)",
                    xy=(0.5, -0.20), xycoords="axes fraction",
                    color=_PAL["dim"], fontsize=6.6, ha="center",
                    family="monospace")
        ax.set_ylim(-0.35, 0.35)
        ax.set_xlim(1.0, N_GAPS + 1.0)
        ax.set_xlabel("ω", color=_PAL["dim"], fontsize=8.2)
        ax.set_ylabel("the two channel functions", color=_PAL["dim"],
                      fontsize=8.2)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper left", fontsize=6.8, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title("each channel is monotone across every gap — so exactly "
                     "two coupled frequencies per gap",
                     color=_PAL["text"], fontsize=9.0, pad=6)

    # ── the headline ────────────────────────────────────────────────────────
    def _draw_complex(self) -> None:
        ax = self.ax_cpx
        lo, hi = 1.1, N_GAPS + 0.9
        ax.axhspan(-0.8, 0.0, color=_PAL["hot"], alpha=0.06, zorder=0)
        ax.annotate("Im ω < 0 — the mode grows",
                    xy=(0.98, -0.70), xycoords=("axes fraction", "data"),
                    color=_PAL["hot"], fontsize=6.6, ha="right",
                    family="monospace")
        ax.axhline(0.0, color=_PAL["free"], lw=1.0, zorder=2)

        herm = MouthPair(SEP, 0.2, -0.13, 0.15 + 0.07j)
        hr = complex_root_search(
            lambda z: complex(np.linalg.det(_krein_complex(z, herm))),
            (lo, hi))
        for r in coupled_spectrum(herm, N_GAPS):
            ax.plot([r["omega"]], [0.0], "o", ms=6.0, color=_PAL["sym"],
                    mec=_PAL["bg"], mew=0.7, zorder=5)
        ax.plot([], [], "o", ms=6.0, color=_PAL["sym"], mec=_PAL["bg"],
                label="self-adjoint A — every root on the axis")

        for kap, mk, alp in ((0.3, "s", 0.95), (1.0, "^", 0.8)):
            ctl = DirectionalThroat(SEP, 1.0, +1, kap)
            got = complex_root_search(ctl.secular, (lo, hi))
            xs = [r.real for r in got["off_axis"]]
            ys = [r.imag for r in got["off_axis"]]
            ax.plot(xs, ys, mk, ms=5.4, color=_PAL["hot"], mec=_PAL["bg"],
                    mew=0.6, alpha=alp, zorder=6,
                    label=f"v55 directional, κ = {kap} — {got['n_off_axis']}"
                          f"/{got['n_roots']} off the axis")

        ax.set_xlim(lo, hi)
        ax.set_ylim(-0.8, 0.8)
        ax.set_xlabel("Re ω", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("Im ω", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper left", fontsize=6.4, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title(f"roots of det M in the complex plane — self-adjoint: "
                     f"{hr['n_off_axis']} off the axis",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # ── the channels ────────────────────────────────────────────────────────
    def _draw_unitary(self) -> None:
        ax = self.ax_uni
        th = np.linspace(0, math.pi / 2, 300)
        ax.plot(np.cos(th) ** 2, np.sin(th) ** 2, lw=1.0, ls=(0, (4, 3)),
                color=_PAL["free"], zorder=2, label="|r|² + |t|² = 1")
        scales = np.geomspace(0.01, 2.0, 60)
        for (a1, a2, b), col, lab, ms, z in (
                ((0.2, 0.2, 0.15), _PAL["sym"], "α = 0.2, β = 0.15", 6.4, 4),
                ((-0.4, 0.07, -0.09 + 0.31j), _PAL["cool"],
                 "α = (−0.4, 0.07), β = −0.09+0.31i", 2.4, 5)):
            p = MouthPair(SEP, a1, a2, b)
            rr, tt = [], []
            for c in scales:
                ch = p.channels(float(c))
                rr.append(abs(ch["reflection_1"]) ** 2)
                tt.append(abs(ch["transmission_21"]) ** 2)
            ax.plot(rr, tt, "o", ms=ms, lw=0, color=col, alpha=0.9,
                    zorder=z, label=lab)
        for kap, col in ((0.3, _PAL["hot"]), (0.6, _PAL["hot"]),
                         (1.0, _PAL["anti"])):
            ax.plot([0.0], [kap ** 2], "x", ms=7.0, mew=1.8, color=col,
                    zorder=6)
            ax.annotate(f"v55, κ={kap}", xy=(0.0, kap ** 2),
                        xytext=(6, -2), textcoords="offset points",
                        color=col, fontsize=6.2, family="monospace")
        ax.set_xlim(-0.04, 1.04)
        ax.set_ylim(-0.04, 1.09)
        ax.set_xlabel("|r|²  — reflection", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("|t|²  — transmission", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper right", fontsize=6.4, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title("the boundary operator is a unitary 2×2 — v55's model is "
                     "off the circle",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # ── switching it off ────────────────────────────────────────────────────
    def _draw_decoupling(self) -> None:
        ax = self.ax_dec
        ts = np.geomspace(1.0, 1e4, 9)
        shifts = []
        for t in ts:
            p = MouthPair(SEP, ALPHA * t, ALPHA * t, BETA * t)
            got = [w for r in spectrum_by_channel(p, 4)["rows"]
                   for w in (r["symmetric"], r["antisymmetric"])
                   if w is not None]
            shifts.append(max(abs(w - round(w)) for w in got))
        inv = 1.0 / ts
        ax.loglog(inv, shifts, "o-", lw=1.4, ms=4.6, color=_PAL["sym"],
                  alpha=0.95, zorder=4, label="worst shift from ω = n+1")
        ref = shifts[-1] * (inv / inv[-1])
        ax.loglog(inv, ref, lw=1.0, ls=(0, (4, 3)), color=_PAL["free"],
                  zorder=3, label="slope 1")
        sl = float(np.polyfit(np.log(inv[-4:]), np.log(shifts[-4:]), 1)[0])
        ax.annotate(f"asymptotic slope  {sl:.4f}", xy=(0.05, 0.88),
                    xycoords="axes fraction", color=_PAL["anti"],
                    fontsize=7.2, family="monospace")
        ax.annotate("off is ‖A‖ → ∞, not A → 0:\nthe diagonal is an INVERSE\n"
                    "scattering length",
                    xy=(0.05, 0.62), xycoords="axes fraction",
                    color=_PAL["dim"], fontsize=6.4, family="monospace",
                    linespacing=1.5)
        ax.set_xlabel("1 / ‖A‖", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("distance to the nearest free frequency",
                      color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="lower right", fontsize=6.6, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], which="both")
        ax.set_title("switch the throat off and the ESU spectrum comes back",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # ── frame ───────────────────────────────────────────────────────────────
    def draw(self) -> None:
        self._draw_secular()
        self._draw_complex()
        self._draw_unitary()
        self._draw_decoupling()
        self.fig.suptitle("v56 — A CONSERVING THROAT CANNOT RING UP",
                          color=_PAL["text"], fontsize=13.2, y=0.962,
                          family="monospace")
        self.fig.text(0.5, 0.908,
                      "a point-supported throat is a SELF-ADJOINT EXTENSION: "
                      "M(ω) = A − Γ(ω) with A Hermitian   ·   flux "
                      "conservation IS A = A†",
                      color=_PAL["dim"], fontsize=8.4, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.062,
                      "v55's off-axis poles were its own directionality — it "
                      "is unstable even at κ = 1, where nothing is lost",
                      color=_PAL["dim"], fontsize=7.8, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.030,
                      "put in: a linear field on a fixed background and the "
                      "boundary matrix itself — four numbers chosen, not "
                      "derived   ·   the throat is POINT-supported, so it has "
                      "no interior and no delay",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")


# ════════════════════════════════════════════════════════════════════════════
def still(path: str) -> None:
    fig = OperatorFigure()
    fig.draw()
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v56.png")
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
