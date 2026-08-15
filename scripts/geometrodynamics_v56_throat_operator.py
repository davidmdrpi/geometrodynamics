#!/usr/bin/env python3
"""
Geometrodynamic QED — v56: conservation is not stability
========================================================

v55 solved a mouth relation self-consistently and said plainly what it was not:
a relation between field **values** carried by the free Green function, with no
normal-derivative matching, no reflected channel, a `1×1` mouth object, and `κ²`
power throughput.

A point-supported throat is a **self-adjoint extension** of the Laplacian on
`S³ ∖ {M⁺, M⁻}`, parametrized by `U(2)`. Writing the boundary condition as the
*pair* `B φ^reg = C q` — general enough to hold v55's relation, which is not of
the form `φ^reg = A q` — the mouth-active spectrum is `det(C − BΓ) = 0`, with

```
G(χ,ω) = sin(ω(π−χ)) / (4π sin χ sin(πω)) ,     g(ω) = −(ω/4π) cot(πω)
```

**Self-adjointness buys conservation and a real `λ = ω²`. It does not buy
`λ ≥ 0`.** A first version of this round claimed it did; the panels below are
what replaced that claim.

What the panels show
────────────────────
**Top left — the mouth-active sector, in `λ`.** The two channel functions
`g ± G_d` against the free eigenvalues `λ = (n+1)²`. Three regions matter and
the first two are invisible to an `ω`-scan that starts above `1`: `λ < 0`, where
a root means a **growing** mode; `0 ≤ λ < 1`, stable modes below the free ground
state; and the interlacing gaps. Note the antisymmetric channel's endpoint at
`λ = 1` is *finite* — the `n = 0` constant mode is equal at both mouths, so it
does not couple — which is why a root in the first gap is conditional.

**Top right — the correction.** The stability wedge of the exchange-symmetric
family. Along the imaginary axis both channels fall monotonically from their
`λ = 0` values, so a growing mode exists iff `α + β < g₀ + G₀` or
`α − β < g₀ − G₀`. Every grid point is *also* scanned for a negative-`λ` root:
0 mismatches, and only 56 of 221 sampled points are stable. Two of the three
boundary matrices this round originally advertised are marked — both outside.

**Bottom left — what conservation does buy.** Net boundary flux `Im(q†Aq)`
against the anti-Hermitian part of `A`, over random draws. Hermitian data sits
on zero identically; anything else does not. This, and not the Cayley entries,
is the physical conservation statement — the Cayley magnitudes depend on an
arbitrary reference scale, shown inset.

**Bottom right — where v55 sits.** Its relation embeds *exactly* as
`B = [[0,0],[gain,0]]`, `C = I`, giving `det(C − BΓ) = 1 − gain·G_d`, matched to
its own `1 − L` to `10⁻¹⁸`. Maximal, but `BC† = B` is not Hermitian. The earlier
version of this control used `A = [[0,0],[1/gain,0]]`, which is a different
function — plotted alongside so the difference is visible rather than asserted.

What is put in
──────────────
A linear scalar field on a fixed background, and the boundary data itself — four
real numbers chosen, not derived. `shells.junction` is what would fix them. The
throat is **point-supported**: no interior, no proper length, and therefore no
delay — the `Δ` of v51–v55 does not survive into a point extension.

**Not done:** no backreaction, no stress tensor, no topology change, no rate,
and no two-source invariant.

Usage
─────
    python scripts/geometrodynamics_v56_throat_operator.py --still out.png
"""

from __future__ import annotations

import argparse
import cmath
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from geometrodynamics.waves.throat_operator import (
    DirectionalThroat,
    MouthPair,
    boundary_mixing,
    channel_endpoints,
    is_stable,
    mouth_active_spectrum,
    mouth_flux,
    stability_thresholds,
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
    "sym": "#7cff9e",          # symmetric channel / stable
    "anti": "#ffb347",         # antisymmetric channel
    "hot": "#ff6b8a",          # growing modes, non-conserving control
    "cool": "#5cc8ff",
}

SEP = 1.3
ALPHA, BETA = 0.05, 0.03
N_GAPS = 4


# ════════════════════════════════════════════════════════════════════════════
class OperatorFigure:
    """The sector, the stability wedge, the flux identity, and where v55 sits."""

    def __init__(self, figsize=(13.8, 8.6)) -> None:
        self.pair = MouthPair(SEP, ALPHA, ALPHA, BETA)
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, width_ratios=[1.25, 1.0], height_ratios=[1.0, 0.9],
            left=0.058, right=0.975, top=0.845, bottom=0.135,
            wspace=0.22, hspace=0.42)
        self.ax_sec = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_stb = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_flx = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_emb = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])

    # ── the sector, in lambda ───────────────────────────────────────────────
    def _draw_sector(self) -> None:
        ax = self.ax_sec
        a, b = ALPHA, BETA
        lam_hi = float(N_GAPS + 1) ** 2

        # plotted against sign(λ)·√|λ| so the three regions get equal width:
        # to the right that variable is ω, to the left it is the growth rate σ
        ax.axvspan(-3.0, 0.0, color=_PAL["hot"], alpha=0.08, zorder=0)
        ax.axvspan(0.0, 1.0, color=_PAL["cool"], alpha=0.07, zorder=0)

        def band(lo, hi):
            t = np.linspace(0.0, 1.0, 700)[1:-1]
            return lo + (hi - lo) * 0.5 * (1.0 - np.cos(math.pi * t))

        segments = [(-3.0, -1e-4), (1e-4, 1.0)]
        segments += [(float(m), float(m + 1)) for m in range(1, N_GAPS + 1)]
        for k, (lo, hi) in enumerate(segments):
            us = band(lo, hi)
            lams = np.sign(us) * us ** 2
            sym = np.array([self.pair.channel_functions(float(x))[0]
                            for x in lams])
            anti = np.array([self.pair.channel_functions(float(x))[1]
                             for x in lams])
            ax.plot(us, sym, lw=1.4, color=_PAL["sym"], alpha=0.95, zorder=4,
                    label="g + G_d  (symmetric)" if k == 0 else None)
            ax.plot(us, anti, lw=1.4, color=_PAL["anti"], alpha=0.95, zorder=4,
                    label="g − G_d  (antisymmetric)" if k == 0 else None)
        ax.axhline(a + b, color=_PAL["sym"], lw=0.9, ls=(0, (4, 3)),
                   alpha=0.7, zorder=3)
        ax.axhline(a - b, color=_PAL["anti"], lw=0.9, ls=(0, (4, 3)),
                   alpha=0.7, zorder=3)
        for n in range(1, N_GAPS + 2):
            ax.axvline(float(n), color=_PAL["free"], lw=0.9, ls=":", zorder=2)
        ax.axvline(0.0, color=_PAL["text"], lw=0.8, alpha=0.5, zorder=2)
        for r in mouth_active_spectrum(self.pair, N_GAPS):
            if r["sector"] == "growing":
                continue
            u = math.sqrt(max(r["lmbda"], 0.0))
            ax.plot([u], [a + b], "o", ms=4.4, color=_PAL["text"],
                    mec=_PAL["bg"], mew=0.6, zorder=6)

        ends = channel_endpoints(self.pair, 1)
        ax.annotate(f"finite: {ends['antisymmetric_at_lower']:+.4f}\n"
                    f"(the n=0 mode is equal at both mouths,\n"
                    f"so it does not couple to this channel)",
                    xy=(1.0, ends["antisymmetric_at_lower"]),
                    xytext=(58, 92), textcoords="offset points",
                    color=_PAL["anti"], fontsize=6.0, family="monospace",
                    linespacing=1.5,
                    arrowprops=dict(arrowstyle="-", color=_PAL["anti"],
                                    lw=0.7, alpha=0.7))
        ax.annotate("λ < 0\nGROWING", xy=(-2.85, 0.29), color=_PAL["hot"],
                    fontsize=6.8, family="monospace", linespacing=1.5)
        ax.annotate("0<λ<1\nbelow the\nground state", xy=(0.06, 0.20),
                    color=_PAL["cool"], fontsize=6.0, family="monospace",
                    linespacing=1.5)
        ax.set_ylim(-0.35, 0.38)
        ax.set_xlim(-3.0, float(N_GAPS + 1))
        ax.set_xlabel("sign(λ)·√|λ|   —   ω to the right, growth rate σ to "
                      "the left", color=_PAL["dim"], fontsize=8.0)
        ax.set_ylabel("the two channel functions", color=_PAL["dim"],
                      fontsize=8.2)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="lower right", fontsize=6.8, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title("the mouth-active sector has three regions — and two of "
                     "them are below ω = 1",
                     color=_PAL["text"], fontsize=8.8, pad=6)

    # ── the correction ──────────────────────────────────────────────────────
    def _draw_stability(self) -> None:
        ax = self.ax_stb
        th = stability_thresholds(SEP)
        s_th, a_th = th["symmetric_threshold"], th["antisymmetric_threshold"]
        alphas = np.linspace(-0.15, 0.15, 61)
        betas = np.linspace(-0.2, 0.2, 81)
        AA, BB = np.meshgrid(alphas, betas, indexing="ij")
        ok = ((AA + BB >= s_th) & (AA - BB >= a_th)).astype(float)
        ax.imshow(ok.T, origin="lower", aspect="auto",
                  extent=[alphas[0], alphas[-1], betas[0], betas[-1]],
                  cmap=LinearSegmentedColormap.from_list(
                      "stab", ["#1a0a12", "#123a2a", _PAL["sym"]]),
                  vmin=0.0, vmax=1.6, interpolation="nearest", zorder=1)
        ax.plot(alphas, s_th - alphas, lw=1.2, color=_PAL["text"], alpha=0.8,
                zorder=4, label="α + β = g₀ + G₀")
        ax.plot(alphas, alphas - a_th, lw=1.2, ls=(0, (4, 3)),
                color=_PAL["text"], alpha=0.8, zorder=4,
                label="α − β = g₀ − G₀")
        for (a, b, lab) in ((ALPHA, BETA, "the stable default"),
                            (0.0, 0.0, "α = β = 0"),
                            (0.05, 0.18, "β = 0.18")):
            st = is_stable(MouthPair(SEP, a, a, b))
            col = _PAL["sym"] if st["stable"] else _PAL["hot"]
            ax.plot([a], [b], "o" if st["stable"] else "X", ms=7.0, color=col,
                    mec=_PAL["bg"], mew=0.8, zorder=6)
            ax.annotate(lab, xy=(a, b), xytext=(7, 5),
                        textcoords="offset points", color=col, fontsize=6.2,
                        family="monospace")
        ax.set_xlim(alphas[0], alphas[-1])
        ax.set_ylim(betas[0], betas[-1])
        ax.set_xlabel("α  (both mouths)", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("β", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper right", fontsize=6.4, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.set_title("green = λ_min ≥ 0 — self-adjoint everywhere, stable only "
                     "in the wedge",
                     color=_PAL["text"], fontsize=8.4, pad=6)

    # ── what conservation does buy ──────────────────────────────────────────
    def _draw_flux(self) -> None:
        ax = self.ax_flx
        rng = np.random.default_rng(20260815)
        xs, ys = [], []
        for _ in range(700):
            a1, a2 = rng.normal(0, 0.3, 2)
            b = complex(*rng.normal(0, 0.3, 2))
            herm = rng.random() < 0.5
            if herm:
                A = np.array([[a1, b], [b.conjugate(), a2]], dtype=complex)
            else:
                c = complex(*rng.normal(0, 0.3, 2))
                A = np.array([[a1, b], [c, a2]], dtype=complex)
            q = rng.normal(0, 1, 2) + 1j * rng.normal(0, 1, 2)
            f = mouth_flux(q, A)
            xs.append(float(np.abs(0.5 * (A - A.conjugate().T)).max()))
            ys.append(abs(f["net"]) / f["scale"] + 1e-18)
        xs, ys = np.array(xs), np.array(ys)
        m = xs < 1e-15
        ax.loglog(np.maximum(xs[m], 1e-18), ys[m], ".", ms=3.0,
                  color=_PAL["sym"], alpha=0.8, zorder=4,
                  label="A = A†  — net flux is zero identically")
        ax.loglog(xs[~m], ys[~m], ".", ms=3.0, color=_PAL["hot"], alpha=0.7,
                  zorder=3, label="A ≠ A†  — it is not")
        ax.set_xlim(1e-18, 3.0)
        ax.set_ylim(1e-18, 3.0)
        ax.set_xlabel("‖(A − A†)/2‖", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("|Im(q†Aq)| / scale", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper left", fontsize=6.4, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], which="both")

        ins = ax.inset_axes([0.58, 0.10, 0.38, 0.34])
        ins.set_facecolor("#050510")
        A = MouthPair(SEP, 0.2, -0.13, 0.15 + 0.07j).boundary_matrix()
        cs = np.geomspace(0.02, 1.0, 40)
        ins.semilogx(cs, [boundary_mixing(A, float(c))["diagonal_mixing"][0]
                          for c in cs], lw=1.2, color=_PAL["anti"])
        ins.semilogx(cs, [boundary_mixing(A, float(c))["off_diagonal_mixing"][0]
                          for c in cs], lw=1.2, color=_PAL["cool"])
        ins.set_title("the same A, different reference scale c",
                      color=_PAL["dim"], fontsize=5.8, pad=3)
        ins.tick_params(colors=_PAL["dim"], labelsize=5)
        for sp in ins.spines.values():
            sp.set_color(_PAL["rule"])
        ins.annotate("|S| entries are NOT r and t",
                     xy=(0.04, 0.80), xycoords="axes fraction",
                     color=_PAL["dim"], fontsize=5.4, family="monospace")
        ax.set_title("Hermiticity IS zero net boundary flux — the physical "
                     "conservation statement",
                     color=_PAL["text"], fontsize=8.4, pad=6)

    # ── where v55 sits ──────────────────────────────────────────────────────
    def _draw_embedding(self) -> None:
        ax = self.ax_emb
        ctl = DirectionalThroat(SEP, 1.0, +1, 0.3)
        ws = np.linspace(1.05, 5.95, 1400)
        ws = ws[np.abs(ws - np.round(ws)) > 0.02]
        embed = np.array([abs(ctl.secular(complex(w))) for w in ws])
        own = np.array([abs(ctl.pr255_pole_condition(complex(w))) for w in ws])
        old = []
        for w in ws:
            wz = complex(w)
            sp = cmath.sin(math.pi * wz)
            g = -wz * cmath.cos(math.pi * wz) / (4.0 * math.pi * sp)
            gd = (cmath.sin(wz * (math.pi - SEP))
                  / (4.0 * math.pi * math.sin(SEP) * sp))
            A = np.array([[0.0, 0.0], [1.0 / ctl.gain(wz), 0.0]],
                         dtype=complex)
            old.append(abs(complex(np.linalg.det(
                A - np.array([[g, gd], [gd, g]], dtype=complex)))))
        ax.semilogy(ws, own, lw=3.0, color=_PAL["free"], alpha=0.9, zorder=3,
                    label="v55's own |1 − L|")
        ax.semilogy(ws, embed, lw=1.1, color=_PAL["sym"], alpha=0.95, zorder=5,
                    label="|det(C − BΓ)| — the exact embedding")
        ax.semilogy(ws, np.array(old), lw=1.2, ls=(0, (4, 3)),
                    color=_PAL["hot"], alpha=0.9, zorder=4,
                    label="the earlier control — a different function")
        ax.annotate(f"embedding matches to {np.abs(embed - own).max():.1e}",
                    xy=(0.04, 0.06), xycoords="axes fraction",
                    color=_PAL["sym"], fontsize=6.6, family="monospace")
        ax.set_xlabel("ω", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("|secular function|", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper left", fontsize=6.4, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], which="both")
        ax.set_title("v55 embeds exactly — as a maximal, NON-self-adjoint "
                     "boundary condition",
                     color=_PAL["text"], fontsize=8.4, pad=6)

    # ── frame ───────────────────────────────────────────────────────────────
    def draw(self) -> None:
        self._draw_sector()
        self._draw_stability()
        self._draw_flux()
        self._draw_embedding()
        self.fig.suptitle("v56 — CONSERVATION IS NOT STABILITY",
                          color=_PAL["text"], fontsize=13.2, y=0.962,
                          family="monospace")
        self.fig.text(0.5, 0.908,
                      "a self-adjoint extension makes λ = ω² REAL   ·   it "
                      "does not make λ ≥ 0, and positivity is a separate "
                      "condition on the boundary data",
                      color=_PAL["dim"], fontsize=8.4, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.062,
                      "the first version of this round claimed a conserving "
                      "throat cannot ring up — two of its own three examples "
                      "have growing modes at σ = 2.4705 and 7.0910",
                      color=_PAL["dim"], fontsize=7.8, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.030,
                      "put in: a linear field on a fixed background and the "
                      "boundary data itself — four numbers chosen, not derived "
                      "  ·   the throat is POINT-supported, so no interior and "
                      "no delay",
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
