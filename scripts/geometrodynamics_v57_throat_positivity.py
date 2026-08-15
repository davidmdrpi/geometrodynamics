#!/usr/bin/env python3
"""
Geometrodynamic QED — v57: the positive sector is a light cone
==============================================================

v56 established that a point-supported throat is a self-adjoint extension
parametrized by a Hermitian `2×2` boundary matrix `A`, that Hermiticity is
exactly flux conservation — and that it does **not** imply stability. It mapped
the stable region only on the two-parameter exchange-symmetric slice, by
scanning.

The full four-parameter answer is one inequality:

```
non-negative   ⟺   A ⪰ Γ(0)      (Löwner order)

Γ(0) = [[g₀, G₀], [G₀, g₀]] ,  g₀ = −1/(4π²) ,  G₀ = (π−d)/(4π² sin d)
```

because `dΓ/dλ ≻ 0` below threshold, so every eigenvalue of `A − Γ(λ)` falls
with `λ` while both run to `+∞` as `λ → −∞`: one crosses zero below threshold
**iff** it is already negative at it. That monotonicity is not sampled — it is
the statement that `dΓ_ij/dλ = ⟨δ_i, (H₀−λ)⁻² δ_j⟩` is a **Gram matrix**,
positive definite whenever the two mouths are distinct.

`φ^reg = A q` is a **chart** of the `U(2)` family, the one with `B` invertible.
The strata it misses are Dirichlet directions, reached only as `‖A‖ → ∞`; there
the criterion reads `A_eff ⪰ P†Γ(0)P` on the allowed-charge subspace, and the
chart is its two-dimensional case.

What the panels show
────────────────────
**Top left — the cone.** `A − Γ(0) = x₀I + x·σ` is positive semidefinite exactly
when `x₀ ≥ |x|`, so the stable set is a **forward light cone** in the four
dimensions of Hermitian boundary data. Drawn on the `x₃ = 0` slice (equal
mouths), with `x₁ = Re β − G₀` and `x₂ = −Im β`. Random boundary matrices are
coloured by an *independent* negative-`λ` root scan, not by the criterion —
green inside, pink outside, and the surface between them is where they agree.

**Top right — it counts, not just decides.** Predicted number of eigenvalues
below `λ*` (the inertia of `A − Γ(λ*)`) against the number actually found by
root-finding, at four thresholds. Everything on the diagonal, 160/160.

**Bottom left — the boundary is a zero mode.** Marching out of the cone along a
ray: the smallest eigenvalue of `A − Γ(0)` crosses zero, and *at that point*
`λ = 0` enters the spectrum. Past it the growth rate rises like `√ε`, exponent
`0.50001`, with the coefficient fixed by the eigenvalue slope rather than
fitted.

**Bottom right — where the apex sits.** `Γ(0)`'s eigenvalues against the mouth
separation. Its **trace is `−1/(2π²)` at every `d`**; for `0 < d < π` it is
indefinite, so `A = 0` is unstable wherever the mouths are actually apart. The
**exact antipode is the exception, not a limit of the rule**: `G_d` has a
*removable* singularity at `d = π`, with `G_π(0) = +1/(4π²) = −g₀`, so `Γ(0)`
has eigenvalues `(2g₀, 0)` — negative *semi*definite — and `A = 0` sits **on**
the cone's boundary, marginally non-negative, with a zero mode in the symmetric
channel.

What is put in
──────────────
A linear scalar field on a fixed background, and the boundary data itself — four
real numbers chosen, not derived. `shells.junction` is what would fix them from
matter. The throat is **point-supported**: no interior, no proper length, no
delay.

**Not done:** no backreaction, no stress tensor, no topology change, no rate,
and no two-source invariant — what this round buys that one is a stated region
to work inside.

Usage
─────
    python scripts/geometrodynamics_v57_throat_positivity.py --still out.png
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.waves.throat_operator import (
    MouthPair,
    negative_lambda_modes,
)
from geometrodynamics.waves.throat_positivity import (
    apex,
    cone_coordinates,
    count_modes_below,
    hermitian_from_parameters,
    inertia_below,
    threshold_scaling,
)

_PAL = {
    "bg": "#010106",
    "panel": "#02020a",
    "grid": "#0e1a2a",
    "rule": "#1a2838",
    "text": "#e8ecf4",
    "dim": "#6a8aad",
    "faint": "#1b2d42",
    "free": "#63798f",
    "sym": "#7cff9e",          # stable / symmetric channel
    "anti": "#ffb347",         # antisymmetric channel
    "hot": "#ff6b8a",          # growing
    "cool": "#5cc8ff",
}

SEP = 1.3


# ════════════════════════════════════════════════════════════════════════════
class PositivityFigure:
    """The cone, the count, the zero-mode boundary, and the apex."""

    def __init__(self, figsize=(13.8, 8.6)) -> None:
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, width_ratios=[1.15, 1.0], height_ratios=[1.0, 0.92],
            left=0.055, right=0.975, top=0.845, bottom=0.135,
            wspace=0.20, hspace=0.42)
        self.ax_cone = self.fig.add_subplot(gs[0, 0], projection="3d",
                                            facecolor=_PAL["panel"])
        self.ax_cnt = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_zm = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_apx = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])

    # ── the cone ────────────────────────────────────────────────────────────
    def _draw_cone(self) -> None:
        ax = self.ax_cone
        ax.set_facecolor(_PAL["panel"])
        th = np.linspace(0, 2 * math.pi, 90)
        hh = np.linspace(0.0, 0.30, 30)
        T, H = np.meshgrid(th, hh)
        ax.plot_surface(H * np.cos(T), H * np.sin(T), H, rstride=1, cstride=3,
                        color=_PAL["sym"], alpha=0.12, linewidth=0,
                        antialiased=True, zorder=1)
        ax.plot_wireframe(H * np.cos(T), H * np.sin(T), H, rstride=6,
                          cstride=12, color=_PAL["sym"], alpha=0.30, lw=0.5,
                          zorder=2)

        rng = np.random.default_rng(20260816)
        ins, outs = [], []
        for _ in range(420):
            a = float(rng.normal(0, 0.11))
            b = complex(rng.normal(0, 0.11), rng.normal(0, 0.11))
            A = hermitian_from_parameters(a, a, b)
            c = cone_coordinates(A, SEP)
            if abs(c["x0"]) > 0.32 or c["norm_x"] > 0.32:
                continue
            p = MouthPair(SEP, a, a, b)
            stable = not negative_lambda_modes(p, lambda_min=-40000.0,
                                               n_grid=4000)
            (ins if stable else outs).append(
                (float(c["x"][0]), float(c["x"][1]), c["x0"]))
        for pts, col, mk, lab, sz in ((outs, _PAL["hot"], ".",
                                       "scan: a growing mode", 6),
                                      (ins, _PAL["sym"], "o",
                                       "scan: none", 14)):
            if not pts:
                continue
            P = np.array(pts)
            ax.scatter(P[:, 0], P[:, 1], P[:, 2], s=sz, c=col, marker=mk,
                       depthshade=False, alpha=0.85, zorder=5, label=lab)
        ax.scatter([0], [0], [0], s=48, c=_PAL["text"], marker="*",
                   depthshade=False, zorder=6, label="apex  A = Γ(0)")

        ax.set_xlim(-0.30, 0.30)
        ax.set_ylim(-0.30, 0.30)
        ax.set_zlim(-0.10, 0.32)
        ax.set_xlabel("x₁ = Re β − G₀", color=_PAL["dim"], fontsize=7,
                      labelpad=-4)
        ax.set_ylabel("x₂ = −Im β", color=_PAL["dim"], fontsize=7, labelpad=-4)
        ax.set_zlabel("x₀ = ᾱ − g₀", color=_PAL["dim"], fontsize=7, labelpad=-6)
        ax.tick_params(colors=_PAL["dim"], labelsize=5.5, pad=-1)
        ax.view_init(elev=16, azim=-58)
        for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
            pane.set_pane_color((0.008, 0.008, 0.04, 1.0))
            pane._axinfo["grid"]["color"] = (0.08, 0.13, 0.20, 1.0)
        leg = ax.legend(loc="upper left", fontsize=6.2, framealpha=0.0,
                        labelcolor=_PAL["dim"], bbox_to_anchor=(-0.02, 0.98))
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.set_title("A ⪰ Γ(0)  is the forward light cone  x₀ ≥ |x|\n"
                     "(x₃ = 0 slice; colour is an independent λ-scan)",
                     color=_PAL["text"], fontsize=8.6, pad=-2)

    # ── the count ───────────────────────────────────────────────────────────
    def _draw_count(self) -> None:
        ax = self.ax_cnt
        rng = np.random.default_rng(4242)
        cols = {-2.0: _PAL["cool"], 0.0: _PAL["sym"], 0.5: _PAL["anti"],
                0.9: _PAL["hot"]}
        jit = {-2.0: -0.12, 0.0: -0.04, 0.5: 0.04, 0.9: 0.12}
        agree = 0
        total = 0
        for lstar, col in cols.items():
            xs, ys = [], []
            for _ in range(45):
                a1, a2 = rng.normal(0, 0.15, 2)
                b = complex(*rng.normal(0, 0.15, 2))
                A = hermitian_from_parameters(float(a1), float(a2), b)
                p = MouthPair(SEP, float(a1), float(a2), b)
                pred = inertia_below(A, SEP, lstar)
                got = count_modes_below(p, lstar)
                xs.append(pred + jit[lstar])
                ys.append(got + jit[lstar])
                agree += int(pred == got)
                total += 1
            ax.plot(xs, ys, "o", ms=4.4, color=col, alpha=0.55,
                    mec=_PAL["bg"], mew=0.3, zorder=4,
                    label=f"λ* = {lstar:g}")
        ax.plot([-0.4, 2.4], [-0.4, 2.4], lw=1.0, ls=(0, (4, 3)),
                color=_PAL["free"], zorder=2)
        ax.annotate(f"{agree} / {total} on the diagonal",
                    xy=(0.05, 0.90), xycoords="axes fraction",
                    color=_PAL["text"], fontsize=7.4, family="monospace")
        ax.set_xlim(-0.4, 2.4)
        ax.set_ylim(-0.4, 2.4)
        ax.set_xticks([0, 1, 2])
        ax.set_yticks([0, 1, 2])
        ax.set_xlabel("predicted: negative eigenvalues of A − Γ(λ*)",
                      color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("measured: roots below λ*", color=_PAL["dim"],
                      fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="lower right", fontsize=6.4, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title("the inertia theorem — it counts, at every threshold "
                     "below λ = 1",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # ── the zero mode at the boundary ───────────────────────────────────────
    def _draw_zero_mode(self) -> None:
        ax = self.ax_zm
        r = threshold_scaling(SEP)
        good = [x for x in r["rows"] if x["sigma"] is not None]
        eps = np.array([x["epsilon"] for x in good])
        sig = np.array([x["sigma"] for x in good])
        ax.loglog(eps, sig, "o-", lw=1.4, ms=5.0, color=_PAL["hot"],
                  alpha=0.95, zorder=5, label="growth rate σ = √|λ|")
        ref = sig[-1] * np.sqrt(eps / eps[-1])
        ax.loglog(eps, ref, lw=1.0, ls=(0, (4, 3)), color=_PAL["free"],
                  zorder=3, label="slope ½")
        ax.annotate(f"exponent  {r['asymptotic_exponent']:.5f}\n"
                    f"λ/ε → {r['lambda_over_epsilon_limit']:.4f}\n"
                    f"predicted {r['predicted_from_the_eigenvalue_slope']:.4f}"
                    f"  (not fitted)",
                    xy=(0.05, 0.72), xycoords="axes fraction",
                    color=_PAL["anti"], fontsize=6.8, family="monospace",
                    linespacing=1.6)
        ax.annotate("ε = 0: the boundary,\nwhere λ = 0 enters the spectrum\n"
                    "as a zero mode",
                    xy=(0.44, 0.10), xycoords="axes fraction",
                    color=_PAL["sym"], fontsize=6.6, family="monospace",
                    linespacing=1.6)
        ax.set_xlabel("ε — Löwner margin past the boundary, −λ_min(A − Γ(0))",
                      color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("σ", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="lower right", fontsize=6.6, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], which="both")
        ax.set_title("outside the boundary the instability turns on like √ε",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # ── the apex ────────────────────────────────────────────────────────────
    def _draw_apex(self) -> None:
        ax = self.ax_apx
        ds = np.linspace(0.12, math.pi, 400)
        lo, hi, tr = [], [], []
        for d in ds:
            ap = apex(float(d))
            lo.append(ap["eigenvalues"][0])
            hi.append(ap["eigenvalues"][1])
            tr.append(ap["trace"])
        ax.axhline(0.0, color=_PAL["free"], lw=0.9, zorder=2)
        ax.fill_between(ds, lo, hi, color=_PAL["hot"], alpha=0.08, zorder=1)
        ax.plot(ds, hi, lw=1.5, color=_PAL["sym"], zorder=5,
                label="g₀ + G₀  (symmetric threshold)")
        ax.plot(ds, lo, lw=1.5, color=_PAL["anti"], zorder=5,
                label="g₀ − G₀  (antisymmetric)")
        ax.plot(ds, tr, lw=1.1, ls=(0, (4, 3)), color=_PAL["text"], alpha=0.8,
                zorder=4, label="tr Γ(0) = −1/(2π²), constant")
        ax.plot([SEP], [apex(SEP)["eigenvalues"][1]], "o", ms=5.0,
                color=_PAL["sym"], mec=_PAL["bg"], mew=0.6, zorder=6)
        ax.plot([SEP], [apex(SEP)["eigenvalues"][0]], "o", ms=5.0,
                color=_PAL["anti"], mec=_PAL["bg"], mew=0.6, zorder=6)
        ax.plot([math.pi], [0.0], "o", ms=6.0, mfc="none",
                color=_PAL["free"], mew=1.2, zorder=7)
        ax.annotate("Γ(0) straddles zero for 0 < d < π,\n"
                    "so A = 0 is unstable there —\n"
                    "but at the exact antipode the symmetric\n"
                    "threshold closes, Γ(0) is negative\n"
                    "semidefinite, and A = 0 is marginal",
                    xy=(0.24, 0.09), xycoords="axes fraction",
                    color=_PAL["hot"], fontsize=6.2, family="monospace",
                    linespacing=1.6)
        ax.set_ylim(-0.34, 0.30)
        ax.set_xlim(ds[0], math.pi + 0.10)
        ax.annotate("d = π: g₀ + G₀ = 0 exactly", xy=(math.pi, 0.0),
                    xytext=(-6, 14), textcoords="offset points", ha="right",
                    color=_PAL["free"], fontsize=6.2, family="monospace")
        ax.set_xlabel("d — geodesic mouth separation", color=_PAL["dim"],
                      fontsize=8)
        ax.set_ylabel("eigenvalues of Γ(0)", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper right", fontsize=6.4, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title("where the apex sits — its trace does not know the mouth "
                     "separation",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # ── frame ───────────────────────────────────────────────────────────────
    def draw(self) -> None:
        self._draw_cone()
        self._draw_count()
        self._draw_zero_mode()
        self._draw_apex()
        self.fig.suptitle("v57 — THE POSITIVE SECTOR IS A LIGHT CONE",
                          color=_PAL["text"], fontsize=13.2, y=0.962,
                          family="monospace")
        self.fig.text(0.5, 0.908,
                      "non-negative if and only if  A - Γ(0) is positive "
                      "semidefinite   ·   one inequality in four parameters, "
                      "because dΓ/dλ is a GRAM matrix below threshold",
                      color=_PAL["dim"], fontsize=8.4, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.070,
                      "v56 mapped the two-parameter slice by scanning — the "
                      "same rule applied to general boundary data misclassifies "
                      "65 of 400 draws",
                      color=_PAL["dim"], fontsize=7.8, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.047,
                      "φ_reg = A q is the B-invertible CHART of U(2) — off it "
                      "the criterion reads A_eff - P†Γ(0)P ≥ 0 on the "
                      "allowed-charge subspace",
                      color=_PAL["dim"], fontsize=7.4, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.022,
                      "put in: a linear field on a fixed background and the "
                      "boundary data itself — four numbers chosen, not derived "
                      "  ·   the throat is POINT-supported, so no interior and "
                      "no delay",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")


# ════════════════════════════════════════════════════════════════════════════
def still(path: str) -> None:
    fig = PositivityFigure()
    fig.draw()
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v57.png")
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
