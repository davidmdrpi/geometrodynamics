#!/usr/bin/env python3
"""
Geometrodynamic QED — v58: static two-source throat tomography
=============================================================

**Not** the roadmap's two-wave collision invariant.  The object here is a
*static* source-interaction kernel at a fixed spectral parameter: it carries no
local null momenta, so it cannot distinguish equal-energy collinear from
counterpropagating waves — the load-bearing control behind
`𝒞 = I_A I_B (k_A·k_B)²`.  The index `(i, j)` labels **mouth channels**, not the
geodesic/winding branches of PRs #253–#255.  The dynamical object is still owed.

What it *is*: PR #253 closed rank counting by naming what it could not supply —
a quantity that **vanishes** when a source is removed rather than merely
becoming underdetermined.  Superposition makes every linear functional additive,
so the object has to be quadratic, and its cross term is

```
𝒞(y_A, y_B) = G(y_A,y_B) + Re Σ_ij G(y_A,c_i) R_ij G(c_j,y_B) ,
R = (A − Γ(λ))⁻¹
```

bilinear in the two source strengths, exactly zero when either is switched off.

What the panels show
────────────────────
**Top left — the mouth-channel matrix**, beside the same matrix for two
*disconnected* scatterers.  They look alike, and that is the point: `Γ` couples
the mouths through the ambient field whatever the boundary data says, so the
off-diagonal is a **cross-mouth** channel and not "through the throat".

**Top right — anisotropy is not the signature.**  Hold the geodesic separation
fixed and move one source over the sphere of that radius.  The free interaction
cannot move at all; the throat's varies by 66% of its mean — and two
**disconnected** scatterers vary by 69%.  A real effect that decides nothing.

**Bottom left — what does discriminate is a parameter count.**  The static
invariant determines three numbers, the entries of `S = Re R`; two independent
scatterers have two knobs, so their image is a *surface* with the exact equation
`S₁₂ = G₀ det S`.  The defect `𝒲 = S₁₂/det S − G₀` is zero on it, and on real
`β` equals **`−β`**.  Plotted against the Löwner margin: the invariant grows as
the cone's boundary is approached and `𝒲` does not move — the answer to PR
#255's caution that a resummed field measures the pole rather than the source.

**Bottom right — the blind family, and the two things that remove it.**  `𝒲 = 0`
has solutions away from `β = 0` only for **complex** `β`.  PR #257's gate
excludes the `Re β > G_d` branch (determinant negative).  And **reality of the
field** excludes the rest: a real scalar needs the self-adjoint domain to be
conjugation-invariant, `A = A*`, hence `β` real — with complex `β` a real static
source produces a *complex* field.  Inside a deliberately time-reversal-breaking
complex extension the family is real, and even there the limitation is the
*protocol*: phase-sensitive complex sources give the full complex `R` and hence
`A = Γ + R⁻¹` at **one** spectral parameter.

The one to remember
───────────────────
At the exact antipode `Γ(0)` is negative semidefinite, so the static response is
singular as `A → 0` and the invariant **diverges** like `1/ε` — while `𝒲` stays
exactly zero through four decades of it.  Size is not evidence.

What is put in
──────────────
A linear field on a fixed background, the mouth positions, and the boundary data
itself — four real numbers chosen, not derived.  Everything is evaluated at a
point **strictly inside** PR #257's cone, margin quoted, and the antipodal
endpoint is tested separately rather than approached.

Usage
─────
    python scripts/geometrodynamics_v58_two_source.py --still out.png
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.waves.throat_operator import MouthPair, gamma_at
from geometrodynamics.waves.throat_positivity import positivity_defect
from geometrodynamics.waves.two_source import (
    WORKING_BOUNDARY,
    WORKING_SEPARATION,
    defect_of_pair,
    free_interaction_energy,
    interaction_energy,
    invisible_partner,
    mouth_channel_invariant,
    random_points,
    recover_boundary,
    ring_points,
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
    "sym": "#7cff9e",          # the discriminator / stable
    "anti": "#ffb347",         # the disconnected null model
    "hot": "#ff6b8a",          # the raw invariant / blind
    "cool": "#5cc8ff",
}

SEP = WORKING_SEPARATION


def _pair() -> MouthPair:
    a1, a2, b = WORKING_BOUNDARY
    return MouthPair(SEP, a1, a2, b)


def _margin(pair: MouthPair) -> float:
    return float(positivity_defect(pair.boundary_matrix(),
                                   pair.separation)["min_eigenvalue"])


# ════════════════════════════════════════════════════════════════════════════
class TwoSourceFigure:
    """The mouth-channel matrix, the false signature, the true one, the
    scope."""

    def __init__(self, figsize=(13.8, 8.6)) -> None:
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, width_ratios=[1.0, 1.05], height_ratios=[1.0, 1.0],
            left=0.065, right=0.975, top=0.845, bottom=0.135,
            wspace=0.24, hspace=0.44)
        self.ax_br = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_iso = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_w = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_bl = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])

    # ── the mouth-channel matrix ────────────────────────────────────────────
    def _draw_channels(self) -> None:
        ax = self.ax_br
        pair = _pair()
        pts = random_points(2, seed=101)
        parts = mouth_channel_invariant(pair, pts[0], pts[1])
        block = np.array(parts["channels"], dtype=float)

        disc = MouthPair(SEP, pair.alpha1, pair.alpha2, 0.0)
        dblock = np.array(
            mouth_channel_invariant(disc, pts[0], pts[1])["channels"],
            dtype=float)
        scale = float(np.abs(np.vstack([block, dblock])).max())

        for k, (mat, name) in enumerate(((block, f"throat  β = "
                                          f"{complex(pair.beta).real:.2f}"),
                                         (dblock, "two DISCONNECTED "
                                          "scatterers  β = 0"))):
            sub = ax.inset_axes([0.06 + 0.50 * k, 0.50, 0.36, 0.32])
            sub.imshow(mat, cmap="magma", vmin=0.0, vmax=scale,
                       interpolation="nearest")
            for i in range(2):
                for j in range(2):
                    sub.text(j, i, f"{mat[i, j]:+.4f}", ha="center",
                             va="center", fontsize=6.4, family="monospace",
                             color=("#101018" if mat[i, j] > 0.55 * scale
                                    else _PAL["text"]))
            sub.set_xticks([0, 1]); sub.set_yticks([0, 1])
            sub.set_xticklabels(["out 1", "out 2"], fontsize=6.0,
                                color=_PAL["dim"])
            sub.set_yticklabels(["in 1", "in 2"], fontsize=6.0,
                                color=_PAL["dim"])
            sub.tick_params(length=0)
            for sp in sub.spines.values():
                sp.set_color(_PAL["rule"])
            sub.set_title(name, color=_PAL["dim"], fontsize=6.4, pad=4)

        ax.axis("off")
        ax.annotate(f"direct channel (∅, ∅):  {parts['direct']:+.4f}"
                    f"   —  {parts['direct'] / parts['throat_total']:.0f}× the "
                    f"whole throat block",
                    xy=(0.03, 0.995), xycoords="axes fraction", va="top",
                    color=_PAL["dim"], fontsize=6.4, family="monospace",
                    linespacing=1.8)
        ax.annotate("off-diagonal = CROSS-MOUTH: in one mouth, out the "
                    "other.\n"
                    "The two blocks look alike, and that is the point:\n"
                    "Γ already couples the mouths through the ambient\n"
                    "field, so an off-diagonal entry is NOT evidence of a\n"
                    "connection.  What separates them is a parameter count —\n"
                    "see the panel below.  (These are mouth channels, NOT\n"
                    "the geodesic branches of PRs #253–255.)",
                    xy=(0.045, 0.40), xycoords="axes fraction", va="top",
                    color=_PAL["hot"], fontsize=6.4, family="monospace",
                    linespacing=1.8)
        ax.set_title("the invariant, resolved on a pair of MOUTH CHANNELS",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # ── anisotropy ──────────────────────────────────────────────────────────
    def _draw_isotropy(self) -> None:
        ax = self.ax_iso
        pair = _pair()
        disc = MouthPair(SEP, pair.alpha1, pair.alpha2, 0.0)
        centre = np.array([math.cos(0.7), 0.0, math.sin(0.7), 0.0])
        ring = ring_points(centre, 1.0, 400, seed=20260816)
        c1 = np.array([1.0, 0.0, 0.0, 0.0])
        xs = np.array([math.acos(float(np.clip(np.dot(p, c1), -1, 1)))
                       for p in ring])
        free = np.array([free_interaction_energy(centre, p) for p in ring])
        thr = np.array([interaction_energy(pair, centre, p) for p in ring])
        dis = np.array([interaction_energy(disc, centre, p) for p in ring])
        ax.scatter(xs, thr, s=7, color=_PAL["sym"], alpha=0.75, zorder=5,
                   label=f"throat  β = {complex(pair.beta).real:.2f}")
        ax.scatter(xs, dis, s=7, color=_PAL["anti"], alpha=0.55, zorder=4,
                   marker="s",
                   label="two DISCONNECTED scatterers  β = 0")
        ax.plot([xs.min(), xs.max()], [free[0], free[0]], lw=1.8,
                color=_PAL["free"], zorder=6,
                label="free field — a single point, by theorem")
        rel = (thr.max() - thr.min()) / abs(thr.mean())
        drel = (dis.max() - dis.min()) / abs(dis.mean())
        ax.annotate(f"spread of the mean:  free "
                    f"{(free.max() - free.min()):.0e},  throat {rel:.0%},  "
                    f"disconnected {drel:.0%}\n"
                    "→ anisotropy detects structure at the mouths, not a "
                    "connection between them",
                    xy=(0.03, 0.055), xycoords="axes fraction",
                    color=_PAL["hot"], fontsize=6.4, family="monospace",
                    linespacing=1.8)
        ax.set_ylim(free[0] - 0.026, max(thr.max(), dis.max()) * 1.10)
        ax.set_xlabel("χ(y_B, mouth 1) — the sphere χ_AB = 1.0 is "
                      "2-dimensional, so this is a band",
                      color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("C", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper right", fontsize=6.4, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title("the same geodesic separation, everywhere on the panel",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # ── the discriminator against the margin ────────────────────────────────
    def _draw_defect(self) -> None:
        ax = self.ax_w
        gam = gamma_at(0.0, SEP).real
        beta = 0.06
        pts = random_points(2, seed=31)
        eps = np.logspace(math.log10(0.4), math.log10(2e-3), 26)
        vals, defs, margins = [], [], []
        for e in eps:
            a = float(gam[0, 0] + abs(beta - gam[0, 1]) + e)
            p = MouthPair(SEP, a, a, beta)
            margins.append(_margin(p))
            vals.append(abs(interaction_energy(p, pts[0], pts[1])))
            defs.append(defect_of_pair(p))
        ax.semilogx(margins, vals, lw=1.6, color=_PAL["hot"], zorder=5,
                    label="|C| — the raw invariant")
        ax.semilogx(margins, [-d for d in defs], lw=1.8, color=_PAL["sym"],
                    zorder=6, label="−W — the discriminator")
        ax.axhline(beta, color=_PAL["dim"], lw=0.9, ls=(0, (3, 3)), zorder=3)
        ax.annotate(f"β = {beta}", xy=(margins[-1], beta),
                    xytext=(6, 6), textcoords="offset points",
                    color=_PAL["dim"], fontsize=6.6, family="monospace")
        drift = max(abs(d + beta) for d in defs)
        ax.annotate(f"driven from margin {margins[0]:.2f} to "
                    f"{margins[-1]:.3f} the invariant grows\n"
                    f"{vals[-1] / vals[0]:.1f}×;  W drifts {drift:.0e} — the "
                    "discriminator does not\nsee the resonance (PR #255's "
                    "caution, answered)",
                    xy=(0.40, 0.86), xycoords="axes fraction", va="top",
                    color=_PAL["dim"], fontsize=6.4, family="monospace",
                    linespacing=1.8)
        ax.set_ylim(0.03, 0.27)
        ax.invert_xaxis()
        ax.set_xlabel("Löwner margin  λ_min(A − Γ(0))  →  the cone's boundary",
                      color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("value", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper left", fontsize=6.6, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], which="both")
        ax.set_title("W = S₁₂/det S − G₀  is exactly −β, and margin-blind",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # ── the blind spot ──────────────────────────────────────────────────────
    def _draw_blind(self) -> None:
        ax = self.ax_bl
        gd = float(gamma_at(0.0, SEP).real[0, 1])
        a1, a2 = 0.30, 0.35
        rbs = np.linspace(-0.16, -0.004, 160)
        ibs, margins = [], []
        for rb in rbs:
            ib = invisible_partner(a1, a2, float(rb), SEP)
            ibs.append(ib if ib is not None else np.nan)
            p = MouthPair(SEP, a1, a2, complex(rb, ib or 0.0))
            margins.append(_margin(p))
        ibs = np.array(ibs)
        margins = np.array(margins)
        ax.plot(rbs, ibs, lw=2.0, color=_PAL["hot"], zorder=6,
                label="W = 0 at λ = 0 — invisible")
        ax.axhline(0.0, color=_PAL["sym"], lw=1.6, zorder=5,
                   label="β real — where W = −β exactly")
        stable = margins > 0.0
        ax.scatter(rbs[stable][::10], ibs[stable][::10], s=18,
                   color=_PAL["sym"], edgecolor=_PAL["bg"], linewidth=0.5,
                   zorder=7, label="…and strictly inside PR #257's cone")
        ax.fill_between(rbs, 0.0, ibs, color=_PAL["hot"], alpha=0.07, zorder=1)
        # what a second frequency does to the same curve
        moved = []
        for rb in rbs:
            ib = invisible_partner(a1, a2, float(rb), SEP, 0.3)
            moved.append(ib if ib is not None else np.nan)
        ax.plot(rbs, moved, lw=1.5, ls=(0, (4, 3)), color=_PAL["cool"],
                zorder=6, label="W = 0 at λ = 0.3 — a different curve")
        # …and at λ = 0.8 the invisibility equation has no real root at all
        # here, so the blind set is not merely moved but empty
        gone = all(invisible_partner(a1, a2, float(rb), SEP, 0.8) is None
                   for rb in rbs)
        blind = MouthPair(SEP, a1, a2,
                          complex(-0.05, invisible_partner(a1, a2, -0.05, SEP)))
        rec = recover_boundary(blind)
        ax.annotate("every point here has Im β ≠ 0, so NONE of it is\n"
                    "compatible with a REAL scalar — a real source there\n"
                    "produces a complex field, and A = A* is exactly what\n"
                    "a conjugation-invariant domain means.\n"
                    f"PR #257 removes the other branch (Re β > G_d = "
                    f"{gd:.4f},\n"
                    "det(A − Γ) < 0);  reality removes this one.\n"
                    "and in a complex extension the limit is the PROTOCOL:\n"
                    "phase-sensitive sources give the full complex R at ONE\n"
                    f"λ, so A = Γ + R⁻¹ to {rec['max_parameter_error']:.0e}."
                    + ("  At λ = 0.8 the set is EMPTY." if gone else ""),
                    xy=(0.035, 0.55), xycoords="axes fraction", va="top",
                    color=_PAL["dim"], fontsize=5.9, family="monospace",
                    linespacing=1.55)
        ax.set_ylim(-0.62, 0.42)
        ax.set_xlabel("Re β", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("Im β", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="lower right", fontsize=6.4, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title("the blind family — stable, and outside the real-field "
                     "sector",
                     color=_PAL["text"], fontsize=8.6, pad=6)

    # ── frame ───────────────────────────────────────────────────────────────
    def draw(self) -> None:
        self._draw_channels()
        self._draw_isotropy()
        self._draw_defect()
        self._draw_blind()
        pair = _pair()
        self.fig.suptitle("v58 — STATIC TWO-SOURCE THROAT TOMOGRAPHY",
                          color=_PAL["text"], fontsize=13.2, y=0.962,
                          family="monospace")
        self.fig.text(0.5, 0.908,
                      "zero without a second source   ·   its DISCONNECTION "
                      "DEFECT is minus the mouth-mixing amplitude, W = −β"
                      "   ·   NOT yet the two-wave invariant",
                      color=_PAL["dim"], fontsize=8.4, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.070,
                      f"evaluated at A = ({pair.alpha1}, {pair.alpha2}, "
                      f"β = {complex(pair.beta).real}) with d = {SEP} — "
                      f"strictly inside PR #257's cone, Löwner margin "
                      f"{_margin(pair):.3f}",
                      color=_PAL["dim"], fontsize=7.6, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.047,
                      "the exact antipode is tested separately: C diverges "
                      "like 1/ε as A → 0 while W stays 0 to 4e-15   ·   size "
                      "is not evidence",
                      color=_PAL["dim"], fontsize=7.4, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.022,
                      "put in: a linear field on a fixed background, the mouth "
                      "positions, and the boundary data — four numbers chosen, "
                      "not derived   ·   the throat is POINT-supported",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")


# ════════════════════════════════════════════════════════════════════════════
def still(path: str) -> None:
    fig = TwoSourceFigure()
    fig.draw()
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v58.png")
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
