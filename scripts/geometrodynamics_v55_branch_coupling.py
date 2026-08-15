#!/usr/bin/env python3
"""
Geometrodynamic QED — v55: the throat solved for, not applied
=============================================================

v54 solved the field and got the strong result: on the Einstein static universe
the conformally coupled retarded Green function has **exact** image support, so
v53's ray branches are the field's branches, with the `1/(4π sin χ)` shell law
and a Maslov sign no path-length ledger could carry.

It did that with the throat on the outside. `φ(M⁺,t) = η φ(M⁻,t+Δ)` was applied
**to the free branches after they were computed** — one traversal, by
construction, because a post-processing step cannot notice that what it re-emits
will come back. Here the identification enters the equation that is solved:

```
a(ω) = ηκ e^{−iωΔ} [ S(ω) + T_d(ω) a(ω) ]   ⟹   a = ηκ e^{−iωΔ} S / (1 − L)
```

and the primitive object is indexed by a **pair of branches**, one per leg:

```
K_ab(ω) = ηκ · s_a A₁ e^{−u ℓ_a} · e^{−iωΔ} · s_b A₂ e^{−u ℓ_b}
```

What the panels show
────────────────────
**Top left — the two waveforms.** The solved field (amber/cyan by sign) against
the one-traversal control (grey). They agree on the primary arrivals and then
the solved field **rings where the control is flat**: echoes at
`ℓ_a + Δ + n(ℓ_c + Δ) + ℓ_b`, on a `κⁿ` ladder, each signed by the product of
every Maslov factor in its word. Those are events at times v54's ledger does not
contain — not corrections to amplitudes it does.

**Bottom left — the primitive.** Band coherence of `K_ab` over the branch pair.
`K_ab` carries the phase `e^{−iω(ℓ_a + Δ + ℓ_b)}`, so v53's closure condition
`ℓ_a + Δ + ℓ_b = 0` is *exactly* the statement that the entry is independent of
`ω`. The bright cells are the closed pairs — an **anti-diagonal**, not a
rectangle. The amplitude factorizes over this index (`K` is rank one); the
condition does not. That is why the pair is the primitive.

**Top right — where the loop lives.** `|L(ω)|`, the round-trip gain, and the
resolvent `|1/(1−L)|`. The peaks sit on the conformal ESU eigenfrequencies
`ω = n+1`, recovered here from the *image* representation: the branch series
sums to `(e^{−uχ} − e^{−u(2π−χ)})/(1 − e^{−2πu})`, whose poles are the spectrum.
`|L|` is also exactly the relative error of the one-traversal answer.

**Bottom right — and where post-processing stops working.** `κ_c = 1/max|T_d|`
against the regulator, log-log, slope `1.000`. As the damping is removed every
coupling is critical at some frequency, and there the one-traversal answer is
the first term of a divergent series rather than a leading term.

What is put in
──────────────
A linear scalar field on a fixed background, with the throat still an
**identification map** carrying a coupling `κ` by hand — `shells.junction`
priced it and the bill is inherited, unpaid. When `Δ + ℓ_c < 0` the loop is
closed in time and `1/(1−L)` is a self-consistency condition rather than a
history sum; it is unique exactly when the branch-resolved loop gain is
subcritical. **Not done:** no backreaction, no topology change, no rate, and no
two-source invariant — the two-throat fringe measured in the module is a
*throat–throat* interference, not that.

Usage
─────
    python scripts/geometrodynamics_v55_branch_coupling.py --still out.png
"""

from __future__ import annotations

import argparse
from typing import List, Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from geometrodynamics.waves.branch_coupling import (
    TWO_PI,
    CoupledThroat,
    coupled_arrivals,
    coupled_waveform,
    critical_coupling,
    leg_branches,
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
    "neg": "#5cc8ff",          # negative field
    "branch": "#7cff9e",
    "control": "#63798f",
    "hot": "#ff6b8a",
}

# coherence runs 0 → 1: dark where a pair dephases, green where it closes
_COH = LinearSegmentedColormap.from_list(
    "coh", ["#02020a", "#0d2033", "#15455a", "#2f8f76", _PAL["branch"]])

CHI_IN, CHI_OUT = 1.2, 0.9
SEP, DELAY, KAPPA = 1.3, 1.0, 0.60
DAMPING, WIDTH = 0.02, 0.05
T_MAX = 11.0


# ════════════════════════════════════════════════════════════════════════════
class CouplingFigure:
    """The solve against the control, the primitive, the loop, and its limit."""

    def __init__(self, figsize=(13.8, 8.4)) -> None:
        self.throat = CoupledThroat(SEP, DELAY, +1, KAPPA)
        self.ts, self.solved = coupled_waveform(
            T_MAX, CHI_IN, CHI_OUT, self.throat, DAMPING, WIDTH)
        _, self.control = coupled_waveform(
            T_MAX, CHI_IN, CHI_OUT, self.throat, DAMPING, WIDTH,
            resolvent=False)
        self.words = coupled_arrivals(CHI_IN, CHI_OUT, self.throat,
                                      T_MAX + 3.0, n_traversal=4, n_k=4,
                                      damping=DAMPING, rel_floor=1e-6)

        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, width_ratios=[1.35, 1.0], height_ratios=[1.0, 0.95],
            left=0.058, right=0.975, top=0.845, bottom=0.135,
            wspace=0.22, hspace=0.38)
        self.ax_wave = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_pair = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_gain = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_scal = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])

    # ── the experiment ──────────────────────────────────────────────────────
    def _draw_wave(self) -> None:
        ax = self.ax_wave
        ts, s, c = self.ts, self.solved, self.control
        # symlog: the echoes are 2% of the primaries and invisible on a linear
        # axis, and their being small is not the point — their existing is
        ax.set_yscale("symlog", linthresh=1e-5, linscale=0.35)
        ax.axhline(0.0, color=_PAL["faint"], lw=0.7)
        ax.plot(ts, c, lw=2.0, color=_PAL["control"], alpha=0.9, zorder=2,
                ls=(0, (5, 2)), label="control — one traversal (v54)")
        ax.plot(ts, s, lw=1.1, color=_PAL["text"], alpha=0.95, zorder=4,
                label="solved — the throat as a boundary condition")
        ax.fill_between(ts, 1e-30, s, where=s > 0, color=_PAL["pos"],
                        alpha=0.30, zorder=3)
        ax.fill_between(ts, -1e-30, s, where=s < 0, color=_PAL["neg"],
                        alpha=0.30, zorder=3)

        one = [w for w in self.words if w["traversals"] == 1]
        drawn: List[float] = []
        for w in self.words:
            if w["t"] > T_MAX or w["traversals"] == 1:
                continue
            if min((abs(w["t"] - o["t"]) for o in one), default=9e9) < 0.4:
                continue
            if any(abs(w["t"] - d) < 0.9 for d in drawn):
                continue
            drawn.append(w["t"])
            ax.axvline(w["t"], color=_PAL["hot"], lw=0.9, ls=(0, (3, 3)),
                       alpha=0.85, zorder=5)
            ax.annotate(f"echo ×{w['traversals']}  "
                        f"{'+' if w['sign'] > 0 else '−'}",
                        xy=(w["t"], 0.90), xycoords=("data", "axes fraction"),
                        color=_PAL["hot"], fontsize=6.6, ha="center",
                        family="monospace", zorder=6)
        for o in one:
            if o["t"] <= T_MAX:
                ax.axvline(o["t"], color=_PAL["faint"], lw=0.8, ls=":",
                           zorder=1)

        ax.set_ylim(-0.4, 0.4)
        ax.set_xlim(0, T_MAX)
        ax.set_xlabel("t", color=_PAL["dim"], fontsize=8.2)
        ax.set_ylabel("φ at the observer  (symlog)", color=_PAL["dim"],
                      fontsize=8.2)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="lower left", fontsize=6.8, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"])
        ax.set_title("the solved field rings where the control is flat — "
                     "arrivals v54's ledger does not contain",
                     color=_PAL["text"], fontsize=9.2, pad=6)

    # ── the primitive ───────────────────────────────────────────────────────
    def _draw_pair(self) -> None:
        ax = self.ax_pair
        n_k = 4
        ins, outs = leg_branches(CHI_IN, n_k), leg_branches(CHI_OUT, n_k)
        delay = -(CHI_IN + CHI_OUT + 2.0 * TWO_PI)   # closes a whole diagonal
        ws = np.linspace(0.5, 12.5, 1600)
        coh = np.zeros((len(ins), len(outs)))
        for i, a in enumerate(ins):
            for j, b in enumerate(outs):
                lag = a["ell"] + delay + b["ell"]
                coh[i, j] = abs(np.exp(-1j * ws * lag).mean())

        ax.imshow(coh, origin="lower", aspect="auto", cmap=_COH, vmin=0.0,
                  vmax=1.0, interpolation="nearest")
        for i, a in enumerate(ins):
            for j, b in enumerate(outs):
                if coh[i, j] > 0.999:
                    ax.add_patch(plt.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1, fill=False,
                        edgecolor=_PAL["pos"], lw=1.4, zorder=4))

        def lab(rs: List[dict]) -> List[str]:
            return [("L" if r["long_way"] else "S") + str(r["winding"])
                    for r in rs]

        ax.set_xticks(range(len(outs)))
        ax.set_xticklabels(lab(outs), fontsize=6.2, family="monospace")
        ax.set_yticks(range(len(ins)))
        ax.set_yticklabels(lab(ins), fontsize=6.2, family="monospace")
        ax.set_xlabel("branch b of the leg  M⁻ → observer", color=_PAL["dim"],
                      fontsize=8)
        ax.set_ylabel("branch a of the leg  source → M⁺", color=_PAL["dim"],
                      fontsize=8)
        ax.tick_params(colors=_PAL["dim"])
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        ax.annotate("bright = band average of e^{-iw(l_a + D + l_b)} is 1,\n"
                    "i.e. l_a + D + l_b = 0  —  3 pairs close,\n"
                    "9 sit in the rectangle they span",
                    xy=(0.985, 0.06), xycoords="axes fraction",
                    color=_PAL["dim"], fontsize=6.2, ha="right", va="bottom",
                    family="monospace", linespacing=1.5)
        ax.set_title("K_ab — closure is broadband coherence, and the closed "
                     "set is a diagonal, not a rectangle",
                     color=_PAL["text"], fontsize=8.4, pad=6)

    # ── the loop ────────────────────────────────────────────────────────────
    def _draw_gain(self) -> None:
        ax = self.ax_gain
        ws = np.linspace(0.05, 8.0, 4000)
        gain = np.array([abs(self.throat.loop_transfer(float(w), DAMPING))
                         for w in ws])
        res = np.array([abs(self.throat.resolvent(float(w), DAMPING))
                        for w in ws])
        ax.axhline(1.0, color=_PAL["hot"], lw=0.9, ls=(0, (4, 3)), alpha=0.8,
                   zorder=2)
        ax.annotate("|L| = 1 — the throat goes critical",
                    xy=(0.34, 1.0), xycoords=("axes fraction", "data"),
                    xytext=(0, -16), textcoords="offset points",
                    color=_PAL["hot"], fontsize=6.6, family="monospace")
        ax.semilogy(ws, gain, lw=1.4, color=_PAL["branch"], alpha=0.95,
                    zorder=4, label="|L(ω)|  round-trip gain")
        ax.semilogy(ws, res, lw=1.2, color=_PAL["pos"], alpha=0.9, zorder=3,
                    label="|1/(1−L)|  the resolvent")
        for n in range(1, 9):
            ax.axvline(float(n), color=_PAL["faint"], lw=0.7, ls=":", zorder=1)
        ax.annotate("ESU eigenfrequencies ω = n+1",
                    xy=(0.98, 0.06), xycoords="axes fraction",
                    color=_PAL["faint"], fontsize=6.4, ha="right",
                    family="monospace")
        ax.set_xlim(0, 8.0)
        ax.set_xlabel("ω", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper left", fontsize=6.6, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], which="both")
        ax.set_title("the loop lives on the spectrum — and |L| IS the error "
                     "of one traversal",
                     color=_PAL["text"], fontsize=8.4, pad=6)

    # ── and its limit ───────────────────────────────────────────────────────
    def _draw_scal(self) -> None:
        ax = self.ax_scal
        gs = np.array([0.16, 0.08, 0.04, 0.02, 0.01, 0.005])
        kc = np.array([critical_coupling(SEP, float(g))["kappa_critical"]
                       for g in gs])
        ax.loglog(gs, kc, "o-", lw=1.4, ms=4.6, color=_PAL["branch"],
                  alpha=0.95, zorder=4, label="κ_c = 1 / max|T_d|")
        ref = kc[0] * (gs / gs[0])
        ax.loglog(gs, ref, lw=1.0, ls=(0, (4, 3)), color=_PAL["control"],
                  zorder=3, label="slope 1")
        slope = float(np.polyfit(np.log(gs), np.log(kc), 1)[0])
        ax.annotate(f"fitted slope  {slope:.4f}",
                    xy=(0.05, 0.10), xycoords="axes fraction",
                    color=_PAL["pos"], fontsize=7.4, family="monospace")
        ax.set_xlabel("γ — the regulator", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("critical coupling κ_c", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        leg = ax.legend(loc="upper left", fontsize=6.6, framealpha=0.0,
                        labelcolor=_PAL["dim"])
        leg.get_frame().set_edgecolor(_PAL["rule"])
        ax.grid(alpha=0.08, color=_PAL["grid"], which="both")
        ax.set_title("as the regulator goes, every coupling is critical "
                     "somewhere",
                     color=_PAL["text"], fontsize=8.4, pad=6)

    # ── frame ───────────────────────────────────────────────────────────────
    def draw(self) -> None:
        self._draw_wave()
        self._draw_pair()
        self._draw_gain()
        self._draw_scal()
        self.fig.suptitle("v55 — THE THROAT SOLVED FOR, NOT APPLIED",
                          color=_PAL["text"], fontsize=13.2, y=0.962,
                          family="monospace")
        self.fig.text(0.5, 0.908,
                      "a(ω) = ηκ e^{−iωΔ}[S + T_d a]  solved, not iterated   ·   "
                      "the primitive K_ab is indexed by a PAIR of branches, one "
                      "per leg",
                      color=_PAL["dim"], fontsize=8.4, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.062,
                      "v54's answer is the n = 0 term of 1/(1−L) = Σ Lⁿ, and "
                      "its relative error is exactly the round-trip gain — the "
                      "solve adds EVENTS, not amplitudes",
                      color=_PAL["dim"], fontsize=7.8, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.030,
                      "put in: a linear field on a fixed background, the throat "
                      "still an identification map with κ by hand   ·   not "
                      "done: backreaction, topology change, and the two-source "
                      "invariant",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")


# ════════════════════════════════════════════════════════════════════════════
def still(path: str) -> None:
    fig = CouplingFigure()
    fig.draw()
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v55.png")
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
