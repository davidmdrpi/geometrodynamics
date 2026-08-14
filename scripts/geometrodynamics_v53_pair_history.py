#!/usr/bin/env python3
"""
Geometrodynamic QED — v53: the interaction event is selected, not inserted
==========================================================================

v51 built one closed history — expanding leg, throat backwards in coordinate
time, collapsing leg. v52 established that pair creation is a threshold on an
invariant and needs **two** independently propagated waves. This round sews two
closed histories at one shared interaction and asks the question that comes
before any attempt at topology change:

    **is that event selected by the closure conditions, or still put in by hand?**

The system
──────────
Unknown: the event `C = (c, t_C)`. Given: two sources with launch times, two
throats with mouths and delays. Every leg is null, so a history closes in
coordinate time exactly on a **geodesic ellipsoid** — the locus whose summed
distance to the two mouths is `|Δ|`. Five equations, five unknowns:

```
|c|² = 1                                normalisation
d(S_A, c) = t_C − τ_A                   C lies on front A
d(S_B, c) = t_C − τ_B                   C lies on front B
d(c, M_A⁺) + d(M_A⁻, c) + Δ_A = 0       history A closes
d(c, M_B⁺) + d(M_B⁻, c) + Δ_B = 0       history B closes
```

What the two panels show
────────────────────────
**Left, both waves:** the two closure ellipsoids and the two null fronts, and
where all four meet — a small number of **isolated** events, each at full
Jacobian rank 5.

**Right, one wave removed:** the same scene with front B dropped. The solutions
do **not** disappear. The rank falls to 4 and they become a **one-parameter
family** — a curve threading the sphere. There is still a locus closing both
histories; there is no longer a *selected* one.

That difference is the whole result, and it is the reason the figure is two
panels rather than one.

What is put in
──────────────
The mouths and the delays: throat data, **given**, not solved for. That is where
the content lives — with `Δ` free, every event on both fronts closes and nothing
is selected. The conjugacy `Q_A + Q_B = 0` is a label, carried and checked.

**Not done:** no action principle, no field equations, no topology change, no
dynamics, no rate, and no worldline. The throats are identification maps on a
fixed round `S³`, and `shells.junction` priced them — two of them here, neither
paid for.

Usage
─────
    python scripts/geometrodynamics_v53_pair_history.py --still out.png
"""

from __future__ import annotations

import argparse
import math
from typing import List, Optional, Sequence, Tuple

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from geometrodynamics.viz.wormhole_ledger import shadow
from geometrodynamics.waves.pair_history import (
    PairHistorySystem,
    _random_system,
)

SEED = 3

_PAL = {
    "bg": "#010106",
    "panel": "#02020a",
    "grid": "#0e1a2a",
    "rule": "#1a2838",
    "text": "#e8ecf4",
    "dim": "#6a8aad",
    "faint": "#1b2d42",
    "ell_a": "#ff9d4d",        # closure ellipsoid of history A
    "ell_b": "#c07af0",        # closure ellipsoid of history B
    "front_a": "#ffe08a",      # null front of wave A
    "front_b": "#7ad4ff",      # null front of wave B
    "sel": "#7cff9e",          # a selected, isolated event
    "fam": "#ff6040",          # the one-parameter family
    "mouth": "#9c7ad0",
    "src": "#cfd8e3",
}


# ════════════════════════════════════════════════════════════════════════════
def _to2d(p3: np.ndarray, az: float = 0.64, el: float = 0.38) -> np.ndarray:
    ca, sa, ce, se = math.cos(az), math.sin(az), math.cos(el), math.sin(el)
    return np.array([p3[0] * ca - p3[1] * sa,
                     (p3[0] * sa + p3[1] * ca) * se + p3[2] * ce])


def _screen(q: Sequence[float]) -> Tuple[np.ndarray, float]:
    p3, depth = shadow(q)
    return _to2d(p3), depth


def _guide() -> List:
    segs = []
    ring = [np.array([math.cos(p), math.sin(p)])
            for p in np.linspace(0, 2 * math.pi, 220)]
    segs += [[ring[i], ring[i + 1]] for i in range(len(ring) - 1)]
    for axis in (2, 1):
        pts = []
        for p in np.linspace(0, 2 * math.pi, 180):
            v = np.zeros(3)
            v[(axis + 1) % 3] = math.cos(p)
            v[(axis + 2) % 3] = math.sin(p)
            pts.append(_to2d(v))
        segs += [[pts[i], pts[i + 1]] for i in range(len(pts) - 1)]
    return segs


def _ellipsoid_cloud(throat, n: int = 700_000, eps: float = 0.006,
                     seed: int = 0) -> np.ndarray:
    """Screen points of the closure locus, by rejection sampling on ``S³``.

    A sampled level set rather than a parametrised surface: the closure locus is
    ``{x : d(x,M⁺) + d(x,M⁻) = |Δ|}``, and drawing it by keeping points whose
    residual is under ``eps`` cannot accidentally draw something else.
    """
    rng = np.random.default_rng(seed)
    xs = rng.normal(size=(n, 4))
    xs /= np.linalg.norm(xs, axis=1)[:, None]
    tot = (np.arccos(np.clip(xs @ throat.m_plus, -1, 1))
           + np.arccos(np.clip(xs @ throat.m_minus, -1, 1)))
    keep = xs[np.abs(tot + throat.delay) < eps]
    return np.array([_screen(x)[0] for x in keep]) if len(keep) else np.zeros(
        (0, 2))


def _front_cloud(centre: np.ndarray, radius: float, n: int = 700_000,
                 eps: float = 0.006, seed: int = 0) -> np.ndarray:
    """Screen points of a null front, sampled the same way as the ellipsoids.

    Drawn as a level set rather than as a ring family on purpose: under the
    orthographic shadow a latitude/longitude mesh collapses edge-on into a comb
    of parallel segments that reads as an artefact rather than a sphere. Every
    locus in this figure is therefore the same kind of object — points whose
    defining residual is under ``eps`` — so none of them can be a drawing
    convention in disguise.
    """
    rng = np.random.default_rng(seed)
    xs = rng.normal(size=(n, 4))
    xs /= np.linalg.norm(xs, axis=1)[:, None]
    d = np.arccos(np.clip(xs @ np.asarray(centre, dtype=float), -1, 1))
    keep = xs[np.abs(d - radius) < eps]
    return np.array([_screen(x)[0] for x in keep]) if len(keep) else np.zeros(
        (0, 2))


# ════════════════════════════════════════════════════════════════════════════
class SelectionFigure:
    """Both waves against one wave: isolation against a family."""

    def __init__(self, seed: int = SEED, figsize=(13.6, 8.2)) -> None:
        rng = np.random.default_rng(seed)
        self.sys: PairHistorySystem = _random_system(rng)
        self.both = self.sys.solve(n_starts=320, seed=seed)
        self.family = self.sys.solve(n_starts=420, seed=seed + 1,
                                     with_b_wave=False, dedupe=2e-2)
        self.guide = _guide()
        self.cloud_a = _ellipsoid_cloud(self.sys.throat_a, seed=seed)
        self.cloud_b = _ellipsoid_cloud(self.sys.throat_b, seed=seed + 7)
        t_ref = self.both[0]["t"] if self.both else 1.0
        self.front_a = _front_cloud(self.sys.source_a, t_ref - self.sys.tau_a,
                                    seed=seed + 11)
        self.front_b = _front_cloud(self.sys.source_b, t_ref - self.sys.tau_b,
                                    seed=seed + 13)

        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, width_ratios=[1.1, 1.0], height_ratios=[1.0, 0.52],
            left=0.025, right=0.975, top=0.845, bottom=0.10,
            wspace=0.10, hspace=0.30)
        self.ax_both = self.fig.add_subplot(gs[:, 0], facecolor=_PAL["panel"])
        self.ax_one = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_book = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])

    # ── one scene ───────────────────────────────────────────────────────────
    def _scene(self, ax, both: bool, dot: float) -> None:
        S = self.sys
        ax.add_collection(LineCollection(self.guide, colors=[_PAL["faint"]],
                                         linewidths=0.4, alpha=0.5, zorder=1))
        for cloud, col in ((self.cloud_a, _PAL["ell_a"]),
                           (self.cloud_b, _PAL["ell_b"])):
            if len(cloud):
                ax.scatter(cloud[:, 0], cloud[:, 1], s=0.7, c=col, alpha=0.32,
                           linewidths=0, zorder=2)

        events = self.both if both else self.family
        clouds = [(self.front_a, _PAL["front_a"])]
        if both:
            clouds.append((self.front_b, _PAL["front_b"]))
        for cloud, col in clouds:
            if len(cloud):
                ax.scatter(cloud[:, 0], cloud[:, 1], s=0.7, c=col, alpha=0.30,
                           linewidths=0, zorder=3)

        for th, lab in ((S.throat_a, "A"), (S.throat_b, "B")):
            for m in (th.m_plus, th.m_minus):
                xy, _ = _screen(m)
                ax.plot([xy[0]], [xy[1]], "s", ms=4.5, color=_PAL["mouth"],
                        zorder=6, alpha=0.9)
        for s, lab in ((S.source_a, "S_A"), (S.source_b, "S_B")):
            xy, _ = _screen(s)
            ax.plot([xy[0]], [xy[1]], "^", ms=6, color=_PAL["src"], zorder=6)

        col = _PAL["sel"] if both else _PAL["fam"]
        for e in events:
            xy, _ = _screen(e["c"])
            ax.plot([xy[0]], [xy[1]], "o", ms=dot, color=col, zorder=8,
                    mec=_PAL["bg"], mew=0.6 if both else 0.0)

        ax.set_xlim(-1.16, 1.16)
        ax.set_ylim(-1.16, 1.16)
        ax.set_aspect("equal")
        ax.axis("off")

    # ── the books ───────────────────────────────────────────────────────────
    def _draw_book(self) -> None:
        ax = self.ax_book
        ax.set_facecolor(_PAL["panel"])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        S = self.sys
        inv = (S.invariant_at(self.both[0], energy=1.0) if self.both
               else {"s": 0.0, "opening_angle": 0.0})
        rows = [
            ("unknowns  (c ∈ S³, t)", "5"),
            ("equations  with both waves", "5"),
            ("selected events", f"{len(self.both)}"),
            ("Jacobian rank at each", f"{sorted({e['rank'] for e in self.both})}"
             if self.both else "—"),
            ("worst residual", f"{max((e['worst_residual'] for e in self.both), default=0):.1e}"),
            ("— with wave B removed —", ""),
            ("equations", "4"),
            ("Jacobian rank", f"{sorted({e['rank'] for e in self.family})}"
             if self.family else "—"),
            ("solutions sampled on the family", f"{len(self.family)}"),
            ("s at first selected event (E = m)", f"{inv['s']:.4f}"),
            ("threshold (2m)²", "4.0000"),
            ("net orientation  ↑ + ↓", f"{S.net_orientation():+d}"),
        ]
        for i, (k, v) in enumerate(rows):
            y = 0.955 - i * 0.079
            c = _PAL["dim"] if v else _PAL["fam"]
            ax.text(0.02, y, k, color=c, fontsize=7.6, family="monospace",
                    va="center")
            ax.text(0.98, y, v, color=_PAL["text"], fontsize=7.6,
                    family="monospace", va="center", ha="right")
        ax.set_title("five equations, five unknowns", color=_PAL["text"],
                     fontsize=8.3, pad=6)

    # ── frame ───────────────────────────────────────────────────────────────
    def draw(self) -> None:
        self._scene(self.ax_both, True, 10.0)
        self.ax_both.set_title(
            f"both waves  —  {len(self.both)} isolated events, Jacobian rank 5",
            color=_PAL["sel"], fontsize=9.6, pad=6)
        self._scene(self.ax_one, False, 3.2)
        self.ax_one.set_title(
            f"wave B removed  —  rank 4: a one-parameter family "
            f"({len(self.family)} sampled)",
            color=_PAL["fam"], fontsize=8.6, pad=6)
        self._draw_book()

        keys = [(_PAL["ell_a"], "closure locus of history A  (geodesic ellipsoid)"),
                (_PAL["ell_b"], "closure locus of history B"),
                (_PAL["front_a"], "null front of wave A"),
                (_PAL["front_b"], "null front of wave B"),
                (_PAL["sel"], "selected event  —  isolated"),
                (_PAL["fam"], "family  —  not selected")]
        for i, (col, lab) in enumerate(keys):
            y = 0.055 - i * 0.030 + 0.16
            self.ax_both.plot([0.02, 0.062], [y, y],
                              transform=self.ax_both.transAxes, lw=2.4,
                              color=col, alpha=0.9, zorder=9)
            self.ax_both.text(0.074, y, lab, transform=self.ax_both.transAxes,
                              color=_PAL["dim"], fontsize=6.6,
                              family="monospace", va="center", zorder=9)

        self.fig.suptitle("v53 — THE INTERACTION EVENT IS SELECTED, "
                          "NOT INSERTED", color=_PAL["text"], fontsize=13.2,
                          y=0.962, family="monospace")
        self.fig.text(0.5, 0.908,
                      "two closed histories sharing one event   ·   each closes "
                      "on a geodesic ellipsoid with its mouths as foci   ·   "
                      "the event is where both ellipsoids meet both null fronts",
                      color=_PAL["dim"], fontsize=8.4, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.042,
                      "removing a wave does not delete the solution — it drops "
                      "the rank from 5 to 4, and the isolated events become a "
                      "one-parameter family",
                      color=_PAL["dim"], fontsize=7.8, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.017,
                      "put in: the mouths and delays are given, not solved for "
                      "— with Δ free every event on both fronts closes   ·   "
                      "the shadow of S³ is 2-to-1, so an event may sit on a "
                      "locus's far sheet: the ranks are the evidence, not the "
                      "apparent coincidences",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")


# ════════════════════════════════════════════════════════════════════════════
def still(path: str, seed: int = SEED) -> None:
    fig = SelectionFigure(seed=seed)
    fig.draw()
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}  ({len(fig.both)} isolated, "
          f"{len(fig.family)} on the family)")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v53.png")
    ap.add_argument("--seed", type=int, default=SEED)
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still, a.seed)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
