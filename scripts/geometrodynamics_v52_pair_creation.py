#!/usr/bin/env python3
"""
Geometrodynamic QED — v52: pair creation is a collision, not a focus
====================================================================

Every earlier wave round in this arc drew **one** wavefront converging on its
antipode and treated the caustic as the interesting event. That was the wrong
object. A caustic is where the *amplitude* gets large; pair creation is a
threshold on an *invariant*:

```
γ γ → e⁺ e⁻      s = 2 E₁E₂ (1 − cos θ) ≥ (2m)²
```

`E` is what focusing raises. `θ` is what a collision supplies, and focusing does
not create one. So this round puts back what the original `S²` picture had and
the later ones dropped: **two** waveforms, one up and one down, from two nearby
sources — and lets the geometry say where they can actually make a pair.

What the geometry says, which is the point
──────────────────────────────────────────
Two sources a geodesic distance `δ` apart, both firing. Their fronts intersect
for `δ/2 ≤ t ≤ π − δ/2`, and at every crossing point

```
1 − cos θ = (1 − cos δ)/sin²t     ⟹     s(t) = 4 E₁E₂ sin²(δ/2) / sin²t
```

an identity of geodesic triangles, so the same on `S²` and `S³` — checked
against embedded tangent vectors that never use the law of cosines. It is
**U-shaped**: head-on (`s = 4E₁E₂`) at *both* ends of the crossing window, and
minimal at the equator where the fronts are largest and nearly collinear.

So a threshold cuts out **two disjoint windows**, never one. The near one sits
on top of the sources, where the fronts have not separated yet and nothing has
propagated independently. **Only the far one is a collision of independent
waves** — each front a half-circumference from home, meeting head-on again.
That is why the second interaction has to be antipodal, and it is derived rather
than staged.

What the picture shows
──────────────────────
Left: the two fronts on `S²`, in **two** colours because they are two waves —
the opposite of v51, where one colour meant one wave. At each crossing point the
two propagation directions are drawn as arrows, so head-on and collinear are
visible rather than asserted. When the invariant clears `(2m)²` the crossing
lights up, and the chord between the two crossing points — through the bulk —
is this program's reading of the created pair.

Right: `s(t)` against the threshold, with both windows shaded. The middle of the
trip is below threshold no matter how bright the waves are.

What is put in
──────────────
The Breit–Wheeler threshold and cross-section are **imported QED**. Treating a
front's null rays as photons is a correspondence, not a derivation. Calling the
chord a throat is interpretation, and `shells.junction` priced that throat — the
bill is inherited, not paid. Fixed round sphere, no backreaction, no rate.

Usage
─────
    python scripts/geometrodynamics_v52_pair_creation.py             # animate
    python scripts/geometrodynamics_v52_pair_creation.py --still out.png
    python scripts/geometrodynamics_v52_pair_creation.py --save out.gif
"""

from __future__ import annotations

import argparse
import math
from typing import List, Optional, Sequence, Tuple

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from geometrodynamics.waves.pair_creation import (
    TWO_M,
    WavePair,
    breit_wheeler_cross_section,
    outgoing_momentum,
)

FPS = 14
DELTA = 0.42
ENERGY = 1.4

_PAL = {
    "bg": "#010106",
    "panel": "#02020a",
    "grid": "#0e1a2a",
    "rule": "#1a2838",
    "text": "#e8ecf4",
    "dim": "#6a8aad",
    "faint": "#1b2d42",
    # TWO colours, because these are two waves — the opposite of v51
    "up": "#ffb347",
    "down": "#5cc8ff",
    "hot": "#7cff9e",          # the crossing, above threshold
    "cold": "#4a5a70",         # the crossing, below threshold
    "throat": "#e050ff",
    "src": "#9aa7b8",
    "anti": "#ff6040",
}


# ════════════════════════════════════════════════════════════════════════════
def _to2d(p3: np.ndarray, az: float = 0.72, el: float = 0.36) -> np.ndarray:
    """A fixed camera, so nothing moves except the physics."""
    ca, sa, ce, se = math.cos(az), math.sin(az), math.cos(el), math.sin(el)
    return np.array([p3[0] * ca - p3[1] * sa,
                     (p3[0] * sa + p3[1] * ca) * se + p3[2] * ce])


def _depth(p3: np.ndarray, az: float = 0.72, el: float = 0.36) -> float:
    """Toward-camera component, for the front/back cue."""
    ca, sa, ce, se = math.cos(az), math.sin(az), math.cos(el), math.sin(el)
    return float((p3[0] * sa + p3[1] * ca) * ce - p3[2] * se)


def _front(centre: np.ndarray, t: float, n: int = 220
           ) -> Tuple[List, np.ndarray]:
    """The geodesic circle of radius ``t`` about ``centre`` on ``S²``.

    The actual level set of geodesic distance, built from an orthonormal frame
    of the plane perpendicular to the centre — not a drawn ellipse.
    """
    c = centre / np.linalg.norm(centre)
    a = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(a, c))) > 0.9:
        a = np.array([0.0, 1.0, 0.0])
    u = a - float(np.dot(a, c)) * c
    u /= np.linalg.norm(u)
    v = np.cross(c, u)
    pts, dep = [], []
    for ph in np.linspace(0.0, 2.0 * math.pi, n):
        q = (math.cos(t) * c
             + math.sin(t) * (math.cos(ph) * u + math.sin(ph) * v))
        pts.append(_to2d(q))
        dep.append(_depth(q))
    segs = [[pts[i], pts[i + 1]] for i in range(len(pts) - 1)]
    mid = np.array([0.5 * (dep[i] + dep[i + 1]) for i in range(len(pts) - 1)])
    return segs, mid


def _sphere_guide() -> List:
    """Silhouette plus two great circles — a scale reference, nothing more."""
    segs = []
    ring = [np.array([math.cos(p), math.sin(p)])
            for p in np.linspace(0, 2 * math.pi, 240)]
    segs += [[ring[i], ring[i + 1]] for i in range(len(ring) - 1)]
    for axis in (2, 1):
        pts = []
        for p in np.linspace(0, 2 * math.pi, 200):
            q = np.zeros(3)
            q[(axis + 1) % 3] = math.cos(p)
            q[(axis + 2) % 3] = math.sin(p)
            pts.append(_to2d(q))
        segs += [[pts[i], pts[i + 1]] for i in range(len(pts) - 1)]
    return segs


# ════════════════════════════════════════════════════════════════════════════
class PairFigure:
    """The S² collision, the invariant against its threshold, and the readout."""

    def __init__(self, pair: Optional[WavePair] = None,
                 figsize=(13.6, 8.2)) -> None:
        self.pair = pair or WavePair(delta=DELTA, energy=ENERGY, dim=3)
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, width_ratios=[1.12, 1.0], height_ratios=[1.0, 0.58],
            left=0.025, right=0.972, top=0.850, bottom=0.095,
            wspace=0.15, hspace=0.36)
        self.ax_s2 = self.fig.add_subplot(gs[:, 0], facecolor=_PAL["panel"])
        self.ax_s = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_book = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])
        self.t0, self.t1 = self.pair.window
        self.guide = _sphere_guide()
        self.inset = None
        self.mom = None
        self.windows = self.pair.windows()

    # ── the S² scene ────────────────────────────────────────────────────────
    def _paint(self, ax, t: float, zoom: bool,
               arrow: float = 0.26) -> None:
        """The scene itself, painted into either the main axes or the inset.

        The same drawing code for both, so the inset cannot drift out of
        agreement with the wide view — it is a crop, not a second picture.
        """
        P = self.pair
        hot = P.above_threshold(t)
        lw = 2.4 if zoom else 1.5

        if not zoom:
            ax.add_collection(LineCollection(
                self.guide, colors=[_PAL["faint"]], linewidths=0.45,
                alpha=0.5, zorder=1))

        for centre, col in ((P.source_a, _PAL["up"]),
                            (P.source_b, _PAL["down"])):
            segs, dep = _front(centre, t, 520 if zoom else 220)
            rgba = np.tile(matplotlib.colors.to_rgba(col), (len(segs), 1))
            rgba[:, 3] = 0.22 + 0.66 * (dep + 1.0) / 2.0
            ax.add_collection(LineCollection(segs, colors=rgba, linewidths=lw,
                                             zorder=3))

        cross_col = _PAL["hot"] if hot else _PAL["cold"]
        try:
            locus = P.locus_at(t, samples=2)
        except ValueError:
            locus = np.zeros((0, 3))

        if len(locus) == 2:
            a2 = _to2d(locus[0] / np.linalg.norm(locus[0]))
            b2 = _to2d(locus[1] / np.linalg.norm(locus[1]))
            # the chord runs THROUGH the interior: this is the "via the bulk"
            # connection, and it is this program's reading, not a result
            ax.plot([a2[0], b2[0]], [a2[1], b2[1]], ls=(0, (2, 3)),
                    lw=1.5 if hot else 0.8,
                    color=_PAL["throat"] if hot else _PAL["faint"],
                    alpha=0.9 if hot else 0.4, zorder=4)

        for x in locus:
            x = x / np.linalg.norm(x)
            xy = _to2d(x)
            ax.plot([xy[0]], [xy[1]], "o", ms=(11 if zoom else 7) if hot
                    else (6 if zoom else 4), color=cross_col, zorder=7)
            # the two propagation directions, so head-on is VISIBLE and not
            # merely asserted by the caption
            for src, col in ((P.source_a, _PAL["up"]),
                             (P.source_b, _PAL["down"])):
                tip = _to2d(x + arrow * outgoing_momentum(src, x, t))
                ax.annotate("", xy=(tip[0], tip[1]), xytext=(xy[0], xy[1]),
                            arrowprops=dict(arrowstyle="-|>", lw=1.8,
                                            color=col, alpha=0.95), zorder=8)

        for pt, col, lab, ms in (
                (P.source_a, _PAL["up"], "source ↑", 6),
                (P.source_b, _PAL["down"], "source ↓", 6),
                (-P.source_a, _PAL["anti"], "antipode ↑", 4),
                (-P.source_b, _PAL["anti"], "antipode ↓", 4)):
            xy = _to2d(pt)
            ax.plot([xy[0]], [xy[1]], "o", ms=ms, color=col, zorder=6,
                    alpha=1.0 if "source" in lab else 0.8)
            if not zoom:
                ax.annotate(lab, xy=(xy[0], xy[1]), xytext=(7, 5),
                            textcoords="offset points", color=col, fontsize=7,
                            family="monospace",
                            alpha=1.0 if "source" in lab else 0.8)

    def _draw_s2(self, t: float) -> None:
        ax = self.ax_s2
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        for old in (self.inset, self.mom):
            if old is not None:
                old.remove()
        self.inset = self.mom = None
        P = self.pair
        hot = P.above_threshold(t)

        self._paint(ax, t, zoom=False, arrow=0.26)
        ax.set_xlim(-1.30, 1.30)
        ax.set_ylim(-1.22, 1.22)
        ax.set_aspect("equal")
        ax.axis("off")

        # A crop on the crossing, because at both threshold windows the fronts
        # are small circles near a pole and the whole event is a few pixels
        # wide in the wide view — which is where the interesting geometry is.
        # Centre on the crossing point facing the camera, at a FIXED
        # magnification.  Framing both points instead makes the zoom collapse
        # to ×1 near the equator — where they are nearly antipodal — and makes
        # it pump through the animation, which reads as motion that is not
        # there.
        half = 0.30
        try:
            locus = P.locus_at(t, samples=2)
            best = max((x / np.linalg.norm(x) for x in locus),
                       key=lambda q: _depth(q))
            c2 = _to2d(best)
        except ValueError:
            c2 = np.zeros(2)
        self.inset = ax.inset_axes([0.615, 0.02, 0.385, 0.385],
                                   facecolor="#04040e")
        # the arrows are drawn in world units, so at ×7 a fixed length walks
        # straight out of the crop — scale it to what the crop actually shows
        self._paint(self.inset, t, zoom=True, arrow=0.55 * half)
        self.inset.set_xlim(c2[0] - half, c2[0] + half)
        self.inset.set_ylim(c2[1] - half, c2[1] + half)
        ax.indicate_inset([c2[0] - half, c2[1] - half, 2 * half, 2 * half],
                          self.inset, edgecolor=_PAL["rule"], alpha=0.75,
                          linewidth=0.7)
        self.inset.set_aspect("equal")
        self.inset.set_xticks([])
        self.inset.set_yticks([])
        for sp in self.inset.spines.values():
            sp.set_color(_PAL["hot"] if hot else _PAL["rule"])
            sp.set_linewidth(1.0)
        self.inset.set_title(
            "the crossing, ×%.0f   —   arrows are projected"
            % (1.30 / half),
            color=_PAL["dim"], fontsize=7, pad=3, family="monospace")

        # a legend, so the inset can carry geometry and nothing else
        keys = [(_PAL["up"], "front from source ↑"),
                (_PAL["down"], "front from source ↓"),
                (_PAL["hot"] if hot else _PAL["cold"],
                 "crossing — pair" if hot else "crossing — below threshold"),
                (_PAL["throat"], "chord through the bulk  (interpretation)")]
        for i, (col, lab) in enumerate(keys):
            y = 0.985 - i * 0.042
            ax.plot([0.02, 0.075], [y, y], transform=ax.transAxes, lw=2.0,
                    color=col, alpha=0.9, zorder=9,
                    ls=(0, (2, 2)) if i == 3 else "-")
            ax.text(0.088, y, lab, transform=ax.transAxes, color=_PAL["dim"],
                    fontsize=6.8, family="monospace", va="center", zorder=9)

        # ── the opening angle, drawn where it is NOT distorted ──────────────
        # The arrows on the sphere are a PROJECTION of the two momenta, and
        # projection does not preserve angles: measured off the picture they
        # are wrong by up to 67°, and by different amounts at the two crossing
        # points, whose true opening angle is identical.  So the angle gets its
        # own axes, drawn in the plane the two momenta actually span.
        theta = P.theta_at(t)
        self.mom = ax.inset_axes([0.015, 0.02, 0.235, 0.235],
                                 facecolor="#04040e")
        m = self.mom
        m.plot([0, 1], [0, 0], lw=2.2, color=_PAL["up"], zorder=3,
               solid_capstyle="round")
        m.plot([0, math.cos(theta)], [0, math.sin(theta)], lw=2.2,
               color=_PAL["down"], zorder=3, solid_capstyle="round")
        m.plot([0, -1], [0, 0], lw=0.8, ls=":", color=_PAL["faint"], zorder=1)
        arc = np.linspace(0.0, theta, 90)
        m.plot(0.34 * np.cos(arc), 0.34 * np.sin(arc), lw=1.0,
               color=_PAL["hot"] if hot else _PAL["cold"], zorder=2)
        m.text(0.44 * math.cos(0.5 * theta), 0.44 * math.sin(0.5 * theta),
               f"{math.degrees(theta):.0f}°", color=_PAL["hot"] if hot
               else _PAL["cold"], fontsize=8, family="monospace",
               ha="center", va="center", zorder=4)
        m.text(-1.0, -0.16, "180° = head-on", color=_PAL["faint"],
               fontsize=6, family="monospace", ha="left", va="top")
        m.plot([0], [0], "o", ms=3, color=_PAL["text"], zorder=4)
        m.set_xlim(-1.15, 1.15)
        m.set_ylim(-0.45, 1.15)
        m.set_aspect("equal")
        m.set_xticks([])
        m.set_yticks([])
        for sp in m.spines.values():
            sp.set_color(_PAL["rule"])
        m.set_title("the opening angle, undistorted", color=_PAL["dim"],
                    fontsize=6.8, pad=3, family="monospace")

        state = "PAIR — s ≥ (2m)²" if hot else "below threshold"
        ax.set_title(f"S²   t = {t:.3f}    θ = {math.degrees(theta):5.1f}°"
                     f"    {state}",
                     color=_PAL["hot"] if hot else _PAL["text"],
                     fontsize=9.8, pad=6)

    # ── the invariant against its threshold ─────────────────────────────────
    def _draw_s(self, t: float) -> None:
        ax = self.ax_s
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        P = self.pair
        ts = np.linspace(self.t0, self.t1, 900)
        s = np.array([P.s_at(x) for x in ts])
        thr = (TWO_M * P.mass) ** 2

        for i, (w0, w1) in enumerate(self.windows):
            ax.axvspan(w0, w1, color=_PAL["hot"], alpha=0.12, zorder=0)
            ax.annotate("independent waves" if i else "source region",
                        xy=(0.5 * (w0 + w1), 0.90), xycoords=("data", "axes fraction"),
                        color=_PAL["hot"] if i else _PAL["dim"], fontsize=6.4,
                        ha="center", va="top", family="monospace")
        ax.axhline(thr, color=_PAL["hot"], lw=1.1, ls="--", alpha=0.8)
        ax.annotate("(2m)²", xy=(self.t1, thr), xytext=(-4, 5),
                    textcoords="offset points", color=_PAL["hot"],
                    fontsize=7, ha="right", family="monospace")
        ax.plot(ts, s, lw=2.0, color=_PAL["text"], alpha=0.9, zorder=4)
        ax.plot([t], [P.s_at(t)], "o", ms=7, zorder=6,
                color=_PAL["hot"] if P.above_threshold(t) else _PAL["cold"])
        ax.axvline(t, color=_PAL["text"], lw=0.8, alpha=0.45)
        ax.axvline(math.pi / 2, color=_PAL["dim"], lw=0.7, ls=":", alpha=0.6)
        ax.annotate("equator — fronts largest,\ninvariant smallest",
                    xy=(math.pi / 2, thr * 0.16), xytext=(6, 0),
                    textcoords="offset points", color=_PAL["dim"],
                    fontsize=6.4, family="monospace", va="center")
        ax.set_yscale("log")
        ax.set_xlim(self.t0, self.t1)
        ax.set_ylim(min(s) * 0.45, max(s) * 2.6)
        ax.set_xlabel("t — how far each front has propagated",
                      color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("s = 4E²sin²(δ/2)/sin²t", color=_PAL["dim"], fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        ax.grid(alpha=0.10, color=_PAL["grid"])
        ax.set_title("head-on at both ends, collinear in the middle — "
                     "so the threshold cuts two windows, never one",
                     color=_PAL["text"], fontsize=8.3, pad=6)

    # ── the readout ─────────────────────────────────────────────────────────
    def _draw_book(self, t: float) -> None:
        ax = self.ax_book
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        P = self.pair
        w = self.windows
        far = w[1] if len(w) == 2 else (None, None)
        rows = [
            ("source separation δ", f"{P.delta:.3f}"),
            ("E / m  (a parameter)", f"{P.energy / P.mass:.3f}"),
            ("opening angle θ now", f"{math.degrees(P.theta_at(t)):.2f}°"),
            ("s now", f"{P.s_at(t):.4f}"),
            ("threshold (2m)²", f"{(TWO_M * P.mass) ** 2:.4f}"),
            ("σ(s) now  [π r_e²]", f"{P.cross_section_at(t):.4f}"),
            ("threshold windows", f"{len(w)}"),
            ("antipodal window opens at",
             "—" if far[0] is None else f"t = {far[0]:.4f}"),
            ("net orientation  ↑ + ↓", f"{P.net_orientation():+d}"),
        ]
        for i, (k, v) in enumerate(rows):
            y = 0.94 - i * 0.104
            ax.text(0.02, y, k, color=_PAL["dim"], fontsize=7.8,
                    family="monospace", va="center")
            ax.text(0.98, y, v, color=_PAL["text"], fontsize=7.8,
                    family="monospace", va="center", ha="right")
        ax.set_title("the collision, counted", color=_PAL["text"],
                     fontsize=8.3, pad=6)

    # ── frame ───────────────────────────────────────────────────────────────
    def draw(self, t: float) -> None:
        self._draw_s2(t)
        self._draw_s(t)
        self._draw_book(t)
        self.fig.suptitle("v52 — PAIR CREATION IS A COLLISION, NOT A FOCUS",
                          color=_PAL["text"], fontsize=13.5, y=0.963,
                          family="monospace")
        self.fig.text(0.5, 0.912,
                      "two waves, one up one down, from sources δ apart   ·   "
                      "s = 2E₁E₂(1 − cos θ):  focusing raises E, only a "
                      "collision supplies θ",
                      color=_PAL["dim"], fontsize=8.6, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.040,
                      "derived: 1 − cos θ = (1 − cos δ)/sin²t, an identity of "
                      "geodesic triangles — head-on at both ends, and only the "
                      "antipodal window is a collision of independently "
                      "propagated waves",
                      color=_PAL["dim"], fontsize=7.6, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.016,
                      "put in: the Breit–Wheeler threshold and cross-section "
                      "are QED; rays-as-photons is a correspondence; the chord "
                      "through the bulk is interpretation, and PR #249 priced "
                      "that throat",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")


# ════════════════════════════════════════════════════════════════════════════
def still(path: str) -> None:
    fig = PairFigure()
    w = fig.windows
    t = 0.5 * (w[1][0] + w[1][1]) if len(w) == 2 else 0.5 * (fig.t0 + fig.t1)
    fig.draw(t)
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def animate(save: Optional[str] = None, frames: int = 150):
    from matplotlib import animation

    holder = PairFigure()
    eps = 1e-4
    times = np.linspace(holder.t0 + eps, holder.t1 - eps, frames)

    def update(i: int):
        holder.draw(float(times[min(i, frames - 1)]))
        return []

    anim = animation.FuncAnimation(holder.fig, update, frames=frames,
                                   interval=1000 // FPS, blit=False)
    if save:
        anim.save(save, fps=FPS, dpi=100,
                  savefig_kwargs={"facecolor": _PAL["bg"]})
        print(f"wrote {save}")
    else:
        plt.show()
    return anim


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG")
    ap.add_argument("--save", metavar="FILE")
    ap.add_argument("--frames", type=int, default=150)
    a = ap.parse_args(argv)
    if a.still:
        matplotlib.use("Agg")
        still(a.still)
        return 0
    if a.save:
        matplotlib.use("Agg")
    animate(save=a.save, frames=a.frames)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
