#!/usr/bin/env python3
"""
Geometrodynamic QED — v51: one conserved wave, seen in pieces
=============================================================

An emitter fires a shell on `S³`. While it expands, a *collapsing* shell sweeps
past it and converges on a receiver, which recoils from the momentum delivered.
Locally that is two unrelated events. Globally it is one object:

```
E ──expand──▶ M_future ──throat, Δt < 0──▶ M_past ──collapse──▶ R
```

The emitted shell reaches the **antipode**, where `S³` refocuses it exactly, and
the antipode carries a wormhole mouth. Through the throat it re-emerges at a
mouth displaced into the *past* — and that emergence **is** the collapsing
shell, the one that appeared to sweep past the emitter and land on the receiver.

What the reference movie could not show
───────────────────────────────────────
`antipodal_crossing` (v10) drew worldlines and their antipodal traces, and had
no notion of the network topology: nothing said that the incoming and outgoing
shells were the *same* wave. Here they are drawn in **one colour**, and the
right-hand panel closes the circuit explicitly — path length against coordinate
time, the one axis on which the two locally-unrelated shells are one continuous
object. The throat is the horizontal jump backwards.

Why the destinations are antipodal
──────────────────────────────────
Not staging, and used **twice**. A geodesic sphere at distance `χ` on `S³` has
area `4π sin²χ`, so an energy-conserving shell has amplitude `∝ 1/sin χ` and
refocuses exactly at `χ = π`. That is why the future mouth sits at the emitter's
antipode, and why the receiver sits at the past mouth's: it is the only place
the arriving shell is genuinely *collapsing* rather than merely arriving —
`dA/dχ = 4π sin 2χ < 0`, checked in `measure_the_receiver_is_struck_by_a_
collapsing_shell` against a displaced receiver, where the same wave is still
expanding when it lands.

The projection, and what it is not evidence for
───────────────────────────────────────────────
The scene is the **orthographic shadow** of `S³`: project `R⁴ → R³` along a
fixed direction `n`, then a fixed camera.

Stereographic projection was the first choice and is wrong for this figure. It
is unbounded **at the pole it projects from**, and a shell launched from a point
sweeps *all* of `S³` — so whatever pole is chosen, the shell crosses it once and
the image blows up there: the radius grows as `2/ε` and never converges. The
first draft projected from `q₃ = +1`, which is the emitter's own position, so
the emitter was a division by zero and never got drawn at all.

The shadow is bounded by the unit ball always, and a shell of radius `χ` has
diameter `≤ 2 sin χ`, so the refocus is *visible*: the shell shrinks to a point.
Better still, its screen extent is proportional to `sin χ` with one constant —
`√(A/4π)` — so the drawn size **is** the area law.

The price is that the shadow is **2-to-1** — the near and far halves of `S³`
land on the same ball, distinguished here only by brightness. **A crossing on
screen is not a crossing on `S³`**, and no claim in this figure rests on one.

What is derived and what is put in
──────────────────────────────────
**Derived:** the antipodal focus, and the self-consistency — a linear wave on
this loop has exactly one amplitude, `A = A_src/(1 − κ)`, matching brute
iteration of 2000 round trips, unique for every `κ ≠ 1`. Nothing is tuned.

**Put in:** the wormhole is an *identification map*, not a solution. `κ` and the
delay `Δ` are parameters. PR #249 priced this throat — a minimal surface has
`σ < 0` identically — so the bill is inherited, not paid.

`Δ` is a dial: it decides whether the receiver is struck before the emitter
fires, and it changes **no conserved quantity at all**. The value used here is
chosen so that both shells are on screen together, which is the informative
movie, not because the bookkeeping prefers it.

Usage
─────
    python scripts/geometrodynamics_v51_wormhole_ledger.py             # animate
    python scripts/geometrodynamics_v51_wormhole_ledger.py --still out.png
    python scripts/geometrodynamics_v51_wormhole_ledger.py --save out.gif
"""

from __future__ import annotations

import argparse
import math
from typing import List, Optional, Sequence, Tuple

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from geometrodynamics.transaction.s3_geometry import nrm4
from geometrodynamics.viz.wormhole_ledger import (
    WormholeLoop,
    geodesic_shell,
    shadow,
    shell_amplitude,
)

FPS = 14
KAPPA = 0.45
# Chosen so both shells are live together for a stretch of the movie: the past
# mouth fires at t = π + Δ < 0, the returning shell sweeps the emitter just
# before it fires, and lands just after.  At the module's default Δ = −12 the
# whole return leg happens before the emission and the ledger does not notice.
DELAY = -5.0

_PAL = {
    "bg": "#010106",
    "panel": "#02020a",
    "grid": "#0e1a2a",
    "rule": "#1a2838",
    "text": "#e8ecf4",
    "dim": "#6a8aad",
    "faint": "#1b2d42",
    # ONE colour for the wave, because it is one wave
    "wave": "#ffd166",
    "emitter": "#70ff38",
    "receiver": "#ff6040",
    "mouth_f": "#38d8ff",
    "mouth_p": "#e050ff",
    "throat": "#e050ff",
    "recoil": "#ff6040",
    "both": "#123048",
}


def _to2d(p3: np.ndarray, az: float = 0.62, el: float = 0.40) -> np.ndarray:
    """A fixed orthographic camera, so nothing moves except the physics."""
    ca, sa, ce, se = math.cos(az), math.sin(az), math.cos(el), math.sin(el)
    return np.array([p3[0] * ca - p3[1] * sa,
                     (p3[0] * sa + p3[1] * ca) * se + p3[2] * ce])


def _screen(q: Sequence[float]) -> Tuple[np.ndarray, float]:
    p3, depth = shadow(q)
    return _to2d(p3), depth


def _shell_segments(centre: np.ndarray, chi: float, n_theta: int = 22,
                    n_phi: int = 44) -> Tuple[List, np.ndarray]:
    """Screen segments of the geodesic 2-sphere of radius ``chi``.

    The *geometry* comes from ``wormhole_ledger.geodesic_shell``, which is where
    it can be measured — and is, to `3e-15` in geodesic distance.  This function
    only projects and stitches, so the drawing cannot quietly disagree with the
    thing the probe checks.
    """
    pts = geodesic_shell(centre, chi, n_theta, n_phi)
    segs, depths = [], []
    for r in range(0, len(pts), n_phi):
        ring = [_screen(p) for p in pts[r:r + n_phi]]
        for i in range(len(ring) - 1):
            segs.append([ring[i][0], ring[i + 1][0]])
            depths.append(0.5 * (ring[i][1] + ring[i + 1][1]))
    return segs, np.asarray(depths)


def _unit_sphere_wire() -> List:
    """The shadow boundary — the great 2-sphere ``q·n = 0``, drawn for scale.

    Deliberately a *silhouette plus two great circles*, not a full wireframe: a
    latitude/longitude globe here is the same drawing idiom as the shells, and
    at the sizes involved the guide and the expanding shell became impossible
    to tell apart.
    """
    segs = []
    ring = [np.array([math.cos(p), math.sin(p)])
            for p in np.linspace(0, 2 * math.pi, 200)]
    segs += [[ring[i], ring[i + 1]] for i in range(len(ring) - 1)]
    for axis in (2, 1):
        pts = []
        for p in np.linspace(0, 2 * math.pi, 160):
            v = np.zeros(3)
            v[(axis + 1) % 3] = math.cos(p)
            v[(axis + 2) % 3] = math.sin(p)
            pts.append(_to2d(v))
        segs += [[pts[i], pts[i + 1]] for i in range(len(pts) - 1)]
    return segs


# ════════════════════════════════════════════════════════════════════════════
class LedgerFigure:
    """The S³ scene, the single wave path, and the books."""

    def __init__(self, loop: Optional[WormholeLoop] = None,
                 figsize=(13.6, 8.2)) -> None:
        self.loop = loop or WormholeLoop(kappa=KAPPA, delay=DELAY)
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, width_ratios=[1.16, 1.0], height_ratios=[1.0, 0.62],
            left=0.025, right=0.972, top=0.855, bottom=0.095,
            wspace=0.14, hspace=0.34)
        self.ax_s3 = self.fig.add_subplot(gs[:, 0], facecolor=_PAL["panel"])
        self.ax_loop = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_book = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])

        t = self.loop.arrival_times()
        self.t_emit = t["emit"]
        self.t_future = t["reach_future_mouth"]
        self.t_past = t["leave_past_mouth"]
        self.t_recv = t["strike_receiver"]
        self.t_sweep = t["sweep_past_emitter"]
        self.t0 = min(self.t_past, self.t_emit) - 0.9
        self.t1 = max(self.t_future, self.t_recv) + 0.9
        self.wire = _unit_sphere_wire()

    # ── the S³ scene ────────────────────────────────────────────────────────
    def _draw_s3(self, t: float) -> None:
        ax = self.ax_s3
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        L = self.loop
        ax.add_collection(LineCollection(self.wire, colors=[_PAL["faint"]],
                                         linewidths=0.45, alpha=0.5, zorder=1))

        live = []
        if self.t_emit <= t <= self.t_future:
            live.append(("emitted", L.emitter.position, t - self.t_emit))
        if self.t_past <= t <= self.t_recv:
            live.append(("returning", L.past_mouth.position, t - self.t_past))

        for _, centre, chi in live:
            if not 1e-3 < chi < math.pi - 1e-3:
                continue
            segs, dep = _shell_segments(centre, chi)
            amp = min(shell_amplitude(chi) / 5.0, 1.0)
            # ONE colour: the outgoing and returning shells are the same wave.
            # Alpha carries the depth q·n, the coordinate the shadow discards —
            # the far half of S³ is drawn faint rather than hidden.
            rgba = np.tile(matplotlib.colors.to_rgba(_PAL["wave"]),
                           (len(segs), 1))
            rgba[:, 3] = (0.13 + 0.42 * amp) * (0.35 + 0.65 * (dep + 1.0) / 2.0)
            ax.add_collection(LineCollection(
                segs, colors=rgba, linewidths=0.55 + 1.5 * amp, zorder=3))
            # name the phase in place, in the wave's own colour: the two labels
            # are the same object, and the sign of dA/dχ is what separates them
            xy, _ = _screen(math.cos(chi) * nrm4(centre))
            word = "expanding" if chi < math.pi / 2 else "collapsing"
            ax.annotate(f"{word}   χ = {chi:.2f}", xy=(xy[0], xy[1]),
                        xytext=(0, -19), textcoords="offset points",
                        color=_PAL["wave"], fontsize=7.0, family="monospace",
                        ha="center", va="center", alpha=0.9, zorder=8)

        for mouth, col, ms in ((L.emitter, _PAL["emitter"], 7),
                               (L.future_mouth, _PAL["mouth_f"], 8),
                               (L.past_mouth, _PAL["mouth_p"], 8),
                               (L.receiver, _PAL["receiver"], 7)):
            xy, _ = _screen(mouth.position)
            ax.plot([xy[0]], [xy[1]], "o", ms=ms, color=col, zorder=6)
            ax.annotate(mouth.label, xy=(xy[0], xy[1]), xytext=(8, 6),
                        textcoords="offset points", color=col, fontsize=7.5,
                        family="monospace")

        # the throat, drawn as what it is: an identification, not a path
        a, _ = _screen(L.future_mouth.position)
        b, _ = _screen(L.past_mouth.position)
        ax.plot([a[0], b[0]], [a[1], b[1]], ls=(0, (2, 3)), lw=1.2,
                color=_PAL["throat"], alpha=0.7, zorder=4)
        mid = 0.35 * a + 0.65 * b
        ax.annotate("throat: an identification, Δt = %+.1f"
                    % L.future_mouth.delay,
                    xy=(mid[0], mid[1]), xytext=(12, -17),
                    textcoords="offset points", color=_PAL["throat"],
                    fontsize=6.8, family="monospace", ha="left", alpha=0.85)

        if abs(t - self.t_sweep) < 0.16:
            xy, _ = _screen(L.emitter.position)
            ax.plot([xy[0]], [xy[1]], "o", ms=17, mfc="none", mew=1.4,
                    color=_PAL["wave"], zorder=7)
        if t >= self.t_recv:
            xy, _ = _screen(L.receiver.position)
            ax.annotate("", xy=(xy[0] + 0.34, xy[1] - 0.24),
                        xytext=(xy[0], xy[1]),
                        arrowprops=dict(arrowstyle="-|>", lw=1.7,
                                        color=_PAL["recoil"]), zorder=7)

        ax.set_xlim(-1.18, 1.18)
        ax.set_ylim(-1.18, 1.18)
        ax.set_aspect("equal")
        ax.axis("off")
        # The phase is read off the SIGN of dA/dχ, not off which leg the wave
        # is on.  An earlier version of this title called the whole return leg
        # "collapsing" because that is what the story says — while the in-scene
        # label, which does read the sign, said "expanding".  Same error the
        # module's receiver placement had, reintroduced in a caption.
        if len(live) == 2:
            phase = "both shells live — and they are one wave"
        elif live:
            name, _, chi = live[0]
            phase = (f"{name} shell — "
                     + ("expanding" if chi < math.pi / 2 else "collapsing"))
        else:
            phase = "—"
        ax.set_title(f"S³, orthographic shadow    t = {t:+.2f}    {phase}",
                     color=_PAL["text"], fontsize=9.6, pad=6)

    # ── the single wave path ────────────────────────────────────────────────
    def _draw_loop(self, t: float) -> None:
        ax = self.ax_loop
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        L = self.loop
        knots = L.wave_path()
        lo, hi = max(self.t_past, self.t_emit), min(self.t_future, self.t_recv)
        if hi > lo:
            ax.axvspan(lo, hi, color=_PAL["both"], alpha=0.55, zorder=0)
            ax.annotate("both shells live", xy=(0.5 * (lo + hi), -0.05),
                        color=_PAL["dim"], fontsize=6.6, ha="center",
                        va="top", family="monospace")

        # the two propagation legs: null, so slope 1 in these coordinates
        ax.plot([knots[0]["t"], knots[1]["t"]], [knots[0]["s"], knots[1]["s"]],
                lw=2.2, color=_PAL["wave"], zorder=4)
        ax.plot([knots[2]["t"], knots[-1]["t"]], [knots[2]["s"], knots[-1]["s"]],
                lw=2.2, color=_PAL["wave"], zorder=4)
        # the throat: a jump backwards in TIME at fixed path length
        ax.annotate("", xy=(knots[2]["t"], knots[2]["s"]),
                    xytext=(knots[1]["t"], knots[1]["s"]),
                    arrowprops=dict(arrowstyle="-|>", lw=1.6, ls="--",
                                    color=_PAL["throat"]), zorder=5)
        offsets = ((-6, 5), (-5, -12), (5, -12), (5, -12), (9, -4))
        aligns = ("right", "right", "left", "left", "left")
        for k, col, off, ha in zip(knots, (_PAL["emitter"], _PAL["mouth_f"],
                                           _PAL["mouth_p"], _PAL["emitter"],
                                           _PAL["receiver"]), offsets, aligns):
            ax.plot([k["t"]], [k["s"]], "o", ms=5.5, color=col, zorder=6)
            ax.annotate(k["label"], xy=(k["t"], k["s"]), xytext=off,
                        textcoords="offset points", color=col, fontsize=6.6,
                        ha=ha, family="monospace", zorder=8)
        ax.axvline(t, color=_PAL["text"], lw=0.9, alpha=0.5, zorder=7)
        ax.axvline(0.0, color=_PAL["dim"], lw=0.7, ls=":", alpha=0.7)
        ax.set_xlim(self.t0, self.t1)
        ax.set_ylim(-0.75, 2.0 * math.pi + 0.5)
        ax.set_yticks([0, math.pi, 2 * math.pi])
        ax.set_yticklabels(["0", "π", "2π"])
        ax.set_xlabel("coordinate time", color=_PAL["dim"], fontsize=8)
        ax.set_ylabel("path length travelled by the wave", color=_PAL["dim"],
                      fontsize=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(_PAL["rule"])
        ax.grid(alpha=0.10, color=_PAL["grid"])
        ax.set_title("the circuit the earlier movie had no concept of — "
                     "path length never resets; the throat is the jump back",
                     color=_PAL["text"], fontsize=8.4, pad=6)

    # ── the books ───────────────────────────────────────────────────────────
    def _draw_book(self, t: float) -> None:
        ax = self.ax_book
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        L = self.loop
        a = L.self_consistent_amplitude()
        it = L.iterate_loop(2000)
        led = L.ledger()
        rows = [
            ("κ (loop transfer, a parameter)", f"{L.kappa.real:+.3f}"),
            ("A = A_src/(1−κ)", f"{a.real:.6f}"),
            ("same, by 2000 round trips", f"{it.real:.6f}"),
            ("|difference|", f"{abs(it - a):.1e}"),
            ("emitter→future mouth (must be π)",
             f"{L.emitter_to_future_mouth:.6f}"),
            ("past mouth→receiver (must be π)",
             f"{L.past_mouth_to_receiver:.6f}"),
            ("flux out of emitter", f"{abs(a) ** 2:.6f}"),
            ("flux into receiver", f"{abs(L.kappa * a) ** 2:.6f}"),
            ("ledger residual", f"{led['residual']:.1e}"),
        ]
        for i, (k, v) in enumerate(rows):
            y = 0.94 - i * 0.104
            ax.text(0.02, y, k, color=_PAL["dim"], fontsize=7.8,
                    family="monospace", va="center")
            ax.text(0.98, y, v, color=_PAL["text"], fontsize=7.8,
                    family="monospace", va="center", ha="right")
        ax.set_title("the ledger — one wave, counted once",
                     color=_PAL["text"], fontsize=8.4, pad=6)

    # ── frame ───────────────────────────────────────────────────────────────
    def draw(self, t: float) -> None:
        self._draw_s3(t)
        self._draw_loop(t)
        self._draw_book(t)
        self.fig.suptitle("v51 — ONE CONSERVED WAVE, SEEN IN PIECES",
                          color=_PAL["text"], fontsize=13.5, y=0.963,
                          family="monospace")
        self.fig.text(0.5, 0.914,
                      "emitter → antipodal mouth → throat (Δt < 0) → past "
                      "mouth → antipodal receiver   ·   both shells are drawn "
                      "in one colour because they are one wave",
                      color=_PAL["dim"], fontsize=8.6, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.040,
                      "derived: the antipodal refocus (4π sin²χ), used twice — "
                      "it is why the arriving shell is collapsing (dA/dχ < 0) "
                      "and not merely arriving   ·   put in: the throat is an "
                      "identification, and PR #249 priced it",
                      color=_PAL["dim"], fontsize=7.6, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.016,
                      "the shadow of S³ is 2-to-1 — brightness is the "
                      "discarded depth — so a crossing on screen is not a "
                      "crossing, and nothing here rests on one",
                      color="#3d5570", fontsize=7.0, ha="center",
                      family="monospace")


# ════════════════════════════════════════════════════════════════════════════
def still(path: str) -> None:
    fig = LedgerFigure()
    fig.draw(fig.t_recv - 0.30)
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def animate(save: Optional[str] = None, frames: int = 150):
    from matplotlib import animation

    holder = LedgerFigure()
    times = np.linspace(holder.t0, holder.t1, frames)

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
