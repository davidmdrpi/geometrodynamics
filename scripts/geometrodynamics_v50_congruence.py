#!/usr/bin/env python3
"""
Geometrodynamic QED — v50: focus the congruence to the pinch, and let it decide
===============================================================================

A singularity is a failure of evolution, not a bright spot. Drawing one as a
glowing dot assumes the answer. The object that does not is a **congruence with
a deforming cross-section**: integrate the geodesic-deviation equation
`F̈ = ½ḧF` from `F(0) = I` and watch the Jacobian `J = det F`, the
cross-sectional area of the bundle. `J → 0` is a caustic of the **map**. It
says nothing about the metric.

The bulk is built from the tidal field itself
─────────────────────────────────────────────
Every line on screen is a principal axis of `F`. For the axisymmetric field the
deformation gradient stays diagonal in the geodesic-polar frame, `F = diag(a,
b)`, so `a` is the principal stretch along `ê_d` — tangent to the slice — and
`b` the one along `ê_ψ`, which is normal to the slice on the sphere and is what
gets extruded into the bulk. The ellipse at each station is the congruence's
cross-section, and its area is `J`. No normals, no invented vector length: the
spin-2 tensor draws its own bulk structure.

What is on screen
─────────────────
* **the congruence** — the slice with a tidal cross-section at every station,
  collapsing as the bundle narrows. The ring underneath is coloured by
  curvature strength `|ḧ|`, the same geometry carrying the other quantity.
* **the spacetime history** — `J(σ, t)`, one spatial cross-section vertically
  and time horizontally, with the `J = 0` locus picked out. This is where a
  topology change would have to show, and the causal bound `t = d − w` is drawn
  on it: nothing moves before the front.
* **the neck** — `√|J|`, the linear scale of the cross-section, for the source
  ring and the converging ring separately.

What the equations decided
──────────────────────────
Of the three outcomes on offer — passage, singular termination, finite-radius
reconnection — this program gives **passage**. At the source ring `J` crosses
zero with slope `−17.877`, converged to `−17.836` under a halved timestep,
plunges to `−471` and stays; a tangency would have driven that slope to zero.
The solver's invariant is unmoved at `2.5e-14`.

The other two were never available, for different reasons worth keeping apart.
Termination needs the geometry to fail, and the background here is a fixed
round `S²` with curvature `1` at every time. Reconnection needs the congruence
to act back on something, and each point's `F` is driven only by the external
`h`, never by its neighbours. So this is not "we looked and did not find them";
it is "this program could not have produced them".

Two rings, and they are not alike
─────────────────────────────────
The source ring closes at peak strain `0.026`, the converging ring at `0.247` —
a factor of ten. Even at `0.25` the antipodal crossing only **grazes** zero, and
the depth of that excursion does not converge. And the neck is never at the
antipode: `h = sin²d·q` vanishes at both poles by spin weight, so the tidal
field is identically zero there and the congruence is never driven. The neck
sits on a ring of radius `0.44 w`, the same ratio across a `3.3×` range of pulse
width. A spin-2 focus has no centre.

Honest scope
────────────
`gain` is a strength dial and is reported as a peak strain. The deviation
equation is exact in `ξ` and linear in `h`; at the strain where the converging
ring closes the field is no longer a weak perturbation, and that is stated
rather than hidden. The drawn cross-sections are **sparse** and their sizes are
clipped for legibility — collapse is resolved, growth saturates — so read the
areas as the story and the panels for the numbers.

Usage
─────
    python scripts/geometrodynamics_v50_congruence.py             # animate
    python scripts/geometrodynamics_v50_congruence.py --still out.png
    python scripts/geometrodynamics_v50_congruence.py --cases out.png
    python scripts/geometrodynamics_v50_congruence.py --save out.gif
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from geometrodynamics.viz.congruence import (
    ANTIPODAL_TIME,
    NECK_CAP,
    TidalCongruence,
)

FPS = 14
GAIN = 60.0
PULSE_WIDTH = 0.18
N_SIGMA = 361
N_RADIAL = 1600
T_END = 1.25 * ANTIPODAL_TIME
STRIDE = 10                    # sparse enough to read as cross-sections
TANGENTIAL = 0.46              # the ê_d axis lies along the slice, so it
                               # is drawn shorter to keep stations apart
R_MID = 1.0
BULK = 0.15                    # radial half-extent given to the ê_ψ axis

_PAL = {
    "bg": "#05070c",
    "panel": "#080b13",
    "text": "#e8ecf4",
    "dim": "#7d8798",
    "faint": "#2a3244",
    "shell": "#39445c",
    "ring": "#ffd166",
    "axis_d": "#7cc4ff",
    "axis_psi": "#ff9f43",
    "pinch": "#ff3b6b",
    "neck": "#c08cff",
    "source": "#8ef0c0",
    "antipode": "#ff8ad0",
}


# ════════════════════════════════════════════════════════════════════════════
def precompute(gain: float = GAIN, frames: int = 150,
               t_end: float = T_END) -> dict:
    """Integrate once and keep the frames; the panels all read from this."""
    c = TidalCongruence(gain=gain, n_sigma=N_SIGMA, n_radial=N_RADIAL,
                        pulse_width=PULSE_WIDTH)
    step = max(1, int((t_end / c.dt) / frames))
    ts, A, B, H, HD = [], [], [], [], []
    k = 0
    while c.t < t_end:
        c.step()
        k += 1
        if k % step == 0:
            ts.append(c.t)
            A.append(c.a.copy())
            B.append(c.b.copy())
            H.append(c.h().copy())
            HD.append(c.h_ddot().copy())
    A, B = np.array(A), np.array(B)
    return {
        "sigma": c.sigma, "distance": c.distance, "t": np.array(ts),
        "a": A, "b": B, "J": A * B, "h": np.array(H), "hdd": np.array(HD),
        "causal_bound": c.causal_bound(c.distance),
        "far": c.far_mask(), "cap": (math.pi - c.distance) < NECK_CAP,
        "gain": gain, "peak_strain": float(np.max(np.abs(np.array(H)))),
        "pulse_width": PULSE_WIDTH,
    }


def _neck_scale(J: np.ndarray) -> np.ndarray:
    """``log₁₀ √|J|`` clipped to a readable band."""
    return np.clip(np.log10(np.sqrt(np.abs(J)) + 1e-6), -3.0, 0.6)


def _compress(v: np.ndarray) -> np.ndarray:
    """``|v|/(1+|v|)``, normalised so an undeformed axis draws at ``BULK``.

    The principal stretches range over decades in both directions, so a hard
    clip turns every station into the same saturated circle and destroys the
    only thing the panel is for.  This map is monotone, bounded, and linear for
    small ``|v|`` — collapse stays proportional, growth saturates gently.
    """
    x = np.abs(v)
    return 2.0 * BULK * x / (1.0 + x)


class CongruenceFigure:
    """The congruence, its worldsheet, and the neck."""

    def __init__(self, data: dict, figsize=(13.6, 8.4)) -> None:
        self.d = data
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            2, 2, width_ratios=[1.15, 1.0], height_ratios=[1.0, 0.78],
            left=0.045, right=0.965, top=0.855, bottom=0.095,
            wspace=0.20, hspace=0.32)
        self.ax_geo = self.fig.add_subplot(gs[:, 0], facecolor=_PAL["bg"])
        self.ax_hist = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_neck = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])

    # ── the embedded geometry ───────────────────────────────────────────────
    def _draw_geometry(self, i: int) -> None:
        ax = self.ax_geo
        ax.clear()
        ax.set_facecolor(_PAL["bg"])
        d = self.d
        sig, J = d["sigma"], d["J"][i]
        a, b, hdd = d["a"][i], d["b"][i], d["hdd"][i]

        # the slice itself, coloured by curvature strength |ḧ|
        n = np.stack([np.cos(sig), np.sin(sig)], axis=1)
        P = R_MID * n
        segs = np.stack([P, np.roll(P, -1, axis=0)], axis=1)
        strength = np.abs(hdd)
        norm = strength / (np.percentile(strength, 99.0) + 1e-12)
        lc = LineCollection(segs, cmap="inferno", linewidths=2.6,
                            array=np.clip(norm, 0, 1), zorder=2)
        ax.add_collection(lc)

        # the tidal cross-sections: every line is a principal axis of F
        idx = np.arange(0, len(sig), STRIDE)
        tang = np.stack([-np.sin(sig), np.cos(sig)], axis=1)
        sa = TANGENTIAL * _compress(a)
        sb = _compress(b)
        scale = _neck_scale(J)
        cmap = plt.get_cmap("magma")
        # bright where the cross-section is collapsing, and never so dark that
        # an ordinary station vanishes into the background
        cnorm = 0.30 + 0.70 * (1.0 - (scale + 3.0) / 3.6)

        self._stations(ax, idx, P, tang, n, sa, sb, J, cmap, cnorm)

        pinched = int(np.sum(J <= 0.0))
        # σ = 0 is the source and σ = ±π the antipode, so they sit on the x axis
        ax.plot([R_MID], [0], marker="o", ms=5, color=_PAL["source"], zorder=6)
        ax.plot([-R_MID], [0], marker="o", ms=5, color=_PAL["antipode"],
                zorder=6)

        # The antipodal cap, where the whole question lives — and it needs its
        # own, finer sampling: the neck ring sits at 0.44 w from the pole, which
        # the sparse stations of the main ring straddle and miss completely.
        inset = ax.inset_axes([0.775, 0.755, 0.225, 0.225],
                              facecolor=_PAL["panel"])
        fine = np.arange(0, len(sig), 2)
        near_pole = fine[(math.pi - d["distance"][fine]) < 0.30]
        self._stations(inset, near_pole, P, tang, n, sa, sb, J, cmap, cnorm,
                       lw=1.0)
        inset.plot([-R_MID], [0], marker="o", ms=4, color=_PAL["antipode"],
                   zorder=6)
        inset.set_xlim(-R_MID - BULK - 0.05, -R_MID + BULK + 0.05)
        inset.set_ylim(-0.30, 0.30)
        inset.set_aspect("equal")
        inset.set_xticks([])
        inset.set_yticks([])
        for sp in inset.spines.values():
            sp.set_color(_PAL["faint"])
        inset.set_title("the antipodal cap", color=_PAL["dim"],
                        fontsize=8.0, pad=3)
        inset.set_xlabel("the neck is a ring; the pole is never driven",
                         color=_PAL["dim"], fontsize=6.2, labelpad=2)

        lim = R_MID + BULK + 0.16
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(
            f"the congruence — every line is a principal axis of F\n"
            f"t = {d['t'][i]:.3f}   min J = {J.min():+.4f}   "
            f"inverted stations: {pinched}",
            color=_PAL["text"], fontsize=10.5, pad=9)

    def _stations(self, ax, idx, P, tang, n, sa, sb, J, cmap, cnorm,
                  lw: float = 1.5) -> None:
        """Draw the principal cross and the cross-section ellipse at each idx."""
        d_segs, psi_segs, cols_d, cols_p = [], [], [], []
        theta = np.linspace(0.0, 2.0 * math.pi, 41)
        ell_segs, ell_cols = [], []
        for k in idx:
            p, tv, nv = P[k], tang[k], n[k]
            la, lb = sa[k], sb[k]
            d_segs.append([p - la * tv, p + la * tv])
            psi_segs.append([p - lb * nv, p + lb * nv])
            inv = J[k] <= 0.0
            cols_d.append(_PAL["pinch"] if inv else _PAL["axis_d"])
            cols_p.append(_PAL["pinch"] if inv else _PAL["axis_psi"])
            e = (p[None, :] + la * np.cos(theta)[:, None] * tv[None, :]
                 + lb * np.sin(theta)[:, None] * nv[None, :])
            ell_segs.append(np.stack([e[:-1], e[1:]], axis=1))
            ell_cols.append(cmap(np.clip(cnorm[k], 0, 1)))

        for e, c in zip(ell_segs, ell_cols):
            ax.add_collection(LineCollection(e, colors=[c], linewidths=0.9,
                                             alpha=0.9, zorder=3))
        ax.add_collection(LineCollection(psi_segs, colors=cols_p,
                                         linewidths=lw, zorder=4))
        ax.add_collection(LineCollection(d_segs, colors=cols_d,
                                         linewidths=lw, zorder=4))

    # ── the spacetime history ───────────────────────────────────────────────
    def _draw_history(self, i: int) -> None:
        ax = self.ax_hist
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        d = self.d
        J = d["J"][: i + 1].T                      # (σ, t)
        t = d["t"][: i + 1]
        # the neck radius on a log scale: J spans many decades and a diverging
        # map centred on zero saturates almost everywhere
        ax.pcolormesh(t, d["sigma"], _neck_scale(J), cmap="magma_r",
                      vmin=-3.0, vmax=0.6, shading="auto", rasterized=True)
        # where the map has turned inside out, drawn as its own layer
        inverted = np.ma.masked_where(J > 0.0, np.ones_like(J))
        ax.pcolormesh(t, d["sigma"], inverted,
                      cmap=matplotlib.colors.ListedColormap([_PAL["pinch"]]),
                      vmin=0.0, vmax=1.0, shading="auto", alpha=0.32,
                      rasterized=True)
        # Nothing may move before the front.  The bound t = |σ| − w has two
        # branches meeting at σ = 0, so it has to be drawn in σ order: sorting
        # by t instead makes one polyline zig-zag between ±σ and fill the panel
        # with a triangular hatch that looks exactly like a light cone.
        ax.plot(d["causal_bound"], d["sigma"], lw=1.1, ls="--",
                color=_PAL["text"], alpha=0.6)
        ax.set_xlim(d["t"][0], d["t"][-1])
        ax.set_ylim(-math.pi, math.pi)
        ax.set_yticks([-math.pi, -math.pi / 2, 0, math.pi / 2, math.pi])
        ax.set_yticklabels(["−π", "−π/2", "0", "π/2", "π"])
        ax.set_xlabel("t", color=_PAL["dim"], fontsize=9)
        ax.set_ylabel("σ", color=_PAL["dim"], fontsize=9)
        ax.tick_params(colors=_PAL["dim"], labelsize=8)
        for s in ax.spines.values():
            s.set_color(_PAL["faint"])
        ax.set_title("spacetime history — neck radius √|J|(σ, t); red is "
                     "J < 0, dashed is t = d − w",
                     color=_PAL["text"], fontsize=9.5, pad=6)

    # ── the neck ────────────────────────────────────────────────────────────
    def _draw_neck(self, i: int) -> None:
        ax = self.ax_neck
        ax.clear()
        ax.set_facecolor(_PAL["panel"])
        d = self.d
        t = d["t"][: i + 1]
        J = d["J"][: i + 1]
        near = np.sqrt(np.abs(J[:, ~d["far"]]).min(axis=1))
        cap = np.sqrt(np.abs(J[:, d["cap"]]).min(axis=1))
        ax.semilogy(t, np.maximum(near, 1e-6), lw=1.6, color=_PAL["source"],
                    label="source ring")
        ax.semilogy(t, np.maximum(cap, 1e-6), lw=1.6, color=_PAL["antipode"],
                    label="converging ring (antipodal cap)")
        ax.axvline(ANTIPODAL_TIME, color=_PAL["dim"], lw=0.9, ls=":")
        ax.set_xlim(d["t"][0], d["t"][-1])
        ax.set_ylim(1e-4, 3.0)
        ax.set_xlabel("t", color=_PAL["dim"], fontsize=9)
        ax.set_ylabel("neck radius  √|J|", color=_PAL["dim"], fontsize=9)
        ax.tick_params(colors=_PAL["dim"], labelsize=8)
        ax.grid(alpha=0.12, color=_PAL["shell"])
        for s in ax.spines.values():
            s.set_color(_PAL["faint"])
        leg = ax.legend(loc="lower left", fontsize=7.6, framealpha=0.0)
        for txt in leg.get_texts():
            txt.set_color(_PAL["dim"])
        ax.set_title("the neck — the linear scale of the cross-section",
                     color=_PAL["text"], fontsize=9.5, pad=6)

    # ── frame ───────────────────────────────────────────────────────────────
    def draw(self, i: int) -> None:
        self._draw_geometry(i)
        self._draw_history(i)
        self._draw_neck(i)
        d = self.d
        self.fig.suptitle(
            "v50 — focus the congruence to the pinch threshold, and let the "
            "equations decide",
            color=_PAL["text"], fontsize=13.5, y=0.965)
        self.fig.text(
            0.5, 0.915,
            f"J = det F from F̈ = ½ḧF   ·   peak strain {d['peak_strain']:.3f}"
            f"   ·   J → 0 is a caustic of the MAP, not of the metric",
            color=_PAL["dim"], fontsize=9.2, ha="center")
        self.fig.text(
            0.5, 0.028,
            "the caustic is a passage: J crosses zero transversally and the "
            "evolution continues   ·   the neck is a ring of radius 0.44 w, "
            "never the antipode — spin weight forces h to vanish at the pole",
            color=_PAL["dim"], fontsize=8.6, ha="center")


# ════════════════════════════════════════════════════════════════════════════
def still(path: str) -> None:
    """The frame where the converging ring's neck is thinnest."""
    data = precompute()
    cap = np.sqrt(np.abs(data["J"][:, data["cap"]]).min(axis=1))
    fig = CongruenceFigure(data)
    fig.draw(int(np.argmin(cap)))
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def cases(path: str) -> None:
    """The three cases side by side, at the amplitudes that separate them."""
    gains = [(4.0, "ordinary focus"), (30.0, "source ring closes"),
             (150.0, "both rings close")]
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 5.0),
                             facecolor=_PAL["bg"])
    for ax, (g, label) in zip(axes, gains):
        d = precompute(gain=g, frames=110)
        ax.set_facecolor(_PAL["panel"])
        t, J = d["t"], d["J"]
        near = np.sqrt(np.abs(J[:, ~d["far"]]).min(axis=1))
        cap = np.sqrt(np.abs(J[:, d["cap"]]).min(axis=1))
        ax.semilogy(t, np.maximum(near, 1e-6), lw=1.7, color=_PAL["source"],
                    label="source ring")
        ax.semilogy(t, np.maximum(cap, 1e-6), lw=1.7, color=_PAL["antipode"],
                    label="converging ring")
        ax.axvline(ANTIPODAL_TIME, color=_PAL["dim"], lw=0.9, ls=":")
        ax.set_ylim(1e-4, 3.0)
        ax.set_xlabel("t", color=_PAL["dim"], fontsize=9)
        ax.set_title(f"{label}\npeak strain {d['peak_strain']:.3f}",
                     color=_PAL["text"], fontsize=10, pad=8)
        ax.tick_params(colors=_PAL["dim"], labelsize=8)
        ax.grid(alpha=0.12, color=_PAL["shell"])
        for s in ax.spines.values():
            s.set_color(_PAL["faint"])
    axes[0].set_ylabel("neck radius  √|J|", color=_PAL["dim"], fontsize=9)
    leg = axes[0].legend(loc="lower left", fontsize=7.6, framealpha=0.0)
    for txt in leg.get_texts():
        txt.set_color(_PAL["dim"])
    fig.suptitle("the three cases — a dip that recovers, one ring closing, "
                 "then both", color=_PAL["text"], fontsize=13, y=0.98)
    fig.text(0.5, 0.02,
             "case 3, a curvature singularity, never appears and cannot: the "
             "background is a fixed round S² with curvature 1 at every time",
             color=_PAL["dim"], fontsize=8.8, ha="center")
    fig.subplots_adjust(left=0.06, right=0.98, top=0.80, bottom=0.16,
                        wspace=0.22)
    fig.savefig(path, dpi=120, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def animate(save: Optional[str] = None, frames: int = 150):
    from matplotlib import animation

    data = precompute(frames=frames)
    holder = CongruenceFigure(data)
    n = len(data["t"])

    def update(i: int):
        holder.draw(min(i, n - 1))
        return []

    anim = animation.FuncAnimation(holder.fig, update, frames=n,
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
    ap.add_argument("--cases", metavar="PNG")
    ap.add_argument("--save", metavar="FILE")
    ap.add_argument("--frames", type=int, default=150)
    a = ap.parse_args(argv)
    if a.still or a.cases:
        matplotlib.use("Agg")
        if a.still:
            still(a.still)
        if a.cases:
            cases(a.cases)
        return 0
    if a.save:
        matplotlib.use("Agg")
    animate(save=a.save, frames=a.frames)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
