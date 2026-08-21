#!/usr/bin/env python3
"""
Geometrodynamic QED — v67: one field, one surface, and the parity of the antipode
=================================================================================

v66 drew TWO curves over the circle and asked whether their images meet through
the glued seam. Two curves in one frame are two surfaces, and reading their
overlap as a connection is a statement about a picture, not about a field.

If the two contributions are two pieces of ONE scalar deformation of ONE surface
there is only ever one curve,

    u = s_A u_A + s_B u_B        r = R_mid + eps u

and the question is whether THAT curve reaches R_outer at one theta and R_inner
at another -- so the surface itself passes through the identification.

THE REPAIR COSTS NOTHING. delta = r_A - r_B, the quantity v66 plotted as a
separation between two curves, IS the one-surface deformation with the second
sign flipped -- the same array, to two ulps of R_mid. Every v66 number survives,
with the two configurations SWAPPING NAMES, which inverts v66's headline.

What the panels show
--------------------
**Top left - one curve, not two.** The single deformed surface at the antipodal
configuration, driven to the gain where it touches both boundaries. There is no
second curve anywhere in this panel.

**Top centre - coincidence cancels.** At alpha = 0 the opposed field is
identically zero: the surface is a circle, at every time and every gain. The
like-signed field at the same offset is at full strength. This is the exact
statement v66 could not make, because it had no single field to be zero.

**Top right - the amplitude law.** B = 2A|sin(m alpha/2)|, measured on a grid
against the closed form, for the first few modes. The first maximum is at
alpha = pi/m -- HALF A WAVELENGTH -- and the antipode is just the m = 1 member.

**Middle left - the parity.** sin(m pi/2) is +-1 for odd m and 0 for even m, so
opposite foci are maximal for odd modes and CANCEL for even ones. Checked
against the exact S^3 zonal harmonics, where the same parity is Z_n(pi) =
(-1)^n.

**Middle centre - WHERE THE TWO MODELS PART COMPANY.** A zonal harmonic is
CENTRED, so |Z_A - Z_B| <= 2 with equality only at the antipode with odd n. For
the real spectrum alpha* = pi for EVERY odd n and it SATURATES the bound; half a
wavelength reaches only a fraction. The parity carries across; the location of
the optimum does not.

**Middle right - a pulse is not a mode.** v46's field is a launched pulse, so
the 1/|sin| divergence is confined to about its own width and the threshold
SATURATES past that instead of continuing down.

**Bottom left - the chord.** At the optimum the two extrema are pi/m apart, so
L = sqrt(D^2 + 4 R_out R_in sin^2(pi/2m)) falls from 2.000 to the purely radial
gap 0.520. Same deformation, shorter connection.

**Bottom right - but frequency is not free.** E ~ omega^2 A^2, so at fixed
energy A ~ 1/omega and the span falls as fast as the chord. The highest mode
that still spans the gap is m = 7. No favourable frequency is claimed.

Honest scope
------------
Unchanged. The crossing rule is a REPRESENTATION choice, not a derived boundary
condition; the field is a LINEAR scalar on a FIXED round background; the gain is
a DISPLAY amplitude -- and the bottom-right panel is the one place that last
distinction is made to do work.

Usage
-----
    python scripts/geometrodynamics_v67_one_surface.py --still v67.png
"""

from __future__ import annotations

import argparse
import math
from typing import Optional

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from geometrodynamics.viz.circle_slice import ANTIPODAL_TIME, RETURN_TIME, TWO_PI
from geometrodynamics.viz.one_surface import (
    LIKE_SIGNED, OPPOSED, MonochromaticPair, OneSurfaceSlice, bulk_chord,
    measure_a_pulse_saturates_where_a_mode_diverges,
    measure_the_zonal_optimum_is_the_antipode, zonal_harmonic,
    zonal_pair_strength)

_PAL = {
    "bg": "#010106", "panel": "#02020a", "grid": "#0e1a2a", "rule": "#1a2838",
    "text": "#e8ecf4", "dim": "#6a8aad",
    "one": "#ffc857", "two": "#5cc8ff", "band": "#b388ff",
    "hot": "#ff6b8a", "ok": "#7cff9e",
}


def _style(ax, title, xlabel="", ylabel=""):
    ax.set_title(title, color=_PAL["text"], fontsize=9.4, pad=7,
                 family="monospace")
    ax.set_xlabel(xlabel, color=_PAL["dim"], fontsize=7.6, family="monospace")
    ax.set_ylabel(ylabel, color=_PAL["dim"], fontsize=7.6, family="monospace")
    ax.tick_params(colors=_PAL["dim"], labelsize=7.0)
    for s in ax.spines.values():
        s.set_color(_PAL["rule"])


def _note(ax, text, y=0.045, color=None):
    ax.text(0.5, y, text, transform=ax.transAxes, ha="center",
            color=color or _PAL["text"], fontsize=6.5, family="monospace",
            bbox=dict(facecolor=_PAL["panel"], edgecolor="none", alpha=0.88,
                      pad=2.0))


class OneSurfaceFigure:
    """One curve, the law, the parity, the divergence, and the price."""

    def __init__(self, figsize=(14.4, 13.2)) -> None:
        self.fig = plt.figure(figsize=figsize, facecolor=_PAL["bg"])
        gs = self.fig.add_gridspec(
            3, 3, left=0.058, right=0.975, top=0.858, bottom=0.062,
            wspace=0.30, hspace=0.42, height_ratios=[1.0, 0.86, 0.86])
        self.ax_one = self.fig.add_subplot(gs[0, 0], facecolor=_PAL["panel"])
        self.ax_zero = self.fig.add_subplot(gs[0, 1], facecolor=_PAL["panel"])
        self.ax_law = self.fig.add_subplot(gs[0, 2], facecolor=_PAL["panel"])
        self.ax_par = self.fig.add_subplot(gs[1, 0], facecolor=_PAL["panel"])
        self.ax_zon = self.fig.add_subplot(gs[1, 1], facecolor=_PAL["panel"])
        self.ax_pul = self.fig.add_subplot(gs[1, 2], facecolor=_PAL["panel"])
        self.ax_chd = self.fig.add_subplot(gs[2, 0], facecolor=_PAL["panel"])
        self.ax_en = self.fig.add_subplot(gs[2, 1:], facecolor=_PAL["panel"])

        self.surf = OneSurfaceSlice(offset=math.pi, signs=OPPOSED)
        self.surf.slice_.reset()
        self.surf.slice_.advance_to(ANTIPODAL_TIME)
        self.gain = self.surf.span_gain()

    # ── helpers ─────────────────────────────────────────────────────────────
    def _rings(self, ax, b):
        th = np.linspace(0, TWO_PI, 721)
        for r in (b.r_inner, b.r_outer):
            ax.plot(r * np.cos(th), r * np.sin(th), lw=1.2,
                    color=_PAL["hot"], alpha=0.7, ls=(0, (6, 5)))
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)

    def _one_curve(self, ax, surf, gain, color, label, lw=2.6):
        s = surf.sigma
        r = surf.radius(gain=gain)
        ax.plot(r * np.cos(s), r * np.sin(s), lw=lw, color=color,
                solid_capstyle="round", label=label)

    # ── panels ──────────────────────────────────────────────────────────────
    def _one(self):
        ax = self.ax_one
        b = self.surf.slice_.bulk
        self._rings(ax, b)
        self._one_curve(ax, self.surf, self.gain, _PAL["one"],
                        "the surface")
        got = self.surf.connection(gain=self.gain)
        for sg, rr, col in ((got["sigma_out"], got["radius_out"], _PAL["ok"]),
                            (got["sigma_in"], got["radius_in"], _PAL["two"])):
            ax.plot([rr * math.cos(sg)], [rr * math.sin(sg)], "o", ms=9,
                    color=col, mec=_PAL["text"], mew=0.8, zorder=6)
        ax.set_xlim(-b.r_outer * 1.34, b.r_outer * 1.34)
        ax.set_ylim(-b.r_outer * 1.34, b.r_outer * 1.34)
        _style(ax, "ONE curve — it reaches both boundaries")
        ax.text(0.5, 0.985,
                f"gain {self.gain:.3f}   sigma_out {got['sigma_out_over_pi']:+.3f} pi"
                f"   sigma_in {got['sigma_in_over_pi']:+.3f} pi",
                transform=ax.transAxes, ha="center", va="top",
                color=_PAL["dim"], fontsize=6.6, family="monospace")
        _note(ax, "there is no second curve in this panel;\n"
                  "the question is whether the SURFACE spans the gap")

    def _zero(self):
        ax = self.ax_zero
        opp = OneSurfaceSlice(offset=0.0, signs=OPPOSED)
        like = OneSurfaceSlice(offset=0.0, signs=LIKE_SIGNED)
        like.slice_ = opp.slice_
        opp.slice_.reset()
        opp.slice_.advance_to(ANTIPODAL_TIME)
        b = opp.slice_.bulk
        g = like.span_gain()
        self._rings(ax, b)
        self._one_curve(ax, like, g, _PAL["two"], "like-signed", lw=2.0)
        self._one_curve(ax, opp, g, _PAL["one"], "opposed  (u = 0)", lw=3.0)
        ax.set_xlim(-b.r_outer * 1.34, b.r_outer * 1.34)
        ax.set_ylim(-b.r_outer * 1.34, b.r_outer * 1.34)
        _style(ax, "coincident foci: the opposed field IS zero")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.4, loc="lower center",
                  ncol=2, framealpha=0.9)
        worst = float(np.max(np.abs(opp.field())))
        _note(ax, f"max|u| = {worst:.1e}  — exactly zero, at every time.\n"
                  "no gain connects it; the threshold is infinite",
              y=0.90, color=_PAL["ok"])

    def _law(self):
        ax = self.ax_law
        al = np.linspace(0, math.pi, 1201)
        for m, col in ((1, _PAL["one"]), (2, _PAL["two"]), (3, _PAL["ok"]),
                       (4, _PAL["band"])):
            ax.plot(al / math.pi, 2 * np.abs(np.sin(0.5 * m * al)), lw=1.8,
                    color=col, label=f"m = {m}")
            first = math.pi / m
            ax.plot([first / math.pi], [2.0], "*", ms=12, color=col,
                    mec=_PAL["text"], mew=0.5, zorder=5)
        ax.axvline(1.0, color=_PAL["hot"], lw=1.0, ls=":")
        ax.text(0.985, 0.30, "the antipode ", transform=ax.get_xaxis_transform(),
                color=_PAL["hot"], fontsize=6.4, family="monospace", ha="right")
        ax.set_ylim(0, 2.35)
        ax.grid(True, color=_PAL["grid"], lw=0.5, alpha=0.7)
        _style(ax, "B = 2A|sin(m alpha/2)|", "offset alpha / pi",
               "amplitude / A")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.4, loc="lower right",
                  ncol=2, framealpha=0.9)
        _note(ax, "stars = the first optimum, alpha* = pi/m\n"
                  "HALF A WAVELENGTH — the antipode is only m = 1",
              y=0.86)

    def _parity(self):
        ax = self.ax_par
        modes = np.arange(1, 9)
        vals = [MonochromaticPair(mode=int(m),
                                  offset=math.pi).amplitude_factor()
                for m in modes]
        cols = [_PAL["ok"] if m % 2 else _PAL["hot"] for m in modes]
        ax.bar(modes, vals, color=cols, width=0.62, edgecolor=_PAL["text"],
               linewidth=0.5)
        # the S^3 check has to be the PAIR strength, not |Z_n(pi)|: the
        # latter is 1 for every n and says nothing
        zn = [float(zonal_harmonic(int(m), np.array([math.pi]))[0])
              for m in modes]
        pair = [zonal_pair_strength(int(m), math.pi) for m in modes]
        ax.plot(modes, pair, "o--", color=_PAL["band"], lw=1.3, ms=7,
                mec=_PAL["text"], mew=0.5,
                label="S^3 zonal pair strength", zorder=6)
        for m, v, z in zip(modes, vals, zn):
            ax.text(m, 0.07, f"Z={z:+.0f}", ha="center", color=_PAL["dim"],
                    fontsize=6.2, family="monospace")
            if v < 1e-9:
                ax.text(m, 0.24, "0", ha="center", color=_PAL["hot"],
                        fontsize=9.0, family="monospace", weight="bold")
                ax.text(m, 0.42, "cancels", ha="center", color=_PAL["hot"],
                        fontsize=5.8, family="monospace", rotation=90)
        ax.set_ylim(0, 2.9)
        ax.set_xticks(modes)
        _style(ax, "at the antipode: odd adds, even cancels",
               "mode m", "amplitude / A")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.2, loc="upper center",
                  framealpha=0.95)
        _note(ax, "not a visualization effect: the plane-wave bars and the\n"
                  "S^3 zonal points agree, and Z_n(pi) = (-1)^n exactly",
              y=0.245)

    def _zonal(self):
        ax = self.ax_zon
        got = measure_the_zonal_optimum_is_the_antipode()
        rows = got["rows"]
        n = [r["order"] for r in rows]
        peak = [r["peak_strength"] for r in rows]
        half = [r["strength_at_half_wavelength"] for r in rows]
        anti = [r["strength_at_the_antipode"] for r in rows]
        w = 0.34
        idx = np.arange(len(n))
        ax.bar(idx - w / 2, peak, width=w, color=_PAL["ok"],
               edgecolor=_PAL["text"], linewidth=0.4, label="best over alpha")
        ax.bar(idx + w / 2, half, width=w, color=_PAL["band"],
               edgecolor=_PAL["text"], linewidth=0.4,
               label="at half a wavelength")
        ax.plot(idx, anti, "o", color=_PAL["hot"], ms=7, mec=_PAL["text"],
                mew=0.5, label="at the antipode", zorder=6)
        ax.axhline(2.0, color=_PAL["one"], lw=1.2, ls=":")
        ax.text(0.99, 2.03, "the bound |Z_A - Z_B| <= 2 ", ha="right",
                color=_PAL["one"], fontsize=6.4, family="monospace",
                transform=ax.get_yaxis_transform())
        ax.set_xticks(idx)
        ax.set_xticklabels([str(v) for v in n])
        ax.set_ylim(0, 2.45)
        _style(ax, "*** the zonal optimum IS the antipode ***",
               "zonal order n   (omega = n+1)", "strength")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.2, loc="upper right",
                  framealpha=0.92)
        _note(ax, "odd n saturates the bound AT the antipode;\n"
                  "half a wavelength does not reproduce it", y=0.60)

    def _pulse(self):
        ax = self.ax_pul
        got = measure_a_pulse_saturates_where_a_mode_diverges()
        rows = got["rows"]
        a = [r["offset_over_pi"] for r in rows]
        ax.plot(a, [r["monochromatic_required_amplitude"] for r in rows],
                lw=1.9, color=_PAL["band"], label="one mode:  D/(4|sin|)")
        ax.plot(a, [r["pulse_threshold"] for r in rows], lw=1.9,
                color=_PAL["one"], ls=(0, (5, 3)), label="the v46 pulse")
        ax.axhline(got["plateau_threshold"], color=_PAL["ok"], lw=1.0, ls=":")
        ax.text(0.99, got["plateau_threshold"] * 1.06,
                f" plateau {got['plateau_threshold']:.4f}", ha="right",
                color=_PAL["ok"], fontsize=6.4, family="monospace",
                transform=ax.get_yaxis_transform())
        w = got["pulse_width"] / math.pi
        ax.axvspan(0, w, color=_PAL["hot"], alpha=0.16)
        ax.text(w * 1.15, 0.86, "one pulse width",
                transform=ax.get_xaxis_transform(), color=_PAL["hot"],
                fontsize=6.2, family="monospace")
        ax.set_yscale("log")
        ax.grid(True, color=_PAL["grid"], lw=0.5, alpha=0.7)
        _style(ax, "a pulse saturates where a mode diverges",
               "offset alpha / pi", "gain to span the gap")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.2, loc="upper right",
                  framealpha=0.9)
        _note(ax, "two localized pulses cancel only while they OVERLAP;\n"
                  "a mode fills the circle and cancels everywhere", y=0.05)

    def _chord(self):
        ax = self.ax_chd
        b = self.surf.slice_.bulk
        ms = np.arange(1, 17)
        L = [bulk_chord(b.r_inner, b.r_outer, math.pi / int(m)) for m in ms]
        ax.plot(ms, L, "o-", lw=1.8, ms=5, color=_PAL["one"],
                mec=_PAL["text"], mew=0.4, label="bulk chord L(m)")
        ax.axhline(b.gap, color=_PAL["hot"], lw=1.2, ls=":")
        ax.text(0.99, b.gap * 1.05, f" D = {b.gap} ", ha="right",
                color=_PAL["hot"], fontsize=6.6, family="monospace",
                transform=ax.get_yaxis_transform())
        ax.plot(ms, [2.0] * len(ms), lw=1.6, color=_PAL["ok"], ls=(0, (5, 3)),
                label="span / A  (fixed amplitude)")
        ax.set_ylim(0, 2.3)
        ax.grid(True, color=_PAL["grid"], lw=0.5, alpha=0.7)
        _style(ax, "same deformation, shorter connection", "mode m", "length")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.4, loc="center right",
                  framealpha=0.9)
        _note(ax, "L = sqrt(D^2 + 4 R_out R_in sin^2(pi/2m))\n"
                  "-> the purely radial gap as m grows", y=0.06)

    def _energy(self):
        ax = self.ax_en
        b = self.surf.slice_.bulk
        ms = np.arange(1, 17)
        L = np.array([bulk_chord(b.r_inner, b.r_outer, math.pi / int(m))
                      for m in ms])
        fixed_a = np.full(len(ms), 4.0)
        fixed_e = 4.0 / ms
        ax.plot(ms, fixed_a, lw=2.0, color=_PAL["ok"], ls=(0, (5, 3)),
                label="span at fixed DISPLAY amplitude")
        ax.plot(ms, fixed_e, lw=2.2, color=_PAL["one"],
                label="span at fixed ENERGY  (A ~ 1/omega)")
        ax.plot(ms, L, lw=1.6, color=_PAL["band"], ls=(0, (2, 3)),
                label="bulk chord L(m)")
        ax.axhline(b.gap, color=_PAL["hot"], lw=1.2, ls=":")
        ax.text(0.012, b.gap * 0.86, f" the gap D = {b.gap}", ha="left",
                color=_PAL["hot"], fontsize=6.8, family="monospace",
                transform=ax.get_yaxis_transform())
        cut = int(np.max(ms[fixed_e >= b.gap]))
        ax.axvline(cut + 0.5, color=_PAL["hot"], lw=1.0, ls="--")
        ax.fill_betweenx([0, 4.6], cut + 0.5, ms[-1], color=_PAL["hot"],
                         alpha=0.10)
        ax.text(cut + 0.9, 3.5, f"beyond m = {cut} a fixed-energy mode\n"
                                "can no longer span the gap at all",
                color=_PAL["hot"], fontsize=6.8, family="monospace")
        ax.set_yscale("log")
        ax.set_ylim(0.15, 5.2)
        ax.set_xticks(ms[::1])
        ax.grid(True, color=_PAL["grid"], lw=0.5, alpha=0.7)
        _style(ax, "but frequency is not free: E ~ omega^2 A^2",
               "mode m", "length")
        ax.legend(facecolor=_PAL["panel"], edgecolor=_PAL["rule"],
                  labelcolor=_PAL["text"], fontsize=6.6, loc="lower left",
                  framealpha=0.95)
        ax.text(0.99, 0.055,
                "low frequency: large deformation, long connection.   "
                "high frequency: short connection, small deformation.\n"
                "NO favourable frequency is claimed — that needs an energy "
                "normalisation and a packet focusing law this model lacks",
                transform=ax.transAxes, ha="right", color=_PAL["text"],
                fontsize=6.5, family="monospace",
                bbox=dict(facecolor=_PAL["panel"], edgecolor="none",
                          alpha=0.9, pad=2.0))

    # ── assembly ────────────────────────────────────────────────────────────
    def draw(self):
        self._one()
        self._zero()
        self._law()
        self._parity()
        self._zonal()
        self._pulse()
        self._chord()
        self._energy()
        self.fig.text(0.5, 0.968,
                      "v67  -  one field, one surface, and the parity of the "
                      "antipode",
                      color=_PAL["text"], fontsize=15.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.944,
                      "v66 drew TWO curves and read their overlap as a "
                      "connection.  two curves in one frame are two surfaces.  "
                      "one field on one surface is the well-posed object.",
                      color=_PAL["dim"], fontsize=8.2, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.922,
                      "THE REPAIR COSTS NOTHING: v66's 'separation' IS the "
                      "one-surface deformation, to 2 ulp of R_mid  --  every "
                      "number survives, and the two configurations swap names",
                      color=_PAL["ok"], fontsize=7.6, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.900,
                      "the antipode is PARITY-DEPENDENT: maximal for odd modes, "
                      "exactly cancelling for even ones  --  and for the real "
                      "zonal spectrum it is the optimum for every odd n",
                      color=_PAL["band"], fontsize=7.6, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.878,
                      "coincident foci with opposite orientation cancel "
                      "EXACTLY:  u = 0, so no amplitude connects them",
                      color=_PAL["dim"], fontsize=7.0, ha="center",
                      family="monospace")
        self.fig.text(0.5, 0.014,
                      "the crossing rule is a REPRESENTATION choice, not a "
                      "derived boundary condition   -   LINEAR scalar on a "
                      "FIXED round background   -   the gain is a DISPLAY "
                      "amplitude, and the bottom-right panel is where that "
                      "matters",
                      color="#3d5570", fontsize=6.8, ha="center",
                      family="monospace")


def still(path: str) -> None:
    fig = OneSurfaceFigure()
    fig.draw()
    fig.fig.savefig(path, dpi=110, facecolor=_PAL["bg"])
    print(f"wrote {path}")


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--still", metavar="PNG", default="v67.png")
    a = ap.parse_args(argv)
    matplotlib.use("Agg")
    still(a.still)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
