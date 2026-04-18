"""
Four-panel geometrodynamics dashboard.

A matplotlib port of the JSX ``GeometrodynamicsViz`` prototype, wired to
the real package primitives instead of duplicating the formulas:

* **Hopf fibration** — ``hopf.connection`` (A, F, holonomy) and
  ``hopf.hopf_circle`` for the stereographically-projected fibres.
* **Wormhole throat** — ``tangherlini.V_tangherlini`` and
  ``tangherlini.solve_radial_modes`` drive the effective potential and
  the live eigenmode shape (cached per ``l``).
* **S³ Green function** — ``transaction.s3_green_potential`` and
  ``transaction.s3_green_field_kernel`` give ``G(ψ)`` and ``|dG/dψ|``.
* **Wheeler–Feynman transaction** — a narrative 4-phase illustration of
  the handshake implemented in ``transaction.handshake``.

Each panel has a pure ``draw_*`` function that renders into an existing
``Figure`` plus a ``plot_*_panel`` convenience wrapper that creates one.
``run_dashboard`` glues them together with a tab selector and a
parameter slider.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Dict, Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from geometrodynamics.constants import DELTA, R_MID, R_OUTER
from geometrodynamics.hopf import (
    hopf_circle,
    hopf_connection,
    hopf_curvature,
    hopf_holonomy,
)
from geometrodynamics.tangherlini import V_tangherlini, solve_radial_modes
from geometrodynamics.transaction import s3_green_field_kernel, s3_green_potential


# ── Palette (matches the JSX) ───────────────────────────────────────────────
_PAL: Dict[str, str] = dict(
    bg="#0a0e17",
    panel="#111827",
    border="#1e293b",
    text="#e2e8f0",
    dim="#64748b",
    accent="#60a5fa",
    accent2="#f472b6",
    accent3="#34d399",
    warn="#fbbf24",
    field="#60a5fa",
    potential="#f472b6",
    mode="#34d399",
    green="#818cf8",
    src="#60a5fa",
    dst="#f472b6",
)
_FIBER_COLORS = ("#818cf8", "#f472b6", "#34d399", "#fbbf24")


def _style_ax(ax: Axes, *, grid: bool = True) -> Axes:
    ax.set_facecolor(_PAL["bg"])
    for spine in ax.spines.values():
        spine.set_color(_PAL["border"])
    ax.tick_params(colors=_PAL["dim"], labelsize=8)
    ax.xaxis.label.set_color(_PAL["dim"])
    ax.yaxis.label.set_color(_PAL["dim"])
    ax.title.set_color(_PAL["text"])
    if grid:
        ax.grid(True, color=_PAL["border"], lw=0.4, alpha=0.5)
    return ax


def _style_fig(fig: Figure) -> Figure:
    fig.patch.set_facecolor(_PAL["bg"])
    return fig


# ═══════════════════════════════════════════════════════════════════════════
# Hopf panel — A(χ), |F(χ)|, fibres
# ═══════════════════════════════════════════════════════════════════════════
def draw_hopf(fig: Figure, chi: float = 0.6, rotate: float = 0.0) -> Dict[str, Axes]:
    """Render the Hopf fibration panel into ``fig``."""
    fig.clear()
    _style_fig(fig)
    gs = fig.add_gridspec(2, 2, width_ratios=(2.0, 1.0), wspace=0.28, hspace=0.35)
    ax3d = fig.add_subplot(gs[:, 0], projection="3d")
    ax_a = fig.add_subplot(gs[0, 1])
    ax_f = fig.add_subplot(gs[1, 1])

    ax3d.set_facecolor(_PAL["bg"])
    ax3d.set_axis_off()
    fibres = [
        (chi, 0.0, _FIBER_COLORS[0]),
        (chi, np.pi, _FIBER_COLORS[1]),
        (np.pi - chi, np.pi / 2, _FIBER_COLORS[2]),
        (np.pi - chi, -np.pi / 2, _FIBER_COLORS[3]),
    ]
    for theta, phi, color in fibres:
        x, y, z = hopf_circle(theta, phi, N=160, twist_offset=rotate, scale=1.6)
        x = np.append(x, x[0])
        y = np.append(y, y[0])
        z = np.append(z, z[0])
        ax3d.plot(x, y, z, color=color, lw=1.8, alpha=0.85)
    lim = 2.6
    ax3d.set_xlim(-lim, lim)
    ax3d.set_ylim(-lim, lim)
    ax3d.set_zlim(-lim, lim)
    ax3d.view_init(elev=22, azim=40 + np.degrees(rotate) * 0.25)
    ax3d.set_title(f"Hopf fibres through S²  (χ = {chi:.2f} rad)", color=_PAL["text"])

    chis = np.linspace(0.0, np.pi, 240)
    _style_ax(ax_a)
    ax_a.plot(chis, hopf_connection(chis), color=_PAL["field"], lw=2.0)
    ax_a.axhline(0.0, color=_PAL["border"], lw=0.7, ls=":")
    ax_a.scatter([chi], [hopf_connection(chi)], color=_PAL["warn"], s=42, zorder=3)
    ax_a.set_xlim(0, np.pi)
    ax_a.set_title("A(χ) = ½ cos(χ)", color=_PAL["field"])
    ax_a.set_xlabel("χ")

    _style_ax(ax_f)
    ax_f.plot(chis, hopf_curvature(chis), color=_PAL["accent2"], lw=2.0)
    ax_f.scatter([chi], [hopf_curvature(chi)], color=_PAL["warn"], s=42, zorder=3)
    ax_f.set_xlim(0, np.pi)
    ax_f.set_title("|F(χ)| = ½ sin(χ)", color=_PAL["accent2"])
    ax_f.set_xlabel("χ")

    hol = float(hopf_holonomy(chi))
    note = (
        f"A(χ) = {float(hopf_connection(chi)):+.3f}   "
        f"|F(χ)| = {float(hopf_curvature(chi)):.3f}   "
        f"∮A = π cos(χ) = {hol:+.3f}   c₁ = 1"
    )
    if abs(chi - np.pi / 2) < 0.15:
        note += "    ← EQUATOR: stable orbit, zero self-energy"
    elif chi < 0.2 or chi > np.pi - 0.2:
        note += "    ← POLE: holonomy = ±π  → spin-½ sign flip"
    fig.text(
        0.02, 0.02, note, color=_PAL["text"], family="monospace", fontsize=9
    )
    return dict(ax3d=ax3d, ax_a=ax_a, ax_f=ax_f)


def plot_hopf_panel(chi: float = 0.6, *, fig: Optional[Figure] = None) -> Figure:
    """Standalone figure-level wrapper for :func:`draw_hopf`."""
    fig = fig if fig is not None else plt.figure(figsize=(11, 5.2))
    draw_hopf(fig, chi=chi)
    return fig


# ═══════════════════════════════════════════════════════════════════════════
# Throat panel — V_eff(r) and eigenmode u(r)
# ═══════════════════════════════════════════════════════════════════════════
@lru_cache(maxsize=8)
def _cached_modes(l: int) -> tuple:
    oms, funcs, _ = solve_radial_modes(l)
    if not funcs:
        return (float("nan"), None, None)
    om0 = float(oms[0])
    r_full = np.asarray(funcs[0]["r_full"])
    u_full = np.asarray(funcs[0]["u_full"])
    return om0, r_full, u_full


def draw_throat(fig: Figure, l: int = 1, phase: float = 0.0) -> Dict[str, Axes]:
    """Render the Tangherlini throat panel into ``fig``."""
    fig.clear()
    _style_fig(fig)
    gs = fig.add_gridspec(2, 1, height_ratios=(1.0, 1.3), hspace=0.35)
    ax_prof = fig.add_subplot(gs[0, 0])
    ax_plot = fig.add_subplot(gs[1, 0])

    # ── Throat profile schematic ────────────────────────────────────────
    ax_prof.set_facecolor(_PAL["bg"])
    ax_prof.set_axis_off()
    s = np.linspace(-1.0, 1.0, 120)
    rn = np.abs(s) ** 2.2
    width = 0.22 + 0.55 * rn
    for side, color in ((+1, _FIBER_COLORS[0]), (-1, _FIBER_COLORS[2])):
        ax_prof.plot(s, side * width, color=color, lw=1.6, alpha=0.85)
    ax_prof.axvline(0.0, color=_PAL["warn"], lw=1.0, ls=":", alpha=0.9)
    ax_prof.text(
        0.0, -0.95, "r = R_MID  (throat)",
        ha="center", color=_PAL["warn"], fontsize=9, family="monospace",
    )
    # Oscillating eigenmode overlay
    env = np.exp(-(s * 2) ** 2)
    osc = env * np.sin(s * np.pi * (1.5 if l == 1 else 2.5 if l == 3 else 3.5))
    amp = 0.22 * np.cos(phase * (1.0 + l * 0.15))
    for side in (+1, -1):
        ax_prof.plot(s, side * (width + osc * amp), color=_PAL["mode"], lw=1.2, alpha=0.7)
    ax_prof.set_xlim(-1.15, 1.15)
    ax_prof.set_ylim(-1.05, 1.05)
    ax_prof.set_title("Wormhole throat profile (schematic)", color=_PAL["text"])

    # ── V_eff(r) and u(r) ───────────────────────────────────────────────
    _style_ax(ax_plot)
    r_grid = np.linspace(R_MID + 5e-4, R_OUTER, 400)
    V = V_tangherlini(r_grid, l)
    V_norm = V / max(np.abs(V).max(), 1e-9)
    ax_plot.plot(r_grid, V_norm, color=_PAL["potential"], lw=2.0, label="V_eff(r)  (norm.)")

    om0, r_full, u_full = _cached_modes(l)
    if r_full is not None and u_full is not None:
        scale = 0.9 * np.cos(phase * (1.0 + l * 0.15))
        ax_plot.plot(
            r_full, scale * u_full, color=_PAL["mode"], lw=2.0,
            label=f"u(r) mode 0   ω₀ = {om0:.3f}",
        )
        ax_plot.axhline(0.0, color=_PAL["border"], lw=0.6, ls=":")
    ax_plot.axvline(R_MID, color=_PAL["warn"], lw=1.0, ls="--", alpha=0.7)
    ax_plot.set_xlabel("r")
    ax_plot.set_title(
        f"V_eff(r) · [l(l+2)/r² + 3rs²/r⁴]     l = {l}     λ = l(l+2) = {l * (l + 2)}",
        color=_PAL["text"],
    )
    ax_plot.set_xlim(R_MID - DELTA * 0.2, R_OUTER + DELTA * 0.1)
    ax_plot.legend(
        loc="upper right", facecolor=_PAL["panel"], edgecolor=_PAL["border"],
        labelcolor=_PAL["text"], fontsize=8,
    )
    fig.text(
        0.02, 0.02,
        "u(R_MID) = 0  (Dirichlet at throat)   "
        "du/dr|_throat → Q → Coulomb law via Maxwell BVP",
        color=_PAL["dim"], family="monospace", fontsize=9,
    )
    return dict(ax_prof=ax_prof, ax_plot=ax_plot)


def plot_throat_panel(l: int = 1, *, fig: Optional[Figure] = None) -> Figure:
    fig = fig if fig is not None else plt.figure(figsize=(11, 5.2))
    draw_throat(fig, l=l)
    return fig


# ═══════════════════════════════════════════════════════════════════════════
# Green panel — G(ψ) and |dG/dψ| on S³
# ═══════════════════════════════════════════════════════════════════════════
def draw_green(
    fig: Figure, src_angle: float = 0.0, phase: float = 0.0
) -> Dict[str, Axes]:
    """Render the S³ Green-function panel into ``fig``."""
    fig.clear()
    _style_fig(fig)
    gs = fig.add_gridspec(2, 2, width_ratios=(1.5, 1.0), hspace=0.35, wspace=0.28)
    ax_s3 = fig.add_subplot(gs[:, 0])
    ax_g = fig.add_subplot(gs[0, 1])
    ax_dg = fig.add_subplot(gs[1, 1])

    # ── S³ refocusing schematic (2D slice) ──────────────────────────────
    ax_s3.set_facecolor(_PAL["bg"])
    ax_s3.set_aspect("equal")
    ax_s3.set_axis_off()
    theta = np.linspace(0, 2 * np.pi, 400)
    ax_s3.plot(np.cos(theta), np.sin(theta), color=_PAL["border"], lw=1.2)

    sx, sy = np.cos(src_angle), np.sin(src_angle)
    ax_s3.scatter(
        [sx], [sy], color=_PAL["src"], s=110, zorder=5,
        edgecolors=_PAL["text"], linewidths=0.7,
    )
    ax_s3.text(sx * 1.12, sy * 1.12, "source", color=_PAL["src"], fontsize=9)

    ax_angle = src_angle + np.pi
    ax, ay = np.cos(ax_angle), np.sin(ax_angle)
    ax_s3.scatter(
        [ax], [ay], color=_PAL["dst"], s=110, zorder=5,
        edgecolors=_PAL["text"], linewidths=0.7,
    )
    ax_s3.text(ax * 1.12, ay * 1.12, "antipode", color=_PAL["dst"], fontsize=9)

    for k in range(5):
        frac = ((phase * 0.5 + k / 5.0) % 1.0)
        r_out = frac * 0.98
        ax_s3.add_patch(plt.Circle(
            (sx, sy), r_out, fill=False,
            edgecolor=_PAL["src"], alpha=(1 - frac) * 0.45, lw=1.4,
        ))
        if frac > 0.6:
            r_in = (1 - (frac - 0.6) / 0.4) * 0.30
            ax_s3.add_patch(plt.Circle(
                (ax, ay), r_in, fill=False,
                edgecolor=_PAL["dst"], alpha=(frac - 0.6) / 0.4 * 0.55, lw=1.8,
            ))
    ax_s3.set_xlim(-1.35, 1.35)
    ax_s3.set_ylim(-1.35, 1.35)
    ax_s3.set_title(
        "wavefront from source refocuses at antipode  (S³ Green fn.)",
        color=_PAL["text"],
    )

    # ── G(ψ) ─────────────────────────────────────────────────────────────
    _style_ax(ax_g)
    psis = np.linspace(0.04, np.pi - 0.04, 400)
    G = np.array([s3_green_potential(p) for p in psis])
    ax_g.plot(psis, G, color=_PAL["green"], lw=2.0)
    ax_g.set_title("G(ψ) = [(π−ψ) cot ψ − ½] / (4π²R)", color=_PAL["green"])
    ax_g.set_xlabel("ψ  (geodesic distance)")
    ax_g.axvline(0.0, color=_PAL["src"], lw=0.8, ls=":")
    ax_g.axvline(np.pi, color=_PAL["dst"], lw=0.8, ls=":")
    ax_g.set_xlim(0, np.pi)

    # ── |dG/dψ| ──────────────────────────────────────────────────────────
    _style_ax(ax_dg)
    dG = np.array([s3_green_field_kernel(p) for p in psis])
    ax_dg.plot(psis, dG, color=_PAL["accent2"], lw=2.0)
    ax_dg.set_title("|dG/dψ|   (field strength, Coulomb peak at ψ → 0)", color=_PAL["accent2"])
    ax_dg.set_xlabel("ψ")
    ax_dg.set_xlim(0, np.pi)
    ax_dg.set_yscale("log")

    fig.text(
        0.02, 0.02,
        "Near source: |E| ~ Q / r²  (Coulomb).   "
        "Shell refocuses at the antipode — Wheeler–Feynman absorber channel.",
        color=_PAL["accent"], family="monospace", fontsize=9,
    )
    return dict(ax_s3=ax_s3, ax_g=ax_g, ax_dg=ax_dg)


def plot_green_panel(
    src_angle: float = 0.0, *, fig: Optional[Figure] = None
) -> Figure:
    fig = fig if fig is not None else plt.figure(figsize=(11, 5.2))
    draw_green(fig, src_angle=src_angle)
    return fig


# ═══════════════════════════════════════════════════════════════════════════
# Transaction panel — narrative 4-step handshake
# ═══════════════════════════════════════════════════════════════════════════
_PHASES = (
    ("① Retarded offer",   "GW shell expands from source",           _PAL["src"]),
    ("② GW hits candidate", "Candidate re-emits toward antipode",     _PAL["accent3"]),
    ("③ Advanced confirm", "Absorber sends time-reversed wave",       _PAL["dst"]),
    ("④ Phase closure",    "φ_offer + φ_confirm ≡ 0 or π (mod 2π)",  _PAL["warn"]),
)


def draw_handshake(fig: Figure, step: float = 0.5) -> Dict[str, Axes]:
    """Render the Wheeler–Feynman transaction panel into ``fig``.

    ``step`` is a continuous cursor in ``[0, 4)`` — integer part selects
    the phase, fractional part advances the animated overlay within it.
    """
    fig.clear()
    _style_fig(fig)
    ax = fig.add_subplot(111)
    ax.set_facecolor(_PAL["bg"])
    ax.set_aspect("equal")
    ax.set_axis_off()

    # S³ great circle
    theta = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(theta), np.sin(theta), color=_PAL["border"], lw=1.2)

    step = float(step) % 4.0
    pi_idx = int(step)
    pf = step - pi_idx

    src_ang = np.pi * 0.85
    sx, sy = np.cos(src_ang), np.sin(src_ang)
    dst_ang = src_ang + np.pi
    dx, dy = np.cos(dst_ang), np.sin(dst_ang)
    cand_ang = src_ang - np.pi * 0.4
    cx, cy = np.cos(cand_ang), np.sin(cand_ang)

    # Phase ①: expanding offer shell
    if pi_idx == 0:
        for k in range(3):
            r = pf * 1.5 - k * 0.09
            if r > 0:
                ax.add_patch(plt.Circle(
                    (sx, sy), r, fill=False,
                    edgecolor=_PAL["src"], alpha=0.30 - k * 0.08, lw=2.0,
                ))
    # Phase ②: src→candidate propagation
    if pi_idx >= 1:
        p = pf if pi_idx == 1 else 1.0
        ex = sx + (cx - sx) * p
        ey = sy + (cy - sy) * p
        ax.plot([sx, ex], [sy, ey], color=_PAL["accent3"], lw=2.0, ls="--", alpha=0.9)
        if p > 0.8:
            ax.add_patch(plt.Circle(
                (cx, cy), (p - 0.8) * 1.5 * 0.30, fill=False,
                edgecolor=_PAL["accent3"], alpha=0.35, lw=1.6,
            ))
    # Phase ③: dst→candidate confirmation (time-reversed)
    if pi_idx >= 2:
        p = pf if pi_idx == 2 else 1.0
        ex = dx + (cx - dx) * p
        ey = dy + (cy - dy) * p
        ax.plot([dx, ex], [dy, ey], color=_PAL["dst"], lw=2.0, ls=":", alpha=0.95)
        if p < 0.6:
            ax.add_patch(plt.Circle(
                (cx, cy), (0.6 - p) * 0.5, fill=False,
                edgecolor=_PAL["dst"], alpha=p * 0.5, lw=1.6,
            ))
    # Phase ④: phase closure flash
    if pi_idx == 3:
        fl = 0.5 + 0.5 * np.sin(pf * np.pi * 4)
        ax.plot(
            [sx, cx, dx], [sy, cy, dy],
            color=_PAL["warn"], lw=3.0, alpha=fl, solid_capstyle="round",
        )
        if pf > 0.5:
            ax.text(0.0, -1.22, "TRANSACTION CONFIRMED",
                    ha="center", color=_PAL["warn"], fontsize=13,
                    family="monospace", fontweight="bold")
            ax.text(0.0, -1.32, "φ_total ≡ 0  (mod 2π)",
                    ha="center", color=_PAL["dim"], fontsize=10,
                    family="monospace")

    for (px, py, color, label, sub) in (
        (sx, sy, _PAL["src"], "source", "q = +Q"),
        (cx, cy, _PAL["accent3"], "candidate", "relay"),
        (dx, dy, _PAL["dst"], "absorber", "q = −Q  (antipodal)"),
    ):
        ax.scatter([px], [py], color=color, s=120, zorder=5,
                   edgecolors=_PAL["text"], linewidths=0.8)
        ax.text(px * 1.12, py * 1.12, label, color=color, fontsize=10,
                family="monospace")
        ax.text(px * 1.12, py * 1.12 - 0.08, sub, color=_PAL["dim"],
                fontsize=8, family="monospace")

    label, sub, color = _PHASES[pi_idx]
    ax.text(-1.45, 1.32, label, color=color, fontsize=14,
            family="monospace", fontweight="bold")
    ax.text(-1.45, 1.22, sub, color=_PAL["dim"], fontsize=10,
            family="monospace")

    for i, (p_label, _, p_color) in enumerate(_PHASES):
        active = (i == pi_idx)
        ax.text(
            -1.35 + i * 0.9, -1.40, p_label,
            color=p_color if active else _PAL["dim"],
            fontsize=10 if active else 9,
            family="monospace",
            fontweight="bold" if active else "normal",
        )

    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.45)
    return dict(ax=ax)


def plot_handshake_panel(
    step: float = 0.5, *, fig: Optional[Figure] = None
) -> Figure:
    fig = fig if fig is not None else plt.figure(figsize=(9, 6.0))
    draw_handshake(fig, step=step)
    return fig


# ═══════════════════════════════════════════════════════════════════════════
# Interactive dashboard
# ═══════════════════════════════════════════════════════════════════════════
_TABS = (
    ("hopf",     "Hopf fibration"),
    ("throat",   "Wormhole throat"),
    ("green",    "S³ Green fn."),
    ("handshake","Transaction"),
)


def run_dashboard(panel: str = "hopf", *, show: bool = True) -> Figure:
    """Launch the interactive 4-panel dashboard with tab + parameter widgets.

    Radio buttons switch panels; a slider controls the active panel's
    parameter (χ for Hopf, l for throat, source angle for Green, step for
    handshake).
    """
    from matplotlib.widgets import RadioButtons, Slider  # local import

    fig = plt.figure(figsize=(12, 6.5))
    _style_fig(fig)
    panel_fig = fig

    state = dict(tab=panel, chi=0.6, l_idx=0, src=0.0, step=0.5)
    l_choices = (1, 3, 5)

    radio_ax = fig.add_axes([0.012, 0.78, 0.11, 0.18], facecolor=_PAL["panel"])
    radio_ax.set_facecolor(_PAL["panel"])
    radio = RadioButtons(
        radio_ax, [lab for _, lab in _TABS],
        active=[key for key, _ in _TABS].index(panel),
        activecolor=_PAL["accent"],
    )
    for label in radio.labels:
        label.set_color(_PAL["text"])
        label.set_fontsize(9)

    slider_ax = fig.add_axes([0.16, 0.03, 0.70, 0.035], facecolor=_PAL["panel"])
    slider = Slider(
        slider_ax, "χ", 0.02, np.pi - 0.02, valinit=state["chi"],
        color=_PAL["accent"],
    )
    slider.label.set_color(_PAL["text"])
    slider.valtext.set_color(_PAL["text"])

    plot_ax_rect = [0.15, 0.10, 0.82, 0.82]

    def _redraw() -> None:
        # Rebuild the panel content using the drawing functions.  We wipe
        # everything that isn't our widget axes each time.
        keep = {radio_ax, slider_ax}
        for a in list(panel_fig.axes):
            if a not in keep:
                panel_fig.delaxes(a)
        if state["tab"] == "hopf":
            draw_hopf(panel_fig, chi=state["chi"])
        elif state["tab"] == "throat":
            draw_throat(panel_fig, l=l_choices[state["l_idx"]])
        elif state["tab"] == "green":
            draw_green(panel_fig, src_angle=state["src"])
        else:
            draw_handshake(panel_fig, step=state["step"])
        panel_fig.canvas.draw_idle()

    def _on_tab(label: str) -> None:
        lookup = {lab: key for key, lab in _TABS}
        state["tab"] = lookup[label]
        if state["tab"] == "hopf":
            slider.label.set_text("χ")
            slider.valmin, slider.valmax = 0.02, np.pi - 0.02
            slider.set_val(state["chi"])
        elif state["tab"] == "throat":
            slider.label.set_text("l idx")
            slider.valmin, slider.valmax = 0, len(l_choices) - 1
            slider.set_val(state["l_idx"])
        elif state["tab"] == "green":
            slider.label.set_text("source angle")
            slider.valmin, slider.valmax = 0.0, 2 * np.pi
            slider.set_val(state["src"])
        else:
            slider.label.set_text("step")
            slider.valmin, slider.valmax = 0.0, 3.999
            slider.set_val(state["step"])
        _redraw()

    def _on_slider(val: float) -> None:
        if state["tab"] == "hopf":
            state["chi"] = float(val)
        elif state["tab"] == "throat":
            state["l_idx"] = int(round(val))
        elif state["tab"] == "green":
            state["src"] = float(val)
        else:
            state["step"] = float(val)
        _redraw()

    radio.on_clicked(_on_tab)
    slider.on_changed(_on_slider)
    _redraw()

    if show:
        plt.show()
    return fig


__all__ = [
    "draw_hopf",
    "draw_throat",
    "draw_green",
    "draw_handshake",
    "plot_hopf_panel",
    "plot_throat_panel",
    "plot_green_panel",
    "plot_handshake_panel",
    "run_dashboard",
]
