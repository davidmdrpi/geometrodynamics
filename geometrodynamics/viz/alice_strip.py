"""Charge without charge built up in four panels, on six Alice handles.

A schematic, not a solved field. One camera and one geometry, with one layer
added per panel:

1. **Wheeler.** Two sheets and six throats in three pairs, as in Wheeler's
   multi-throat picture, in one row facing the camera. Sourceless electric flux runs down one throat of each
   pair and up the other. Gauss's law, ∮E·dA = ±4πq, gives each mouth an
   apparent charge: red (+) where flux leaves onto a sheet, blue (−) where it
   enters.
2. **+ Kaluza–Klein.** A zoom on one waist shows the Hopf fibre χ over the
   throat. The matter wave is the winding mode Ψ = ψ(x)e^{ikχ}, and the charge
   is the winding k (``docs/guidance_law_from_5d.md``). Flux through a throat is
   fibre momentum, not a pasted-on source.
3. **+ Compact antipodal closure.** The two sheets are hemispheres of one closed
   surface, joined at one equator, and the bottom point (x, y) is the top point
   (−x, −y). Dashed rings show the six bottom mouths where they sit on that one
   surface. A pulse from a mouth p crosses the equator and refocuses at −p on
   the cover (``antipodal_focusing``), and −p is p itself. Focusing alone
   creates nothing (#329).
4. **+ Alice seams.** Every handle is the seam type s_t = s_f = −1. It reverses
   time orientation and the fibre, χ → −χ, so q → −q. That is the only
   time-reversing type that lets trapped flux through its tube
   (``docs/klein_handle_electric.md``, E2). Each flux loop is one closed
   worldline, gold where proper time runs with the circulation and cyan where
   it runs against it.

Where the seam sits along a throat is a choice of chart. That each loop
crosses two seams is not.
"""
from __future__ import annotations

import numpy as np

from geometrodynamics.viz.alice_circulation import (
    BLUE, CYAN, DEPTH, FRONT, GOLD, GREEN, HEIGHT, PAIR, RED, RIM, SEAM, VIEW, Point, _arrow, column_radius,
    identify,
)

SHEET_RADIUS = 2.9
# Six throats in one row facing the camera, three pairs as in Wheeler's figure. The row
# sits in front of the sheet centre, so the antipodal images form a separate back row.
_ACROSS = np.array([-VIEW[1], VIEW[0]]) / np.linalg.norm(VIEW[:2])
_TOWARD = VIEW[:2] / np.linalg.norm(VIEW[:2])
ROW = (-2.05, -1.3, -0.375, 0.375, 1.3, 2.05)
FORWARD = 1.0
COLUMNS = tuple(t * _ACROSS + FORWARD * _TOWARD for t in ROW)
CIRCUITS = ((0, 1), (3, 2), (4, 5))          # (down handle, up handle) of each flux loop
WAVE = "#d36be8"
WAVE_SOURCE = 2                              # the top mouth the panel-3 pulse starts from
CHORD_HANDLE = 0                             # the bottom mouth paired with its image in panel 3
ZOOM_HANDLE = 1                              # the waist shown in the panel-2 inset


# ---- geometry ---------------------------------------------------------------------------------

def mouths_on_single_surface():
    """Each handle's two mouths as points of the one surface (top-sheet coordinates)."""
    return [(c.copy(), identify(c)) for c in COLUMNS]


def sheet_height(xy, top: bool) -> np.ndarray:
    """Embedding height of a sheet, with a funnel into each of the six mouths."""
    xy = np.atleast_2d(xy)
    dip = np.zeros(len(xy))
    for c in COLUMNS:
        rho = np.maximum(np.linalg.norm(xy - c, axis=1), RIM)
        dip += DEPTH * (RIM / rho) ** 3.5
    z = HEIGHT / 2 - dip
    return z if top else -z


def rim_point(c):
    return c + RIM * FRONT


def mouth_height(i: int) -> float:
    """Height of handle i's top mouth: where the funnel meets the column."""
    return float(sheet_height(rim_point(COLUMNS[i]), True)[0])


def throat_radius(i: int, z):
    """Column radius of handle i: the finite-mouth neck, scaled to this mouth's height."""
    return column_radius(np.asarray(z) * (HEIGHT / 2 - DEPTH) / mouth_height(i))


def mouth_charges():
    """Apparent Gauss charge of every drawn mouth, keyed by (sheet, handle index)."""
    q = {}
    for down, up in CIRCUITS:
        q[("top", down)], q[("bottom", down)] = -1, +1     # flux enters at the top, leaves at the bottom
        q[("top", up)], q[("bottom", up)] = +1, -1
    return q


def seam_signs():
    """(s_t, s_f) of every handle: each reverses time orientation and the fibre."""
    return {i: (-1, -1) for i in range(len(COLUMNS))}


def winding_after_seam(k: int, s_f: int) -> int:
    """The fibre winding (the charge) after crossing a seam that maps χ → s_f χ."""
    return s_f * k


# ---- the flux loops as closed worldlines -------------------------------------------------------

def _sheet_path(p, q, t, top, bulge):
    a, b = rim_point(p), rim_point(q)
    mid = 0.5 * (p + q)
    out = mid / (np.linalg.norm(mid) + 1e-12)
    xy = (1 - t) * a + t * b + bulge * np.sin(np.pi * t) * out
    for c in (p, q):
        d = xy - c
        rho = np.linalg.norm(d)
        if rho < RIM:
            xy = c + d / (rho + 1e-12) * RIM
    return np.array([xy[0], xy[1], sheet_height(xy, top)[0]])


def _throat_path(i, t, down):
    zr = mouth_height(i)
    z = zr - 2 * zr * t if down else -zr + 2 * zr * t
    xy = COLUMNS[i] + throat_radius(i, z) * FRONT
    return np.array([xy[0], xy[1], z])


SEGMENTS = (("top", 0.30), ("down", 0.20), ("bottom", 0.30), ("up", 0.20))
BULGE = 0.35


def worldline(k: int, s: float) -> Point:
    """Flux loop k at parameter s ∈ [0, 1): top sheet, down throat, bottom sheet, up throat."""
    down, up = CIRCUITS[k]
    s = s % 1.0
    edges = np.cumsum([0.0] + [w for _, w in SEGMENTS])
    j = int(np.searchsorted(edges, s, side="right") - 1)
    t = (s - edges[j]) / SEGMENTS[j][1]
    name = SEGMENTS[j][0]

    def at(tt):
        if name == "top":
            return _sheet_path(COLUMNS[up], COLUMNS[down], tt, True, BULGE)
        if name == "down":
            return _throat_path(down, tt, True)
        if name == "bottom":
            return _sheet_path(COLUMNS[down], COLUMNS[up], tt, False, BULGE)
        return _throat_path(up, tt, False)

    h = 1e-4
    xyz = at(t)
    tangent = at(min(t + h, 1.0)) - at(max(t - h, 0.0))
    tangent /= np.linalg.norm(tangent) + 1e-15
    if name == "top":
        orient = +1
    elif name == "bottom":
        orient = -1
    elif name == "down":
        orient = +1 if t < 0.5 else -1
    else:
        orient = -1 if t < 0.5 else +1
    return Point(xyz, tangent, orient, name if name in ("top", "bottom") else "throat")


# ---- the exterior pulse on the closed surface ---------------------------------------------------

def _on_sphere(xy, top):
    """Read a sheet point as a point of the sphere of radius SHEET_RADIUS (top: upper hemisphere)."""
    zs = np.sqrt(max(SHEET_RADIUS ** 2 - float(np.dot(xy, xy)), 0.0))
    return np.array([xy[0], xy[1], zs if top else -zs])


def wavefront(theta: float, n: int = 361):
    """The geodesic circle at angle θ about the source mouth, as sheet points and a top/bottom mask.

    On the sphere the pulse reaches the equator and refocuses at the antipode −p
    at θ = π. That antipode lies on the bottom sheet at −c, which the
    identification makes the source point itself.
    """
    p = _on_sphere(COLUMNS[WAVE_SOURCE], True) / SHEET_RADIUS
    e1 = np.cross(p, [0.0, 0.0, 1.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(p, e1)
    phi = np.linspace(0, 2 * np.pi, n)
    X = SHEET_RADIUS * (np.cos(theta) * p[None, :]
                        + np.sin(theta) * (np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2))
    return X[:, :2], X[:, 2] >= 0


# ---- drawing -----------------------------------------------------------------------------------

TH = np.linspace(0, 2 * np.pi, 121)


def _draw_base(ax, labels):
    for top in (True, False):
        for R in np.linspace(0.3, SHEET_RADIUS, 16):
            xy = np.column_stack([R * np.cos(TH), R * np.sin(TH)])
            ax.plot(xy[:, 0], xy[:, 1], sheet_height(xy, top), color=GREEN, lw=0.3, alpha=0.5)
        for phi in np.linspace(0, 2 * np.pi, 36, endpoint=False):
            R = np.linspace(0.05, SHEET_RADIUS, 90)
            xy = np.column_stack([R * np.cos(phi), R * np.sin(phi)])
            ax.plot(xy[:, 0], xy[:, 1], sheet_height(xy, top), color=GREEN, lw=0.22, alpha=0.32)
        for c in COLUMNS:
            for rho in np.geomspace(RIM * 1.02, 0.62, 6):
                xy = c + rho * np.column_stack([np.cos(TH), np.sin(TH)])
                ax.plot(xy[:, 0], xy[:, 1], sheet_height(xy, top), color=GREEN, lw=0.35, alpha=0.55)
    q = mouth_charges()
    for i, c in enumerate(COLUMNS):
        zr = mouth_height(i)
        z = np.linspace(-zr, zr, 90)
        r = throat_radius(i, z)
        for phi in np.linspace(0, 2 * np.pi, 16, endpoint=False):
            ax.plot(c[0] + r * np.cos(phi), c[1] + r * np.sin(phi), z, color=GREEN, lw=0.3, alpha=0.45)
        for zz in np.linspace(-zr, zr, 7):
            rr = throat_radius(i, zz)
            ax.plot(c[0] + rr * np.cos(TH), c[1] + rr * np.sin(TH), np.full_like(TH, zz), color=GREEN, lw=0.28, alpha=0.38)
        for sheet, zz in (("top", zr), ("bottom", -zr)):
            col = RED if q[(sheet, i)] > 0 else BLUE
            ax.plot(c[0] + RIM * np.cos(TH), c[1] + RIM * np.sin(TH), np.full_like(TH, zz), color=col, lw=2.0)
            if labels:
                ax.text(c[0] + 1.3 * RIM * FRONT[0], c[1] - 1.3 * RIM * FRONT[1] - 0.1,
                        zz + (0.2 if sheet == "top" else -0.3), "+" if q[(sheet, i)] > 0 else "−",
                        color=col, fontsize=10, weight="bold")
    # sourceless flux of each loop: a bundle on each sheet and through each throat
    for down, up in CIRCUITS:
        for top in (True, False):
            a, b = (up, down) if top else (down, up)
            for bulge in np.linspace(-0.3, 0.6, 7):
                t = np.linspace(0.02, 0.98, 80)
                P = np.array([_sheet_path(COLUMNS[a], COLUMNS[b], tt, top, bulge) for tt in t])
                ax.plot(P[:, 0], P[:, 1], P[:, 2] + 0.01, color=GOLD, lw=0.6, ls=(0, (4, 3)), alpha=0.75)
        for i in (down, up):
            c, zr = COLUMNS[i], mouth_height(i)
            z = np.linspace(-zr, zr, 60)
            r = 0.55 * throat_radius(i, z)
            for phi in np.linspace(0, 2 * np.pi, 4, endpoint=False):
                ax.plot(c[0] + r * np.cos(phi), c[1] + r * np.sin(phi), z, color=GOLD, lw=0.6, ls=(0, (4, 3)), alpha=0.7)


def _draw_identification(ax, alpha, chord):
    q = mouth_charges()
    for i, c in enumerate(COLUMNS):
        xy = identify(c) + RIM * np.column_stack([np.cos(TH), np.sin(TH)])
        col = RED if q[("bottom", i)] > 0 else BLUE
        ax.plot(xy[:, 0], xy[:, 1], sheet_height(xy, True) + 0.02, color=col, lw=1.5, ls=(0, (3, 2)), alpha=alpha)
    if chord:
        c = COLUMNS[CHORD_HANDLE]
        zb = -mouth_height(CHORD_HANDLE)
        img = identify(c)
        zi = float(sheet_height(img, True)[0])
        ax.plot([c[0], img[0]], [c[1], img[1]], [zb, zi], color=PAIR, lw=0.9, ls=(0, (1, 3)), alpha=0.75)
        _arrow(ax, np.array([c[0], c[1], zb]) + 0.12 * VIEW, np.array([0, 0, -0.4]), PAIR)
        _arrow(ax, np.array([img[0], img[1], zi]), np.array([0, 0, 0.4]), PAIR)
        ax.text(img[0], img[1], zi + 0.55, "same point", color=PAIR, fontsize=7)


def _draw_closure(ax):
    # one equator: the two rims are the same circle
    perp = np.array([-VIEW[1], VIEW[0]]) / np.linalg.norm(VIEW[:2])
    t = np.linspace(0, 1, 60)
    for sgn in (+1, -1):
        xy = sgn * SHEET_RADIUS * perp[None, :] * (1 + 0.04 * np.sin(np.pi * t))[:, None]
        z = HEIGHT / 2 * np.cos(np.pi * t)
        ax.plot(xy[:, 0], xy[:, 1], z, color=PAIR, lw=0.9, ls=(0, (1, 3)), alpha=0.7)
    ax.text(*(SHEET_RADIUS * 0.86 * perp), 0.0, "one\nequator", color=PAIR, fontsize=7, ha="center")
    # the pulse from p: out across the top sheet, over the equator, in to −p on the bottom
    for theta in (0.45, 1.05, 2.1, np.pi - 0.45):
        xy, top = wavefront(theta)
        for sheet in (True, False):
            m = top == sheet
            if not m.any():
                continue
            P = np.column_stack([xy, sheet_height(xy, sheet) + (0.03 if sheet else 0.03)])
            P[~m] = np.nan
            for c in COLUMNS:
                P[np.linalg.norm(xy - c, axis=1) < RIM * 1.05] = np.nan
            ax.plot(P[:, 0], P[:, 1], P[:, 2], color=WAVE, lw=0.9, alpha=0.6)
    src = COLUMNS[WAVE_SOURCE]
    foc = identify(src)
    zf = float(sheet_height(foc, False)[0])
    ax.scatter([foc[0]], [foc[1]], [zf + 0.03], s=60, marker="*", color=WAVE, depthshade=False)
    ax.text(foc[0] + 0.5 * perp[0], foc[1] + 0.5 * perp[1], zf + 0.12, "−p ≡ p", color=WAVE, fontsize=8)
    ax.text(src[0] + 0.2, src[1] + 0.2, mouth_height(WAVE_SOURCE) + 0.35, "p", color=WAVE, fontsize=9)


def _draw_alice(ax):
    for k in range(len(CIRCUITS)):
        s = np.linspace(0, 1, 1201)
        W = [worldline(k, x) for x in s]
        P = np.array([w.xyz for w in W]) + 0.02 * np.r_[0, 0, 1]
        O = np.array([w.orientation for w in W])
        for sign, col in ((+1, GOLD), (-1, CYAN)):
            Q = P.copy()
            Q[O != sign] = np.nan
            ax.plot(Q[:, 0], Q[:, 1], Q[:, 2], color=col, lw=1.5, alpha=0.95)
    for i, c in enumerate(COLUMNS):
        r = 1.2 * throat_radius(i, 0.0)
        ax.plot(c[0] + r * np.cos(TH), c[1] + r * np.sin(TH), np.zeros_like(TH), color=SEAM, lw=1.2, alpha=0.95)
        u = np.linspace(-1, 1, 30)
        ang = np.arctan2(FRONT[1], FRONT[0]) + 0.9 * u
        ax.plot(c[0] + 1.5 * r * np.cos(ang), c[1] + 1.5 * r * np.sin(ang), 0.12 * np.tanh(2.5 * u),
                color=SEAM, lw=1.0, alpha=0.95)
    c = COLUMNS[5]
    ax.text(c[0] - 0.2, c[1] + 0.55, 0.0, "Alice seam", color=SEAM, fontsize=7)
    balls, ghosts = [], []
    for k in range(len(CIRCUITS)):
        for s0 in (0.12 + 0.11 * k, 0.47 + 0.11 * k):
            p = worldline(k, s0)
            col = GOLD if p.orientation > 0 else CYAN
            base = p.xyz + np.array([0, 0, 0.03])
            balls.append((base, col))
            _arrow(ax, base + 0.09 * VIEW, 0.36 * p.orientation * p.tangent, col)
            if p.sheet == "bottom":
                g = identify(p.xyz[:2])
                ghosts.append([g[0], g[1], float(sheet_height(g, True)[0]) + 0.03])
    B = np.array([b for b, _ in balls])
    ax.scatter(B[:, 0], B[:, 1], B[:, 2], s=40, c=[c for _, c in balls], depthshade=False)
    if ghosts:
        G = np.array(ghosts)
        ax.scatter(G[:, 0], G[:, 1], G[:, 2], s=40, facecolors="none", edgecolors=CYAN, linewidths=1.2, depthshade=False)


def _kk_inset(fig, rect):
    ax = fig.add_axes(rect, facecolor="black")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.62, 0.62)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#666666")
    for y in (0.3, -0.3):
        ax.plot([0.05, 0.95], [y, y], color=GREEN, lw=0.9, alpha=0.8)
    th = np.linspace(0, 2 * np.pi, 80)
    for x in np.linspace(0.1, 0.9, 5):
        front = np.cos(th) >= 0
        ex, ey = x + 0.035 * np.cos(th), 0.3 * np.sin(th)
        ax.plot(np.where(front, ex, np.nan), ey, color=GREEN, lw=0.8, alpha=0.8)
        ax.plot(np.where(~front, ex, np.nan), ey, color=GREEN, lw=0.6, ls=":", alpha=0.6)
    t = np.linspace(0, 1, 400)
    x = 0.07 + 0.84 * t
    ang = 2 * np.pi * 3.5 * t
    y, depth = 0.3 * np.cos(ang), np.sin(ang)
    ax.plot(np.where(depth >= 0, x, np.nan), y, color=GOLD, lw=1.6)
    ax.plot(np.where(depth < 0, x, np.nan), y, color=GOLD, lw=0.9, ls=":")
    ax.annotate("", xy=(x[-1], y[-1]), xytext=(x[-8], y[-8]),
                arrowprops=dict(arrowstyle="-|>", color=GOLD, lw=1.4))
    ax.text(0.5, 0.45, "fibre S¹ (χ) over each point of the waist", color="#cccccc", fontsize=6.5, ha="center")
    ax.text(0.5, -0.55, "winding k  =  charge q", color=GOLD, fontsize=7, ha="center")
    return ax


TITLES = (
    "1  Wheeler: sourceless flux through six throats",
    "2  + Kaluza–Klein: charge is fibre winding",
    "3  + closure: top ≡ bottom, antipodal refocus",
    "4  + Alice seams: both time orientations",
)


def render(path_png, dpi: int = 110):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from mpl_toolkits.mplot3d import proj3d

    fig = plt.figure(figsize=(20, 7.4), facecolor="black")
    axes = []
    for i in range(4):
        x0 = 0.0025 + 0.25 * i
        ax = fig.add_axes([x0, 0.15, 0.245, 0.8], projection="3d", facecolor="none")
        ax.set_axis_off()
        ax.view_init(elev=13, azim=45)
        lim = 0.71 * SHEET_RADIUS
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-HEIGHT / 2 - 0.15, HEIGHT / 2 + 0.15)
        ax.set_box_aspect((1, 1, 1.05))
        _draw_base(ax, labels=False)
        fig.text(x0 + 0.122, 0.945, TITLES[i], color="white", ha="center", fontsize=11.5)
        axes.append(ax)
    # 1: Wheeler
    fig.text(0.125, 0.115, r"$\oint \mathbf{E}\cdot d\mathbf{A} = \pm 4\pi q$", color="white", ha="center", fontsize=17)
    fig.text(0.125, 0.045, "three flux loops, no source anywhere\nred +: flux leaves a mouth   ·   blue −: flux enters",
             color="#bbbbbb", ha="center", fontsize=8.5, linespacing=1.5)
    # 2: Kaluza–Klein
    c, zw = COLUMNS[ZOOM_HANDLE], 0.0
    r = 1.6 * throat_radius(ZOOM_HANDLE, zw)
    axes[1].plot(c[0] + r * np.cos(TH), c[1] + r * np.sin(TH), np.full_like(TH, zw), color="white", lw=1.0)
    inset = _kk_inset(fig, [0.285, 0.035, 0.18, 0.12])
    fig.text(0.375, 0.165, r"$\Psi = \psi(x)\,e^{ik\chi}$:  flux through a throat is fibre momentum",
             color="#bbbbbb", ha="center", fontsize=8.5)
    # 3: closure
    _draw_identification(axes[2], alpha=1.0, chord=True)
    _draw_closure(axes[2])
    fig.text(0.625, 0.045,
             "(x, y)↓ ~ (−x, −y)↑: dashed rings are the bottom mouths on the one surface\n"
             "magenta: a pulse from p crosses the equator and refocuses at −p, which is p\n"
             "focusing alone creates nothing (#329)",
             color="#bbbbbb", ha="center", fontsize=8.5, linespacing=1.5)
    # 4: Alice
    _draw_identification(axes[3], alpha=0.45, chord=False)
    _draw_alice(axes[3])
    fig.text(0.875, 0.045,
             "gold: proper time with the circulation   ·   cyan: reversed\n"
             "white band: Alice seam, s_t = s_f = −1:  t → −t,  χ → −χ,  q → −q\n"
             "hollow: a bottom-sheet ball at its image on the one surface",
             color="#bbbbbb", ha="center", fontsize=8.5, linespacing=1.5)
    fig.text(0.998, 0.005, "schematic, no field equation solved   ·   geometrodynamics/viz/alice_strip.py",
             color="#777777", ha="right", fontsize=7)
    # connect the zoom ring to the inset
    fig.canvas.draw()
    x2, y2, _ = proj3d.proj_transform(c[0] + r * FRONT[0], c[1] + r * FRONT[1], zw, axes[1].get_proj())
    fx, fy = fig.transFigure.inverted().transform(axes[1].transData.transform((x2, y2)))
    fig.add_artist(Line2D([fx, 0.375], [fy, 0.155], transform=fig.transFigure, color="white", lw=0.7, ls=":"))
    fig.savefig(path_png, dpi=dpi, facecolor="black")
    plt.close(fig)
    return inset


if __name__ == "__main__":
    from pathlib import Path
    out = Path(__file__).resolve().parents[2] / "docs" / "images" / "alice_circulation"
    out.mkdir(parents=True, exist_ok=True)
    render(out / "strip.png", dpi=130)
    print("wrote", out / "strip.png")
