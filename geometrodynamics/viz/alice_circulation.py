"""Charge without charge on one antipodally closed sheet, with two Alice handles.

A schematic, not a solved field. It draws Wheeler's sourceless electric flux in
the BAM setting of #338/#340: the two sheets shown are **one** surface, with the
bottom point (x, y) identified with the top point (−x, −y). Two vertical
columns are the two Alice handles. Column 1 stands at (d, 0) and column 2 at
(0, d), so after the identification each joins a pair of antipodal points of
the single surface: four mouths in two antipodal pairs.

Flux circulates top sheet → throat 1 → bottom sheet → throat 2 → top sheet with
no source. Gauss's law, ∮E·dA = ±4πq, then gives each mouth an apparent charge:
positive (red) where flux leaves onto a sheet, negative (blue) where it enters.

Balls circulate on the same closed worldline. Both handles reverse time
orientation, so a ball's proper-time arrow runs with the circulation on the top
sheet and against it on the bottom sheet. One circular worldline carries both
timelike directions. A ball on the bottom sheet is also shown as a hollow ghost
at its antipode on the top sheet, because they are the same surface.

Each handle is the time-reversing, charge-conjugating seam type, s_t = s_f = −1.
That is the only time-reversing type that admits flux through the tube
(``docs/klein_handle_electric.md``, E2: trapped flux needs s_t·s_f = +1). The
flip is drawn as a seam at each waist. On a non-orientable space where to cut
is a choice of chart; that the circuit crosses two such seams is not.

The throat profile is the finite-mouth neck r(s) = sqrt(s² + b²) of
``bulk.finite_mouth`` with a = 0.3, stretched vertically for display.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

A_MOUTH = 0.3
RIM = np.sin(A_MOUTH)                 # mouth radius, sin a
WAIST = np.sin(A_MOUTH) ** 2          # neck waist, b = sin^2 a
HEIGHT = 3.0                          # sheet separation (display only)
DEPTH = 0.45                          # funnel depth below each sheet
OFFSET = 1.45                         # column distance from the sheet centre
SHEET_RADIUS = 2.7
COLUMNS = (np.array([OFFSET, 0.0]), np.array([0.0, OFFSET]))   # handle 1, handle 2


def identify(xy_bottom: np.ndarray) -> np.ndarray:
    """The top-sheet point equal to a bottom-sheet point: (x, y)↓ ~ (−x, −y)↑."""
    return -np.asarray(xy_bottom, dtype=float)


def mouths_on_single_surface():
    """The four mouths as points of the one surface (top-sheet coordinates), by handle."""
    return [(c.copy(), identify(c)) for c in COLUMNS]


def sheet_height(xy: np.ndarray, top: bool) -> np.ndarray:
    """Embedding height of a sheet, with a funnel dipping into every mouth on it."""
    xy = np.atleast_2d(xy)
    dip = np.zeros(len(xy))
    for c in COLUMNS:
        rho = np.maximum(np.linalg.norm(xy - c, axis=1), RIM)
        dip += DEPTH * (RIM / rho) ** 2.2
    z = HEIGHT / 2 - dip
    return z if top else -z


def column_radius(z: np.ndarray) -> np.ndarray:
    """Throat radius: the hyperbolic neck profile, rim radius at the funnel bottom."""
    zr = HEIGHT / 2 - DEPTH
    k = np.sqrt(RIM ** 2 - WAIST ** 2) / zr
    return np.sqrt(WAIST ** 2 + (k * np.asarray(z)) ** 2)


# ---- flux and charges -----------------------------------------------------------------------

def mouth_charges():
    """Apparent charge sign of each drawn mouth from the circulating flux (Gauss's law).

    The flux runs down throat 1 and up throat 2. On the top sheet it leaves
    mouth 2 (+) and enters mouth 1 (−); on the bottom sheet it leaves mouth 1 (+)
    and enters mouth 2 (−).
    """
    return {("top", 1): -1, ("top", 2): +1, ("bottom", 1): +1, ("bottom", 2): -1}


def seam_signs():
    """(s_t, s_f) of each Alice handle: time orientation and fibre (charge) reversal."""
    return {1: (-1, -1), 2: (-1, -1)}


def flux_through_tube_allowed(s_t: int, s_f: int) -> bool:
    """E2 of the Klein-handle electric test: trapped flux threads the tube iff s_t·s_f = +1."""
    return s_t * s_f == +1


# ---- the closed worldline -------------------------------------------------------------------

@dataclass
class Point:
    xyz: np.ndarray
    tangent: np.ndarray
    orientation: int          # +1: proper time along the circulation, −1: reversed
    sheet: str                # "top", "bottom" or "throat"


FRONT = np.array([0.7071, 0.7071])          # unit vector of the generator facing the camera


def rim_point(c: np.ndarray) -> np.ndarray:
    """Where the circuit crosses a mouth's rim: on the camera-facing generator."""
    return c + RIM * FRONT


def _sheet_path(p: np.ndarray, q: np.ndarray, t: float, top: bool, bulge: float) -> np.ndarray:
    """A curve on a sheet from mouth p to mouth q, rim to rim, bowed away from the sheet centre."""
    a, b = rim_point(p), rim_point(q)
    mid = 0.5 * (p + q)
    out = mid / (np.linalg.norm(mid) + 1e-12)
    xy = (1 - t) * a + t * b + bulge * np.sin(np.pi * t) * out
    for c in (p, q):                                   # stay on the sheet, outside both mouths
        d = xy - c
        rho = np.linalg.norm(d)
        if rho < RIM:
            xy = c + d / (rho + 1e-12) * RIM
    z = sheet_height(xy, top)[0]
    return np.array([xy[0], xy[1], z])


def _throat_path(c: np.ndarray, t: float, down: bool) -> np.ndarray:
    """Down or up the camera-facing generator of a column, between the two sheets' rim points."""
    zr = sheet_height(rim_point(c), True)[0]
    z = zr - 2 * zr * t if down else -zr + 2 * zr * t
    xy = c + column_radius(z) * FRONT
    return np.array([xy[0], xy[1], z])


SEGMENTS = (("top", 0.30), ("down1", 0.20), ("bottom", 0.30), ("up2", 0.20))


def worldline(s: float) -> Point:
    """The closed circuit at parameter s ∈ [0, 1): top sheet, throat 1, bottom sheet, throat 2."""
    s = s % 1.0
    edges = np.cumsum([0.0] + [w for _, w in SEGMENTS])
    k = int(np.searchsorted(edges, s, side="right") - 1)
    t = (s - edges[k]) / SEGMENTS[k][1]
    m1, m2 = COLUMNS
    name = SEGMENTS[k][0]

    def at(tt):
        if name == "top":
            return _sheet_path(m2, m1, tt, True, 0.55)
        if name == "down1":
            return _throat_path(m1, tt, True)
        if name == "bottom":
            return _sheet_path(m1, m2, tt, False, 0.55)
        return _throat_path(m2, tt, False)

    h = 1e-4
    xyz = at(t)
    tangent = at(min(t + h, 1.0)) - at(max(t - h, 0.0))
    tangent /= np.linalg.norm(tangent) + 1e-15
    # Each Alice handle reverses time orientation at its waist (t = 1/2 along the throat).
    if name == "top":
        orient = +1
    elif name == "bottom":
        orient = -1
    elif name == "down1":
        orient = +1 if t < 0.5 else -1
    else:
        orient = -1 if t < 0.5 else +1
    sheet = name if name in ("top", "bottom") else "throat"
    return Point(xyz, tangent, orient, sheet)


def orientation_flips(n: int = 4000) -> int:
    """How many times the proper-time orientation flips around one circuit."""
    o = np.array([worldline(i / n).orientation for i in range(n)])
    return int(np.count_nonzero(o != np.roll(o, 1)))


# ---- rendering ------------------------------------------------------------------------------

GOLD, CYAN, RED, BLUE, GREEN = "#f5c542", "#55d6f0", "#f0504a", "#4a7df0", "#3fbf7f"
SEAM, PAIR = "#f2f2f2", "#d9d9d9"
VIEW = np.array([np.cos(np.radians(13)) * np.cos(np.radians(45)),
                 np.cos(np.radians(13)) * np.sin(np.radians(45)), np.sin(np.radians(13))])


def _draw_static(ax):
    th = np.linspace(0, 2 * np.pi, 181)
    # sheets: rings about the sheet centre and radial lines
    for top in (True, False):
        for R in np.linspace(0.3, SHEET_RADIUS, 18):
            xy = np.column_stack([R * np.cos(th), R * np.sin(th)])
            ax.plot(xy[:, 0], xy[:, 1], sheet_height(xy, top), color=GREEN, lw=0.35, alpha=0.55)
        for phi in np.linspace(0, 2 * np.pi, 40, endpoint=False):
            R = np.linspace(0.05, SHEET_RADIUS, 80)
            xy = np.column_stack([R * np.cos(phi), R * np.sin(phi)])
            ax.plot(xy[:, 0], xy[:, 1], sheet_height(xy, top), color=GREEN, lw=0.25, alpha=0.35)
        # funnel rings around each mouth
        for c in COLUMNS:
            for rho in np.geomspace(RIM * 1.02, 0.95, 9):
                xy = c + rho * np.column_stack([np.cos(th), np.sin(th)])
                ax.plot(xy[:, 0], xy[:, 1], sheet_height(xy, top), color=GREEN, lw=0.4, alpha=0.6)
    # columns
    zr = HEIGHT / 2 - DEPTH
    z = np.linspace(-zr, zr, 120)
    for c in COLUMNS:
        for phi in np.linspace(0, 2 * np.pi, 28, endpoint=False):
            r = column_radius(z)
            ax.plot(c[0] + r * np.cos(phi), c[1] + r * np.sin(phi), z, color=GREEN, lw=0.35, alpha=0.5)
        for zz in np.linspace(-zr, zr, 9):
            r = column_radius(zz)
            ax.plot(c[0] + r * np.cos(th), c[1] + r * np.sin(th), np.full_like(th, zz), color=GREEN, lw=0.3, alpha=0.4)
    # mouths, coloured by apparent charge
    q = mouth_charges()
    for i, c in enumerate(COLUMNS, start=1):
        for top, zz in (("top", zr), ("bottom", -zr)):
            col = RED if q[(top, i)] > 0 else BLUE
            ax.plot(c[0] + RIM * np.cos(th), c[1] + RIM * np.sin(th), np.full_like(th, zz), color=col, lw=2.2)
            ax.text(c[0] + 1.25 * RIM, c[1] - 1.25 * RIM, zz + (0.18 if top == "top" else -0.25),
                    "+" if q[(top, i)] > 0 else "−", color=col, fontsize=12, weight="bold")
    # where the bottom mouths sit on the single surface: dashed images on the top sheet
    for i, c in enumerate(COLUMNS, start=1):
        p = identify(c)
        col = RED if q[("bottom", i)] > 0 else BLUE
        xy = p + RIM * np.column_stack([np.cos(th), np.sin(th)])
        ax.plot(xy[:, 0], xy[:, 1], sheet_height(xy, True) + 0.02, color=col, lw=1.7, ls=(0, (3, 2)), alpha=1.0)
    # sourceless flux: dashed lines on each sheet between its two mouths, and through each throat
    m1, m2 = COLUMNS
    for top in (True, False):
        for bulge in np.linspace(-1.3, 1.6, 11):
            t = np.linspace(0.02, 0.98, 120)
            xy = np.array([_sheet_path(m2, m1, tt, top, bulge)[:2] for tt in t])
            ax.plot(xy[:, 0], xy[:, 1], sheet_height(xy, top) + 0.01, color=GOLD, lw=0.7, ls=(0, (4, 3)), alpha=0.8)
    for c in COLUMNS:
        for phi in np.linspace(0, 2 * np.pi, 6, endpoint=False):
            r = 0.6 * column_radius(z)
            ax.plot(c[0] + r * np.cos(phi), c[1] + r * np.sin(phi), z, color=GOLD, lw=0.7, ls=(0, (4, 3)), alpha=0.75)
    # the worldline itself, gold where proper time runs with the circulation, cyan where reversed
    s = np.linspace(0, 1, 1601)
    W = [worldline(x) for x in s]
    P = np.array([w.xyz for w in W]) + 0.02 * np.r_[0, 0, 1]
    O = np.array([w.orientation for w in W])
    for sign, col in ((+1, GOLD), (-1, CYAN)):
        Q = P.copy()
        Q[O != sign] = np.nan
        ax.plot(Q[:, 0], Q[:, 1], Q[:, 2], color=col, lw=1.6, alpha=0.95)
    # the Alice seams: a band at each waist with a half-twist tick
    for c in COLUMNS:
        r = 1.18 * column_radius(0.0)
        ax.plot(c[0] + r * np.cos(th), c[1] + r * np.sin(th), np.zeros_like(th), color=SEAM, lw=1.3, alpha=0.95)
        u = np.linspace(-1, 1, 40)
        ang = np.arctan2(FRONT[1], FRONT[0]) + 0.9 * u
        ax.plot(c[0] + 1.5 * r * np.cos(ang), c[1] + 1.5 * r * np.sin(ang), 0.12 * np.tanh(2.5 * u),
                color=SEAM, lw=1.1, alpha=0.95)
    c2 = COLUMNS[1]
    ax.text(c2[0] + 0.75, c2[1] + 0.75, -0.05, "Alice seam\nt → −t,  q → −q", color=SEAM, fontsize=7.5)
    # one antipodal pairing: bottom mouth 1 is the top-sheet point −c1, normals ↓ and ↑
    c1 = COLUMNS[0]
    zb = -(HEIGHT / 2 - DEPTH)
    img = identify(c1)
    zi = sheet_height(img, True)[0]
    ax.plot([c1[0], img[0]], [c1[1], img[1]], [zb, zi], color=PAIR, lw=0.9, ls=(0, (1, 3)), alpha=0.7)
    ax.quiver(c1[0], c1[1], zb, 0, 0, -0.35, color=PAIR, lw=1.2, arrow_length_ratio=0.35)
    ax.quiver(img[0], img[1], zi, 0, 0, 0.35, color=PAIR, lw=1.2, arrow_length_ratio=0.35)
    ax.text(img[0] + 0.12, img[1] + 0.12, zi + 0.4, "same point", color=PAIR, fontsize=7.5)


def _arrow(ax, base, d, color):
    """A shaft with a chevron head drawn in the screen plane (3D quiver loses the head on vertical arrows)."""
    tip = base + d
    u = d / (np.linalg.norm(d) + 1e-15)
    right = np.cross(VIEW, [0.0, 0.0, 1.0])
    right /= np.linalg.norm(right)
    w = np.cross(VIEW, u)
    if np.linalg.norm(w) < 1e-6:
        w = right
    w /= np.linalg.norm(w)
    h = 0.42 * np.linalg.norm(d)
    lines = ax.plot(*np.array([base, tip]).T, color=color, lw=2.0)
    for sgn in (+1, -1):
        wing = tip - h * u + sgn * 0.55 * h * w
        lines += ax.plot(*np.array([tip, wing]).T, color=color, lw=2.0)
    return lines


def render(path_gif=None, path_png=None, frames: int = 54, balls: int = 6, dpi: int = 90):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import animation

    fig = plt.figure(figsize=(7.2, 8.6), facecolor="black")
    ax = fig.add_axes([-0.06, 0.09, 1.12, 0.95], projection="3d", facecolor="black")
    ax.set_axis_off()
    ax.view_init(elev=13, azim=45)
    lim = SHEET_RADIUS * 0.74
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-HEIGHT / 2 - 0.15, HEIGHT / 2 + 0.15)
    ax.set_box_aspect((1, 1, 1.05))
    _draw_static(ax)
    fig.text(0.5, 0.075, r"$\oint \mathbf{E}\cdot d\mathbf{A} = \pm 4\pi q$", color="white", ha="center", fontsize=20)
    fig.text(0.5, 0.035, "top ≡ bottom:  (x, y)↓ ~ (−x, −y)↑   ·   two Alice handles, four antipodal mouths",
             color="#cccccc", ha="center", fontsize=9)
    fig.text(0.5, 0.012, "gold: proper time along the circuit · cyan: reversed · white band: Alice seam · dotted: same point · schematic",
             color="#999999", ha="center", fontsize=7.5)

    dots = ax.scatter([], [], [], s=60, depthshade=False)
    ghosts = ax.scatter([], [], [], s=60, depthshade=False, facecolors="none", linewidths=1.4)
    arrows = []

    def update(f):
        nonlocal arrows
        for a in arrows:
            a.remove()
        arrows = []
        pts = [worldline(f / frames + j / balls) for j in range(balls)]
        xyz = np.array([p.xyz for p in pts])
        cols = [GOLD if p.orientation > 0 else CYAN for p in pts]
        dots._offsets3d = (xyz[:, 0], xyz[:, 1], xyz[:, 2] + 0.03)
        dots.set_color(cols)
        g = [p for p in pts if p.sheet == "bottom"]
        if g:
            gxy = np.array([identify(p.xyz[:2]) for p in g])
            gz = sheet_height(gxy, True) + 0.03
            ghosts._offsets3d = (gxy[:, 0], gxy[:, 1], gz)
            ghosts.set_edgecolor(CYAN)
        else:
            ghosts._offsets3d = ([], [], [])
        for p, c in zip(pts, cols):
            d = 0.32 * p.orientation * p.tangent
            arrows += _arrow(ax, p.xyz + np.array([0, 0, 0.03]) + 0.09 * VIEW, d, c)
        return [dots, ghosts, *arrows]

    if path_png is not None:
        update(int(frames * 0.12))
        fig.savefig(path_png, dpi=dpi * 2, facecolor="black")
    if path_gif is not None:
        anim = animation.FuncAnimation(fig, update, frames=frames, blit=False)
        anim.save(path_gif, writer=animation.PillowWriter(fps=15), dpi=dpi)
    plt.close(fig)


if __name__ == "__main__":
    from pathlib import Path
    out = Path(__file__).resolve().parents[2] / "docs" / "images" / "alice_circulation"
    out.mkdir(parents=True, exist_ok=True)
    render(out / "circulation.gif", out / "circulation.png")
    print("wrote", out)
