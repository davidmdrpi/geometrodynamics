"""The sphere with Alice handles, drawn from ``waves/wave_mouth_sphere``.

The 3D object is the one of ``viz/klein_handles`` (#338) moved by
``viz/moving_mouth_pairs`` (#340's relief): a unit sphere with two Klein-type handles.
Here every frame is a state of one solved Hamiltonian:

* sphere colour and radial relief: the solved real scalar, sampled at the
  moved sphere vertices;
* handle positions: the solved pair axes;
* rim colour: the solved throat coordinate of that mouth;
* the strip under the frames: the solved energy ledger.

What is NOT solved: nothing propagates along the drawn necks (they are grey
for that reason), the handle shape does not change when a throat tips, and
a Klein handle's passage through the sphere wall is a drawing artifact.
Renderer: matplotlib, painter-sorted, so no GL is needed.
"""
from __future__ import annotations

import numpy as np

from ..waves import wave_mouth_sphere as dynamics
from .klein_handles import HandleSpec, build_surface, default_handles
from .moving_mouth_pairs import MovingPairGeometry, axis_rotation

VIEW = np.array([0.55, -0.62, 0.56]) / np.linalg.norm([0.55, -0.62, 0.56])
COLLAR = 0.34                 # > mouth radius; 2*COLLAR < closest approach
FRAME_TIMES = (0.0, 3.0, 4.4, 6.0, 12.0, 30.0)
INK, MUTED, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"
LEDGER = (("wave", "scalar field on the sphere", "#2a78d6"),
          ("shape", "throat support and motion", "#eda100"),
          ("coupling", "coupling", "#9e9e9e"),
          ("kinetic", "mouth-pair kinetic", "#e34948"))


class SphereScene:
    """Builds the #338 handled sphere at the solved pair axes, frame by frame.

    #340's ``MovingPairGeometry.coordinates`` warps one fixed mesh and is only
    valid for small displacements; here a pair moves by more than a radian,
    so the surface is rebuilt with its handles at the solved axes. The rims
    are then exact at every frame. ``MovingPairGeometry.deform`` is still
    used for the radial relief and its taper into the necks.
    """

    def __init__(self, model, n_sphere=2600, n_ring=24):
        self.model = model
        self.n_sphere, self.n_ring = n_sphere, n_ring
        self.reference = default_handles(2)
        self.axes0 = np.array([h.axis for h in self.reference], float)

    def build(self, axes):
        handles = []
        for spec, start, end in zip(self.reference, self.axes0, axes):
            end = np.asarray(end, float) / np.linalg.norm(end)
            side = axis_rotation(start, end) @ np.asarray(spec.side, float)
            handles.append(HandleSpec(axis=tuple(end), side=tuple(side),
                                      klein=spec.klein))
        surface = build_surface(handles=handles, n_sphere=self.n_sphere,
                                n_ring=self.n_ring)
        return surface, MovingPairGeometry(surface, collar=COLLAR)

    def state(self, a, axes, scale, excursion=0.10):
        """Surface, drawn points and the field on sphere vertices for one state."""
        surface, geometry = self.build(axes)
        field = self.model.sample(a, surface.points[geometry.sphere])
        drawn = geometry.deform(surface.points, field, scale, excursion)
        return surface, geometry, drawn, field


def throat_colour(q, loaded, top):
    """Rim colour: open, tipped past the barrier (loaded only), or swinging."""
    if loaded and q < top:
        return "#e34948"
    return INK


def draw_state(ax, scene, a, axes, q, scale, view=VIEW):
    from matplotlib import cm, colors
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    model = scene.model
    surf, geometry, drawn, field = scene.state(a, axes, scale)
    values = np.zeros(len(drawn))
    values[geometry.sphere] = field
    sphere_tri = np.zeros(len(surf.tris), bool)
    sphere_tri[surf.sphere_tris] = True
    tris = surf.tris
    tri_pts = drawn[tris]
    centroid = tri_pts.mean(axis=1)
    normal = np.cross(tri_pts[:, 1] - tri_pts[:, 0], tri_pts[:, 2] - tri_pts[:, 0])
    normal /= np.linalg.norm(normal, axis=1)[:, None] + 1e-300
    outward = np.sign(np.einsum("ij,ij->i", normal, centroid))
    facing = np.einsum("ij,j->i", normal * outward[:, None], view)
    keep = ~sphere_tri | (centroid @ view > -0.05)
    norm = colors.Normalize(-scale, scale)
    face = cm.RdBu_r(norm(values[tris].mean(axis=1)))
    face[~sphere_tri] = colors.to_rgba("#b9b8b2")
    shade = 0.62 + 0.38 * np.clip(np.abs(facing), 0, 1)
    face[:, :3] *= shade[:, None]
    mesh = Poly3DCollection(tri_pts[keep], facecolors=face[keep],
                            edgecolors=face[keep], linewidths=0.2)
    ax.add_collection3d(mesh)
    top = model.barrier_top()
    for h, pair in enumerate(surf.rims):
        for k, rim in enumerate(pair):
            ring = drawn[np.r_[rim, rim[:1]]]
            if ring.mean(axis=0) @ view < -0.05:
                continue
            loaded = model.stored[h, k] > 0
            ax.plot(*(ring * 1.004).T, color=throat_colour(q[h, k], loaded, top),
                    lw=2.6 if loaded else 1.3, zorder=10)
    lim = 0.98
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim), zlim=(-lim, lim))
    ax.set_box_aspect((1, 1, 1))
    elev = np.degrees(np.arcsin(view[2]))
    azim = np.degrees(np.arctan2(view[1], view[0]))
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()


def nearest(times, t):
    return int(np.argmin(np.abs(np.asarray(times) - t)))


def field_scale(model, data):
    """One fixed colour and relief scale for every frame of a run."""
    pts = build_surface(n_sphere=400, n_ring=8).points[:300]
    pts = pts / np.linalg.norm(pts, axis=1)[:, None]
    peak = max(np.abs(model.sample(a, pts)).max() for a in data["a"][::10])
    return 0.5 * float(peak)


def summary(model, data):
    q, x = data["q"], data["x"]
    top = model.barrier_top()
    total = data["total"]
    tipped = []
    for h in range(2):
        for k in range(2):
            below = q[:, h, k] < top
            if model.stored[h, k] <= 0 or not below[-1]:
                tipped.append(None)
                continue
            above = np.flatnonzero(~below)
            tipped.append(float(data["t"][above[-1] + 1 if above.size else 0]))
    sep = np.arccos(np.clip(np.abs(np.sum(x[:, 0] * x[:, 1], axis=1)), 0, 1))
    return {
        "energy_initial": float(total[0]),
        "energy_max_relative_error": float(np.abs(total - total[0]).max()
                                           / max(abs(total[0]), 1e-300)),
        "angular_momentum_drift": float(data["angular_momentum_drift"]),
        "field_energy_initial": float(data["wave"][0]),
        "field_energy_final": float(data["wave"][-1]),
        "barrier_top": top,
        "transition_times": tipped,
        "min_coordinate": q.min(axis=0).ravel().tolist(),
        "final_coordinates": q[-1].ravel().tolist(),
        "pair_displacements": np.arccos(np.clip(
            np.sum(x[-1] * x[0], axis=1), -1, 1)).tolist(),
        "max_speed": float(model.speed(data["angular"]).max()),
        "min_pair_separation": float(sep.min()),
        "stored_released": float(np.sum(model.support(q[0]) - model.support(q[-1]))),
        # Share of the final field energy in the uniform l = 0 mode.
        "field_energy_final_l0": float(model.wave_energy_by_degree(
            data["a"][-1], data["p"][-1])[0]),
    }


def _ledger(ax, data, marks=(), now=None):
    t = data["t"]
    parts = {"wave": data["wave"], "shape": data["shape"].sum(axis=(1, 2)),
             "coupling": data["coupling"].sum(axis=(1, 2)),
             "kinetic": data["kinetic"].sum(axis=1)}
    ax.stackplot(t, *[parts[k] for k, _, _ in LEDGER],
                 colors=[c for _, _, c in LEDGER],
                 labels=[lab for _, lab, _ in LEDGER],
                 edgecolor=SURFACE, linewidth=0.6)
    ax.plot(t, data["total"], color=INK, lw=1.2, label="total")
    for m in marks:
        ax.axvline(m, color=MUTED, lw=0.7, ls=":")
    if now is not None:
        ax.axvline(now, color=INK, lw=1.4)
    ax.set_xlim(t[0], t[-1])
    ax.set_xlabel("time")
    ax.set_ylabel("energy")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def render(model, data, path, title=None, times=FRAME_TIMES, scene=None,
           scale=None):
    """Frames of the handled sphere over the solved ledger."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    scene = scene or SphereScene(model)
    scale = scale or field_scale(model, data) or 1.0
    s = summary(model, data)
    n = len(times)
    fig = plt.figure(figsize=(2.9 * n, 7.2), facecolor=SURFACE)
    grid = fig.add_gridspec(2, n, height_ratios=(1.35, 1.0), hspace=0.12,
                            wspace=0.0, left=0.05, right=0.99, top=0.86,
                            bottom=0.17)
    shown = []
    for i, target in enumerate(times):
        j = nearest(data["t"], target)
        shown.append(float(data["t"][j]))
        ax = fig.add_subplot(grid[0, i], projection="3d", facecolor=SURFACE)
        draw_state(ax, scene, data["a"][j], data["x"][j], data["q"][j], scale)
        ax.set_title(f"t = {data['t'][j]:.1f}", fontsize=10, color=INK, y=1.1)
    ax = fig.add_subplot(grid[1, :], facecolor=SURFACE)
    _ledger(ax, data, marks=shown)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=5,
              fontsize=8.5, frameon=False)
    ax.set_title(f"Ledger: total conserved to {s['energy_max_relative_error']:.1e}"
                 " relative; dotted lines are the frames above",
                 fontsize=10, loc="left", color=INK)
    fig.suptitle(title or "A solved wave on the sphere with Alice handles: it "
                 "moves the mouths and tips the loaded throat", fontsize=12,
                 color=INK, x=0.05, ha="left")
    fig.text(0.05, 0.012, "Colour and relief: the solved scalar on the round "
             "sphere. Thick rim: the loaded mouth, red once past its support "
             "barrier. Necks are grey: nothing is solved on them. Comparator "
             "model; support assumed; not a GR solution.", fontsize=8,
             color=MUTED)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return dict(s, frame_times=shown, field_scale=scale)


def animate(model, data, path, frames=72, fps=12, scene=None):
    """The same drawing as a GIF, with a moving cursor on the ledger."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    scene = scene or SphereScene(model, n_sphere=1800, n_ring=20)
    scale = field_scale(model, data)
    idx = np.unique(np.linspace(0, len(data["t"]) - 1, frames).round().astype(int))
    fig = plt.figure(figsize=(5.2, 6.4), facecolor=SURFACE)
    grid = fig.add_gridspec(2, 1, height_ratios=(2.6, 1.0), hspace=0.0,
                            left=0.13, right=0.97, top=0.97, bottom=0.1)
    ax3 = fig.add_subplot(grid[0], projection="3d", facecolor=SURFACE)
    axl = fig.add_subplot(grid[1], facecolor=SURFACE)

    def draw(j):
        ax3.cla()
        axl.cla()
        draw_state(ax3, scene, data["a"][j], data["x"][j], data["q"][j], scale)
        ax3.set_title(f"t = {data['t'][j]:.1f}", fontsize=10, color=INK, y=0.95)
        _ledger(axl, data, now=data["t"][j])

    anim = FuncAnimation(fig, draw, frames=list(idx))
    anim.save(path, writer=PillowWriter(fps=fps), dpi=80)
    plt.close(fig)
    return len(idx)
