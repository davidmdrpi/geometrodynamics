#!/usr/bin/env python
"""Momentum-transfer movie for PR #339's antipodal-kick model.

Run from a checkout of PR #339 (needs geometrodynamics.waves.antipodal_kick):

    PYTHONPATH=. python render_momentum_transfer.py --output momentum_transfer.gif

The model is used unchanged, except the attachment mass, which is set to
M = 0.1 (PR default 0.3) so the solved motion is visible without magnifying
positions. Horizon 3*pi + 1.2: first focus at pi, coasting, next focus at 3*pi.
Field shading/height uses a fixed display gain. Bead positions, trails and
arrows come from the solved X(t) and L(t).
"""
import argparse
from dataclasses import replace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont

from geometrodynamics.waves.antipodal_kick import AntipodalKick, Parameters, harmonics

BG = "#0f1722"
CASES = [("Symmetric focus", 0.0, "#6edbd5"), ("Directional focus", 0.7, "#ffcc66"),
         ("Mirrored directional", -0.7, "#8ec5ff")]
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
BOLD = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"


def solve(mass, end, out_dt):
    model = AntipodalKick(replace(Parameters(), mass=mass))
    runs = []
    for _, eta, _ in CASES:
        arrays, stats = model.evolve(model.initial(eta), dt=0.002, end=end, output_dt=out_dt)
        runs.append(arrays)
        assert stats["max_energy_relative"] < 1e-4 and stats["max_angular_relative"] < 1e-8
    return model, runs


def phase(t):
    if t < np.pi - 0.6:
        return "incoming wave converging on the defect"
    if t < np.pi + 0.6:
        return "FOCUS: momentum exchange"
    if t < 3 * np.pi - 0.6:
        return "defect coasts with the momentum it received"
    if t < 3 * np.pi + 0.6:
        return "NEXT FOCUS (3π): second exchange"
    return "coasting"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--output", default="momentum_transfer.gif")
    ap.add_argument("--mass", type=float, default=0.1)
    ap.add_argument("--stride", type=int, default=2)
    args = ap.parse_args()
    end, out_dt = 3 * np.pi + 1.2, 0.05
    model, runs = solve(args.mass, end, out_dt)
    t = runs[1]["t"]
    ref = float(model.sample(model.target_state(0.7)[0], [1.0, 0.0, 0.0]))

    base = pv.Sphere(theta_resolution=72, phi_resolution=48)
    pts = base.points.astype(float)
    pts /= np.linalg.norm(pts, axis=1)[:, None]
    basis = harmonics(model.ell, model.order, pts)
    fields = [(basis @ r["a"].T).real.T for r in runs]

    W, HV, TOP = 1080, 420, 92
    pl = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(W, HV), border=False)
    pl.set_background(BG)
    meshes, beads, trails, arrows = [], [], [], []
    for j, (name, _, color) in enumerate(CASES):
        pl.subplot(0, j)
        mesh = base.copy()
        mesh.points = pts.copy()
        mesh["u"] = fields[j][0]
        actor = pl.add_mesh(mesh, scalars="u", clim=(-ref, ref), cmap="coolwarm", show_scalar_bar=False,
                            smooth_shading=True, specular=0.2)
        meshes.append(actor.mapper.dataset)
        ring = np.column_stack([np.full(100, 1.004), 0.09 * np.cos(np.linspace(0, 2 * np.pi, 100)),
                                0.09 * np.sin(np.linspace(0, 2 * np.pi, 100))])
        pl.add_mesh(pv.lines_from_points(ring), color="white", line_width=1.5)
        beads.append(color)
        pl.camera_position = [(5.2, 0.0, 1.1), (0.25, 0.0, 0.0), (0, 0, 1)]
        pl.camera.parallel_projection = True
        pl.camera.parallel_scale = 1.22
    pl.show(auto_close=False, interactive=False)

    # Charts: momentum ledger for the directional case and positions for all cases.
    fig, axes = plt.subplots(1, 2, figsize=(W / 100, 2.75), dpi=100, facecolor=BG)
    for ax in axes:
        ax.set_facecolor(BG)
        ax.tick_params(colors="#d6e1e8", labelsize=8)
        for sp in ax.spines.values():
            sp.set_color("#506070")
        for c in (np.pi, 3 * np.pi):
            ax.axvspan(c - 0.6, c + 0.6, color="white", alpha=0.07)
        ax.set_xlim(0, end)
        ax.set_xlabel("time (R/c)", color="#d6e1e8", fontsize=8)
        ax.grid(alpha=0.1)
    d = runs[1]
    Jf, La = d["wave_angular"][:, 2], d["angular"][:, 2]
    axes[0].plot(t, Jf, color="#ff8a80", lw=2, label="wave field")
    axes[0].plot(t, La, color="#ffcc66", lw=2, label="defect")
    axes[0].plot(t, Jf + La, color="white", lw=1.2, ls="--", label="total (conserved)")
    axes[0].set_title("Directional case: angular momentum moves from wave to defect", color="white", fontsize=9)
    axes[0].legend(facecolor=BG, labelcolor="white", fontsize=7, loc="center right", framealpha=0.6)
    for r, (name, _, color) in zip(runs, CASES):
        axes[1].plot(t, np.degrees(np.unwrap(np.arctan2(r["x"][:, 1], r["x"][:, 0]))), color=color, lw=2, label=name)
    axes[1].set_title("Defect position along the equator (degrees, not magnified)", color="white", fontsize=9)
    axes[1].legend(facecolor=BG, labelcolor="white", fontsize=7, loc="upper left", framealpha=0.6)
    cursors = [ax.axvline(0, color="white", lw=1) for ax in axes]
    fig.tight_layout(pad=1.2)

    f_title, f_text, f_small = (ImageFont.truetype(BOLD, 20), ImageFont.truetype(FONT, 15),
                                ImageFont.truetype(FONT, 12))
    lscale = 0.38 / max(np.abs(r["angular"]).max() for r in runs)
    frames = []
    for i in range(0, len(t), args.stride):
        for j, r in enumerate(runs):
            field = fields[j][i]
            meshes[j].points = pts * (1 + 0.16 * np.tanh(field / ref))[:, None]
            meshes[j].point_data["u"][:] = field
            meshes[j].Modified()
            pl.subplot(0, j)
            x = r["x"][i]
            lift = lambda p: p * (1 + 0.16 * np.tanh(float(model.sample(r["a"][i], p)) / ref) + 0.045)
            centre = lift(x)
            radius = 0.06 * (1 + 2 * float(r["q"][i]))
            pl.add_mesh(pv.Sphere(radius=radius, center=centre, theta_resolution=20, phi_resolution=14),
                        color=CASES[j][2], smooth_shading=True, name=f"bead{j}")
            path = np.array([p * 1.03 for p in r["x"][: i + 1]])
            if len(path) > 1 and np.linalg.norm(path[-1] - path[0]) > 1e-3:
                pl.add_mesh(pv.lines_from_points(path), color=CASES[j][2], line_width=4, name=f"trail{j}")
            else:
                pl.remove_actor(f"trail{j}")
            L = r["angular"][i]
            v = np.cross(L, x)
            if np.linalg.norm(L) * lscale > 0.04:
                direction = v / np.linalg.norm(v)
                length = np.linalg.norm(L) * lscale
                arrow = pv.Arrow(start=(0.0, 0.0, 0.0), direction=direction, scale=length,
                                 tip_length=0.35, tip_radius=0.14, shaft_radius=0.05)
                arrow.translate(centre + direction * radius, inplace=True)
                pl.add_mesh(arrow, color="white", name=f"arrow{j}")
            else:
                pl.remove_actor(f"arrow{j}")
        pl.render()
        view = Image.fromarray(pl.screenshot(return_img=True)).convert("RGB")
        for c in cursors:
            c.set_xdata([t[i], t[i]])
        fig.canvas.draw()
        chart = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3])
        canvas = Image.new("RGB", (W, TOP + HV + 20 + chart.height + 44), BG)
        dr = ImageDraw.Draw(canvas)
        dr.text((18, 10), "Momentum transfer at antipodal refocusing", font=f_title, fill="white")
        dr.text((18, 38), phase(t[i]), font=f_text, fill="#ffcc66" if "FOCUS" in phase(t[i]) else "#aab6c4")
        dr.text((W - 18, 12), f"t = {t[i]:5.2f}", font=f_title, fill="white", anchor="ra")
        dr.text((W - 18, 40), "focus at π and 3π", font=f_small, fill="#aab6c4", anchor="ra")
        canvas.paste(view, (0, TOP))
        for j, (name, _, color) in enumerate(CASES):
            dr.text((j * W // 3 + W // 6, 70), name, font=f_text, fill=color, anchor="ma")
        dr.text((18, TOP + HV + 2), "bead = defect at its solved position · line = its path · arrow = its momentum · ring = start",
                font=f_small, fill="#d6e1e8")
        canvas.paste(chart, (0, TOP + HV + 20))
        dr.text((18, canvas.height - 38), "PR #339 model unchanged except defect mass M = 0.1 (PR default 0.3) so motion is visible "
                "without magnifying positions.", font=f_small, fill="#aab6c4")
        dr.text((18, canvas.height - 20), "Conservative scalar on S² (−Δ + ¼). Field height/colour use a fixed display gain. "
                "No GR or quantization claim.", font=f_small, fill="#aab6c4")
        frames.append(canvas)
    sw = Image.new("RGB", (8 * 180, ((len(frames) + 7) // 8) * 140), BG)
    for k, f in enumerate(frames):
        sw.paste(f.resize((180, 140)), ((k % 8) * 180, (k // 8) * 140))
    pal = sw.quantize(colors=255)
    small = [f.resize((f.width * 9 // 10, f.height * 9 // 10)).quantize(palette=pal, dither=Image.Dither.NONE)
             for f in frames]
    durations = [140 if abs(t[min(k * args.stride, len(t) - 1)] - c) < 0.7 else 70
                 for k in range(len(small)) for c in [min((np.pi, 3 * np.pi), key=lambda c: abs(t[min(k * args.stride, len(t) - 1)] - c))]]
    durations[-1] = 1500
    small[0].save(args.output, save_all=True, append_images=small[1:], duration=durations, loop=0, optimize=True)
    for label, target in (("focus", np.pi + 0.3), ("coast", 2 * np.pi), ("end", end)):
        frames[int(np.argmin(abs(t[:: args.stride] - target)))].save(args.output.replace(".gif", f"_{label}.png"))
    pl.close()
    plt.close(fig)
    print("wrote", args.output, len(frames), "frames")


if __name__ == "__main__":
    main()
