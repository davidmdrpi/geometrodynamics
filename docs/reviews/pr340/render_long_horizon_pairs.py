#!/usr/bin/env python
"""Long-horizon movie for PR #340's two moving antipodal mouth pairs.

Run from a checkout of PR #340 (needs geometrodynamics.waves.moving_mouth_pairs):

    PYTHONPATH=. xvfb-run -a python render_long_horizon_pairs.py --output long_horizon

The PR's MovingMouthPairs model and parameters are used unchanged (M = 0.1); only
the horizon is extended from pi + 1.2 to 4 pi + 1.2. The symmetric run is stopped
when mouths of different handles reach the PR protocol's disjointness limit
(2*0.3 + 0.05 rad). Drawing differs from the PR renderer in two ways, both needed
for large motion: the sphere mesh is fixed in space and the moving mouths are cut
out of it each frame (the PR's collar warp folds the mesh beyond about 12 degrees), and the
scalar field is shown as colour only, with no radial display, so the only geometric
motion on screen is the solved motion of the mouths.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont

from geometrodynamics.viz.klein_handles import build_surface
from geometrodynamics.viz.moving_mouth_pairs import axis_rotation
from geometrodynamics.waves.moving_mouth_pairs import MovingMouthPairs, harmonics

BG = "#101923"
COLORS = ["#55d9d1", "#ffc568"]
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
BOLD = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
CASES = [("Directional packets", 0.7), ("Symmetric packets", 0.0)]


def faces(tris):
    return np.column_stack([np.full(len(tris), 3), tris]).ravel()


def solve(model, end, out_dt, cache):
    if cache.exists():
        data = np.load(cache)
        return [{k[len(f"c{j}_"):]: data[k] for k in data.files if k.startswith(f"c{j}_")}
                for j in range(len(CASES))]
    runs, store = [], {}
    for j, (_, eta) in enumerate(CASES):
        arrays, stats = model.evolve(model.initial(eta), dt=0.002, end=end, output_dt=out_dt)
        assert stats["max_energy_relative"] < 1e-4 and stats["max_angular_relative"] < 1e-8
        run = {k: arrays[k] for k in ["t", "a", "x", "angular", "wave_angular", "separation"]}
        runs.append(run)
        store.update({f"c{j}_{k}": v for k, v in run.items()})
        print(f"solved {CASES[j][0]}: energy err {stats['max_energy_relative']:.1e}, "
              f"J err {stats['max_angular_relative']:.1e}", flush=True)
    np.savez_compressed(cache, **store)
    return runs


def phase(t):
    k = int(np.round(t/np.pi))
    if k >= 1 and abs(t-k*np.pi) < 0.6:
        side = "+ mouths (target side)" if k % 2 else "− mouths (source side)"
        return f"FOCUS at {k}π on the {side}", True
    if t < 0.6:
        return "launch: packets leave the − mouths", True
    return "mouths coast with the momentum they took", False


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--output", type=Path, default=Path("long_horizon"))
    ap.add_argument("--out-dt", type=float, default=0.1)
    ap.add_argument("--preview", type=int, nargs="*", help="render only these frame indices as PNGs")
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    model = MovingMouthPairs()
    end = 4*np.pi+1.2
    runs = solve(model, end, args.out_dt, args.output/"long_runs.npz")
    t = runs[0]["t"]
    limit = 2*0.3+0.05

    # Symmetric case: stop at the protocol's disjoint-mouth limit.
    stops = []
    for r in runs:
        cross = np.minimum(r["separation"], np.pi-r["separation"])
        bad = np.flatnonzero(cross <= limit)
        stops.append(bad[0]-1 if len(bad) else len(t)-1)
    print("last valid frame per case:", [(s, round(float(t[s]), 2)) for s in stops], flush=True)

    # Fixed free-target amplitude scale, as in the PR renderer.
    target = model.initial()
    model.drift(target, np.pi)
    reference = max(abs(model.sample(target.a, [[1., 0., 0.], [0., 0., 1.]])))

    surf = build_surface(n_sphere=1800, n_ring=32, ds_target=.008)
    a_mouth = surf.a
    ref_axes = np.array([d.c for d in surf.drawn])
    assert np.allclose(ref_axes, runs[0]["x"][0])
    tube_faces = [faces(surf.tris[surf.tri_tube == h]) for h in range(2)]

    # Eulerian sphere: fixed points in space, field sampled there, moving holes cut out.
    sphere = pv.Icosphere(radius=1.0, nsub=5)
    pts = np.asarray(sphere.points, float)
    pts /= np.linalg.norm(pts, axis=1)[:, None]
    sphere.points = pts
    basis = harmonics(model.ell, model.order, pts).astype(np.complex64)

    W, HV = 1200, 470
    pl = pv.Plotter(shape=(1, 2), off_screen=True, window_size=(W, HV), border=False)
    pl.set_background(BG)
    circle = np.linspace(0, 2*np.pi, 97)
    for j, (name, _) in enumerate(CASES):
        pl.subplot(0, j)
        for h in range(2):
            for sgn in (1, -1):
                c = sgn*ref_axes[h]
                e1 = np.cross(c, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1)
                e2 = np.cross(c, e1)
                ring = 1.004*(np.cos(a_mouth)*c+np.sin(a_mouth)*(np.outer(np.cos(circle), e1)+np.outer(np.sin(circle), e2)))
                pl.add_mesh(pv.lines_from_points(ring), color="white", opacity=0.7, line_width=2)
        pl.add_point_labels(np.eye(3)*1.16, ["+x", "+y", "+z"], font_size=13, text_color="white",
                            point_color="white", point_size=3, shape=None, always_visible=True)
        # Look down the (x+z) diagonal, where both pairs' +mouths start and move face-on.
        pl.camera_position = [(5.0, -2.2, 5.0), (0.15, 0.0, 0.15), (0, 0, 1)]
        pl.camera.parallel_projection = True
        pl.camera.parallel_scale = 1.95
    pl.show(auto_close=False, interactive=False)

    # Charts.
    fig, axes = plt.subplots(1, 3, figsize=(W/100, 2.7), dpi=100, facecolor=BG)
    for ax in axes:
        ax.set_facecolor(BG)
        ax.tick_params(colors="#d5e0ec", labelsize=8)
        for sp in ax.spines.values():
            sp.set_color("#52657b")
        for k in range(5):
            ax.axvspan(k*np.pi-.6, k*np.pi+.6, color="white", alpha=.08, lw=0)
        ax.axvline(np.pi+1.2, color="#aabacd", lw=.8, ls=":")
        ax.set_xlim(0, end)
        ax.set_xticks([k*np.pi for k in range(5)], ["0", "π", "2π", "3π", "4π"])
        ax.grid(alpha=.1)
    for j, r in enumerate(runs):
        s = stops[j]+1
        ls = "-" if j == 0 else "--"
        x0 = r["x"][0]
        disp = np.degrees(np.arccos(np.clip(np.sum(r["x"]*x0[None], -1), -1, 1)))
        for h in range(2):
            axes[0].plot(t[:s], disp[:s, h], color=COLORS[h], ls=ls, lw=2,
                         label=f"pair {h+1}, {'directional' if j == 0 else 'symmetric'}")
            axes[2].plot(t[:s], np.linalg.norm(r["angular"][:s, h], axis=1), color=COLORS[h], ls=ls, lw=2)
        axes[1].plot(t[:s], np.degrees(r["separation"][:s])-90, color="#ad9cff", ls=ls, lw=2,
                     label="directional" if j == 0 else "symmetric")
    axes[1].axhline(np.degrees(limit)-90, color="#ef839a", lw=1)
    axes[1].text(0.2, np.degrees(limit)-90+2.5, "disjoint-mouth limit", color="#ef839a", fontsize=7)
    axes[0].text(np.pi+1.3, 1, "PR #340 horizon", color="#aabacd", fontsize=7, rotation=90, va="bottom")
    axes[0].set_title("Axis displacement (true degrees)", color="white", fontsize=10)
    axes[1].set_title("Change in angle between pair axes (degrees)", color="white", fontsize=10)
    axes[2].set_title("Pair angular momentum |L| (taken from field)", color="white", fontsize=10)
    axes[0].legend(fontsize=7, facecolor=BG, labelcolor="white", framealpha=.5, loc="upper left")
    axes[1].legend(fontsize=7, facecolor=BG, labelcolor="white", framealpha=.5, loc="center left")
    cursors = [ax.axvline(0, color="white", lw=.8) for ax in axes]
    fig.tight_layout(pad=1.2)

    f_title, f_text, f_small = (ImageFont.truetype(BOLD, 20), ImageFont.truetype(FONT, 15),
                                ImageFont.truetype(FONT, 12))
    TOP = 78
    frames = []
    for i in (args.preview if args.preview else range(len(t))):
        for j, r in enumerate(runs):
            k = min(i, stops[j])
            X = r["x"][k]
            pl.subplot(0, j)
            mouths = np.vstack([X, -X])
            d = np.arccos(np.clip(pts@mouths.T, -1, 1)).min(axis=1)
            mesh = sphere.copy()
            mesh["u"] = (basis@r["a"][k].astype(np.complex64)).real
            mesh["d"] = d
            cut = mesh.clip_scalar(scalars="d", value=a_mouth, invert=False)
            cut.set_active_scalars("u")
            pl.add_mesh(cut, scalars="u", clim=(-reference/4, reference/4), cmap="coolwarm",
                        show_scalar_bar=False, smooth_shading=True, opacity=.86, name=f"sphere{j}")
            for h in range(2):
                R = axis_rotation(ref_axes[h], X[h])
                moved = surf.points@R.T
                pl.add_mesh(pv.PolyData(moved, tube_faces[h]), color="#8a969d", smooth_shading=True,
                            name=f"tube{j}{h}")
                for q, rim in enumerate(surf.rims[h]):
                    pl.add_mesh(pv.lines_from_points(1.004*moved[np.r_[rim, rim[0]]]), color=COLORS[h],
                                line_width=5, name=f"rim{j}{h}{q}")
                for q, sgn in enumerate((1, -1)):
                    path = sgn*r["x"][:k+1, h]
                    path = 1.01*path[path@(sgn*X[h]) < np.cos(a_mouth)]  # stop the trail at the rim
                    if len(path) > 1 and np.linalg.norm(path[-1]-path[0]) > 1e-3:
                        pl.add_mesh(pv.lines_from_points(path), color=COLORS[h], line_width=6,
                                    name=f"trail{j}{h}{q}")
                    else:
                        pl.remove_actor(f"trail{j}{h}{q}")
        pl.render()
        view = Image.fromarray(pl.screenshot(return_img=True)).convert("RGB")
        for c in cursors:
            c.set_xdata([t[i], t[i]])
        fig.canvas.draw()
        chart = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3])
        canvas = Image.new("RGB", (W, TOP+HV+24+chart.height+42), BG)
        dr = ImageDraw.Draw(canvas)
        text, hot = phase(t[i])
        dr.text((18, 10), "Two antipodal mouth pairs moving across the sphere", font=f_title, fill="white")
        dr.text((18, 38), text, font=f_text, fill="#ffcc66" if hot else "#aabacd")
        dr.text((W-18, 12), f"t = {t[i]:5.2f}", font=f_title, fill="white", anchor="ra")
        dr.text((W-18, 40), "PR #340 model unchanged, run to 4π + 1.2", font=f_small, fill="#aabacd", anchor="ra")
        canvas.paste(view, (0, TOP))
        for j, (name, _) in enumerate(CASES):
            dr.text((j*W//2+18, TOP+4), name, font=f_text, fill="white")
            if i > stops[j]:
                dr.text((j*W//2+18, TOP+26), f"stopped at t = {t[stops[j]]:.2f}: mouths of the two",
                        font=f_small, fill="#ef839a")
                dr.text((j*W//2+18, TOP+42), "handles reached the 0.65 rad disjointness limit",
                        font=f_small, fill="#ef839a")
        dr.line([(W//2, TOP+8), (W//2, TOP+HV-8)], fill="#52657b", width=1)
        dr.text((18, TOP+HV+4), "Coloured rims = mouths at their solved positions · coloured lines = paths of the four "
                "mouth centres · white rings = starting rims", font=f_small, fill="#d5e0ec")
        canvas.paste(chart, (0, TOP+HV+24))
        dr.text((18, canvas.height-36), "Field drawn as colour only on a sphere fixed in space (no radial display): the "
                "only geometric motion shown is the mouths'. Necks are neutral drawings (#338), rotated rigidly with",
                font=f_small, fill="#aabacd")
        dr.text((18, canvas.height-19), "their pair axis; each Klein handle passes through the sphere wall only because "
                "it cannot be drawn in 3D without crossing. Round-S² scalar, M = 0.1. No GR or quantization claim.",
                font=f_small, fill="#aabacd")
        frames.append(canvas)
        if args.preview:
            canvas.save(args.output/f"preview_{i:03d}.png")
            continue
        if i % 20 == 0:
            print(f"frame {i+1}/{len(t)}", flush=True)

    if args.preview:
        return
    sw = Image.new("RGB", (8*150, ((len(frames)+7)//8)*104), BG)
    for k, f in enumerate(frames):
        sw.paste(f.resize((150, 104)), ((k % 8)*150, (k//8)*104))
    pal = sw.quantize(colors=255)
    small = [f.resize((f.width*4//5, f.height*4//5)).quantize(palette=pal, dither=Image.Dither.NONE) for f in frames]
    durations = [130 if phase(tt)[1] else 75 for tt in t]
    durations[-1] = 2000
    small[0].save(args.output/"long_horizon_pairs.gif", save_all=True, append_images=small[1:],
                  duration=durations, loop=0, optimize=True)
    for label, target in [("pr_horizon", np.pi+1.2), ("after_2pi", 2*np.pi+.6), ("after_3pi", 3*np.pi+.6),
                          ("end", end)]:
        frames[int(np.argmin(abs(t-target)))].save(args.output/f"{label}.png", optimize=True)
    pl.close()
    plt.close(fig)
    print("wrote", args.output/"long_horizon_pairs.gif", len(frames), "frames", flush=True)


if __name__ == "__main__":
    main()
