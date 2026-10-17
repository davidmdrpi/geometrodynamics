#!/usr/bin/env python
"""Render solved pair motion on #338's two-handle geometry, with fixed camera.

Scalar samples live on the round sphere; neutral necks have NO solved field.
Fibers follow the coordinate map and do not generate the motion.
"""
import argparse
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import pyvista as pv
from restore_moving_mouth_assets import load_asset

from geometrodynamics.viz.klein_handles import build_surface, is_orientable, euler_characteristic
from geometrodynamics.viz.moving_mouth_pairs import MovingPairGeometry, axis_rotation
from geometrodynamics.waves.antipodal_kick import Parameters
from geometrodynamics.waves.moving_mouth_pairs import MovingMouthPairs, harmonics

BG = '#101923'
COLORS = ['#55d9d1', '#ffc568']


def faces(tris):
    return np.column_stack([np.full(len(tris), 3), tris]).ravel()


def render(data, output):
    output.mkdir(parents=True, exist_ok=True)
    archive = np.load(io.BytesIO(load_asset(data/'render_states.npz')))
    summary = json.loads((data/'summary.json').read_text())
    model = MovingMouthPairs(Parameters(**summary['parameters']))
    surf = build_surface(n_sphere=1800, n_ring=32, ds_target=.008)
    geom = MovingPairGeometry(surf)
    # Fixed free-target amplitude scale for every case/frame.
    target = model.initial()
    model.drift(target, np.pi)
    reference = max(abs(model.sample(target.a, [[1., 0., 0.], [0., 0., 1.]])))
    times = archive['directional_t']
    cases = ['directional', 'symmetric']
    plotter = pv.Plotter(shape=(1, 2), off_screen=True, window_size=(1200, 450))
    plotter.set_background(BG)
    sphere_actors, tube_actors, rim_actors, fiber_actors = [], [], [], []
    glyph_colors = ['#55d9d1', '#ffc568', '#ad9cff', '#ef839a']
    for j, name in enumerate(cases):
        plotter.subplot(0, j)
        mesh = pv.PolyData(surf.points.copy(), faces(surf.tris[surf.sphere_tris]))
        mesh['u'] = np.zeros(len(surf.points))
        actor = plotter.add_mesh(mesh, scalars='u', clim=(-reference/4, reference/4), cmap='coolwarm',
                                 show_scalar_bar=False, smooth_shading=True, opacity=.82)
        sphere_actors.append(actor.mapper.dataset)
        tubes, rims, fibers = [], [], []
        for h, ids in enumerate(surf.tube_ids):
            tube = pv.PolyData(surf.points.copy(), faces(surf.tris[surf.tri_tube == h]))
            actor = plotter.add_mesh(tube, color='#8a969d', smooth_shading=True)
            tubes.append(actor.mapper.dataset)
            for rim in surf.rims[h]:
                reference_ring = surf.points[np.r_[rim, rim[0]]]*1.004
                plotter.add_mesh(pv.lines_from_points(reference_ring), color='white', opacity=.5, line_width=1)
                actor = plotter.add_mesh(pv.lines_from_points(reference_ring), color=COLORS[h], line_width=4)
                rims.append((actor.mapper.dataset, np.r_[rim, rim[0]]))
            for k, col in enumerate(range(0, surf.n_ring, 8)):
                path = ids[:, col]
                actor = plotter.add_mesh(pv.lines_from_points(surf.points[path]), color=glyph_colors[k], line_width=2)
                fibers.append((actor.mapper.dataset, path))
        tube_actors.append(tubes)
        rim_actors.append(rims)
        fiber_actors.append(fibers)
        # Fixed coordinate landmarks, not the moving mesh's material labels.
        plotter.add_point_labels(np.eye(3)*1.14, ['+x', '+y', '+z'], font_size=13,
                                  text_color='white', point_color='white', point_size=3,
                                  shape=None, always_visible=True)
        plotter.add_text('Directional packets' if j == 0 else 'Symmetric packets',
                         position='upper_left', font_size=15, color='white')
        plotter.camera_position = [(4.7, -6., 3.3), (0, 0, .2), (0, 0, 1)]
        plotter.camera.parallel_projection = True
        plotter.camera.parallel_scale = 1.86
    plotter.show(auto_close=False, interactive=False)
    trace = np.genfromtxt(data/'directional.csv', delimiter=',', names=True)
    projection = np.array([1., 0., 1.])/np.sqrt(2)
    fig, axes = plt.subplots(1, 3, figsize=(12, 2.6), dpi=100, facecolor=BG)
    for ax in axes:
        ax.set_facecolor(BG)
        ax.tick_params(colors='#d5e0ec', labelsize=8)
        for spine in ax.spines.values():
            spine.set_color('#52657b')
        ax.axvspan(0, .6, color='#5b8bbc', alpha=.14)
        ax.axvspan(np.pi-.6, np.pi+.6, color='white', alpha=.09)
        ax.axvline(np.pi, color='white', alpha=.4, lw=.7)
        ax.set_xlim(0, times[-1])
        ax.set_xlabel('Time (R/c)', color='white', fontsize=9)
        ax.grid(alpha=.1)
    momenta = []
    for h in range(2):
        x = np.column_stack([trace[f'x_{h}_{i}'] for i in range(3)])
        initial_axis = geom.axes[h]
        tilt = np.degrees(np.arctan2(np.linalg.norm(np.cross(initial_axis, x), axis=1), x@initial_axis))
        axes[0].plot(trace['t'], tilt, color=COLORS[h], label=f'Pair {h+1}')
        angular = np.column_stack([trace[f'angular_{h}_{i}'] for i in range(3)])@projection
        axes[1].plot(trace['t'], angular, color=COLORS[h], label=f'Pair {h+1}')
        momenta.append(angular)
    wave = np.column_stack([trace[f'wave_angular_{i}'] for i in range(3)])@projection
    axes[1].plot(trace['t'], sum(momenta), color='white', lw=2, label='Pair sum')
    axes[1].plot(trace['t'], wave[0]-wave, color='#ef839a', ls='--', label='Field loss')
    axes[2].plot(trace['t'], np.degrees(trace['separation'])-90, color='#ad9cff')
    axes[0].set_title('Magnified axis displacement (degrees)', color='white', fontsize=10)
    axes[1].set_title('Angular momentum along (x+z)/sqrt(2)', color='white', fontsize=10)
    axes[2].set_title('Change in pair separation (degrees)', color='white', fontsize=10)
    for ax in axes[:2]:
        ax.legend(fontsize=7, facecolor=BG, labelcolor='white', framealpha=.5)
    cursors = [ax.axvline(0, color='white', lw=.8) for ax in axes]
    fig.tight_layout(pad=1.6)
    font = ImageFont.truetype('DejaVuSans.ttf', 19)
    small = ImageFont.truetype('DejaVuSans.ttf', 14)
    frames = []
    checks = dict(minimum_sphere_orientation_ratio=1., maximum_rim_error=0.,
                  maximum_far_landmark_error=0., euler_characteristic=euler_characteristic(surf.tris),
                  orientable=is_orientable(surf.tris), frames=len(times), cases=cases,
                  fixed_field_scale=float(reference), display_gain=4., sphere_vertices=len(geom.sphere),
                  minimum_collar_gap=float('inf'),
                  note='Necks are neutral: no neck scalar or moving handled-metric solver is implied.')
    far = geom.sphere[np.all(geom.weights == 0, axis=1)]
    for i, time in enumerate(times):
        for j, name in enumerate(cases):
            axes_now = archive[f'{name}_x'][i]
            collar_gap = np.arccos(np.clip(abs(axes_now[0]@axes_now[1]), 0, 1))-2*.68
            checks['minimum_collar_gap'] = min(checks['minimum_collar_gap'], float(collar_gap))
            if collar_gap <= 0:
                raise ValueError('Moving drawing collars overlap')
            points = geom.coordinates(axes_now)
            ratio = geom.minimum_sphere_orientation_ratio(points)
            checks['minimum_sphere_orientation_ratio'] = min(checks['minimum_sphere_orientation_ratio'], ratio)
            if ratio <= 0:
                raise ValueError('Sphere map folded: refuse a misleading moving-mouth render')
            checks['maximum_far_landmark_error'] = max(checks['maximum_far_landmark_error'],
                    float(np.max(abs(points[far]-surf.points[far]))))
            for h, rims in enumerate(surf.rims):
                rotation = axis_rotation(geom.axes[h], axes_now[h])
                for rim in rims:
                    checks['maximum_rim_error'] = max(checks['maximum_rim_error'],
                            float(np.max(abs(points[rim]-surf.points[rim]@rotation.T))))
            # Important: Eulerian sampling at moved sphere points, NOT cached
            # values carried by moving mesh labels. No field is sampled on necks.
            field = (harmonics(model.ell, model.order, points[geom.sphere])@archive[f'{name}_a'][i]).real
            deformed = geom.deform(points, 4*field, reference)
            mesh = sphere_actors[j]
            mesh.points = deformed
            mesh.point_data['u'][geom.sphere] = field
            mesh.compute_normals(cell_normals=False, point_normals=True, inplace=True)
            mesh.Modified()
            for mesh in tube_actors[j]:
                mesh.points = deformed
                mesh.compute_normals(cell_normals=False, point_normals=True, inplace=True)
                mesh.Modified()
            for mesh, path in rim_actors[j]+fiber_actors[j]:
                mesh.points = deformed[path]
                mesh.Modified()
        plotter.render()
        picture = Image.fromarray(plotter.screenshot(return_img=True)).convert('RGB')
        for cursor in cursors:
            cursor.set_xdata([time, time])
        fig.canvas.draw()
        chart = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3])
        canvas = Image.new('RGB', (1200, 820), BG)
        draw = ImageDraw.Draw(canvas)
        draw.text((20, 12), f'TWO MOVING ANTIPODAL PAIRS  |  t = {time:.2f}  |  free refocus at pi', font=font, fill='white')
        canvas.paste(picture, (0, 45))
        draw.text((20, 495), 'Actual solved positions, fixed camera. White rims = initial locations. Fibers follow the moving necks.',
                  font=small, fill='#d5e0ec')
        draw.text((20, 517), 'Scalar deformation gain 4. Neutral necks: no computed throat wave or fluid flow. Mouth positions are unscaled.',
                  font=small, fill='#d5e0ec')
        canvas.paste(chart, (0, 542))
        draw.text((20, 798), 'Directional-case diagnostics. Blue window = launch; pale window = refocus. No action quantization or GR claim.',
                  font=small, fill='#aabacd')
        frames.append(canvas)
        if i % 10 == 0:
            print(f'Rendered {i+1}/{len(times)} frames', flush=True)
    for label, target in [('launch', 0), ('transit', np.pi/2), ('focus', np.pi), ('after', times[-1])]:
        frames[int(np.argmin(abs(times-target)))].save(output/f'{label}.png', optimize=True)
    swatches = Image.new('RGB', (8*180, 8*123), BG)
    for i, frame in enumerate(frames):
        swatches.paste(frame.resize((180, 123)), ((i%8)*180, (i//8)*123))
    palette = swatches.quantize(colors=128)
    smaller = [f.resize((960, 656)).quantize(palette=palette, dither=Image.Dither.NONE) for f in frames]
    durations = [max(30, int(1400*(t-times[i]))) for i, t in enumerate(times[1:])]+[1600]
    smaller[0].save(output/'moving_mouth_pairs.gif', save_all=True, append_images=smaller[1:],
                    loop=0, duration=durations, optimize=True)
    (output/'geometry_checks.json').write_text(json.dumps(checks, indent=2)+'\n')
    plotter.close()
    plt.close(fig)
    print(json.dumps(checks), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path,
                        default=Path('experiments/closure_ledger/runs/20261016_moving_mouth_pairs'))
    parser.add_argument('--output', type=Path, default=Path('docs/images/moving_mouth_pairs'))
    args = parser.parse_args()
    render(args.data, args.output)
