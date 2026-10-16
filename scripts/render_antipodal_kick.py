#!/usr/bin/env python
"""Render archived spectral states, never prescribe a focus kick or trajectory.

Requires the viz extra and a working VTK backend (EGL works headlessly).
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import pyvista as pv

from geometrodynamics.waves.antipodal_kick import AntipodalKick, harmonics


NAMES = ['symmetric', 'directional', 'dephased']
COLORS = ['#6edbd5', '#ffcc66', '#c5a1ff']
BG = '#101820'


def diagnostics(data, output):
    """Plot torque localization and the actual field/attachment momentum ledger."""
    names = NAMES+['mirror']
    traces = {n: np.genfromtxt(data/f'{n}.csv', delimiter=',', names=True) for n in names}
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), dpi=140, constrained_layout=True)
    for ax in axes.flat:
        ax.axvspan(np.pi-.6, np.pi+.6, color='#808080', alpha=.12)
        ax.axvline(np.pi, color='#666666', lw=.8, ls=':')
        ax.set_xlabel('Time (R/c)')
        ax.grid(alpha=.2)
    for name, color in zip(names, ['#148c87', '#b47a00', '#8553ba', '#287cb4']):
        tr = traces[name]
        axes[0, 0].plot(tr['t'], tr['torque_z'], color=color, label=name)
        axes[0, 1].plot(tr['t'], tr['angular_z'], color=color, label=name)
    axes[0, 0].set_title('Torque concentrates near refocus')
    axes[0, 0].set_ylabel('Torque about z')
    axes[0, 1].set_title('Reversing the incoming bias reverses recoil')
    axes[0, 1].set_ylabel('Attachment angular momentum Lz')
    axes[0, 1].legend(fontsize=8)
    tr = traces['directional']
    axes[1, 0].plot(tr['t'], tr['angular_z'], color='#b47a00', lw=2.5, label='Attachment gain')
    axes[1, 0].plot(tr['t'], tr['wave_angular_z'][0]-tr['wave_angular_z'],
                    color='#203d5b', ls='--', label='Field loss')
    axes[1, 0].set_title('Reciprocal momentum exchange (directional case)')
    axes[1, 0].set_ylabel('Change in angular momentum')
    axes[1, 0].legend(fontsize=8)
    for filename, label in [('directional.csv', 'dt <= 0.002'),
                             ('convergence_dt_half.csv', 'dt <= 0.001')]:
        tr = np.genfromtxt(data/filename, delimiter=',', names=True)
        axes[1, 1].plot(tr['t'], (tr['energy']/tr['energy'][0]-1)*1e6, label=label)
    axes[1, 1].set_title('Energy error under time refinement')
    axes[1, 1].set_ylabel('Relative energy error (parts per million)')
    axes[1, 1].legend(fontsize=8)
    fig.savefig(output/'diagnostics.png')
    plt.close(fig)


def render(data, output):
    output.mkdir(parents=True, exist_ok=True)
    archive = np.load(data/'render_states.npz')
    traces = {name: np.genfromtxt(data/f'{name}.csv', delimiter=',', names=True)
              for name in NAMES}
    model = AntipodalKick()
    base = pv.Sphere(theta_resolution=64, phi_resolution=40)
    points = base.points.astype(float)
    points /= np.linalg.norm(points, axis=1)[:, None]
    # Reference is the prescribed free target profile, fixed for every frame/case.
    ref = float(model.sample(model.target_state(.7)[0], [1., 0., 0.]))
    basis = harmonics(model.ell, model.order, points)
    fields = {name: (basis@archive[f'{name}_a'].T).real.T for name in NAMES}
    times = archive['directional_t']
    meshes, bubbles = [], []
    plotter = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(960, 340))
    plotter.set_background(BG)
    for j, name in enumerate(NAMES):
        plotter.subplot(0, j)
        mesh = base.copy()
        mesh.points = points.copy()
        mesh['u'] = fields[name][0]
        actor = plotter.add_mesh(mesh, scalars='u', clim=(-ref, ref), cmap='coolwarm',
                                 show_scalar_bar=False, smooth_shading=True)
        meshes.append(actor.mapper.dataset)
        bubble = pv.Sphere(radius=.10, center=(1.08, 0, 0), theta_resolution=24, phi_resolution=16)
        actor = plotter.add_mesh(bubble, color=COLORS[j], smooth_shading=True)
        bubbles.append(actor.mapper.dataset)
        ring_angle = np.linspace(0, 2*np.pi, 100)
        ring = np.column_stack([np.ones(100)*1.008, .14*np.cos(ring_angle), .14*np.sin(ring_angle)])
        plotter.add_mesh(pv.lines_from_points(ring), color='white', line_width=1)
        plotter.add_text(['Symmetric focus', 'Directional focus', 'Dephased control'][j],
                         position='upper_left', font_size=12, color=COLORS[j])
        plotter.camera_position = [(5., -3., 1.8), (0, 0, 0), (0, 0, 1)]
        plotter.camera.parallel_projection = True
        plotter.camera.parallel_scale = 1.35
    plotter.show(auto_close=False, interactive=False)
    # Curves use full archived output, not movie-decimated coefficients.
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 2.55), dpi=100, facecolor=BG)
    for ax in axes:
        ax.set_facecolor(BG)
        ax.tick_params(colors='#d6e1e8', labelsize=8)
        for spine in ax.spines.values():
            spine.set_color('#506070')
        ax.axvspan(np.pi-.6, np.pi+.6, color='white', alpha=.07)
        ax.axvline(np.pi, color='white', alpha=.4, lw=.8, ls=':')
        ax.set_xlim(0, times[-1])
        ax.set_xlabel('Time  (R/c)', color='#d6e1e8', fontsize=9)
        ax.grid(alpha=.1)
    axes[0].set_title('Attachment angular momentum Lz', color='white', fontsize=10)
    axes[1].set_title('Position along surface (degrees)', color='white', fontsize=10)
    for name, color in zip(NAMES, COLORS):
        tr = traces[name]
        axes[0].plot(tr['t'], tr['angular_z'], color=color, lw=1.8, label=name)
        axes[1].plot(tr['t'], np.degrees(np.arctan2(tr['x_y'], tr['x_x'])), color=color, lw=1.8)
    axes[0].legend(facecolor=BG, labelcolor='white', fontsize=7, loc='upper left', framealpha=.5)
    cursors = [ax.axvline(0, color='white', lw=1) for ax in axes]
    fig.tight_layout(pad=1.5)
    frames = []
    font = ImageFont.truetype('DejaVuSans.ttf', 16)
    small = ImageFont.truetype('DejaVuSans.ttf', 12)
    for i, time in enumerate(times):
        for j, name in enumerate(NAMES):
            a, x, q = [archive[f'{name}_{key}'][i] for key in ['a', 'x', 'q']]
            field = fields[name][i]
            meshes[j].points = points*(1+.18*np.tanh(field/ref))[:, None]
            meshes[j].compute_normals(cell_normals=False, point_normals=True, inplace=True)
            meshes[j].point_data['u'][:] = field
            meshes[j].Modified()
            radius = .10*(1+2*q)
            centre = x*(1+.18*np.tanh(float(model.sample(a, x))/ref)+.8*radius)
            fresh = pv.Sphere(radius=radius, center=centre, theta_resolution=24, phi_resolution=16)
            bubbles[j].points = fresh.points
        plotter.render()
        pic = Image.fromarray(plotter.screenshot(return_img=True)).convert('RGB')
        for cursor in cursors:
            cursor.set_xdata([time, time])
        fig.canvas.draw()
        chart = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3])
        canvas = Image.new('RGB', (960, 670), BG)
        draw = ImageDraw.Draw(canvas)
        draw.text((20, 12), f'ANTIPODAL REFOCUSING  |  t = {time:.2f}   |   focus at pi = 3.14',
                  font=font, fill='white')
        canvas.paste(pic, (0, 42))
        canvas.paste(chart, (0, 390))
        draw.text((20, 365), 'Solved scalar displacement; fixed radial display gain. Attachment position is not magnified.',
                  font=small, fill='#d6e1e8')
        draw.text((20, 647), 'Conservative scalar toy on S2. Shading = field sign; pale ring = original target. No GR or quantization claim.',
                  font=small, fill='#aabac8')
        frames.append(canvas)
    for label, target in [('incoming', np.pi-.35), ('focus', np.pi), ('after', times[-1])]:
        frames[int(np.argmin(abs(times-target)))].save(output/f'{label}.png', optimize=True)
    swatches = Image.new('RGB', (8*160, 8*112), BG)
    for i, frame in enumerate(frames):
        swatches.paste(frame.resize((160, 112)), ((i%8)*160, (i//8)*112))
    palette = swatches.quantize(colors=128)
    smaller = [frame.resize((768, 536)).quantize(palette=palette, dither=Image.Dither.NONE)
               for frame in frames]
    durations = [max(30, int(1500*(times[i+1]-t))) for i, t in enumerate(times[:-1])]+[1400]
    smaller[0].save(output/'antipodal_kick.gif', save_all=True, append_images=smaller[1:],
                    duration=durations, loop=0, optimize=True)
    plotter.close()
    plt.close(fig)
    diagnostics(data, output)
    print(f'Rendered {len(frames)} archived frames to {output}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path,
        default=Path('experiments/closure_ledger/runs/20261016_antipodal_kick'))
    parser.add_argument('--output', type=Path, default=Path('docs/images/antipodal_kick'))
    args = parser.parse_args()
    render(args.data, args.output)
