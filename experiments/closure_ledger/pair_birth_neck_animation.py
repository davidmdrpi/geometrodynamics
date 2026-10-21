#!/usr/bin/env python
"""Animation of the neck channel: a lopsided pair birth with the history clock.

    python experiments/closure_ledger/pair_birth_neck_animation.py [end] [frame_dt]

Defaults: end = 17 (birth, first flight, contact, release, refocusing, re-birth and the
start of the second flight), frame_dt = 0.05.

Panels:
- the sphere-side transverse field psi (B = psi_t) on an equal-area map with the birth
  circle as equator, with the mouths' rims and recent tracks; after contact, the released
  deformation pulse is drawn as contours;
- the neck field psi = eta_0 + eta_1 cos(theta) + eta_2 sin(theta), drawn on the catenoid
  unrolled at its true circumference, on the same colour scale;
- the clock offset A = tau_A - tau_B and the clock momentum mu;
- what the clock changes: separation with the history (and reversed) clock minus the
  clock-off neck, through the first flight.

Writes docs/images/pair_birth/neck_channel.mp4 (every frame, 1280 x 720, 20 frames per
unit time; needs imageio-ffmpeg) and a lighter GIF preview, neck_channel.gif (every fourth
frame, 720 px wide). Model: geometrodynamics/waves/pair_birth.py (NeckChannel, defaults).
About 5 minutes on 3 cores.
"""
import os
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path

for _var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_var, '1')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

from geometrodynamics.waves.moving_mouth_pairs import harmonics
from geometrodynamics.waves.pair_birth import NeckChannel, PairBirth, Parameters, birth_positions

IMAGES = Path(__file__).resolve().parents[2]/'docs'/'images'/'pair_birth'
PLANE = ((0., 0., 1.), (1., 0., 0.))
P, T = (np.array(e) for e in PLANE)
N = np.cross(P, T)
LON = np.linspace(-np.pi, np.pi, 121)
LAT = np.linspace(-np.pi/2, np.pi/2, 61)
A_RIM = .3


def lopsided():
    _, vel, _, _ = birth_positions(A_RIM, 1., PLANE)
    return np.array([.45*vel[0], .15*vel[1]])


def grid_points():
    lon, lat = np.meshgrid(LON, LAT)
    return (np.cos(lat)[..., None]*(np.cos(lon)[..., None]*P + np.sin(lon)[..., None]*T)
            + np.sin(lat)[..., None]*N).reshape(-1, 3)


def lonlat(x):
    x = np.asarray(x)
    return np.arctan2(x @ T, x @ P), np.arcsin(np.clip(x @ N, -1, 1))


def simulate(args):
    clock, end, frame_dt = args
    m = PairBirth(Parameters(), NeckChannel(clock=clock))
    frames = []
    if clock == 'history':
        pts = grid_points()
        Y = harmonics(m.ell, m.order, pts)
        uY = harmonics(m.u_ell, m.u_order, pts)

        def observe(s):
            psi = (Y @ (s.P.sum(0)/m.omega)).real
            deform = m.field_density(s, Y, uY)[1] if s.release is not None else np.zeros(len(pts))
            frames.append(dict(t=s.t, x=s.x.copy(), alive=s.alive, cycle=s.cycle, A=s.tau[0]-s.tau[1],
                               mu=s.mu, eta=s.eta.copy(), d=m.separation(s.x),
                               psi=psi.reshape(len(LAT), len(LON)).astype(np.float32),
                               deform=np.asarray(deform).reshape(len(LAT), len(LON)).astype(np.float32),
                               rebirth=None if s.release is None else s.release['time'] + np.pi))
    else:
        observe = None
    h, info = m.run_history(m.initial(velocities=lopsided()), end, cycles=2, dt=.004, output_dt=frame_dt,
                            observer=observe)
    contacts = [e['time'] for e in info['events'] if e['kind'] == 'contact']
    return clock, dict(t=h['t'], d=h['d'], alive=h['alive'], clock=h['clock'], mu=h['mu'], contacts=contacts,
                       neck=dict(length=m.neck_length, nodes=m.neck_nodes, b=np.sin(A_RIM)**2), frames=frames)


def rim_line(x, n=72):
    """Lon/lat of a mouth's rim (geodesic circle of radius a), broken where it wraps."""
    e1 = np.cross(x, N if abs(x @ N) < .9 else P)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(x, e1)
    th = np.linspace(0, 2*np.pi, n)
    ring = np.cos(A_RIM)*x + np.sin(A_RIM)*(np.cos(th)[:, None]*e1 + np.sin(th)[:, None]*e2)
    lon, lat = lonlat(ring)
    jump = np.r_[False, np.abs(np.diff(lon)) > np.pi]
    lon, lat = lon.astype(float), lat.astype(float)
    lon[jump], lat[jump] = np.nan, np.nan
    return lon, lat


def phase(f, contact):
    if f['alive']:
        return f"flight {f['cycle'] + 1}: mouths connected by the neck"
    return (f"rims touched at t = {contact:.2f}: neck removed, energy released as a pulse (contours); "
            f"re-birth at its antipodal focus, t = {f['rebirth']:.2f}")


def render(runs, path_gif, path_mp4=None):
    hist, off, rev = runs['history'], runs['off'], runs['reversed']
    frames = hist['frames']
    vmax = .6*float(np.percentile(np.abs(np.stack([f['psi'] for f in frames])), 99.5))    # extremes saturate
    neck = hist['neck']
    L, sigma, b = neck['length'], neck['nodes'], neck['b']
    radius = np.sqrt((sigma - L/2)**2 + b**2)
    theta = np.linspace(0, 2*np.pi, 49)
    tube_x = (theta[None, :] - np.pi)*radius[:, None]
    tube_y = np.broadcast_to(sigma[:, None], tube_x.shape)
    first = hist['contacts'][0]
    # separation differences through the first flight (common output times)
    diffs = {}
    for name, other in (('history', hist), ('reversed', rev)):
        n = min(np.argmax(~other['alive']), np.argmax(~off['alive']))      # first records after contact
        same = np.abs(other['t'][:n] - off['t'][:n]) < 1e-9
        diffs[name] = (other['t'][:n][same], (other['d'][:n] - off['d'][:n])[same])

    fig = plt.figure(figsize=(12.8, 7.2), dpi=100)
    gs = fig.add_gridspec(3, 2, width_ratios=[2.05, 1], height_ratios=[1.25, 1, 1], hspace=.6, wspace=.13,
                          left=.01, right=.95, top=.88, bottom=.08)
    ax_map = fig.add_subplot(gs[:, 0], projection='mollweide')
    ax_neck = fig.add_subplot(gs[0, 1])
    ax_clock = fig.add_subplot(gs[1, 1])
    ax_diff = fig.add_subplot(gs[2, 1])
    tw = ax_clock.twinx()

    alive = hist['alive']
    ax_clock.plot(hist['t'], np.where(alive, hist['clock'], np.nan), color='k', lw=1.2)
    tw.plot(hist['t'], np.where(alive, hist['mu'], np.nan), color='#c0392b', lw=1)
    ax_clock.set_ylabel('offset A', fontsize=8)
    tw.set_ylabel('μ', color='#c0392b', fontsize=8)
    tw.tick_params(axis='y', colors='#c0392b', labelsize=7)
    ax_clock.tick_params(labelsize=7)
    ax_clock.set_xlim(0, hist['t'][-1])
    ax_clock.set_title('history offset A = τ_A − τ_B (black) and clock momentum μ (red):\n'
                       'rest energies M + μ and M − μ', fontsize=8)
    cur_clock = ax_clock.axvline(0, color='#1f5fa8', lw=1)
    for name, col, lab in (('history', '#c0392b', 'history clock − clock off'),
                           ('reversed', '#e67e22', 'reversed clock − clock off')):
        ax_diff.plot(*diffs[name], color=col, lw=1.2, label=lab)
    ax_diff.axhline(0, color='k', lw=.4)
    ax_diff.set_xlim(0, hist['t'][-1])
    ax_diff.axvline(first, color='k', lw=.5, ls=':')
    ax_diff.text(first, ax_diff.get_ylim()[0], f'first contact {first - off["contacts"][0]:+.4f}\n'
                 f'vs clock off  ', fontsize=7, va='bottom', ha='right')
    ax_diff.legend(fontsize=7, loc='upper left')
    ax_diff.tick_params(labelsize=7)
    ax_diff.set_xlabel('time t', fontsize=8)
    ax_diff.set_title('what the clock changes: separation difference (first flight)', fontsize=8)
    cur_diff = ax_diff.axvline(0, color='#1f5fa8', lw=1)

    ax_map.grid(False)
    ax_map.set_xticklabels([])
    ax_map.set_yticklabels([])
    glon, glat = np.meshgrid(LON, LAT)
    mesh = ax_map.pcolormesh(glon, glat, frames[0]['psi'], cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='gouraud',
                             rasterized=True)
    ax_map.plot([0], [0], '+', color='k', ms=8)
    ax_map.text(0, -.12, 'P', ha='center', fontsize=9)
    ax_map.text(np.pi*.985, 0, '−P', ha='right', va='center', fontsize=9)
    ax_map.axhline(0, color='k', lw=.3, alpha=.4)
    cb = fig.colorbar(mesh, ax=ax_map, orientation='horizontal', fraction=.035, pad=.02, shrink=.55, extend='both')
    cb.set_label('transverse field ψ (B = ψ_t), sphere and neck', fontsize=8)
    cb.ax.tick_params(labelsize=7)
    colors = ('#117a3a', '#7d2fa5')
    rims = [ax_map.plot([], [], color=c, lw=2)[0] for c in colors]
    tracks = [ax_map.plot([], [], color=c, lw=1.2, alpha=.8)[0] for c in colors]
    labels = [ax_map.text(0, 0, k, fontsize=10, fontweight='bold', color=c, ha='center')
              for k, c in zip('AB', colors)]
    contours = []

    tube = ax_neck.pcolormesh(tube_x, tube_y, np.zeros(tube_x.shape), cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                              shading='gouraud', rasterized=True)
    ax_neck.plot(-np.pi*radius, sigma, color='k', lw=.8)
    ax_neck.plot(np.pi*radius, sigma, color='k', lw=.8)
    ax_neck.set_aspect('equal')
    ax_neck.set_xlim(-1.05*np.pi*np.sin(A_RIM), 1.05*np.pi*np.sin(A_RIM))
    ax_neck.set_ylim(-.04, L + .04)
    ax_neck.set_xticks([])
    ax_neck.set_yticks([0, L/2, L], ['mouth A', 'waist', 'mouth B'], fontsize=7)
    for tick, c in zip(ax_neck.get_yticklabels(), (colors[0], 'k', colors[1])):
        tick.set_color(c)
    neck_text = ax_neck.text(.01, .5, '', transform=ax_neck.transAxes, fontsize=7, va='center', ha='left')
    mu_text = ax_neck.text(.99, .5, '', transform=ax_neck.transAxes, fontsize=7, va='center', ha='right')
    ax_neck.set_title('neck interior, unrolled at its true circumference (#338 catenoid)', fontsize=8)
    title = fig.text(.5, .965, '', ha='center', fontsize=10.5)
    subtitle = fig.text(.5, .925, '', ha='center', fontsize=9, color='#444')

    canvas_frames = []
    trail = []
    for f in frames:
        mesh.set_array(f['psi'])
        for c in contours:
            c.remove()
        contours = []
        if f['deform'].max() > 0:
            dm = f['deform']
            levels = dm.max()*np.array([.05, .2, .5])
            lon, lat = np.meshgrid(LON, LAT)
            contours = [ax_map.contour(lon, lat, dm, levels=levels, colors='k', linewidths=[.5, .8, 1.2])]
        if f['alive']:
            trail.append(f['x'])
            trail = trail[-60:]
            for k in range(2):
                lon, lat = rim_line(f['x'][k])
                rims[k].set_data(lon, lat)
                tl, tb = lonlat(np.array([x[k] for x in trail]))
                tl = tl.astype(float)
                tb = tb.astype(float)
                jump = np.r_[False, np.abs(np.diff(tl)) > np.pi]
                tl[jump] = np.nan
                tracks[k].set_data(tl, tb)
                cl, cb_ = lonlat(f['x'][k])
                labels[k].set_position((cl, cb_ + .32))
                labels[k].set_visible(True)
            eta = f['eta']
            psi_neck = eta[0][:, None] + eta[1][:, None]*np.cos(theta)[None, :] + eta[2][:, None]*np.sin(theta)[None, :]
            tube.set_array(psi_neck)
            A = f['A']
            neck_text.set_text(f'offset A = {A:+.3f}\nA → B: L + A = {L + A:.3f}\nB → A: L − A = {L - A:.3f}')
            mu_text.set_text(f'clock momentum\nμ = {f["mu"]:+.1e}\nrest energies M ± μ')
        else:
            trail = []
            for k in range(2):
                rims[k].set_data([], [])
                tracks[k].set_data([], [])
                labels[k].set_visible(False)
            tube.set_array(np.zeros(tube_x.shape))
            neck_text.set_text('no handle\nuntil re-birth')
            mu_text.set_text('')
        cur_clock.set_xdata([f['t']])
        cur_diff.set_xdata([f['t']])
        title.set_text(f't = {f["t"]:5.2f}      lopsided pair birth (speeds 0.45 / 0.15) with the neck channel '
                       f'and the history clock' + (f'      separation d = {f["d"]:.3f}' if f['alive'] else ''))
        subtitle.set_text(phase(f, first))
        fig.canvas.draw()
        canvas_frames.append(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy())
    plt.close(fig)

    try:
        import imageio_ffmpeg
    except ImportError:
        imageio_ffmpeg = None
    if imageio_ffmpeg is not None and path_mp4 is not None:
        h, w = canvas_frames[0].shape[:2]
        writer = imageio_ffmpeg.write_frames(str(path_mp4), (w, h), fps=20, codec='libx264', pix_fmt_out='yuv420p',
                                             quality=8, macro_block_size=8)
        writer.send(None)
        for a in canvas_frames:
            writer.send(np.ascontiguousarray(a))
        writer.close()
        # GIF preview for the docs: every fourth frame at 720 px wide, with a palette fitted by ffmpeg.
        subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(), '-loglevel', 'error', '-y', '-i', str(path_mp4), '-vf',
                        "select='not(mod(n,4))',setpts=N/5/TB,scale=720:-1:flags=lanczos,split[a][b];"
                        "[a]palettegen=max_colors=128:stats_mode=diff[p];"
                        "[b][p]paletteuse=dither=bayer:bayer_scale=4:diff_mode=rectangle", '-r', '5', str(path_gif)],
                       check=True)
        return
    images = [Image.fromarray(a).resize((720, 405), Image.LANCZOS).convert('P', palette=Image.ADAPTIVE, colors=128)
              for a in canvas_frames[::4]]
    images[0].save(path_gif, save_all=True, append_images=images[1:], duration=200, loop=0, optimize=True)


def main(end=17., frame_dt=.05):
    with Pool(3) as pool:
        runs = dict(pool.map(simulate, [(c, end, frame_dt) for c in ('history', 'off', 'reversed')]))
    IMAGES.mkdir(parents=True, exist_ok=True)
    render(runs, IMAGES/'neck_channel.gif', IMAGES/'neck_channel.mp4')


if __name__ == '__main__':
    main(*map(float, sys.argv[1:3]))
