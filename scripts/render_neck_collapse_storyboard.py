#!/usr/bin/env python
"""Storyboard of the neck-collapse protocol (docs/neck_collapse_protocol.md) on #340 geometry.

ILLUSTRATIVE, NOT A RESULT. The protocol is not frozen and no coupled run has been made,
so the waist trajectory r(t) is scripted. Everything else shown is computed:

- the field: exact free evolution (kappa = 0) of #340's symmetric x-packet, started
  pi + delta before its target focus (delta = 1.2, the protocol's start convention);
- the released pulse: exact free evolution of an impulsive deposit at the two former mouth
  sites, axisymmetric footprints with the baseline surface release 0.0233 (Stage 0);
- the neck shapes and V(r): the Stage-0 catenoid family at h/a = 0.60, Delta V = 0.01;
- the geometry: #338's Klein handles on #340's pair axes, drawn on a fixed sphere with the
  mouths cut out, field as colour only.

    PYTHONPATH=. xvfb-run -a python scripts/render_neck_collapse_storyboard.py
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'experiments/closure_ledger'))
import neck_stage0_check as stage0  # noqa: E402
from geometrodynamics.viz.klein_handles import build_surface  # noqa: E402
from geometrodynamics.waves.moving_mouth_pairs import MovingMouthPairs, harmonics  # noqa: E402

BG = '#101923'
ENERGY = LinearSegmentedColormap.from_list('energy', ['#26313d', '#3b1f5c', '#a2306b', '#f0703b', '#fbe36a'])
TEAL, AMBER, RED = '#55d9d1', '#ffc568', '#ef839a'
FONT = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
BOLD = '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'
DELTA, H, DV, R_S, E_RELEASE = 1.2, stage0.H_BASE, .01, stage0.R_SURGERY, .0233
T_SOURCE, T_TARGET = DELTA, np.pi+DELTA
RETRACT = .5


T_COLLAPSE = T_TARGET+.45   # scripted


def waist_ratio(t, rs_ratio):
    """SCRIPTED waist r(t)/r0: rest, ringing after the source focus, driven past the
    barrier at the target focus, then falling to the surgery radius at T_COLLAPSE."""
    if t <= T_SOURCE:
        return 1.
    ring = lambda u: 1-.10*np.sin(4*(u-T_SOURCE))
    ta = T_TARGET-.3
    if t <= ta:
        return ring(t)
    v0 = ring(ta)
    return v0+(rs_ratio-v0)*min((t-ta)/(T_COLLAPSE-ta), 1.)**2


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--output', type=Path, default=Path('docs/images/neck_collapse'))
    ap.add_argument('--dt', type=float, default=.06)
    ap.add_argument('--preview', type=int, nargs='*')
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    # Stage-0 neck: potential and shapes (computed).
    rb, r0 = stage0.catenoid_waists(H)
    a0, ab = stage0.area(r0, H), stage0.area(rb, H)
    gamma = DV/(ab-a0)
    rgrid = np.linspace(R_S, stage0.R_OUTER, 160)
    vgrid = gamma*(np.array([stage0.area(r, H) for r in rgrid])-a0)
    v_caps = gamma*(2*np.pi-a0)
    rb_ratio, rs_ratio = rb/r0, R_S/r0
    t_c = T_COLLAPSE

    # Field (computed, free evolution): #340's symmetric x-packet, own normalisation 0.1.
    m = MovingMouthPairs()
    f = m.ell*(m.ell+1)*np.exp(-m.par.sigma**2*m.ell*(m.ell+1)/2)*harmonics(m.ell, m.order, [1., 0., 0.]).conj()
    f = f*np.sqrt(m.par.energy/m.wave_energy(f, 0*f))

    def packet(t):  # a(T) = f, p(T) = 0 at the target focus T
        return f*np.cos(m.omega*(t-T_TARGET))

    def packet_p(t):
        return -m.omega*f*np.sin(m.omega*(t-T_TARGET))

    sites = np.array([[1., 0., 0.], [-1., 0., 0.]])
    w = m.smoothing*harmonics(m.ell, m.order, sites).conj()
    dp = w.sum(axis=0)
    dp *= np.sqrt(2*E_RELEASE)/np.linalg.norm(dp)      # (1/2)|dp|^2 = released energy

    def deposit(t):
        return np.zeros_like(dp) if t < t_c else dp*np.sin(m.omega*(t-t_c))/m.omega

    def deposit_p(t):
        return np.zeros_like(dp) if t < t_c else dp*np.cos(m.omega*(t-t_c))

    # Geometry: #338 handles on #340's initial pair axes (pair 1 at z, pair 2 at x).
    surf = build_surface(n_sphere=1800, n_ring=32, ds_target=.008)
    ids = surf.tube_ids[1]
    s_station = surf.s[ids[:, 0]]/surf.neck.S
    centre = surf.points[ids].mean(axis=1)
    bump = (1-s_station**2)**2
    faces = lambda tris: np.column_stack([np.full(len(tris), 3), tris]).ravel()
    tris2 = surf.tris[surf.tri_tube == 1]
    tube1 = pv.PolyData(surf.points.copy(), faces(surf.tris[surf.tri_tube == 0]))

    sphere = pv.Icosphere(radius=1., nsub=5)
    pts = np.asarray(sphere.points, float)
    pts /= np.linalg.norm(pts, axis=1)[:, None]
    sphere.points = pts
    basis = harmonics(m.ell, m.order, pts).astype(np.complex64)
    caps = [np.arccos(np.clip(pts@s, -1, 1)) < .3 for s in sites]
    area_w = 4*np.pi/len(pts)

    def energy_density(a, p):
        u = (basis@a.astype(np.complex64)).real
        ut = (basis@p.astype(np.complex64)).real
        g2 = sum(((basis@(d@a).astype(np.complex64)).real)**2 for d in m.generators)
        return .5*(ut**2+g2+u**2/4)

    def cap_energy(t):
        e = energy_density(deposit(t), deposit_p(t))
        return float(sum(e[c].sum() for c in caps)*area_w)

    e_scale = .5*float(energy_density(deposit(t_c+.02), deposit_p(t_c+.02)).max())
    e_field = .5*float(energy_density(packet(T_TARGET), packet_p(T_TARGET)).max())

    t_end = t_c+np.pi+.6
    times = np.unique(np.r_[np.arange(0, t_end, args.dt), T_SOURCE, T_TARGET, t_c, t_c+np.pi/2, t_c+np.pi, t_end])
    cap_e = np.array([cap_energy(t) for t in times])
    mouth_e = np.array([[float(energy_density(packet(t), packet_p(t))[c].sum()*area_w) for c in caps] for t in times])

    W, HV, TOP = 1200, 640, 78
    cells = [(row, col) for row in range(2) for col in range(2)]   # rows: field, released energy
    views = [(6., -3.3, 2.7), (-6., -3.3, 2.7)]                     # columns: toward +x, toward -x
    pl = pv.Plotter(shape=(2, 2), off_screen=True, window_size=(W, HV), border=False)
    pl.set_background(BG)
    for row, col in cells:
        pl.subplot(row, col)
        pl.add_mesh(tube1, color='#8a969d', smooth_shading=True)
        for rim in surf.rims[0]:
            pl.add_mesh(pv.lines_from_points(1.004*surf.points[np.r_[rim, rim[0]]]), color=TEAL, line_width=4)
        pl.add_point_labels(np.vstack([np.eye(3), -np.eye(3)[:1]])*1.16, ['+x', '+y', '+z', '−x'], font_size=12,
                            text_color='white', point_color='white', point_size=3, shape=None, always_visible=False)
        pl.camera_position = [views[col], (0., 0., .05), (0, 0, 1)]
        pl.camera.parallel_projection = True
        pl.camera.parallel_scale = 1.9
    pl.show(auto_close=False, interactive=False)

    fig, axes = plt.subplots(1, 4, figsize=(W/100, 2.75), dpi=100, facecolor=BG,
                             gridspec_kw=dict(width_ratios=[1.25, .8, 1.15, 1.15]))
    for ax in axes:
        ax.set_facecolor(BG)
        ax.tick_params(colors='#d5e0ec', labelsize=8)
        for sp in ax.spines.values():
            sp.set_color('#52657b')
        ax.grid(alpha=.1)
    ax = axes[0]
    ax.plot(rgrid/r0, vgrid, color='white', lw=1.6)
    ax.axhline(v_caps, color=RED, lw=1, ls='--')
    ax.text(.08, v_caps+.004, 'two closed caps', color=RED, fontsize=7)
    for x, lab in [(1, 'r₀'), (rb_ratio, 'barrier r_b'), (1/r0, 'r = a'), (rs_ratio, 'r_s')]:
        ax.axvline(x, color='#52657b', lw=.8)
        ax.text(x+.02, .135, lab, color='#aabacd', fontsize=7, rotation=90, va='top')
    marker, = ax.plot([1], [0], 'o', color=AMBER, ms=8)
    ax.set_xlim(0, stage0.R_OUTER/r0)
    ax.set_ylim(-.05, .14)
    ax.set_xlabel('waist r / r₀', color='white', fontsize=8)
    ax.set_title('Neck potential V(r) (Stage 0, ΔV = 0.01)', color='white', fontsize=9)
    ax = axes[1]
    zz = np.linspace(-H, H, 200)
    prof_hi, = ax.plot([], [], color=AMBER, lw=2)
    prof_lo, = ax.plot([], [], color=AMBER, lw=2)
    discs = [ax.plot([s*H, s*H], [-1, 1], color=AMBER, lw=3, alpha=0)[0] for s in (-1, 1)]
    ax.set_xlim(-.75, .75)
    ax.set_ylim(-1.45, 1.45)
    ax.set_aspect('equal')
    ax.set_xlabel('z (rims at ±h, h/a = 0.60)', color='white', fontsize=8)
    ax.set_title('Neck profile ρ(z)', color='white', fontsize=9)
    ax = axes[2]
    ax.plot(times, mouth_e[:, 0]*1e2, color=AMBER, lw=2, label='+x mouth cap')
    ax.plot(times, mouth_e[:, 1]*1e2, color=AMBER, lw=1.4, ls='--', label='−x mouth cap')
    for x, lab in [(T_SOURCE, 'source focus'), (T_TARGET, 'target focus'), (t_c, 't_c')]:
        ax.axvline(x, color='#52657b', lw=.8)
        ax.text(x+.08, mouth_e.max()*1e2*.35, lab, color='#aabacd', fontsize=7, rotation=90, va='bottom')
    ax.set_xlim(times[0], times[-1])
    ax.set_xlabel('time', color='white', fontsize=8)
    ax.set_title('Packet energy at pair 2\'s mouths (10⁻²)', color='white', fontsize=9)
    ax.legend(fontsize=7, facecolor=BG, labelcolor='white', framealpha=.5, loc='upper right')
    cursor_drive = ax.axvline(0, color='white', lw=.8)
    ax = axes[3]
    ax.plot(times-t_c, cap_e*1e3, color=TEAL, lw=2)
    ymax = cap_e.max()*1e3*1.15
    ax.set_ylim(0, ymax)
    for x, lab in [(0, 't_c'), (np.pi/2, 't_c + π/2'), (np.pi, 't_c + π')]:
        ax.axvline(x, color='#52657b', lw=.8)
        ax.text(x+.05, ymax*.9, lab, color='#aabacd', fontsize=7)
    ax.set_xlim(times[0]-t_c, times[-1]-t_c)
    ax.set_xlabel('time after collapse', color='white', fontsize=8)
    ax.set_title('Released pulse at former sites (10⁻³)', color='white', fontsize=9)
    cursor = ax.axvline(0, color='white', lw=.8)
    fig.tight_layout(pad=1.)

    f_title, f_text, f_small = ImageFont.truetype(BOLD, 20), ImageFont.truetype(FONT, 15), ImageFont.truetype(FONT, 12)
    frames, frame_t = [], []
    for i, t in enumerate(times):
        if args.preview and i not in args.preview:
            continue
        ratio = waist_ratio(t, rs_ratio) if t < t_c else rs_ratio
        pinch = t >= t_c
        s_cut = min((t-t_c)/RETRACT, 1.) if pinch else -1.
        closed = pinch and s_cut >= 1.
        for row, col in cells:
            pl.subplot(row, col)
            # The field domain is the full sphere (as in #340), so no mouths are cut out here.
            mesh = sphere.copy()
            # Energy density, not displacement: at a focal instant the energy can be all in u_t.
            if row == 0:
                mesh['e'] = energy_density(packet(t)+deposit(t), packet_p(t)+deposit_p(t))
                pl.add_mesh(mesh, scalars='e', clim=(0, e_field), cmap=ENERGY, show_scalar_bar=False,
                            smooth_shading=True, name='sphere')
            else:
                mesh['e'] = energy_density(deposit(t), deposit_p(t))
                pl.add_mesh(mesh, scalars='e', clim=(0, e_scale), cmap=ENERGY, show_scalar_bar=False,
                            smooth_shading=True, name='sphere')
            # Pair 2's neck: scripted waist, scaled about each ring's centroid; split after pinch.
            pts2 = surf.points.copy()
            factor = 1+(ratio-1)*bump
            pts2[ids] = centre[:, None]+factor[:, None, None]*(surf.points[ids]-centre[:, None])
            keep = tris2
            if pinch:
                ok = np.abs(surf.s/surf.neck.S) > s_cut
                keep = tris2[np.all(ok[tris2], axis=1)]
            if len(keep):
                pl.add_mesh(pv.PolyData(pts2, faces(keep)), color='#c9a76a', smooth_shading=True, opacity=.5,
                            name='tube2')
            else:
                pl.remove_actor('tube2')
            for q, rim in enumerate(surf.rims[1]):
                if closed:
                    pl.remove_actor(f'rim2{q}')
                else:
                    pl.add_mesh(pv.lines_from_points(1.004*surf.points[np.r_[rim, rim[0]]]), color=AMBER,
                                line_width=4, name=f'rim2{q}')
        pl.render()
        view = Image.fromarray(pl.screenshot(return_img=True)).convert('RGB')

        marker.set_data([ratio], [np.interp(ratio*r0, rgrid, vgrid)])
        marker.set_alpha(0 if pinch else 1)
        r = ratio*r0
        if pinch:
            prof_hi.set_data([], [])
            prof_lo.set_data([], [])
            for d in discs:
                d.set_alpha(1)
        else:
            if r <= 1:
                beta = np.arccosh(1/r)/H
                rho = r*np.cosh(beta*zz)
            else:
                beta = np.arccos(1/r)/H
                rho = r*np.cos(beta*zz)
            prof_hi.set_data(zz, rho)
            prof_lo.set_data(zz, -rho)
            for d in discs:
                d.set_alpha(0)
        cursor.set_xdata([t-t_c, t-t_c])
        cursor_drive.set_xdata([t, t])
        fig.canvas.draw()
        chart = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3])

        canvas = Image.new('RGB', (W, TOP+HV+24+chart.height+42), BG)
        dr = ImageDraw.Draw(canvas)
        if t < T_SOURCE-.6:
            phase, hot = 'symmetric packet converging on the source-side mouth (−x)', False
        elif t < T_SOURCE+.6:
            phase, hot = 'SOURCE FOCUS on pair 2\'s −x mouth: neck starts ringing', True
        elif t < T_TARGET-.6:
            phase, hot = 'packet in transit; neck rings about r₀', False
        elif t < t_c:
            phase, hot = 'TARGET FOCUS on +x: waist driven past the barrier', True
        elif t < t_c+RETRACT:
            phase, hot = 'PINCH (imposed surgery rule): halves retract, energy released', True
        elif abs(t-(t_c+np.pi)) < .6:
            phase, hot = 'released pulse RE-CONVERGES on the two former mouth sites (t_c + π)', True
        else:
            phase, hot = 'released ring travels; caps closed', False
        dr.text((18, 10), 'Neck-collapse protocol (#342): storyboard on #340 geometry', font=f_title, fill='white')
        dr.text((18, 38), phase, font=f_text, fill='#ffcc66' if hot else '#aabacd')
        dr.text((W-18, 12), f't = {t:5.2f}', font=f_title, fill='white', anchor='ra')
        dr.text((W-18, 40), 'ILLUSTRATIVE — waist r(t) scripted, not computed', font=f_small, fill=RED, anchor='ra')
        canvas.paste(view, (0, TOP))
        dr.text((18, TOP+4), 'looking at +x (target side)', font=f_text, fill='white')
        dr.text((W//2+18, TOP+4), 'looking at −x (source side)', font=f_text, fill='white')
        dr.text((18, TOP+28), 'field energy density: free packet + released pulse', font=f_small, fill='#aabacd')
        dr.text((18, TOP+HV//2+6), 'released pulse alone: energy density (difference field, own scale)', font=f_small,
                fill='#aabacd')
        dr.line([(W//2, TOP+8), (W//2, TOP+HV-8)], fill='#52657b', width=1)
        dr.line([(8, TOP+HV//2), (W-8, TOP+HV//2)], fill='#2b3a4a', width=1)
        dr.text((18, TOP+HV+4), 'Amber: pair 2 (±x), the collapsing neck · teal: pair 1 (±z), spectator · '
                'field on the full fixed sphere (the model\'s field domain)', font=f_small, fill='#d5e0ec')
        canvas.paste(chart, (0, TOP+HV+24))
        dr.text((18, canvas.height-36), 'Computed: free field (κ = 0) of #340\'s symmetric x-packet, started π + 1.2 before its '
                'target focus; released pulse = impulsive deposit of the Stage-0 surface release 0.0233', font=f_small,
                fill='#aabacd')
        dr.text((18, canvas.height-19), 'at the two former mouth sites; V(r) and ρ(z) from the Stage-0 catenoid family. '
                'The protocol is not frozen and no coupled run was made. No GR or quantization claim.', font=f_small,
                fill='#aabacd')
        frames.append(canvas)
        frame_t.append(t)
        if args.preview:
            canvas.save(args.output/f'preview_{i:03d}.png')
    pl.close()
    plt.close(fig)
    if args.preview:
        return
    frame_t = np.array(frame_t)
    for label, target in [('focus', T_TARGET), ('pinch', t_c+.25), ('reconvergence', t_c+np.pi)]:
        frames[int(np.argmin(abs(frame_t-target)))].save(args.output/f'{label}.png', optimize=True)
    sw = Image.new('RGB', (8*150, ((len(frames)+7)//8)*110), BG)
    for k, fr in enumerate(frames):
        sw.paste(fr.resize((150, 110)), ((k % 8)*150, (k//8)*110))
    pal = sw.quantize(colors=255)
    small = [fr.resize((fr.width*4//5, fr.height*4//5)).quantize(palette=pal, dither=Image.Dither.NONE) for fr in frames]
    hot = [abs(t-T_SOURCE) < .6 or abs(t-T_TARGET) < .6 or 0 <= t-t_c < RETRACT+.1 or abs(t-t_c-np.pi) < .4
           for t in frame_t]
    durations = [140 if h else 70 for h in hot]
    durations[-1] = 2000
    small[0].save(args.output/'storyboard.gif', save_all=True, append_images=small[1:], duration=durations,
                  loop=0, optimize=True)
    i_half, i_full = (int(np.argmin(abs(times-(t_c+x)))) for x in (np.pi/2, np.pi))
    print(f'collapse t_c = {t_c:.3f}; released-pulse cap energy at t_c+pi/2 = {cap_e[i_half]:.3e}, '
          f'at t_c+pi = {cap_e[i_full]:.3e} (ratio {cap_e[i_full]/cap_e[i_half]:.1f}); frames {len(frames)}')


if __name__ == '__main__':
    main()
