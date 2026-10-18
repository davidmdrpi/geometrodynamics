#!/usr/bin/env python
"""Storyboard of the neck-collapse protocol (docs/neck_collapse_protocol.md) on #340 geometry.

ILLUSTRATIVE, NOT A RESULT. The protocol is not frozen and no coupled run has been made.
The field is uncoupled (kappa = 0): it cannot cause anything the neck does here.

Scripted or illustrative (drawing choices, not dynamics):
- the waist trajectory r(t): ringing, contraction, pinch time;
- the overview handle morph: #338's neck scaled about its ring centroids by a polynomial bump;
- the retraction: in the overview, the neck mesh is removed progressively from the middle;
  in the close-up, the two halves shorten into caps that end as the two-disc state.

Computed:
- the field: exact free evolution of #340's symmetric x-packet, started pi + delta before its
  target focus (delta = 1.2, the protocol's start convention);
- the released pulse: exact free evolution of a PRESCRIBED impulsive velocity pulse at the two
  former mouth sites, whose self-energy equals the Stage-0 surface release (0.0233). The
  storyboard does not enforce the protocol's surgery ledger: added to the packet, the pulse
  also carries interference energy Re<p_packet, dp>, printed at run time;
- the neck shapes in the close-up and V(r): the Stage-0 catenoid family at h/a = 0.60,
  Delta V = 0.01, evaluated at the scripted r(t); the end state is the two-disc solution;
- cap energies: Gauss-Legendre (polar) x uniform (azimuth) quadrature over 0.3-rad caps.

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
DELTA, H, DV, R_S, CAP = 1.2, stage0.H_BASE, .01, stage0.R_SURGERY, .3
T_SOURCE, T_TARGET = DELTA, np.pi+DELTA
T_COLLAPSE = T_TARGET+.45   # scripted
RETRACT = .5                # scripted
FIELD_LOOK = dict(ambient=.6, diffuse=.4, specular=0.)


def waist_ratio(t, rs_ratio):
    """SCRIPTED waist r(t)/r0: rest, ringing after the source focus, contraction from the target
    focus, reaching the surgery radius at T_COLLAPSE."""
    if t <= T_SOURCE:
        return 1.
    ring = lambda u: 1-.10*np.sin(4*(u-T_SOURCE))
    ta = T_TARGET-.3
    if t <= ta:
        return ring(t)
    v0 = ring(ta)
    return v0+(rs_ratio-v0)*min((t-ta)/(T_COLLAPSE-ta), 1.)**2


def profile(r, z):
    """Stage-0 family: r cosh(beta z) for r <= a, its barrel continuation for r > a."""
    if r <= 1:
        return r*np.cosh(np.arccosh(1/r)/H*z)
    return r*np.cos(np.arccos(1/r)/H*z)


def cap_quadrature(axis, radius=CAP, n_mu=24, n_phi=48):
    """Points and weights covering a spherical cap of given angular radius around axis."""
    x, wx = np.polynomial.legendre.leggauss(n_mu)
    lo = np.cos(radius)
    mu = lo+(1-lo)*(x+1)/2
    phi = 2*np.pi*np.arange(n_phi)/n_phi
    st = np.sqrt(1-mu**2)
    local = np.stack([np.outer(st, np.cos(phi)), np.outer(st, np.sin(phi)), np.outer(mu, np.ones(n_phi))], -1)
    axis = np.asarray(axis, float)
    e1 = np.cross(axis, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    pts = local[..., 0, None]*e1+local[..., 1, None]*e2+local[..., 2, None]*axis
    wts = np.outer(wx*(1-lo)/2, np.full(n_phi, 2*np.pi/n_phi))
    return pts.reshape(-1, 3), wts.ravel()


def revolution(zs, rhos, phi0, phi1, n_phi=48):
    """Surface of revolution about the x-axis (close-up frame): x = z, radius rho(z)."""
    phi = np.linspace(phi0, phi1, n_phi)
    zz, pp = np.meshgrid(zs, phi, indexing='ij')
    rr = np.broadcast_to(np.asarray(rhos)[:, None], zz.shape)
    grid = pv.StructuredGrid(zz, rr*np.cos(pp), rr*np.sin(pp))
    return grid.extract_surface(algorithm='dataset_surface')


def disc(x, radius, n=48):
    return pv.Disc(center=(x, 0, 0), inner=0., outer=max(radius, 1e-3), normal=(1, 0, 0), r_res=2, c_res=n)


def colorbar(draw, x, y, height, vmax, label, font):
    for k in range(height):
        rgb = ENERGY(1-k/(height-1))[:3]
        draw.line([(x, y+k), (x+10, y+k)], fill=tuple(int(255*c) for c in rgb))
    draw.rectangle([x, y, x+10, y+height], outline='#52657b')
    draw.text((x+14, y-2), f'{vmax:.2g}', font=font, fill='#d5e0ec')
    draw.text((x+14, y+height-12), '0', font=font, fill='#d5e0ec')
    draw.text((x, y-16), label, font=font, fill='#aabacd')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--output', type=Path, default=Path('docs/images/neck_collapse'))
    ap.add_argument('--dt', type=float, default=.06)
    ap.add_argument('--preview', type=int, nargs='*')
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    # Stage-0 neck (computed).
    rb, r0 = stage0.catenoid_waists(H)
    a0, ab = stage0.area(r0, H), stage0.area(rb, H)
    gamma = DV/(ab-a0)
    rgrid = np.linspace(R_S, stage0.R_OUTER, 160)
    vgrid = gamma*(np.array([stage0.area(r, H) for r in rgrid])-a0)
    v_caps = gamma*(2*np.pi-a0)
    release = gamma*(stage0.area(R_S, H)-2*np.pi)
    rb_ratio, rs_ratio, t_c = rb/r0, R_S/r0, T_COLLAPSE

    # Field (computed, free evolution): #340's symmetric x-packet, own normalisation 0.1.
    m = MovingMouthPairs()
    f = m.ell*(m.ell+1)*np.exp(-m.par.sigma**2*m.ell*(m.ell+1)/2)*harmonics(m.ell, m.order, [1., 0., 0.]).conj()
    f = f*np.sqrt(m.par.energy/m.wave_energy(f, 0*f))
    packet = lambda t: f*np.cos(m.omega*(t-T_TARGET))            # a(T) = f, p(T) = 0
    packet_p = lambda t: -m.omega*f*np.sin(m.omega*(t-T_TARGET))
    sites = np.array([[1., 0., 0.], [-1., 0., 0.]])
    dp = (m.smoothing*harmonics(m.ell, m.order, sites).conj()).sum(axis=0)
    dp *= np.sqrt(2*release)/np.linalg.norm(dp)                    # prescribed self-energy = release
    deposit = lambda t: np.zeros_like(dp) if t < t_c else dp*np.sin(m.omega*(t-t_c))/m.omega
    deposit_p = lambda t: np.zeros_like(dp) if t < t_c else dp*np.cos(m.omega*(t-t_c))
    interference = float(np.vdot(packet_p(t_c), dp).real)

    def density(basis, a, p):
        u = (basis@a.astype(np.complex64)).real
        ut = (basis@p.astype(np.complex64)).real
        g2 = sum(((basis@(d@a).astype(np.complex64)).real)**2 for d in m.generators)
        return .5*(ut**2+g2+u**2/4)

    quad = [cap_quadrature(s) for s in sites]
    qbasis = [harmonics(m.ell, m.order, q[0]).astype(np.complex64) for q in quad]
    cap = lambda a, p: [float(density(b, a, p)@q[1]) for b, q in zip(qbasis, quad)]

    t_end = t_c+np.pi+.6
    times = np.unique(np.r_[np.arange(0, t_end, args.dt), T_SOURCE, T_TARGET, t_c, t_c+np.pi/2, t_c+np.pi, t_end])
    pulse_cap = np.array([sum(cap(deposit(t), deposit_p(t))) for t in times])
    mouth_e = np.array([cap(packet(t), packet_p(t)) for t in times])
    i_half, i_full = (int(np.argmin(abs(times-(t_c+x)))) for x in (np.pi/2, np.pi))
    print(f'prescribed pulse self-energy {release:.10f}; interference with packet at t_c {interference:+.14f}; '
          f'field-energy jump {release+interference:.14f}')
    print(f'pulse energy in the two {CAP}-rad caps: t_c+pi/2 {pulse_cap[i_half]:.10e}, t_c+pi {pulse_cap[i_full]:.10e}, '
          f'ratio {pulse_cap[i_full]/pulse_cap[i_half]:.4f}')

    # Sphere display and signed meridian profile.
    sphere = pv.Icosphere(radius=1., nsub=5)
    pts = np.asarray(sphere.points, float)
    pts /= np.linalg.norm(pts, axis=1)[:, None]
    sphere.points = pts
    basis = harmonics(m.ell, m.order, pts).astype(np.complex64)
    theta = np.linspace(0, np.pi, 241)
    line_basis = harmonics(m.ell, m.order, np.stack([np.cos(theta), 0*theta, np.sin(theta)], -1)).astype(np.complex64)
    u_ref = float(np.abs(line_basis@f.astype(np.complex64)).max())
    e_field = .5*float(density(basis, packet(T_TARGET), packet_p(T_TARGET)).max())
    e_pulse = .5*float(density(basis, deposit(t_c+.02), deposit_p(t_c+.02)).max())

    # Geometry: #338 handles on #340's initial pair axes (pair 1 at z, pair 2 at x).
    surf = build_surface(n_sphere=1800, n_ring=32, ds_target=.008)
    ids = surf.tube_ids[1]
    s_station = surf.s[ids[:, 0]]/surf.neck.S
    centre = surf.points[ids].mean(axis=1)
    bump = (1-s_station**2)**2
    faces = lambda tris: np.column_stack([np.full(len(tris), 3), tris]).ravel()
    tris2 = surf.tris[surf.tri_tube == 1]
    tube1 = pv.PolyData(surf.points.copy(), faces(surf.tris[surf.tri_tube == 0]))

    W, HV, TOP = 1200, 600, 78
    pl = pv.Plotter(shape=(2, 3), groups=[(np.s_[:], 2)], off_screen=True, window_size=(W, HV), border=False)
    pl.set_background(BG)
    views = [(6., -3.3, 2.7), (-6., -3.3, 2.7)]
    for row in range(2):
        for col in range(2):
            pl.subplot(row, col)
            pl.add_mesh(tube1, color='#8a969d', smooth_shading=True)
            for rim in surf.rims[0]:
                pl.add_mesh(pv.lines_from_points(1.004*surf.points[np.r_[rim, rim[0]]]), color=TEAL, line_width=3)
            pl.add_point_labels(np.array([[1.16, 0, 0], [0, 1.16, 0], [-1.16, 0, 0]]), ['+x', '+y', '−x'],
                                font_size=11, text_color='white', point_color='white', point_size=3, shape=None)
            pl.camera_position = [views[col], (0., 0., .05), (0, 0, 1)]
            pl.camera.parallel_projection = True
            pl.camera.parallel_scale = 1.9
    pl.subplot(0, 2)
    for x in (-H, H):
        ring = np.stack([np.full(97, x), np.cos(np.linspace(0, 2*np.pi, 97)), np.sin(np.linspace(0, 2*np.pi, 97))], -1)
        pl.add_mesh(pv.lines_from_points(ring), color=AMBER, line_width=4)
    pl.add_mesh(pv.Line((-1.1, 0, 0), (1.1, 0, 0)), color='#52657b', line_width=1)

    def closeup_camera():
        pl.subplot(0, 2)
        pl.camera_position = [(1.0, -4.0, 1.6), (0., 0., 0.), (0, 0, 1)]
        pl.camera.parallel_projection = True
        pl.camera.parallel_scale = 1.45

    closeup_camera()
    pl.show(auto_close=False, interactive=False)

    fig, axes = plt.subplots(1, 4, figsize=(W/100, 2.75), dpi=100, facecolor=BG,
                             gridspec_kw=dict(width_ratios=[1.15, 1.2, 1.1, 1.1]))
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
    ax.set_xlabel('waist r / r₀ (marker: scripted)', color='white', fontsize=8)
    ax.set_title('Neck potential V(r), Stage 0', color='white', fontsize=9)
    ax = axes[1]
    line_total, = ax.plot([], [], color='white', lw=1.6, label='packet + pulse')
    line_pulse, = ax.plot([], [], color=TEAL, lw=1.4, label='pulse alone')
    ax.axhline(0, color='#52657b', lw=.8)
    ax.set_xlim(0, np.pi)
    ax.set_ylim(-1.05*u_ref, 1.05*u_ref)
    ax.set_xticks([0, np.pi/2, np.pi], ['+x', '+z', '−x'])
    ax.set_xlabel('meridian +x → +z → −x', color='white', fontsize=8)
    ax.set_title('Signed scalar u (display, not membrane shape)', color='white', fontsize=9)
    ax.legend(fontsize=7, facecolor=BG, labelcolor='white', framealpha=.5, loc='lower left')
    ax = axes[2]
    ax.plot(times, mouth_e[:, 0]*1e2, color=AMBER, lw=2, label='+x cap')
    ax.plot(times, mouth_e[:, 1]*1e2, color=AMBER, lw=1.4, ls='--', label='−x cap')
    for x, lab in [(T_SOURCE, 'source focus'), (T_TARGET, 'target focus'), (t_c, 't_c (scripted)')]:
        ax.axvline(x, color='#52657b', lw=.8)
        ax.text(x+.08, mouth_e.max()*1e2*.35, lab, color='#aabacd', fontsize=7, rotation=90, va='bottom')
    ax.set_xlim(times[0], times[-1])
    ax.set_xlabel('time', color='white', fontsize=8)
    ax.set_title('Packet energy at pair 2\'s mouths (10⁻²)', color='white', fontsize=9)
    ax.legend(fontsize=7, facecolor=BG, labelcolor='white', framealpha=.5, loc='upper right')
    cursor_drive = ax.axvline(0, color='white', lw=.8)
    ax = axes[3]
    ax.plot(times-t_c, pulse_cap*1e3, color=TEAL, lw=2)
    ymax = pulse_cap.max()*1e3*1.15
    ax.set_ylim(0, ymax)
    for x, lab in [(0, 't_c'), (np.pi/2, '+π/2'), (np.pi, '+π')]:
        ax.axvline(x, color='#52657b', lw=.8)
        ax.text(x+.05, ymax*.9, lab, color='#aabacd', fontsize=7)
    ax.set_xlim(times[0]-t_c, times[-1]-t_c)
    ax.set_xlabel('time after scripted collapse', color='white', fontsize=8)
    ax.set_title('Pulse alone at former sites (10⁻³)', color='white', fontsize=9)
    cursor = ax.axvline(0, color='white', lw=.8)
    fig.tight_layout(pad=1.)

    f_title, f_text = ImageFont.truetype(BOLD, 20), ImageFont.truetype(FONT, 15)
    f_small, f_tiny = ImageFont.truetype(FONT, 12), ImageFont.truetype(FONT, 11)
    zs = np.linspace(-H, H, 81)
    frames, frame_t = [], []
    for i, t in enumerate(times):
        if args.preview and i not in args.preview:
            continue
        ratio = waist_ratio(t, rs_ratio) if t < t_c else rs_ratio
        pinch = t >= t_c
        s_cut = min((t-t_c)/RETRACT, 1.) if pinch else -1.
        closed = pinch and s_cut >= 1.
        for row in range(2):
            for col in range(2):
                pl.subplot(row, col)
                mesh = sphere.copy()   # field domain: the full sphere, as in #340
                if row == 0:
                    mesh['e'] = density(basis, packet(t)+deposit(t), packet_p(t)+deposit_p(t))
                    clim = (0, e_field)
                else:
                    mesh['e'] = density(basis, deposit(t), deposit_p(t))
                    clim = (0, e_pulse)
                pl.add_mesh(mesh, scalars='e', clim=clim, cmap=ENERGY, show_scalar_bar=False,
                            smooth_shading=True, name='sphere', **FIELD_LOOK)
                pts2 = surf.points.copy()   # illustrative morph of #338's neck, not the catenoid
                factor = 1+(ratio-1)*bump
                pts2[ids] = centre[:, None]+factor[:, None, None]*(surf.points[ids]-centre[:, None])
                keep = tris2
                if pinch:   # illustrative retraction: progressive removal from the middle
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
                                    line_width=3, name=f'rim2{q}')
        # Close-up: Stage-0 shapes at the scripted waist; cutaway shows the inside.
        pl.subplot(0, 2)
        r = ratio*r0
        if not pinch:
            rho = np.array([profile(r, z) for z in zs])
            # cutaway at the plane y = 0: keep the half away from the camera (y >= 0) opaque
            far, near = revolution(zs, rho, -np.pi/2, np.pi/2), revolution(zs, rho, np.pi/2, 3*np.pi/2)
            pl.add_mesh(far, color='#c9a76a', smooth_shading=True, name='neck_far')
            pl.add_mesh(near, color='#c9a76a', smooth_shading=True, opacity=.12, name='neck_near')
            for nm, sgn in [('edge_a', 1), ('edge_b', -1)]:
                pl.add_mesh(pv.lines_from_points(np.stack([zs, 0*zs, sgn*rho], -1)), color=AMBER, line_width=4, name=nm)
            for nm in ('half_l', 'half_r', 'cap_l', 'cap_r'):
                pl.remove_actor(nm)
        else:
            for nm in ('neck_far', 'neck_near', 'edge_a', 'edge_b'):
                pl.remove_actor(nm)
            z_cut = H*s_cut   # illustrative retraction: halves shorten into caps
            for side, (half, cp) in [(1, ('half_r', 'cap_r')), (-1, ('half_l', 'cap_l'))]:
                if z_cut < H-1e-6:
                    zz = side*np.linspace(z_cut, H, 30)
                    rr = np.array([profile(r, z) for z in zz])
                    pl.add_mesh(revolution(zz, rr, 0, 2*np.pi), color='#c9a76a', smooth_shading=True, name=half)
                    pl.add_mesh(disc(side*z_cut, profile(r, z_cut)), color='#c9a76a', name=cp)
                else:
                    pl.remove_actor(half)
                    pl.add_mesh(disc(side*H, 1.), color='#c9a76a', opacity=.85, name=cp)   # two-disc end state
        closeup_camera()
        pl.render()
        view = Image.fromarray(pl.screenshot(return_img=True)).convert('RGB')

        marker.set_data([ratio], [np.interp(r, rgrid, vgrid)])
        marker.set_alpha(0 if pinch else 1)
        line_total.set_data(theta, (line_basis@(packet(t)+deposit(t)).astype(np.complex64)).real)
        line_pulse.set_data(theta, (line_basis@deposit(t).astype(np.complex64)).real)
        cursor.set_xdata([t-t_c, t-t_c])
        cursor_drive.set_xdata([t, t])
        fig.canvas.draw()
        chart = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3])

        canvas = Image.new('RGB', (W, TOP+HV+24+chart.height+62), BG)
        dr = ImageDraw.Draw(canvas)
        crossed = (not pinch) and ratio < rb_ratio
        if t < T_SOURCE-.6:
            phase = 'packet approaching the −x mouth (computed free field)'
        elif t < T_SOURCE+.6:
            phase = 'source focus at −x (computed) · neck: scripted ringing'
        elif t < T_TARGET-.6:
            phase = 'packet in transit (computed) · neck: scripted ringing about r₀'
        elif crossed:
            phase = 'scripted waist past the barrier r_b · contracting to r_s'
        elif not pinch:
            phase = 'target focus at +x (computed) · neck: scripted contraction'
        elif not closed:
            phase = 'pinch at r_s → two caps (imposed surgery rule; retraction drawn illustratively)'
        elif abs(t-(t_c+np.pi)) < .6:
            phase = 'prescribed pulse re-converges on both former sites (computed free field)'
        else:
            phase = 'prescribed pulse travelling (computed free field) · caps closed'
        hot = abs(t-T_SOURCE) < .6 or abs(t-T_TARGET) < .6 or crossed or (pinch and not closed) or \
            abs(t-(t_c+np.pi)) < .6
        dr.text((18, 10), 'Neck-collapse protocol (#342): storyboard on #340 geometry', font=f_title, fill='white')
        dr.text((18, 38), phase, font=f_text, fill='#ffcc66' if hot else '#aabacd')
        dr.text((W-18, 12), f't = {t:5.2f}', font=f_title, fill='white', anchor='ra')
        dr.text((W-18, 40), 'ILLUSTRATIVE: neck motion scripted, field uncoupled (κ = 0)', font=f_small, fill=RED,
                anchor='ra')
        canvas.paste(view, (0, TOP))
        dr.text((52, TOP+4), 'looking at +x (target side)', font=f_small, fill='white')
        dr.text((452, TOP+4), 'looking at −x (source side)', font=f_small, fill='white')
        dr.text((818, TOP+4), 'Neck close-up: Stage-0 shapes at the scripted waist', font=f_small, fill='white')
        dr.text((818, TOP+22), 'cutaway; after the pinch, retraction is illustrative', font=f_tiny, fill='#aabacd')
        colorbar(dr, 12, TOP+44, HV//2-80, e_field, 'field energy (packet + pulse)', f_tiny)
        colorbar(dr, 12, TOP+HV//2+44, HV//2-80, e_pulse, 'prescribed pulse alone', f_tiny)
        dr.line([(8, TOP+HV//2), (800, TOP+HV//2)], fill='#2b3a4a', width=1)
        dr.line([(800, TOP+8), (800, TOP+HV-8)], fill='#52657b', width=1)
        canvas.paste(chart, (0, TOP+HV+24))
        dr.text((18, TOP+HV+4), 'Amber: pair 2 (±x), the neck in the storyboard · teal: pair 1 (±z), spectator · '
                'colour: field energy density on the full fixed sphere', font=f_small, fill='#d5e0ec')
        foot = [
            'Scripted/illustrative: waist r(t) and pinch time; overview neck morph (polynomial bump about ring centroids); '
            'retraction (mesh removal, cap shortening).',
            'Computed: free field (κ = 0) of #340\'s symmetric x-packet; a PRESCRIBED velocity pulse at the former sites with '
            f'self-energy {release:.4f} (Stage-0 release); Stage-0 shapes and V(r).',
            f'Not enforced: the protocol\'s surgery ledger. Added to the packet, the pulse raises field energy by '
            f'{release+interference:.4f} at t_c (interference {interference:+.4f}). No coupled run; protocol not frozen.',
        ]
        for k, line in enumerate(foot):
            dr.text((18, canvas.height-58+k*18), line, font=f_tiny, fill='#aabacd')
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
    slow = [abs(t-T_SOURCE) < .6 or abs(t-T_TARGET) < .6 or 0 <= t-t_c < RETRACT+.1 or abs(t-t_c-np.pi) < .4
            for t in frame_t]
    durations = [150 if s else 70 for s in slow]
    durations[-1] = 2000
    small[0].save(args.output/'storyboard.gif', save_all=True, append_images=small[1:], duration=durations,
                  loop=0, optimize=True)
    print(f'frames {len(frames)}')


if __name__ == '__main__':
    main()
