#!/usr/bin/env python
"""Tables and figures for docs/mouth_kick_experiments.md from the rim-flux probe runs.

J is the angular momentum the sphere band hands to each rim (rim momentum flux), not a
body impulse; "a" and "b" are the ambient directions it would move A and B if each rim's
flux acted on its mouth as a rigid body.

    python experiments/closure_ledger/mouth_kick_report.py
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

RUNS = Path(__file__).resolve().parent/'runs'/'20261019_mouth_kicks'
IMAGES = Path(__file__).resolve().parents[2]/'docs'/'images'/'mouth_kicks'
AXIS = np.array([1., 0., 0.])
T1, T2 = np.pi+.5, 2*np.pi+.5
KC, TC, GC = '#c0392b', '#1f5fa8', '#138d75'


def load(part):
    path = RUNS/f'{part}.json'
    if not path.exists():
        return []
    runs = json.loads(path.read_text())['runs']
    for r in runs:
        r['t'] = np.array(r['t'])
        J = np.array(r['J'])
        if r['mode'] == 'out':                          # A emitted the beam: its launch recoil goes to A
            J[:, 0] -= np.array(r['L0'])
        e1, e2 = np.array(r['e1']), np.array(r['e2'])
        r['a'], r['b'] = frame2(np.cross(J[:, 0], AXIS), e1, e2), frame2(-np.cross(J[:, 1], AXIS), e1, e2)
        if 'Jx' in r:
            Jx = np.array(r['Jx'])
            if r['mode'] == 'out':
                Jx[:, 0] -= np.array(r['Lx0'])
            r['ax'], r['bx'] = frame2(np.cross(Jx[:, 0], AXIS), e1, e2), frame2(-np.cross(Jx[:, 1], AXIS), e1, e2)
    return runs


def frame2(v, e1, e2):
    return np.stack([v @ e1, v @ e2], -1)


def pick(runs, **kw):
    out = [r for r in runs if all((abs(r[k]-v) < 1e-9) if isinstance(v, float) else r[k] == v for k, v in kw.items())]
    if len(out) != 1:
        raise KeyError(f'{kw}: {len(out)} matches')
    return out[0]


def at(r, key, T):
    return r[key][int(np.argmin(np.abs(r['t']-T)))]


def window(r, key, t0, t1):
    return at(r, key, t1)-at(r, key, t0)


PLATEAU = (3.0, 3.7)        # after the first passes over both mouths, before the second


def plateau(r, key):
    """Values of a rim-flux history inside the first-pass plateau window."""
    w = (r['t'] >= PLATEAU[0]-1e-9) & (r['t'] <= PLATEAU[1]+1e-9)
    return r[key][w]


def follow(a, b):
    return float(b @ a/(a @ a))


def side(psi):
    return np.array([np.cos(psi), np.sin(psi)])


def scan_table(scan):
    lines = ['| launch | handle | entry side | A rim flux | B rim flux along A | follow b·a/\\|a\\|² (first-pass plateau) | follow (one loop, T = 2π + ½) |',
             '|---|---|---|---|---|---|---|']
    for mode in ('release', 'out', 'in'):
        for h in ('klein', 'torus'):
            for p in (0., np.pi/4, np.pi/2):
                r = pick(scan, mode=mode, handle=h, psi=p)
                a, b = plateau(r, 'a').mean(0), plateau(r, 'b').mean(0)
                a2, b2 = at(r, 'a', T2), at(r, 'b', T2)
                lines.append(f"| {mode} | {h} | {np.degrees(p):.0f}° | {1e3*np.linalg.norm(a):.1f} | {1e3*(b @ a)/np.linalg.norm(a):+.1f} | "
                             f"{follow(a, b):+.3f} | {follow(a2, b2):+.3f} |")
    return '\n'.join(lines)


def throat_ratio(r):
    u = side(r['psi'])
    return float((window(r, 'b', 1., 2.2) @ u)/(window(r, 'a', 0., 1.) @ u))


def figure_exchange(scan, path):
    get = lambda h, m, p: pick(scan, handle=h, mode=m, psi=p)
    fig = plt.figure(figsize=(15, 9.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 3)
    ax = fig.add_subplot(gs[0, 0]); ax.set_aspect('equal'); ax.axis('off')
    ax.add_patch(Circle((0, 0), 1, fill=False, lw=1.5, color='0.3'))
    for y, lab in [(1, 'A'), (-1, 'B')]:
        ax.plot([-.3, .3], [y*.955, y*.955], color='k', lw=4); ax.text(0, y*1.13, f'mouth {lab}', ha='center', va='center', fontsize=11)
    th = np.linspace(np.pi/2-1.0, -np.pi/2+.3, 60)
    ax.plot(1.06*np.cos(th), 1.06*np.sin(th), color=GC, lw=2.5)
    ax.annotate('', xy=(1.06*np.cos(th[-1]), 1.06*np.sin(th[-1])), xytext=(1.06*np.cos(th[-3]), 1.06*np.sin(th[-3])),
                arrowprops=dict(arrowstyle='-|>', color=GC, lw=2))
    ax.text(1.15, 0, 'wave travels\nover the +u side', color=GC, fontsize=9.5, va='center')
    ax.plot(np.sin(1.0), np.cos(1.0), 'o', color=GC, ms=9); ax.text(np.sin(1.)+.06, np.cos(1.)+.05, 'release\n(side +u)', fontsize=9, color=GC)
    for y in (1, -1):
        ax.add_patch(FancyArrowPatch((0, y*.955), (-.55, y*.83), arrowstyle='-|>', mutation_scale=18, color=KC, lw=2.5))
    ax.text(-.62, .62, 'flux at A: −u', color=KC, ha='right', fontsize=10); ax.text(-.62, -.62, 'flux at B: −u', color=KC, ha='right', fontsize=10)
    ax.plot(-1, 0, '*', color=KC, ms=14); ax.text(-1.08, 0, 'common mode:\nwould close in here', ha='right', va='center', fontsize=9.5, color=KC)
    ax.set(xlim=(-1.9, 1.9), ylim=(-1.3, 1.3)); ax.set_title('(a) antipodal exchange on a closed sphere: the\nemitter- and receiver-side rim fluxes point the same way', fontsize=11)
    ax = fig.add_subplot(gs[0, 1:])
    for h, ls in [('klein', '-'), ('torus', '--')]:
        for p, lw in [(0., 2.2), (np.pi/2, 1.2)]:
            r = get(h, 'release', p); u = side(p)
            ax.plot(r['t'], -1e3*r['a'] @ u, color=KC, ls=ls, lw=lw, label=f'A, {h}, ψ₀ = {np.degrees(p):.0f}°')
            ax.plot(r['t'], -1e3*r['b'] @ u, color=TC, ls=ls, lw=lw, label=f'B, {h}, ψ₀ = {np.degrees(p):.0f}°')
    ax.axvspan(*PLATEAU, color='0.85', alpha=.6); ax.text(PLATEAU[0]+.03, 92, 'first-pass\nplateau', fontsize=9, color='0.4')
    ax.axhline(0, color='k', lw=.6)
    ax.set(xlabel='t', ylabel='cumulative rim flux toward −u (×10⁻³)',
           title='(b) one release next to A, nothing assigned: the rim fluxes at A and B point the same way')
    ax.legend(fontsize=7.5, ncol=2, loc='lower left', framealpha=.92); ax.grid(alpha=.3)
    ax = fig.add_subplot(gs[1, 0])
    ps = np.arange(8)*np.pi/4
    for h, c, mk in [('klein', KC, 'o'), ('torus', TC, 's')]:
        ax.plot(np.degrees(ps), [throat_ratio(get(h, 'in', p)) for p in ps], mk+'-', color=c, label=h)
    pp = np.linspace(0, 2*np.pi, 200)
    ax.plot(np.degrees(pp), throat_ratio(get('torus', 'in', 0.))*np.cos(2*pp), ':', color=KC, lw=1, label='Klein: m = 1 parity, ∝ cos 2ψ₀')
    ax.axhline(0, color='k', lw=.6)
    ax.set(xlabel='entry side ψ₀ (deg, from the handle mirror axis e₁)', ylabel="B's exit flux / A's entry flux (along −u)",
           title='(c) through the throat only: ~7.5% of the sideways flux;\nKlein sign set by entry side, torus always aligned')
    ax.set_xticks(range(0, 361, 90)); ax.legend(fontsize=8.5); ax.grid(alpha=.3)
    ax = fig.add_subplot(gs[1, 1:])
    labels, vals, cols = [], [], []
    groups = [('release', 'release next to A (both channels)'), ('out', 'A emits a beam over the sphere'), ('in', 'wave sent into A (throat channel)')]
    for m, name in groups:
        for h in ('klein', 'torus'):
            for p in (0., np.pi/4, np.pi/2):
                r = get(h, m, p); a, b = plateau(r, 'a').mean(0), plateau(r, 'b').mean(0)
                labels.append(f"{'K' if h == 'klein' else 'T'}{np.degrees(p):.0f}°"); vals.append(follow(a, b)); cols.append(KC if h == 'klein' else TC)
    x = np.arange(len(vals)); ax.bar(x, vals, color=cols)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8.5)
    for k, (m, name) in enumerate(groups):
        ax.text(6*k+2.5, 1.07, name, ha='center', fontsize=9.5)
        if k:
            ax.axvline(6*k-.5, color='0.7', lw=.6)
    ax.axhline(0, color='k', lw=.6); ax.set_ylim(-.3, 1.15)
    ax.text(17.4, -.27, 'K = Klein, T = torus; angle = entry side ψ₀', ha='right', fontsize=8.5, color='0.35')
    ax.set(ylabel="B's rim flux along A's, b·a/|a|²", title='(d) first-pass plateau: positive = the two rim fluxes point the same way; negative = opposite')
    ax.grid(alpha=.3, axis='y')
    fig.suptitle("Rim momentum flux from one classical wave exchange between antipodal mouths (#338's glued sphere, −Δ + ¼)\n"
                 "J = angular momentum leaving the sphere band at each rim: momentum into the throat, not a body impulse", fontsize=11.5)
    fig.savefig(path, dpi=100)


def flux_rows(flux):
    """Per scenario, over the first-pass plateau: quadratic and cross rim fluxes along A's quadratic flux.

    Unit wave energy and unit flux energy. The quadratic fluxes are steady in the window; the
    cross term is not, so its mean and range are both reported, and the follow values are
    given as the mean and range over the window.
    """
    rows = []
    for mode in ('release', 'out', 'in'):
        for h in ('klein', 'torus'):
            for p in (0., np.pi/4, np.pi/2):
                r = pick(flux, mode=mode, handle=h, psi=p)
                A, B, AX, BX = (plateau(r, k) for k in ('a', 'b', 'ax', 'bx'))
                ah = A.mean(0)/np.linalg.norm(A.mean(0))
                fp = np.array([follow(a+x, b+y) for a, b, x, y in zip(A, B, AX, BX)])
                fm = np.array([follow(a-x, b-y) for a, b, x, y in zip(A, B, AX, BX)])
                rows.append(dict(mode=mode, handle=h, psi=p, aq=float((A @ ah).mean()), bq=float((B @ ah).mean()),
                                 axa=float((AX @ ah).mean()), bxa=float((BX @ ah).mean()),
                                 bx_range=[float((BX @ ah).min()), float((BX @ ah).max())],
                                 follow_plus=float(fp.mean()), follow_plus_range=[float(fp.min()), float(fp.max())],
                                 follow_minus=float(fm.mean()), follow_minus_range=[float(fm.min()), float(fm.max())]))
    return rows


def flux_table(rows):
    lines = ["| launch | handle | side | quadratic A, B | cross A | cross B (range) | follow, signs aligned (range) | follow, signs opposed (range) |",
             '|---|---|---|---|---|---|---|---|']
    for r in rows:
        lines.append(f"| {r['mode']} | {r['handle']} | {np.degrees(r['psi']):.0f}° | {1e3*r['aq']:.1f}, {1e3*r['bq']:+.1f} | "
                     f"{1e3*r['axa']:+.1f} | {1e3*r['bxa']:+.1f} ({1e3*r['bx_range'][0]:+.0f} to {1e3*r['bx_range'][1]:+.0f}) | "
                     f"{r['follow_plus']:+.2f} ({r['follow_plus_range'][0]:+.2f} to {r['follow_plus_range'][1]:+.2f}) | "
                     f"{r['follow_minus']:+.2f} ({r['follow_minus_range'][0]:+.2f} to {r['follow_minus_range'][1]:+.2f}) |")
    return '\n'.join(lines)


def figure_flux(flux, path):
    rows = flux_rows(flux)
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.2), constrained_layout=True)
    ax = axs[0]
    r = pick(flux, mode='release', handle='klein', psi=np.pi/2); u = side(np.pi/2)
    for key, c, lab in [('a', KC, 'A, quadratic'), ('b', TC, 'B, quadratic'), ('ax', KC, 'A, flux cross term'), ('bx', TC, 'B, flux cross term')]:
        ax.plot(r['t'], -1e3*r[key] @ u, color=c, ls='-' if len(key) == 1 else '--', lw=2 if len(key) == 1 else 1.5, label=lab)
    ax.axvspan(*PLATEAU, color='0.85', alpha=.6); ax.axhline(0, color='k', lw=.6)
    ax.text(PLATEAU[0]+.03, -88, 'first-pass\nplateau', fontsize=8.5, color='0.35')
    ax.set(xlabel='t', ylabel='cumulative rim flux toward −u (×10⁻³)', title='(a) release next to A, Klein, ψ₀ = 90°, massless field:\nquadratic rim flux vs the term linear in the throat flux (unit energies)')
    ax.legend(fontsize=8.5); ax.grid(alpha=.3)
    ax = axs[1]
    x = np.arange(len(rows)); w = .38
    for off, key, col, lab in [(-w/2, 'follow_plus', GC, 'wave and flux signs aligned'), (w/2, 'follow_minus', '#8e44ad', 'signs opposed')]:
        mean = np.array([r[key] for r in rows]); rng = np.array([r[key+'_range'] for r in rows]).T
        ax.bar(x+off, mean, w, color=col, label=lab, yerr=[mean-rng[0], rng[1]-mean], capsize=2, error_kw=dict(lw=.8, ecolor='0.3'))
    ax.set_xticks(x); ax.set_xticklabels([f"{r['mode'][:3]} {'K' if r['handle'] == 'klein' else 'T'}{np.degrees(r['psi']):.0f}" for r in rows], rotation=60, fontsize=8)
    ax.axhline(0, color='k', lw=.6)
    ax.set(ylabel='follow = b·a/|a|², plateau mean (whiskers: range)', title='(b) equal wave and flux energy: the relative sign decides whether the rim fluxes align')
    ax.legend(fontsize=9); ax.grid(alpha=.3, axis='y')
    fig.savefig(path, dpi=100)
    return rows


NECK = 2*np.sin(.3)*np.cos(.3)            # geometric neck length 2S (null transit through the throat)
RETURN = np.pi-2*.3                        # rim-to-rim geodesic over the sphere
LOOP = NECK+RETURN                         # shortest null loop; the chronology horizon is offset = -LOOP


def offset_rows(offs):
    rows = []
    for r in offs:
        a, b = at(r, 'a', T1), at(r, 'b', T1)
        rows.append(dict(handle=r['handle'], psi=r['psi'], offset=r['offset'], tilt=r['tilt'], arrival_A=r['arrival_rim_A'],
                         arrival_B=r['arrival_rim_B'], follow=follow(a, b), approach=float(np.linalg.norm(a+b)), drift=r['energy_drift']))
    return sorted(rows, key=lambda r: (r['handle'], r['psi'], r['offset']))


def offset_table(rows):
    lines = ['| handle | side | offset Δ | clock tilt | rim A, first | rim B, first | B − A (threshold) | geometric transit 2S + Δ | follow, T = π + ½ |',
             '|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        lines.append(f"| {r['handle']} | {np.degrees(r['psi']):.0f}° | {r['offset']:+.2f} | {r['tilt']:.2f} | {r['arrival_A']:.3f} | "
                     f"{r['arrival_B']:+.3f} | {r['arrival_B']-r['arrival_A']:+.3f} | {NECK+r['offset']:+.3f} | {r['follow']:+.2f} |")
    return '\n'.join(lines)


def figure_offsets(offs, path):
    rows = offset_rows(offs)
    fig, axs = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    ax = axs[0]
    cmap = plt.cm.viridis
    sel = sorted([r for r in offs if r['handle'] == 'klein' and abs(r['psi']) < 1e-9], key=lambda r: r['offset'])
    for k, r in enumerate(sel):
        u = side(r['psi']); c = cmap(k/(len(sel)-1))
        ax.plot(r['t'], -1e3*r['b'] @ u, color=c, lw=1.8, label=f"Δ = {r['offset']:+.2f}")
    ax.axvline(0, color='k', lw=.6); ax.axhline(0, color='k', lw=.6)
    ax.set(xlabel='external time t', ylabel="B's cumulative rim flux toward −u (×10⁻³)", xlim=(-1.7, 4.5),
           title="(a) B's rim flux in external time, Klein, ψ₀ = 0:\nthe throat part moves earlier with Δ < 0")
    ax.legend(fontsize=8, ncol=2); ax.grid(alpha=.3)
    ax = axs[1]
    for h, p, c, mk in [('klein', 0., KC, 'o'), ('klein', np.pi/2, KC, 's'), ('torus', 0., TC, 'o'), ('torus', np.pi/2, TC, 's')]:
        rs = [r for r in rows if r['handle'] == h and abs(r['psi']-p) < 1e-9]
        ax.plot([r['offset'] for r in rs], [r['arrival_B'] for r in rs], mk+'-', color=c, label=f'rim B, {h} {np.degrees(p):.0f}°', ms=5)
        ax.plot([r['offset'] for r in rs], [r['arrival_A'] for r in rs], mk+':', color=c, alpha=.6, ms=4)
    r0 = [r for r in rows if r['handle'] == 'klein' and r['psi'] == 0 and r['offset'] == 0][0]
    dd = np.linspace(-2, 1, 10)
    ax.plot(dd, r0['arrival_B']+dd, 'k--', lw=.8, label='slope 1')
    ax.axhline(0, color='0.6', lw=.6)
    ax.set(xlabel='clock offset Δ (A → B)', ylabel='first external time the field passes threshold at the rim',
           title=f"(b) threshold arrival at rim B tracks Δ; rim A (dotted) does not.\nBelow Δ = −{r0['arrival_B']:.2f} it leaves B before the release at t = 0")
    ax.legend(fontsize=8); ax.grid(alpha=.3)
    ax = axs[2]
    rs = [r for r in rows if r['handle'] == 'klein' and r['psi'] == 0 and r['offset'] <= .5]
    d = np.array([r['offset'] for r in rs]); one = np.array([r['arrival_B']-r['arrival_A'] for r in rs])
    lag0 = float(one[d == 0][0])
    dd = np.linspace(-3.4, .6, 50)
    ax.plot(d, one, 'o', color=KC, label=f'measured one-way lag, threshold detection ({lag0:.3f} + Δ)')
    ax.plot(dd, NECK+dd, color='k', lw=1.2, label=f'geometric null transit of the neck, 2S + Δ (2S = {NECK:.3f})')
    ax.plot(dd, LOOP+dd, '--', color=TC, lw=1.5, label=f'geometric null loop, neck + return over the sphere ({LOOP:.3f} + Δ)')
    ax.axvline(-LOOP, color='k', lw=1, ls=':'); ax.text(-LOOP+.05, 2.0, f'chronology horizon\nΔ = −{LOOP:.3f}\n(null loop closes)', fontsize=8.5)
    ax.axhline(0, color='k', lw=.6)
    ax.set(xlabel='clock offset Δ', ylabel='external-time interval', xlim=(-3.4, .6),
           title='(c) measured packet lag vs geometric null times.\nOnly the geometric null loop is drawn to the horizon')
    ax.legend(fontsize=7.5, loc='upper left'); ax.grid(alpha=.3)
    fig.savefig(path, dpi=100)
    return rows


def ab_from(rec, e1, e2):
    J = np.array(rec['J'])
    return frame2(np.cross(J[:, 0], AXIS), e1, e2), frame2(-np.cross(J[:, 1], AXIS), e1, e2), np.array(rec['t'])


def offset_checks_summary(offs, ochecks):
    """Frame-spacing convergence, rim-weight independence and the direct comparison, per offset."""
    out = {}
    for c in ochecks:
        e1, e2 = np.array(c['e1']), np.array(c['e2'])
        main = pick(offs, handle=c['handle'], psi=c['psi'], offset=c['offset'])
        ref = [ab_from(rec, e1, e2) for rec in c['by_weight']]
        row = {}
        for T in (2.2, T1):
            i = lambda tt: int(np.argmin(np.abs(tt-T)))
            a0, b0, t0 = ref[0]
            scale = max(np.linalg.norm(a0[i(t0)]), np.linalg.norm(b0[i(t0)]))
            diff = lambda a, b, t: max(np.linalg.norm(a[i(t)]-a0[i(t0)]), np.linalg.norm(b[i(t)]-b0[i(t0)]))/scale
            row[f'T={T:.2f}'] = dict(frame_dt_0025_vs_00125=diff(main['a'], main['b'], main['t']),
                                     weights={str(tuple(w)): diff(*ref[k]) for k, w in enumerate(c['weights'][1:], 1)},
                                     direct_vs_frames=diff(*ab_from(c['direct'], e1, e2)) if 'direct' in c else None,
                                     scale=scale)
        out[f"offset={c['offset']:+.2f}"] = row
    return out


def checks_summary(checks, scan):
    out = {}
    if not checks:
        return out
    bare = pick(checks, tag='bare')
    out['bare_flux_over_L0'] = float(np.abs(np.array(bare['J'])).max()/np.linalg.norm(bare['L0']))
    ring = pick(checks, tag='ring'); rel = pick(scan, mode='release', handle='klein', psi=0.)
    out['ring_over_release'] = float(max(np.linalg.norm(at(ring, 'a', T1)), np.linalg.norm(at(ring, 'b', T1)))/np.linalg.norm(at(rel, 'a', T1)))
    ww = pick(checks, tag='wide_weight'); base = pick(scan, mode='release', handle='klein', psi=np.pi/4)
    out['weight_dependence'] = float(max(np.linalg.norm(at(ww, k, T1)-at(base, k, T1)) for k in 'ab')/np.linalg.norm(at(base, 'a', T1)))
    sf = pick(checks, tag='sign_flip'); b2 = pick(scan, mode='release', handle='klein', psi=np.pi/2)
    out['sign_flip_max_difference'] = float(max(np.abs(at(sf, k, T1)-at(b2, k, T1)).max() for k in 'ab'))
    mesh = []
    for r in [c for c in checks if c['tag'] == 'mesh48k']:
        b = pick(scan, mode=r['mode'], handle=r['handle'], psi=r['psi'])
        mesh.append(max(np.linalg.norm(at(r, k, T1)-at(b, k, T1)) for k in 'ab')/np.linalg.norm(at(b, 'a', T1)))
    out['mesh_24k_vs_48k'] = float(max(mesh))
    out['throat_ratio_mesh48k'] = {f"{r['handle']}_{np.degrees(r['psi']):.0f}": throat_ratio(r) for r in checks if r['tag'] == 'mesh48k' and r['mode'] == 'in'}
    out['throat_ratio_wide_packet'] = {f"{r['handle']}_{np.degrees(r['psi']):.0f}": float((window(r, 'b', 1.2, 2.6) @ side(r['psi']))/(window(r, 'a', 0., 1.2) @ side(r['psi'])))
                                       for r in checks if r['tag'] == 'wide_packet'}
    ml = []
    for r in [c for c in checks if c['tag'] == 'massless']:
        b = pick(scan, mode=r['mode'], handle=r['handle'], psi=r['psi'])
        ml.append(max(np.linalg.norm(at(r, k, T1)-at(b, k, T1)) for k in 'ab')/np.linalg.norm(at(b, 'a', T1)))
    out['massless_vs_conformal'] = float(max(ml))
    out['max_energy_drift'] = float(max(r['energy_drift'] for r in checks+scan if 'energy_drift' in r))
    return out


def main():
    IMAGES.mkdir(parents=True, exist_ok=True)
    scan, flux, offs, checks, ochecks = load('scan'), load('flux'), load('offsets'), load('checks'), load('offset_checks')
    report = {}
    if scan:
        figure_exchange(scan, IMAGES/'exchange_kicks.png')
        print(scan_table(scan), '\n')
        report['throat_ratio'] = {f"{h}_{np.degrees(p):.0f}": throat_ratio(pick(scan, mode='in', handle=h, psi=p))
                                  for h in ('klein', 'torus') for p in np.arange(8)*np.pi/4}
    if flux:
        rows = figure_flux(flux, IMAGES/'flux_cross.png')
        print(flux_table(rows), '\n')
        report['flux'] = [{k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in r.items()} for r in rows]
    if offs:
        rows = figure_offsets(offs, IMAGES/'clock_offsets.png')
        print(offset_table(rows), '\n')
        report['offsets'] = rows
    if checks and scan:
        report['checks'] = checks_summary(checks, scan)
        print(json.dumps(report['checks'], indent=2))
    if offs and ochecks:
        report['offset_checks'] = offset_checks_summary(offs, ochecks)
        print(json.dumps(report['offset_checks'], indent=2))
    (RUNS/'summary.json').write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    files = sorted(list(RUNS.glob('*.json'))+list(IMAGES.glob('*.png')))
    (RUNS/'sha256.json').write_text(json.dumps({p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files if p.name != 'sha256.json'},
                                               indent=2, sort_keys=True)+'\n')


if __name__ == '__main__':
    main()
