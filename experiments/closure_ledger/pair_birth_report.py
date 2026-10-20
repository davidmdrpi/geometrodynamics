#!/usr/bin/env python
"""Tables and figures for docs/pair_birth_timeline.md from the pair-birth probe runs.

    python experiments/closure_ledger/pair_birth_report.py
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from geometrodynamics.waves.pair_birth import PairBirth, Parameters, bare_timeline, front_crossings

RUNS = Path(__file__).resolve().parent/'runs'/'20261019_pair_birth'
IMAGES = Path(__file__).resolve().parents[2]/'docs'/'images'/'pair_birth'
A = .3
WINDOW = .35
MODES = ('full', 'dressed', 'instant')
COLORS = dict(full='#c0392b', dressed='#1f5fa8', instant='#7f8c8d', free='#138d75')
LABELS = dict(full='full field (birth flash)', dressed='no flash (dressed birth)', instant='instantaneous Coulomb only',
              free='no field (q = 0)')


def load(part):
    runs = json.loads((RUNS/f'{part}.json').read_text())['runs']
    for r in runs:
        for k in ('t', 'phi', 'd', 'alive', 'force_coulomb', 'force_exchange', 'force_own', 'force_flash', 'J', 'x'):
            if k in r:
                r[k] = np.array(r[k])
        r['energy'] = {k: np.array(v) for k, v in r['energy'].items()}
    return runs


def pick(runs, **kw):
    out = [r for r in runs if all((abs(r[k]-v) < 1e-9) if isinstance(v, float) else r[k] == v for k, v in kw.items())]
    if len(out) != 1:
        raise KeyError(f'{kw}: {len(out)} matches')
    return out[0]


def side_marker(side):
    return {'P': 'o', '-P': '^', None: 'x'}[side]


# ---- events -------------------------------------------------------------------------------------
def window_impulse(r, c):
    t, f = r['t'], r['force_flash'][:, 0 if c['source'] == 'own' else 1]
    sel = np.abs(t-c['time']) <= WINDOW
    return float(np.trapezoid(f[sel], t[sel]))


def first_return_events(r, crossings, lo):
    """Crossings of mouth A by the flash after one or two antipodal foci, with the flash impulse
    of the matching source in a window around each, and whether it points along the ray."""
    hi = (r['contact_time'] or r['t'][-1]) - WINDOW
    rows = []
    for c in crossings:
        if c['mouth'] == 0 and c['foci'] in (1, 2) and lo <= c['time'] <= hi:
            imp = window_impulse(r, c)
            rows.append(dict(time=c['time'], source=c['source'], direction=c['direction'], foci=c['foci'],
                             impulse=imp, along_ray=bool(np.sign(imp) == c['direction'])))
    return rows


def direct_events(r):
    return [dict(time=c['time'], source=c['source'], direction=c['direction'], impulse=window_impulse(r, c))
            for c in r['crossings'] if c['mouth'] == 0 and c['foci'] == 0 and c['source'] == 'partner']


def event_summary(featured, retimed):
    rows = []
    for r in featured:
        rows += [dict(run=f"natural v={r['v']:.3f}", prediction='flash timing', **e)
                 for e in first_return_events(r, r['crossings'], 1.5)]
    for r in retimed:
        td, shift = r['retime']
        if shift is None:
            continue
        al = r['alive']
        unshifted = [c for c in front_crossings(r['phi'][al], r['t'][al]) if c['time'] >= td]
        name = f"shift {shift:+.1f} v={r['v']:.3f}"
        rows += [dict(run=name, prediction='flash timing', **e)
                 for e in first_return_events(r, [c for c in r['crossings'] if c['time'] >= td], td)]
        rows += [dict(run=name, prediction='unshifted timing (wrong)', **e)
                 for e in first_return_events(r, unshifted, td)]
    stats = {}
    for pred in ('flash timing', 'unshifted timing (wrong)'):
        sel = [e for e in rows if e['prediction'] == pred]
        stats[pred] = dict(n=len(sel), along_ray=int(sum(e['along_ray'] for e in sel)),
                           mean_abs_impulse=float(np.mean([abs(e['impulse']) for e in sel])) if sel else None)
    return rows, stats


# ---- tables -------------------------------------------------------------------------------------
def outcome_table(scan):
    lines = ['| v | free (q = 0) | full field | no flash | instantaneous Coulomb |', '|---|---|---|---|---|']
    for v in sorted({r['v'] for r in scan}):
        cells = [f'{(np.pi-2*A)/v:.2f} at −P']
        for mode in MODES:
            r = pick(scan, v=v, mode=mode)
            cells.append('none by 34' if r['contact_time'] is None else
                         f"{r['contact_time']:.2f} at {'P' if r['contact_side'] == 'P' else '−P'}")
        lines.append(f'| {v:.3f} | ' + ' | '.join(cells) + ' |')
    return '\n'.join(lines)


def threshold(scan, mode):
    vs = sorted({r['v'] for r in scan})
    back = [v for v in vs if pick(scan, v=v, mode=mode)['contact_side'] == 'P']
    over = [v for v in vs if pick(scan, v=v, mode=mode)['contact_side'] == '-P']
    return max(back), min(over)


def contact_ledger(r):
    i = np.flatnonzero((np.abs(r['t']-r['contact_time']) < 1e-5) & r['alive'])[-1]     # just before surgery
    e = {k: float(v[i]) for k, v in r['energy'].items()}
    e0 = {k: float(v[0]) for k, v in r['energy'].items()}
    return dict(kinetic0=e0['kinetic'], coulomb0=e0['coulomb'], kinetic=e['kinetic'], coulomb=e['coulomb'],
                transverse=e['transverse'], released=r['released_energy'],
                speed=float(np.linalg.norm(r['contact_speed'][0])))


def flash_energy(scan, v):
    """Energy of the birth flash: a free wave, so it equals the dressing it removes."""
    return float(pick(scan, v=v, mode='dressed')['energy']['transverse'][0])


def ledger_table(scan):
    lines = ['| v | kinetic at birth | birth flash | kinetic at contact | transverse field at contact | released |',
             '|---|---|---|---|---|---|']
    for v in sorted({r['v'] for r in scan}):
        r = pick(scan, v=v, mode='full')
        if r['contact_time'] is None:
            continue
        c = contact_ledger(r)
        lines.append(f"| {v:.3f} | {c['kinetic0']:.4f} | {flash_energy(scan, v):.4f} | {c['kinetic']:.4f} | "
                     f"{c['transverse']:.4f} | {c['released']:.4f} |")
    return '\n'.join(lines)


def retimed_table(retimed, featured):
    lines = ['| v | flash after t = 3 | returns shift | energy change at t = 3 | contact |', '|---|---|---|---|---|']
    for v in (.30, 1/3):
        r = pick(featured, v=v)
        lines.append(f"| {v:.3f} | natural | 0 | 0 | {r['contact_time']:.2f} at P |")
        for shift in (None, -1., -.5, .5, 1.):
            r = [q for q in retimed if abs(q['v']-v) < 1e-9 and q['retime'][1] == shift][0]
            what = 'absorbed' if shift is None else 're-timed'
            sh = '—' if shift is None else f'{-shift:+.1f} (later)' if shift < 0 else f'{-shift:+.1f} (earlier)'
            lines.append(f"| {v:.3f} | {what} | {sh} | {r['retime_energy']:+.4f} | "
                         f"{r['contact_time']:.2f} at {'P' if r['contact_side'] == 'P' else '−P'} |")
    return '\n'.join(lines)


# ---- figures ------------------------------------------------------------------------------------
def spacetime_bare(ax, operator, v=1/3, end=4*np.pi):
    th, t, e, _, _ = bare_timeline(operator=operator, t=np.linspace(0, end, 641))
    full = np.concatenate([e[:, ::-1], e[:, 1:]], 1)
    phi = np.concatenate([-th[::-1], th[1:]])
    ax.imshow(np.log10(full+1e-3), origin='lower', aspect='auto', cmap='magma', vmin=-2, vmax=2,
              extent=[phi[0], phi[-1], t[0], t[-1]])
    tt = np.linspace(0, (np.pi-2*A)/v, 200)
    for s in (1, -1):
        ax.plot(s*(A+v*tt), tt, color='w', ls='--', lw=1.2)
    for k in range(1, 4):
        ax.axhline(k*np.pi, color='w', lw=.4, alpha=.5)
    ax.text(0, 2*np.pi+.25, 'back at P, sign flipped' if operator == 'conformal' else 'back near P', color='w',
            ha='center', fontsize=8)
    ax.text(-np.pi+.1, np.pi+.25, 'focus at −P', color='w', fontsize=8)
    ax.set_xlabel('angle along the birth circle from P')
    ax.set_ylabel('time t')
    ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi], ['−P', '', 'P', '', '−P'])


def spacetime_run(ax, r, end=34.):
    """Ring-averaged energy density around P, mirrored to the two sides of the birth circle."""
    mp = r['map']
    t = np.array(mp['t'])
    mirror = lambda a: np.concatenate([np.array(a)[:, ::-1], np.array(a)[:, 1:]], 1)
    tr, de = mirror(mp['log_transverse']), mirror(mp['log_deformation'])
    tc = r['contact_time']
    before, after = t < tc, t >= tc
    ax.imshow(np.where(before[:, None], tr, np.nan), origin='lower', aspect='auto', cmap='magma', vmin=-4.5, vmax=-1.5,
              extent=[-np.pi, np.pi, t[0], t[-1]], interpolation='nearest')
    ax.imshow(np.where(after[:, None], de, np.nan), origin='lower', aspect='auto', cmap='viridis', vmin=-2, vmax=2,
              extent=[-np.pi, np.pi, t[0], t[-1]], interpolation='nearest')
    al = r['alive']
    for j in range(2):
        ax.plot(r['phi'][al, j], r['t'][al], color='w', lw=1.6)
    tt = np.linspace(0, r['free_contact'], 200)
    for s in (1, -1):
        ax.plot(s*(A+r['v']*tt), tt, color='w', ls=':', lw=1)
    for c in r['crossings']:
        if c['foci'] in (1, 2):
            ax.plot(r['phi'][al, c['mouth']][np.argmin(np.abs(r['t'][al]-c['time']))], c['time'], 'o', ms=4,
                    mfc='none', mec='#7fdbff')
    ax.plot(0, tc, '*', color='#ffd166', ms=13, mec='k')
    for k, lab in ((1, 'released pulse focuses at −P'), (2, 'back at P, sign flipped')):
        if tc+k*np.pi < end:
            ax.axhline(tc+k*np.pi, color='w', lw=.4, alpha=.6)
            ax.text(-np.pi+.1, tc+k*np.pi+.25, lab, color='w', fontsize=7)
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(0, end)
    ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi], ['−P', '', 'P', '', '−P'])
    ax.set_xlabel('distance from P (mirrored: the two sides of the birth circle)')
    ax.set_ylabel('time t')


def figure_timeline(featured, path):
    fig, axes = plt.subplots(1, 3, figsize=(15, 7.4), gridspec_kw=dict(wspace=.28))
    spacetime_bare(axes[0], 'conformal')
    axes[0].set_title('(a) No mouths: a pulse released at P\n'
                      '(conformal operator; dashed: free mouths, v = 1/3)', fontsize=10)
    for ax, v, lab in ((axes[1], .30, 'b'), (axes[2], 1/3, 'c')):
        r = pick(featured, v=v)
        spacetime_run(ax, r)
        ax.set_title(f'({lab}) Mouths born touching at P, v = {v:.3f}\n'
                     f'(contact at t = {r["contact_time"]:.2f}, then the released pulse)', fontsize=10)
    fig.text(.5, .0, 'Colour: ring-averaged energy density around P, log scale. Magma: the transverse (flux) field '
             'before contact; green: the released deformation pulse after it.\nWhite: mouth worldlines (solid) and '
             'free motion (dotted). Rings: first returns of the birth flash after one or two antipodal foci. '
             'Star: rim contact and release.', ha='center', va='top', fontsize=9)
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def figure_outcomes(scan, checks, path):
    m = PairBirth(Parameters())
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), gridspec_kw=dict(wspace=.3))
    ax = axes[0]
    d = np.linspace(2*A, np.pi, 300)
    ax.plot(d, m.coulomb_interaction(d)-m.coulomb_interaction(2*A), color=COLORS['full'],
            label='flux mouths (opposite Wheeler charges)')
    ax.plot(d, m.scalar_interaction(d)-m.scalar_interaction(2*A), color='#8e44ad', ls='--',
            label='same charges as scalar sources')
    ax.axhline(0, color='k', lw=.5)
    ax.set_xlabel('separation d')
    ax.set_ylabel('static interaction energy, from contact')
    ax.set_title('(a) Sign of the static force\nflux mouths attract; the pull vanishes at d = π', fontsize=10)
    ax.legend(fontsize=8)
    ax = axes[1]
    for v, ls in ((.25, '-'), (.30, '-'), (1/3, '-'), (.36, '-')):
        for mode in ('full', 'instant'):
            r = pick(scan, v=v, mode=mode)
            al = r['alive']
            ax.plot(r['t'][al], r['d'][al], color=COLORS[mode], lw=1.4 if mode == 'full' else .9,
                    ls=ls if mode == 'full' else '--')
            if r['contact_time'] is not None:
                ax.plot(r['contact_time'], 2*A, side_marker(r['contact_side']), color=COLORS[mode], ms=6)
        ax.text(pick(scan, v=v, mode='full')['t_d_max'], pick(scan, v=v, mode='full')['d_max']+.04, f'{v:.2f}',
                fontsize=7, color=COLORS['full'], ha='center')
    ax.axhline(np.pi, color='k', lw=.5)
    ax.set_xlabel('time t')
    ax.set_ylabel('separation d')
    ax.set_ylim(0, 3.3)
    ax.set_title('(b) Separation: full field (solid), instantaneous Coulomb (dashed)\n'
                 'o: contact back at P,  ▲: contact at −P', fontsize=10)
    ax = axes[2]
    vs = np.array(sorted({r['v'] for r in scan}))
    ax.plot(vs, (np.pi-2*A)/vs, color=COLORS['free'], marker='^', mfc='none', ls=':', label=LABELS['free'])
    for mode in MODES:
        for v in vs:
            r = pick(scan, v=v, mode=mode)
            if r['contact_time'] is not None:
                ax.plot(v, r['contact_time'], side_marker(r['contact_side']), color=COLORS[mode],
                        mfc=COLORS[mode] if r['contact_side'] == 'P' else 'none', ms=7)
        tc = [pick(scan, v=v, mode=mode)['contact_time'] for v in vs]
        ax.plot(vs, [np.nan if c is None else c for c in tc], color=COLORS[mode], lw=.8, label=LABELS[mode])
    ax.set_xlabel('birth speed v of each mouth')
    ax.set_ylabel('contact time')
    ax.set_ylim(0, 35)
    ax.set_title('(c) Annihilation time and place\nfilled o: back at P; open ▲: at −P', fontsize=10)
    ax.legend(fontsize=7)
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def force_panel(ax, r, title):
    al = r['alive'] & (r['t'] <= r['contact_time'] - .05)
    t = r['t'][al]
    flash = r['force_flash'][al]
    transverse = (r['force_exchange'] + r['force_own'])[al]
    ax.plot(t, r['force_coulomb'][al], color='k', lw=1, label='Coulomb (static part, from the partner)')
    ax.plot(t, flash[:, 0], color='#2e86de', lw=1.1, label="birth flash: mouth's own")
    ax.plot(t, flash[:, 1], color='#e67e22', lw=1.1, label="birth flash: partner's")
    ax.plot(t, transverse - flash.sum(1), color='#95a5a6', lw=.8, label='rest of the transverse field')
    for c in r['crossings']:
        if c['mouth'] == 0 and c['foci'] >= 1:
            ax.axvline(c['time'], color='#2e86de' if c['source'] == 'own' else '#e67e22',
                       ls='-' if c['direction'] < 0 else ':', lw=.8, alpha=.8)
    ax.axhline(0, color='k', lw=.4)
    ax.set_ylim(-.08, .06)
    ax.set_ylabel('outward force on A')
    ax.set_xlabel('time t')
    ax.set_title(title, fontsize=10)
    tw = ax.twinx()
    tw.plot(r['t'][al], r['d'][al], color=COLORS['free'], ls='--', lw=.9)
    tw.set_ylim(0, np.pi)
    tw.set_ylabel('separation d (green dashed)', fontsize=8)


def figure_flash(featured, retimed, scan, path):
    fig = plt.figure(figsize=(15, 9.2))
    gs = fig.add_gridspec(2, 2, hspace=.42, wspace=.25)
    force_panel(fig.add_subplot(gs[0, 0]), pick(featured, v=.30),
                '(a) v = 0.30: forces on A; lines: flash rays crossing A\n(solid: ray heading to P, dotted: heading away)')
    force_panel(fig.add_subplot(gs[0, 1]), pick(featured, v=1/3),
                '(b) v = 1/3: four returns while the pair lingers near antipodes')
    fig.axes[0].legend(fontsize=7, loc='lower right')
    for k, v in enumerate((.30, 1/3)):
        ax = fig.add_subplot(gs[1, k])
        r = pick(featured, v=v)
        curves = [('natural flash', r, 'k', '-')]
        for shift, col in ((None, '#7f8c8d'), (-1., '#8e44ad'), (-.5, '#c39bd3'), (.5, '#f5b041'), (1., '#d35400')):
            q = [x for x in retimed if abs(x['v']-v) < 1e-9 and x['retime'][1] == shift][0]
            lab = 'flash absorbed at t = 3' if shift is None else f'returns {abs(shift):.1f} {"later" if shift < 0 else "earlier"}'
            curves.append((lab, q, col, '-' if shift is None else '--'))
        curves.append(('no flash (dressed birth)', pick(scan, v=v, mode='dressed'), COLORS['dressed'], ':'))
        for lab, q, col, ls in curves:
            al = q['alive']
            ax.plot(q['t'][al], q['d'][al], color=col, ls=ls, lw=1.2, label=lab)
            if q['contact_time'] is not None:
                ax.plot(q['contact_time'], 2*A, side_marker(q['contact_side']), color=col, ms=7,
                        mfc=col if q['contact_side'] == 'P' else 'none')
        ax.axvline(3., color='k', lw=.5, ls=':')
        ax.axhline(np.pi, color='k', lw=.4)
        ax.set_ylim(0, 3.3)
        ax.set_xlim(0, 35)
        ax.set_xlabel('time t')
        ax.set_ylabel('separation d')
        ax.set_title(f'({"cd"[k]}) v = {v:.3f}: same birth, flash re-timed at t = 3\n'
                     'filled o: back at P; open ▲: at −P', fontsize=10)
        ax.legend(fontsize=7, loc='lower center' if v < .32 else 'lower right')
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def figure_ledger(featured, checks, path):
    r = pick(featured, v=1/3)
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), gridspec_kw=dict(wspace=.3))
    t, E = r['t'], r['energy']
    ax = axes[0]
    al = r['alive']
    ax.plot(t[al], E['kinetic'][al], color=COLORS['full'], label='mouth kinetic')
    ax.plot(t[al], E['coulomb'][al]-E['coulomb'][0], color='k', label='Coulomb (change from birth)')
    ax.plot(t, E['transverse'], color='#e67e22', label='transverse field (waves)')
    total = E['kinetic']+E['coulomb']-E['coulomb'][0]+E['transverse']
    ax.plot(t[al], total[al], color='#7f8c8d', ls='--', label='sum')
    ax.axvline(r['contact_time'], color='k', lw=.5, ls=':')
    ax.set_xlabel('time t')
    ax.set_ylabel('energy')
    ax.set_title('(a) v = 1/3: exchange between mouths and field', fontsize=10)
    ax.legend(fontsize=7)
    ax = axes[1]
    parts = [('rest (neck and flux, 2M)', E['rest'], '#5d6d7e'), ('kinetic', E['kinetic'], COLORS['full']),
             ('Coulomb', E['coulomb'], 'k'), ('transverse', E['transverse'], '#e67e22'),
             ('deformation (released)', E['deformation'], '#27ae60')]
    ax.stackplot(t, *[p[1] for p in parts], labels=[p[0] for p in parts], colors=[p[2] for p in parts], alpha=.85)
    ax.plot(t, r['energy']['rest']*0 + sum(p[1] for p in parts), color='w', lw=.8)
    ax.set_xlabel('time t')
    ax.set_ylabel('energy')
    ax.set_title('(b) Full ledger through contact and release\nclock support: none (no clock offset)', fontsize=10)
    ax.legend(fontsize=7, loc='lower left')
    ax = axes[2]
    q = pick(checks, tag='independent')
    x = np.array(q['x'])
    lon = np.arctan2(x[..., 0], x[..., 2])
    lat = np.arcsin(np.clip(x[..., 1], -1, 1))
    al = q['alive']
    for j, col in ((0, COLORS['full']), (1, COLORS['dressed'])):
        ax.plot(lon[al, j], lat[al, j], color=col, label=f'mouth {"AB"[j]}')
        ax.plot(lon[0, j], lat[0, j], 'o', color=col)
    k = np.flatnonzero(al)[-1]
    ax.plot(lon[k].mean(), lat[k].mean(), '*', color='#ffd166', ms=13, mec='k', label='contact')
    ax.set_xlabel('angle from P along the birth circle')
    ax.set_ylabel('latitude off the circle')
    ax.set_title(f'(c) Independent mouths (unequal speeds, off-circle kick)\n'
                 f'energy drift {q["energy_drift"]:.1e}, angular momentum drift {q["angular_drift"]:.1e}', fontsize=10)
    ax.legend(fontsize=7)
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def load_history():
    path = RUNS/'history.json'
    if not path.exists():
        return []
    runs = json.loads(path.read_text())['runs']
    for r in runs:
        for k in ('t', 'd', 'tau', 'alive', 'cycle', 'J', 'x', 'clock', 'external', 'first_order', 'chronology_margin'):
            r[k] = np.array([np.nan if v is None else v for v in r[k]] if k in
                            ('clock', 'external', 'first_order', 'chronology_margin') else r[k], dtype=float
                            if k != 'alive' else bool)
        r['energy'] = {k: np.array(v) for k, v in r['energy'].items()}
    return runs


def history_name(r):
    return dict(history_symmetric=f'symmetric, v = {r["v"]:.2f}', history_independent='unequal + off-circle kick',
                history_lopsided='lopsided (0.45 / 0.15)')[r['tag']]


def history_rows(runs):
    rows = []
    for r in runs:
        births = [dict(time=0., speed=None)] + [e for e in r['events'] if e['kind'] == 'birth']
        contacts = [e for e in r['events'] if e['kind'] == 'contact']
        for c, (b, k) in enumerate(zip(births, contacts)):
            sel = r['alive'] & (r['cycle'] == c)
            ext = np.abs(r['external'][sel])
            rows.append(dict(run=history_name(r), cycle=c, birth=b['time'], contact=k['time'],
                             flight=k['time']-b['time'], contact_point=k['point'], contact_speed=k['speed'],
                             max_offset=float(np.nanmax(ext)) if np.isfinite(ext).any() else 0.,
                             offset_at_contact=float(k['clocks'][0]-k['clocks'][1]),
                             min_margin=float(np.nanmin(r['chronology_margin'][sel])),
                             released=k['energy'], rebirth_correction=None if b['speed'] is None else b['correction']))
    return rows


def history_table(rows):
    lines = ['| run | cycle | born | contact | place | max throat offset | offset at contact | min chronology margin |',
             '|---|---|---|---|---|---|---|---|']
    for q in rows:
        x, y, z = q['contact_point']
        place = ('P' if z > .999 else '−P' if z < -.999 else
                 f'circle, {np.degrees(np.arctan2(x, z)):.0f}° from P' if abs(y) < 1e-6 else 'off the circle')
        lines.append(f"| {q['run']} | {q['cycle']+1} | {q['birth']:.2f} | {q['contact']:.2f} | {place} | "
                     f"{q['max_offset']:.4f} | {q['offset_at_contact']:+.5f} | {q['min_margin']:.3f} |")
    return '\n'.join(lines)


def figure_history(runs, path):
    fig, axes = plt.subplots(2, 2, figsize=(15, 9), gridspec_kw=dict(hspace=.38, wspace=.22))
    ax = axes[0, 0]
    for r, col in ((pick(runs, tag='history_symmetric', v=.30), COLORS['full']),
                   (pick(runs, tag='history_symmetric', v=1/3), COLORS['dressed'])):
        al = r['alive']
        ax.plot(np.where(al, r['t'], np.nan), np.where(al, r['d'], np.nan), color=col, lw=1.3, label=history_name(r))
        for e in r['events']:
            ax.plot(e['time'], 2*A, '*' if e['kind'] == 'contact' else '^', color=col, ms=9 if e['kind'] == 'contact'
                    else 6, mec='k', mew=.4)
            if e['kind'] == 'contact':
                ax.text(e['time'], 2*A-.22, 'P' if e['point'][2] > 0 else '−P', ha='center', fontsize=7, color=col)
    ax.axhline(np.pi, color='k', lw=.4)
    ax.set_ylim(0, 3.3)
    ax.set_xlabel('time t')
    ax.set_ylabel('separation d')
    ax.set_title('(a) Repeated cycles: born touching (▲), contact (★), re-born at the antipode after π',
                 fontsize=10)
    ax.legend(fontsize=8, loc='upper right')
    ax = axes[0, 1]
    for r, col in ((pick(runs, tag='history_lopsided'), '#8e44ad'), (pick(runs, tag='history_independent'), '#e67e22'),
                   (pick(runs, tag='history_symmetric', v=.30), COLORS['full'])):
        ax.plot(r['t'], r['external'], color=col, lw=1.2, label=history_name(r))
        for e in r['events']:
            if e['kind'] == 'birth':
                ax.axvline(e['time'], color=col, lw=.5, ls=':')
    ax.axhline(0, color='k', lw=.4)
    ax.set_xlabel('time t')
    ax.set_ylabel('throat clock offset (external time)')
    ax.set_title('(b) History offset: zero at every birth, set by the velocity history, back near zero at contact',
                 fontsize=10)
    ax.legend(fontsize=8)
    ax = axes[1, 0]
    for r, col in ((pick(runs, tag='history_lopsided'), '#8e44ad'), (pick(runs, tag='history_independent'), '#e67e22')):
        al = r['alive']
        ax.plot(np.where(al, r['t'], np.nan), np.where(al, r['d'] + r['neck_length'], np.nan), color=col, lw=1,
                ls='--', label=f'{history_name(r)}: d + L_neck')
        ax.plot(r['t'], np.abs(r['external']), color=col, lw=1.4, label=f'{history_name(r)}: |offset|')
    ax.set_xlabel('time t')
    ax.set_ylabel('external time')
    ax.set_title('(c) Chronology check: a closed loop needs |offset| ≥ d + L_neck; the gap never closes', fontsize=10)
    ax.legend(fontsize=7, loc='upper right')
    ax = axes[1, 1]
    r = pick(runs, tag='history_symmetric', v=.30)
    E = r['energy']
    ax.plot(r['t'], E['kinetic'], color=COLORS['full'], label='mouth kinetic')
    ax.plot(r['t'], E['transverse'], color='#e67e22', label='transverse field (waves)')
    ax.plot(r['t'], np.where(r['alive'], E['coulomb']-E['coulomb'][0], np.nan), color='k', lw=.9,
            label='Coulomb (change from birth)')
    births = [e for e in r['events'] if e['kind'] == 'birth']
    for e in births:
        ax.axvline(e['time'], color='k', lw=.4, ls=':')
    ax.set_xlabel('time t')
    ax.set_ylabel('energy')
    ax.set_title('(d) Symmetric v = 0.30: each flight leaves energy in the waves, so each birth is slower',
                 fontsize=10)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def main():
    scan, featured, retimed, checks = (load(p) for p in ('scan', 'featured', 'retimed', 'checks'))
    history = load_history()
    IMAGES.mkdir(parents=True, exist_ok=True)
    if history:
        figure_history(history, IMAGES/'history.png')
    figure_timeline(featured, IMAGES/'timeline.png')
    figure_outcomes(scan, checks, IMAGES/'outcomes.png')
    figure_flash(featured, retimed, scan, IMAGES/'flash_returns.png')
    figure_ledger(featured, checks, IMAGES/'ledger.png')
    events, stats = event_summary(featured, retimed)
    m = PairBirth(Parameters())
    report = dict(
        thresholds={mode: threshold(scan, mode) for mode in MODES},
        coulomb_barrier=float(m.coulomb_interaction(np.pi)-m.coulomb_interaction(2*A)),
        flash_energy={f'{v:.3f}': flash_energy(scan, v) for v in sorted({r['v'] for r in scan})},
        direct_exchange=[dict(v=r['v'], **e) for r in featured for e in direct_events(r)],
        first_returns=stats, first_return_events=events,
        featured={f"{r['v']:.3f}": dict(contact_time=r['contact_time'], contact_side=r['contact_side'],
                                        d_max=r['d_max'], t_d_max=r['t_d_max'], **contact_ledger(r),
                                        energy_drift=r['energy_drift'], angular_drift=r['angular_drift'])
                  for r in featured},
        retimed=[dict(v=r['v'], shift=r['retime'][1], retime_energy=r['retime_energy'],
                      contact_time=r['contact_time'], contact_side=r['contact_side']) for r in retimed],
        checks={f"{r['tag']} v={r['v']:.3f}": dict(contact_time=r['contact_time'], contact_side=r['contact_side'],
                                                    energy_drift=r['energy_drift'], angular_drift=r['angular_drift'],
                                                    released=r['released_energy'],
                                                    unreleased_angular=r['unreleased_angular'])
                for r in checks},
        max_energy_drift=max(r['energy_drift'] for part in (scan, featured, retimed, checks) for r in part),
        max_mirror_asymmetry=max(max(r['asymmetry_coulomb'], r['asymmetry_exchange'], r['asymmetry_own'])
                                 for part in (scan, featured, retimed) for r in part),
        history=dict(rows=history_rows(history),
                     runs={f"{r['tag']} v={r['v']:.3f}": dict(energy_drift=r['energy_drift'],
                                                             angular_drift=r['angular_drift'],
                                                             max_external_offset=r['max_external_offset'],
                                                             max_clock_difference=r['max_clock_difference'],
                                                             min_chronology_margin=r['min_chronology_margin'],
                                                             births=sum(e['kind'] == 'birth' for e in r['events'])+1)
                           for r in history}) if history else None,
        tables=dict(outcomes=outcome_table(scan), ledger=ledger_table(scan), retimed=retimed_table(retimed, featured),
                    history=history_table(history_rows(history)) if history else None))
    (RUNS/'summary.json').write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    files = sorted(RUNS.glob('*.json')) + sorted(IMAGES.glob('*.png'))
    (RUNS/'sha256.json').write_text(json.dumps({p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                                for p in files if p.name != 'sha256.json'}, indent=2)+'\n')
    for k in ('outcomes', 'ledger', 'retimed', 'history'):
        print(report['tables'][k], '\n')
    print(json.dumps(dict(thresholds=report['thresholds'], barrier=report['coulomb_barrier'],
                          first_returns=stats, direct=report['direct_exchange'],
                          max_energy_drift=report['max_energy_drift'],
                          max_mirror_asymmetry=report['max_mirror_asymmetry']), indent=1))


if __name__ == '__main__':
    main()
