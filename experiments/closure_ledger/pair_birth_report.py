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
    return load_history_runs(json.loads(path.read_text())['runs'])


def load_history_runs(runs):
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


NECK_ARRAYS = ('masses', 'mu', 'spin', 'ends', 'throat_force', 'port_force', 'exchange_force')
NECK_COLORS = dict(none='#7f8c8d', off='#1f5fa8', history='#c0392b', reversed='#e67e22')
NECK_LABELS = dict(none='no neck', off='neck, clock off', history='neck, history clock',
                   reversed='neck, reversed clock')


def load_neck(part):
    path = RUNS/f'{part}.json'
    if not path.exists():
        return []
    runs = load_history_runs(json.loads(path.read_text())['runs'])
    for r in runs:
        for k in NECK_ARRAYS:
            if r.get(k) is not None:
                r[k] = np.array(r[k], float)
        r['birth'], r['clock_mode'] = neck_tag(r['tag'])
    return runs


def neck_tag(tag):
    """('lopsided', 'history') from 'neck_lopsided_history'; the torus runs read 'lopsided_torus'."""
    *birth, clock = tag.split('_')[1:]
    return '_'.join(birth), clock


def contacts(r):
    return [e for e in r['events'] if e['kind'] == 'contact']


def neck_rows(runs):
    """Per birth and clock mode: contacts, offsets, clock momentum, neck share and conservation."""
    rows = []
    for r in runs:
        al = r['alive']
        e = r['energy']
        neck = e.get('neck', np.zeros_like(r['t'])) + e.get('port', np.zeros_like(r['t']))
        rows.append(dict(
            birth=r['birth'], clock=r['clock_mode'],
            contacts=[round(c['time'], 5) for c in contacts(r)],
            places=[c['point'] for c in contacts(r)],
            offsets_at_contact=[c['clocks'][0] - c['clocks'][1] for c in contacts(r)],
            max_offset=float(np.nanmax(np.abs(r['clock'][al]))),
            max_mu=float(np.abs(r['mu']).max()) if r.get('mu') is not None else 0.,
            max_neck=float(neck.max()), mean_neck_share=float(np.mean((neck/e['transverse'].clip(1e-12))[al])),
            max_throat_force=float(r['throat_force'].max()) if r.get('throat_force') is not None else 0.,
            max_exchange_force=float(r['exchange_force'].max()) if r.get('exchange_force') is not None else None,
            energy_drift=r['energy_drift'], angular_drift=r['angular_drift'],
            min_margin=r['min_chronology_margin']))
    return rows


def feedback(runs, birth):
    """Contact-time shifts per cycle: clock feedback (history - off, reversed - off) and the
    neck's own effect (off - none)."""
    by = {r['clock_mode']: r for r in runs if r['birth'] == birth}
    t = {k: np.array([c['time'] for c in contacts(v)]) for k, v in by.items()}
    n = min(len(v) for v in t.values())
    return dict(history=(t['history'][:n] - t['off'][:n]).tolist(),
                reversed=(t['reversed'][:n] - t['off'][:n]).tolist() if 'reversed' in t else None,
                neck=(t['off'][:n] - t['none'][:n]).tolist() if 'none' in t else None)


def neck_table(rows):
    lines = ['| birth | clock | contacts | max \\|offset\\| | offset at contact | max \\|μ\\| | neck + port, max | '
             'neck / waves | energy drift |', '|---|---|---|---|---|---|---|---|---|']
    for q in rows:
        lines.append(f"| {q['birth']} | {q['clock']} | {', '.join(f'{c:.3f}' for c in q['contacts'])} | "
                     f"{q['max_offset']:.4f} | {max(abs(o) for o in q['offsets_at_contact']):.1e} | "
                     f"{q['max_mu']:.2e} | {q['max_neck']:.2e} | {q['mean_neck_share']:.3f} | {q['energy_drift']:.1e} |")
    return '\n'.join(lines)


def calibration_rows(checks):
    rows = {}
    for r in checks:
        key = r['tag'].rsplit('_', 1)[0]
        rows.setdefault(key, {})[r['tag'].rsplit('_', 1)[1]] = r
    out = []
    for key, pair in sorted(rows.items()):
        h, o = pair['history'], pair['off']
        ch, co = contacts(h)[0], contacts(o)[0]
        e = h['energy']
        out.append(dict(run=key, port=h['neck']['port'], nodes=h['neck'].get('nodes', 65), dt=h['dt'],
                        contact_shift=ch['time'] - co['time'], mu_at_contact=float(h['mu'][h['alive']][-1]),
                        mean_neck=float(np.mean((e['neck'] + e['port'])[h['alive']])),
                        offset_shift=(ch['clocks'][0] - ch['clocks'][1]) - (co['clocks'][0] - co['clocks'][1]),
                        energy_drift=max(h['energy_drift'], o['energy_drift'])))
    return out


def calibration_table(rows):
    lines = ['| port g | nodes | dt | mean neck + port energy | μ at contact | contact shift, history − off | '
             'energy drift |', '|---|---|---|---|---|---|---|']
    for q in sorted(rows, key=lambda q: (q['port'], q['nodes'], -q['dt'])):
        lines.append(f"| {q['port']:g} | {q['nodes']} | {q['dt']:g} | {q['mean_neck']:.2e} | {q['mu_at_contact']:+.2e} | "
                     f"{q['contact_shift']:+.2e} | {q['energy_drift']:.1e} |")
    return '\n'.join(lines)


def figure_neck(runs, checks, path):
    fig, axes = plt.subplots(2, 2, figsize=(16, 9.5), gridspec_kw=dict(hspace=.42, wspace=.42))
    lop = {r['clock_mode']: r for r in runs if r['birth'] == 'lopsided'}
    first = lambda r: r['alive'] & (r['cycle'] == 0)
    ax = axes[0, 0]
    r = lop['history']
    ax.plot(r['t'], r['clock'], color='k', lw=1.2, label='offset A = τ_A − τ_B (history)')
    ax.set_xlabel('time t')
    ax.set_ylabel('offset A')
    tw = ax.twinx()
    for k in ('history', 'reversed'):
        tw.plot(lop[k]['t'], lop[k]['mu'], color=NECK_COLORS[k], lw=1.2, label=f'μ, {NECK_LABELS[k]}')
    tw.set_ylabel('clock momentum μ')
    for e in contacts(r):
        ax.axvline(e['time'], color='k', lw=.4, ls=':')
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = tw.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7, loc='lower left')
    ax.set_title('(a) Lopsided birth: the offset loads the clock momentum μ\n'
                 '(rest energies M ± μ); reversing the clock flips its sign', fontsize=10)
    ax = axes[0, 1]
    for k, ref, lab in (('history', 'off', 'history − off (clock feedback)'),
                        ('reversed', 'off', 'reversed − off'), ('off', 'none', 'off − no neck (neck itself)')):
        a, b = lop[k], lop[ref]
        n = min(np.flatnonzero(first(a))[-1], np.flatnonzero(first(b))[-1])
        same = np.abs(a['t'][:n] - b['t'][:n]) < 1e-9                     # common output times before contact
        ax.plot(a['t'][:n][same], (a['d'][:n] - b['d'][:n])[same], color=NECK_COLORS[k if k != 'off' else 'none'],
                lw=1.2, label=lab)
    ax.axhline(0, color='k', lw=.4)
    ax.set_xlabel('time t (first flight)')
    ax.set_ylabel('separation difference')
    ax.set_title('(b) Lopsided birth, first flight: what the clock and the neck\nchange in the separation',
                 fontsize=10)
    ax.legend(fontsize=7)
    ax = axes[1, 0]
    r = lop['history']
    e = r['energy']
    ax.plot(r['t'], e['transverse'], color='#7f8c8d', lw=1, label='transverse field (sphere)')
    ax.plot(r['t'], e['neck'], color=NECK_COLORS['history'], lw=1.2, label='neck field')
    ax.plot(r['t'], e['port'], color=NECK_COLORS['reversed'], lw=1, label='ports')
    ax.set_yscale('log')
    ax.set_ylim(1e-6, 1e-1)
    ax.set_xlabel('time t')
    ax.set_ylabel('energy')
    ax.legend(fontsize=7, loc='lower right')
    ax.set_title('(c) Lopsided birth, history clock: the neck holds a few per cent\nof the wave energy; '
                 'it is emptied at each contact', fontsize=10)
    ax = axes[1, 1]
    rows = calibration_rows(checks)
    for nodes, mk in ((65, 'o'), (129, 's')):
        sel = sorted((q for q in rows if q['nodes'] == nodes and q['dt'] > .003), key=lambda q: q['port'])
        ax.plot([q['port'] for q in sel], [q['contact_shift'] for q in sel], mk + '-', color='#c0392b',
                label=f'contact shift, {nodes} nodes')
        ax.plot([q['port'] for q in sel], [q['mu_at_contact'] for q in sel], mk + '--', color='#1f5fa8',
                label=f'μ at contact, {nodes} nodes')
    ax.axhline(0, color='k', lw=.4)
    ax.set_xscale('log')
    ax.set_xticks([25, 100, 400], ['25', '100', '400'])
    ax.minorticks_off()
    ax.set_xlabel('port stiffness g')
    ax.set_ylabel('time shift, or μ')
    ax.set_title('(d) Calibration: first lopsided contact, history clock against\nclock off, '
                 'at two neck resolutions', fontsize=10)
    ax.legend(fontsize=7)
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


ANTIPODAL_NAMES = dict(symmetric='symmetric, v = 0.30', symmetric_third='symmetric, v = 1/3',
                      independent='unequal + off-circle kick', lopsided='lopsided (0.45 / 0.15)',
                      fast_lopsided='fast lopsided (0.85 / 0.15)')


def load_antipodal(part='antipodal'):
    path = RUNS/f'{part}.json'
    if not path.exists():
        return []
    runs = json.loads(path.read_text())['runs']
    for r in runs:
        for k in ('t', 'd', 'alive', 'cycle', 'clock', 'total'):
            r[k] = np.array([np.nan if v is None else v for v in r[k]], dtype=float if k != 'alive' else bool)
        r['x'] = np.array(r['x'])
        r['energy'] = {k: np.array(v, float) for k, v in r['energy'].items()}
        if part == 'antipodal':
            r['birth_name'] = r['tag'][len('antipodal_'):].rsplit('_', 2 if r['tag'].endswith('torus') else 1)[0]
        else:                       # klein_<birth>_<handle> or klein_sweep_<birth>_<handle>_x<scale>
            name = r['tag'][len('klein_'):]
            r['sweep'] = name.startswith('sweep_')
            r['birth_name'] = name[len('sweep_'):].split(f"_{r['handle']}")[0] if r['sweep'] else \
                name[:-len(r['handle'])-1]
    return runs


def place(point):
    p = np.asarray(point)
    if p[2] > .999:
        return 'P'
    if p[2] < -.999:
        return '−P'
    if abs(p[1]) < 1e-6:
        return f'{np.degrees(np.arctan2(p[0], p[2])):.0f}°'
    return 'off circle'


def antipodal_rows(runs):
    rows = []
    for r in runs:
        ev = r['events']
        rows.append(dict(
            tag=r['tag'], birth=r['birth_name'], static=r['static'], handle=r['handle'],
            speeds=r['birth_speed'], contacts=[(e['time'], place(e['point'])) for e in ev if e['kind'] == 'contact'],
            captures=sum(e['kind'] == 'capture' for e in ev), exits=sum(e['kind'] == 'throat_exit' for e in ev),
            advanced=sum(e['kind'] == 'throat_exit' and e['advanced'] for e in ev),
            causality=sum(e['kind'] == 'capture' and e.get('exit') == 'causality_boundary' for e in ev),
            reflections=sum(e['kind'] == 'reflection' for e in ev),
            captured_energy=sum(e['into_throat'] for e in ev if e['kind'] == 'capture'),
            max_offset=float(np.nanmax(np.abs(r['clock']))), min_margin=r['min_chronology_margin'],
            energy_drift=r['energy_drift'], angular_drift=r['angular_drift']))
    return rows


def antipodal_table(runs, history):
    lines = ['| birth | rule | birth speeds | contacts (t, place) | captures | throat exits (advanced) | '
             'max \\|offset\\| | min chronology margin |', '|---|---|---|---|---|---|---|---|']
    cont = {'symmetric': pick(history, tag='history_symmetric', v=.30), 'symmetric_third':
            pick(history, tag='history_symmetric', v=1/3), 'independent': pick(history, tag='history_independent'),
            'lopsided': pick(history, tag='history_lopsided')} if history else {}
    order = ['symmetric', 'symmetric_third', 'independent', 'lopsided', 'fast_lopsided']
    rows = sorted(antipodal_rows(runs), key=lambda q: (order.index(q['birth']), q['static'], q['handle']))
    for name in order:
        label = ANTIPODAL_NAMES[name]
        if name in cont:
            h = cont[name]
            cs = ', '.join(f"{e['time']:.2f} {place(e['point'])}" for e in h['events'] if e['kind'] == 'contact')
            ext = np.abs(h['clock'][h['alive']])
            lines.append(f"| {label} | continuous (#345 §7) | — | {cs} | — | — | "
                         f"{np.nanmax(ext):.3f} | {h['min_chronology_margin']:.3f} |")
            label = ''
        for q in (q for q in rows if q['birth'] == name):
            rule = ('antipodal, strict' if q['static'] == 'none' else 'antipodal + static field') + \
                (', torus' if q['handle'] == 'torus' else '')
            cs = ', '.join(f'{t:.2f} {p}' for t, p in q['contacts'])
            lines.append(f"| {label} | {rule} | {q['speeds'][0]:.3f} / {q['speeds'][1]:.3f} | {cs} | {q['captures']} | "
                         f"{q['exits']} ({q['advanced']}) | {q['max_offset']:.3f} | {q['min_margin']:.3f} |")
            label = ''
    return '\n'.join(lines)


def figure_antipodal(runs, history, path):
    fig, axes = plt.subplots(2, 2, figsize=(15, 9), gridspec_kw=dict(hspace=.42, wspace=.22))
    by = {r['tag']: r for r in runs}
    styles = dict(continuous=('#7f8c8d', 'continuous (#345)'), none=('#c0392b', 'antipodal, strict'),
                  coulomb=('#1f5fa8', 'antipodal + static field'))

    def alive_line(ax, t, y, alive, **kw):
        ax.plot(np.where(alive, t, np.nan), np.where(alive, y, np.nan), **kw)

    ax = axes[0, 0]
    if history:
        h = pick(history, tag='history_symmetric', v=.30)
        alive_line(ax, h['t'], h['d'], h['alive'], color=styles['continuous'][0], lw=1.2, label=styles['continuous'][1])
    for static in ('none', 'coulomb'):
        r = by[f'antipodal_symmetric_{static}']
        alive_line(ax, r['t'], r['d'], r['alive'], color=styles[static][0], lw=1.2, label=styles[static][1])
        for e in r['events']:
            if e['kind'] == 'contact':
                ax.text(e['time'], 2*A - (.25 if static == 'none' else .45), place(e['point']), ha='center',
                        fontsize=7, color=styles[static][0])
    ax.axhline(np.pi, color='k', lw=.4)
    ax.set_ylim(0, 3.3)
    ax.set_xlabel('time t')
    ax.set_ylabel('separation d')
    ax.set_title('(a) Symmetric v = 0.30: under the strict rule the pair coasts to −P every cycle;\n'
                 'the static field alone restores #345\'s P, −P, P, ...', fontsize=10)
    ax.legend(fontsize=8, loc='upper right')
    ax = axes[0, 1]
    if history:
        h = pick(history, tag='history_lopsided')
        ax.plot(h['t'], h['clock'], color=styles['continuous'][0], lw=1.2, label=styles['continuous'][1])
    for static in ('none', 'coulomb'):
        r = by[f'antipodal_lopsided_{static}']
        ax.plot(r['t'], r['clock'], color=styles[static][0], lw=1.2, label=styles[static][1])
    L = by['antipodal_lopsided_none']['length']
    for s_ in (1, -1):
        ax.axhline(s_*L, color='k', lw=.6, ls='--')
    ax.text(.5, -L - .05, '−L: beyond it a throat exit lands in the past', fontsize=7, va='top')
    ax.axhline(0, color='k', lw=.4)
    ax.set_xlabel('time t')
    ax.set_ylabel('clock offset A = τ_A − τ_B')
    ax.set_title('(b) Lopsided birth: without the attraction nothing swaps the roles,\nso the offset grows past L',
                 fontsize=10)
    ax.legend(fontsize=8, loc='upper right')
    ax = axes[1, 0]
    r = by['antipodal_lopsided_none']
    rows = dict(reflection=0, capture=1, throat_exit=2, contact=3, birth=4)
    for e in r['events']:
        y = rows[e['kind']]
        if e['kind'] == 'reflection':
            ax.plot(e['time'], y, '|', color='#7f8c8d', ms=7)
        elif e['kind'] == 'capture':
            ax.plot(e['time'], y, 'o', color=('#117a3a', '#7d2fa5')[e['mouth']], ms=6)
        elif e['kind'] == 'throat_exit':
            ax.plot(e['time'], y, '^', color=('#117a3a', '#7d2fa5')[e['mouth']], ms=7,
                    mfc='white' if e['advanced'] else None, mec='#c0392b' if e['advanced'] else None)
            if e['advanced']:
                ax.annotate('', xy=(e['time'], y + .15), xytext=(e['t_in'], rows['capture'] + .15),
                            arrowprops=dict(arrowstyle='->', color='#c0392b', lw=.8))
        else:
            ax.plot(e['time'], y, '*' if e['kind'] == 'contact' else 'v', color='k', ms=9)
    ax.set_yticks(list(rows.values()), ['reflection at an empty focus', 'capture at a focus (A green, B purple)',
                                        'throat exit (hollow: advanced)', 'rim contact, release', 'birth'], fontsize=8)
    ax.set_ylim(-.6, 4.6)
    ax.set_xlabel('time t')
    ax.set_title('(c) Lopsided birth, strict rule: every interaction is an antipodal focus;\n'
                 'red arrows join a capture to its exit in the past', fontsize=10)
    ax = axes[1, 1]
    for static in ('none', 'coulomb'):
        r = by[f'antipodal_fast_lopsided_{static}']
        margin = np.where(r['alive'], r['d'] + r['length'] - np.abs(r['clock']), np.nan)
        ax.plot(r['t'], margin, color=styles[static][0], lw=1.2, label=f'{styles[static][1]}: d + L − |A|')
    ax.axhline(0, color='k', lw=.6)
    ax.fill_between([0, HISTORY_END], -3, 0, color='#c0392b', alpha=.07)
    ax.set_ylim(-2, 2)
    ax.set_xlim(0, HISTORY_END)
    ax.set_xlabel('time t')
    ax.set_ylabel('chronology margin')
    ax.set_title('(d) Fast lopsided birth (0.85 / 0.15): the offset outgrows d + L;\n'
                 'below zero a closed null loop through the throat is possible', fontsize=10)
    ax.legend(fontsize=8, loc='upper right')
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


KLEIN_COLORS = dict(klein='#8e44ad', torus='#16a085')
KLEIN_ORDER = ['symmetric', 'symmetric_third', 'independent', 'lopsided', 'fast_lopsided']


def klein_rows(runs):
    rows = []
    for r in runs:
        ev = r['events']
        caps = [e for e in ev if e['kind'] == 'capture']
        exits = [e for e in ev if e['kind'] == 'throat_exit']
        i = int(np.searchsorted(r['t'], .1))
        rows.append(dict(
            tag=r['tag'], birth=r['birth_name'], handle=r['handle'], scale=r['flash_scale'], sweep=r['sweep'],
            end=r['end'], speeds=r['birth_speed'],
            contacts=[(e['time'], place(e['point'])) for e in ev if e['kind'] == 'contact'],
            push=sum(e['coupling'] > 0 for e in caps), pull=sum(e['coupling'] < 0 for e in caps),
            quadrature=sum(e['coupling'] == 0 for e in caps), exits=len(exits),
            capped=sum(bool(e.get('saturated')) for e in caps + exits),
            capture_drate=float(sum(e['drate'] for e in caps)), exit_drate=float(sum(e['drate'] for e in exits)),
            birth_rate=float((r['d'][i] - r['d'][0])/(r['t'][i] - r['t'][0])),
            flash=float(r['energy']['fronts'][0]),
            max_kick=float(max([np.linalg.norm(e['kick']) for e in caps] or [0])),
            energy_drift=r['energy_drift'], angular_drift=r['angular_drift']))
    return rows


def klein_table(runs):
    lines = ['| birth | throat | birth speeds | contacts (t, place) | captures: push / pull / quadrature | '
             'Σ Δḋ at captures | Σ Δḋ at exits | net Σ Δḋ | ḋ at birth |',
             '|---|---|---|---|---|---|---|---|---|']
    rows = [q for q in klein_rows(runs) if not q['sweep']]
    for name in KLEIN_ORDER:
        label = ANTIPODAL_NAMES[name]
        for q in sorted((q for q in rows if q['birth'] == name), key=lambda q: q['handle']):
            cs = ', '.join(f'{t:.2f} {p}' for t, p in q['contacts'])
            capped = f" ({q['capped']} capped)" if q['capped'] else ''
            lines.append(f"| {label} | {q['handle']} | {q['speeds'][0]:.3f} / {q['speeds'][1]:.3f} | {cs} | "
                         f"{q['push']} / {q['pull']} / {q['quadrature']}{capped} | {q['capture_drate']:+.4f} | "
                         f"{q['exit_drate']:+.4f} | {q['capture_drate'] + q['exit_drate']:+.4f} | "
                         f"{q['birth_rate']:.3f} |")
            label = ''
    return '\n'.join(lines)


def klein_sweep_table(runs):
    lines = ['| birth | flash × | throat | birth speeds | flash energy | contacts to t = 30 | '
             'captures: push / pull / quadrature | net Σ Δḋ |', '|---|---|---|---|---|---|---|---|']
    rows = sorted((q for q in klein_rows(runs) if q['sweep']), key=lambda q: (q['birth'], q['scale'], q['handle']))
    for q in rows:
        cs = ', '.join(f'{t:.2f} {p}' for t, p in q['contacts'])
        lines.append(f"| {ANTIPODAL_NAMES[q['birth']]} | {q['scale']:g} | {q['handle']} | {q['speeds'][0]:.3f} / "
                     f"{q['speeds'][1]:.3f} | {q['flash']:.4f} | {cs} | "
                     f"{q['push']} / {q['pull']} / {q['quadrature']} | "
                     f"{q['capture_drate'] + q['exit_drate']:+.4f} |")
    return '\n'.join(lines)


def handle_turns():
    """A front in quadrature captured at the south pole and leaving the north pole, for each handle
    and each component of its angular momentum: the factor the throat applies (cf. the tests)."""
    from geometrodynamics.waves.antipodal_pair_birth import AntipodalPairBirth, AntipodalState, Front
    from geometrodynamics.waves.pair_birth import birth_frames
    out = {}
    for handle in ('klein', 'torus'):
        m = AntipodalPairBirth(Parameters(cutoff=12, footprint=.2), handle=handle, exchange='signed')
        factors = []
        for axis in range(3):
            x = np.array([[0., 0., 1.], [0., 0., -1.]])
            s = AntipodalState(0., x, np.zeros((2, 3)), np.zeros(2), birth_frames(x, np.array([0., 1., 0.])),
                               np.zeros(2))
            J = .002*np.eye(3)[axis]
            s.fronts = [Front(x[0].copy(), 0., .01, J.copy(), -1j, 'probe', 1, 1.)]
            _, info = m.run(s, np.pi + m.length + .05, dt=.01)
            f = next(f for f in info['final_state'].fronts if f.kind.startswith('exit'))
            factors.append(float(f.J @ J/(J @ J)))
        out[handle] = factors + [float(f.sign)]
    return out


def figure_klein(runs, path):
    fig, axes = plt.subplots(2, 2, figsize=(15, 9), gridspec_kw=dict(hspace=.45, wspace=.22))
    by = {r['tag']: r for r in runs}

    def cumulative(r):
        ev = [e for e in r['events'] if e['kind'] in ('capture', 'throat_exit')]
        ev = sorted(ev, key=lambda e: e['time'])
        t = np.r_[0., [e['time'] for e in ev], r['end']]
        y = np.r_[0., np.cumsum([e['drate'] for e in ev]), 0.]
        y[-1] = y[-2]
        return t, y

    ax = axes[0, 0]
    for name, ls in (('symmetric', '-'), ('symmetric_third', '--')):
        for handle in ('klein', 'torus'):
            r = by[f'klein_{name}_{handle}']
            t, y = cumulative(r)
            ax.step(t, y, where='post', color=KLEIN_COLORS[handle], ls=ls, lw=1.4,
                    label=f"{handle}, {ANTIPODAL_NAMES[name]}")
            for e in r['events']:
                if e['kind'] == 'contact':
                    ax.plot(e['time'], .012 if handle == 'klein' else .016, 'v', color=KLEIN_COLORS[handle],
                            ms=5, alpha=.8 if ls == '-' else .4)
    ax.axhline(0, color='k', lw=.5)
    ax.set_ylim(-.06, .02)
    ax.set_xlabel('time t')
    ax.set_ylabel('accumulated change of separation rate Σ Δḋ')
    ax.set_title('(a) Symmetric births: Klein exchange accumulates approach, torus about none;\n'
                 'stopping the pair would need Σ Δḋ ≈ −0.57 (triangles: contacts at −P)', fontsize=10)
    ax.legend(fontsize=8, loc='lower left')

    ax = axes[0, 1]
    r_k, r_t = by['klein_symmetric_klein'], by['klein_symmetric_torus']
    for r, dx in ((r_k, -.25), (r_t, .25)):
        c = KLEIN_COLORS[r['handle']]
        for e in r['events']:
            if e['kind'] not in ('capture', 'throat_exit') or e.get('mouth') != 0:
                continue                                  # the run is mirror-symmetric: show mouth A's events
            t = e['time'] + dx
            ax.plot([t, t], [0, e['drate']], color=c, lw=1)
            if e['kind'] == 'capture':
                face = c if e['coupling'] > 0 else 'white'
                ax.plot(t, e['drate'], 'o' if e['coupling'] != 0 else 'x', color=c, mfc=face, ms=6)
            else:
                ax.plot(t, e['drate'], '^', color=c, ms=7)
    ax.axhline(0, color='k', lw=.5)
    ax.annotate('the first exit front comes back:\ntorus: same sign, pushes\nKlein: flipped, pulls, and the\n'
                'throat then holds 2J', xy=(25.8, -.0061), xytext=(11.5, -.0068), fontsize=8, va='top',
                arrowprops=dict(arrowstyle='->', lw=.7))
    ax.plot([], [], 'o', color='k', ms=6, label='capture, push')
    ax.plot([], [], 'o', color='k', mfc='white', ms=6, label='capture, pull')
    ax.plot([], [], 'x', color='k', ms=6, label='capture in quadrature')
    ax.plot([], [], '^', color='k', ms=7, label='throat exit (recoil)')
    for handle in ('klein', 'torus'):
        ax.plot([], [], color=KLEIN_COLORS[handle], lw=3, label=handle)
    ax.set_ylim(-.0125, .0036)
    ax.set_xlabel('time t (Klein drawn left, torus right)')
    ax.set_ylabel('Δḋ at the event')
    ax.set_title('(b) Symmetric v = 0.30, each exchange at mouth A (B mirrors it)', fontsize=10)
    ax.legend(fontsize=7, loc='lower right', ncol=2)

    ax = axes[1, 0]
    turns = handle_turns()
    labels = ['orbital,\nmotion along e1\n(the birth circle)', 'orbital,\nmotion along e2',
              'spin about\nthe mouth axis', 'amplitude sign\n(compression)']
    xs = np.arange(4)
    for i, handle in enumerate(('klein', 'torus')):
        ax.bar(xs + (i - .5)*.36, turns[handle], .34, color=KLEIN_COLORS[handle], label=handle)
    ax.axhline(0, color='k', lw=.5)
    ax.set_xticks(xs, labels, fontsize=8)
    ax.set_yticks([-1, 0, 1], ['reversed', '', 'kept'])
    ax.set_title('(c) What one passage does (mouths antipodal, measured): each handle turns J by π;\n'
                 'Klein flips spin and sign, torus the sideways orbital part', fontsize=10)
    ax.legend(fontsize=8, loc='lower right')

    ax = axes[1, 1]
    rows = [q for q in klein_rows(runs) if q['sweep']]
    for name, marker in (('symmetric', 'o'), ('lopsided', 's')):
        for handle, dx in (('klein', -.04), ('torus', .04)):
            q = sorted((q for q in rows if q['birth'] == name and q['handle'] == handle), key=lambda q: q['scale'])
            sc = np.array([q_['scale'] for q_ in q])
            ax.plot(sc*(1 + dx), [q_['capture_drate'] + q_['exit_drate'] for q_ in q], marker=marker, ls='-',
                    color=KLEIN_COLORS[handle], mfc=KLEIN_COLORS[handle] if name == 'symmetric' else 'white',
                    label=f'{handle}, {ANTIPODAL_NAMES[name]}')
    ax.axhline(0, color='k', lw=.5)
    ax.set_xscale('log', base=2)
    ax.set_xticks([1, 2, 4], ['1', '2', '4'])
    ax.set_xlabel('flash scale (the pair pays for it: birth speed falls)')
    ax.set_ylabel('net Σ Δḋ to t = 30')
    ax.set_title('(d) A stronger flash: no trend, no binding; the two handles differ\n'
                 'only where a returning exit front is captured in phase', fontsize=10)
    ax.legend(fontsize=7, loc='upper left')
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def load_elliptic():
    path = RUNS/'elliptic.json'
    return json.loads(path.read_text()) if path.exists() else None


def elliptic_rows(data):
    rows = []
    for c in data['cycles']:
        caps = {}
        for e in c['captures']:
            i = sum(b['time'] <= e['time'] + 1e-9 for b in c['births']) - 1
            caps.setdefault(i, []).append(e)
        rows.append(dict(label=c['label'], kick=c['kick'], speeds=[b['speed'] for b in c['births']],
                         flash=[b['flash'] for b in c['births']], contacts=c['contacts'],
                         fronts_after_contact=c['fronts_after_contact'],
                         captures_before_contact=[sum(e['time'] <= tc + 1e-9 and e['time'] >= (c['births'][i]['time'])
                                                      for e in caps.get(i, [])) for i, tc in enumerate(c['contacts'])],
                         energy_drift=c['energy_drift'], angular_drift=c['angular_drift']))
    return rows


ELLIPTIC_NAMES = dict(resonant_k1='resonant, k = 1 (u₁ less 0.02%)', resonant_k2='resonant, k = 2 (u₂ less 0.02%)',
                      off_resonance_v030='off resonance (kick 0.30)')


def elliptic_table(data):
    lines = ['| cycle | birth speeds, cycles 1–4 | flash at each birth | contacts | flash captures before each '
             'contact | fronts in flight after each contact |', '|---|---|---|---|---|---|']
    for q in elliptic_rows(data):
        lines.append(f"| {ELLIPTIC_NAMES[q['label']]} | {', '.join(f'{v:.5f}' for v in q['speeds'])} | "
                     f"{', '.join(f'{f:.5f}' for f in q['flash'])} | {', '.join(f'{t:.3f}' for t in q['contacts'])} | "
                     f"{', '.join(str(n) for n in q['captures_before_contact'])} | "
                     f"{', '.join(str(n) for n in q['fronts_after_contact'])} |")
    return '\n'.join(lines)


def figure_elliptic(data, path):
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), gridspec_kw=dict(hspace=.38, wspace=.22))
    by = {c['label']: c for c in data['cycles']}

    def theta(x):
        x = np.asarray(x)
        return np.arctan2(x[..., 0], x[..., 2])

    def flight(ax, c, title):
        h = c['history']
        t, x, alive = np.asarray(h['t']), np.asarray(h['x']), np.asarray(h['alive'])
        tc = c['contacts'][0]
        keep = (t <= tc + 1e-9) & alive
        for k, (col, name) in enumerate((('#117a3a', 'A (+q)'), ('#7d2fa5', 'B (−q)'))):
            ax.plot(theta(x[keep, k]), t[keep], color=col, lw=2.6, label=f'mouth {name}')
        for k in range(2):                                      # each mouth's image lies on the other
            ax.plot(np.angle(np.exp(1j*(theta(x[keep, k]) + np.pi))), tc - t[keep], color='k', lw=1., ls='--',
                    label='image of the other mouth under (n, t, q) → (−n, t_c − t, −q)' if k == 0 else None)
        for n in range(1, int(tc/np.pi) + 2):
            for k, col in enumerate(('#117a3a', '#7d2fa5')):
                th0 = theta(x[0, k]) + (np.pi if n % 2 else 0.)
                th0 = np.angle(np.exp(1j*th0))
                ax.plot(th0, n*np.pi, 'o', mfc='white', mec=col, ms=8, mew=1.6)
        for e in c['captures']:
            if e['time'] <= tc + 1e-9:
                ax.plot(theta(x[keep, e['mouth']])[-1], e['time'], 'o', color=('#117a3a', '#7d2fa5')[1 - e['mouth']],
                        mec='k', ms=10, zorder=6)
        ax.axhline(tc/2, color='#b9770e', lw=.8, ls=':')
        ax.text(-3.05, tc/2 + .1, 'mid-flight: mouths antipodal, the involution’s fixed slice', fontsize=7, color='#b9770e')
        ax.axhline(tc, color='k', lw=.6)
        ax.text(-3.05, tc + .1, f'contact at −P, t_c = {tc:.3f} = {tc/np.pi:.3f} π', fontsize=7)
        ax.set_xlim(-np.pi, np.pi)
        ax.set_ylim(0, max(tc, 3*np.pi) + 2.4)
        ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi], ['−P', '−π/2', 'P', 'π/2', '−P'])
        ax.set_xlabel('position on the birth circle')
        ax.set_ylabel('time t')
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7, loc='upper center')

    flight(axes[0, 0], by['off_resonance_v030'],
           '(a) Off resonance (birth speed 0.283): the mouths are each other’s image,\n'
           'but each flash focuses (open circles: A’s green, B’s purple) away from contact')
    flight(axes[0, 1], by['resonant_k1'],
           '(b) Resonant u₁: each flash focuses on the partner at contact (filled, the flash’s colour):\n'
           'emission at birth and absorption at contact are images of each other')

    ax = axes[1, 0]
    scan = data['scan']
    u = np.array([r['speed'] for r in scan])
    mis = np.array([r['flash_mismatch'] for r in scan])
    early = np.array([any(abs(c['time'] - r['contact']) > 1e-3 for c in r['captures']) for r in scan])
    ax.plot(u, mis, '-', color='#7f8c8d', lw=1)
    ax.plot(u[~early], mis[~early], 'o', color='#1f5fa8', ms=4, label='no flash capture before contact')
    ax.plot(u[early], mis[early], 's', color='#c0392b', ms=4, label='a flash captured before contact, not at it')
    from geometrodynamics.waves.antipodal_pair_birth import closure_band
    for k, uk in enumerate(data['meta']['resonant_speeds']):
        if uk > u.min() - .01:
            low, _ = closure_band(A, k)
            ax.axvspan(low, uk, color='#117a3a', alpha=.18, lw=0,
                       label='cycle closes: the caught flash is still in the throat at contact' if k == 0 else None)
            ax.axvline(uk, color='k', lw=.6, ls='--')
            ax.text(uk, 3.25, f'u{"₀₁₂₃"[k]}', ha='center', fontsize=9)
    ax.set_ylim(0, 3.45)
    ax.set_xlabel('birth speed u (after paying for the flash)')
    ax.set_ylabel('|t_c − nearest odd multiple of π|')
    ax.set_title('(c) Emission and absorption are exact images only at u_k = (π − 2a)/((2k + 1)π);\n'
                 'cycles close just below it (shaded); trajectories are invariant at every speed', fontsize=10)
    ax.legend(fontsize=8, loc='upper center', bbox_to_anchor=(.5, -.13), ncol=2)

    ax = axes[1, 1]
    st = data['static']
    E = np.geomspace(1e-3, 3, 200)
    us = st['birth_speed']
    for slope, style, label in ((1., '-', 'one front absorbed, fully directed: Δp ≤ E'),
                                (2., '--', 'one in, one out, both directed (a mirror): Δp ≤ 2E'),
                                (2*us, ':', f'one in, one out, emitted isotropically at u = {us:.3f}: Δp ≈ 2uE')):
        ax.plot(E, slope*E, style, color='#1f5fa8', lw=1.2, label=label)
    need = st['momentum_change_per_mouth']
    ax.axhline(need, color='#c0392b', lw=1)
    ax.text(1.1e-3, need*1.15, f'turning a mouth around: Δp = 2γMu = {need:.3f}', fontsize=8, color='#c0392b')
    ax.plot(st['energy_exchanged_per_mouth'], need, '*', color='#c0392b', ms=14,
            label=f"static field: Δp = {need:.3f} for ΔE = {st['energy_exchanged_per_mouth']:.4f}")
    ef = st['flash_energy_per_mouth']
    ax.plot([ef, ef], [2*us*ef, 2*ef], color='#117a3a', lw=3, solid_capstyle='butt')
    ax.plot(ef, 2*ef, 'o', color='#117a3a', ms=7, label=f'the birth flash, {ef:.4f} per mouth: at most {2*ef:.4f}')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-3, 3)
    ax.set_ylim(1e-4, 3)
    ax.set_xlabel('energy exchanged per mouth per flight')
    ax.set_ylabel('momentum moved per mouth')
    ax.set_title('(d) Why fronts cannot do what the static field does: it moves 14× more momentum\n'
                 'than energy; fronts move at most 2× (two-way), or 2u× if emitted at rest', fontsize=10)
    ax.legend(fontsize=7, loc='lower right')
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def load_s3():
    path = RUNS/'s3.json'
    return json.loads(path.read_text()) if path.exists() else None


def s3_tables(data):
    bands = ['| k | u_k | S2 closure band (neck sin 2a) | S3 closure band (3-D neck) |', '|---|---|---|---|']
    for b in data['bands']:
        bands.append(f"| {b['k']} | {b['s2'][1]:.5f} | ({b['s2'][0]:.5f}, {b['s2'][1]:.5f}] | "
                     f"({b['s3'][0]:.5f}, {b['s3'][1]:.5f}] |")
    tally = ['| space | per focus | handle | kick | birth speed | Σ Δḋ at captures | Σ Δḋ at exits (exits) | net | '
             'contacts |', '|---|---|---|---|---|---|---|---|---|']
    for r in data['tally']:
        sign = r.get('per_focus', '-1' if r['space'] == 'S3' else '-i').replace('-', '−')
        tally.append(f"| {r['space']} | {sign} | {r['handle']} | {r['kick']:.3f} | {r['birth_speed']:.3f} | "
                     f"{r['capture_drate']:+.4f} | {r['exit_drate']:+.4f} ({r['exits']}) | "
                     f"{r['capture_drate'] + r['exit_drate']:+.4f} | {', '.join(f'{t:.2f}' for t in r['contacts'])} |")
    signs = ['| data | energy off U(T) z = K z, S = (−1)^m | with S = −(−1)^m | momentum density at the foci / energy |',
             '|---|---|---|---|']
    for name in dict.fromkeys(r['data'] for r in data.get('signs', [])):
        rows = [r for r in data['signs'] if r['data'] == name]
        signs.append(f"| {name} | {max(r['same'] for r in rows):.2g} | {max(r['opposite'] for r in rows):.2g} | "
                     f"{max(r['focal_momentum'] for r in rows):.2g} |")
    return '\n'.join(bands), '\n'.join(tally), '\n'.join(signs)


def figure_s3(data, path):
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), gridspec_kw=dict(hspace=.38, wspace=.22))
    pr = data['profiles']
    psi = np.array(pr['psi'])
    ax = axes[0, 0]
    u0 = np.array(pr['u0'])
    ax.plot(psi, u0, color='#7f8c8d', lw=1.5, label='flash at t = 0 (released at rest, centred on the emitter)')
    ax.plot(psi, -u0[::-1], color='k', lw=4, alpha=.25, label='its antipodal image, sign flipped')
    ax.plot(psi, pr['s3_exact'], color='#1f5fa8', lw=1.5, label='S³ at t = π, exact')
    ax.plot(psi[::8], np.array(pr['s3_fd'])[::8], 'o', color='#1f5fa8', ms=3, label='S³ at t = π, finite volumes')
    ax.plot(psi, pr['s2_u'], color='#c0392b', lw=1.5, label='S² (conformal) at t = π: displacement')
    p2 = np.array(pr['s2_p'])
    ax.plot(psi, p2/np.abs(p2).max(), color='#c0392b', lw=1, ls='--', label='S² at t = π: velocity (scaled)')
    ax.set_xticks([0, np.pi/2, np.pi], ['emitter', 'equator', 'antipode'])
    ax.set_xlabel('distance from the emitter ψ (zonal flash)')
    ax.set_ylabel('field')
    ax.set_title('(a) After π the S³ flash is its own mirror image, sign flipped and still at rest;\n'
                 'on S² the focus turns all displacement into velocity (the quarter phase)', fontsize=10)
    ax.legend(fontsize=7, loc='lower left')

    ax = axes[0, 1]
    k = np.array([c['k'] for c in data['closure']])
    s3 = np.maximum([c['s3'] for c in data['closure']], 1e-32)
    s2 = [c['s2_conformal'] for c in data['closure']]
    ax.semilogy(k, s3, 'o-', color='#1f5fa8', label='S³ conformal, released at rest, S = (−1)^m')
    if data.get('signs'):
        for name, marker, ls, label in (('at rest (P = 0)', 'o', ':', 'S³, released at rest, opposite sign'),
                                        ('birth phase 45 degrees', 'D', '-.', 'S³, 45° birth phase, either sign')):
            rows = sorted((r for r in data['signs'] if r['data'] == name), key=lambda r: r['k'])
            ax.semilogy([r['k'] for r in rows], [r['opposite'] for r in rows], marker + ls, color='#1f5fa8',
                        mfc='white', label=label)
    ax.semilogy(k, s2, 's-', color='#c0392b', label='S² conformal, released at rest')
    ax.semilogy(k, [.5923, .2427, .8686, .0939], 'v--', color='#c0392b', mfc='white',
                label='S² Maxwell flash, released at rest (#349)')
    ax.semilogy(k, [.0017, .0152, .0417, .0800], '^--', color='#b9770e', mfc='white',
                label='S² Maxwell flash, best birth phase (#349 review)')
    ax.set_xticks(k, [f'k = {i}\n(T = {2*i + 1}π)' for i in k])
    ax.set_ylim(1e-32, 3)
    ax.set_ylabel('energy fraction to change for U(T) z = K z')
    ax.set_title('(b) Time-reversal closure of the flash: on S³ exact or wholly off, by the sign of the\n'
                 'identification (not yet derived); on S² half off (conformal) or dispersed (Maxwell)', fontsize=10)
    ax.legend(fontsize=7, loc='center right')

    ab = data['absorption']
    a = data['meta']['mouth_angle']

    def series(**kw):
        rows = sorted((r for r in ab if all(r[k] == v for k, v in kw.items())), key=lambda r: r['delta'])
        return np.array([r['delta'] for r in rows])/a, rows

    ax = axes[1, 0]
    for sphere, col in ((3, '#1f5fa8'), (2, '#c0392b')):
        x, rows = series(sphere=sphere, speed=0., offset='along')
        ax.plot(x, [r['energy_fraction'] for r in rows], color=col, lw=1.6, label=f'S{"³" if sphere == 3 else "²"}: absorbing mouth')
    ax.plot([0, 1, 1, x.max()], [1, 1, 0, 0], color='k', lw=1, ls=':', label='strict rule: capture iff the focus is inside')
    ax.set_xlabel('offset of the focus from the mouth centre, δ/a (a = 0.3)')
    ax.set_ylabel('fraction of the flash the mouth takes')
    ax.set_title('(c) A mouth near the focus takes rays on both sides of it: no jump at the rim;\n'
                 'beyond it S³ takes 1 − cos α, sin α = sin a / sin δ', fontsize=10)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    for sphere, offset, col, ls in ((3, 'along', '#1f5fa8', '-'), (3, 'across', '#1f5fa8', '--'),
                                    (2, 'along', '#c0392b', '-'), (2, 'across', '#c0392b', '--')):
        x, rows = series(sphere=sphere, speed=.27, offset=offset)
        ax.plot(x, [r['along'] for r in rows], color=col, ls=ls, lw=1.5,
                label=f'S{"³" if sphere == 3 else "²"}, offset {offset} the line of motion')
    x, rows = series(sphere=3, speed=0., offset='along')
    ax.plot(x, [r['along'] for r in rows], color='k', lw=2.5, label='emitter at rest (S³ and S², either offset): exactly 0')
    ax.axvline(1, color='k', lw=.6, ls=':')
    ax.set_xlabel('offset of the focus from the mouth centre, δ/a')
    ax.set_ylabel('momentum taken / γE₀u')
    ax.set_title('(d) The mouth takes the emitter’s momentum, never more: all of it while it holds the focus\n'
                 '(across: cos δ as orbital motion, the rest as spin); none from an emitter at rest', fontsize=10)
    ax.legend(fontsize=7, loc='upper right')
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


HISTORY_END = 60.


def main():
    scan, featured, retimed, checks = (load(p) for p in ('scan', 'featured', 'retimed', 'checks'))
    history = load_history()
    neck, neck_checks = load_neck('neck'), load_neck('neck_checks')
    if neck and neck_checks:
        figure_neck(neck, neck_checks, IMAGES/'neck.png')
    antipodal = load_antipodal()
    if antipodal:
        figure_antipodal(antipodal, history, IMAGES/'antipodal.png')
    klein = load_antipodal('klein')
    if klein:
        figure_klein(klein, IMAGES/'klein.png')
    elliptic = load_elliptic()
    if elliptic:
        figure_elliptic(elliptic, IMAGES/'elliptic.png')
    s3 = load_s3()
    if s3:
        figure_s3(s3, IMAGES/'s3.png')
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
        neck=dict(rows=neck_rows(neck), feedback={b: feedback(neck, b) for b in sorted({r['birth'] for r in neck})
                                                  if not b.endswith('torus')},
                  torus=dict(history=[c['time'] for c in contacts(pick(neck, tag='neck_lopsided_torus_history'))],
                             off=[c['time'] for c in contacts(pick(neck, tag='neck_lopsided_torus_off'))]),
                  calibration=calibration_rows(neck_checks)) if neck and neck_checks else None,
        antipodal=dict(rows=antipodal_rows(antipodal)) if antipodal else None,
        klein=dict(rows=klein_rows(klein), turns=handle_turns()) if klein else None,
        elliptic=dict(rows=elliptic_rows(elliptic), static=elliptic['static'],
                      resonant_speeds=elliptic['meta']['resonant_speeds'],
                      max_trajectory_mismatch=max(r['trajectory_mismatch'] for r in elliptic['scan']))
        if elliptic else None,
        s3=dict(refocus=s3['refocus'], closure=s3['closure'], signs=s3.get('signs'), bands=s3['bands'],
                router=s3['router'], couplings=s3['couplings'], tally=s3['tally']) if s3 else None,
        tables=dict(outcomes=outcome_table(scan), ledger=ledger_table(scan), retimed=retimed_table(retimed, featured),
                    history=history_table(history_rows(history)) if history else None,
                    neck=neck_table(neck_rows(neck)) if neck else None,
                    calibration=calibration_table(calibration_rows(neck_checks)) if neck_checks else None,
                    antipodal=antipodal_table(antipodal, history) if antipodal else None,
                    klein=klein_table(klein) if klein else None,
                    klein_sweep=klein_sweep_table(klein) if klein else None,
                    elliptic=elliptic_table(elliptic) if elliptic else None,
                    s3_bands=s3_tables(s3)[0] if s3 else None, s3_tally=s3_tables(s3)[1] if s3 else None,
                    s3_signs=s3_tables(s3)[2] if s3 else None))
    (RUNS/'summary.json').write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    files = sorted(RUNS.glob('*.json')) + sorted(p for ext in ('png', 'gif', 'mp4') for p in IMAGES.glob(f'*.{ext}'))
    (RUNS/'sha256.json').write_text(json.dumps({p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                                for p in files if p.name != 'sha256.json'}, indent=2)+'\n')
    for k in ('outcomes', 'ledger', 'retimed', 'history', 'neck', 'calibration', 'antipodal', 'klein',
              'klein_sweep', 'elliptic', 's3_bands', 's3_tally'):
        print(report['tables'][k], '\n')
    if report['neck']:
        print(json.dumps(dict(feedback=report['neck']['feedback'], torus=report['neck']['torus']), indent=1))
    print(json.dumps(dict(thresholds=report['thresholds'], barrier=report['coulomb_barrier'],
                          first_returns=stats, direct=report['direct_exchange'],
                          max_energy_drift=report['max_energy_drift'],
                          max_mirror_asymmetry=report['max_mirror_asymmetry']), indent=1))


if __name__ == '__main__':
    main()
