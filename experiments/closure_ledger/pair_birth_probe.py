#!/usr/bin/env python
"""Pair birth to annihilation for two Wheeler-flux mouths: speed scan, re-timed flashes, checks.

    python experiments/closure_ledger/pair_birth_probe.py scan       # speed x {full, dressed, instant}
    python experiments/closure_ledger/pair_birth_probe.py featured   # v = 0.30, 1/3 full runs with field maps
    python experiments/closure_ledger/pair_birth_probe.py retimed    # same birth, flash removed or re-timed at t = 3
    python experiments/closure_ledger/pair_birth_probe.py checks     # step, cutoff, free motion, independent mouths
    python experiments/closure_ledger/pair_birth_probe.py history    # repeated cycles with history clocks
    python experiments/closure_ledger/pair_birth_probe.py neck       # neck channel: clock history, off, reversed
    python experiments/closure_ledger/pair_birth_probe.py neck_checks    # port stiffness, nodes, step
    python experiments/closure_ledger/pair_birth_probe.py antipodal  # the history under the antipodal-only rule
    python experiments/closure_ledger/pair_birth_probe.py klein      # signed exchange: Klein against torus throats
    python experiments/closure_ledger/pair_birth_probe.py elliptic   # antipodal time reversal, resonant speeds
    python experiments/closure_ledger/pair_birth_probe.py s3         # the S3 benchmark: focusing, closure, absorption

Each part writes JSON to experiments/closure_ledger/runs/20261019_pair_birth/.
Model: geometrodynamics/waves/pair_birth.py. Write-up: docs/pair_birth_timeline.md.
Runs use four worker processes, each with single-threaded BLAS.
"""
import json
import os
import sys
import time

for _var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_var, '1')
from dataclasses import replace
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from geometrodynamics.waves.antipodal_pair_birth import AntipodalPairBirth, elliptic_mismatch, resonant_speeds
from geometrodynamics.waves.pair_birth import (NeckChannel, PairBirth, Parameters, birth_positions, circle_angle,
                                               front_crossings)

OUT = Path(__file__).resolve().parent/'runs'/'20261019_pair_birth'
END = 34.
DT = .004
OUTPUT_DT = .02
SPEEDS = (.20, .25, .28, .30, .32, 1/3, .34, .36, .40, .45, .50)
RETIME = 3.0                   # the birth flash is near its antipodal foci, far from the mouths
SHIFTS = (None, -1., -.5, .5, 1.)
PLANE = ((0., 0., 1.), (1., 0., 0.))
P, T = (np.array(e) for e in PLANE)
N = np.cross(P, T)


def outward(h, key):
    """Component of a per-mouth force along each mouth's direction away from P on the circle."""
    e = np.stack([np.cross(N, h['x'][:, 0]), -np.cross(N, h['x'][:, 1])], 1)
    return np.einsum('tji,tji->tj', h[key], e)


RING_THETA = np.linspace(0, np.pi, 61)
RING_AZIMUTHS = 24


def ring_points():
    """Rings of constant distance theta from P, for ring-averaged field maps.

    The birth flash is dipolar and weak along the birth circle itself (its current axis), so
    the maps average each ring around P; the mouths sit on the ring at their distance from P.
    """
    az = 2*np.pi*np.arange(RING_AZIMUTHS)/RING_AZIMUTHS
    th, az = np.meshgrid(RING_THETA, az, indexing='ij')
    return (np.cos(th)[..., None]*P + np.sin(th)[..., None]*(np.cos(az)[..., None]*T + np.sin(az)[..., None]*N)
            ).reshape(-1, 3)


def thinning(h, info, every):
    keep = np.zeros(len(h['t']), bool)
    keep[::every] = True
    keep[-1] = True
    if info['contact_time'] is not None:
        keep[np.abs(h['t']-info['contact_time']) < 1e-12] = True      # just before and just after surgery
    return keep


def summarise(h, info, extra, every):
    """Thin a history to every `every` outputs and attach event and ledger data."""
    alive = h['energy_parts']['rest'] > 0
    phi = circle_angle(h['x'], P, T)
    keep = thinning(h, info, every)
    r = lambda a, k=7: np.round(np.asarray(a)[keep], k).tolist()
    out = dict(t=r(h['t'], 6), phi=r(phi, 6), d=r(h['d'], 6), alive=np.asarray(alive)[keep].tolist(),
               energy={k: r(v, 9) for k, v in h['energy_parts'].items()}, J=r(h['J'], 10))
    for key in ('coulomb', 'exchange', 'own'):
        f = outward(h, key)
        out[f'force_{key}'] = r(f[:, 0], 8)
        out[f'asymmetry_{key}'] = float(np.abs(f[:, 0]-f[:, 1]).max())
    if alive.sum() > 3:
        cr = front_crossings(phi[alive], h['t'][alive])
        if extra.get('retime') is not None:
            td, shift = extra['retime']
            cr = [c for c in cr if c['time'] < td]
            if shift is not None:
                cr += [c for c in front_crossings(phi[alive], h['t'][alive], advance=shift) if c['time'] >= td]
        out['crossings'] = cr
    info = {k: v for k, v in info.items() if k != 'final_state'}
    i = int(np.argmax(np.where(alive, h['d'], -1)))
    out.update(info, d_max=float(h['d'][i]), t_d_max=float(h['t'][i]),
               contact_side=None if info['contact_point'] is None
               else ('P' if np.dot(info['contact_point'], P) > 0 else '-P'), **extra)
    return out


def flash_force(m, h, s0, retime=None):
    """Post hoc outward force of the birth flash on mouth A, by source copy (own, partner).

    The flash is a free wave, so its force along the recorded trajectory follows from its
    initial data alone. With a re-timed run the flash acting after t_d is the shifted one.
    """
    Q0, P0 = m.flash(s0)
    e = np.cross(N, h['x'][:, 0])
    f = np.einsum('tsi,ti->ts', m.free_wave_force(h, Q0, P0)[:, 0], e)        # (time, source)
    if retime is not None:
        td, shift = retime
        late = h['t'] >= td - 1e-9
        f[late] = 0.
        if shift is not None:
            f2 = np.einsum('tsi,ti->ts', m.free_wave_force(h, *m.free_wave(Q0, P0, shift))[:, 0], e)
            f[late] = f2[late]
    return f


def job(args):
    tag, v, mode, kw = args
    par = Parameters(**kw.get('parameters', {}))
    if mode == 'instant':
        par = replace(par, transverse=False)
    m = PairBirth(par)
    retime = kw.get('retime')
    if kw.get('velocities') is not None:
        s0 = m.initial(velocities=kw['velocities'])
    else:
        s0 = m.initial(v, dressed=(mode == 'dressed'))
    run_kw = dict(dt=kw.get('dt', DT), output_dt=OUTPUT_DT, retime=retime)
    if kw.get('field_map'):
        run_kw.update(field_points=ring_points(), field_dt=.1)
    t0 = time.time()
    h, info = m.run(s0, kw.get('end', END), **run_kw)
    extra = dict(tag=tag, v=v, mode=mode, retime=retime, parameters=par.__dict__,
                 free_contact=(np.pi-2*par.mouth_angle)/v if v else None, seconds=round(time.time()-t0, 1))
    every = kw.get('every', 5)
    if kw.get('velocities') is not None:
        extra['velocities'] = np.asarray(kw['velocities']).tolist()
        extra['x'] = np.round(h['x'][thinning(h, info, every)], 6).tolist()
    out = summarise(h, info, extra, every)
    if mode == 'full':
        out['force_flash'] = np.round(flash_force(m, h, s0, retime)[thinning(h, info, every)], 8).tolist()
    if kw.get('field_map'):
        mp = h['maps']
        floor = 1e-7
        ring = lambda a: np.asarray(a).reshape(len(a), len(RING_THETA), RING_AZIMUTHS).mean(2)
        out['map'] = dict(t=np.round(mp['t'], 4).tolist(), theta=np.round(RING_THETA, 6).tolist(),
                          log_transverse=np.round(np.log10(ring(mp['transverse'])+floor), 2).tolist(),
                          log_deformation=np.round(np.log10(ring(mp['deformation'])+floor), 2).tolist(), floor=floor)
    return out


HISTORY_END = 60.
NECK = dict(port=100., nodes=65)    # port stiffness and neck nodes (calibration: neck_checks)


def history_job(args):
    """Several birth-to-annihilation cycles, with each mouth's proper time recorded. With a
    neck (a dict of NeckChannel fields) the offset acts back through the neck channel."""
    tag, v, cycles, velocities, every, *more = args
    neck, end, dt = (more + [None, HISTORY_END, DT][len(more):])[:3]
    m = PairBirth(Parameters(), None if neck is None else NeckChannel(**neck))
    s0 = m.initial(v) if velocities is None else m.initial(velocities=velocities)
    t0 = time.time()
    h, info = m.run_history(s0, end, cycles=cycles, dt=dt, output_dt=OUTPUT_DT)
    keep = np.zeros(len(h['t']), bool)
    keep[::every] = True
    keep[-1] = True
    for e in info['events']:                       # both sides of every surgery
        keep[np.abs(h['t']-e['time']) < 1e-9] = True
    def r(a, k=7):
        a = np.round(np.asarray(a, float)[keep], k)
        return [None if np.isnan(v_) else v_ for v_ in a.tolist()] if a.ndim == 1 else a.tolist()
    out = dict(tag=tag, v=v, mode='history' if neck is None else 'neck', retime=None, cycles=cycles, end=end,
               velocities=None if velocities is None else np.asarray(velocities).tolist(),
               t=r(h['t'], 6), d=r(h['d'], 6), x=np.round(h['x'][keep], 6).tolist(), tau=r(h['tau'], 8),
               alive=h['alive'][keep].tolist(), cycle=h['cycle'][keep].tolist(),
               energy={k: r(v_, 9) for k, v_ in h['energy_parts'].items()}, J=r(h['J'], 10),
               clock=r(h['clock'], 10), external=r(h['external'], 10), first_order=r(h['first_order'], 10),
               chronology_margin=r(h['chronology_margin'], 8),
               events=[{k: (np.asarray(v_).tolist() if isinstance(v_, np.ndarray) else v_) for k, v_ in e.items()}
                       for e in info['events']],
               seconds=round(time.time()-t0, 1),
               **{k: v_ for k, v_ in info.items() if k not in ('events', 'final_state')})
    out.update(contact_time=None, contact_side=None, energy_drift=info['energy_drift'])
    if neck is not None:
        norm = lambda a: np.linalg.norm(a, axis=-1)
        out.update(neck=neck, dt=dt, masses=r(h['masses'], 9), mu=r(h['mu'], 10), spin=r(h['spin'], 10),
                   ends=np.round(h['ends'][keep], 8).tolist(), throat_force=r(norm(h['throat']), 9),
                   port_force=r(norm(h['port']), 9), exchange_force=r(norm(h['exchange']), 9))
    return out


def jobs_for(part):
    if part == 'scan':
        return [('scan', v, mode, dict(every=10)) for v in SPEEDS for mode in ('full', 'dressed', 'instant')]
    if part == 'featured':
        return [('featured', v, 'full', dict(every=2, field_map=True)) for v in (.30, 1/3)]
    if part == 'retimed':
        return [('retimed', v, 'full', dict(every=3, retime=(RETIME, a))) for v in (.30, 1/3) for a in SHIFTS]
    if part == 'checks':
        x, vel, _, _ = birth_positions(.3, .3, PLANE)
        independent = vel.copy()
        independent[0] *= 1.15
        independent[1] = .85*vel[1] + .08*N
        return [('dt_half', 1/3, 'full', dict(dt=DT/2)),
                ('cutoff40', 1/3, 'full', dict(parameters=dict(cutoff=40))),
                ('dt_half', .30, 'full', dict(dt=DT/2)),
                ('cutoff40', .30, 'full', dict(parameters=dict(cutoff=40))),
                ('free', 1/3, 'full', dict(parameters=dict(charge=0.), end=12.)),
                ('independent', .3, 'full', dict(velocities=independent.tolist()))]
    if part == 'history':
        x, vel, _, _ = birth_positions(.3, .3, PLANE)
        independent = vel.copy()
        independent[0] *= 1.15
        independent[1] = .85*vel[1] + .08*N
        _, slow, _, _ = birth_positions(.3, .15, PLANE)
        _, fast, _, _ = birth_positions(.3, .45, PLANE)
        lopsided = np.array([fast[0], slow[1]])      # one mouth travels, the other nearly stays
        return [('history_symmetric', .30, 5, None, 5), ('history_symmetric', 1/3, 3, None, 5),
                ('history_independent', .30, 5, independent.tolist(), 5),
                ('history_lopsided', .30, 5, lopsided.tolist(), 5)]
    if part == 'antipodal':
        births = dict(history_births(), symmetric_third=None, fast_lopsided=None)
        _, slow, _, _ = birth_positions(.3, .15, PLANE)
        _, fast, _, _ = birth_positions(.3, .85, PLANE)
        births['fast_lopsided'] = np.array([fast[0], slow[1]]).tolist()
        out = []
        for name, vel in births.items():
            v = 1/3 if name == 'symmetric_third' else .30
            for static in ('none', 'coulomb'):
                out.append((f'antipodal_{name}_{static}', v, vel, static, 'klein'))
        return out + [('antipodal_lopsided_none_torus', .30, births['lopsided'], 'none', 'torus')]
    if part == 'klein':
        births = dict(history_births(), symmetric_third=None)
        _, slow, _, _ = birth_positions(.3, .15, PLANE)
        _, fast, _, _ = birth_positions(.3, .85, PLANE)
        births['fast_lopsided'] = np.array([fast[0], slow[1]]).tolist()
        out = [(f'klein_{name}_{handle}', 1/3 if name == 'symmetric_third' else .30, vel, 'none', handle, 'signed')
               for name, vel in births.items() for handle in ('klein', 'torus')]
        return out + [(f'klein_sweep_{name}_{handle}_x{scale:g}', .30, births[name], 'none', handle, 'signed', scale, 30.)
                      for name in ('symmetric', 'lopsided') for scale in (1., 2., 4.) for handle in ('klein', 'torus')]
    if part in ('neck', 'neck_checks'):
        births = history_births()
        if part == 'neck':
            out = [(f'neck_{name}_{clock}', .30, 5, vel, 5, None if clock == 'none' else dict(NECK, clock=clock))
                   for name, vel in births.items() for clock in ('none', 'off', 'history', 'reversed')]
            return out + [(f'neck_lopsided_torus_{clock}', .30, 5, births['lopsided'], 5,
                           dict(NECK, clock=clock, handle='torus')) for clock in ('off', 'history')]
        first = 14.                                    # through the lopsided pair's first contact
        lop = births['lopsided']
        out = [(f'check_port{g:g}_nodes{n}_{clock}', .30, 1, lop, 2, dict(port=g, nodes=n, clock=clock), first)
               for g in (25., 100., 400.) for n in (65, 129) for clock in ('off', 'history')]
        return out + [(f'check_dt_half_{clock}', .30, 1, lop, 2, dict(NECK, clock=clock), first, DT/2)
                      for clock in ('off', 'history')]
    raise SystemExit(f'unknown part {part!r}')


def antipodal_job(args):
    """One history under the antipodal-only rule (strict, or with the static field kept; the
    exchange at rest, or signed with Klein or torus reversal)."""
    tag, v, velocities, static, handle, *more = args
    exchange, scale, end = (more + ['rest', 1., HISTORY_END][len(more):])[:3]
    m = AntipodalPairBirth(Parameters(), static=static, handle=handle, exchange=exchange, flash_scale=scale)
    s0 = m.initial(v) if velocities is None else m.initial(velocities=velocities)
    t0 = time.time()
    h, info = m.run(s0, end, cycles=5, dt=DT, output_dt=OUTPUT_DT)
    keep = np.zeros(len(h['t']), bool)
    keep[::5] = True
    keep[-1] = True
    for e in info['events']:
        keep[np.abs(h['t'] - e['time']) < 1e-9] = True
    r = lambda a, k=7: [None if not np.isfinite(x) else x for x in np.round(np.asarray(a, float)[keep], k).tolist()]

    def clean(e):
        return {k: (np.asarray(x).tolist() if isinstance(x, np.ndarray) else float(x) if isinstance(x, np.floating)
                    else x) for k, x in e.items()}
    return dict(tag=tag, v=v, mode='antipodal', static=static, handle=handle, exchange=exchange, flash_scale=scale,
                retime=None, end=end,
                velocities=None if velocities is None else np.asarray(velocities).tolist(),
                birth_speed=np.linalg.norm(m.velocity(s0), axis=1).tolist(),
                t=r(h['t'], 6), d=r(h['d'], 6), x=np.round(h['x'][keep], 6).tolist(), tau=np.round(h['tau'][keep], 8).tolist(),
                alive=h['alive'][keep].tolist(), cycle=h['cycle'][keep].tolist(), clock=r(h['clock'], 9),
                energy={k: r(v_, 10) for k, v_ in h['ledger'].items()}, total=r(h['energy'], 12),
                J=np.round(h['J'][keep], 12).tolist(), n_fronts=h['n_fronts'][keep].tolist(),
                events=[clean(e) for e in info['events']], energy_drift=info['energy_drift'],
                angular_drift=info['angular_drift'], min_chronology_margin=info['min_chronology_margin'],
                length=info['length'], seconds=round(time.time()-t0, 1), contact_time=None, contact_side=None)


def history_births():
    """Birth momenta of the history runs: symmetric (v = 0.30), independent, lopsided."""
    x, vel, _, _ = birth_positions(.3, .3, PLANE)
    independent = vel.copy()
    independent[0] *= 1.15
    independent[1] = .85*vel[1] + .08*N
    _, slow, _, _ = birth_positions(.3, .15, PLANE)
    _, fast, _, _ = birth_positions(.3, .45, PLANE)
    return dict(symmetric=None, independent=independent.tolist(), lopsided=np.array([fast[0], slow[1]]).tolist())


def main(part):
    jobs = jobs_for(part)
    OUT.mkdir(parents=True, exist_ok=True)
    t0, runs = time.time(), []
    with Pool(4) as pool:
        worker = (job if part in ('scan', 'featured', 'retimed', 'checks') else
                  antipodal_job if part in ('antipodal', 'klein') else history_job)
        for r in pool.imap_unordered(worker, jobs):
            runs.append(r)
            print(f"{len(runs)}/{len(jobs)} {r['tag']} v={r['v']:.3f} {r['mode']} retime={r['retime']} "
                  f"contact={r['contact_time']} {r['contact_side']} drift={r['energy_drift']:.1e} "
                  f"{time.time()-t0:.0f}s", flush=True)
    runs.sort(key=lambda r: (r['tag'], r['mode'], r['v'], str(r['retime'])))
    meta = dict(part=part, END=END, DT=DT, output_dt=OUTPUT_DT, plane=PLANE, numpy=np.__version__)
    (OUT/f'{part}.json').write_text(json.dumps(dict(meta=meta, runs=runs), separators=(',', ':'))+'\n')


def elliptic_part():
    """The strict history under (n, t, q) -> (-n, t_c - t, -q): a speed scan of the first flight,
    resonant against off-resonant cycles, and the momentum the static field moves per unit energy."""
    from scipy.optimize import brentq
    m = AntipodalPairBirth(Parameters())
    a = m.par.mouth_angle

    def birth_speed(v):
        return float(np.linalg.norm(m.velocity(m.initial(v)), axis=1)[0])

    def kick_for(u):
        return brentq(lambda v: birth_speed(v) - u, .02, .95, xtol=1e-14)

    def compact(h, keep_every=5):
        idx = np.unique(np.r_[np.arange(0, len(h['t']), keep_every), len(h['t']) - 1])
        return dict(t=np.round(h['t'][idx], 6).tolist(), x=np.round(h['x'][idx], 7).tolist(),
                    alive=h['alive'][idx].tolist())

    resonant = resonant_speeds(a, 4)
    speeds = sorted(set(np.round(np.arange(.10, .845, .01), 4).tolist() + (resonant*(1 - 1e-4)).tolist()))
    scan = []
    for u in speeds:
        v = kick_for(u)
        s = m.initial(v)
        h, info = m.run(s, (np.pi - 2*a)/u + .02, dt=DT, output_dt=OUTPUT_DT)
        trajectory, flash, contact, captures = elliptic_mismatch(h, info)
        scan.append(dict(kick=v, speed=u, contact=contact, trajectory_mismatch=trajectory, flash_mismatch=flash,
                         captures=[dict(time=c['time'], front=c['front'], mouth=c['mouth']) for c in captures],
                         energy_drift=info['energy_drift']))
        print(f'u={u:.5f} t_c={contact:.4f} trajectory {trajectory:.1e} flash {flash:.4f} '
              f'captures {[(round(c["time"], 4), c["front"]) for c in captures]}', flush=True)
    cycles = []
    for label, v in [('resonant_k1', kick_for(resonant[1]*(1 - 2e-4))),
                     ('resonant_k2', kick_for(resonant[2]*(1 - 2e-4))), ('off_resonance_v030', .30)]:
        s = m.initial(v)
        end = 4*((np.pi - 2*a)/birth_speed(v) + np.pi) - .5
        h, info = m.run(s, end, cycles=4, dt=DT, output_dt=OUTPUT_DT)
        ev = info['events']
        births = [dict(time=0., speed=birth_speed(v), flash=float(sum(f.energy for f in s.fronts)))] + \
            [dict(time=e['time'], speed=e['speed'][0], flash=float(sum(e['flash'])))
             for e in ev if e['kind'] == 'birth']
        contacts = [e['time'] for e in ev if e['kind'] == 'contact']
        n_fronts = np.asarray(h['n_fronts'])
        left = [int(n_fronts[np.flatnonzero(np.abs(np.asarray(h['t']) - tc) < 1e-12)[-1]]) for tc in contacts]
        cycles.append(dict(label=label, kick=v, births=births, contacts=contacts, fronts_after_contact=left,
                           captures=[dict(time=e['time'], front=e['front'], mouth=e['mouth']) for e in ev
                                     if e['kind'] == 'capture'],
                           energy_drift=info['energy_drift'], angular_drift=info['angular_drift'], history=compact(h)))
        print(label, [round(b['speed'], 6) for b in births], 'fronts left', left, flush=True)
    # The static field: momentum moved per unit energy exchanged, first flight of a v = 0.30 kick.
    ms = AntipodalPairBirth(Parameters(), static='coulomb')
    s = ms.initial(.30)
    h, info = ms.run(s, 12., dt=DT, output_dt=OUTPUT_DT)
    contact = next(e for e in info['events'] if e['kind'] == 'contact')
    alive = np.asarray(h['alive'])
    L = np.asarray(h['L'])[alive]
    p_birth, p_contact = np.linalg.norm(L[0], axis=1), np.linalg.norm(L[-1], axis=1)
    kinetic = np.asarray(h['ledger']['kinetic'])[alive]
    v_birth = ms.velocity(s)
    turned = bool(np.all(np.sum(np.cross(L[-1], np.asarray(h['x'])[alive][-1])*v_birth, axis=1) < 0))
    flash = float(sum(f.energy for f in s.fronts))/2
    static = dict(kick=.30, contact=contact['time'], momentum_birth=p_birth.tolist(),
                  momentum_contact=p_contact.tolist(),
                  reversed=turned, momentum_change_per_mouth=float(np.mean(p_birth + p_contact)),
                  energy_exchanged_per_mouth=float((kinetic.max() - kinetic.min())/2),
                  flash_energy_per_mouth=flash, birth_speed=float(np.linalg.norm(v_birth[0])))
    print('static field', static, flush=True)
    meta = dict(part='elliptic', DT=DT, output_dt=OUTPUT_DT, resonant_speeds=resonant.tolist(), numpy=np.__version__)
    (OUT/'elliptic.json').write_text(json.dumps(dict(meta=meta, scan=scan, cycles=cycles, static=static),
                                                separators=(',', ':'))+'\n')


def s3_part():
    """S3 against S2: exact refocusing, closure of a flash released at rest, what a finite absorbing
    mouth takes from a flash, and the strict router with the 3-D focus phase and neck."""
    from scipy.optimize import brentq
    from geometrodynamics.waves.antipodal_pair_birth import closure_band
    from geometrodynamics.waves.s3_focusing import (absorbed_by_mouth, closure_correction, midpoint_grid,
                                                    neck_length_s3, s3_zonal_evolve, zonal_energy, zonal_fd,
                                                    zonal_modes_s2, zonal_modes_s3)
    a = Parameters().mouth_angle
    refocus, profiles = [], None
    for n in (200, 400, 800):
        psi = midpoint_grid(n)
        u0, p0 = np.exp(-(psi/.25)**2), np.zeros(n)
        s = np.sin(psi)
        ue, pe = s3_zonal_evolve(u0, p0, np.pi)
        uf, pf = zonal_fd(u0, p0, np.pi, sphere=3)
        u2, p2 = zonal_fd(u0, p0, np.pi, sphere=2)
        norm3, norm2 = np.sqrt(np.sum(u0**2*s**2)), np.sqrt(np.sum(u0**2*s))
        refocus.append(dict(cells=n, exact_mirror_error=float(np.sqrt(np.sum((ue + u0[::-1])**2*s**2))/norm3),
                            fd_mirror_error=float(np.sqrt(np.sum((uf + u0[::-1])**2*s**2))/norm3),
                            fd_energy_drift=abs(zonal_energy(uf, pf) - zonal_energy(u0, p0))/zonal_energy(u0, p0),
                            s2_displacement=float(np.sqrt(np.sum(u2**2*s))/norm2),
                            s2_mirror_distance=float(np.sqrt(np.sum((u2 + u0[::-1])**2*s))/norm2),
                            s2_energy_drift=abs(zonal_energy(u2, p2, 2) - zonal_energy(u0, p0, 2))/zonal_energy(u0, p0, 2)))
        if n == 400:
            profiles = dict(psi=psi.tolist(), u0=u0.tolist(), s3_exact=ue.tolist(), s3_fd=uf.tolist(),
                            s2_u=u2.tolist(), s2_p=p2.tolist())
        print('refocus', refocus[-1], flush=True)
    psi = midpoint_grid(400)
    u0, p0 = np.exp(-(psi/.25)**2), np.zeros(400)
    m3, q3, v3 = zonal_modes_s3(u0, p0)
    l2, q2, v2 = zonal_modes_s2(u0, p0, lmax=150)
    closure = [dict(k=k, s3=closure_correction(m3, (-1.)**m3, q3, v3, (2*k + 1)*np.pi),
                    s2_conformal=closure_correction(l2 + .5, (-1.)**(l2 + 1), q2, v2, (2*k + 1)*np.pi))
               for k in range(4)]
    print('closure', closure, flush=True)
    deltas = np.round(np.linspace(0, 1, 41), 4)
    absorption = [dict(sphere=sphere, speed=speed, offset='across' if across else 'along', delta=float(d),
                       **absorbed_by_mouth(a, float(d), speed, across=across, sphere=sphere, nodes=160))
                  for sphere in (3, 2) for speed in (0., .27) for across in (False, True) for d in deltas]
    L2, L3 = float(np.sin(2*a)), neck_length_s3(a)
    bands = [dict(k=k, s2=closure_band(a, k), s3=closure_band(a, k, length=L3)) for k in range(4)]
    m = AntipodalPairBirth(Parameters(), reflection_phase=-1, neck_length=L3)

    def kick_for(u):
        return brentq(lambda v: np.linalg.norm(m.velocity(m.initial(v)), axis=1)[0] - u, .02, .95, xtol=1e-14)

    router = []
    for k in (1, 2):
        low, uk = closure_band(a, k, length=L3)
        for label, u in (('inside', low*(1 + 3e-3)), ('below', low*(1 - 3e-3)), ('at_u_k', uk)):
            h, info = m.run(m.initial(kick_for(u)), 3*((np.pi - 2*a)/u + np.pi) - .5, cycles=3, dt=DT,
                            output_dt=OUTPUT_DT)
            births = [e['speed'][0] for e in info['events'] if e['kind'] == 'birth']
            router.append(dict(k=k, case=label, speed=u, speed_drift=float(max(abs(b - u) for b in births)),
                               fronts_left=len(info['final_state'].fronts), energy_drift=info['energy_drift']))
            print('router', router[-1], flush=True)
    couplings = {}
    for label, phase, length in (('S3: -1 per focus, 3-D neck', -1., L3), ('S2: -i per focus', -1j, None)):
        ms = AntipodalPairBirth(Parameters(), reflection_phase=phase, neck_length=length, exchange='signed')
        h, info = ms.run(ms.initial(.30), HISTORY_END, cycles=5, dt=DT, output_dt=OUTPUT_DT)
        caps = [e['coupling'] for e in info['events'] if e['kind'] == 'capture']
        couplings[label] = dict(push=sum(c > 0 for c in caps), pull=sum(c < 0 for c in caps),
                                quadrature=sum(c == 0 for c in caps), energy_drift=info['energy_drift'],
                                angular_drift=info['angular_drift'])
    print('couplings', couplings, flush=True)
    tally = []
    for label, phase, length in (('S3', -1., L3), ('S2', -1j, None)):
        for handle in ('klein', 'torus'):
            for v in (.30, 1/3):
                ms = AntipodalPairBirth(Parameters(), reflection_phase=phase, neck_length=length, exchange='signed',
                                        handle=handle)
                s = ms.initial(v)
                h, info = ms.run(s, HISTORY_END, cycles=5, dt=DT, output_dt=OUTPUT_DT)
                ev = info['events']
                tally.append(dict(space=label, handle=handle, kick=v,
                                  birth_speed=float(np.linalg.norm(ms.velocity(s), axis=1)[0]),
                                  capture_drate=float(sum(e['drate'] for e in ev if e['kind'] == 'capture')),
                                  exit_drate=float(sum(e['drate'] for e in ev if e['kind'] == 'throat_exit')),
                                  exits=sum(e['kind'] == 'throat_exit' for e in ev),
                                  contacts=[e['time'] for e in ev if e['kind'] == 'contact'],
                                  energy_drift=info['energy_drift'], angular_drift=info['angular_drift']))
                print('tally', {k: tally[-1][k] for k in ('space', 'handle', 'kick', 'capture_drate', 'exit_drate')},
                      flush=True)
    meta = dict(part='s3', mouth_angle=a, neck_s2=L2, neck_s3=L3, numpy=np.__version__,
                flash='the S2 Maxwell flash energies of PairBirth; only the routing and neck are S3')
    (OUT/'s3.json').write_text(json.dumps(dict(meta=meta, refocus=refocus, profiles=profiles, closure=closure,
                                               absorption=absorption, bands=bands, router=router,
                                               couplings=couplings, tally=tally), separators=(',', ':'))+'\n')


if __name__ == '__main__':
    part = sys.argv[1] if len(sys.argv) > 1 else 'scan'
    {'elliptic': elliptic_part, 's3': s3_part}.get(part, lambda: main(part))()
