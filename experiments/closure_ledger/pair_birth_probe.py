#!/usr/bin/env python
"""Pair birth to annihilation for two Wheeler-flux mouths: speed scan, re-timed flashes, checks.

    python experiments/closure_ledger/pair_birth_probe.py scan       # speed x {full, dressed, instant}
    python experiments/closure_ledger/pair_birth_probe.py featured   # v = 0.30, 1/3 full runs with field maps
    python experiments/closure_ledger/pair_birth_probe.py retimed    # same birth, flash removed or re-timed at t = 3
    python experiments/closure_ledger/pair_birth_probe.py checks     # step, cutoff, free motion, independent mouths

Each part writes JSON to experiments/closure_ledger/runs/20261019_pair_birth/.
Model: geometrodynamics/waves/pair_birth.py. Write-up: docs/pair_birth_timeline.md.
Runs use four worker processes.
"""
import json
import sys
import time
from dataclasses import replace
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from geometrodynamics.waves.pair_birth import (PairBirth, Parameters, birth_positions, circle_angle,
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
    raise SystemExit(f'unknown part {part!r}')


def main(part):
    jobs = jobs_for(part)
    OUT.mkdir(parents=True, exist_ok=True)
    t0, runs = time.time(), []
    with Pool(4) as pool:
        for r in pool.imap_unordered(job, jobs):
            runs.append(r)
            print(f"{len(runs)}/{len(jobs)} {r['tag']} v={r['v']:.3f} {r['mode']} retime={r['retime']} "
                  f"contact={r['contact_time']} {r['contact_side']} drift={r['energy_drift']:.1e} "
                  f"{time.time()-t0:.0f}s", flush=True)
    runs.sort(key=lambda r: (r['tag'], r['mode'], r['v'], str(r['retime'])))
    meta = dict(part=part, END=END, DT=DT, output_dt=OUTPUT_DT, plane=PLANE, numpy=np.__version__)
    (OUT/f'{part}.json').write_text(json.dumps(dict(meta=meta, runs=runs), separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'scan')
