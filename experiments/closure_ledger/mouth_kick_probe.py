#!/usr/bin/env python
"""Single-exchange kicks on two antipodal mouths: scan, checks, throat flux, clock offsets.

    python experiments/closure_ledger/mouth_kick_probe.py scan      # handle x launch x entry side
    python experiments/closure_ledger/mouth_kick_probe.py checks    # conservation, mesh, controls
    python experiments/closure_ledger/mouth_kick_probe.py flux      # static throat flux x wave cross term
    python experiments/closure_ledger/mouth_kick_probe.py offsets   # mouth clock offsets (MTY)

Each part writes JSON to experiments/closure_ledger/runs/20261019_mouth_kicks/.
Model and bookkeeping: geometrodynamics/waves/rim_kicks.py. Write-up:
docs/mouth_kick_experiments.md. Runs use four worker processes.
"""
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from geometrodynamics.waves.rim_kicks import (FluxField, GluedWave, glued_surface, offset_history,
                                              release_source)

OUT = Path(__file__).resolve().parent/'runs'/'20261019_mouth_kicks'
T_END = 2*np.pi+.5
D0 = 1.0
SAVE_DT = .04
_surfaces = {}


def surface(handle, n_sphere):
    key = (handle, n_sphere)
    if key not in _surfaces:
        _surfaces[key] = glued_surface(handle, n_sphere=n_sphere)
    return _surfaces[key]


def thin(h, extra=None):
    """Keep every SAVE_DT of a history, as plain lists."""
    step = max(1, int(round(SAVE_DT/(h.t[1]-h.t[0]))))
    keep = slice(0, None, step)
    out = dict(t=np.round(h.t[keep], 6).tolist(), J=np.round(h.J[keep], 10).tolist(), L0=h.L[0].tolist())
    if h.Jx is not None:
        out.update(Jx=np.round(h.Jx[keep], 10).tolist(), Lx0=h.Lx[0].tolist())
    if h.energy is not None:
        out['energy_drift'] = float(np.abs(h.energy-h.energy[0]).max()/h.energy[0])
    out.update(extra or {})
    return out


def direct(args):
    """One direct run: (tag, handle, mode, psi, n_sphere, width, d0, sign, m2, weight, flux)."""
    tag, handle, mode, psi, n_sphere, width, d0, sign, m2, weight, flux = args
    s = surface(handle, n_sphere)
    g = GluedWave(s, m2=m2, weight=weight)
    if mode == 'ring':
        u, v = g.ring(d0, width)
    else:
        u, v, _ = g.packet(d0, psi, mode, width)
    static = None
    if flux:
        F = FluxField.unit_energy()
        static = F.band_gradient(g.kicks)/np.sqrt(F.energy)
    h = g.run(sign*u, sign*v, T_END, static=static)
    e1, e2 = g.e1, g.e2
    return thin(h, dict(tag=tag, handle=handle, mode=mode, psi=psi, n_sphere=n_sphere, width=width, d0=d0, sign=sign,
                        m2=m2, weight=list(weight), flux=flux, e1=e1.tolist(), e2=e2.tolist()))


def offset_run(args):
    handle, psi, offset, n_sphere = args
    s = surface(handle, n_sphere)
    h, info = offset_history(s, release_source(psi, D0, s), (offset,) if offset else (), t_max=T_END,
                             flux=FluxField.unit_energy())
    from geometrodynamics.waves.rim_kicks import mouth_frame
    e1, e2 = mouth_frame(s)
    return thin(h, dict(tag='offset', handle=handle, mode='release', psi=psi, offset=offset, n_sphere=n_sphere,
                        e1=e1.tolist(), e2=e2.tolist(), **info))


def jobs_for(part):
    w = (.35, .9)
    psis = [k*np.pi/4 for k in range(8)]
    if part == 'scan':
        return direct, [('scan', h, m, p, 24000, .25, D0, 1, .25, w, False)
                        for h in ('klein', 'torus') for m in ('release', 'out', 'in') for p in psis]
    if part == 'checks':
        j = [('bare', 'none', 'out', .7, 24000, .25, D0, 1, .25, w, False),
             ('ring', 'klein', 'ring', 0., 24000, .25, .8, 1, .25, w, False),
             ('wide_weight', 'klein', 'release', np.pi/4, 24000, .25, D0, 1, .25, (.5, 1.3), False),
             ('sign_flip', 'klein', 'release', np.pi/2, 24000, .25, D0, -1, .25, w, False)]
        j += [('mesh48k', h, m, p, 48000, .25, D0, 1, .25, w, False)
              for h in ('klein', 'torus') for m in ('release', 'in') for p in (0., np.pi/2)]
        j += [('wide_packet', h, 'in', p, 24000, .35, 1.2, 1, .25, w, False) for h in ('klein', 'torus') for p in (0., np.pi/2)]
        j += [('massless', h, m, np.pi/2, 24000, .25, D0, 1, 0., w, False) for h in ('klein', 'torus') for m in ('release', 'in')]
        return direct, j
    if part == 'flux':
        return direct, [('flux', h, m, p, 24000, .25, D0, 1, 0., w, True)
                        for h in ('klein', 'torus') for m in ('release', 'out', 'in') for p in (0., np.pi/4, np.pi/2)]
    if part == 'offsets':
        return offset_run, [(h, p, o, 12000) for h in ('klein', 'torus') for p in (0., np.pi/2)
                            for o in (-2., -1.5, -1., -.5, 0., .5, 1.)]
    raise SystemExit(f'unknown part {part!r}')


def main(part):
    fn, jobs = jobs_for(part)
    OUT.mkdir(parents=True, exist_ok=True)
    t0, runs = time.time(), []
    with Pool(4) as pool:
        for r in pool.imap_unordered(fn, jobs):
            runs.append(r)
            print(f"{len(runs)}/{len(jobs)} {r['tag']} {r['handle']} {r['mode']} psi={r['psi']:.3f} "
                  f"{'offset=' + str(r.get('offset')) if 'offset' in r else ''} {time.time()-t0:.0f}s", flush=True)
    runs.sort(key=lambda r: (r['tag'], r['handle'], r['mode'], r['psi'], r.get('offset', 0), r.get('n_sphere', 0)))
    meta = dict(part=part, T_END=T_END, D0=D0, save_dt=SAVE_DT, mouth_angle=.3, numpy=np.__version__)
    (OUT/f'{part}.json').write_text(json.dumps(dict(meta=meta, runs=runs), separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'scan')
