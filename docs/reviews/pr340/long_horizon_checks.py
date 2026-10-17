"""Convergence of the long-horizon symmetric run with the PR's own evolve(), plus a
displaced-pair stability test for a single symmetric packet.

Run from a checkout of PR #340:  PYTHONPATH=. python long_horizon_checks.py
"""
import json
from dataclasses import replace

import numpy as np

from geometrodynamics.waves.moving_mouth_pairs import MovingMouthPairs, Parameters

end = 4*np.pi+1.2
res = {}
for tag, par, dt in [('base', Parameters(mass=.1), .002), ('dt_half', Parameters(mass=.1), .001),
                     ('L40', Parameters(mass=.1, cutoff=40), .002)]:
    m = MovingMouthPairs(par)
    for eta in (0., .7):
        arr, st = m.evolve(m.initial(eta), dt=dt, end=end, output_dt=.04)
        sep = np.degrees(arr['separation'])
        cross = np.degrees(np.minimum(arr['separation'], np.pi-arr['separation']))
        bad = np.flatnonzero(cross < np.degrees(.65))
        res[f'{tag}_eta{eta}'] = dict(final_L=arr['angular'][-1].tolist(), Eerr=st['max_energy_relative'],
                                      Jerr=st['max_angular_relative'], sep_change_final=float(sep[-1]-90),
                                      min_cross_deg=float(cross.min()),
                                      first_violation_t=float(arr['t'][bad[0]]) if len(bad) else None,
                                      disp=st['final_displacement'])
        print(tag, eta, json.dumps(res[f'{tag}_eta{eta}']), flush=True)

# Stability of a symmetric focus for a slightly displaced pair: z-packet only, pair 1 tilted
# by delta toward +x (pair 2 left at x, far from the z-packet's foci).
m = MovingMouthPairs()
for delta in (1., 3., 6.):
    s = m.initial(0.)
    # keep only the z-packet: rebuild from the PR's preparation with x-packet weight 0
    from geometrodynamics.waves.moving_mouth_pairs import harmonics
    prof = [m.ell*(m.ell+1)*np.exp(-m.par.sigma**2*m.ell*(m.ell+1)/2)*harmonics(m.ell, m.order, x).conj()
            for x in [[1., 0., 0.], [0., 0., 1.]]]
    f_full = prof[0]+.65*prof[1]
    scale = np.sqrt(m.par.energy/m.wave_energy(f_full, 0*f_full))
    par = (-1.)**m.ell
    s.a[:] = 0
    s.p = par*m.omega*.65*prof[1]*scale
    d = np.radians(delta)
    s.x[0] = [np.sin(d), 0., np.cos(d)]
    arr, st = m.evolve(s, end=end, output_dt=.04)
    t = arr['t']
    tilt = np.degrees(np.arccos(np.clip(arr['x'][:, 0, 2], -1, 1)))
    out = {f'{k}pi+0.6': float(tilt[np.argmin(abs(t-(k*np.pi+.6)))]) for k in range(1, 5)}
    res[f'stability_{delta}'] = dict(tilt_deg=out, L1=arr['angular'][-1, 0].tolist())
    print('displaced pair', delta, json.dumps(res[f'stability_{delta}']), flush=True)
json.dump(res, open('long_checks.json', 'w'), indent=1)
