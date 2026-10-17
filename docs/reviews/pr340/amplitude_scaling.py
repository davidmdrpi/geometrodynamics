"""Weak-field scaling of the single-packet nonadditivity residual (PR #340 follow-up review).

Scales the PR's prepared field (both packets together, and each alone with the joint
normalization) by lambda, i.e. energy 0.1*lambda^2, and divides by lambda^2.
At fixed mouth positions the torque is quadratic in the field, so a residual that is
mainly the bilinear packet-packet term is constant after division by lambda^2, while
trajectory/backreaction corrections vanish as lambda -> 0.

Run from a checkout of PR #340 at 650e825:  PYTHONPATH=. python amplitude_scaling.py
"""
import json

import numpy as np

from geometrodynamics.waves.moving_mouth_pairs import MovingMouthPairs, harmonics

m = MovingMouthPairs()


def states(eta, lam):
    prof = [m.ell*(m.ell+1)*np.exp(-m.par.sigma**2*m.ell*(m.ell+1)/2)*harmonics(m.ell, m.order, p).conj()
            for p in [[1., 0., 0.], [0., 0., 1.]]]
    f = [prof[0], .65*prof[1]]
    v = [-eta*(m.generators[2]@f[0]), -eta*(m.generators[0]@f[1])]
    scale = lam*np.sqrt(m.par.energy/m.wave_energy(sum(f), sum(v)))
    out = {}
    for name, w in [('both', (1, 1)), ('x_only', (1, 0)), ('z_only', (0, 1))]:
        s = m.initial(eta)
        fi = w[0]*f[0]+w[1]*f[1]
        vi = w[0]*v[0]+w[1]*v[1]
        s.a = -(-1.)**m.ell*vi*scale/m.omega
        s.p = (-1.)**m.ell*m.omega*fi*scale
        out[name] = s
    return out


def angle(a, b):
    return np.degrees(np.arccos(np.clip(np.sum(a*b, -1), -1, 1)))


res = {}
for eta, tag in [(.7, 'directional'), (0., 'symmetric')]:
    for lam in (1., .5, .25):
        runs = {k: m.evolve(s)[0] for k, s in states(eta, lam).items()}
        L = {k: r['angular'][-1] for k, r in runs.items()}
        resid = L['both']-L['x_only']-L['z_only']
        shared = angle(runs['both']['x'][-1, 0], runs['both']['x'][-1, 1])-90
        composite = angle(runs['z_only']['x'][-1, 0], runs['x_only']['x'][-1, 1])-90
        r = dict(lam=lam, resid_over_lam2=(resid/lam**2).tolist(),
                 diag_over_lam2=[float(L['z_only'][0, 0]/lam**2), float(L['x_only'][1, 2]/lam**2)],
                 shared_sep_over_lam2=shared/lam**2, composite_sep_over_lam2=composite/lam**2,
                 difference_over_lam2=(shared-composite)/lam**2, shared_sep=shared, composite_sep=composite)
        res[f'{tag}_{lam}'] = r
        print(tag, json.dumps(r), flush=True)
json.dump(res, open('amplitude_scaling.json', 'w'), indent=1)
