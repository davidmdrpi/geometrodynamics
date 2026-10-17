"""Direct small-amplitude residual for PR #341 (review): lambda = 1/32 and 1/64, no extrapolation.

Run with a checkout of PR #341 on PYTHONPATH.
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
    out = []
    for wgt in [(1, 1), (1, 0), (0, 1)]:
        s = m.initial(eta)
        s.a = -(-1.)**m.ell*(wgt[0]*v[0]+wgt[1]*v[1])*scale/m.omega
        s.p = (-1.)**m.ell*m.omega*(wgt[0]*f[0]+wgt[1]*f[1])*scale
        out.append(s)
    return out


res = {}
for eta, tag in [(.7, 'directional'), (0., 'symmetric')]:
    for lam in (1/32, 1/64):
        L = [m.evolve(s)[0]['angular'][-1] for s in states(eta, lam)]
        r = (L[0]-L[1]-L[2])/lam**2
        res[f'{tag}_{lam}'] = dict(C1y=float(r[0, 1]), C2y=float(r[1, 1]), full=r.tolist())
        print(tag, lam, r[0, 1], r[1, 1], flush=True)
json.dump(res, open('small_lambda.json', 'w'), indent=1)
