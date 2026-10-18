"""Cross-check of start_offset.py with the PR's moving model: directional, delta = 1.2, lambda = 1/64.

Run with a checkout of PR #341 on PYTHONPATH.
"""
import numpy as np

from geometrodynamics.waves.moving_mouth_pairs import MovingMouthPairs, harmonics

m = MovingMouthPairs()
eta, lam, delta = .7, 1/64, 1.2
prof = [m.ell*(m.ell+1)*np.exp(-m.par.sigma**2*m.ell*(m.ell+1)/2)*harmonics(m.ell, m.order, p).conj()
        for p in [[1., 0., 0.], [0., 0., 1.]]]
f = [prof[0], .65*prof[1]]
v = [-eta*(m.generators[2]@f[0]), -eta*(m.generators[0]@f[1])]
scale = lam*np.sqrt(m.par.energy/m.wave_energy(sum(f), sum(v)))
L = []
for wgt in [(1, 1), (1, 0), (0, 1)]:
    s = m.initial(eta)
    a = -(-1.)**m.ell*(wgt[0]*v[0]+wgt[1]*v[1])*scale/m.omega
    p = (-1.)**m.ell*m.omega*(wgt[0]*f[0]+wgt[1]*f[1])*scale
    co, si = np.cos(m.omega*delta), np.sin(m.omega*delta)
    s.a, s.p = a*co-p*si/m.omega, p*co+m.omega*a*si
    arr, st = m.evolve(s, end=np.pi+delta+1.2)
    L.append(arr['angular'][-1])
r = (L[0]-L[1]-L[2])/lam**2
print('moving model, delta=1.2, lambda=1/64: C1y', r[0, 1], 'C2y', r[1, 1])
