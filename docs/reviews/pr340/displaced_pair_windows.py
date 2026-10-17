"""Displaced-pair test (PR #340 review): one symmetric z-packet, pair 1 tilted 1 degree
toward +x and at rest. Prints pair 1's L_y gained in each focus window k*pi +- 0.6.

Run from a checkout of PR #340:  PYTHONPATH=. python displaced_pair_windows.py
"""
import numpy as np

from geometrodynamics.waves.moving_mouth_pairs import MovingMouthPairs, harmonics

m = MovingMouthPairs()
prof = [m.ell*(m.ell+1)*np.exp(-m.par.sigma**2*m.ell*(m.ell+1)/2)*harmonics(m.ell, m.order, x).conj()
        for x in [[1., 0., 0.], [0., 0., 1.]]]
scale = np.sqrt(m.par.energy/m.wave_energy(prof[0]+.65*prof[1], 0*prof[0]))
s = m.initial(0.)
s.a[:] = 0
s.p = (-1.)**m.ell*m.omega*.65*prof[1]*scale
d = np.radians(1.)
s.x[0] = [np.sin(d), 0, np.cos(d)]
end = 4*np.pi+1.2
marks = [0]+[k*np.pi+e for k in range(1, 5) for e in (-.6, .6)]+[end]
t, prev = 0., 0.
for b in marks[1:]:
    n = int(np.ceil((b-t)/.002))
    h = (b-t)/n
    for _ in range(n):
        m.step(s, h)
    t = b
    tilt = np.degrees(np.arccos(np.clip(s.x[0, 2], -1, 1)))
    print(f't={t:6.3f} L1y={s.angular[0, 1]:+.3e} (step {s.angular[0, 1]-prev:+.2e}) tilt {tilt:5.2f} deg')
    prev = s.angular[0, 1]
