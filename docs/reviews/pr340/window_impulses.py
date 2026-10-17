"""Per-mouth impulse in each focus window of the long directional run (PR #340 review).

Reads long_directional.npz written by review_controls.py. Mouth impulses are
projected on each pair's main axis (pair 1: x, pair 2: z).
"""
import numpy as np

o = np.load('long_directional.npz')
t, x, imp = o['t'], o['x'], o['imp']
idx = lambda tt: int(np.argmin(abs(t-tt)))
names = ['P1+', 'P1-', 'P2+', 'P2-']
marks = [0, .6]+[k*np.pi+d for k in range(1, 5) for d in (-.6, .6)]+[t[-1]]
for a, b in zip(marks[:-1], marks[1:]):
    d = imp[idx(b)]-imp[idx(a)]
    row = [d[0, 0, 0], d[0, 1, 0], d[1, 0, 2], d[1, 1, 2]]
    print(f'[{a:5.2f},{b:5.2f}]', ' '.join(f'{n} {r:+.2e}' for n, r in zip(names, row)))
for k in range(1, 5):
    i = idx(k*np.pi)
    print(f'{k}pi: pair 1 +z mouth is {np.degrees(np.arccos(x[i, 0]@[0, 0, 1])):.1f} deg from +z; '
          f'pair 2 +x mouth is {np.degrees(np.arccos(x[i, 1]@[1, 0, 0])):.1f} deg from +x')
