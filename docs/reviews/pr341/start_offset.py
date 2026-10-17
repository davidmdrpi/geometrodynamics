"""PR #341 review control: start the prepared packets delta earlier (free backpropagation by
pi + delta), with resting, unexcited mouths at the new t = 0, and run to the same horizon
relative to the packets (pi + delta + 1.2). delta = 0 is the PR's preparation, which switches
the coupling on in the middle of the source-side focus. Fixed-axis leading-order mixed y.

Run from this folder (imports fixed_axis_mechanism.py), PR #341 checkout on PYTHONPATH.
"""
import json

import numpy as np

import fixed_axis_mechanism as mech


def run(eta, delta, dt=.001):
    x, z = mech.components(eta)
    for st in (x, z):  # free backpropagation by delta (exact)
        co, si = np.cos(mech.m.omega*delta), np.sin(mech.m.omega*delta)
        st['a'], st['p'] = st['a']*co-st['p']*si/mech.m.omega, st['p']*co+mech.m.omega*st['a']*si
    end = np.pi+delta+1.2
    n = int(np.ceil(end/dt)); h = end/n
    cum = np.zeros((4, 3)); prev = None
    win = {}
    marks = {'source focus end': delta+.6, 'target window start': np.pi+delta-.6, 'target window end': np.pi+delta+.6}
    for i in range(n+1):
        ex, gx = mech.fields(x); ez, gz = mech.fields(z)
        T = mech.k*(ez[:, None]*gx+ex[:, None]*gz)
        if prev is not None:
            cum = cum+(T+prev)/2*h
        prev = T
        for name, tm in marks.items():
            if name not in win and i*h >= tm-1e-12:
                win[name] = float(cum[0, 1]+cum[1, 1])
        if i < n:
            for st in (x, z):
                mech.kick(st, h/2); mech.drift(st, h); mech.kick(st, h/2)
    return dict(C1y=float(cum[0, 1]+cum[1, 1]), C2y=float(cum[2, 1]+cum[3, 1]), cumulative_at=win)


out = {}
for eta, tag in [(.7, 'directional'), (0., 'symmetric')]:
    for delta in (0., .3, .6, .9, 1.2):
        r = run(eta, delta)
        out[f'{tag}_{delta}'] = r
        print(tag, delta, json.dumps(r), flush=True)
json.dump(out, open('start_offset.json', 'w'), indent=1)
