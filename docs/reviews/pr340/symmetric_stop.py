"""Symmetric long run of PR #340's model, stopped at the first integration step where
min(separation, pi - separation) <= 0.65 rad (the protocol's disjoint-mouth limit).

Run from a checkout of PR #340:  PYTHONPATH=. python symmetric_stop.py
"""
import json

import numpy as np

from geometrodynamics.waves.moving_mouth_pairs import MovingMouthPairs

m = MovingMouthPairs()
out = {}
for dt in (.002, .001):
    s = m.initial(0.)
    x0 = s.x.copy()
    e0 = m.energy(s)
    j0 = m.wave_angular(s)+s.angular.sum(0)
    t, prev = 0., None
    while True:
        cross = lambda st: min(np.arccos(np.clip(st.x[0]@st.x[1], -1, 1)), np.pi-np.arccos(np.clip(st.x[0]@st.x[1], -1, 1)))
        if cross(s) <= .65:
            break
        prev = (t, s.copy())
        m.step(s, dt)
        t += dt
    tl, sl = prev
    sep = np.degrees(np.arccos(np.clip(sl.x[0]@sl.x[1], -1, 1)))-90
    disp = np.degrees(np.arccos(np.clip(np.sum(sl.x*x0, 1), -1, 1)))
    out[dt] = dict(last_admissible_t=tl, first_crossing_t=t, sep_change_deg=sep, disp_deg=disp.tolist(),
                   L=sl.angular.tolist(), Eerr=abs(m.energy(sl)-e0)/e0,
                   Jerr=float(np.linalg.norm(m.wave_angular(sl)+sl.angular.sum(0)-j0)/e0))
    print(dt, json.dumps(out[dt]), flush=True)
