"""Independent controls for PR #340 using the PR's MovingMouthPairs unchanged.

Run from a checkout of PR #340:  PYTHONPATH=. python review_controls.py
Writes controls.json and long_{directional,symmetric}.npz (used by the other scripts).

1. single-packet decomposition (x-packet only, z-packet only, both) at the
   combined state's normalization, so the parts superpose linearly in the field;
2. long horizon 4*pi + 1.2 with exact checkpoints at every k*pi +- 0.6;
3. drawing-frame twist: integrated body rotation vs the renderer's shortest rotation.
"""
import json
import sys

import numpy as np

from geometrodynamics.waves.moving_mouth_pairs import MovingMouthPairs, PairState, harmonics

m = MovingMouthPairs()


def initial(eta, wx=1., wz=.65):
    prof = [m.ell*(m.ell+1)*np.exp(-m.par.sigma**2*m.ell*(m.ell+1)/2)*harmonics(m.ell, m.order, x).conj()
            for x in [[1., 0., 0.], [0., 0., 1.]]]
    # normalization of the full two-packet launch, as in MovingMouthPairs.initial
    f_full = prof[0]+.65*prof[1]
    v_full = -eta*(m.generators[2]@prof[0]+.65*m.generators[0]@prof[1])
    scale = np.sqrt(m.par.energy/m.wave_energy(f_full, v_full))
    f = wx*prof[0]+wz*prof[1]
    v = -eta*(wx*m.generators[2]@prof[0]+wz*m.generators[0]@prof[1])
    par = (-1.)**m.ell
    return PairState(-par*v*scale/m.omega, par*m.omega*f*scale, np.array([[0., 0., 1.], [1., 0., 0.]]),
                     np.zeros((2, 3)), np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2, 3)))


def run(s, end, dt=.002, checkpoints=()):
    s = s.copy()
    times = np.unique(np.r_[np.arange(0, end, .02), [c for c in checkpoints if c <= end], end])
    R = np.repeat(np.eye(3)[None], 2, axis=0)
    rec = dict(t=[0.], x=[s.x.copy()], L=[s.angular.copy()], imp=[s.mouth_impulse.copy()],
               J=[m.wave_angular(s)], R=[R.copy()], E=[m.energy(s)])
    j0 = m.wave_angular(s)+s.angular.sum(0)
    t = 0.
    for target in times[1:]:
        n = max(1, int(np.ceil((target-t)/dt-1e-10)))
        h = (target-t)/n
        for _ in range(n):
            # body rotation integrated with the same per-pair geodesic flow the drift uses
            m.kick(s, h/2)
            for k in range(2):
                ang = s.angular[k]
                ln = np.linalg.norm(ang)
                if ln:
                    ax = ang/ln
                    th = ln*h/m.par.mass
                    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
                    R[k] = (np.eye(3)+np.sin(th)*K+(1-np.cos(th))*K@K)@R[k]
            m.drift(s, h)
            m.kick(s, h/2)
        t = target
        for k, v in [('t', t), ('x', s.x.copy()), ('L', s.angular.copy()), ('imp', s.mouth_impulse.copy()),
                     ('J', m.wave_angular(s)), ('R', R.copy()), ('E', m.energy(s))]:
            rec[k].append(v)
    out = {k: np.asarray(v) for k, v in rec.items()}
    # normalized as in the PR: by max(initial energy, |J0|)
    out['Jerr'] = np.abs(out['J']+out['L'].sum(1)-j0).max()/max(out['E'][0], np.linalg.norm(j0), 1e-12)
    out['Eerr'] = np.abs(out['E']-out['E'][0]).max()/out['E'][0]
    return out


def deg(x0, x):
    return np.degrees(np.arccos(np.clip(np.sum(x0*x, -1), -1, 1)))


def shortest(start, end):
    v = np.cross(start, end)
    c = start@end
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3)+K+K@K/(1+c)


def twist(o):
    """Angle about x(t) between integrated body frame and the shortest-rotation frame."""
    out = []
    for k in range(2):
        x0 = o['x'][0, k]
        e = np.cross(x0, [0.3, 0.5, 0.8]); e /= np.linalg.norm(e)
        tw = []
        for R, x in zip(o['R'][:, k], o['x'][:, k]):
            a, b = R@e, shortest(x0, x)@e
            tw.append(np.degrees(np.arctan2(np.cross(b, a)@x, a@b)))
        out.append(np.array(tw))
    return np.array(out)


res = {}
# 1. Single-packet decomposition over the PR horizon
end = np.pi+1.2
for eta in (.7, 0.):
    tag = 'directional' if eta else 'symmetric'
    parts = {name: run(initial(eta, wx, wz), end) for name, (wx, wz) in
             dict(both=(1, .65), x_only=(1, 0), z_only=(0, .65)).items()}
    L = {k: v['L'][-1] for k, v in parts.items()}
    cross = L['both']-L['x_only']-L['z_only']
    J0 = {k: v['J'][0] for k, v in parts.items()}
    res[tag] = dict(J0={k: v.tolist() for k, v in J0.items()}, L={k: v.tolist() for k, v in L.items()},
                    cross=cross.tolist(), Jerr=max(v['Jerr'] for v in parts.values()),
                    Eerr=max(v['Eerr'] for v in parts.values()),
                    sep_change={k: float(deg(v['x'][-1, 0], v['x'][-1, 1])-90) for k, v in parts.items()},
                    disp={k: deg(v['x'][0], v['x'][-1]).tolist() for k, v in parts.items()})
    print(tag, json.dumps(res[tag], indent=1), flush=True)

# 2. Long horizon with exact focus windows
end = 4*np.pi+1.2
cps = sorted({.6} | {k*np.pi+d for k in range(1, 5) for d in (-.6, 0, .6)})
for eta in (.7, 0.):
    tag = 'directional' if eta else 'symmetric'
    o = run(initial(eta), end, checkpoints=cps)
    np.savez(f'long_{tag}.npz', **{k: v for k, v in o.items() if k not in ('Jerr', 'Eerr')})
    idx = lambda tt: int(np.argmin(abs(o['t']-tt)))
    windows = {}
    prev = 0
    for k in range(1, 5):
        lo, hi = idx(k*np.pi-.6), idx(k*np.pi+.6)
        windows[f'{k}pi'] = dict(between=(o['imp'][lo]-o['imp'][prev]).tolist(),
                                 window=(o['imp'][hi]-o['imp'][lo]).tolist())
        prev = hi
    tw = twist(o)
    res[f'long_{tag}'] = dict(final_L=o['L'][-1].tolist(), J0=o['J'][0].tolist(),
                              disp=deg(o['x'][0], o['x'][-1]).tolist(),
                              sep_change=float(deg(o['x'][-1, 0], o['x'][-1, 1])-90),
                              Jerr=o['Jerr'], Eerr=o['Eerr'], windows=windows,
                              max_twist_deg=np.abs(tw).max(1).tolist(),
                              L_at={f'{k}pi+0.6': o['L'][idx(k*np.pi+.6)].tolist() for k in range(1, 5)})
    print(f'long_{tag}', json.dumps({k: v for k, v in res[f'long_{tag}'].items() if k != 'windows'}, indent=1), flush=True)
json.dump(res, open('controls.json', 'w'), indent=1)
