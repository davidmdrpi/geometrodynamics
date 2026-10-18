"""Independent fixed-axis linear response for PR #341 (review): Strang integrator, not DOP853.

Propagates the x- and z-packet components (joint normalization) through the coupled
field/shape equations with the mouth axes held fixed. Splits pair 1's mixed y torque
by mouth (+z, -z) and by term: e^z grad U^x (own-packet shape mismatch x other packet's
gradient) versus e^x grad U^z. Integrates with the trapezoid rule.

Run from this folder with a checkout of PR #341 on PYTHONPATH:
    PYTHONPATH=<pr341 checkout> python fixed_axis_mechanism.py   # writes mechanism.json, mechanism_*.npz
"""
import json

import numpy as np

from geometrodynamics.waves.moving_mouth_pairs import MovingMouthPairs, harmonics

m = MovingMouthPairs()
k, mu, Om = m.par.coupling, m.par.shape_mass, m.par.shape_frequency


def components(eta):
    prof = [m.ell*(m.ell+1)*np.exp(-m.par.sigma**2*m.ell*(m.ell+1)/2)*harmonics(m.ell, m.order, p).conj()
            for p in [[1., 0., 0.], [0., 0., 1.]]]
    f = [prof[0], .65*prof[1]]
    v = [-eta*(m.generators[2]@f[0]), -eta*(m.generators[0]@f[1])]
    scale = np.sqrt(m.par.energy/m.wave_energy(sum(f), sum(v)))
    par = (-1.)**m.ell
    return [dict(a=-par*vi*scale/m.omega, p=par*m.omega*fi*scale, q=np.zeros(4), pq=np.zeros(4)) for fi, vi in zip(f, v)]


s0 = m.initial()
w = m.interaction(s0)[0].reshape(4, -1)           # mouths: (pair1 +z, pair1 -z, pair2 +x, pair2 -x)


def fields(st):
    e = st['q']-(w.conj()@st['a']).real
    g = np.stack([(w.conj()@(d@st['a'])).real for d in m.generators], axis=-1)
    return e, g


def kick(st, h):
    e, _ = fields(st)
    st['p'] = st['p']+h*k*(e@w)
    st['pq'] = st['pq']-h*k*e


def drift(st, h):
    co, si = np.cos(m.omega*h), np.sin(m.omega*h)
    st['a'], st['p'] = st['a']*co+st['p']*si/m.omega, st['p']*co-m.omega*st['a']*si
    co, si = np.cos(Om*h), np.sin(Om*h)
    st['q'], st['pq'] = st['q']*co+st['pq']*si/(mu*Om), st['pq']*co-mu*Om*st['q']*si


def run(eta, dt):
    x, z = components(eta)
    end = np.pi+1.2
    n = int(np.ceil(end/dt))
    h = end/n
    t = np.arange(n+1)*h
    T = np.zeros((n+1, 2, 4, 3))   # term (0: e^z grad U^x, 1: e^x grad U^z), mouth, component
    for i in range(n+1):
        ex, gx = fields(x)
        ez, gz = fields(z)
        T[i, 0] = k*ez[:, None]*gx
        T[i, 1] = k*ex[:, None]*gz
        if i < n:
            for st in (x, z):
                kick(st, h/2); drift(st, h); kick(st, h/2)
    cum = np.concatenate([np.zeros((1, 2, 4, 3)), np.cumsum((T[1:]+T[:-1])/2*h, axis=0)])
    return t, T, cum



if __name__ == '__main__':
    out = {}
    for eta, tag in [(.7, 'directional'), (0., 'symmetric')]:
        res = {}
        for dt in (.002, .001):
            t, T, cum = run(eta, dt)
            pair = cum[-1].sum(axis=0)          # sum over terms -> (mouth, comp)
            res[dt] = dict(C1y=float(pair[0, 1]+pair[1, 1]), C2y=float(pair[2, 1]+pair[3, 1]))
        idx = lambda tt: int(np.argmin(abs(t-tt)))
        win = [0, .6, np.pi-.6, np.pi+.6, t[-1]]
        names = ['source half [0,.6]', 'between [.6,pi-.6]', 'target [pi-.6,pi+.6]', 'after']
        split = {}
        for term, tn in enumerate(['e^z gradU^x', 'e^x gradU^z']):
            for mouth, mn in [(0, '+z'), (1, '-z')]:
                split[f'pair1 {mn} {tn}'] = [float(cum[idx(b), term, mouth, 1]-cum[idx(a), term, mouth, 1]) for a, b in zip(win[:-1], win[1:])]
        # time of the peak |torque| for the dominant term, summed over pair 1's mouths
        tau = T[:, :, 0:2, 1].sum(axis=2)
        out[tag] = dict(C=res, windows=names, split=split,
                        peak_time={tn: float(t[np.argmax(abs(tau[:, j]))]) for j, tn in enumerate(['e^z gradU^x', 'e^x gradU^z'])})
        np.savez(f'mechanism_{tag}.npz', t=t, T=T, cum=cum)
        print(tag, json.dumps(out[tag], indent=1), flush=True)
    json.dump(out, open('mechanism.json', 'w'), indent=1)
