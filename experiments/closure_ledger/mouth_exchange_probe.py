#!/usr/bin/env python
"""Independent fixed-axis bilinear response and moving-mouth amplitude test."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from geometrodynamics.waves.antipodal_kick import Parameters
from geometrodynamics.waves.moving_mouth_pairs import MovingMouthPairs
from moving_mouth_review_probe import components

EVENTS = np.array([0., .6, np.pi-.6, np.pi, np.pi+.6, np.pi+1.2])


def fixed_reference(model, eta, rtol=2e-10, atol=2e-12):
    """Linear field/shape ODE, plus passive DIRECT mixed torque quadrature.

    No moving-mouth integrator or subtraction of finite-amplitude torques.
    Both independently propagated packets see all four fixed shape oscillators.
    """
    states = components(model, eta)
    n = len(model.ell)
    width = 2*n+8
    w = model.interaction(states[0])[0].reshape(4, n)
    k = model.par.coupling
    mu, freq = model.par.shape_mass, model.par.shape_frequency
    initial = np.concatenate([np.r_[s.a, s.p, s.q.ravel(), s.pq.ravel()] for s in states]+[np.zeros(12)])

    def unpack(z):
        a, p = z[:n], z[n:2*n]
        q, pq = z[2*n:2*n+4], z[2*n+4:width]
        e = q-(w.conj()@a).real
        grad = np.stack([(w.conj()@(d@a)).real for d in model.generators], axis=-1)
        return a, p, q, pq, e, grad

    def rhs(t, y):
        values = [unpack(y[i*width:(i+1)*width]) for i in range(2)]
        dy = []
        for a, p, q, pq, e, grad in values:
            dy.append(np.r_[p, -model.omega**2*a+k*(e@w), pq/mu, -mu*freq**2*q-k*e])
        e0, g0 = values[0][-2:]
        e1, g1 = values[1][-2:]
        dy.append((k*(e0[:, None]*g1+e1[:, None]*g0)).ravel())
        return np.concatenate(dy)

    sol = solve_ivp(rhs, (0., EVENTS[-1]), initial, method='DOP853',
                    rtol=rtol, atol=atol, t_eval=EVENTS)
    if not sol.success:
        raise RuntimeError(sol.message)
    impulses = sol.y[2*width:].real.T.reshape(-1, 2, 2, 3)
    field = []
    for y in sol.y.T:
        a0,p0,*_ = unpack(y[:width])
        a1,p1,*_ = unpack(y[width:2*width])
        field.append([-np.real(np.vdot(p0, d@a1)+np.vdot(p1, d@a0)) for d in model.generators])
    field = np.array(field)
    ledger = field-field[0]+impulses.sum(axis=(1, 2))
    return dict(t=EVENTS.tolist(), mouth_mixed_impulse=impulses.tolist(),
                pair_mixed_impulse=impulses.sum(axis=2).tolist(), mixed_field=field.tolist(),
                max_mixed_ledger=float(np.max(np.linalg.norm(ledger, axis=1))), nfev=sol.nfev)


def moving(model, state, dt):
    s = state.copy()
    e0 = model.energy(s)
    j0 = model.wave_angular(s)+s.angular.sum(axis=0)
    escale = max(e0, 1e-12)
    jscale = max(e0, np.linalg.norm(j0), 1e-12)
    stats = dict(max_energy_relative=0., max_angular_relative=0., minimum_distance=np.pi/2, stop=None)
    records = []
    def record(t):
        energy = model.energy_components(s)
        records.append(dict(t=float(t), angular=s.angular.tolist(), x=s.x.tolist(),
            mouth_impulse=s.mouth_impulse.tolist(), field=model.wave_angular(s).tolist(),
            energy={k: np.asarray(v).tolist() for k,v in energy.items()}))
    record(0.)
    t = 0.
    for target in EVENTS[1:]:
        steps = int(np.ceil((target-t)/dt-1e-10))
        h = (target-t)/steps
        for _ in range(steps):
            before = s.copy()
            model.step(s, h)
            sep = np.arccos(np.clip(s.x[0]@s.x[1], -1, 1))
            distance = min(sep, np.pi-sep)
            stats['minimum_distance'] = min(stats['minimum_distance'], float(distance))
            stats['max_energy_relative'] = max(stats['max_energy_relative'], abs(model.energy(s)-e0)/escale)
            stats['max_angular_relative'] = max(stats['max_angular_relative'],
                float(np.linalg.norm(model.wave_angular(s)+s.angular.sum(axis=0)-j0)/jscale))
            if distance <= .65:
                stats['stop'] = [float(t), float(t+h)]
                s = before
                record(t)
                return dict(records=records, stats=stats)
            t += h
        t = float(target)
        record(t)
    return dict(records=records, stats=stats)


def triple(model, eta, lam, dt):
    singles = components(model, eta)
    states = [model.initial(eta)]+singles
    runs = []
    for s in states:
        s.a *= lam
        s.p *= lam
        runs.append(moving(model, s, dt))
    if any(r['stats']['stop'] is not None for r in runs):
        return dict(runs=runs, incomplete=True)
    a = np.array([[r['angular'] for r in run['records']] for run in runs])
    f = np.array([[r['field'] for r in run['records']] for run in runs])
    x = np.array([[r['x'] for r in run['records']] for run in runs])
    mixed = a[0]-a[1]-a[2]
    mixed_f = f[0]-f[1]-f[2]
    composite = np.arccos(np.clip(np.sum(x[2,:,0]*x[1,:,1], axis=1), -1, 1))-np.pi/2
    return dict(runs=runs, incomplete=False, mixed=mixed.tolist(), mixed_field=mixed_f.tolist(),
                normalized_final=(mixed[-1]/lam**2).tolist(), composite_change=composite.tolist())


def job(eta, cutoff, dt, lam):
    model = MovingMouthPairs(Parameters(mass=.1, cutoff=cutoff))
    return fixed_reference(model, eta) if lam is None else triple(model, eta, lam, dt)


def run(out):
    out.mkdir(parents=True, exist_ok=True)
    protocol = hashlib.sha256(Path('docs/mouth_exchange_protocol.md').read_bytes()).hexdigest()
    result = json.loads((out/'results.json').read_text()) if (out/'results.json').exists() else dict(protocol_sha256=protocol, cases={})
    if result['protocol_sha256'] != protocol:
        raise ValueError('Checkpoint belongs to a different protocol')
    with ProcessPoolExecutor(max_workers=4) as pool:
        jobs = {}
        for name, eta in [('directional', .7), ('symmetric', 0.)]:
            case = result['cases'].setdefault(name, dict(reference={}, moving={}))
            for cutoff in [32, 40]:
                if str(cutoff) not in case['reference']:
                    jobs[pool.submit(job, eta, cutoff, None, None)] = (name, 'reference', str(cutoff))
                settings = [(.002, [1., .5, .25, .125]), (.001, [.25, .125])] if cutoff == 32 else [(.002, [.25, .125])]
                for dt, amplitudes in settings:
                    for lam in amplitudes:
                        key = f'L{cutoff}_dt{dt}_a{lam}'
                        if key not in case['moving']:
                            jobs[pool.submit(job, eta, cutoff, dt, lam)] = (name, 'moving', key)
        for future in as_completed(jobs):
            name, kind, key = jobs[future]
            result['cases'][name][kind][key] = future.result()
            (out/'results.json').write_text(json.dumps(result, indent=2)+'\n')
            print(name, kind, key, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('experiments/closure_ledger/runs/20261017_mouth_exchange'))
    run(parser.parse_args().output)
