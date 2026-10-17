#!/usr/bin/env python
"""Run the frozen two-pair comparison, recording both mouth and field budgets."""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from geometrodynamics.waves.antipodal_kick import Parameters
from geometrodynamics.waves.moving_mouth_pairs import MovingMouthPairs


def save_trace(path, data):
    names, columns = [], []
    for key, values in data.items():
        if key == 'a':
            continue
        if values.ndim == 1:
            names.append(key)
            columns.append(values)
        else:
            for idx in np.ndindex(values.shape[1:]):
                names.append(key+'_'+'_'.join(map(str, idx)))
                columns.append(values[(slice(None),)+idx])
    np.savetxt(path, np.column_stack(columns), delimiter=',', header=','.join(names),
               comments='', fmt='%.12e')


def run(output):
    output.mkdir(parents=True, exist_ok=True)
    par = Parameters(mass=.10)
    protocol = Path(__file__).resolve().parents[2]/'docs/moving_mouth_pairs_protocol.md'
    summary = dict(parameters=asdict(par), dt=.002, horizon=float(np.pi+1.2),
                   protocol_sha256=hashlib.sha256(protocol.read_bytes()).hexdigest(),
                   runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                   cases={}, convergence={})
    saved, render = {}, {}
    cases = [('directional', .7, 1., False, 32, .002),
             ('symmetric', 0., 1., False, 32, .002),
             ('mirror', -.7, 1., False, 32, .002),
             ('uncoupled', .7, 0., False, 32, .002),
             ('empty', .7, 1., True, 32, .002),
             ('L40', .7, 1., False, 40, .002),
             ('dt_half', .7, 1., False, 32, .001)]
    for name, eta, coupling, empty, cutoff, dt in cases:
        model = MovingMouthPairs(replace(par, cutoff=cutoff, coupling=coupling))
        data, stats = model.evolve(model.initial(eta, empty=empty), dt=dt)
        summary['convergence' if name in ['L40', 'dt_half'] else 'cases'][name] = stats
        save_trace(output/f'{name}.csv', data)
        saved[name] = {k: data[k] for k in ['angular', 'x', 'q', 'mouth_impulse']}
        if name in ['directional', 'symmetric']:
            ids = np.unique(np.r_[np.arange(0, len(data['t']), 2),
                                 np.argmin(abs(data['t']-np.pi)), len(data['t'])-1])
            for key in ['t', 'a', 'x', 'angular', 'q']:
                render[f'{name}_{key}'] = data[key][ids]
        print(name, json.dumps(stats), flush=True)
    np.savez_compressed(output/'render_states.npz', **render)
    base = saved['directional']['angular'][-1]
    changes = {name: float(np.linalg.norm(saved[name]['angular'][-1]-base)) for name in ['L40', 'dt_half']}
    stats = list(summary['cases'].values())+list(summary['convergence'].values())
    mirror_errors = {
        'position': float(np.max(abs(saved['mirror']['x']-saved['directional']['x']*[1, -1, 1]))),
        'angular': float(np.max(abs(saved['mirror']['angular']-saved['directional']['angular']*[-1, 1, -1]))),
        'shape': float(np.max(abs(saved['mirror']['q']-saved['directional']['q']))),
        'mouth_impulse': float(np.max(abs(saved['mirror']['mouth_impulse']-saved['directional']['mouth_impulse']*[-1, 1, -1]))),
    }
    gates = dict(energy=all(s['max_energy_relative'] < 1e-4 for s in stats),
                 angular=all(s['max_angular_relative'] < 1e-8 for s in stats),
                 constraints=all(max(s['max_sphere_error'], s['max_tangency_error']) < 1e-10 for s in stats),
                 space_convergence=changes['L40'] < max(.02*np.linalg.norm(base), 1e-7),
                 time_convergence=changes['dt_half'] < max(.02*np.linalg.norm(base), 1e-7),
                 mirror=max(mirror_errors.values()) < 1e-8,
                 null_controls=all(summary['cases'][n]['max_pair_angular'] < 1e-13 for n in ['empty', 'uncoupled']))
    numerics = all(gates.values())
    gates.update(both_pairs_move=bool(np.all(np.linalg.norm(base, axis=1) > 1e-6)),
                 relative_motion=summary['cases']['directional']['max_separation_change'] > 1e-5,
                 disjoint_mouths=all(s['minimum_cross_mouth_distance'] > .65 for s in stats))
    summary.update(gates={k: bool(v) for k, v in gates.items()},
                   final_impulse_difference=changes, mirror_max_errors=mirror_errors,
                   verdict=('UNRESOLVED_NUMERICALLY' if not numerics else
                            'SUPPORTED_FOR_THIS_EFFECTIVE_FAMILY' if all(gates.values()) else 'NOT_SUPPORTED'))
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    files = sorted(output.glob('*.csv'))+[output/'render_states.npz', output/'summary.json']
    (output/'sha256.json').write_text(json.dumps({p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                                for p in files}, indent=2)+'\n')
    print(summary['verdict'], summary['gates'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('experiments/closure_ledger/runs/20261016_moving_mouth_pairs'))
    run(parser.parse_args().output)
