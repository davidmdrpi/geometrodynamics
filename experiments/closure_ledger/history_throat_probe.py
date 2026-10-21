"""Archive the frozen history-coupled control matrix, including failures."""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path

import numpy as np

from geometrodynamics.waves.legacy_history_throat import HistoryThroat, Parameters


ROOT = Path('experiments/closure_ledger/runs/20261020_history_throat')
CASES = dict(baseline={}, asymmetric=dict(birth_a=.45, birth_b=-.15),
    no_throat=dict(coupling=0.), blocked=dict(blocked=True),
    spatial_even=dict(spatial_parity=1), reversed=dict(time_orientation=-1),
    reversed_anchor9=dict(time_orientation=-1, neck_horizon=9.),
    absorbing_load=dict(endpoint_load=1.),
    tight_tolerance=dict(rtol=2e-10, fixed_point_tolerance=2e-8),
    higher_resolution=dict(lmax=12, nodes=25, spacing=.01),
    high_no_throat=dict(lmax=12, nodes=25, spacing=.01, coupling=0.),
    high_blocked=dict(lmax=12, nodes=25, spacing=.01, blocked=True),
    high_reversed=dict(lmax=12, nodes=25, spacing=.01, time_orientation=-1),
    high_asymmetric=dict(lmax=12, nodes=25, spacing=.01, birth_a=.45, birth_b=-.15))


def archive(name, root=ROOT):
    root.mkdir(parents=True, exist_ok=True)
    parameters = replace(Parameters(), **CASES[name])
    model = HistoryThroat(parameters)
    print(f'CASE {name}', flush=True)
    result = model.solve(verbose=True)
    p = parameters; rows = []
    for t, state, energy, angular in zip(result['t'], result['states'],
                                        result['exterior_energy'], result['angular']):
        u, momentum, q, mechanical, clocks, clock_energy = model.unpack(state)
        f, fq, gap, _, ports, mu, body_energy = model.quantities(state, result['endpoint_history'])
        velocity = mechanical/body_energy
        port_force = mu/body_energy*p.coupling*model.eta*gap*(fq@u)
        coulomb, gradient = model.coulomb(q)
        field_energy = .5*(momentum@momentum+np.dot(model.omega*u, model.omega*u))
        rows.append([t, *q, *mechanical, *clocks, *clock_energy, *velocity, *mu,
            np.arccos(np.cos(q[0]-q[1])), energy, angular, field_energy, coulomb,
            *ports, *port_force, *(-gradient), *gap])
    rows = np.asarray(rows)
    columns = ('t,qa,qb,pa,pb,sa,sb,ka,kb,va,vb,mua,mub,separation,exterior_energy,'
        'angular,scalar_energy,coulomb_energy,port_a,port_b,port_force_a,port_force_b,'
        'coulomb_force_a,coulomb_force_b,gap_a,gap_b')
    np.savetxt(root/f'{name}.csv', rows, delimiter=',', header=columns, comments='', fmt='%.12g')
    np.savez_compressed(root/f'{name}.npz', t=result['t'], states=result['states'],
        sigma=result['sigma'], neck_states=result['neck_states'], neck_energy=result['neck_energy'])
    events = []
    for source in ([0, 1] if p.coupling and not p.blocked else []):
        for emission in np.arange(0., result['t'][-1], .5):
            for traversals in range(1, 7):
                destination = 1-source if traversals % 2 else source
                event = model.transfer_event(result, float(emission), source, destination,
                                              traversals*p.length)
                if event is not None:
                    event['traversals'] = traversals
                    event['interpretation'] = 'geometric path prediction, not a tagged packet measurement'
                    events.append(event)
    summary = dict(case=name, parameters=asdict(p), **result['summary'])
    summary.update(initial_separation=float(rows[0, 13]), final_separation=float(rows[-1, 13]),
        maximum_separation=float(np.max(rows[:, 13])), actual_birth_velocities=rows[0, 9:11].tolist(),
        final_positions=rows[-1, 1:3].tolist(), final_momenta=rows[-1, 3:5].tolist(),
        minimum_effective_mass=float(np.min(rows[:, 11:13])),
        maximum_simultaneous_clock_difference=float(np.max(abs(rows[:, 6]-rows[:, 5]))),
        max_neck_energy=float(np.max(result['neck_energy'])),
        port_impulses=np.trapezoid(rows[:, 20:22], rows[:, 0], axis=0).tolist(),
        geometric_supported_events=len(events),
        geometric_advanced_events=sum(e['elapsed_exterior'] < 0 for e in events),
        initial_data_note='Same specified field and nominal momenta; different port Hamiltonians need not have equal initial energies.',
        scope='Constitutive scalar plus Coulomb model; fixed neck, meridional mouths, dynamical clock energies. No GR/QED or topology-change claim.')
    (root/f'{name}.json').write_text(json.dumps(summary, indent=2)+'\n')
    (root/f'{name}_events.json').write_text(json.dumps(events, separators=(',', ':'))+'\n')
    print(json.dumps({k: summary[k] for k in ['case', 'converged', 'iterations', 'stop',
        'final_separation', 'energy_relative_residual', 'angular_absolute_residual',
        'neck_clock_work_residual', 'ledger_pass']}), flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=['all', *CASES], default='all')
    parser.add_argument('--out', type=Path, default=ROOT)
    parser.add_argument('--resume', action='store_true', help='Skip only completed case JSONs')
    args = parser.parse_args()
    for name in CASES if args.case == 'all' else [args.case]:
        if args.resume and (args.out/f'{name}.json').exists():
            print(f'SKIP completed {name}', flush=True)
            continue
        try:
            archive(name, args.out)
        except (ValueError, RuntimeError) as exc:
            failure = dict(case=name, status='FAILED', error=str(exc),
                           parameters=asdict(replace(Parameters(), **CASES[name])))
            args.out.mkdir(parents=True, exist_ok=True)
            (args.out/f'{name}_failure.json').write_text(json.dumps(failure, indent=2)+'\n')
            print(json.dumps(failure), flush=True)
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(args.out.glob('*')) if p.is_file() and p.name != 'sha256.json'}
    (args.out/'sha256.json').write_text(json.dumps(hashes, indent=2)+'\n')


if __name__ == '__main__':
    main()
