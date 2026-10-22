"""Replay the locally frozen scalar candidate; write all controls and failures."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path

import numpy as np

from geometrodynamics.waves.coupled_mouth_throat import CoupledMouthThroat, Parameters


CASES = dict(baseline={}, no_wave={'energy': 0.}, no_offset={'offset': 0.},
    mistimed={'offset': -.6}, fixed_mouths={'fixed_mouths': True},
    fixed_radius={'fixed_radius': True}, dephased={}, tight_tolerance={},
    higher_resolution={'lmax': 14, 'nodes': 25})
DEFAULT_ROOT = Path('experiments/closure_ledger/runs/20261020_coupled_mouth_throat_replay')


def run_case(name, root=DEFAULT_ROOT):
    root.mkdir(parents=True, exist_ok=True)
    model = CoupledMouthThroat(replace(Parameters(), **CASES[name]))
    rtol = 2e-10 if name == 'tight_tolerance' else 2e-8
    t, states, stop = model.evolve(rtol=rtol, dephased=name == 'dephased')
    rows, totals, angular = [], [], []
    for ti, state in zip(t, states):
        _, _, q, p = model.unpack(state)
        _, energy = model.evaluate(state, True)
        sep = np.arccos(np.cos(q[1]-q[0]))
        speed = np.linalg.norm(p[:2])/model.par.mass
        rows.append([ti, *q, *p, sep, speed, *energy.values(), state[-1]])
        totals.append(sum(energy.values()))
        angular.append(model.angular(state)+state[-1])
    rows, totals, angular = map(np.asarray, (rows, totals, angular))
    columns = 'tau,qa,qb,radius,Pa,Pb,Pr,separation,speed,'+','.join(energy)+',clock_constraint_angular_reaction'
    np.savetxt(root/f'{name}.csv', rows, delimiter=',', header=columns, comments='', fmt='%.12g')
    keep = np.unique(np.r_[np.arange(0, len(t), 4), len(t)-1])
    np.savez_compressed(root/f'{name}.npz', t=t[keep], states=states[keep], initial=states[0], final=states[-1])
    windows = []
    for j in range(3):
        start, end = 2*np.pi*j, 2*np.pi*(j+1)
        if t[-1] < end-1e-10:
            break
        a, b = np.argmin(abs(t-start)), np.argmin(abs(t-end))
        ds, dv = rows[b, 7]-rows[a, 7], rows[b, 8]-rows[a, 8]
        windows.append(dict(start=start, end=end, separation_change=float(ds), speed_change=float(dv),
            passes=bool(ds < -.001 and dv > .001)))
    relative = float(np.max(abs(totals-totals[0]))/max(abs(totals[0]), 1e-30))
    result = dict(case=name, parameters=asdict(model.par), rtol=rtol, stop=stop,
        final_time=float(t[-1]), initial_separation=float(rows[0, 7]), final_separation=float(rows[-1, 7]),
        minimum_radius_ratio=float(rows[:, 3].min()/model.par.radius0), maximum_speed=float(rows[:, 8].max()),
        maximum_energy_relative_residual=relative,
        maximum_angular_balance_residual=float(np.max(abs(angular-angular[0]))),
        maximum_throat_energy_fraction=float(rows[:, 10].max()/max(totals[0], 1e-30)),
        initial_energy=float(totals[0]), final_clock_constraint_reaction=float(states[-1, -1]),
        return_windows=windows,
        acceleration_gate=bool(len(windows) == 3 and all(w['passes'] for w in windows)),
        shrinkage_gate=bool(rows[:, 3].min() <= .75*model.par.radius0),
        cutoff_reached=stop == 'radius_cutoff', energy_gate=relative < 1e-5,
        complete_mechanism=False,
        scope='Unholed sphere with constitutive finite ports; not glued MTY worldtubes or a GR solution.')
    (root/f'{name}.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('case', 'stop', 'final_separation', 'minimum_radius_ratio',
        'maximum_energy_relative_residual', 'maximum_angular_balance_residual')}), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=list(CASES)+['all'], default='all')
    parser.add_argument('--out', type=Path, default=DEFAULT_ROOT)
    args = parser.parse_args()
    for case in CASES if args.case == 'all' else [args.case]:
        run_case(case, args.out)
