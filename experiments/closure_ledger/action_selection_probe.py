"""Frozen nonlinear canonical-action test with a separate receiver-readiness verdict."""
import argparse
import base64
import gzip
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from geometrodynamics.waves import action_selection as a

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'experiments/closure_ledger/runs/20260928_action_selection'


def encode_states(raw):
    data = json.dumps(raw, separators=(',', ':'), allow_nan=False).encode()
    return base64.b64encode(gzip.compress(data, mtime=0)).decode()+'\n'


def decode_states(path):
    return json.loads(gzip.decompress(base64.b64decode(Path(path).read_text())))


def validate_schedule(raw):
    """A replay cannot silently relabel times, amplitudes or initial states."""
    if raw.get('nodes') != 64 or not np.array_equal(raw.get('times'), a.TIMES):
        raise ValueError('raw archive does not match frozen node/time schedule')
    if tuple(row['epsilon'] for row in raw.get('cases', [])) != a.AMPLITUDES:
        raise ValueError('raw archive does not match frozen amplitude schedule')
    for epsilon, values in [(row['epsilon'], row['states']) for row in raw['cases']]+[(.08, raw['rk45_states'])]:
        states = np.asarray(values, dtype=float)
        if states.shape != (64, len(a.TIMES), 29) or not np.isfinite(states).all():
            raise ValueError('raw state dimensions or values invalid')
        if not np.allclose(states[:, 0], a.initial_loop(epsilon), rtol=1e-12, atol=1e-14):
            raise ValueError('raw initial states differ from frozen preparation')


def summarize(raw):
    validate_schedule(raw)
    result = dict(freeze=a.FREEZE, observable='closed preparation-loop circulation, not absorbed action',
                  times=a.TIMES.tolist(), cases=[], receiver_selection=a.receiver_readiness(),
                  numerical_gates={})
    gates = dict(constraints=True, chart=True, initial_action=True, full_circulation=True,
                 zero_action=True, loop_refinement=True, integrator_refinement=True)
    for item in raw['cases']:
        eps, states = item['epsilon'], np.array(item['states'])
        pred = a.initial_prediction(eps)
        rows = []
        for nodes in (16, 32, 64):
            action = a.circulation(states[::64//nodes])
            rows.append(dict(nodes=nodes, action={k:v.tolist() for k,v in action.items()},
                             change={k:(v-v[0]).tolist() for k,v in action.items()}))
        d = a.diagnostics(states)
        fine = rows[-1]['action']
        err = [[max(abs(np.array(rows[j+1]['action'][k])-rows[j]['action'][k])) for k in a.SECTORS]
               for j in (0,1)]
        initial_err = abs(fine['total'][0]-pred)
        drift = float(np.max(abs(np.array(fine['total'])-fine['total'][0])))
        gates['constraints'] &= d['sampled_constraint_max'] < 1e-9
        gates['chart'] &= d['sampled_det_error'] < 1e-8 and d['sampled_symmetry_error'] < 1e-8 and min(d['sampled_minimum_A_H_M']) > 0
        if eps:
            gates['initial_action'] &= initial_err/pred < 1e-8
            gates['full_circulation'] &= drift/pred < 1e-6
            gates['loop_refinement'] &= max(err[-1])/pred < 1e-5
        else:
            gates['zero_action'] &= max(abs(v) for k in a.SECTORS for v in fine[k]) < 1e-12
        result['cases'].append(dict(epsilon=eps, prediction=pred, diagnostics=d, quadratures=rows,
                                    initial_absolute_error=initial_err, circulation_absolute_drift=drift,
                                    sector_order=list(a.SECTORS), quadrature_max_differences=err))
    reference = result['cases'][-1]['quadratures'][-1]['action']
    other = a.circulation(np.array(raw['rk45_states']))
    difference = {k:float(np.max(abs(other[k]-reference[k]))) for k in a.SECTORS}
    rk_diag = a.diagnostics(np.array(raw['rk45_states']))
    gates['constraints'] &= rk_diag['sampled_constraint_max'] < 1e-9
    gates['chart'] &= rk_diag['sampled_det_error'] < 1e-8 and rk_diag['sampled_symmetry_error'] < 1e-8 and min(rk_diag['sampled_minimum_A_H_M']) > 0
    gates['integrator_refinement'] = max(difference.values())/a.initial_prediction(.08) < 1e-5
    result['integrator_control'] = dict(epsilon=.08, method='RK45', diagnostics=rk_diag,
                                       action={k:v.tolist() for k,v in other.items()}, max_differences=difference)
    # Signed-sector slopes are omitted near zero or sign changes, never clipped.
    slopes = []
    for left, right in zip(result['cases'][1:-1], result['cases'][2:]):
        row = dict(amplitudes=[left['epsilon'],right['epsilon']], log_slopes={})
        for key in a.SECTORS:
            low = np.array(left['quadratures'][-1]['action'][key])
            high = np.array(right['quadratures'][-1]['action'][key])
            floor = 1e-8*max(left['prediction'],right['prediction'])
            row['log_slopes'][key] = [float(np.log(abs(h/l))/np.log(2)) if abs(l)>floor and abs(h)>floor and l*h>0 else None for l,h in zip(low,high)]
        slopes.append(row)
    result['amplitude_log_slopes'] = slopes
    result['numerical_gates'] = {k:bool(v) for k,v in gates.items()}
    result['numerical_verdict'] = 'PASS_CANONICAL_CONTINUUM_CONTROL' if all(gates.values()) else 'REGISTERED_NUMERICAL_FAILURE'
    result['physical_scope'] = 'Continuous preparation-loop action in this sector; receiver action selection remains untested.'
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=RUN)
    parser.add_argument('--replay', type=Path, help='recompute diagnostics from a .json.gz.b64 raw-state archive')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.replay:
        raw = decode_states(args.replay)
    else:
        raw = dict(times=a.TIMES.tolist(), nodes=64, cases=[])
        for epsilon in a.AMPLITUDES:
            print('Evolving epsilon', epsilon, flush=True)
            states = a.evolve_loop(epsilon)
            raw['cases'].append(dict(epsilon=epsilon, states=states.tolist()))
        print('Independent RK45 control', flush=True)
        raw['rk45_states'] = a.evolve_loop(.08, method='RK45').tolist()
    archive = args.output_dir/'states.json.gz.b64'
    archive.write_text(encode_states(raw))
    result = summarize(raw)
    result['environment'] = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)
    result['source_sha256'] = {}
    for p in ['geometrodynamics/waves/action_selection.py', 'geometrodynamics/waves/nonlinear_supported_tt.py',
              'experiments/closure_ledger/action_selection_probe.py','docs/dynamical_action_selection_prereg.md']:
        result['source_sha256'][p] = hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
    result['archive_sha256'] = hashlib.sha256(archive.read_bytes()).hexdigest()
    (args.output_dir/'action.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('numerical_verdict','numerical_gates','receiver_selection')}, indent=2))
    return 0 if result['numerical_verdict']=='PASS_CANONICAL_CONTINUUM_CONTROL' else 1


if __name__ == '__main__':
    raise SystemExit(main())
