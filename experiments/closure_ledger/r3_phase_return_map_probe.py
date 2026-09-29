"""Frozen cubic return-map experiment; no long trajectory or state resetting."""
import argparse
import base64
import gzip
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from geometrodynamics.waves import r3_phase_return_map as rm
from geometrodynamics.waves import esu_floquet as fl
from geometrodynamics.waves.taylor_jets import Jet, ring

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'experiments/closure_ledger/runs/20260929_r3_phase_return_map'
RADII = [.008, .004, .002]
ANGLES = np.arange(8)*np.pi/4
SOURCES = ['geometrodynamics/waves/taylor_jets.py', 'geometrodynamics/waves/r3_phase_return_map.py',
           'geometrodynamics/waves/nonlinear_supported_tt.py', 'geometrodynamics/waves/esu_floquet.py',
           'experiments/closure_ledger/r3_phase_return_map_probe.py', 'docs/r3_return_map_prereg.md']


def serial(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def source_hashes():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def embedding(reduction):
    h = [Jet(ring(2), c) for c in reduction['graph']]
    w = [ring(2).variable(i) for i in range(2)]
    return h+rm.transform(reduction['elliptic_basis'], w)


def measure():
    raw = dict(freeze=rm.FREEZE, jets=[], validation=[])
    for method in rm.METHODS:
        print('Integrating cubic variational map:', method, flush=True)
        raw['jets'].append(rm.integrate_jets(method))
    reduction = rm.reduce_map(raw['jets'][0]['coefficients'])
    W = embedding(reduction)
    for angle in ANGLES:
        print('Independent full-state validation angle:', angle, flush=True)
        for radius in RADII:
            w = radius*np.array([np.cos(angle), np.sin(angle)])
            z = np.array([p.evaluate(w) for p in W])
            raw['validation'].append(dict(angle=float(angle), radius=radius, z=z,
                                          history=rm.direct_return(z)))
    return raw


def analyze(raw):
    if raw.get('freeze') != rm.FREEZE or [r['method'] for r in raw['jets']] != list(rm.METHODS):
        raise ValueError('wrong freeze or integrator schedule')
    for row in raw['jets']:
        c = np.asarray(row['coefficients'])
        if c.shape != (4, ring(4).size) or not np.isfinite(c).all():
            raise ValueError('invalid cubic coefficients')
        t = np.asarray(row['return_time'])
        if t.shape != (ring(4).size,) or not np.isfinite(t).all():
            raise ValueError('invalid return-time coefficients')
    reductions = [rm.reduce_map(r['coefficients']) for r in raw['jets']]
    diagnostics = []
    tensor_trace = float(np.trace(fl.monodromy('T', 2)))
    for row, reduced in zip(raw['jets'], reductions):
        L = reduced['linear']
        absolute, scaled = rm.symplectic_residual(row['coefficients'])
        diagnostics.append(dict(background_error=float(np.max(abs(np.asarray(row['coefficients'])[:, 0]))),
                                background_time_error=float(abs(row['return_time'][0]-np.pi)),
                                offblock=float(max(np.max(abs(L[:2, 2:])), np.max(abs(L[2:, :2])))),
                                linear_symplectic=float(np.max(abs(L.T@rm.J4@L-rm.J4))),
                                linear_trace_error=float(abs(np.trace(L[2:, 2:])-tensor_trace)),
                                symplectic_absolute=absolute, symplectic_scaled=scaled))
    primary = reductions[0]
    F = rm.map_polynomials(raw['jets'][0]['coefficients'])
    W = embedding(primary)
    center = [Jet(ring(2), c) for c in primary['center']]
    expected = [(float(a), r) for a in ANGLES for r in RADII]
    if [(r['angle'], r['radius']) for r in raw['validation']] != expected:
        raise ValueError('missing, duplicate or relabelled validation points')
    rows = []
    for record in raw['validation']:
        angle, radius, history = record['angle'], record['radius'], record['history']
        w = radius*np.array([np.cos(angle), np.sin(angle)])
        z = np.array([p.evaluate(w) for p in W])
        if not np.allclose(record['z'], z, atol=1e-12, rtol=1e-12):
            raise ValueError('incorrect centre preparation')
        y0, returned = np.asarray(history['initial']), np.asarray(history['returned'])
        samples, times = np.asarray(history['samples']), np.asarray(history['sample_times'])
        if (samples.shape != (257, 29) or times.shape != (257,) or
                not np.isfinite(samples).all() or not np.isfinite(times).all()):
            raise ValueError('invalid validation history')
        if not np.allclose(y0, rm.initial_full(z), rtol=1e-11, atol=1e-13):
            raise ValueError('changed full initial data')
        if not np.allclose(samples[[0, -1]], [y0, returned], rtol=1e-12, atol=1e-13):
            raise ValueError('endpoint mismatch')
        if not np.allclose(times, np.linspace(0., history['return_time'], 257), atol=1e-13, rtol=0):
            raise ValueError('changed validation time schedule')
        exact = rm.canonical_state(returned)
        if abs(returned[3]) > 1e-10 or returned[7] >= 0:
            raise ValueError('wrong return section')
        residual = max(float(np.max(rm.full.constraints(y)['normalized'])) for y in samples)
        chart = min(min(y[0], rm.full.ingredients(y)[0], np.linalg.eigvalsh(rm.full.unpack(y)[4]).min()) for y in samples)
        formulation = max(np.linalg.norm(exact-history['phase_canonical']), abs(history['return_time']-history['phase_time']))
        predicted = np.array([p.evaluate(z) for p in F])
        center_at = np.linalg.solve(primary['elliptic_basis'], exact[2:])
        on_graph = np.array([p.evaluate(center_at) for p in W[:2]])
        rows.append(dict(angle=angle, radius=radius, map_error=float(np.linalg.norm(exact-predicted)),
                         graph_error=float(np.linalg.norm(exact[:2]-on_graph)),
                         formulation_error=float(formulation), constraint_max=residual, chart_min=float(chart)))
    remainders = []
    for i, angle in enumerate(ANGLES):
        group = rows[3*i:3*i+3]
        for key in ['map_error', 'graph_error']:
            values = [r[key] for r in group]
            slopes = [float(np.log2(a/b)) if a > 0 and b > 0 else None for a, b in zip(values, values[1:])]
            passed = values[-1] < 1e-6 and all(b < 1e-10 or (s is not None and s >= 3.5)
                                              for b, s in zip(values[1:], slopes))
            remainders.append(dict(angle=float(angle), observable=key, errors=values,
                                   orders=slopes, passes=passed))
    error = abs(reductions[0]['nu']-reductions[1]['nu'])
    gates = dict(
        G1=all(d['background_error'] < 1e-9 and d['background_time_error'] < 1e-9 and d['offblock'] < 1e-9
               and d['linear_trace_error'] < 1e-9 and d['linear_symplectic'] < 1e-8
               and abs(r['rho0']-1.484666408416) < 1e-9 for d, r in zip(diagnostics, reductions)),
        G2=all(d['symplectic_absolute'] < 1e-4 and d['symplectic_scaled'] < 1e-9 for d in diagnostics),
        G3=all(max(r['homological_conditions']) < 1e10 and max(r['homological_residuals']) < 1e-7
               and r['graph_residual'] < 1e-7 and abs(r['radial_resonant']) < 1e-7 for r in reductions),
        G4=error < 1e-5*max(1., abs(primary['nu'])) and abs(primary['nu']) > 100*max(error, 1e-9),
        G5=all(r['formulation_error'] < 1e-9 and r['constraint_max'] < 1e-9 and r['chart_min'] > 0 for r in rows)
           and all(r['passes'] for r in remainders))
    verdict = 'UNRESOLVED'
    if all(gates.values()):
        verdict = ('SHIFT_TOWARD_TARGET' if np.sign(primary['nu']) == np.sign(1.5-primary['rho0'])
                   else 'SHIFT_AWAY_FROM_TARGET')
    return dict(freeze=rm.FREEZE, nu=primary['nu'], rho0=primary['rho0'], gates=gates, verdict=verdict,
                integrator_nu_difference=float(error),
                full_jet_difference=float(np.max(abs(np.asarray(raw['jets'][0]['coefficients'])-raw['jets'][1]['coefficients']))),
                reductions=reductions, diagnostics=diagnostics, validation=rows, remainders=remainders,
                resonance_crossing='NOT_TESTED', full_state_closure='NOT_TESTED', action_selection='NOT_ESTABLISHED')


def load_raw(path):
    return json.loads(gzip.decompress(base64.b64decode(Path(path).read_bytes())))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=RUN)
    parser.add_argument('--replay', type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rawpath, reportpath = args.output_dir/'raw.json.gz.b64', args.output_dir/'report.json'
    if rawpath.exists() or reportpath.exists():
        raise ValueError('output exists: preserve evidence and choose another directory')
    sources = source_hashes()
    raw = load_raw(args.replay) if args.replay else measure()
    if source_hashes() != sources:
        raise RuntimeError('sources changed during measurement')
    rawbytes = (base64.b64encode(gzip.compress(json.dumps(raw, default=serial, separators=(',', ':'), allow_nan=False).encode(), mtime=0))+b'\n')
    rawpath.write_bytes(rawbytes)
    report = analyze(raw)
    report.update(source_sha256=sources, raw_sha256=hashlib.sha256(rawbytes).hexdigest(),
                  environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__))
    reportpath.write_text(json.dumps(report, default=serial, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: report[k] for k in ['nu', 'rho0', 'gates', 'verdict']}, indent=2))
    return int(report['verdict'] == 'UNRESOLVED')


if __name__ == '__main__':
    raise SystemExit(main())
