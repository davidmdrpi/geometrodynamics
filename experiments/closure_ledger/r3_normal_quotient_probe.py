"""Registered normal-quotient experiment and deterministic evidence replay."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from geometrodynamics.waves import r3_normal_quotient as nq, r3_family as rf
from geometrodynamics.waves import nonlinear_supported_tt as d
from experiments.closure_ledger import r3_family_probe as old
from experiments.closure_ledger import r3_family_replay as parent

ROOT = old.ROOT
RUN = ROOT/'experiments/closure_ledger/runs/20261003_r3_normal_quotient'
FREEZE = 'c5c8c6be168408a82b86c9208c5182983c22bffc'
SOURCES = tuple(sorted(set(old.SOURCES + (
    'geometrodynamics/waves/r3_normal_quotient.py',
    'experiments/closure_ledger/r3_normal_quotient_probe.py',
    'experiments/closure_ledger/r3_family_replay.py',
    'docs/r3_normal_quotient_prereg.md'))))
EPS = (2e-6, 1e-6)
STEPS = ((1024, 2048), (2048, 4096))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources():
    return {p: digest(ROOT/p) for p in SOURCES}


def inputs():
    for name, sha in parent.SHA256.items():
        if digest(old.RUN_DIR/name) != sha:
            raise ValueError('parent archive fingerprint mismatch')
    F = json.loads((old.RUN_DIR/'stage_F.json').read_text())
    S = json.loads((old.RUN_DIR/'stage_S.json').read_text())
    if F['sources'] != old.sources() or S['sources'] != old.sources():
        raise ValueError('parent source mismatch')
    parent.validate_s(S, F)
    pts, _ = parent.schedule(F)
    return S['samples'], pts


def args_for(sample, pts):
    j = sample['index']
    return np.r_[sample['v'][:6], np.zeros(6)], pts[(j+1) % len(pts)][:6]-pts[j-1][:6]


def compensated(z):
    y = rf.to_state(z); grav = d.constraints(y)['gravitational_momentum']
    def state(x):
        yy = y.copy(); yy[4:7] = x[:3]; yy[7] = x[3]
        return yy
    sol = root(lambda x: d.constraints(state(x))['residual'], np.r_[-grav/y[7], y[7]], tol=1e-11)
    out = state(sol.x)
    if not np.isfinite(out).all() or np.max(abs(d.constraints(out)['residual'])) > 1e-8 or out[7] >= 0:
        raise ArithmeticError('constraint completion failed')
    return out


def evolve(y):
    times = []
    for _ in range(2):
        s = solve_ivp(lambda t, u: d.conformal_rhs(u), (0., np.pi+.8), y,
                      events=rf._clock, method='DOP853', rtol=1e-12, atol=1e-14, max_step=.025)
        hits = [i for i, t in enumerate(s.t_events[0]) if t > 1.]
        if not s.success or not hits:
            raise ArithmeticError('full return failed')
        y = s.y_events[0][hits[0]]; times.append(float(s.t_events[0][hits[0]]))
    return y, times


def score(raw):
    """Recompute every diagnostic; never accept stored gate booleans."""
    if raw['freeze'] != FREEZE or raw['sources'] != sources() or raw['parent_hashes'] != parent.SHA256:
        raise ValueError('provenance mismatch')
    samples, pts = inputs()
    if len(raw['samples']) != 12 or len(raw['refined']) != 2 or len(raw['perturbations']) != 12:
        raise ValueError('incomplete evidence schedule')
    reduced = []
    for rec, smp in zip(raw['samples'], samples):
        M = nq.finite(rec['M2'], (12, 12))
        if rec['index'] != smp['index'] or not np.array_equal(M, smp['M2']):
            raise ValueError('sample differs from authenticated parent')
        reduced.append(nq.reduce_sample(M, *args_for(smp, pts)))
    z, chord = args_for(samples[0], pts)
    base = reduced[0]; M0 = np.array(samples[0]['M2'])
    refinements = []
    for rec, steps in zip(raw['refined'], STEPS):
        if rec['steps'] != list(steps):
            raise ValueError('wrong refinement schedule')
        M = nq.finite(rec['M2'], (12, 12)); rr = nq.reduce_sample(M, z, chord)
        diff = float(np.linalg.norm(M-M0, 2)/np.linalg.norm(M0, 2))
        a, b = rr['modal'], base['modal']
        trace = abs(a.get('elliptic_trace', 1e20)-b.get('elliptic_trace', -1e20))
        power = abs(max(a.get('power_norms', [1e20]))/max(b.get('power_norms', [1.]))-1)
        refinements.append(dict(steps=list(steps), matrix_relative_error=diff, trace_difference=trace,
                                power_relative_difference=power, reduction=rr,
                                ok=bool(all(rr['checks'].values()) and diff <= 1e-7 and trace <= 1e-6 and power <= 1e-4)))
    Q, U, B = (np.array(base[k]) for k in ('Q', 'U', 'B'))
    errors = []; constraints = []
    for rec, (eps, k) in zip(raw['perturbations'], [(e, k) for e in EPS for k in range(6)]):
        if rec['epsilon'] != eps or rec['direction'] != k:
            raise ValueError('wrong perturbation schedule')
        finals = []
        for side, sign in (('minus', -1), ('plus', 1)):
            ini = nq.finite(rec[side]['initial'], (29,)); fin = nq.finite(rec[side]['final'], (29,))
            times = nq.finite(rec[side]['return_times'], (2,))
            if np.any(times <= 1) or abs(fin[3]) > 1e-8 or fin[7] >= 0 or abs(fin[2]-ini[2]-sum(times)) > 1e-8:
                raise ValueError('invalid return endpoint')
            # Initial compensation shifts section momenta only at second order.
            target = z+sign*eps*(Q@U[:, k])
            expected = compensated(target)
            if np.max(abs(ini-expected)) > 1e-12:
                raise ValueError('initial state is not the scheduled physical perturbation')
            constraints.extend(float(np.max(abs(d.constraints(y)['residual']))) for y in (ini, fin))
            finals.append(rf.to_section(fin))
        response = Q.T@(finals[1]-finals[0])/(2*eps)
        error = float(np.linalg.norm(response-B@U[:, k])/max(1., np.linalg.norm(B@U[:, k])))
        errors.append(dict(epsilon=eps, direction=k, relative_error=error))
    checks = {key: bool(all(r['checks'][key] for r in reduced)) for key in ('Q1', 'Q2', 'Q3', 'Q4')}
    checks['Q5'] = all(r['ok'] for r in refinements)
    checks['Q6'] = max(e['relative_error'] for e in errors) <= 1e-4 and max(constraints) <= 1e-8
    return dict(checks=checks, verdict='BOUNDED_CENTER_QUOTIENT_NUMERICALLY' if all(checks.values()) else 'NORMAL_RESPONSE_UNRESOLVED',
                scope='conditional linear normal response; no nonlinear or action-selection claim',
                unreduced='UNREDUCED_HYPERBOLIC_PAIR_PRESENT' if all(r['checks']['Q3'] for r in reduced) else 'UNRESOLVED',
                reduced=reduced, refinements=refinements, perturbation_errors=errors,
                max_constraint=max(constraints))


def produce(directory=RUN):
    if directory.exists():
        raise FileExistsError('append-only output directory exists')
    parent.replay()
    src = sources(); samples, pts = inputs()
    raw = dict(freeze=FREEZE, sources=src, parent_hashes=parent.SHA256,
               samples=[dict(index=s['index'], M2=s['M2']) for s in samples], refined=[], perturbations=[])
    z, chord = args_for(samples[0], pts)
    base = nq.reduce_sample(samples[0]['M2'], z, chord)
    Q, U = np.array(base['Q']), np.array(base['U'])
    if U.shape[1] != 6:
        raise ArithmeticError('centre dimension differs; cannot run registered six-direction schedule')
    for steps in STEPS:
        js = [rf.DP(np.r_[samples[0]['v'][i:i+6], np.zeros(6)], 12, steps=steps)[1] for i in (0, 6)]
        raw['refined'].append(dict(steps=list(steps), M2=(js[1]@js[0]).tolist()))
        print('completed refinement', steps, flush=True)
    for eps in EPS:
        for k in range(6):
            rec = dict(epsilon=eps, direction=k)
            for side, sign in (('minus', -1), ('plus', 1)):
                ini = compensated(z+sign*eps*(Q@U[:, k])); fin, times = evolve(ini)
                rec[side] = dict(initial=ini.tolist(), final=fin.tolist(), return_times=times)
            raw['perturbations'].append(rec)
    result = score(raw)
    if sources() != src:
        raise ValueError('sources changed during measurement')
    directory.mkdir(parents=True)
    for name, obj in (('raw.json', raw), ('result.json', result)):
        (directory/name).write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')
    manifest = {n: digest(directory/n) for n in ('raw.json', 'result.json')}
    (directory/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('checks', 'verdict', 'max_constraint')}, indent=2))


def replay(directory=RUN, manifest_sha=None):
    if manifest_sha is None:
        raise ValueError('a pinned manifest digest is required')
    if digest(directory/'manifest.json') != manifest_sha:
        raise ValueError('manifest fingerprint mismatch')
    manifest = json.loads((directory/'manifest.json').read_text())
    if set(manifest) != {'raw.json', 'result.json'}:
        raise ValueError('wrong evidence inventory')
    for name, sha in manifest.items():
        if digest(directory/name) != sha:
            raise ValueError('evidence fingerprint mismatch: '+name)
    raw = json.loads((directory/'raw.json').read_text()); fresh = score(raw)
    saved = json.loads((directory/'result.json').read_text())
    if fresh['checks'] != saved['checks'] or fresh['verdict'] != saved['verdict']:
        raise ValueError('decision does not reproduce')
    return fresh


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('mode', choices=['run', 'replay'])
    ap.add_argument('--manifest-sha'); args = ap.parse_args()
    if args.mode == 'run':
        produce()
    else:
        out = replay(manifest_sha=args.manifest_sha)
        print(json.dumps({k: out[k] for k in ('checks', 'verdict', 'max_constraint')}, indent=2))
