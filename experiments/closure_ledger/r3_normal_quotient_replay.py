"""Portable authenticated replay, added after the registered measurement.

The producer recomputed a nonunique Schur frame when matching perturbations.
Different SciPy builds can rotate that frame. This replay uses the authenticated
recorded frame after checking orthonormality and projector agreement with the
independently recomputed quotient and centre. It also removes the producer's
unit floor in the response-error denominator to enforce a strictly relative
Q6 error. Thresholds and recorded decisions are unchanged.
"""
import argparse
import json
import numpy as np
from geometrodynamics.waves import r3_normal_quotient as nq, r3_family as rf
from geometrodynamics.waves import nonlinear_supported_tt as d
from experiments.closure_ledger.r3_normal_quotient_probe import (
    FREEZE, RUN, STEPS, EPS, sources, inputs, args_for, compensated, digest, parent)

MANIFEST_SHA = '132a99d9f051d7e96f96c5b8d8673d7e9deeffe045342ea1563c66ebe1a2f163'

def score(raw, saved):
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
    # A Schur basis in a repeated neutral eigenspace is not unique. Authenticate
    # the saved frame, then verify its subspace against this environment's frame.
    Q = nq.finite(saved['reduced'][0]['Q'], (12, 8))
    U = nq.finite(saved['reduced'][0]['U'], (8, 6))
    Qnew, Unew = (np.array(base[k]) for k in ('Q', 'U'))
    lift = Q@U
    for defect in (np.linalg.norm(Q.T@Q-np.eye(8), 2),
                   np.linalg.norm(U.T@U-np.eye(6), 2),
                   np.linalg.norm(Q@Q.T-Qnew@Qnew.T, 2),
                   np.linalg.norm(lift@lift.T-(Qnew@Unew)@(Qnew@Unew).T, 2)):
        if defect > 1e-7:
            raise ValueError('saved perturbation frame is not the verified quotient centre')
    B = Q.T@M0@Q
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
        scale = np.linalg.norm(B@U[:, k])
        if scale <= 1e-12:
            raise ValueError('zero predicted response')
        error = float(np.linalg.norm(response-B@U[:, k])/scale)
        errors.append(dict(epsilon=eps, direction=k, relative_error=error))
    checks = {key: bool(all(r['checks'][key] for r in reduced)) for key in ('Q1', 'Q2', 'Q3', 'Q4')}
    checks['Q5'] = all(r['ok'] for r in refinements)
    checks['Q6'] = max(e['relative_error'] for e in errors) <= 1e-4 and max(constraints) <= 1e-8
    return dict(checks=checks, verdict='BOUNDED_CENTER_QUOTIENT_NUMERICALLY' if all(checks.values()) else 'NORMAL_RESPONSE_UNRESOLVED',
                scope='conditional linear normal response; no nonlinear or action-selection claim',
                unreduced='UNREDUCED_HYPERBOLIC_PAIR_PRESENT' if all(r['checks']['Q3'] for r in reduced) else 'UNRESOLVED',
                reduced=reduced, refinements=refinements, perturbation_errors=errors,
                max_constraint=max(constraints))


def replay(directory=RUN):
    if digest(directory/'manifest.json') != MANIFEST_SHA:
        raise ValueError('manifest fingerprint mismatch')
    manifest = json.loads((directory/'manifest.json').read_text())
    if set(manifest) != {'raw.json', 'result.json'}:
        raise ValueError('wrong evidence inventory')
    for name, sha in manifest.items():
        if digest(directory/name) != sha:
            raise ValueError('evidence fingerprint mismatch: '+name)
    raw = json.loads((directory/'raw.json').read_text())
    saved = json.loads((directory/'result.json').read_text())
    fresh = score(raw, saved)
    if fresh['checks'] != saved['checks'] or fresh['verdict'] != saved['verdict'] or fresh['unreduced'] != saved['unreduced']:
        raise ValueError('decision does not reproduce')
    for a, b in zip(fresh['reduced'], saved['reduced']):
        # Basis-free diagnostics; do not compare arbitrary SVD/eigenvector signs.
        if not np.allclose(a['modal']['power_norms'], b['modal']['power_norms'], rtol=1e-6, atol=1e-8):
            raise ValueError('power envelope does not reproduce')
        if abs(a['modal']['elliptic_trace']-b['modal']['elliptic_trace']) > 1e-6:
            raise ValueError('elliptic trace does not reproduce')
    return fresh


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); ap.parse_args()
    r = replay()
    print(json.dumps({k: r[k] for k in ('checks', 'verdict', 'max_constraint')}, indent=2))
