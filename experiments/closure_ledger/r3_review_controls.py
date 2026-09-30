"""Retrospective checks of 2e984ac; no nonlinear selection experiment.

The historical producer and archive are unchanged. Replay authenticates and
re-scores the saved diagnostics; it cannot establish shadowing without states.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from geometrodynamics.waves import esu_floquet as fl
from geometrodynamics.waves import r3_resonance as r3
from experiments.closure_ledger import r3_resonance_probe as probe

ARCHIVE = probe.RUN_DIR/'r3_resonance.json'
ARCHIVE_SHA256 = '902c309fb71c6469d85370bc8ec089c256943edce6ae1f9fe118bc511fe8a337'


def exact_event_time(epsilon, eigenvalue):
    """Static ansatz R(tau)x+epsilon*dphi(x)=0 near tau=0, NOT GR evolution.

    R(tau)=-(sqrt(3)/2)sin(2 tau); dphi(x)=eigenvalue*x.
    Require a simple time crossing on the central branch.
    """
    z = 2*np.asarray(epsilon)*np.asarray(eigenvalue)/np.sqrt(3)
    if not np.isfinite(z).all() or np.any(abs(z) >= 1):
        raise ValueError('no simple event on the central branch')
    return .5*np.arcsin(z)


def validate_raw(raw):
    """Reject missing, duplicate, malformed and nonfinite diagnostic records."""
    expected = {(e, p, n) for e in probe.LADDER for p in probe.POLS for n in r3.TOLERANCES}
    rows = raw.get('runs', [])
    keys = [(r['eps'], r['pol'], r['integrator']) for r in rows]
    if len(keys) != len(expected) or set(keys) != expected:
        raise ValueError('incomplete, duplicate or relabelled run schedule')
    trace = raw.get('trace_T2')
    if not np.isfinite(trace) or abs(trace) >= 2:
        raise ValueError('invalid elliptic linear trace')
    for row in rows:
        if row.get('error') is not None:
            raise ValueError('failed run has no complete replay evidence')
        for key, shape in [('increments', (48,)), ('residuals', (48,)),
                           ('A_dev', (48,)), ('radii', (48, 2))]:
            values = np.asarray(row.get(key), dtype=float)
            if values.shape != shape or not np.isfinite(values).all():
                raise ValueError('invalid '+key)
            if key != 'increments' and np.any(values < 0):
                raise ValueError('negative '+key)
        radii = np.asarray(row['radii'])
        if np.any(radii[:, 0] > radii[:, 1]):
            raise ValueError('reversed radius bounds')
        windows = row.get('windows', [])
        if len(windows) != 12:
            raise ValueError('incomplete restart ledger')
        for w in windows:
            if not np.isfinite(w['kick']) or type(w['widenings']) is not int or w['widenings'] < 0:
                raise ValueError('invalid restart ledger')


def replay(path=ARCHIVE):
    # Keep the standalone kinematic helpers free of the symbolic audit imports.
    from experiments.closure_ledger.esu_floquet_probe import close

    blob = Path(path).read_bytes()
    if hashlib.sha256(blob).hexdigest() != ARCHIVE_SHA256:
        raise ValueError('historical archive fingerprint mismatch')
    rec = json.loads(blob)
    if (rec['freeze'], rec['correction']) != (probe.FREEZE, probe.CORRECTION):
        raise ValueError('wrong historical freezes')
    if rec['sources'] != probe.sources():
        raise ValueError('historical producer changed')
    validate_raw(rec['raw'])
    fresh = json.loads(json.dumps(probe.score(rec['raw']), allow_nan=False))
    # Archive bytes stay exact. Recomputed floats can differ with BLAS order;
    # close retains exact structure, labels, booleans and None throughout.
    if fresh['verdict'] != rec['result']['verdict'] or not close(fresh, rec['result'], 1e-9):
        raise ValueError('saved verdict does not reproduce')
    restarts = [abs(w['kick']) for r in rec['raw']['runs'] for w in r['windows'][1:]]
    return dict(archive_sha256=ARCHIVE_SHA256, registered_verdict=fresh['verdict'],
                diagnostic_replay='VERIFIED', trajectory_shadowing='NOT_ESTABLISHED',
                reason='No pre/post restart states or variational sensitivity ledger archived',
                largest_later_A_kick=max(restarts))


def linear_control(method='DOP853', samples=256):
    """Same clock phase/angle/K as the nonlinear run, but only the linear ODE.

    The unwrapped linear path selects the Floquet branch independently of any
    nonlinear amplitude. No zero-amplitude nonlinear angle is evaluated.
    """
    start = np.pi/4
    ts = start+np.arange(48*samples+1)*(np.pi/samples)
    sol = solve_ivp(lambda t, y: fl.rhs_tensor(t, y, 2), (ts[0], ts[-1]), [1., 0.],
                    method=method, t_eval=ts, rtol=1e-12, atol=1e-14)
    if not sol.success or sol.y.shape != (2, len(ts)) or not np.isfinite(sol.y).all():
        raise ArithmeticError('linear control failed')
    f = fl.background(ts)[2]
    wrapped = np.arctan2(-sol.y[1]/(3*f), sol.y[0])
    theta = np.unwrap(wrapped)
    if np.max(abs(np.diff(theta))) >= np.pi/2:
        raise ArithmeticError('angle sampling unresolved')
    increments = np.diff(theta[::samples])
    trace = float(np.trace(fl.monodromy('T', 2)))
    a = np.arccos(trace/2)/(2*np.pi)
    candidates = np.array([1+a, 2-a])
    turns = (theta[-1]-theta[0])/(2*np.pi*48)
    rho0 = float(candidates[np.argmin(abs(candidates-turns))])
    estimates = {str(k): r3.birkhoff(increments[:k]) for k in (24, 48)}
    # A different degree-one phase coordinate: same infinite-time rotation,
    # generally different finite-window estimates.
    theta_alt = theta+.2*np.sin(2*theta)
    alternative = r3.birkhoff(np.diff(theta_alt[::samples]))
    return dict(method=method, samples_per_period=samples, rho0=rho0,
                rho=estimates, bias48=estimates['48']-rho0,
                half_window_difference=abs(estimates['48']-estimates['24']),
                alternative_angle_rho48=alternative,
                increments=increments.tolist())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = dict(scope='Retrospective archive and linear controls only', replay=replay(),
                  linear=[linear_control(), linear_control('RK45', 512)],
                  exact_static_event_times={str(e): float(exact_event_time(e, 1.))
                                            for e in (.001, .01, .05)})
    result['source_sha256'] = {str(Path(__file__).relative_to(probe.ROOT)):
                               hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    output = json.dumps(result, indent=2, allow_nan=False)+'\n'
    if args.output:
        if args.output.resolve() == ARCHIVE.resolve():
            raise ValueError('cannot overwrite the historical archive')
        args.output.parent.mkdir(parents=True, exist_ok=True)
        # Review records are also append-only: choose a new path for a rerun.
        with args.output.open('x') as stream:
            stream.write(output)
    print(output)


if __name__ == '__main__':
    main()
