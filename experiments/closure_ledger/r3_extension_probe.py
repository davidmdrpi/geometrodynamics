"""R3 extension probe (docs/r3_extension_prereg.md).

part B: leading twist for all five homogeneous n=2 polarisations.
part A: LRS invariant-circle family continued to large amplitude.
score : per-part labels plus the closure-condition table (part D).
Archives are append-only and bind their sources.
"""
import argparse
import hashlib
import json
import time
from fractions import Fraction
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from geometrodynamics.waves import esu_floquet as fl
from geometrodynamics.waves import r3_extension as rx
from geometrodynamics.waves import r3_return_map as rm

ROOT = Path(__file__).resolve().parents[2]
RUN_DIR = ROOT/'experiments/closure_ledger/runs/20260929_r3_extension'
SOURCES = ('geometrodynamics/waves/jets.py', 'geometrodynamics/waves/r3_return_map.py',
           'geometrodynamics/waves/r3_extension.py', 'geometrodynamics/waves/nonlinear_supported_tt.py',
           'geometrodynamics/waves/esu_floquet.py', 'experiments/closure_ledger/r3_extension_probe.py')
NU_316 = -0.9501180968942993
STEPS = (1024, 2048, 4096)
A_START, A_FACTOR, A_CAP, GRID = .004, 2**.25, 2., 63
Q_MAX = 12


def sources():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def _write(path, rec):
    if path.exists():
        raise FileExistsError('archives are append-only: '+str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec, indent=1, allow_nan=False)+'\n')


def _cplx(T):
    return {'re': np.real(T).tolist(), 'im': np.imag(T).tolist()}


def measure_b():
    runs = {n: rx.jet_return_map(3, n) for n in STEPS}
    R12 = rm.richardson(runs[1024]['P'], runs[2048]['P'])
    R24 = rm.richardson(runs[2048]['P'], runs[4096]['P'])
    z0 = np.r_[1., np.zeros(11)]
    nf, nf12 = rx.multimode_normal_form(R24, z0), rx.multimode_normal_form(R12, z0)
    plain = {k: ([v.real, v.imag] if isinstance(v, complex) else v) for k, v in nf.items() if k != 'T'}
    return dict(normal_form=plain, T=_cplx(nf['T']),
                T_coarse=_cplx(nf12['T']), c_coarse=nf12['c'],
                constraint_max=float(abs(runs[4096]['constraint'].c).max()),
                section_q_max=float(abs(runs[4096]['q'].c).max()),
                coef_scale=max(float(abs(p.c).max()) for p in R24))


def score_b(rec, theta0):
    nf = rec['normal_form']
    T = np.array(rec['T']['re'])+1j*np.array(rec['T']['im'])
    Tc = np.array(rec['T_coarse']['re'])+1j*np.array(rec['T_coarse']['im'])
    c = nf['c']
    sup, sup_c = rx.sup_nu(T, c), rx.sup_nu(Tc, rec['c_coarse'])
    rng = np.random.default_rng(29)
    cov = 0.
    for _ in range(5):
        Rq, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        Rq *= np.sign(np.linalg.det(Rq))
        D = rx.so3_rep(Rq)
        for _ in range(5):
            w = rng.normal(size=5)+1j*rng.normal(size=5)
            cov = max(cov, abs(rx.nu_of(T, D @ w, c)-rx.nu_of(T, w, c)))
    lin = [rx.nu_of(T, rng.normal(size=5), c) for _ in range(20)]
    nu_lrs = rx.nu_of(T, np.eye(5)[0], c)
    scale = max(abs(sup['max']), abs(sup['min']))
    conv = max(abs(sup['max']-sup_c['max']), abs(sup['min']-sup_c['min']), abs(nu_lrs-rx.nu_of(Tc, np.eye(5)[0], rec['c_coarse'])))
    checks = dict(
        B1_fixed_point=nf['fixed_error'] <= 1e-10,
        B2_linear=abs(nf['theta']-theta0) <= 1e-9 and nf['off_block'] <= 1e-9 and nf['block_spread'] <= 1e-9,
        B3_dissipative=sup['dissipative_max'] <= 1e-6*sup['twist_scale'],
        B4_so3_covariance=cov <= 1e-7*scale,
        B5_lrs_matches_316=abs(nu_lrs-NU_316) <= 1e-8,
        B6_linear_polarisations_equal=max(lin)-min(lin) <= 1e-7*scale,
        B7_convergence=conv <= 1e-7*scale,
        B8_nonresonance=nf['min_relative_divisor'] >= .1 and min(nf['resonance_distances']) >= .1,
        B9_constraint=rec['constraint_max'] <= 1e-8*max(1., rec['coef_scale']) and rec['section_q_max'] <= 1e-8*max(1., rec['coef_scale']),
    )
    fails = [k for k, v in checks.items() if not v]
    resolved = abs(sup['max']) >= 100*max(conv, 1e-12)
    if fails or not resolved:
        label = 'UNRESOLVED'
    elif sup['max'] < 0:
        label = 'ALL_POLARISATIONS_SHIFT_AWAY'
    else:
        label = 'SOME_POLARISATION_SHIFTS_TOWARD'
    return dict(label=label, checks=checks, failures=fails, nu_lrs=nu_lrs, nu_max=sup['max'], nu_min=sup['min'],
                argmax=sup['max_polarisation'], argmin=sup['min_polarisation'],
                sym2_upper_bound=sup['sym2_upper_bound'], sym2_lower_bound=sup['sym2_lower_bound'],
                certified_by_sym2_bound=bool(sup['sym2_upper_bound'] < 0), convergence=conv,
                so3_defect=cov, linear_polarisation_spread=max(lin)-min(lin), dissipative_max=sup['dissipative_max'])


def measure_a():
    circles, status = [], 'CAP_REACHED'
    with Pool(4) as pool:
        J = rm.fd_jacobian(rx._P, rm.ZSTAR)
        K, omega = rm.linear_circle(J, A_START, GRID)
        a_prev, a = None, A_START
        while a <= A_CAP*(1+1e-12):
            K0 = K if a_prev is None else rm.ZSTAR+(K-rm.ZSTAR)*(a/a_prev)
            rec, err = None, None
            t = time.time()
            try:
                c = rx.invariant_circle_parallel(a, K0, omega, GRID, pool)
                res = pool.map(rm.esu_map, list(c['K']))
                cmax = max(r[1] for r in res)
                ok = c['residual'] <= 1e-10 and c['fourier_tail'] <= 1e-10 and cmax <= 1e-10
                rec = dict(a=a, omega=c['omega'], action=c['action'], residual=c['residual'],
                           fourier_tail=c['fourier_tail'], constraint_max=cmax, iterations=c['iterations'],
                           history=c['history'], seconds=time.time()-t, K=c['K'].tolist(), ok=bool(ok))
            except (ArithmeticError, ValueError, np.linalg.LinAlgError) as e:
                ok, err = False, f'{type(e).__name__}: {e}'
            if ok:
                circles.append(rec)
                K, omega, a_prev = c['K'], c['omega'], a
                if omega >= np.pi-1e-9:
                    status = 'REACHED_PI'
                    break
                a = a*A_FACTOR
                continue
            circles.append(dict(a=a, ok=False, error=err, **({k: rec[k] for k in ('omega', 'action', 'residual', 'fourier_tail', 'constraint_max', 'iterations')} if rec else {})))
            if a_prev is not None and abs(a/a_prev-A_FACTOR) < 1e-12:
                a = a_prev*2**.125            # one registered retry at the half step
                continue
            status = 'FAMILY_ENDED'
            break
    return dict(circles=circles, status=status)


def score_a(rec, theta0):
    good = [c for c in rec['circles'] if c['ok']]
    w = np.array([c['omega'] for c in good])
    I = np.array([c['action'] for c in good])
    order = np.argsort(I)
    w, I = w[order], I[order]
    steps = np.diff(w)
    turn = bool(np.any(steps > 1e-9))
    reached = bool(np.any(w >= np.pi-1e-9))
    if len(good) < 8:
        label = 'UNRESOLVED'
    elif reached:
        label = 'CROSSING_IN_FAMILY'
    elif turn:
        label = 'TURN_IN_FAMILY'
    else:
        label = 'NO_TURN_IN_FAMILY'
    last_fail = next((c for c in reversed(rec['circles']) if not c['ok']), None)
    frac = w/(2*np.pi)
    lo, hi = float(frac.min()), float(theta0/(2*np.pi))
    rationals = sorted({Fraction(p, q) for q in range(1, Q_MAX+1) for p in range(0, q+1) if lo <= p/q <= hi})
    low_order = [str(f) for f in (Fraction(0), Fraction(1, 4), Fraction(1, 3), Fraction(1, 2), Fraction(2, 3), Fraction(3, 4), Fraction(1))]
    return dict(label=label, status=rec['status'], n_circles=len(good), a_max=float(max(c['a'] for c in good)),
                I_max=float(I.max()), frac_range=[lo, hi], max_omega_increase=float(steps.max()) if len(steps) else None,
                end_reason=last_fail.get('error') if last_fail else None,
                closure_rationals_q_le_12=[str(f) for f in rationals],
                low_order_closures_in_range=[s for s in low_order if lo <= float(Fraction(s)) <= hi])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('part', choices=['A', 'B', 'score'])
    args = ap.parse_args()
    theta0 = float(np.arccos(np.trace(fl.monodromy('T', 2))/2))
    src = sources()
    if args.part in ('A', 'B'):
        rec = measure_b() if args.part == 'B' else measure_a()
        if sources() != src:
            raise RuntimeError('sources changed during the run')
        _write(RUN_DIR/f'part_{args.part}.json', dict(spec='docs/r3_extension_prereg.md', sources=src, **rec))
        return
    A = json.loads((RUN_DIR/'part_A.json').read_text())
    B = json.loads((RUN_DIR/'part_B.json').read_text())
    if A['sources'] != src or B['sources'] != src:
        raise RuntimeError('archive sources do not match the committed code')
    result = dict(theta0=theta0, A=score_a(A, theta0), B=score_b(B, theta0),
                  C='NOT_TESTED: inhomogeneous n>=3 tensor, vector and scalar sectors (specification section 5)')
    _write(RUN_DIR/'result.json', result)
    print(json.dumps(result, indent=1))


if __name__ == '__main__':
    main()
