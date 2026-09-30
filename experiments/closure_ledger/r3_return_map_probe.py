"""Constrained R3 return map: leading frequency-shift coefficient nu.

Specification: docs/r3_return_map_prereg.md. Method 1 (jets + Birkhoff normal
form) and Method 2 (invariant circles) are independent evaluations of the same
clock-section map. Writes an archive binding its sources.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from geometrodynamics.waves import esu_floquet as fl
from geometrodynamics.waves import r3_return_map as rm

ROOT = Path(__file__).resolve().parents[2]
RUN_DIR = ROOT/'experiments/closure_ledger/runs/20260929_r3_return_map'
SOURCES = ('geometrodynamics/waves/jets.py', 'geometrodynamics/waves/r3_return_map.py',
           'geometrodynamics/waves/nonlinear_supported_tt.py', 'geometrodynamics/waves/esu_floquet.py',
           'experiments/closure_ledger/r3_return_map_probe.py')
STEPS = (1024, 2048, 4096)
LADDER = (.004, .008, .016, .032)
GRID = 31
RADAU_CHECK = .008


def sources():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def theta0_310():
    return float(np.arccos(np.trace(fl.monodromy('T', 2))/2))


def measure_m1():
    runs = {n: rm.jet_return_map(3, n) for n in STEPS}
    R12 = rm.richardson(runs[1024]['P'], runs[2048]['P'])
    R24 = rm.richardson(runs[2048]['P'], runs[4096]['P'])
    nf, nf12 = rm.normal_form(R24), rm.normal_form(R12)
    scale = max(float(abs(p.c).max()) for p in R24)
    ablation = rm.normal_form(rm.jet_return_map(3, 4096, fixed_time=True)['P'])
    tjet = runs[4096]['time']
    return dict(normal_form=nf, normal_form_coarse=nf12, coef_scale=scale,
                symplectic_defect=rm.symplectic_defect(R24),
                constraint_max=float(abs(runs[4096]['constraint'].c).max()),
                section_q_max=float(abs(runs[4096]['q'].c).max()),
                return_time_quadratic={str(tuple(int(v) for v in e)): float(tjet.c[i])
                                       for i, e in enumerate(tjet.s['exps']) if tjet.s['deg'][i] == 2},
                ablation_fixed_time_nu=ablation['nu'],
                coefficients=[p.c.tolist() for p in R24])


def measure_m2():
    P = lambda z: rm.esu_map(z)[0]
    fixed, fixed_res, _ = rm.esu_map(rm.ZSTAR)
    J = rm.fd_jacobian(P, rm.ZSTAR)
    K, omega = rm.linear_circle(J, LADDER[0], GRID)
    circles = []
    for a in LADDER:
        K = rm.ZSTAR+(K-rm.ZSTAR)*(a/(2*abs(np.fft.fft(K[:, 2])[1])/GRID))
        t = time.time()
        c = rm.invariant_circle(P, a, K, omega, GRID)
        K, omega = c['K'], c['omega']
        cres = max(rm.esu_map(k)[1] for k in K)
        rec = dict(a=a, omega=c['omega'], action=abs(c['action']), action_signed=c['action'],
                   residual=c['residual'], iterations=c['iterations'], history=c['history'],
                   fourier_tail=c['fourier_tail'], constraint_max=cres, K=K.tolist(), seconds=time.time()-t)
        if a == RADAU_CHECK:
            T, _ = rm.shift_matrix(GRID, omega)
            PR = np.array([rm.esu_map(k, method='Radau')[0] for k in K])
            rec['radau_invariance_residual'] = float(np.max(abs(PR-T @ K)))
        circles.append(rec)
    return dict(fixed_point_error=float(np.max(abs(fixed-rm.ZSTAR))), fixed_constraint=fixed_res,
                fd_jacobian=J.tolist(), fd_symplectic_defect=float(np.max(abs(J.T@rm.OMEGA@J-rm.OMEGA))),
                fd_eigenvalues=[[float(v.real), float(v.imag)] for v in np.linalg.eigvals(J)], circles=circles)


def fit_m2(circles, theta0):
    I = np.array([c['action'] for c in circles])
    w = np.array([c['omega'] for c in circles])
    y = (w-theta0)/(2*np.pi*I)
    quad = np.polyfit(I, y, 2)
    lin = np.polyfit(I[:3], y[:3], 1)
    return float(quad[-1]), float(abs(quad[-1]-lin[-1])), y.tolist()


def score(m1, m2, theta0):
    nf, fails = m1['normal_form'], []
    nu1 = nf['nu']
    conv = abs(nu1-m1['normal_form_coarse']['nu'])
    chk = dict(
        C1_fixed_point=nf['fixed_error'] <= 1e-10 and m2['fixed_point_error'] <= 1e-10,
        C2_linear=abs(nf['theta']-theta0) <= 1e-9 and abs(np.angle(complex(*max(m2['fd_eigenvalues'], key=lambda v: v[1])))-theta0) <= 1e-6,
        C3_symplectic=m1['symplectic_defect'] <= 1e-10*max(1., m1['coef_scale']) and m2['fd_symplectic_defect'] <= 1e-4,
        C4_dissipative=abs(nf['dissipative']) <= 1e-5*abs(nf['twist']),
        C5_constraint=m1['constraint_max'] <= 1e-8*max(1., m1['coef_scale']) and m1['section_q_max'] <= 1e-8*max(1., m1['coef_scale'])
        and m2['fixed_constraint'] <= 1e-10 and all(c['constraint_max'] <= 1e-10 for c in m2['circles']),
        C6_nonresonance=min(nf['resonance_distances']) >= .1 and nf['min_quadratic_divisor'] >= .1,
        C7_m1_convergence=conv <= 1e-6*abs(nu1)+1e-9,
        D1_invariance=all(c['residual'] <= 1e-10 for c in m2['circles']),
        D2_fourier=all(c['fourier_tail'] <= 1e-10 for c in m2['circles']),
        D3_radau=any('radau_invariance_residual' in c for c in m2['circles'])
        and all(c.get('radau_invariance_residual', 0) <= 1e-9 for c in m2['circles']),
    )
    nu2, unc2, ratios = fit_m2(m2['circles'], theta0)
    chk['A1_agreement'] = unc2 <= 1e-3*abs(nu1) and abs(nu1-nu2) <= max(1e-4*abs(nu1), 3*unc2)
    I = np.array([c['action'] for c in m2['circles']])
    rem = np.array([c['omega']-theta0-2*np.pi*nu1*c['action'] for c in m2['circles']])
    slopes = [float(np.log(abs(rem[k+1]/rem[k]))/np.log(I[k+1]/I[k])) for k in range(len(I)-1)
              if min(abs(rem[k]), abs(rem[k+1])) >= 1e-8]
    chk['A2_remainder'] = bool(slopes) and all(s >= 1.8 for s in slopes)
    fails = [k for k, v in chk.items() if not v]
    resolved = abs(nu1) >= 10*max(conv, abs(nu1-nu2), unc2)
    if fails or not resolved:
        label = 'UNRESOLVED'
    else:
        label = 'SHIFT_TOWARD_TARGET' if nu1 > 0 else 'SHIFT_AWAY_FROM_TARGET'
    return dict(label=label, checks={k: bool(v) for k, v in chk.items()}, failures=fails,
                nu_m1=nu1, nu_m1_convergence=conv, nu_m2=nu2, nu_m2_uncertainty=unc2,
                m2_ratios=ratios, actions=I.tolist(), remainders=rem.tolist(), remainder_slopes=slopes,
                theta0=theta0, frac_rho0=theta0/(2*np.pi), target_frac=.5,
                ablation_fixed_time_nu=m1['ablation_fixed_time_nu'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, default=RUN_DIR/'return_map.json')
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError('archives are append-only; choose a new path')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    src = sources()
    theta0 = theta0_310()
    m1 = measure_m1()
    m2 = measure_m2()
    if sources() != src:
        raise RuntimeError('sources changed during the run')
    result = score(m1, m2, theta0)
    blob = json.dumps(dict(spec='docs/r3_return_map_prereg.md', sources=src, m1=m1, m2=m2, result=result),
                      indent=1, allow_nan=False, default=lambda o: [o.real, o.imag] if isinstance(o, complex) else str(o))
    args.output.write_text(blob+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('m2_ratios',)}, indent=1))


if __name__ == '__main__':
    main()
