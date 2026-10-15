"""Post-hoc diagnostics of the breaking scan (written after the registered labels were fixed).

Is the coherent harmonic-10 structure of lambda(phi) in the main scan a property
of the flow, or a systematic error of the map evaluation? Any smooth map error
summed over five nodes spaced 2 pi/5 apart also contributes only harmonics that
are multiples of 5, so lambda is re-solved (chord Newton, archived Jacobian) with:
  T   DOP853 at tightened tolerances 1e-13/1e-15, all 60 phases
  R   fixed-step RK4 on the reduced state with a Newton-solved return time,
      Richardson-extrapolated over (2048, 4096); no solve_ivp, no event location.
      12 phases (every fifth).
Output: runs/20261005_r3_breaking/posthoc.json. Not registered; descriptive only.
"""
import json
from multiprocessing import Pool
import numpy as np
from geometrodynamics.waves import r3_breaking as b
from geometrodynamics.waves import r3_extension as rx
from geometrodynamics.waves import r3_family as f
from geometrodynamics.waves import r3_return_map as rm
from experiments.closure_ledger import r3_breaking_probe as probe

LEAD = 16


def P_tight(z):
    return rm.esu_map(np.asarray(z, float), tol=(1e-13, 1e-15))[0]


def _pack(u):
    return [u[0], u[1], u[2], u[3], [list(u[4:7]), list(u[7:10]), list(u[10:13])],
            [list(u[13:16]), list(u[16:19]), list(u[19:22])]]


def _unpack(y):
    return np.array([y[0], y[1], y[2], y[3]]+[e for row in y[4] for e in row]+[e for row in y[5] for e in row], float)


def _rk4_return(r0, T, steps):
    h = T/steps
    y = rx.rk4(_pack(r0), h, steps-LEAD)
    s = LEAD*h
    for _ in range(6):
        yc = rx.rk4(y, s/LEAD, LEAD)
        s = s-yc[2]/yc[3]
    return _unpack(rx.rk4(y, s/LEAD, LEAD))


def P_rk4(z):
    z = np.asarray(z, float)
    zf = np.r_[z, np.zeros(8)]
    T = rm.esu_map(z)[2]
    r0 = f._reduced(f.to_state(zf))
    a, c = _rk4_return(r0, T, 2048), _rk4_return(r0, T, 4096)
    return f.to_section(f._full(c+(c-a)/15))[:4]


def _one(args):
    kind, j = args
    pt = probe._read('scan_main.json')['points'][j]
    r = b.chord_resolve(P_tight if kind == 'T' else P_rk4, pt, tol=1e-13, maxit=4)
    return dict(kind=kind, j=j, lam=pt['lam'], lam_alt=r['lam'], residual=r['residual'])


def main():
    jobs = [('T', j) for j in range(probe.N_PHASE)]+[('R', j) for j in range(0, probe.N_PHASE, 5)]
    with Pool(4) as pool:
        rows = pool.map(_one, jobs)
    out = {}
    for kind in ('T', 'R'):
        rr = [r for r in rows if r['kind'] == kind]
        lam, alt = np.array([r['lam'] for r in rr]), np.array([r['lam_alt'] for r in rr])
        out[kind] = dict(rows=rr, max_abs_lam=float(np.abs(alt).max()), max_diff=float(np.abs(alt-lam).max()),
                         max_residual=float(max(r['residual'] for r in rr)),
                         correlation=float(np.corrcoef(lam, alt)[0, 1]))
    lam = np.array([r['lam_alt'] for r in rows if r['kind'] == 'T'])
    out['T']['harmonics'] = [float(x) for x in np.abs(np.fft.rfft(lam))[:16]/len(lam)]
    probe.RUN_DIR.mkdir(parents=True, exist_ok=True)
    (probe.RUN_DIR/'posthoc.json').write_text(json.dumps(dict(sources=probe.sources(), **out), indent=1))
    print(json.dumps({k: {kk: v for kk, v in out[k].items() if kk != 'rows'} for k in out}, indent=1))


if __name__ == '__main__':
    main()
