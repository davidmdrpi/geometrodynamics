"""Two-return family: closed-loop action and transverse stability (docs/r3_family_prereg.md).

stage F: trace the family of two-return solutions (P^2(z)=z, diagonal subsystem)
         by pseudo-arclength continuation around its full loop; loop action;
         bracketing circles for the interpolated resonant action.
stage S: at 12 loop points, exact-jet DP^2 in all 12 homogeneous section
         coordinates: family-tangent check, diagonal transverse pair,
         off-diagonal pairs, direct-perturbation and Radau checks.
score  : registered labels.
"""
import argparse
import hashlib
import json
import time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from geometrodynamics.waves import r3_family as rf

ROOT = Path(__file__).resolve().parents[2]
RUN_DIR = ROOT/'experiments/closure_ledger/runs/20261001_r3_family'
SOURCES = ('geometrodynamics/waves/r3_family.py', 'geometrodynamics/waves/r3_extension.py',
           'geometrodynamics/waves/jets.py', 'geometrodynamics/waves/nonlinear_supported_tt.py',
           'experiments/closure_ledger/r3_family_probe.py')
# #318 seed-0 two-return nodes (published in docs/diagonal_bianchi.md), used only as a Newton start
START = np.array([[1.000503231205612, -0.000025712699647, 0.077190977191945, 0.000802940697183,
                   0.000364201695803, -0.254452104496236],
                  [1.000436620468819, 0.000025726421217, -0.083326512003473, -0.001163768798227,
                   -0.000312623403504, 0.211950511588411]])
DS, MAX_STEPS, N_SAMPLE = .02, 400, 12
CIRCLE_A = (.06400, .07610925536, .09050966799)


def sources():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def _P6(z):
    return rf.P(np.r_[z, np.zeros(6)])[:6]


def _fdcol(args):
    z, k, h = args
    e = np.zeros(6)
    e[k] = h
    return (_P6(z+e)-_P6(z-e))/(2*h)


def jac2(v, pool, h=1e-6):
    cols = pool.map(_fdcol, [(v[:6], k, h) for k in range(6)]+[(v[6:], k, h) for k in range(6)])
    J0, J1 = np.array(cols[:6]).T, np.array(cols[6:]).T
    return np.block([[J0, -np.eye(6)], [-np.eye(6), J1]])


def resid(v, pool):
    p = pool.map(_P6, [v[:6], v[6:]])
    return np.r_[p[0]-v[6:], p[1]-v[:6]]


def corrector(v, pool, t=None, vp=None, tol=1e-12, maxit=10):
    for it in range(maxit):
        F = resid(v, pool)
        if t is not None:
            F = np.r_[F, t @ (v-vp)]
        if np.abs(F).max() < tol:
            break
        J = jac2(v, pool)
        if t is not None:
            J = np.vstack([J, t])
        v = v+np.linalg.lstsq(J, -F, rcond=None)[0]
    return v, float(np.abs(resid(v, pool)).max()), it+1


def tangent(v, pool, prev=None):
    _, s, Vt = np.linalg.svd(jac2(v, pool))
    t = Vt[-1]
    if prev is not None and t @ prev < 0:
        t = -t
    return t, float(s[-1]), float(s[-2])


def stage_f():
    rec = dict(points=[], stop=None)
    with Pool(4) as pool:
        v, r, it = corrector(START.ravel(), pool)
        v_start = v.copy()
        rec['start'] = dict(v=v.tolist(), residual=r, iterations=it)
        t, s1, s2 = tangent(v, pool)
        travelled = 0.
        for k in range(MAX_STEPS):
            vp = v+DS*t
            v_new, r, it = corrector(vp, pool, t, vp)
            travelled += np.linalg.norm(v_new-v)
            v = v_new
            t, s1, s2 = tangent(v, pool, t)
            rec['points'].append(dict(v=v.tolist(), residual=r, iterations=it, smin=s1, s2=s2))
            if travelled > 10*DS and np.linalg.norm(v-v_start) < 1.5*DS:
                # land exactly on the start along the start's own tangent row
                t0 = tangent(v_start, pool)[0]
                v_end, r_end, _ = corrector(v, pool, t0, v_start)
                rec['stop'] = dict(reason='RETURNED_TO_START', landing_distance=float(np.linalg.norm(v_end-v_start)),
                                   landing_residual=r_end, steps=k+1)
                break
        if rec['stop'] is None:
            rec['stop'] = dict(reason='MAX_STEPS', steps=MAX_STEPS)
        # bracketing circles (independent of the loop), 63 nodes
        rec['circles'] = circles(pool)
    return rec


def circles(pool):
    from geometrodynamics.waves.r3_return_map import shift_matrix
    M = 63
    J = rf._fd(_P6, np.r_[1., np.zeros(5)], 1e-6)
    lam, vec = np.linalg.eig(J[2:4, 2:4])
    i = int(np.argmax(lam.imag))
    q = vec[:, i]/vec[0, i]
    th = 2*np.pi*np.arange(M)/M
    base = np.r_[1., np.zeros(5)]
    out = []
    K, om, a_prev = None, float(np.angle(lam[i])), None
    ladder = [.004*2**(k/4) for k in range(18) if .004*2**(k/4) < CIRCLE_A[0]]+list(CIRCLE_A)
    for a in ladder:
        if K is None:
            K = base+np.c_[np.zeros((M, 2)), a*np.real(q[None]*np.exp(1j*th)[:, None]),
                           a*np.real(1j*q[None]*np.exp(1j*th)[:, None])]
        else:
            K = base+(K-base)*(a/a_prev)
        e1 = np.exp(-1j*th)/M
        for it in range(20):
            PK = np.array(pool.map(_P6, list(K)))
            T, dT = shift_matrix(M, om)
            R = PK-T @ K
            c1 = e1 @ K[:, 2]
            res = np.r_[R.ravel(), c1.real-a/2, c1.imag]
            if np.abs(res).max() < 5e-11:
                break
            Jb = pool.map(_fdpt, list(K))
            Jm = np.zeros((6*M+2, 6*M+1))
            for j in range(M):
                Jm[6*j:6*j+6, 6*j:6*j+6] = Jb[j]
            Jm[:6*M, :6*M] -= np.kron(T, np.eye(6))
            Jm[:6*M, 6*M] = -(dT @ K).ravel()
            Jm[6*M, 2:6*M:6] = e1.real
            Jm[6*M+1, 2:6*M:6] = e1.imag
            st = np.linalg.lstsq(Jm, -res, rcond=None)[0]
            K = K+st[:6*M].reshape(M, 6)
            om += st[6*M]
        a_prev = a
        if a in CIRCLE_A:
            out.append(dict(a=a, omega=float(om), action=rf.loop_action(K), residual=float(np.abs(res).max())))
    return out


def _fdpt(z):
    return rf._fd(_P6, z, 1e-6)


def stage_s(F):
    pts = [np.array(F['start']['v'])]+[np.array(p['v']) for p in F['points']]
    Z = np.array([p[:6] for p in pts])
    s = np.r_[0., np.cumsum(np.linalg.norm(np.diff(Z, axis=0), axis=1))]
    picks = [int(np.argmin(abs(s-x))) for x in np.linspace(0, s[-1], N_SAMPLE, endpoint=False)]
    out = []
    with Pool(4) as pool:
        for j in picks:
            v = pts[j]
            z0, z1 = v[:6], v[6:]
            p0, J0 = rf.DP(np.r_[z0, np.zeros(6)], 12)
            p1, J1 = rf.DP(np.r_[z1, np.zeros(6)], 12)
            M2 = J1 @ J0
            Jn = np.block([[J0[:6, :6], -np.eye(6)], [-np.eye(6), J1[:6, :6]]])      # exact-jet two-node Jacobian
            _, sv, Vt = np.linalg.svd(Jn)
            t = Vt[-1][:6]/np.linalg.norm(Vt[-1][:6])
            radau = max(np.abs(rf.P(np.r_[rf.P(np.r_[z0, np.zeros(6)], method='Radau')[:6], np.zeros(6)],
                                    method='Radau')[:6]-z0))
            out.append(dict(index=j, v=v.tolist(), M2=M2.tolist(), closure_jet=float(max(np.abs(p0[:6]-z1).max(), np.abs(p1[:6]-z0).max())),
                            radau_closure=float(radau), tangent_defect=float(np.linalg.norm((M2[:6, :6]-np.eye(6)) @ t)),
                            smin_exact=float(sv[-1]), s2_exact=float(sv[-2])))
        # direct-perturbation check at the first sample, along each off-diagonal eigenvector (real part)
        v = pts[picks[0]]
        M2 = np.array(out[0]['M2'])
        lam, V = np.linalg.eig(M2[6:, 6:])
        direct = []
        for k in range(6):
            dv = np.r_[np.zeros(6), np.real(V[:, k])]
            dv /= np.linalg.norm(dv)
            eps = 1e-6
            z = np.r_[v[:6], np.zeros(6)]
            P2 = lambda w: rf.P(rf.P(w))
            lin = (P2(z+eps*dv)-P2(z-eps*dv))/(2*eps)
            direct.append(float(np.linalg.norm(lin-M2 @ dv)/np.linalg.norm(M2 @ dv)))
    return dict(samples=out, direct_perturbation=direct)


def score(F, S):
    pts = [F['start']]+F['points']
    res_max = max(p['residual'] for p in pts)
    samples = S['samples']
    Z = np.array([np.array(p['v'])[:6] for p in pts])
    I_all, I_half = rf.loop_action(Z), rf.loop_action(Z[::2])
    circ = F['circles']
    Ic, wc = np.array([c['action'] for c in circ]), np.array([c['omega'] for c in circ])
    quad, lin = np.polyfit(wc, Ic, 2), np.polyfit(wc[1:], Ic[1:], 1)
    I_star_q, I_star_l = float(np.polyval(quad, np.pi)), float(np.polyval(lin, np.pi))
    diag, off, esu, decouple, dets = [], [], [], 0., 0.
    for smp in samples:
        M2 = np.array(smp['M2'])
        decouple = max(decouple, np.abs(M2[:6, 6:]).max()/np.abs(M2).max(), np.abs(M2[6:, :6]).max()/np.abs(M2).max())
        dets = max(dets, abs(np.linalg.det(M2)-1))
        pd = rf.classify_pairs(np.linalg.eigvals(M2[:6, :6]))
        e = max(pd, key=lambda p: abs(p['trace']))                  # Einstein-static pair
        rest = [p for p in pd if p is not e]
        unit = min(rest, key=lambda p: abs(p['trace']-2))            # family (unit) pair
        esu.append(e)
        diag.append([p for p in rest if p is not unit][0])          # diagonal transverse pair
        off.append(rf.classify_pairs(np.linalg.eigvals(M2[6:, 6:])))
    checks = dict(
        F1_loop_residuals=res_max <= 1e-10,
        F2_loop_closed=F['stop']['reason'] == 'RETURNED_TO_START' and F['stop']['landing_distance'] <= 1e-8,
        F3_tangent_is_unit_eigenvector=max(s['tangent_defect'] for s in samples) <= 1e-6,
        F4_radau_closure=max(s['radau_closure'] for s in samples) <= 1e-9,
        F5_action_converged=abs(I_all-I_half) <= 1e-6*I_all,
        S1_block_decoupling=decouple <= 1e-8,
        S2_determinant=dets <= 1e-6,
        S3_direct_perturbation=max(S['direct_perturbation']) <= 1e-3,
        S4_jet_closure=max(s['closure_jet'] for s in samples) <= 1e-9,
    )
    fam_ok = all(checks[k] for k in checks if k.startswith('F'))
    stab_ok = all(checks[k] for k in checks if k.startswith('S'))
    kinds_d = {p['kind'] for p in diag}
    kinds_o = {p['kind'] for pairs in off for p in pairs}
    return dict(
        checks=checks,
        FAMILY='CLOSED_FAMILY_LOOP_NUMERICALLY' if fam_ok else 'FAMILY_UNRESOLVED',
        loop_action=I_all, loop_action_uncertainty=abs(I_all-I_half),
        circle_interpolated_action=dict(quadratic=I_star_q, linear=I_star_l),
        ACTION_CONSISTENCY='CONSISTENT' if abs(I_all-I_star_q) <= 3*(abs(I_star_q-I_star_l)+abs(I_all-I_half)) else 'INCONSISTENT',
        DIAGONAL_TRANSVERSE=('UNRESOLVED' if not stab_ok else 'ELLIPTIC' if kinds_d == {'ELLIPTIC'} else
                             'HYPERBOLIC_PRESENT' if 'HYPERBOLIC' in kinds_d else 'MARGINAL'),
        OFF_DIAGONAL=('UNRESOLVED' if not stab_ok else 'ALL_ELLIPTIC' if kinds_o == {'ELLIPTIC'} else
                      'HYPERBOLIC_PRESENT' if 'HYPERBOLIC' in kinds_o else 'MARGINAL'),
        esu_pairs=[p['trace'] for p in esu], diagonal_pairs=diag, off_diagonal_pairs=off,
        loop_points=len(pts), max_loop_residual=res_max, landing=F['stop'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('stage', choices=['F', 'S', 'score'])
    args = ap.parse_args()
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    src = sources()
    path = RUN_DIR/f'stage_{args.stage}.json' if args.stage != 'score' else RUN_DIR/'result.json'
    if path.exists():
        raise FileExistsError('append-only: '+str(path))
    if args.stage == 'F':
        rec = stage_f()
    elif args.stage == 'S':
        rec = stage_s(json.loads((RUN_DIR/'stage_F.json').read_text()))
    else:
        F = json.loads((RUN_DIR/'stage_F.json').read_text())
        S = json.loads((RUN_DIR/'stage_S.json').read_text())
        if F['sources'] != src or S['sources'] != src:
            raise RuntimeError('sources changed')
        rec = score(F, S)
        print(json.dumps({k: v for k, v in rec.items() if k not in ('diagonal_pairs', 'off_diagonal_pairs')}, indent=1))
    if sources() != src:
        raise RuntimeError('sources changed during the run')
    path.write_text(json.dumps(dict(sources=src, **rec), indent=1, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
