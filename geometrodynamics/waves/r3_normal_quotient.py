"""Linear normal quotient of the #319 family; no nonlinear reduction is claimed."""
import numpy as np
from scipy.linalg import schur
from .r3_extension import E

POWERS = 2**np.arange(11)


def finite(x, shape):
    x = np.asarray(x, float)
    if x.shape != shape or not np.isfinite(x).all():
        raise ValueError(f'expected finite {shape}')
    return x


def rotation_generators(z):
    z = finite(z, (12,))
    beta = np.einsum('k,kij->ij', z[2::2], E)
    momentum = np.einsum('k,kij->ij', z[3::2], E)
    out = np.zeros((12, 3))
    for k, (i, j) in enumerate(((0, 1), (0, 2), (1, 2))):
        O = np.zeros((3, 3)); O[i, j] = 1.; O[j, i] = -1.
        out[2::2, k] = np.einsum('kij,ij->k', E, O@beta-beta@O)
        out[3::2, k] = np.einsum('kij,ij->k', E, O@momentum-momentum@O)
    return out


def quotient(M, generators):
    """Euclidean quotient; invariance is returned as a gate, not assumed."""
    n, k = generators.shape
    U, s, _ = np.linalg.svd(generators, full_matrices=True)
    if s[-1] <= 1e-6:
        raise ValueError('dependent removed directions')
    R, Q = U[:, :k], U[:, k:]
    return Q.T@M@Q, Q, R, float(np.linalg.norm(Q.T@M@R, 2))


def modal_check(C, neutral=4):
    """Numerical similarity to neutral identity plus a planar rotation."""
    n = neutral+2
    C = finite(C, (n, n))
    _, sv, Vt = np.linalg.svd(C-np.eye(n))
    ev, V = np.linalg.eig(C)
    candidates = np.flatnonzero(ev.imag > 1e-3)
    ans = dict(singular_values=sv.tolist(), eigenvalues=np.c_[ev.real, ev.imag].tolist(),
               neutral_dimension=int(np.sum(sv <= 1e-7)), modal_ok=False)
    if len(candidates) != 1 or ans['neutral_dimension'] != neutral or sv[1] <= 1e-3:
        return ans
    j = candidates[0]; v = V[:, j]
    W = np.column_stack((Vt[-neutral:].T, v.real, v.imag))
    cond = float(np.linalg.cond(W))
    if not np.isfinite(cond) or cond > 1e4:
        ans['condition'] = cond if np.isfinite(cond) else None
        return ans
    unit = ev[j]/abs(ev[j])
    D = np.eye(n); D[-2:, -2:] = [[unit.real, unit.imag], [-unit.imag, unit.real]]
    Wi = np.linalg.inv(W); G = Wi.T@Wi
    similarity = float(np.linalg.norm(C@W-W@D, 2)/np.linalg.norm(W, 2))
    metric = float(np.linalg.norm(C.T@G@C-G, 2)/np.linalg.norm(G, 2))
    norms = [float(np.linalg.norm(np.linalg.matrix_power(C, int(k)), 2)) for k in POWERS]
    err = float(abs(abs(ev[j])-1))
    ans.update(condition=cond, similarity_defect=similarity, metric_defect=metric,
               metric=G.tolist(), similarity_basis=W.tolist(), elliptic_trace=float(2*ev[j].real),
               elliptic_modulus_error=err, power_norms=norms,
               modal_ok=bool(err <= 1e-7 and similarity <= 1e-7 and metric <= 1e-7
                             and max(norms) <= 1.01*cond))
    return ans


def reduce_sample(M, z, chord):
    M = finite(M, (12, 12)); z = finite(z, (12,)); chord = finite(chord, (6,))
    if np.linalg.norm(chord) == 0:
        raise ValueError('zero family chord')
    _, fs, ft = np.linalg.svd(M[:6, :6]-np.eye(6))
    tangent = ft[-1]
    chord = chord/np.linalg.norm(chord)
    sine = float(np.sqrt(max(0., 1-min(1., abs(tangent@chord))**2)))
    rot = rotation_generators(z)
    rs = np.linalg.svd(rot, compute_uv=False)
    B, Q, R, leak = quotient(M, np.column_stack((rot, np.r_[tangent, np.zeros(6)])))
    T, Uall, count = schur(B, output='real', sort=lambda re, im: .1 < np.hypot(re, im) < 10.)
    U = Uall[:, :count]; C = U.T@B@U
    excluded = np.linalg.eigvals(T[count:, count:])
    centre_leak = float(np.linalg.norm(B@U-U@C, 2))
    checks = dict(Q1=bool(rs[-1] > 1e-6 and fs[-1] <= 1e-7 and sine <= 1e-2),
                  Q2=bool(leak <= 1e-7 and centre_leak <= 1e-7),
                  Q3=bool(count == 6 and len(excluded) == 2
                          and min(abs(excluded)) < .1 and max(abs(excluded)) > 10), Q4=False)
    modal = modal_check(C) if count == 6 else dict(modal_ok=False)
    checks['Q3'] &= modal.get('neutral_dimension') == 4
    checks['Q4'] = modal['modal_ok']
    return dict(Q=Q.tolist(), R=R.tolist(), B=B.tolist(), U=U.tolist(), C=C.tolist(),
                rotation_singular_values=rs.tolist(), family_singular_values=fs.tolist(),
                family_chord_sine=sine, removed_leakage=leak, centre_leakage=centre_leak,
                excluded_eigenvalues=np.c_[excluded.real, excluded.imag].tolist(), modal=modal, checks=checks)
