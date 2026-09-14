"""Classical operator controls for the public field-to-apparatus freeze.

The graph, intrinsic circle and imposed potential have DIFFERENT actions.
No function here constructs a bulk-to-throat interface or event frequencies.
"""
from functools import lru_cache
import copy
import itertools
import numpy as np
from . import nonlinear_supported_tt as bulk

BASELINE = '4cd86541d3b836b35561b0c4a3a54629d28851cd'
PREREG = 'e528c847d172a4629dc1bcd737fc46979b259a2f'
SEED = 2026091401
GRIDS = (32, 64, 128)
PHASES = tuple(float(j*np.pi/4) for j in range(8))
MODES = np.array([1, -1, 2, -2])


def encode(z):
    z = np.asarray(z)
    return np.stack((z.real, z.imag), axis=-1).tolist()


def decode(z):
    z = np.asarray(z, dtype=float)
    if z.shape[-1:] != (2,) or not np.isfinite(z).all():
        raise ValueError('finite complex pairs required')
    return z[..., 0] + 1j*z[..., 1]


def relative(x, y):
    x, y = np.asarray(x), np.asarray(y)
    if x.shape != y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():
        return float('inf')
    return float(np.linalg.norm(x-y)/max(1., np.linalg.norm(x), np.linalg.norm(y)))


def grid(N):
    return 2*np.pi*np.arange(N)/N


def graph(epsilon, phase, m, N=8, derivative=False):
    """Independent incidence assembly, full winding space and wrap retained."""
    chi = grid(N) + np.pi/N
    f = np.cos(m*chi + phase)
    if np.min(1+epsilon*f) <= 0:
        raise ValueError('positive link radii required')
    weights = -2*f if derivative else (1+epsilon*f)**-2
    incidence = np.zeros((N, N))
    for j in range(N):
        incidence[j, j] = -1
        incidence[j, (j+1) % N] = 1
    return incidence.T @ (weights[:, None]*incidence)


def graph_block(H):
    chi = grid(len(H))
    modes = np.exp(1j*chi[:, None]*np.array([1, -1]))/np.sqrt(len(H))
    return modes.conj().T @ H @ modes


def quartet(q, A, circle, N):
    chi = grid(N)
    x = np.zeros((N, 4))
    x[:, circle[0]], x[:, circle[1]] = np.cos(chi), np.sin(chi)
    fields = x @ bulk.B(q).T/A
    intensity = np.sum(fields**2, axis=1)
    return dict(mean=float(np.mean(intensity)),
                fourier=encode(np.mean(intensity[:, None]*np.exp(-1j*chi[:, None]*np.arange(1, 4)), axis=0)))


def circle(R0, epsilon, m, phase, N):
    """Two independent coordinate expressions for the arclength eigenmodes.

    The coordinate flux derivative retains R' terms until cancellation.
    K1 and W1 are differentiated stiffness and kinetic forms in a fixed,
    W0-normalized +/-1 basis; they are individually nonzero at m=2.
    """
    chi = grid(N)
    f = np.cos(m*chi+phase)
    R = R0*(1+epsilon*f)
    Rp = -R0*epsilon*m*np.sin(m*chi+phase)
    if R0 <= 0 or np.min(R) <= 0:
        raise ValueError('positive metric required')
    s = R0*(chi+epsilon*(np.sin(m*chi+phase)-np.sin(phase))/m)
    lam = MODES**2/R0**2
    u = np.exp(1j*MODES[:, None]*s/R0)
    up = 1j*MODES[:, None]*R*u/R0
    upp = (1j*MODES[:, None]*Rp/R0-MODES[:, None]**2*R**2/R0**2)*u
    lhs = Rp*up/R**2-upp/R
    rhs = lam[:, None]*R*u
    # Each mode normalized with the true kinetic measure R dchi.
    gram = np.einsum('in,jn,n->ij', u.conj(), u, R)/(N*R0)
    ks = np.array([1, -1])
    modes = np.exp(1j*ks[:, None]*chi)
    W1 = np.einsum('in,jn,n->ij', modes.conj(), modes, f)/N
    K1 = -ks[:, None]*ks[None, :]*W1/R0**2
    return dict(lhs=encode(lhs), rhs=encode(rhs), kinetic_gram=encode(gram),
                K1=encode(K1), W1=encode(W1),
                circumference=float(2*np.pi*np.mean(R)))


def potential(v, m, phase, N):
    chi = grid(N)
    modes = np.exp(1j*np.array([1, -1])[:, None]*chi)
    V = v*np.cos(m*chi+phase)
    return np.einsum('in,jn,n->ij', modes.conj(), modes, V)/N


@lru_cache(None)
def _certificate():
    import sympy as s
    q = s.symbols('q0:4', real=True)
    # Construct the symbolic matrix from the inherited, integer S_i convention.
    matrices = [s.eye(4)] + [s.Matrix(v.astype(int).tolist()) for v in bulk.S]
    B = sum((q[i]*matrices[i] for i in range(4)), s.zeros(4))
    gram = (B.T*B-s.eye(4)*sum(v*v for v in q)).applyfunc(s.expand)
    x = s.symbols('x0:4', real=True)
    intensity = s.expand(((B*s.Matrix(x)).T*(B*s.Matrix(x)))[0]
                         -sum(v*v for v in q)*sum(v*v for v in x))
    c, d, r, z = s.symbols('c d r z', nonzero=True, real=True)
    blocks = []
    for m in (1, 2, 3):
        f = ((c+s.I*d)*z**m+(c-s.I*d)*z**(-m))/2
        W = s.Matrix(2, 2, lambda i, j: s.expand(f).coeff(z, (1, -1)[i]-(1, -1)[j]))
        K = s.Matrix(2, 2, lambda i, j: -(1, -1)[i]*(1, -1)[j]*W[i, j]/r**2)
        residual = (K-W/r**2).applyfunc(s.simplify)
        blocks.append(dict(m=m, W1=[str(v) for v in W], K1=[str(v) for v in K],
                           generalized=[str(v) for v in residual]))
    # Fourier potential selection is derived from the Laurent polynomial,
    # not from a requested coefficient passed to a numerical routine.
    selections = {str(m):str(s.expand(((c+s.I*d)*z**m+(c-s.I*d)*z**(-m))/2).coeff(z, 2)) for m in (1, 2, 3)}
    t, e, R0, k, phase = s.symbols('chi epsilon R0 k phase', real=True)
    mapped = []
    for m in (1, 2, 3):
        R = R0*(1+e*s.cos(m*t+phase))
        arc = R0*(t+e*(s.sin(m*t+phase)-s.sin(phase))/m)
        u = s.exp(s.I*k*arc/R0)
        residual = s.simplify(-s.diff(s.diff(u,t)/R,t)-k*k/R0**2*R*u)
        mapped.append(str(residual))
    return dict(quaternion_gram=[str(v) for v in gram], intensity=str(intensity),
                circle_blocks=blocks, mapped_equation=mapped, potential_selection=selections)


def certificate():
    return copy.deepcopy(_certificate())


def quartet_cases():
    qs = np.random.default_rng(SEED).normal(size=(64, 4))
    for i, q in enumerate(qs):
        for A, ij, N in itertools.product((1., 2.), itertools.combinations(range(4), 2), GRIDS):
            yield dict(draw=i, q=q.tolist(), A=A, circle=list(ij), N=N)


def circle_cases():
    for R0, epsilon, m, phase, N in itertools.product((1.,2.), (0.,-.05,.05,-.1,.1), (1,2,3), (0.,np.pi/4,np.pi/2), GRIDS):
        yield dict(R0=R0, epsilon=epsilon, m=m, phase=phase, N=N)
