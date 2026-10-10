"""High-precision clock-section return map of the LRS sector by the Taylor-series method.

Same system as r3_return_map (kappa = a = 1, conformal time), written for
Taylor recurrences in the variables y = (A, V = A', q, W = q', x, P = H x'):

    A' = V,   V' = -A (r + ell)/6 + A^3,
    q' = W,   W' = -(trinv + (r + ell)/6) q,
    x' = X = P/H,   P' = (q^2 - 4H) fb_inv - 4 H fb_sq,

with H = A^2 - q^2/6, ell = X^2, and, writing E_c = exp(c x/sqrt 6),
    trinv = 2 E_-2 + E_4,  r = 2(2 trinv - tr M^2) = 8 E_-2 - 2 E_-8,
    fb_inv = (2 E_-2 - 2 E_4)/sqrt 6,  fb_sq = (2 E_4 - 2 E_-8)/sqrt 6.
P' follows from (H x')' = H x'' + H' x' and r3_return_map.rhs. The section is
q = 0, q' < 0, first crossing after t = 1 (as esu_map), with the constraint
solved for q' exactly as r3_return_map.section_to_state. Section coordinates
z = (A, p_A, x, p_x) = (A, -6V, x, P).

Arithmetic is gmpy2 mpfr at `bits` precision; each step takes the order-N
Taylor polynomial with step h = rho eps^(1/N), rho the radius estimated from
the last two coefficients; the crossing is found by Newton on the q polynomial
inside the step. Pure Python; deterministic given (bits, N, eps).
"""
import gmpy2
from gmpy2 import mpfr

C1 = dict(bits=160, order=40, eps=1e-40)
C2 = dict(bits=200, order=50, eps=1e-50)


def _conv(u, v, k):
    return gmpy2.fsum(u[j]*v[k-j] for j in range(k+1))


def _exp_next(x, E, c, k):
    """k-th coefficient of exp(c x) from the recurrence E' = c x' E."""
    return c*gmpy2.fsum(j*x[j]*E[k-j] for j in range(1, k+1))/k


def coefficients(y, N):
    """Normalized Taylor coefficients (orders 0..N) of the six state variables at the current point."""
    s6 = gmpy2.sqrt(mpfr(6))
    c2, c4, c8 = mpfr(-2)/s6, mpfr(4)/s6, mpfr(-8)/s6
    A, V, q, W, x, P = ([v] for v in y)
    E2, E4, E8 = [gmpy2.exp(c2*x[0])], [gmpy2.exp(c4*x[0])], [gmpy2.exp(c8*x[0])]
    Asq, Q, H, X, ell, S, G, fbi, fbs, QmH = [], [], [], [], [], [], [], [], [], []
    for k in range(N):
        if k:
            E2.append(_exp_next(x, E2, c2, k))
            E4.append(_exp_next(x, E4, c4, k))
            E8.append(_exp_next(x, E8, c8, k))
        Asq.append(_conv(A, A, k))
        Q.append(_conv(q, q, k))
        H.append(Asq[k]-Q[k]/6)
        X.append((P[k]-gmpy2.fsum(H[j]*X[k-j] for j in range(1, k+1)))/H[0])
        ell.append(_conv(X, X, k))
        S.append(8*E2[k]-2*E8[k]+ell[k])
        G.append(2*E2[k]+E4[k]+S[k]/6)
        fbi.append((2*E2[k]-2*E4[k])/s6)
        fbs.append((2*E4[k]-2*E8[k])/s6)
        QmH.append(Q[k]-4*H[k])
        Vp = -_conv(A, S, k)/6+_conv(Asq, A, k)
        Wp = -_conv(G, q, k)
        Pp = _conv(QmH, fbi, k)-4*_conv(H, fbs, k)
        A.append(V[k]/(k+1))
        V.append(Vp/(k+1))
        q.append(W[k]/(k+1))
        W.append(Wp/(k+1))
        x.append(X[k]/(k+1))
        P.append(Pp/(k+1))
    return [A, V, q, W, x, P]


def _horner(c, s):
    out = c[-1]
    for a in reversed(c[:-1]):
        out = out*s+a
    return out


def _dhorner(c, s):
    out = len(c[1:])*c[-1]
    for k in range(len(c)-2, 0, -1):
        out = out*s+k*c[k]
    return out


def _radius(cs, N):
    rho = None
    for c in cs:
        for k in (N-1, N):
            if c[k] != 0:
                r = abs(c[k])**(-mpfr(1)/k)
                rho = r if rho is None else min(rho, r)
    return rho


def constraint(y):
    A, V, q, W, x, P = y
    s6 = gmpy2.sqrt(mpfr(6))
    E2, E4, E8 = gmpy2.exp(-2*x/s6), gmpy2.exp(4*x/s6), gmpy2.exp(-8*x/s6)
    H = A*A-q*q/6
    trinv, r = 2*E2+E4, 8*E2-2*E8
    X = P/H
    return -3*V*V+W*W/2+H*X*X/2-H*r/2+q*q*trinv/2+3*A**4/2


def section_to_state(z):
    A, pA, x, px = (mpfr(v) for v in z)
    s6 = gmpy2.sqrt(mpfr(6))
    E2, E8 = gmpy2.exp(-2*x/s6), gmpy2.exp(-8*x/s6)
    H = A*A
    V, X = -pA/6, px/H
    r = 8*E2-2*E8
    w2 = 2*(3*V*V-H*X*X/2+H*r/2-3*A**4/2)
    if w2 <= 0:
        raise ValueError('no real clock velocity')
    return [A, V, mpfr(0), -gmpy2.sqrt(w2), x, px]


def flow_to_section(y, t_min=1, t_max=5, order=40, eps=1e-40, direction=-1):
    """Integrate from state y to the first q = 0 crossing with sign(q') = direction after t_min.
    Returns (state, time, steps)."""
    N = order
    e = mpfr(eps)**(mpfr(1)/N)
    t, steps = mpfr(0), 0
    while t < t_max:
        cs = coefficients(y, N)
        h = _radius(cs, N)*e
        q0, q1 = cs[2][0], _horner(cs[2], h)
        if t+h > t_min and direction*q0 < 0 and direction*q1 >= 0:
            s = h*q0/(q0-q1)
            for _ in range(200):
                ds = _horner(cs[2], s)/_dhorner(cs[2], s)
                s -= ds
                if abs(ds) <= abs(s)*mpfr(2)**(-gmpy2.get_context().precision+8):
                    break
            else:
                raise ArithmeticError('crossing Newton did not converge')
            if not 0 <= s <= h:
                raise ArithmeticError('crossing outside step')
            return [_horner(c, s) for c in cs], t+s, steps+1
        y = [_horner(c, h) for c in cs]
        t += h
        steps += 1
    raise ArithmeticError('no clock section')


def hp_map(z, bits=160, order=40, eps=1e-40):
    """P(z) at `bits` precision. Returns (z' as mpfr list, |constraint| at return, return time)."""
    with gmpy2.context(gmpy2.get_context(), precision=bits):
        y0 = section_to_state(z)
        y, T, _ = flow_to_section(y0, order=order, eps=eps)
        A, V, q, W, x, P = y
        return [+A, -6*V, +x, +P], abs(constraint(y)), T


def half_map(z, bits=160, order=40, eps=1e-40):
    """h(z): from the downward clock section to the next upward crossing, read in z = (A, p_A, x, p_x).
    The flow commutes with S: (q, q') -> (-q, -q'), which maps upward to downward crossings and fixes z,
    so P = h o h."""
    with gmpy2.context(gmpy2.get_context(), precision=bits):
        y, T, _ = flow_to_section(section_to_state(z), t_min=.5, order=order, eps=eps, direction=1)
        A, V, q, W, x, P = y
        return [+A, -6*V, +x, +P], abs(constraint(y)), T


def to_mpfr(v, bits):
    with gmpy2.context(gmpy2.get_context(), precision=bits):
        return [mpfr(a) if not isinstance(a, str) else mpfr(a) for a in v]
