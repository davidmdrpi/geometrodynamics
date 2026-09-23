"""Independent arbitrary-precision coordinate jets of saved diagonal data.

No constraint ODE or conformal-curvature identity is used here. The engine
supports diagonal metrics and diagonal covariant extrinsic curvature.
"""
from bisect import bisect_right
from functools import lru_cache
import mpmath as mp


def norm(v):
    return mp.sqrt(sum(x*x for x in v))


class Polynomial:
    def __init__(self, record, component=None):
        self.x = [mp.mpf(v) for v in record['x']]
        self.c = [[mp.mpf(v if component is None else v[component])
                   for v in row] for row in record['c']]

    def __call__(self, s):
        i = min(len(self.x)-2, max(0, bisect_right(self.x, s)-1))
        z = s-self.x[i]
        value = mp.mpf(0)
        for row in self.c:
            value = value*z+row[i]
        return value


class SavedData:
    def __init__(self, record, profiles):
        self.psi = Polynomial(record['solution'], 0)
        self.theta = Polynomial(profiles['theta'])
        self.tensor = Polynomial(profiles['tensor'])
        self.L, self.eta = mp.mpf(record['L']), mp.mpf(record['eta'])

    def step(self, point):
        distance = min(abs(point[0]-k) for poly in
                       (self.psi, self.theta, self.tensor) for k in poly.x)
        if distance == 0:
            raise ValueError('physical point lies on a polynomial knot')
        return min(mp.mpf('.008'), distance/8)

    def values(self, s, t, azimuth):
        f = mp.mpf(7)/8
        q = mp.sqrt(3)/2
        psi, theta = self.psi(s), self.theta(s)
        a = self.eta*self.tensor(s)
        p = self.eta/mp.cosh(self.L)**3*mp.sin(mp.pi*s/(2*self.L))
        angular = [mp.sin(t)*mp.cos(azimuth), mp.sin(t)*mp.sin(azimuth), mp.cos(t)]
        phi = [q*mp.sin(theta)]+[q*mp.cos(theta)*v for v in angular]
        tangent = [q*mp.cos(theta)]+[-q*mp.sin(theta)*v for v in angular]
        g = [psi**4/f, psi**4/f, psi**4/f*mp.sin(t)**2]
        K = [2*a/psi**2/mp.sqrt(f), -a/psi**2/mp.sqrt(f),
             -a/psi**2/mp.sqrt(f)*mp.sin(t)**2]
        Pi = [mp.sqrt(f)*psi**-6*p*v for v in tangent]
        return g+K+phi+Pi


def coordinate_check(values, point, h=None, f=None, cosmological=None):
    """Direct jets if h=None; otherwise fourth-order value-only differences.

    The values function returns diagonal g, diagonal K, four fields, four Pi.
    All derivatives and contractions retain the active mpmath precision.
    """
    f = mp.mpf(7)/8 if f is None else mp.mpf(f)
    cosmological = mp.mpf(3)/2 if cosmological is None else mp.mpf(cosmological)
    x = tuple(mp.mpf(v) for v in point)
    # mp.diff temporarily raises precision; precision belongs in the cache key.
    @lru_cache(None)
    def cached(precision, args):
        return values(*args)

    def sample(args):
        return cached(mp.mp.prec, tuple(args))

    center = sample(x)
    g, K, phi, Pi = center[:3], center[3:6], center[6:10], center[10:14]
    if any(v <= 0 for v in g):
        raise ValueError('nonpositive metric')
    first = [[mp.mpf(0) for _ in range(10)] for _ in range(3)]
    second = [[[mp.mpf(0) for _ in range(3)] for _ in range(3)] for _ in range(3)]
    if h is None:
        for a in range(3):
            orders = tuple(int(j == a) for j in range(3))
            for k in range(10):
                first[a][k] = mp.diff(lambda *z: sample(z)[k], x, orders)
            for b in range(a+1):
                orders = tuple(int(j == a)+int(j == b) for j in range(3))
                for k in range(3):
                    second[a][b][k] = second[b][a][k] = mp.diff(
                        lambda *z: sample(z)[k], x, orders)
    else:
        h = mp.mpf(h)
        if h <= 0:
            raise ValueError('nonpositive coordinate step')
        offsets = (-2, -1, 1, 2)
        w1 = [mp.mpf(v)/12 for v in (1, -8, 8, -1)]
        w2 = [mp.mpf(v)/12 for v in (-1, 16, 16, -1)]
        def shift(a, o, b=None, u=0):
            z = list(x)
            z[a] += o*h
            if b is not None:
                z[b] += u*h
            return sample(z)
        for a in range(3):
            samples = [shift(a, o) for o in offsets]
            first[a] = [sum(w*v[k] for w, v in zip(w1, samples))/h for k in range(10)]
            second[a][a] = [(sum(w*v[k] for w, v in zip(w2, samples))
                             -mp.mpf(5)/2*g[k])/h**2 for k in range(3)]
            for b in range(a):
                mixed = [sum(wa*wb*shift(a, oa, b, ob)[k]
                             for oa, wa in zip(offsets, w1)
                             for ob, wb in zip(offsets, w1))/h**2 for k in range(3)]
                second[a][b] = second[b][a] = mixed

    def numerator(a, b, c):
        return ((first[b][a] if a == c else 0)
                +(first[c][a] if a == b else 0)
                -(first[a][b] if b == c else 0))
    gamma = [[[numerator(a, b, c)/(2*g[a]) for c in range(3)]
              for b in range(3)] for a in range(3)]
    def dgamma(d, a, b, c):
        dn = ((second[d][b][a] if a == c else 0)
              +(second[d][c][a] if a == b else 0)
              -(second[d][a][b] if b == c else 0))
        return dn/(2*g[a])-first[d][a]/g[a]*gamma[a][b][c]
    ricci = [sum(dgamma(a, a, i, i)-dgamma(i, a, i, a)
                 +sum(gamma[a][a][b]*gamma[b][i][i]
                      -gamma[a][i][b]*gamma[b][i][a] for b in range(3))
                 for a in range(3)) for i in range(3)]
    R = sum(ricci[i]/g[i] for i in range(3))
    mixed = [K[i]/g[i] for i in range(3)]
    dmixed = [[first[a][3+i]/g[i]-K[i]*first[a][i]/g[i]**2
               for i in range(3)] for a in range(3)]
    terms = [[dmixed[i][i]-sum(dmixed[i]) for i in range(3)],
             [sum(gamma[j][j][i]*mixed[i] for j in range(3)) for i in range(3)],
             [-sum(gamma[j][j][i]*mixed[j] for j in range(3)) for i in range(3)]]
    divergence = [sum(row[i] for row in terms) for i in range(3)]
    trace, K2 = sum(mixed), sum(v*v for v in mixed)
    gradient = [[first[i][6+A] for i in range(3)] for A in range(4)]
    kinetic = sum(v*v for v in Pi)
    spatial = sum(gradient[A][i]**2/g[i] for A in range(4) for i in range(3))
    current = [sum(Pi[A]*gradient[A][i] for A in range(4)) for i in range(3)]
    geom = R+trace**2-K2
    H = f*geom-kinetic-spatial-2*cosmological
    M = [f*divergence[i]+current[i] for i in range(3)]
    hs = max(mp.mpf(1), f*(abs(R)+trace**2+K2)+kinetic+spatial+2*abs(cosmological))
    ms = max(mp.mpf(1), f*sum(norm(row) for row in terms)+norm(current))
    return dict(metric=g, K=K, fields=phi, Pi=Pi, gradient=gradient,
                R=R, trace_K=trace, K2=K2, momentum_terms=terms,
                divergence=divergence, kinetic=kinetic, spatial=spatial,
                current=current, H=H, M=M, Hscale=hs, Mscale=ms,
                wrong_f_H=(geom-kinetic-spatial-2*cosmological)/hs)


def calibration_values(sphere=False, tensor=False):
    def values(s, t, azimuth):
        factor = 1/mp.cosh(s)**2 if sphere else mp.mpf(1)
        return ([factor, factor, factor*mp.sin(t)**2]
                +([mp.mpf(0), s, s*mp.sin(t)**2] if tensor else [mp.mpf(0)]*3)
                +[mp.mpf(0)]*8)
    return values
