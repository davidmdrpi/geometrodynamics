"""Truncated multivariate power series (jet transport).

A Jet is a polynomial in `nvar` variables truncated at total degree `order`,
with real or complex coefficients. Arithmetic, analytic functions and
composition are exact up to the truncation order. Used to Taylor-expand
Poincare return maps without finite differences.
"""
from functools import lru_cache
from itertools import combinations_with_replacement
from math import factorial
import numpy as np


@lru_cache(maxsize=None)
def space(nvar, order):
    mons = [()]
    for deg in range(1, order+1):
        mons += list(combinations_with_replacement(range(nvar), deg))
    exps = np.zeros((len(mons), nvar), dtype=int)
    for i, m in enumerate(mons):
        for v in m:
            exps[i, v] += 1
    index = {tuple(e): i for i, e in enumerate(exps)}
    deg = exps.sum(1)
    I, J, T = [], [], []
    for i in range(len(mons)):
        for j in range(len(mons)):
            if deg[i]+deg[j] <= order:
                I.append(i)
                J.append(j)
                T.append(index[tuple(exps[i]+exps[j])])
    return dict(nvar=nvar, order=order, exps=exps, index=index, deg=deg,
                I=np.array(I), J=np.array(J), T=np.array(T), n=len(mons))


class Jet:
    __array_priority__ = 1000

    def __init__(self, c, nvar, order):
        self.s = space(nvar, order)
        self.c = np.asarray(c)

    # construction -----------------------------------------------------------
    @classmethod
    def const(cls, value, nvar, order, dtype=float):
        c = np.zeros(space(nvar, order)['n'], dtype=dtype)
        c[0] = value
        return cls(c, nvar, order)

    @classmethod
    def var(cls, k, value, nvar, order, dtype=float):
        j = cls.const(value, nvar, order, dtype)
        e = [0]*nvar
        e[k] = 1
        j.c[j.s['index'][tuple(e)]] = 1
        return j

    def _new(self, c):
        return Jet(c, self.s['nvar'], self.s['order'])

    def _lift(self, other):
        if isinstance(other, Jet):
            return other
        c = np.zeros(self.s['n'], dtype=np.result_type(self.c, other))
        c[0] = other
        return self._new(c)

    @property
    def value(self):
        return self.c[0]

    def coef(self, exps):
        return self.c[self.s['index'][tuple(exps)]]

    # arithmetic -------------------------------------------------------------
    def __add__(self, o):
        o = self._lift(o)
        return self._new(self.c+o.c)
    __radd__ = __add__

    def __neg__(self):
        return self._new(-self.c)

    def __sub__(self, o):
        return self+(-self._lift(o))

    def __rsub__(self, o):
        return self._lift(o)-self

    def __mul__(self, o):
        if not isinstance(o, Jet):
            return self._new(self.c*o)
        s = self.s
        w = self.c[s['I']]*o.c[s['J']]
        if np.iscomplexobj(w):
            c = np.bincount(s['T'], w.real, s['n'])+1j*np.bincount(s['T'], w.imag, s['n'])
        else:
            c = np.bincount(s['T'], w, s['n'])
        return self._new(c)
    __rmul__ = __mul__

    def __truediv__(self, o):
        if not isinstance(o, Jet):
            return self._new(self.c/o)
        return self*o.reciprocal()

    def __rtruediv__(self, o):
        return self._lift(o)*self.reciprocal()

    def __pow__(self, k):
        if not (isinstance(k, int) and k >= 0):
            raise ValueError('nonnegative integer powers only')
        out = self._lift(1.)
        for _ in range(k):
            out = out*self
        return out

    # analytic functions: f(c+p) = sum f^(k)(c)/k! p^k, p nilpotent ----------
    def _series(self, derivs):
        p = self-self.value
        out = self._lift(derivs[0])
        pk = self._lift(1.)
        for k in range(1, self.s['order']+1):
            pk = pk*p
            out = out+pk*(derivs[k]/factorial(k))
        return out

    def reciprocal(self):
        a = self.value
        return self._series([(-1)**k*factorial(k)*a**(-k-1) for k in range(self.s['order']+1)])

    def sqrt(self):
        a = self.value
        d, coef = [], 1.
        for k in range(self.s['order']+1):
            d.append(coef*a**(.5-k))
            coef *= (.5-k)
        return self._series(d)

    def exp(self):
        return self._series([np.exp(self.value)]*(self.s['order']+1))

    def sin(self):
        a = self.value
        cyc = [np.sin(a), np.cos(a), -np.sin(a), -np.cos(a)]
        return self._series([cyc[k % 4] for k in range(self.s['order']+1)])

    def cos(self):
        a = self.value
        cyc = [np.cos(a), -np.sin(a), -np.cos(a), np.sin(a)]
        return self._series([cyc[k % 4] for k in range(self.s['order']+1)])

    # calculus and composition ----------------------------------------------
    def deriv(self, k):
        s = self.s
        c = np.zeros_like(self.c)
        for i, e in enumerate(s['exps']):
            if e[k] > 0:
                e2 = e.copy()
                e2[k] -= 1
                c[s['index'][tuple(e2)]] += e[k]*self.c[i]
        return self._new(c)

    def truncate(self, order):
        c = self.c.copy()
        c[self.s['deg'] > order] = 0
        return self._new(c)

    def homogeneous(self, degree):
        c = np.zeros_like(self.c)
        mask = self.s['deg'] == degree
        c[mask] = self.c[mask]
        return self._new(c)

    def compose(self, args):
        """self(args), args a list of Jets (any space) with the same count as nvar."""
        s = self.s
        one = args[0]._lift(1.)
        powers = []
        for a in args:
            pw = [one]
            for _ in range(s['order']):
                pw.append(pw[-1]*a)
            powers.append(pw)
        out = one*0
        for i, e in enumerate(s['exps']):
            if self.c[i] == 0:
                continue
            term = one*self.c[i]
            for v, k in enumerate(e):
                if k:
                    term = term*powers[v][k]
            out = out+term
        return out

    def __call__(self, x):
        x = np.asarray(x)
        return sum(self.c[i]*np.prod(x**e) for i, e in enumerate(self.s['exps']))


def variables(point, order, dtype=float):
    n = len(point)
    return [Jet.var(k, point[k], n, order, dtype) for k in range(n)]


def linear_part(jets):
    n = jets[0].s['nvar']
    return np.array([[j.coef(np.eye(n, dtype=int)[k]) for k in range(n)] for j in jets])


def jacobian(jets):
    n = jets[0].s['nvar']
    return [[j.deriv(k) for k in range(n)] for j in jets]
