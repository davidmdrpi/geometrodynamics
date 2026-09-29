"""Small ordinary-coefficient multivariate Taylor algebra (total degree <= 3)."""
from functools import lru_cache
from itertools import product
import math
import numpy as np


@lru_cache(None)
def ring(n, order=3):
    return Ring(n, order)


class Ring:
    def __init__(self, n, order):
        self.n, self.order = n, order
        self.powers = sorted((p for p in product(range(order+1), repeat=n)
                              if sum(p) <= order), key=lambda p: (sum(p), p))
        self.index = {p: i for i, p in enumerate(self.powers)}
        self.degrees = np.array([sum(p) for p in self.powers])
        self.size = len(self.powers)
        triples = [(i, j, self.index[tuple(a+b for a, b in zip(p, q))])
                   for i, p in enumerate(self.powers) for j, q in enumerate(self.powers)
                   if sum(p)+sum(q) <= order]
        self.left, self.right, self.dest = np.array(triples).T

    def constant(self, x=0.):
        c = np.zeros(self.size, dtype=np.result_type(x, float))
        c[0] = x
        return Jet(self, c)

    def variable(self, j):
        c = self.constant()
        p = tuple(int(i == j) for i in range(self.n))
        c.c[self.index[p]] = 1.
        return c


class Jet:
    __array_priority__ = 1000

    def __init__(self, algebra, coefficients):
        self.r = algebra
        self.c = np.asarray(coefficients)

    def other(self, value):
        if isinstance(value, Jet):
            if self.r is not value.r:
                raise ValueError('different Taylor rings')
            return value
        return self.r.constant(value)

    def __add__(self, other):
        return Jet(self.r, self.c+self.other(other).c)

    __radd__ = __add__

    def __neg__(self):
        return Jet(self.r, -self.c)

    def __sub__(self, other):
        return self+-self.other(other)

    def __rsub__(self, other):
        return self.other(other)+-self

    def __mul__(self, other):
        if not isinstance(other, Jet):
            return Jet(self.r, self.c*other)
        other = self.other(other)
        weights = self.c[self.r.left]*other.c[self.r.right]
        if np.iscomplexobj(weights):
            c = np.bincount(self.r.dest, weights=weights.real, minlength=self.r.size)
            c = c+1j*np.bincount(self.r.dest, weights=weights.imag, minlength=self.r.size)
        else:
            c = np.bincount(self.r.dest, weights=weights, minlength=self.r.size)
        return Jet(self.r, c)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        if isinstance(exponent, int) and exponent >= 0:
            out = self.r.constant(1.)
            for _ in range(exponent):
                out = out*self
            return out
        if self.c[0] == 0:
            raise ValueError('noninteger/negative power at zero constant')
        delta = (self-self.c[0])/self.c[0]
        out, power, coefficient = self.r.constant(1.), self.r.constant(1.), 1.
        for k in range(1, self.r.order+1):
            power = power*delta
            coefficient *= (exponent-k+1)/k
            out = out+coefficient*power
        return self.c[0]**exponent*out

    def __truediv__(self, other):
        return self*(other**-1) if isinstance(other, Jet) else self*(1/other)

    def exp(self):
        delta = self-self.c[0]
        out, power = self.r.constant(1.), self.r.constant(1.)
        for k in range(1, self.r.order+1):
            power = power*delta
            out = out+power/math.factorial(k)
        return np.exp(self.c[0])*out

    def derivative(self, j):
        out = self.r.constant(0.*self.c[0])
        for i, p in enumerate(self.r.powers):
            if p[j]:
                q = list(p)
                q[j] -= 1
                out.c[self.r.index[tuple(q)]] += p[j]*self.c[i]
        return out

    def homogeneous(self, degree):
        return Jet(self.r, np.where(self.r.degrees == degree, self.c, 0))

    def compose(self, variables):
        target = variables[0].r
        powers = [[v**k for k in range(self.r.order+1)] for v in variables]
        out = target.constant(0.*self.c[0])
        for coefficient, indices in zip(self.c, self.r.powers):
            if coefficient == 0:
                continue
            term = target.constant(coefficient)
            for j, k in enumerate(indices):
                term = term*powers[j][k]
            out = out+term
        return out

    def evaluate(self, point):
        return sum(c*np.prod(np.asarray(point)**p) for c, p in zip(self.c, self.r.powers))
