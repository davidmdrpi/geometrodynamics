"""Uniform-port response on an exact round S³ exterior with four geodesic holes."""
import json
import numpy as np
from scipy.optimize import brentq
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh, splu

from .receiver import rim_area, solve


def chamber(n, radius):
    """Q1 forms on 0 <= u <= 1, 0 <= alpha <= pi/2 around the +e0 mouth."""
    i, j = np.indices((n, n)).reshape(2, -1)
    corner = i * (n + 1) + j
    nodes = corner[:, None] + [0, n + 1, n + 2, 1]
    angular_step = np.pi / (2 * n)
    quadrature, weights = np.polynomial.legendre.leggauss(3)
    shape = (len(i), 4, 4)
    stiffness = np.zeros(shape)
    mass = np.zeros(shape)
    for xi, wx in zip((quadrature + 1) / 2, weights / 2):
        u = (i + xi) / n
        for eta, wy in zip((quadrature + 1) / 2, weights / 2):
            alpha = (j + eta) * angular_step
            limit = np.arctan2(1, np.cos(alpha))
            dr_du = limit - radius
            dr_dalpha = u * np.sin(alpha) / (1 + np.cos(alpha)**2)
            r = radius + u * dr_du
            basis = np.array([(1-xi)*(1-eta), xi*(1-eta), xi*eta, (1-xi)*eta])
            du = n * np.array([eta-1, 1-eta, eta, -eta])
            da = np.array([xi-1, -xi, xi, 1-xi]) / angular_step
            radial = du / dr_du[:, None]
            angular = (da - dr_dalpha[:, None] * radial) / np.sin(r)[:, None]
            # Integrate psi, reflect x0 and x1, divide by two for the quotient.
            volume = 4*np.pi * np.sin(r)**2 * np.sin(alpha) * dr_du * angular_step/n * wx*wy
            stiffness += volume[:, None, None] * (
                radial[:, :, None] * radial[:, None, :] + angular[:, :, None] * angular[:, None, :])
            mass += volume[:, None, None] * np.outer(basis, basis)
    rows = np.broadcast_to(nodes[:, :, None], shape).ravel()
    cols = np.broadcast_to(nodes[:, None, :], shape).ravel()
    return tuple(coo_matrix((a.ravel(), (rows, cols)), shape=((n+1)**2,)*2).tocsc()
                 for a in (stiffness, mass))


class RoundExterior:
    """The second chamber is the first reflected through x0=x1."""
    def __init__(self, n, area):
        self.radius = np.arcsin(np.sqrt(area / (4 * np.pi)))
        stiffness, mass = chamber(n, self.radius)
        self.volume = 2 * mass.sum()
        self.sectors = []
        rim = np.arange(n+1)
        for end in ((n+1)**2, n*(n+1)):
            # Even under chamber exchange: natural outer boundary. Odd: zero trace.
            interior = np.arange(n+1, end)
            blocks = []
            for form in (stiffness, mass):
                load = np.asarray(form[interior][:, rim].sum(axis=1)).ravel()
                blocks.append((form[interior][:, interior], load, form[rim][:, rim].sum()))
            self.sectors.append(blocks)

    def poles(self, count=3):
        values = []
        for (k, _, _), (m, _, _) in self.sectors:
            lam = eigsh(k, M=m, k=count, sigma=0, which='LM',
                        v0=np.linspace(1, 2, k.shape[0]), return_eigenvectors=False)
            values.append(np.sqrt(np.sort(lam)).tolist())
        return values

    def ports(self, w):
        eigenvalues = []
        for (k, kload, kb), (m, mload, mb) in self.sectors:
            load = kload - w*w*mload
            response = splu(k-w*w*m).solve(load)
            eigenvalues.append(kb-w*w*mb-load @ response)
        even, odd = eigenvalues
        return np.array([[even+odd, even-odd], [even-odd, even+odd]]) / 2


def tube_ports(area, w, length=0.95):
    return area*w/np.sin(length*w) * np.array([[np.cos(length*w), -1], [-1, np.cos(length*w)]])


def measure(n, area):
    ext = RoundExterior(n, area)
    r = ext.radius
    expected_volume = np.pi**2 - 4*np.pi*(r-np.sin(r)*np.cos(r))
    frequencies = [0.2, 1., 4.599769873349991]
    rows = []
    for w in frequencies:
        h, e = tube_ports(area, w), ext.ports(w)
        tau = np.pi/w if w == frequencies[-1] else 0.2375
        rows.append(dict(solve(h, e, w, tau), exterior_cross=e[0, 1],
                         zero_clock_cross=abs(e[0, 1]+h[0, 1]),
                         shifted_cross=abs(e[0, 1]+np.exp(1j*w*tau)*h[0, 1])))
    poles = ext.poles(4)
    # Fixed after the 32/64 fits and before the 128/256 dark-port calculations.
    prediction = 5.879137825138902 + (5.879137825138902 - 5.88861095780784) / 3
    bracket = (5.75, 6.0)
    h = tube_ports(area, prediction)
    e = ext.ports(prediction)
    dark = solve(h, e, prediction, np.pi / prediction)
    bright = solve(h, e, prediction, 0.)
    def mismatch(w):
        return ext.ports(w)[0, 1] - tube_ports(area, w)[0, 1]
    # The fourth eigenvalue in each ordered symmetry sector exceeds the bracket.
    assert min(sector[-1] for sector in poles) > bracket[1]
    assert all(not bracket[0] <= w <= bracket[1] for sector in poles for w in sector)
    assert int(bracket[0]*.95/np.pi) == int(bracket[1]*.95/np.pi)
    signs = [mismatch(w) for w in bracket]
    root = brentq(mismatch, *bracket) if signs[0]*signs[1] < 0 else None
    switch = dict(bracket=list(bracket), endpoint_mismatches=signs, fitted_root=root,
                  fixed_frequency=prediction, fixed_shift=np.pi/prediction,
                  fixed_dark_power=dark["receiver_power"], zero_shift_power=bright["receiver_power"],
                  fixed_power_ratio=dark["receiver_power"]/bright["receiver_power"])
    return dict(subdivisions=n, radius=r, switch=switch, quotient_exterior_volume=ext.volume,
                exact_volume=expected_volume, symmetry_sector_poles=poles, responses=rows)


if __name__ == "__main__":
    area = rim_area()
    print(json.dumps([measure(n, area) for n in (8, 16, 32, 64, 128, 256)], indent=2))
