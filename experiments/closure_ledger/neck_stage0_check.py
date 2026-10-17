#!/usr/bin/env python
"""Stage-0 numbers for docs/neck_collapse_protocol.md: catenoid neck between two rims.

Checks the reduced cosh-profile family against Goldschmidt's catenoid results and
prints the baseline neck parameters quoted in the protocol. No field or mouth dynamics.

    python experiments/closure_ledger/neck_stage0_check.py
"""
import numpy as np
from scipy import integrate, optimize

A_RIM = 1.0              # rim radius a (neck units)
H_BASE = 0.60            # baseline half-separation h/a
MU, OMEGA = .05, 5.      # #340 mouth shape mass and frequency
BARRIERS = (.0025, .005, .01, .02)


def area(r, h, a=A_RIM):
    """Area of rho(z) = r cosh(beta z), beta = arccosh(a/r)/h, for |z| <= h."""
    beta = np.arccosh(a/r)/h
    f = lambda z: r*np.cosh(beta*z)*np.sqrt(1+(r*beta*np.sinh(beta*z))**2)
    return 4*np.pi*integrate.quad(f, 0, h, epsabs=1e-13, epsrel=1e-13, limit=400)[0]


def catenoid_waists(h, a=A_RIM):
    """Waists r with r cosh(h/r) = a: (thin barrier, thick stable)."""
    rmin = optimize.minimize_scalar(lambda r: r*np.cosh(h/r), bounds=(.02, a),
                                    method='bounded', options={'xatol': 1e-14}).x
    g = lambda r: r*np.cosh(h/r)-a
    return optimize.brentq(g, .02, rmin, xtol=1e-15), optimize.brentq(g, rmin, a, xtol=1e-15)


def critical_ratio():
    """Largest h/a with a catenoid: min_r r cosh(h/r) = a."""
    def excess(h):
        return optimize.minimize_scalar(lambda r: r*np.cosh(h/r), bounds=(.02, A_RIM), method='bounded',
                                        options={'xatol': 1e-14}).fun-A_RIM
    return optimize.brentq(excess, .5, .7, xtol=1e-14)


def crossing_ratio():
    """h/a where the stable catenoid's area equals two discs, 2 pi a^2."""
    return optimize.brentq(lambda h: area(catenoid_waists(h)[1], h)-2*np.pi*A_RIM**2, .45, .65, xtol=1e-12)


def family_stationary_points(h, n=400, step=1e-6):
    d = lambda r: (area(r+step, h)-area(r-step, h))/(2*step)
    rs = np.linspace(.02, .999, n)
    v = np.array([d(r) for r in rs])
    return [optimize.brentq(d, rs[i], rs[i+1], xtol=1e-12) for i in range(n-1) if v[i]*v[i+1] < 0]


def main():
    hc, hx = critical_ratio(), crossing_ratio()
    print(f'critical h/a = {hc:.5f} (Goldschmidt 0.6627); disc crossing h/a = {hx:.5f} (0.5277)')
    assert abs(hc-.66274) < 1e-4 and abs(hx-.52770) < 1e-4
    for h in (.55, .60, .65):
        rb, r0 = catenoid_waists(h)
        crit = family_stationary_points(h)
        assert len(crit) == 2 and abs(crit[0]-rb) < 1e-6 and abs(crit[1]-r0) < 1e-6
        a0, ab = area(r0, h), area(rb, h)
        print(f'h/a = {h:.2f}: r0 = {r0:.6f}, rb = {rb:.6f}, family stationary points {crit[0]:.6f}, {crit[1]:.6f}; '
              f'barrier {(ab-a0)/a0:.4f} A0, release to two discs {(a0-2*np.pi)/a0:.4f} A0')
    rb, r0 = catenoid_waists(H_BASE)
    a0, ab = area(r0, H_BASE), area(rb, H_BASE)
    e = 1e-4
    a2 = (area(r0+e, H_BASE)-2*a0+area(r0-e, H_BASE))/e**2
    print(f'baseline h/a = {H_BASE}: A(r0) = {a0:.6f}, A(rb) = {ab:.6f}, two discs = {2*np.pi:.6f}, '
          f'A(1e-4) = {area(1e-4, H_BASE):.6f}, A\'\'(r0) = {a2:.4f}')
    print('dV       gamma     zeta    Q_b      harmonic  waist inertia')
    for dv in BARRIERS:
        gamma = dv/(ab-a0)
        zeta = np.sqrt(2*MU*OMEGA**2/(gamma*a2))
        print(f'{dv:<8} {gamma:.5f}  {zeta:.3f}  {(rb-r0)/zeta:+.3f}   {-np.sqrt(dv/(MU*OMEGA**2)):+.3f}    {2*MU/zeta**2:.4f}')


if __name__ == '__main__':
    main()
