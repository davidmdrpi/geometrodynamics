"""Exact spherical-fibre clock experiment. Run with Python 3, SciPy, NumPy."""
from math import isqrt
import json
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import solve_ivp

L = np.pi / np.sqrt(5)
T = 2 * np.pi
b = 1 / np.sqrt(2)


def modes(limit):
    result = []
    for m in range(limit + 1):
        n = isqrt(4 + 5 * m*m)
        if n*n == 4 + 5*m*m and (n - m) % 2 == 0:
            result.append([n, m])
    return result


def field(t, x, amplitude=1):
    # Same real, normalized axial l=1 vector harmonic in both terms.
    k = np.pi/L
    u = np.cos(2*t) + np.cos(3*t)*np.cos(k*x)
    ut = -2*np.sin(2*t) - 3*np.sin(3*t)*np.cos(k*x)
    ux = -k*np.cos(3*t)*np.sin(k*x)
    return [amplitude*v for v in (u, ut, ux)]


def clock_control(quarters):
    """Compose deck words in units L=1, T=4. The angular slot acts on Ω as a whole."""
    from fractions import Fraction

    identity = np.eye(4, dtype=int)
    p = np.diag([-1, 1, -1, 1])
    d = np.diag([1, 1, -1, 1])
    d[:2, 3] = [1, quarters]

    def word(j, k):
        # C and D commute; use integer translations for their signed powers.
        value = identity.copy()
        value[:2, 3] = [k, 4*j + quarters*k]
        value[2, 2] = -1 if k % 2 else 1
        return value

    # Conjugation must preserve both deck generators; P is its own inverse.
    words = [word(j, k) for j in range(-2, 3) for k in (-1, 0, 1)]
    descends = all(any(np.array_equal(p @ g @ p, w) for w in words)
                   for g in (word(1, 0), d))
    # On the fundamental strip |x| <= 1/2, a fixed point of P C^j D^k
    # requires |k| <= 1. For these four offsets the listed j exhaust time matches.
    fixed = None
    if descends:
        fixed = False
        for w in words:
            action = p @ w
            x = Fraction(int(action[0, 3]), 1 - int(action[0, 0]))
            point = np.array([x, 0, 1, 1], dtype=object)
            if abs(x) <= Fraction(1, 2) and np.array_equal(action @ point, point):
                fixed = True
                break
    return dict(tau_over_T=str(Fraction(quarters, 4)), descends=descends, fixed_orbit=fixed)


def verify():
    t = np.arange(2048)*T/2048
    u, ut, ux = field(t, L/2)
    j = -ut*ux
    rng = np.random.default_rng(91)
    times, xs = rng.uniform(-6, 6, (2, 130))
    seam = max(np.max(np.abs(a-bv)) for a,bv in zip(
        field(times, xs),field(times+T/2, xs+L)))
    base = field(times, xs)
    reflected = field(times, -xs)
    parity_error = max(np.max(np.abs(reflected[k] - sign*base[k]))
                       for k,sign in enumerate([1,1,-1]))
    # Independent second-order ND Sturm-Liouville solve on [0,L/2],
    # cell-centred mesh, even ghost at 0 and odd ghost at L/2.
    refinement = []
    for count in (16,32,64,128):
        dx = L/(2*count)
        diag = np.full(count,2.)
        diag[0],diag[-1] = 1.,3.
        ev = eigh_tridiagonal(diag,-np.ones(count-1),select='i',
                              select_range=(0,0),eigvals_only=True)[0]
        refinement.append(float(abs(np.sqrt(4+ev/dx**2)-3)))
    # Independent two-port propagation across both mouth sections, one full
    # optical period. Monodromy must be -I at 3 and must fail off resonance.
    def transfer(omega):
        def rhs(x, y):
            f = y.reshape(2,2)
            return np.vstack((f[1], (4-omega**2)*f[0])).ravel()
        sol = solve_ivp(rhs,[0,L],np.eye(2).ravel(),rtol=1e-12,atol=1e-14,
                        method='DOP853')
        assert sol.success
        return float(np.linalg.norm(sol.y[:,-1].reshape(2,2)+np.eye(2)))
    # Directed flux balance on a finite interval, independent time and space
    # integration. Use a fixed lifted interval, with both mouths physically
    # distinguished on the cover. Its ledger is not an inventory on the quotient.
    from numpy.polynomial.legendre import leggauss
    z,w = leggauss(48)
    left,right = L*.25,L*.5
    def energy(time, amp=1):
        f,ft,fx = field(time,(right-left)*z/2+(right+left)/2,amp)
        return (right-left)/4 * (w @ (ft*ft + fx*fx + 4*f*f))
    end = .8
    ts = (z+1)*end/2
    def flux(x):
        _,ft,fx = field(ts,x)
        return -ft*fx
    balance = abs(energy(end)-energy(0)-end/2*(w @ (flux(left)-flux(right))))
    rows = modes(20000)
    clocks = [clock_control(num) for num in range(4)]
    assert [(r['descends'],r['fixed_orbit']) for r in clocks]==[(True,True),(False,None),(True,False),(False,None)]
    answer = dict(clock_controls=clocks,modes_n_m=rows, omega3_grid_errors=refinement,
                  seam_residual=seam, antipodal_residual=parity_error,
                  peak_abs_flux=float(np.max(np.abs(j))),
                  mean_flux=float(j.mean()),
                  flux_half_period_residual=float(np.max(np.abs(j+np.roll(j,1024)))),
                  interval_energy_balance=balance,
                  energy_ratio_at_0p37=energy(.3,.37)/energy(.3),
                  monodromy_at_3=transfer(3), monodromy_at_3p01=transfer(3.01),
                  omega_after_one_percent_longer=float(np.sqrt(4+5/1.01**2)))
    assert rows[:4] == [[2,0],[3,1],[7,3],[18,8]]
    assert all(3.9 < p/q < 4.1 for p,q in zip(refinement,refinement[1:]))
    assert max(seam,parity_error,balance,abs(j.mean()),answer['flux_half_period_residual']) < 1e-11
    assert answer['peak_abs_flux'] > 4 and answer['monodromy_at_3'] < 1e-10
    assert answer['monodromy_at_3p01'] > .01
    assert abs(answer['energy_ratio_at_0p37']-.37**2) < 1e-13
    return answer


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
