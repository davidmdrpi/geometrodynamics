"""Continuum test of the registered ground-doublet splitting, with independent roots."""
import json
from math import cos, cosh, gamma, log, pi, sin, sqrt, tan

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def prediction(a, m):
    norm = sqrt(pi)*gamma(m+1)/(2*gamma(m+1.5))
    prefactor = 2*4**m/((m+1)*norm)
    return prefactor*tan(a/2)**(4*m)


def measure(a, m):
    length = -log(tan(a/2))
    b = sin(a)**2

    def propagate(rhs, initial):
        state = initial
        for radius in (lambda t: b*cosh(t), lambda t: 1/cosh(length-t)):
            sol = solve_ivp(lambda t, y: rhs(radius(t)**2,y), (0,length),state,
                            method="DOP853",rtol=1e-12,atol=2e-14)
            if not sol.success:
                raise RuntimeError(sol.message)
            state = sol.y[:,-1]
        return state

    def ground(initial):
        def residual(lam):
            return propagate(lambda r,y: [y[1],(m*m-lam*r)*y[0]],initial)[1]
        lam0 = m*(m+1)
        return brentq(residual,.7*lam0,1.3*lam0,xtol=1e-14)

    ln, ld = ground([1.,0.]), ground([0.,1.])

    # The unwrapped phase atan2(u,u') crosses each interior node upwards
    # through a multiple of pi. Its first Neumann endpoint is pi/2.
    def endpoint_phase(lam, start):
        def phase(r, y):
            angle = y[0]
            return [cos(angle)**2 + (lam*r - m*m)*sin(angle)**2]
        end, = propagate(phase, [start])
        return end

    if max(abs(endpoint_phase(ln, pi/2)-pi/2), abs(endpoint_phase(ld, 0.)-pi/2)) > 1e-8:
        raise ArithmeticError("shooting selected an excited parity level")
    if endpoint_phase(ld, pi/2) >= pi:
        raise ArithmeticError("a lower equator-odd antiperiodic mode exists")

    def overlap(r,y):
        n,p,d,q,_ = y
        return [p,(m*m-ln*r)*n,q,(m*m-ld*r)*d,r*n*d]

    _,_,_,_,integral = propagate(overlap,[1.,0.,0.,1.,0.])
    split = float(2/(integral*(sqrt(ld)+sqrt(ln))**2))
    pred = prediction(a,m)
    return dict(a=a,m=m,lambda_N=ln,lambda_D=ld,split=split,prediction=pred,
                relative_error=pred/split-1)


def independent(a, m, digits=45):
    """Exact hemisphere hypergeometric solution matched to neck Taylor propagation."""
    with mp.workdps(digits):
        a = mp.mpf(str(a))
        b = mp.sin(a)**2
        length = -mp.log(mp.tan(a/2))
        z = mp.cos(a)
        steps = int(mp.ceil(length*3))
        h = length/steps
        order = digits+12
        # Taylor coefficients of rho²=b²(1+cosh(2y))/2 in each neck step.
        coeffs = []
        for j in range(steps):
            even, odd = mp.cosh(2*j*h), mp.sinh(2*j*h)
            row = [b*b*(even+1)/2]
            for k in range(1,order):
                row.append(b*b*(odd if k%2 else even)*2**(k-1)/mp.factorial(k))
            coeffs.append(row)

        def mismatch(nu, odd):
            lam = nu*(nu+1)
            value, slope = (0,1) if odd else (1,0)
            for rho in coeffs:
                v = [mp.mpf(value),mp.mpf(slope)]
                for k in range(order-1):
                    conv = mp.fdot(rho[:k+1], v[k::-1])
                    v.append((m*m*v[k]-lam*conv)/((k+1)*(k+2)))
                value = mp.polyval(v[::-1],h)
                slope = mp.polyval([k*v[k] for k in range(order,0,-1)],h)
            aa, bb = (m-nu)/2,(m+nu+1)/2
            g = mp.hyp2f1(aa,bb,mp.mpf('.5'),z*z)
            deriv = 2*aa*bb*mp.hyp2f1(aa+1,bb+1,mp.mpf('1.5'),z*z)
            # Logarithmic matching with opposite outward coordinates; use a
            # determinant so accidental zeros of the rim value cause no poles.
            return slope*g+value*(-m*z*g+2*z*(1-z*z)*deriv)

        roots = [mp.findroot(lambda v: mismatch(v, odd), (m,mp.mpf(m)+mp.mpf('.0001')))
                 for odd in (False,True)]
        wn,wd = [mp.sqrt(nu*(nu+1)) for nu in roots]
        return float(2*(wd-wn)/(wd+wn))


def probe():
    primary = [measure(a,m) for a in (.08,.12,.18,.24) for m in (1,2,3,4)]
    stress = [measure(a,m) for a in (.4,.6) for m in (2,4)]
    checks = []
    for a,m in ((.08,4),(.24,1),(.4,2)):
        result = independent(a,m)
        row = next(r for r in primary+stress if (r["a"],r["m"])==(a,m))
        checks.append(dict(a=a,m=m,direct_split=result,relative_difference=result/row["split"]-1))
    by_cell = {(r["a"],r["m"]):abs(r["relative_error"]) for r in primary}
    gates = dict(each_within_three_percent=max(by_cell.values())<=.03,
                 endpoints_improve=all(by_cell[.08,m]<by_cell[.24,m] for m in range(1,5)),
                 independent_agreement=max(abs(r["relative_difference"]) for r in checks)<1e-7)
    return dict(primary=primary,stress=stress,independent=checks,gates=gates,passed=all(gates.values()))


if __name__ == "__main__":
    print(json.dumps(probe(),indent=2,allow_nan=False))
