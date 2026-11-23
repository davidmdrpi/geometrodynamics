"""Five-dimensional radial shear symbol for R + f(phi) GB - (dphi)^2/2."""
import json
from itertools import combinations

import numpy as np
import sympy as sp

SIGN = np.array([-1, 1, 1, 1, 1])
METRIC = np.diag(SIGN)


def put_plane(tensor, a, b, c, d, value):
    """Add a symmetric bilinear form on antisymmetric index pairs."""
    for i, j, s in ((a, b, 1), (b, a, -1)):
        for k, l, t in ((c, d, 1), (d, c, -1)):
            tensor[i, j, k, l] += s*t*value
            if (a, b) != (c, d):
                tensor[k, l, i, j] += s*t*value


def curvature(x, y, z):
    riem = np.zeros((5,)*4)
    for a, b in combinations(range(5), 2):
        if (a, b) == (0, 1):
            section = x
        elif a == 0:
            section = 0
        elif a == 1:
            section = y
        else:
            section = z
        put_plane(riem, a, b, a, b, SIGN[a]*SIGN[b]*section)
    return riem


def wedge(a, b):
    return np.einsum("ac,bd->abcd", a, b) - np.einsum("ad,bc->abcd", a, b)


def field_equations(riem, hessian, f, p):
    """Metric equation (both indices down) and GB from full covariant contractions."""
    ric = np.einsum("a,abad->bd", SIGN, riem)
    scalar = np.sum(SIGN*np.diag(ric))
    raised = ric*SIGN[:, None]*SIGN
    factors = np.einsum("a,b,c,d->abcd", SIGN, SIGN, SIGN, SIGN)
    gb = np.sum(riem*riem*factors) - 4*np.sum(ric*raised) + scalar**2
    lanczos = 2*(scalar*ric - 2*(ric*SIGN) @ ric
                 - 2*np.einsum("cd,acbd->ab", raised, riem)
                 + np.einsum("acde,bcde,c,d,e->ab", riem, riem, SIGN, SIGN, SIGN))
    lanczos -= METRIC*gb/2
    p_tensor = (riem - wedge(METRIC, ric) - wedge(ric, METRIC)
                + scalar*wedge(METRIC, METRIC)/2)
    coupling = 4*np.einsum("acbd,cd->ab", p_tensor, hessian)
    grad = np.array([0., p, 0., 0., 0.])
    return (ric - scalar*METRIC/2 + f*lanczos + coupling
            - np.outer(grad, grad)/2 + METRIC*p**2/4), gb


def shear_symbol(riem, hessian, f, p, freq, wave, angular=0., polarization=(3, 4)):
    """Complex-step principal variation of the full equations for a transverse shear."""
    variation = np.zeros_like(riem)
    xi = [freq, wave, angular, 0., 0.]
    a, b = polarization
    if a == b or xi[a] or xi[b]:
        raise ValueError("choose distinct transverse polarization indices")
    for i in np.flatnonzero(xi):
        for j in np.flatnonzero(xi):
            put_plane(variation, i, a, j, b, xi[i]*xi[j]/2)
    step = 1e-20
    metric, gb = field_equations(riem+1j*step*variation, hessian, f, p)
    return metric.imag/step, gb.imag/step


def exact_reduction():
    # Orthogonal sectional curvatures X=-N''/N, Y=-r''/r and Z=1/r².
    y, z, f, p, d0, d1 = sp.symbols("Y Z f p D0 D1", nonzero=True)
    radial = -3*z + 12*z*d0 - p**2/4
    temporal = -3*(y+z) - 12*f*y*z + 12*z*d1 + p**2/4
    slopes = sp.solve([radial, temporal], [d0, d1])
    ct = sp.factor((1+4*f*y-4*d1).subs(slopes))
    cr = sp.factor((1-4*d0).subs(slopes))
    return dict(time_coefficient=str(ct), radial_coefficient=str(cr),
                radial_speed_squared=str(sp.factor(cr/ct)))


def audit():
    # Exact on-shell jet. Set Z=1, Y=-1/2, p=f_phi=1, f=0; solve the two
    # constraints, then the angular and scalar equations for X and phi''.
    y, z, p, f = -0.5, 1., 1., 0.
    d0 = (3*z+p*p/4)/(12*z)
    d1 = (3*(y+z)+12*f*y*z-p*p/4)/(12*z)
    hessian = np.diag([-d0, d1, 0., 0., 0.])  # both indices raised
    def angular(x):
        return field_equations(curvature(x,y,z), hessian, f, p)[0][2, 2]
    x = -angular(0.)/(angular(1.)-angular(0.))
    riem = curvature(x,y,z)
    eq, gb = field_equations(riem,hessian,f,p)
    lapse_slope = d0  # N=1, f_phi=p=1
    phi_second = -lapse_slope*p-gb
    f_second = (d1-phi_second)/p**2
    # Check against independently contracted off-shell jets including f != 0,
    # with mixed temporal/radial covectors. The scalar GB variation must vanish.
    rng = np.random.default_rng(236)
    worst = 0.
    for _ in range(32):
        xx, yy, ff, pp, h0, h1, om, k = rng.normal(size=8)
        zz = rng.uniform(0.3, 2.)
        rr = curvature(xx,yy,zz)
        hh = np.diag([-h0,h1,0.,0.,0.])
        tensor, scalar_source = shear_symbol(rr,hh,ff,pp,om,k)
        prediction = (-(1+4*ff*yy-4*h1)*om**2 + (1-4*h0)*k**2)/2
        expected = np.zeros((5,5))
        expected[3,4] = expected[4,3] = prediction
        background, _ = field_equations(rr,hh,ff,pp)
        temporal = -3*(yy+zz)-12*ff*yy*zz+12*zz*h1+pp**2/4
        radial = -3*zz+12*zz*h0-pp**2/4
        worst = max(worst, float(abs(tensor-expected).max()), abs(scalar_source),
                    abs(-background[0,0]-temporal), abs(background[1,1]-radial))
    tpart, _ = shear_symbol(riem,hessian,f,p,1.,0.)
    spart, _ = shear_symbol(riem,hessian,f,p,0.,1.)
    return dict(symbolic=exact_reduction(), random_contraction_error=float(worst),
                jet=dict(X=x, Y=y, Z=z, phi_prime=p, f=f, f_phi=1.,
                         N_log_prime=lapse_slope, phi_second=phi_second, f_phi_phi=f_second,
                         metric_residual=float(abs(eq).max()),
                         scalar_residual=abs(phi_second+lapse_slope*p+gb),
                         radial_speed_squared=float(-spart[3,4]/tpart[3,4])))


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
