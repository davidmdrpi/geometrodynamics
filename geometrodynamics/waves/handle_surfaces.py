"""Null expansions and oriented constraint fluxes, not particle momenta.

K_ij = -L_n gamma_ij/2; future null normals are n +/- e_s (dot=-2).
The state order is A, B, k, l, u, v, P, Q in the Einstein frame.
"""
import numpy as np


def expansions(state, spatial_derivative):
    y = np.asarray(state, dtype=float)
    dy = np.asarray(spatial_derivative, dtype=float)
    if y.ndim != 2 or y.shape[0] != 8 or dy.shape != y.shape:
        raise ValueError('expected matching (8,N) state and derivative')
    if not np.isfinite(y).all() or not np.isfinite(dy).all():
        raise ValueError('nonfinite surface data')
    A, B, _, ell, u, v, P, Q = y
    f = 1-(u*u+v*v)/6
    if np.min(A) <= 0 or np.min(B) <= 0 or np.min(f) <= 0:
        raise ValueError('positive metric factors and f required')
    ft = -(u*P+v*Q)/3
    fs = -(u*dy[4]+v*dy[5])/3
    h = 2*dy[1]/(A*B)
    einstein = np.array([-2*ell+h, -2*ell-h])
    jordan = np.sqrt(f)*np.array([einstein[0]-ft/f-fs/(A*f),
                                 einstein[1]-ft/f+fs/(A*f)])
    return dict(Einstein=einstein, Jordan=jordan)


def classify(pair, tolerance=1e-8):
    """Tolerance labels near-zero data unresolved, never strictly trapped."""
    p = np.asarray(pair, dtype=float)
    if p.ndim != 2 or p.shape[0] != 2 or not np.isfinite(p).all():
        raise ValueError('expected finite (2,N) expansions')
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError('invalid tolerance')
    return dict(future_trapped=np.all(p < -tolerance, axis=0),
                past_trapped=np.all(p > tolerance, axis=0),
                untrapped=(np.min(p, axis=0) < -tolerance)
                          & (np.max(p, axis=0) > tolerance),
                near_marginal=np.any(abs(p) <= tolerance, axis=0))


def section_flux(state, index, normal=1):
    """Integral (K_ij-trK gamma_ij) X^i nu^j dA, X=partial_s.

Unnormalised geometric constraint flux; no 1/(8 pi G) prefactor,
no asymptotic frame, and no claim this is an object's linear momentum.
"""
    y = np.asarray(state, dtype=float)
    if y.ndim != 2 or y.shape[0] != 8 or not np.isfinite(y).all():
        raise ValueError('expected finite (8,N) state')
    if normal not in (-1, 1) or np.min(y[:2]) <= 0:
        raise ValueError('unit normal sign and positive metric required')
    return float(-normal*8*np.pi*y[0,index]*y[1,index]**2*y[3,index])


def periodic_derivative(values, dx):
    """Independent sixth-order derivative for periodic geometric scalars."""
    v = np.asarray(values)
    return sum(c*np.roll(v, -j, axis=-1) for j,c in
               [(-3,-1),(-2,9),(-1,-45),(1,45),(2,-9),(3,1)])/(60*dx)


def area_expansions(state, dx):
    """Independent area-rate route: theta=2(nR +/- e_s R)/R.

Differentiate each frame's areal radius directly, without using the
conformal-expansion formula or the evolution spatial stencil.
"""
    A, B, _, ell, u, v, P, Q = np.asarray(state)
    f = 1-(u*u+v*v)/6
    ft = -(u*P+v*Q)/3
    R = B/np.sqrt(f)
    nt = -B*ell-B*ft/(2*f)
    radial_e = periodic_derivative(B,dx)/A
    radial_j = np.sqrt(f)*periodic_derivative(R,dx)/A
    return dict(Einstein=2*np.array([-B*ell+radial_e,-B*ell-radial_e])/B,
                Jordan=2*np.array([nt+radial_j,nt-radial_j])/R)
