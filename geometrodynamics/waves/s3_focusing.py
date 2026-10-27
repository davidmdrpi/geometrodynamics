"""Antipodal focusing on the round three-sphere, against the two-sphere used so far.

The strict pair model (`antipodal_pair_birth`) routes fronts on S2 with the 2-D caustic phase
-i per focus. BAM's space is S3. This module provides the S3 counterparts of three things that
model and #349 rely on, and checks them independently:

- **Exact refocusing.** For the conformal wave operator on the unit S3, -Delta + 1, a degree-n
  harmonic has frequency n + 1 and antipodal parity (-1)^n, so free evolution for time pi is
  exactly minus the antipodal map, U(pi) = -A, on every mode. A flash released at rest arrives at
  the antipode as a sign-flipped mirror image, with no dispersion and no quarter-period phase.
  On S2 the conformal operator, -Delta + 1/4, has frequencies l + 1/2, so U(pi) maps displacement
  wholly into velocity (the 2-D Gouy phase -i). Zonal (axisymmetric) data make both checkable: on
  S3, v = u sin(psi) obeys the 1-D wave equation on [0, pi] with fixed ends.
- **Time-reversal closure of a flash.** With K(Q, P) = (S Q, -S P), S the charge-conjugated
  antipodal parity, a history closes when U(T) z = K z. `closure_correction` measures how far
  given data are from that, as the energy-orthogonal projection #349 uses.
- **What a finite absorbing mouth receives.** A flash emitted isotropically in its emitter's rest
  frame travels on great circles through the emitter x0 and its antipode, energy per direction
  D^d with D = 1/(gamma (1 - w.n)) (d = 3 on S3, 2 on S2). An absorbing mouth (a geodesic ball)
  takes every ray that meets it. `absorbed_by_mouth` integrates the energy and momentum it takes.
  Every great circle through x0 that meets the ball carries two rays, n and -n, with opposite
  angular momentum, and the ball takes both, one before the focus and one after. So a resting
  emitter gives an absorbing mouth no momentum at any offset; only the emitter's motion does.

Geometric-optics transport is exact for thin zonal fronts on S3, where Huygens' principle holds.
Nothing here models partial reflection at the mouth, throat deformation or backreaction.
"""
import numpy as np
from scipy.special import eval_legendre

from .physical_throat import VacuumThroat


# ---- exact zonal propagation ----------------------------------------------------------------
def midpoint_grid(n):
    """Cell centres psi_j = (j + 1/2) pi/n on (0, pi); the antipodal map is j -> n - 1 - j."""
    return (np.arange(n) + .5)*np.pi/n


def s3_zonal_evolve(u0, p0, t):
    """Exact free evolution of zonal data for the conformal operator on S3, from the sine series
    of v = u sin(psi) (frequencies m = n + 1, m >= 1). Data on `midpoint_grid(len(u0))`."""
    u0, p0 = np.asarray(u0, float), np.asarray(p0, float)
    n = len(u0)
    psi = midpoint_grid(n)
    s = np.sin(psi)
    m = np.arange(1, n + 1)
    basis = np.sin(np.outer(psi, m))                   # discrete sine transform (type II) basis
    weights = np.full(n, 2./n)
    weights[-1] = 1./n
    a = weights*(basis.T @ (u0*s))
    b = weights*(basis.T @ (p0*s))
    v = basis @ (a*np.cos(m*t) + b*np.sin(m*t)/m)
    w = basis @ (-a*m*np.sin(m*t) + b*np.cos(m*t))
    return v/s, w/s


def zonal_fd(u0, p0, t, sphere=3, courant=.2):
    """Independent finite-volume evolution of the zonal conformal wave equation on S2 or S3,
        u_tt = (1/sin^k) (sin^k u')' - mu u,   k = sphere - 1, mu = 1 (S3) or 1/4 (S2),
    with exact cell volumes and the pole faces closed (sin = 0). Velocity Verlet."""
    u, p = np.array(u0, float), np.array(p0, float)
    n = len(u)
    h = np.pi/n
    faces = np.arange(n + 1)*h
    k = sphere - 1
    mu = 1. if sphere == 3 else .25
    flux = np.sin(faces[1:-1])**k
    if k == 2:
        volume = (h - (np.sin(2*faces[1:]) - np.sin(2*faces[:-1]))/2)/2
    else:
        volume = np.cos(faces[:-1]) - np.cos(faces[1:])

    def accel(u):
        f = flux*(u[1:] - u[:-1])/h
        div = np.zeros(n)
        div[:-1] += f
        div[1:] -= f
        return div/volume - mu*u

    steps = int(np.ceil(t/(courant*h)))
    dt = t/steps
    a = accel(u)
    for _ in range(steps):
        p += .5*dt*a
        u += dt*p
        a = accel(u)
        p += .5*dt*a
    return u, p


def zonal_energy(u, p, sphere=3):
    """Energy of zonal data on the midpoint grid: (1/2) int (p^2 + u'^2 + mu u^2) sin^k dpsi
    (the solid-angle factor is dropped)."""
    n = len(u)
    h = np.pi/n
    psi = midpoint_grid(n)
    k = sphere - 1
    mu = 1. if sphere == 3 else .25
    faces = np.arange(1, n)*h
    grad = (u[1:] - u[:-1])/h
    return .5*h*(np.sum((p**2 + mu*u**2)*np.sin(psi)**k) + np.sum(grad**2*np.sin(faces)**k))


# ---- time-reversal closure in modes ---------------------------------------------------------
def closure_correction(omega, parity, q, p, duration):
    """Energy fraction by which data (q, p) must change to satisfy U(T) z = K z, with
    K(Q, P) = (S Q, -S P): the projection z+ = (z + K U(T) z)/2 of #349 (K U(T) is an
    energy-orthogonal involution). Mode by mode: omega, parity S, coefficients q, p."""
    omega, parity, q, p = (np.asarray(x) for x in (omega, parity, q, p))
    c, s = np.cos(omega*duration), np.sin(omega*duration)
    qt, pt = q*c + p*s/omega, p*c - omega*q*s
    qi, pi = parity*qt, -parity*pt
    dq, dp = (qi - q)/2, (pi - p)/2
    weight = np.ones_like(omega) if q.ndim == 1 else 1.
    energy = np.sum(weight*(omega**2*np.abs(q)**2 + np.abs(p)**2))
    return float(np.sum(weight*(omega**2*np.abs(dq)**2 + np.abs(dp)**2))/energy)


def zonal_modes_s3(u0, p0):
    """Sine-series coefficients of v = u sin(psi) and v_t: modes m = n + 1 with frequency m,
    antipodal parity (-1)^n = (-1)^(m+1) and energy weight pi/2 (dropped: uniform)."""
    n = len(u0)
    psi = midpoint_grid(n)
    basis = np.sin(np.outer(psi, np.arange(1, n + 1)))
    weights = np.full(n, 2./n)
    weights[-1] = 1./n
    m = np.arange(1, n + 1)
    return m.astype(float), weights*(basis.T @ (u0*np.sin(psi))), weights*(basis.T @ (p0*np.sin(psi)))


def zonal_modes_s2(u0, p0, lmax=None):
    """Legendre coefficients of zonal data on S2 (Gauss-Legendre in cos psi), with the energy
    normalisation folded in so that sum(omega^2 q^2 + p^2) is the energy up to a constant."""
    n = len(u0)
    lmax = n - 1 if lmax is None else lmax
    x, w = np.polynomial.legendre.leggauss(n)
    psi = np.arccos(x)[::-1]
    grid = midpoint_grid(n)
    q0 = np.interp(psi, grid, u0)
    pp = np.interp(psi, grid, p0)
    ls = np.arange(lmax + 1)
    P = np.array([eval_legendre(l, np.cos(psi)) for l in ls])
    w = w[::-1]
    norm = np.sqrt((2*ls + 1)/2)
    return ls, norm*(P @ (w*q0)), norm*(P @ (w*pp))


# ---- a finite absorbing mouth ---------------------------------------------------------------
def absorbed_by_mouth(a, delta, speed, across=False, sphere=3, nodes=400):
    """Energy and momentum a geodesic ball of radius a takes from one flash.

    The emitter sits at e0 and moves with velocity speed*e1; it emits isotropically in its rest
    frame, total rest energy 1, so the lab flash has energy gamma and momentum gamma*speed. The
    flash focuses at -e0. The ball is centred delta from -e0, along e1 (the line of motion) or,
    with across=True, along e2. A ray along unit tangent n is the great circle cos s e0 + sin s n;
    it meets the ball iff |n.d| >= cos(alpha), sin(alpha) = sin(a)/sin(delta) (every ray if
    delta <= a), so the ball takes two antipodal caps of directions. Each ray carries angular
    momentum dE e0 ^ n, which gives the ball linear momentum dE (n (e0.c) - e0 (n.c)).

    Returns dict(energy fraction of the lab flash; momentum along the emitter's motion carried to
    the ball, `along`, and the rest, `across`, in units of the flash momentum gamma*speed when
    speed > 0, else absolute)."""
    if not 0 < a < np.pi/2 or not 0 <= delta < np.pi - a or not 0 <= speed < 1:
        raise ValueError('Need 0 < a < pi/2, 0 <= delta < pi - a and 0 <= speed < 1')
    dim = sphere + 1
    e = np.eye(dim)
    d = e[2] if across else e[1]
    c = -np.cos(delta)*e[0] + np.sin(delta)*d
    # The emitter's direction e1, carried to the ball: along the circle of motion it arrives as
    # -(sin(delta) e0 + cos(delta) e1); across, e1 is normal to the e0-e2 circle and stays -e1 at
    # the focus side. At delta = 0 both are -e1, where a full capture receives -gamma*speed*e1.
    forward = -(e[1] if across else np.sin(delta)*e[0] + np.cos(delta)*e[1])
    gamma = 1/np.sqrt(1 - speed**2)
    cos_alpha = 0. if delta <= a else np.sqrt(max(0., 1 - (np.sin(a)/np.sin(delta))**2))
    energy, momentum = 0., np.zeros(dim)
    if sphere == 3:
        others = [v for v in (e[1], e[2], e[3]) if not np.allclose(v, d)]
        x, wx = np.polynomial.legendre.leggauss(nodes)               # x = cos(beta) about +d
        phi = 2*np.pi*(np.arange(2*nodes) + .5)/(2*nodes)
        for sign in (1., -1.):                                        # the caps about +d and -d
            cb = cos_alpha + (1 - cos_alpha)*(x + 1)/2
            wb = wx*(1 - cos_alpha)/2
            sb = np.sqrt(1 - cb**2)
            for cbi, sbi, wbi in zip(cb, sb, wb):
                n = (sign*cbi*d[None, :] + sbi*(np.cos(phi)[:, None]*others[0] + np.sin(phi)[:, None]*others[1]))
                D = 1/(gamma*(1 - speed*n[:, 1]))
                dE = wbi*(2*np.pi/len(phi))*D**3/(4*np.pi)
                energy += dE.sum()
                momentum += (dE[:, None]*(n*(e[0] @ c) - np.outer(n @ c, e[0]))).sum(0)
    else:
        alpha = np.arccos(cos_alpha)
        base = np.arctan2(d[2], d[1])
        t, wt = np.polynomial.legendre.leggauss(nodes)
        for centre in (base, base + np.pi):
            ang = centre + alpha*t
            n = np.cos(ang)[:, None]*e[1] + np.sin(ang)[:, None]*e[2]
            D = 1/(gamma*(1 - speed*n[:, 1]))
            dE = alpha*wt*D**2/(2*np.pi)
            energy += dE.sum()
            momentum += (dE[:, None]*(n*(e[0] @ c) - np.outer(n @ c, e[0]))).sum(0)
    unit = gamma*speed if speed > 0 else 1.
    perpendicular = momentum - (momentum @ forward)*forward
    return dict(energy_fraction=float(energy/gamma), along=float(momentum @ forward)/unit,
                across=float(np.linalg.norm(perpendicular))/unit, tangent_check=float(abs(momentum @ c)))


def neck_length_s3(mouth_angle):
    """Proper length of the repository's three-dimensional vacuum bridge (`VacuumThroat`): neck
    radius sin^3 a, glued at the rim sin a. The S2 model uses sin 2a."""
    return float(VacuumThroat(mouth_angle).length())
