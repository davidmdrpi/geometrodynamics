"""Nonlinear two-region initial data and fixed-window snapshot observables.

This is not a spacetime evolution. K=0 and scalar normal velocity=0 initially.
The quartet is assumed matter, not an emergent vacuum degree of freedom.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.special import eval_jacobi

FREEZE = '27bf9c6d5840506bd1c786b4cddbf1d0ee2f79a6'
VOLUME = 2*np.pi**2
Q0 = np.sqrt(3)/2
F0 = 7/8
PSI0 = F0**.25
CENTERS = np.array([[1., 0., 0., 0.], [np.cos(1.2), np.sin(1.2), 0., 0.]])
ROTATIONS = [(a, b) for a in range(4) for b in range(a+1, 4)]


class DiskBasis:
    """SO(2)-invariant S3 harmonics; even degrees enforce antipodal symmetry.

    S3 projects to the (x0,x1) unit disk with uniform measure 2 pi dx0 dx1.
    Columns have mean-square one. No spherical or homogeneous truncation.
    """
    def __init__(self, degree):
        if degree < 0 or degree % 2:
            raise ValueError('degree must be a nonnegative even integer')
        self.modes = [(n, m, s) for n in range(0, degree+1, 2)
                      for m in range(0, n+1, 2) for s in range(1 if m == 0 else 2)]
        self.lambdas = np.array([n*(n+2) for n, _, _ in self.modes])

    def evaluate(self, xy, gradient=False):
        xy = np.asarray(xy)
        r = np.linalg.norm(xy, axis=1)
        if np.any(r == 0) or np.any(r > 1+1e-14):
            raise ValueError('evaluation requires 0 < disk radius <= 1')
        theta = np.arctan2(xy[:, 1], xy[:, 0])
        cr, sr = xy[:, 0]/r, xy[:, 1]/r
        columns, gradients = [], []
        for n, m, s in self.modes:
            k = (n-m)//2
            p = eval_jacobi(k, 0, m, 2*r*r-1)
            norm = np.sqrt((n+1)*(1 if m == 0 else 2))
            radial = norm*r**m*p
            trig = np.sin(m*theta) if s else np.cos(m*theta)
            columns.append(radial*trig)
            if gradient:
                dp = np.zeros_like(r) if k == 0 else 2*r*(k+m+1)*eval_jacobi(k-1, 1, m+1, 2*r*r-1)
                dr = norm*((m*r**(m-1)*p if m else 0)+r**m*dp)*trig
                dt = radial*m*(np.cos(m*theta) if s else -np.sin(m*theta))
                gradients.append(np.stack([cr*dr-sr*dt/r, sr*dr+cr*dt/r], axis=1))
        Y = np.array(columns).T
        return (Y, np.stack(gradients, axis=-1)) if gradient else Y


def disk_grid(radial, angular):
    z, w = np.polynomial.legendre.leggauss(radial)
    r = np.sqrt((z+1)/2)
    theta = 2*np.pi*np.arange(angular)/angular
    xy = np.stack([r[:, None]*np.cos(theta), r[:, None]*np.sin(theta)], axis=-1).reshape(-1, 2)
    weights = np.repeat(np.pi*(w/2)*(2*np.pi/angular), angular)
    return xy, weights


def sources(xy, amplitudes):
    """Analytic field source, ambient tangent gradient, and round Laplacian."""
    amplitudes = np.asarray(amplitudes, dtype=float)
    if amplitudes.shape != (2,) or not np.all(np.isfinite(amplitudes)):
        raise ValueError('two finite independently specified amplitudes required')
    z = xy @ CENTERS[:, :2].T
    b = np.cosh(8*z)/np.cosh(8)
    db = 8*np.sinh(8*z)/np.cosh(8)
    ddb = 64*b
    q = Q0*(1+b @ amplitudes)
    dq = Q0*(db*amplitudes) @ CENTERS[:, :2]
    lap = Q0*((1-z*z)*ddb-3*z*db) @ amplitudes
    norm = np.sum(dq*dq, axis=1)-np.sum(xy*dq, axis=1)**2
    f = 1-q*q/6
    if np.min(f) <= .1:
        raise ArithmeticError('registered f domain stop')
    S = norm/f**2+3*q*q/f
    U = 1.5/f**2
    return dict(q=q, dq=dq, lap=lap, norm=norm, f=f, S=S, U=U)


@dataclass
class InitialData:
    degree: int
    amplitudes: tuple
    coefficients: np.ndarray
    iterations: int
    projected_residual: float

    def fields(self, xy):
        basis = DiskBasis(self.degree)
        Y, dY = basis.evaluate(xy, gradient=True)
        psi = Y @ self.coefficients
        if np.min(psi) <= .1 or not np.all(np.isfinite(psi)):
            raise ArithmeticError('registered psi domain stop')
        grad = dY @ self.coefficients
        lap = -Y @ (basis.lambdas*self.coefficients)
        s = sources(xy, self.amplitudes)
        rho = .5*psi**-4*s['S']+s['U']
        curvature = psi**-5*(-8*lap+6*psi)
        H = curvature-2*rho
        return dict(**s, psi=psi, dpsi=grad, rho=rho, curvature=curvature,
                    H=H, Hn=H/np.maximum(1, abs(curvature)+2*abs(rho)))


def solve(degree=20, radial=28, angular=88, amplitudes=(.02, .03)):
    basis = DiskBasis(degree)
    xy, weights = disk_grid(radial, angular)
    Y = basis.evaluate(xy)
    project = Y.T*(weights/VOLUME)
    if np.max(abs(project @ Y-np.eye(len(basis.modes)))) > 1e-10:
        raise ValueError('quadrature under-resolves retained modes')
    s = sources(xy, amplitudes)
    c = np.zeros(len(basis.modes)); c[0] = PSI0
    def residual(coefficients):
        p = Y @ coefficients
        return 8*basis.lambdas*coefficients+project @ ((6-s['S'])*p-2*s['U']*p**5)
    for iteration in range(31):
        psi = Y @ c
        F = residual(c)
        error = float(np.max(abs(F)))
        if error < 1e-11:
            return InitialData(degree, tuple(amplitudes), c, iteration, error)
        if iteration == 30:
            break
        J = np.diag(8*basis.lambdas)+(project*(6-s['S']-10*s['U']*psi**4)) @ Y
        dc = np.linalg.solve(J, -F)
        for power in range(30):
            trial = c+dc*2.**-power
            if np.min(Y @ trial) > .1 and np.max(abs(residual(trial))) < error:
                c = trial
                break
        else:
            raise ArithmeticError('Newton line search failed')
    raise ArithmeticError('Newton iteration limit')


def window(points, center, radius=.45):
    """Paired caps. Ambient derivatives are projected by contractions later."""
    center = np.asarray(center)
    if center.shape != (4,) or not np.isclose(center @ center, 1):
        raise ValueError('unit ambient center required')
    if not 0 < radius < np.pi/2:
        raise ValueError('cap radius must lie in (0, pi/2)')
    z = points @ center
    h = np.maximum(0, (z*z-np.cos(radius)**2)/np.sin(radius)**2)
    return h**4, (8*h**3*z/np.sin(radius)**2)[:, None]*center


def track_snapshot(points, weights, psi, rho, current, centers=CENTERS):
    """Separate fixed-window measurements in a shared ambient rotation frame.

    weights integrate S3, current is the round-metric dual of the physical
    covector j_i. Half-cover integrals report twisted-RP3 inventories. For
    this quotient interpretation the supplied scalar densities must be even,
    current odd, and the snapshot must include the full cover quadrature.
    This routine does not identify a horizon or evolve/move the windows.
    """
    points, weights = np.asarray(points), np.asarray(weights)
    psi, rho, current = map(np.asarray, (psi, rho, current))
    if points.ndim != 2 or points.shape[1] != 4 or weights.shape != (len(points),):
        raise ValueError('full S3 point/weight arrays required')
    if psi.shape != weights.shape or rho.shape != weights.shape or current.shape != points.shape:
        raise ValueError('snapshot arrays have inconsistent shapes')
    if not all(np.all(np.isfinite(a)) for a in (points, weights, psi, rho, current)):
        raise ValueError('finite snapshot required')
    if np.min(psi) <= 0 or np.min(weights) <= 0 or np.min(rho) < 0:
        raise ValueError('positive metric, weights and nonnegative energy required')
    if np.max(abs(np.sum(points*points, axis=1)-1)) > 1e-10:
        raise ValueError('points must lie on S3')
    if np.max(abs(np.sum(points*current, axis=1))) > 1e-10*max(1., np.max(abs(current))):
        raise ValueError('current must be tangent')
    centers = np.asarray(centers)
    if centers.shape != (2, 4):
        raise ValueError('two region centers required')
    if np.arccos(np.clip(abs(centers[0] @ centers[1]), 0, 1)) <= .9:
        raise ValueError('region caps overlap on the quotient')
    measure = weights*psi**6/2
    charges = np.stack([points[:, a]*current[:, b]-points[:, b]*current[:, a]
                        for a, b in ROTATIONS], axis=1)
    out = []
    for name, center in zip(('A', 'B'), centers):
        w, _ = window(points, center)
        energy = float((measure*w) @ rho)
        moment = np.sum((measure*w*rho*np.sign(points @ center))[:, None]*points, axis=0)
        length = np.linalg.norm(moment)
        if energy <= 0 or length <= 1e-14:
            raise ValueError('region has no resolved energy centroid')
        out.append(dict(id=name, volume=float(measure @ w), energy=energy,
                        centroid=(moment/length).tolist(),
                        momentum=(measure*w @ charges).tolist()))
    return out


def diagnostics(data, radial=64, angular=192):
    """Initial matter momentum rate and independent smooth-window balance.

    SO(2) symmetry makes five rotation rates vanish. Explicit orbit expansion
    below measures all six inventories and the four centroid components.
    """
    xy, weights = disk_grid(radial, angular)
    f = data.fields(xy)
    q, psi, target_f = f['q'], f['psi'], f['f']
    loggrad = f['dpsi']/psi[:, None]
    def dot(a, b):
        return np.sum(a*b, axis=1)-np.sum(xy*a, axis=1)*np.sum(xy*b, axis=1)
    X = np.stack([-xy[:, 1], xy[:, 0]], axis=1)
    Xq, Xlog = np.sum(X*f['dq'], axis=1), np.sum(X*loggrad, axis=1)
    # Sigma-model normal acceleration Pdot=A_rad x+A_tangent.
    arad = psi**-4*(f['lap']-3*q+2*dot(loggrad, f['dq'])+q*f['norm']/(3*target_f))-q/target_f
    atan_X = psi**-4*((2+q*q/(3*target_f))*Xq+2*q*Xlog)
    jdot_X = -arad*Xq/target_f**2-q/target_f*atan_X
    trace_stress = psi**-4*f['S']-3*f['rho']
    measure = weights*psi**6/2
    windows = []
    for center in CENTERS:
        points = np.c_[xy, np.sqrt(1-np.sum(xy*xy, axis=1)), np.zeros(len(xy))]
        w, dw = window(points, center)
        windows.append((w, dw[:, :2]))
    windows.append((1-windows[0][0]-windows[1][0], -windows[0][1]-windows[1][1]))
    ledger = []
    for name, (w, dw) in zip(('A', 'B', 'exterior'), windows):
        Xw = np.sum(X*dw, axis=1)
        stress_flux = psi**-4*(Xq*dot(f['dq'], dw)/target_f**2+q*q*Xw/target_f)-f['rho']*Xw
        metric_work = 2*w*Xlog*trace_stress
        direct = float(measure @ (w*jdot_X))
        flux, work = float(measure @ stress_flux), float(measure @ metric_work)
        ledger.append(dict(id=name, direct_rate_01=direct, stress_flux_01=flux,
                           metric_work_01=work, balance_rate_01=flux+work,
                           balance_error=direct-flux-work,
                           other_five_rates=0., other_rates_reason='SO(2) symmetry'))
    # Four equally weighted orbit nodes integrate each centroid/rotation exactly.
    phi = np.arange(4)*np.pi/2
    transverse = np.sqrt(1-np.sum(xy*xy, axis=1))
    points = np.concatenate([np.c_[xy, transverse*np.cos(p), transverse*np.sin(p)] for p in phi])
    regions = track_snapshot(points, np.tile(weights/4, 4), np.tile(psi, 4),
                             np.tile(f['rho'], 4), np.zeros_like(points))
    rho0 = .5*PSI0**-4*(3*Q0**2/F0)+1.5/F0**2
    for row, (w, _) in zip(regions, windows):
        row['background_energy'] = float(weights @ w/2*PSI0**6*rho0)
        row['excess_energy'] = row['energy']-row['background_energy']
    return dict(regions=regions, ledger=ledger,
                total_direct_rate_01=float(measure @ jdot_X),
                total_metric_work_01=float(measure @ (2*Xlog*trace_stress)),
                H_max=float(np.max(abs(f['H']))), H_normalized_max=float(np.max(abs(f['Hn']))),
                psi_min=float(np.min(psi)), f_min=float(np.min(target_f)))


def coordinate_curvature(data, coordinate, step):
    """Ricci scalar from finite differences of the physical coordinate metric.

    Independent of the spectral Laplacian and the source/constraint formula.
    Hopf coordinates (r,theta,phi), away from either coordinate axis.
    """
    basis = DiskBasis(data.degree)
    def metric(z):
        r, theta, _ = z
        xy = np.array([[r*np.cos(theta), r*np.sin(theta)]])
        p = float((basis.evaluate(xy) @ data.coefficients)[0])
        return np.diag([1/(1-r*r), r*r, 1-r*r])*p**4
    eye = np.eye(3)*step
    def christoffel(z):
        g = metric(z)
        dg = np.array([(metric(z+e)-metric(z-e))/(2*step) for e in eye])
        inv = np.linalg.inv(g)
        return np.array([[[sum(inv[k, l]*(dg[i, l, j]+dg[j, l, i]-dg[l, i, j])/2
                               for l in range(3)) for j in range(3)] for i in range(3)] for k in range(3)])
    z = np.asarray(coordinate, dtype=float)
    G = christoffel(z)
    dG = np.array([(christoffel(z+e)-christoffel(z-e))/(2*step) for e in eye])
    Ric = np.array([[sum(dG[k, k, i, j]-dG[j, k, i, k]
                         +sum(G[k, i, j]*G[l, k, l]-G[l, i, k]*G[k, j, l] for l in range(3))
                         for k in range(3)) for j in range(3)] for i in range(3)])
    return float(np.sum(np.linalg.inv(metric(z))*Ric))
