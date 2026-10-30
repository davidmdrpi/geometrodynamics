"""Exact dipole-sector Maxwell fields and transparent worldtubes on unit S3.

This is a 3+1-dimensional Einstein-cylinder reference, not the S2 router and
not a wormhole scattering solution. An observation surface does not absorb
the field. Its stress flux must be combined with stored field momentum before
interpreting it as a force. Units: R=c=epsilon0=mu0=1.
"""
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import eval_gegenbauer
from scipy.linalg import eigh_tridiagonal


def unit(v):
    v = np.asarray(v, float)
    if not np.isfinite(v).all() or np.linalg.norm(v) == 0:
        raise ValueError("A finite nonzero vector is required")
    return v / np.linalg.norm(v)


def quadrature(start, stop, order):
    x, w = leggauss(order)
    return start + (stop-start)*(x+1)/2, w*(stop-start)/2


def radial_basis(n, chi):
    """f_n=sin(chi)^2 C_(n-2)^2(cos chi)/sqrt(pi*(n^2-1)/8).

    Integral f_n f_m dchi = delta_nm. The TE dipole radial operator is
    -d_chi^2 + 2 csc(chi)^2, with eigenvalue n^2, n>=2.
    Return f/sin^2 and f'/sin, which stay regular at both poles.
    """
    n = np.asarray(n, int)
    chi = np.atleast_1d(chi)
    c, s = np.cos(chi), np.sin(chi)
    norm = np.sqrt(np.pi*(n*n-1)/8)
    h = np.array([eval_gegenbauer(j-2, 2, c) for j in n]) / norm[:, None]
    dh = np.array([4*eval_gegenbauer(j-3, 3, c) if j > 2 else np.zeros_like(c)
                   for j in n]) / norm[:, None]
    return h, 2*c*h-s*s*dh


@dataclass
class DipolePulse:
    """One exact TE l=1,m=0 sector; duality supplies the TM control.

    Initial radial potential is projected from sin^2(chi)*exp((cos chi-1)/w^2).
    Radiation is prepared data. It is NOT derived from charged-mouth birth.
    A common canonical phase rotates (n*q,p) while preserving total energy.
    """
    cutoff: int = 48
    width: float = .2
    phase: float = 0.
    axis: tuple = (0., 0., 1.)
    duality: float = 0.
    energy: float = 1.

    def __post_init__(self):
        if int(self.cutoff) != self.cutoff or self.cutoff < 2:
            raise ValueError("Integer cutoff >=2 required")
        if not 0 < self.width < 1 or not np.isfinite(self.energy) or self.energy <= 0:
            raise ValueError("Require 0<width<1 and finite positive energy")
        if not np.isfinite([self.phase, self.duality]).all():
            raise ValueError("Finite phases required")
        self.axis = unit(self.axis)
        if self.axis.shape != (3,):
            raise ValueError("The polarization axis has three components")
        self.n = np.arange(2, self.cutoff+1)
        chi, weights = quadrature(0, np.pi, max(128, 4*self.cutoff))
        h, _ = radial_basis(self.n, chi)
        f = h*np.sin(chi)**2
        q = f @ (weights*np.sin(chi)**2*np.exp((np.cos(chi)-1)/self.width**2))
        q *= np.sqrt(self.energy/(.5*np.sum((self.n*q)**2)))
        self.q0, self.p0 = q*np.cos(self.phase), -self.n*q*np.sin(self.phase)

    def state(self, time):
        c, s = np.cos(self.n*time), np.sin(self.n*time)
        return self.q0*c+self.p0*s/self.n, self.p0*c-self.n*self.q0*s

    def prepare(self, points):
        """Cache spatial factors for many time evaluations; points are in R4."""
        x = np.asarray(points, float)
        if x.ndim != 2 or x.shape[1] != 4 or not np.isfinite(x).all():
            raise ValueError("Require finite (N,4) points on unit S3")
        if np.max(abs(np.sum(x*x, axis=1)-1)) > 1e-10:
            raise ValueError("Points must lie on unit S3")
        s = np.linalg.norm(x[:, 1:], axis=1)
        c = x[:, 0]
        direction = np.divide(x[:, 1:], s[:, None], out=np.zeros_like(x[:, 1:]),
                              where=s[:, None] > 1e-14)
        direction[s <= 1e-14] = (0, 0, 1)
        chi = np.arctan2(s, c)
        h, derivative = radial_basis(self.n, chi)
        mu = direction @ self.axis
        radial = np.column_stack([-s, c[:, None]*direction])
        polar = np.column_stack([np.zeros(len(x)), self.axis-mu[:, None]*direction])
        azimuthal = np.column_stack([np.zeros(len(x)), np.cross(self.axis, direction)])
        return dict(h=h, derivative=derivative, s=s, mu=mu,
                    radial=radial, polar=polar, azimuthal=azimuthal)

    def fields(self, time, prepared):
        q, p = self.state(time)
        c = np.sqrt(3/(8*np.pi))
        alpha = -c*prepared['s']*(p @ prepared['h'])
        beta = 2*c*(q @ prepared['h'])
        gamma = c*(q @ prepared['derivative'])
        e = alpha[:, None]*prepared['azimuthal']
        b = (beta*prepared['mu'])[:, None]*prepared['radial'] + gamma[:, None]*prepared['polar']
        # Exact E cross B in the oriented orthonormal (chi,theta,phi) frame.
        momentum = ((alpha*gamma*(1-prepared['mu']**2))[:, None]*prepared['radial']
                    -(alpha*beta*prepared['mu'])[:, None]*prepared['polar'])
        ct, st = np.cos(self.duality), np.sin(self.duality)
        return ct*e+st*b, ct*b-st*e, momentum


@dataclass
class Worldtube:
    """Stationary closed S2 observation surface around an antipodal receiver.

    Offset moves the centre away from the south pole in the x0,x1 plane.
    r(u)=radius*(1+deformation*P2(u dot axis)) defines a shape perturbation.
    Neither parameter deforms the spacetime metric or installs a throat.
    Killing momentum uses xi(x)=h(c dot x)-c(h dot x), xi(c)=h.
    """
    radius: float = .3
    offset: float = 0.
    deformation: float = 0.
    angular_order: int = 24
    shape_axis: tuple = (0., 0., 1.)

    def __post_init__(self):
        if not (0 < self.radius < .7 and abs(self.offset) < .8 and abs(self.deformation) < .5):
            raise ValueError("Require 0<radius<.7, |offset|<.8, |deformation|<.5")
        if int(self.angular_order) != self.angular_order or self.angular_order < 4:
            raise ValueError("Integer angular order >=4 required")
        self.center = np.array([-np.cos(self.offset), np.sin(self.offset), 0., 0.])
        self.heading = np.array([np.sin(self.offset), np.cos(self.offset), 0., 0.])
        self.frame = np.array([self.heading, [0, 0, 1, 0], [0, 0, 0, 1]])
        mu, wm = leggauss(self.angular_order)
        phi = 2*np.pi*(np.arange(2*self.angular_order)+.5)/(2*self.angular_order)
        m, ph = np.meshgrid(mu, phi, indexing='ij')
        self.u = np.column_stack([(np.sqrt(1-m*m)*np.cos(ph)).ravel(),
                                  (np.sqrt(1-m*m)*np.sin(ph)).ravel(), m.ravel()])
        self.weights = np.repeat(wm, len(phi))*2*np.pi/len(phi)
        axis = unit(self.shape_axis)
        if axis.shape != (3,):
            raise ValueError("Shape axis has three components")
        z = self.u @ axis
        self.r = self.radius*(1+self.deformation*(3*z*z-1)/2)
        grad = 3*self.radius*self.deformation*z[:, None]*(axis-z[:, None]*self.u)
        direction = self.u @ self.frame
        self.points = np.cos(self.r)[:, None]*self.center+np.sin(self.r)[:, None]*direction
        radial = -np.sin(self.r)[:, None]*self.center+np.cos(self.r)[:, None]*direction
        stretch = np.sqrt(1+np.sum(grad*grad, axis=1)/np.sin(self.r)**2)
        self.normal = (radial-(grad @ self.frame)/np.sin(self.r)[:, None])/stretch[:, None]
        self.area = self.weights*np.sin(self.r)**2*stretch
        self.xi = self.killing(self.points)
        target = np.array([-1., 0., 0., 0.])
        distance = np.arccos(np.clip(self.center @ target, -1, 1))
        if distance < 1e-14:
            self.focus_inside = True
        else:
            direction_to_focus = (target-np.cos(distance)*self.center)/np.sin(distance)
            zfocus = (direction_to_focus @ self.frame.T) @ axis
            self.focus_inside = bool(distance < self.radius*(1+self.deformation*(3*zfocus*zfocus-1)/2))

    def killing(self, x):
        return (x @ self.center)[:, None]*self.heading-(x @ self.heading)[:, None]*self.center

    def volume(self, order=24):
        v, w = quadrature(0, 1, order)
        r = self.r[:, None]*v
        direction = self.u @ self.frame
        x = (np.cos(r)[:, :, None]*self.center
             +np.sin(r)[:, :, None]*direction[:, None, :]).reshape(-1, 4)
        weights = (self.weights[:, None]*self.r[:, None]*w*np.sin(r)**2).ravel()
        return x, weights


def surface_rates(pulse, tube, time, prepared=None):
    """Outward energy/Killing-momentum flux; positive radial stress is compression.

    incoming_characteristic is a diagnostic flux norm, not absorbed energy.
    Nonzero vacuum fields on both sides of a transparent surface exert no
    mechanical force. Flux plus stored-momentum change is the relevant audit.
    """
    if prepared is None:
        prepared = pulse.prepare(tube.points)
    e, b, s = pulse.fields(time, prepared)
    rho = .5*np.sum(e*e+b*b, axis=1)
    en = np.sum(e*tube.normal, axis=1)
    bn = np.sum(b*tube.normal, axis=1)
    sn = np.sum(s*tube.normal, axis=1)
    tnormal = rho[:, None]*tube.normal-en[:, None]*e-bn[:, None]*b
    tangential_square = 2*rho-en*en-bn*bn
    return dict(energy_flux=float(tube.area @ sn),
                momentum_flux=float(tube.area @ np.sum(tnormal*tube.xi, axis=1)),
                normal_stress=float(tube.area @ (rho-en*en-bn*bn)),
                incoming_characteristic=float(tube.area @ ((tangential_square-2*sn)/4)),
                outgoing_characteristic=float(tube.area @ ((tangential_square+2*sn)/4)))


def stored(pulse, tube, times, radial_order=24):
    x, w = tube.volume(radial_order)
    xi = tube.killing(x)
    # Avoid allocating a large (mode x volume) cache at high angular orders.
    result = np.zeros((len(times), 2))
    for start in range(0, len(x), 4096):
        end = start+4096
        prepared = pulse.prepare(x[start:end])
        for j, time in enumerate(times):
            e, b, s = pulse.fields(time, prepared)
            result[j, 0] += w[start:end] @ (.5*np.sum(e*e+b*b, axis=1))
            result[j, 1] += w[start:end] @ np.sum(s*xi[start:end], axis=1)
    return result


def encounter(pulse, tube, start, stop, time_order=96, radial_order=24):
    if not np.isfinite([start, stop]).all() or stop <= start:
        raise ValueError("Finite increasing observation window required")
    times, weights = quadrature(start, stop, time_order)
    prepared = pulse.prepare(tube.points)
    rates = [surface_rates(pulse, tube, t, prepared) for t in times]
    integral = {key: float(weights @ np.array([r[key] for r in rates])) for key in rates[0]}
    endpoints = stored(pulse, tube, [start, stop], radial_order)
    delta_energy, delta_momentum = endpoints[1]-endpoints[0]
    # J_matter = -Delta P_field - integral outward T.n; exactly zero in vacuum.
    return dict(start=float(start), stop=float(stop), integrals=integral,
                stored_initial=endpoints[0].tolist(), stored_final=endpoints[1].tolist(),
                energy_balance_residual=float(delta_energy+integral['energy_flux']),
                transparent_matter_impulse=float(-delta_momentum-integral['momentum_flux']),
                incident_only_impulse=float(-integral['momentum_flux']),
                focus_inside=tube.focus_inside,
                scope='Vacuum observation worldtube, no absorption/reflection law or mouth recoil')


class ConductingReceiver:
    """Centered, fixed perfect-conductor S2 boundary: a reflection comparator.

    Exterior domain 0<chi<pi-radius. TE radial f vanishes at both ends.
    The operator -d_chi^2+2/sin(chi)^2 is discretized by a symmetric second-order
    finite difference matrix and evolved by exact discrete normal modes.
    Tangential force cancels by symmetry. The radius load -f'(boundary)^2/2
    follows from the exterior Maxwell stress, not from an inserted force law.
    This is an ideal conducting shell, not a Klein throat or an evolving mouth.
    """
    def __init__(self, radius=.3, width=.2, cells=512):
        if not 0 < radius < .7 or not 0 < width < 1 or int(cells) != cells or cells < 32:
            raise ValueError('Require 0<radius<.7, 0<width<1, integer cells>=32')
        self.radius, self.width, self.cells = radius, width, cells
        self.h = (np.pi-radius)/(cells+1)
        self.chi = self.h*np.arange(1, cells+1)
        self.diagonal = 2/self.h**2+2/np.sin(self.chi)**2
        self.off = np.full(cells-1, -1/self.h**2)
        eigenvalues, self.modes = eigh_tridiagonal(self.diagonal, self.off)
        self.omega = np.sqrt(eigenvalues)
        initial = np.sin(self.chi)**2*np.exp((np.cos(self.chi)-1)/width**2)
        coefficients = self.modes.T @ (initial*np.sqrt(self.h))
        self.coefficients = coefficients/np.sqrt(.5*np.sum(eigenvalues*coefficients**2))
        self.boundary_derivative = (-4*self.modes[-1]+self.modes[-2])/(2*self.h*np.sqrt(self.h))

    def radius_force(self, times):
        times = np.atleast_1d(times)
        slope = (np.cos(times[:, None]*self.omega)*self.coefficients) @ self.boundary_derivative
        return -.5*slope*slope

    def encounter(self, start, stop, order=192):
        times, weights = quadrature(start, stop, order)
        forces = self.radius_force(times)
        return dict(radius=self.radius, width=self.width, cells=self.cells,
                    start=float(start), stop=float(stop), quadrature_order=order,
                    radius_impulse=float(weights @ forces),
                    peak_inward_radius_force=float(-np.min(forces)),
                    prepared_energy=float(.5*np.sum((self.omega*self.coefficients)**2)),
                    vector_impulse=0., vector_impulse_reason='Exact angular symmetry',
                    scope='Fixed ideal conductor reflection; generalized radius load, not translational attraction')
