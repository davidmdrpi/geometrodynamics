"""Autonomous scalar port model; not a glued moving-mouth GR solution.

State order is (u, p, q_A, q_B, r, P_A, P_B, P_r, clock reaction).
The last coordinate integrates the reaction from the prescribed clock and
static constraints. It does not feed back on the dynamics.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import cho_factor, cho_solve
from scipy.special import sph_harm_y


def real_basis(lmax, theta, phi, derivatives=False):
    """Orthonormal real harmonics: m=0, then sqrt(2) Re/Im for m>0."""
    theta, phi = np.broadcast_arrays(theta, phi)
    values, dt, dp, ell = [], [], [], []
    for l in range(lmax + 1):
        for m in range(l + 1):
            if not derivatives:
                y = sph_harm_y(l, m, theta, phi)
                values.extend([y.real] if m == 0 else [np.sqrt(2)*y.real, np.sqrt(2)*y.imag])
                continue
            y, dy = sph_harm_y(l, m, theta, phi, diff_n=1)
            parts = [(y.real, dy[..., 0].real, dy[..., 1].real)] if m == 0 else [
                (np.sqrt(2)*y.real, np.sqrt(2)*dy[..., 0].real, np.sqrt(2)*dy[..., 1].real),
                (np.sqrt(2)*y.imag, np.sqrt(2)*dy[..., 0].imag, np.sqrt(2)*dy[..., 1].imag)]
            for a, b, c in parts:
                values.append(a); dt.append(b); dp.append(c); ell.append(l)
    y = np.stack(values, axis=-1)
    if derivatives:
        return y, np.stack(dt, axis=-1), np.stack(dp, axis=-1)
    return y


@dataclass(frozen=True)
class Parameters:
    lmax: int = 10
    nodes: int = 17
    sigma: float = .30
    footprint: float = .20
    length: float = .565
    radius0: float = float(np.sin(.3)**2)
    mass: float = .3
    radial_mass: float = .2
    stiffness: float = 10.
    coupling: float = 1.
    energy: float = .08
    qa: float = .35
    qb: float = float(np.pi-.35)
    source: float = 1.1
    offset: float = -1.2
    fixed_mouths: bool = False
    fixed_radius: bool = False


class CoupledMouthThroat:
    def __init__(self, parameters=Parameters()):
        self.par = p = parameters
        if (p.lmax < 1 or p.nodes < 3 or min(p.length, p.radius0, p.mass,
                p.radial_mass, p.sigma, p.footprint) <= 0 or
                min(p.energy, p.coupling, p.stiffness) < 0):
            raise ValueError('Invalid resolution, scale or energy')
        self.kappa = -p.offset/(np.cos(p.qa)-np.cos(p.qb))
        if abs(self.kappa) >= .9:
            raise ValueError('Initial clock tilt exceeds the registered margin')
        self.ns = (p.lmax+1)**2
        self.n = self.ns+p.nodes
        self.ell = np.concatenate([np.full(2*l+1, l) for l in range(p.lmax+1)])
        self.ks = (self.ell+.5)**2
        self.smoothing = np.exp(-.5*p.footprint**2*self.ell*(self.ell+1))
        z, wt = np.polynomial.legendre.leggauss(2*p.lmax+5)
        phi = np.arange(4*p.lmax+8)*2*np.pi/(4*p.lmax+8)
        theta, ph = np.broadcast_arrays(np.arccos(z)[:, None], phi[None, :])
        weights = np.broadcast_to(wt[:, None]*2*np.pi/len(phi), theta.shape).ravel()
        y, yt, yp = real_basis(p.lmax, theta, ph, True)
        y, yt, yp = [a.reshape(-1, self.ns) for a in (y, yt, yp)]
        st, ct = np.sin(theta).ravel(), np.cos(theta).ravel()
        self.ms = y.T @ ((weights*(1-self.kappa**2*st**2))[:, None]*y)
        self.bs = y.T @ ((-weights*self.kappa*st)[:, None]*yt)
        self.ms = (self.ms+self.ms.T)/2
        self.factor = cho_factor(self.ms)
        dy = np.cos(ph).ravel()[:, None]*yt-(ct/st*np.sin(ph).ravel())[:, None]*yp
        d = y.T @ (weights[:, None]*dy)
        self.d = (d-d.T)/2
        self.mb = self.ms@self.d-self.d@self.ms
        self.bb = self.bs@self.d-self.d@self.bs
        h = 1/(p.nodes-1)
        self.mn = np.full(p.nodes, h); self.mn[[0, -1]] *= .5
        self.bn = np.zeros((p.nodes, p.nodes))
        self.kn = np.zeros_like(self.bn)
        for j in range(p.nodes-1):
            self.bn[j:j+2, j:j+2] += np.array([[-1., 1.], [-1., 1.]])/2
            self.kn[j:j+2, j:j+2] += np.array([[1., -1.], [-1., 1.]])/h

    def unpack(self, state):
        return state[:self.n], state[self.n:2*self.n], state[2*self.n:2*self.n+3], state[2*self.n+3:2*self.n+6]

    def neck(self, q):
        p = self.par
        qa, qb, r = q
        a = self.kappa*(np.cos(qb)-np.cos(qa))-p.offset
        aq = self.kappa*np.array([np.sin(qa), -np.sin(qb)])
        w = (r*r+(.1*p.radius0)**2)/(1.01*p.radius0**2)
        wr = 2*r/(1.01*p.radius0**2)
        m = w*p.length*(1-(a/p.length)**2)*self.mn
        return a, aq, w, wr, m

    def ports(self, q):
        # Positions are on the x-z great circle. Using the rotation generator
        # avoids the polar-coordinate derivative ambiguity at either pole.
        theta = np.arccos(np.cos(q[:2]))
        phi = np.where(np.sin(q[:2]) >= 0, 0., np.pi)
        f = real_basis(self.par.lmax, theta, phi)*self.smoothing
        return f, f@self.d

    def evaluate(self, state, diagnostics=False):
        p = self.par
        u, mom, q, mech = self.unpack(state)
        us, un = u[:self.ns], u[self.ns:]
        ps, pn = mom[:self.ns], mom[self.ns:]
        a, aq, w, wr, mn = self.neck(q)
        if np.min(mn) <= 0:
            raise ValueError('Non-spacelike neck clock slice')
        vs = cho_solve(self.factor, ps-self.bs@us)
        vn = (pn-w*a*(self.bn@un))/mn
        f, fq = self.ports(q)
        gap = un[[0, -1]]-f@us
        g = p.coupling*(q[2]/p.radius0)**2
        gr = 2*p.coupling*q[2]/p.radius0**2
        pdots = self.bs.T@vs-self.ks*us+g*(f.T@gap)
        pdotn = w*a*(self.bn.T@vn)-w/p.length*(self.kn@un)
        pdotn[[0, -1]] -= g*gap
        fa = -w*a/p.length*np.dot(vn*self.mn, vn)+w*np.dot(vn, self.bn@un)
        forces = np.empty(3)
        forces[:2] = fa*aq+g*gap*(fq@us)
        forces[2] = (.5*wr*p.length*(1-(a/p.length)**2)*np.dot(vn*self.mn, vn)
                     +wr*a*np.dot(vn, self.bn@un)
                     -.5*wr/p.length*np.dot(un, self.kn@un)
                     -.5*gr*np.dot(gap, gap)-p.stiffness*(q[2]-p.radius0))
        qdot = mech/np.array([p.mass, p.mass, p.radial_mass])
        clock_torque = (-.5*vs@self.mb@vs-vs@self.bb@us+fa*aq.sum())
        if p.fixed_mouths:
            clock_torque -= forces[:2].sum()
            forces[:2] = 0.; qdot[:2] = 0.
        if p.fixed_radius:
            forces[2] = 0.; qdot[2] = 0.
        rhs = np.r_[vs, vn, pdots, pdotn, qdot, forces, -clock_torque]
        if not diagnostics:
            return rhs
        energies = dict(sphere=.5*(vs@self.ms@vs+np.dot(us*self.ks, us)),
            throat=.5*(np.dot(vn*mn, vn)+w/p.length*(un@self.kn@un)),
            port=.5*g*np.dot(gap, gap), mouth_kinetic=np.dot(mech[:2], mech[:2])/(2*p.mass),
            radial_kinetic=mech[2]**2/(2*p.radial_mass),
            neck_potential=.5*p.stiffness*(q[2]-p.radius0)**2)
        return rhs, energies

    def initial(self, dephased=False):
        p = self.par
        state = np.zeros(2*self.n+7)
        state[2*self.n:2*self.n+3] = [p.qa, p.qb, p.radius0]
        lam = self.ell*(self.ell+1)
        us = lam*np.exp(-.5*p.sigma**2*lam)*real_basis(p.lmax, p.source, 0.)
        vs = .7*(self.d@us)
        if dephased:
            # Fixed mode phases; a control, not a search over timing.
            phases = np.array([0.0, -0.6556919081789173, 1.8800633058144132, -1.9943377210282789, 1.6258161186045408, 3.042546614319629, -2.9666243919137143, 1.3529129533658857, 0.8540670043415368, 2.693507565988488, 1.5256417717610584])
            if p.lmax > 10:
                raise ValueError("Recovered dephased control is specified only through l=10")
            phase = phases[self.ell]
            us, vs = (us*np.cos(phase)+vs/(self.ell+.5)*np.sin(phase),
                      vs*np.cos(phase)-us*(self.ell+.5)*np.sin(phase))
        state[:self.ns] = us
        state[self.n:self.n+self.ns] = self.ms@vs+self.bs@us
        energy = sum(self.evaluate(state, True)[1].values())
        state[:2*self.n] *= np.sqrt(p.energy/energy)
        return state

    def angular(self, state):
        u, mom, _, mech = self.unpack(state)
        return mech[:2].sum()-mom[:self.ns]@self.d@u[:self.ns]

    def evolve(self, horizon=8*np.pi, rtol=2e-8, initial=None, dephased=False):
        state = self.initial(dephased) if initial is None else np.array(initial, float)
        p = self.par
        def cutoff(t, y):
            return y[2*self.n+2]-.25*p.radius0
        def chart(t, y):
            a = self.neck(y[2*self.n:2*self.n+3])[0]
            return .9-max(abs(self.kappa), abs(a/p.length))
        def overlap(t, y):
            q = y[2*self.n:2*self.n+2]
            return np.arccos(np.cos(q[1]-q[0]))-.3
        for event in (cutoff, chart, overlap):
            event.terminal = True; event.direction = -1
        sol = solve_ivp(lambda t, y: self.evaluate(y), (0., horizon), state,
            method='DOP853', rtol=rtol, atol=rtol*.01, max_step=.05,
            dense_output=True, events=(cutoff, chart, overlap))
        if not sol.success:
            raise RuntimeError(sol.message)
        times = np.unique(np.r_[np.arange(0., sol.t[-1], .04),
            np.arange(0., sol.t[-1], 2*np.pi), sol.t[-1]])
        return times, sol.sol(times).T, next((name for name, ts in zip(
            ('radius_cutoff', 'clock_margin', 'mouth_overlap'), sol.t_events) if len(ts)), 'horizon')
