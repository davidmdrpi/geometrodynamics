"""Pair birth, separation, return and contact for two Wheeler-flux mouths on the round S².

Two independently moving mouths carry equal and opposite flux q through one handle:
Wheeler's charge without charge. On the sphere side Gauss's law makes the two
charges opposite whatever the gluing, so the handle's orientation does not enter.

**Field.** The sphere-side field is 2+1 Maxwell on the unit sphere, in Coulomb gauge.
- The longitudinal part is fixed instantaneously by the mouth positions.
- The transverse part is A = n x grad chi, carried by modes with omega_l = sqrt(l(l+1)).
- In 2+1 dimensions this is one scalar wave field psi, with E = n x grad psi and
  B = psi_t, that winds by q around each rim. The mouths are its vortex and
  antivortex.

**Mouths.** The mouths are relativistic point particles of rest energy M, with
Gaussian footprints of width f, as #340's.

**Action.** Everything follows from one Lagrangian:

    L = sum_k [-M sqrt(1 - v_k^2) + q_k v_k . A(x_k)] - U(x_A, x_B)
        + 1/2 sum_lm l(l+1) (|c_lm'|^2 - l(l+1) |c_lm|^2),

where chi = sum c_lm Y_lm and U = (1/2) int |grad phi|^2 is the Coulomb energy.
- **Forces are two-sided.** Each force on a mouth is reaction to a change of field
  momentum, so total energy and total canonical angular momentum are conserved.
- **The field is split by source.** It is kept as one copy per mouth, each driven
  by that mouth's current alone. The force on A then splits exactly into:
  - the field B made: Coulomb plus transverse, the retarded exchange;
  - A's own transverse field: local reaction, plus A's own emission returning
    after antipodal refocusing.
- **Regularity law.** The glued surface is regular while the two rim discs are
  disjoint, d > 2a.
  - When the dynamics brings the rims into contact, the handle can no longer be
    glued. The mouths, their neck and their Coulomb field are then removed.
  - Their energy is released as a pulse of the conformal deformation field
    (#339/#344 operator, omega = l + 1/2) at the contact point.
  - The contact time is an outcome of the dynamics. The surgery itself is a
    stated rule, not something the dynamics derives.

**Scope.** Classical throughout: nothing is quantized. The neck shape and the
rest energy M (neck surface energy plus the flux energy inside the neck) are
held fixed.
"""
from dataclasses import dataclass

import numpy as np
from numpy.polynomial import legendre
from scipy import optimize, sparse

from .moving_mouth_pairs import harmonics


@dataclass(frozen=True)
class Parameters:
    cutoff: int = 32
    footprint: float = .12       # Gaussian smoothing of each mouth's charge and current (#340's footprint)
    charge: float = .75          # flux through the handle; mouth A carries +q, mouth B -q
    mass: float = 1.             # rest energy of each mouth, held fixed
    mouth_angle: float = .3      # rim radius a (#338); the surface is regular while d > 2a
    release_width: float = .25
    transverse: bool = True      # False: instantaneous Coulomb only, no waves


def rotate(v, axis, angle):
    """Rotate vectors v (..., 3) about unit axes (..., 3) by angles (...)."""
    angle = np.asarray(angle)[..., None]
    co, si = np.cos(angle), np.sin(angle)
    return v*co + np.cross(axis, v)*si + axis*np.sum(axis*v, -1, keepdims=True)*(1-co)


def birth_positions(a, v, plane=((0., 0., 1.), (1., 0., 0.))):
    """Mouths touching at P (rims in contact) and moving apart along one great circle.

    Returns (x, velocity, P, tangent): A at angle +a from P towards the tangent,
    B at -a, each with speed v away from P.
    """
    P, T = (np.asarray(e, float) for e in plane)
    x = np.array([np.cos(a)*P + np.sin(a)*T, np.cos(a)*P - np.sin(a)*T])
    vel = v*np.array([-np.sin(a)*P + np.cos(a)*T, -np.sin(a)*P - np.cos(a)*T])
    return x, vel, P, T


def rotation_generators(pairs):
    """D_i on Y_lm coefficients, for the listed (l, m): D_i c are the coefficients of
    (e_i x n).grad of the field, so e_i.(n x grad f) = Y . D_i c."""
    index = {lm: j for j, lm in enumerate(pairs)}
    rows, cols, vals = [], [], []
    for j, (l, m) in enumerate(pairs):
        if m < l:
            rows.append(index[(l, m+1)])
            cols.append(j)
            vals.append(np.sqrt((l-m)*(l+m+1)))
    n = len(pairs)
    plus = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
    minus = plus.T.tocsr()
    order = np.array([m for _, m in pairs])
    return 0.5j*(plus+minus), 0.5*(plus-minus), sparse.diags(1j*order, format='csr')


@dataclass
class PairState:
    t: float
    x: np.ndarray            # (2, 3) mouth positions
    L: np.ndarray            # (2, 3) mechanical angular momenta x cross p
    Q: np.ndarray            # (2, n) transverse modes, one copy per source mouth
    P: np.ndarray            # (2, n)
    alive: bool = True
    release: dict = None     # set at contact: time, point, energy and pulse coefficients

    def copy(self):
        return PairState(self.t, self.x.copy(), self.L.copy(), self.Q.copy(), self.P.copy(), self.alive,
                         None if self.release is None else dict(self.release))


class PairBirth:
    def __init__(self, parameters=Parameters()):
        p = self.par = parameters
        if (p.cutoff < 2 or p.footprint <= 0 or p.mass <= 0 or not 0 < p.mouth_angle < np.pi/4
                or p.release_width <= 0):
            raise ValueError('Positive footprint, mass and release width, and 0 < a < pi/4 required')
        pairs = [(l, m) for l in range(1, p.cutoff+1) for m in range(-l, l+1)]
        self.ell, self.order = np.array(pairs).T
        self.lam = self.ell*(self.ell+1.)
        self.omega = np.sqrt(self.lam)
        self.smoothing = np.exp(-p.footprint**2*self.lam/2)
        self.n = len(pairs)
        self.generators = rotation_generators(pairs)
        self._stacked = sparse.vstack(self.generators).tocsr()
        self._cache = (None, None)
        self.charges = np.array([p.charge, -p.charge])
        # Coulomb energy U = 1/2 sum_jk q_j q_k S(x_j . x_k), S = sum_l s_l^2 (2l+1) P_l / (4 pi l(l+1)).
        ls = np.arange(p.cutoff+1)
        coef = np.zeros(p.cutoff+1)
        coef[1:] = np.exp(-p.footprint**2*ls[1:]*(ls[1:]+1))*(2*ls[1:]+1)/(4*np.pi*ls[1:]*(ls[1:]+1))
        self.S, self.dS = coef, legendre.legder(coef)
        # Deformation field for the release: conformal operator -Delta + 1/4 (exact antipodal refocusing).
        upairs = [(l, m) for l in range(p.cutoff+1) for m in range(-l, l+1)]
        self.u_ell, self.u_order = np.array(upairs).T
        self.u_omega = self.u_ell + .5
        self.u_generators = rotation_generators(upairs)

    # ---- field pieces -------------------------------------------------------------------------
    def footprints(self, x):
        key = np.asarray(x).tobytes()
        if self._cache[0] != key:
            self._cache = (key, self.smoothing*harmonics(self.ell, self.order, x).conj())      # (2, n)
        return self._cache[1]

    def _derivatives(self, c):
        """D_i c for i = 0, 1, 2, shape (3, n)."""
        return (self._stacked @ c).reshape(3, self.n)

    def _rot(self, w, c):
        """Ambient vectors n x grad chi at each mouth for coefficients c, shape (2, 3)."""
        return (w.conj() @ self._derivatives(c).T).real

    def potentials(self, s, w=None):
        """Per source copy k: A[k, j], E_T[k, j] (ambient vectors) and B[k, j] at mouth j."""
        w = self.footprints(s.x) if w is None else w
        A = np.stack([self._rot(w, s.Q[k]/self.omega) for k in range(2)])
        E = np.stack([-self._rot(w, s.P[k]/self.omega) for k in range(2)])
        B = np.stack([-(w.conj() @ (self.omega*s.Q[k])).real for k in range(2)])
        return A, E, B

    def source(self, x, vel, w=None):
        """Generalised force on each copy's modes from its own mouth's current: u_k / sqrt(lambda)."""
        w = self.footprints(x) if w is None else w
        out = np.empty((2, self.n), complex)
        for k in range(2):
            out[k] = -self.charges[k]*(vel[k] @ self._derivatives(w[k]))/self.omega
        return out

    def coulomb(self, x):
        """Coulomb energy (self plus interaction) and forces on the two mouths."""
        mu = float(np.clip(x[0] @ x[1], -1, 1))
        qa, qb = self.charges
        energy = .5*(qa**2+qb**2)*legendre.legval(1., self.S) + qa*qb*legendre.legval(mu, self.S)
        g = -qa*qb*legendre.legval(mu, self.dS)
        force = np.array([g*(x[1]-mu*x[0]), g*(x[0]-mu*x[1])])
        return energy, force

    def coulomb_interaction(self, d):
        """Interaction part of U at separation d (vectorised); the self part is constant."""
        return self.charges[0]*self.charges[1]*legendre.legval(np.cos(d), self.S)

    def scalar_interaction(self, d):
        """Static interaction of the same two charges as SOURCES of the conformal scalar
        (-Delta + 1/4, #339): U = -q_A q_B G(d). Scalar exchange pulls like charges together
        and pushes a source-sink pair apart, the opposite sign to the flux mouths."""
        ls = np.arange(self.par.cutoff+1)
        coef = np.exp(-self.par.footprint**2*ls*(ls+1))*(2*ls+1)/(4*np.pi*(ls+.5)**2)
        return -self.charges[0]*self.charges[1]*legendre.legval(np.cos(d), coef)

    # ---- mechanics -----------------------------------------------------------------------------
    def gamma(self, L):
        return np.sqrt(1+np.sum(L*L, -1)/self.par.mass**2)

    def velocity(self, s):
        return np.cross(s.L, s.x)/(self.par.mass*self.gamma(s.L))[:, None]

    def separation(self, x):
        return float(np.arccos(np.clip(x[0] @ x[1], -1, 1)))

    def separation_rate(self, s):
        v = self.velocity(s)
        d = self.separation(s.x)
        return -(v[0] @ s.x[1] + v[1] @ s.x[0])/max(np.sin(d), 1e-300)

    def initial(self, v=1/3, dressed=False, plane=((0., 0., 1.), (1., 0., 0.)), velocities=None):
        """Mouths in rim contact at P, kicked apart with speed v each.

        - Default: the field is the static Coulomb field of the touching pair and the
          transverse field is zero. The sudden start of the current then radiates a birth
          flash. By linearity the flash is exactly the free wave launched by -Q_qs, where
          Q_qs is the quasi-static (magnetostatic) field of the initial currents
          (see `flash`).
        - dressed=True: the transverse field starts as Q_qs. This removes the flash and is
          the no-flash control.
        - velocities overrides the symmetric start (independent mouths).
        """
        x, vel, _, _ = birth_positions(self.par.mouth_angle, v, plane)
        if velocities is not None:
            vel = np.asarray(velocities, float)
            vel = vel - np.sum(vel*x, -1, keepdims=True)*x
        speed = np.linalg.norm(vel, axis=1)
        if np.any(speed >= 1):
            raise ValueError('Mouth speeds must be below 1')
        p = self.par.mass*vel/np.sqrt(1-speed**2)[:, None]
        Q = np.zeros((2, self.n), complex)
        if dressed and self.par.transverse:
            Q = self.quasi_static(x, vel)
        return PairState(0., x, np.cross(x, p), Q, np.zeros((2, self.n), complex))

    def quasi_static(self, x, vel):
        """Magnetostatic transverse field of the currents q_k v_k at x_k, one copy per mouth."""
        return self.source(x, vel)/self.lam

    def free_wave(self, Q, P, tau):
        """Free transverse evolution of (Q, P) by tau (any sign)."""
        co, si = np.cos(self.omega*tau), np.sin(self.omega*tau)
        return Q*co + P*si/self.omega, P*co - self.omega*Q*si

    def flash(self, s0):
        """The birth flash of a default start: the free wave launched by -Q_qs at t = 0."""
        return -self.quasi_static(s0.x, np.cross(s0.L, s0.x)/(self.par.mass*self.gamma(s0.L))[:, None]), \
            np.zeros((2, self.n), complex)

    def free_wave_force(self, history, Q0, P0, t0=0.):
        """Force of a free transverse wave (Q0, P0 at time t0, one copy per source) on each
        mouth along a recorded trajectory: shape (time, 2 mouths, 2 sources, 3)."""
        out = np.zeros((len(history['t']), 2, 2, 3))
        g = np.sqrt(1 + np.sum(history['L']**2, -1)/self.par.mass**2)
        for i, (t, x, L) in enumerate(zip(history['t'], history['x'], history['L'])):
            if history['energy_parts']['rest'][i] == 0:
                continue
            Q, P = self.free_wave(Q0, P0, t - t0)
            w = self.footprints(x)
            v = np.cross(L, x)/(self.par.mass*g[i])[:, None]
            for k in range(2):
                E = -self._rot(w, P[k]/self.omega)
                B = -(w.conj() @ (self.omega*Q[k])).real
                for j in range(2):
                    out[i, j, k] = self.charges[j]*(E[j] + B[j]*np.cross(v[j], x[j]))
        return out

    # ---- ledgers -------------------------------------------------------------------------------
    def transverse_energy(self, s):
        Q, P = s.Q.sum(0), s.P.sum(0)
        return float(.5*(np.vdot(P, P).real + np.vdot(self.omega*Q, self.omega*Q).real))

    def deformation_energy(self, s):
        if s.release is None:
            return 0.
        a, p = s.release['a0'], s.release['p0']
        return float(.5*(np.vdot(p, p).real + np.vdot(self.u_omega*a, self.u_omega*a).real))

    def deformation_angular(self, s):
        if s.release is None:
            return np.zeros(3)
        a, p = s.release['a0'], s.release['p0']
        return -np.array([np.vdot(p, g @ a).real for g in self.u_generators])

    def energy_components(self, s):
        out = dict(rest=0., kinetic=0., coulomb=0., transverse=self.transverse_energy(s),
                   deformation=self.deformation_energy(s))
        if s.alive:
            g = self.gamma(s.L)
            out.update(rest=2*self.par.mass, kinetic=float(self.par.mass*np.sum(g-1)),
                       coulomb=float(self.coulomb(s.x)[0]))
        return out

    def field_angular(self, Q, P):
        return -np.array([np.vdot(P, g @ Q).real for g in self.generators])

    def angular_momentum(self, s):
        """Total canonical angular momentum: mouths x cross (p + qA) plus field."""
        J = self.field_angular(s.Q.sum(0), s.P.sum(0)) + self.deformation_angular(s)
        if s.alive:
            A, _, _ = self.potentials(s)
            p = np.cross(s.L, s.x) + self.charges[:, None]*A.sum(0)
            J = J + np.cross(s.x, p).sum(0)
        return J

    def forces(self, s):
        """Instantaneous forces on each mouth, by source: Coulomb from the partner,
        transverse from the partner's copy, transverse from the mouth's own copy."""
        if not s.alive:
            return dict(coulomb=np.zeros((2, 3)), exchange=np.zeros((2, 3)), own=np.zeros((2, 3)))
        _, fc = self.coulomb(s.x)
        out = dict(coulomb=fc, exchange=np.zeros((2, 3)), own=np.zeros((2, 3)))
        if self.par.transverse:
            _, E, B = self.potentials(s)
            v = self.velocity(s)
            for j in range(2):
                for k in range(2):
                    f = self.charges[j]*(E[k, j] + B[k, j]*np.cross(v[j], s.x[j]))
                    out['own' if k == j else 'exchange'][j] += f
        return out

    # ---- integrator ----------------------------------------------------------------------------
    def _field_flow(self, s, tau):
        """Free transverse field for tau with mouths frozen: p_mech changes by -q dA (E_T kick)."""
        if not self.par.transverse:
            return
        w = self.footprints(s.x) if s.alive else None
        if s.alive:
            A0 = self._rot(w, s.Q.sum(0)/self.omega)
        co, si = np.cos(self.omega*tau), np.sin(self.omega*tau)
        s.Q, s.P = s.Q*co + s.P*si/self.omega, s.P*co - self.omega*s.Q*si
        if s.alive:
            A1 = self._rot(w, s.Q.sum(0)/self.omega)
            s.L -= self.charges[:, None]*np.cross(s.x, A1-A0)

    def _magnetic(self, s, tau):
        """Rotation of each mouth's momentum about its normal by the magnetic field."""
        if not self.par.transverse:
            return
        _, _, B = self.potentials(s)
        angle = -self.charges*B.sum(0)*tau/(self.par.mass*self.gamma(s.L))
        s.L = rotate(s.L, s.x, angle)

    def _coulomb_kick(self, s, tau):
        s.L += tau*np.cross(s.x, self.coulomb(s.x)[1])

    def _particle_flow(self, s, dt):
        """Mouths in the frozen transverse field and Coulomb force; the field modes take
        the mouths' currents as a source (trapezoid along the drift)."""
        self._coulomb_kick(s, dt/2)
        self._magnetic(s, dt/2)
        if self.par.transverse:
            s.P += dt/2*self.source(s.x, self.velocity(s))
        length = np.linalg.norm(s.L, axis=1)
        g = self.gamma(s.L)
        for j in range(2):
            if length[j] > 0:
                s.x[j] = rotate(s.x[j], s.L[j]/length[j], length[j]*dt/(self.par.mass*g[j]))
        s.x /= np.linalg.norm(s.x, axis=1, keepdims=True)
        if self.par.transverse:
            s.P += dt/2*self.source(s.x, self.velocity(s))
        self._magnetic(s, dt/2)
        self._coulomb_kick(s, dt/2)

    def step(self, s, dt):
        if s.alive:
            self._field_flow(s, dt/2)
            self._particle_flow(s, dt)
            self._field_flow(s, dt/2)
        else:
            self._field_flow(s, dt)
        s.t += dt

    # ---- regularity law and release ------------------------------------------------------------
    def contact(self, s):
        return s.alive and self.separation(s.x) <= 2*self.par.mouth_angle and self.separation_rate(s) < 0

    def release(self, s):
        """Remove the mouths at rim contact and release their energy and angular momentum
        as a deformation pulse at the contact point c.

        - Energy released: rest + kinetic + the whole Coulomb (longitudinal) field energy,
          which must vanish with the charges. The transverse field continues freely.
        - Pulse: u_t = alpha s + sum_i gamma_i D_i s and u = sum_i beta_i D_i s, with
          s = exp(-w^2 l(l+1)/2) Y*(c). The beta and gamma make the pulse carry the mouths'
          canonical angular momentum (gamma is needed for the spin about c, which an
          axisymmetric shape alone cannot hold), and alpha then fixes the energy. They are
          solved by least squares with the smallest beta and gamma.
        - Returns the angular momentum left unaccounted (zero to solver tolerance).
        """
        e = self.energy_components(s)
        released = e['rest'] + e['kinetic'] + e['coulomb']
        J_before = self.angular_momentum(s)
        c = s.x[0] + s.x[1]
        c /= np.linalg.norm(c)
        s.alive = False
        J_mouths = J_before - self.angular_momentum(s)
        shape = np.exp(-self.par.release_width**2*self.u_ell*(self.u_ell+1)/2)*harmonics(
            self.u_ell, self.u_order, c).conj()
        Ds = np.array([g @ shape for g in self.u_generators])

        def pulse(z):
            alpha, beta, gamma = z[0], z[1:4], z[4:7]
            return beta @ Ds, alpha*shape + gamma @ Ds

        def residual(z):
            a, p = pulse(z)
            J = -np.array([np.vdot(p, g @ a).real for g in self.u_generators])
            E = .5*(np.vdot(p, p).real + np.vdot(self.u_omega*a, self.u_omega*a).real)
            return np.r_[J - J_mouths, (E - released)/max(released, 1e-12), 1e-6*z[1:]]

        z0 = np.r_[np.sqrt(2*released/np.vdot(shape, shape).real), np.zeros(6)]
        z = optimize.least_squares(residual, z0, xtol=1e-15, ftol=1e-15, gtol=1e-15).x
        a0, p0 = pulse(z)
        s.release = dict(time=s.t, point=c, energy=released, p0=p0, a0=a0)
        return J_before - self.angular_momentum(s)

    def deformation(self, s, points, t=None):
        """Deformation field u and u_t at points (unit vectors) at time t (default s.t)."""
        if s.release is None:
            z = np.zeros(np.shape(points)[:-1])
            return z, z
        tau = (s.t if t is None else t) - s.release['time']
        Y = harmonics(self.u_ell, self.u_order, points)
        p0 = s.release['p0']
        a0 = s.release['a0']
        co, si = np.cos(self.u_omega*tau), np.sin(self.u_omega*tau)
        u = (Y @ (a0*co + p0*si/self.u_omega)).real
        ut = (Y @ (p0*co - self.u_omega*a0*si)).real
        return u, ut

    def field_density(self, s, Y, uY=None):
        """Energy densities at points with harmonics Y = harmonics(ell, order, points)
        (and uY for the deformation modes): transverse (E_T^2 + B^2)/2 and deformation
        (u_t^2 + |grad u|^2 + u^2/4)/2."""
        Q, P = s.Q.sum(0), s.P.sum(0)
        E = np.stack([(Y @ d).real for d in self._derivatives(P/self.omega)])
        B = (Y @ (self.omega*Q)).real
        transverse = .5*(np.sum(E**2, 0) + B**2)
        deformation = np.zeros_like(B)
        if s.release is not None and uY is not None:
            tau = s.t - s.release['time']
            co, si = np.cos(self.u_omega*tau), np.sin(self.u_omega*tau)
            a = s.release['a0']*co + s.release['p0']*si/self.u_omega
            p = s.release['p0']*co - self.u_omega*s.release['a0']*si
            grad = np.stack([(uY @ (g @ a)).real for g in self.u_generators])
            u = (uY @ a).real
            deformation = .5*((uY @ p).real**2 + np.sum(grad**2, 0) + .25*u**2)
        return transverse, deformation

    # ---- runs ----------------------------------------------------------------------------------
    def run(self, s, end, dt=.004, output_dt=.02, refine_contact=True, field_points=None, field_dt=None,
            retime=None):
        """Evolve to `end`, recording ledgers and the forces on the mouths by source.

        At rim contact the step is refined by secant iteration to the contact separation 2a,
        and the surgery rule is applied there.

        retime=(t_d, shift) is the timing control, for a default (flash-emitting) start. At
        t_d the birth flash, a known free wave, is taken out of the field and, unless shift
        is None, put back advanced by shift. The mouths keep the same birth and lose the same
        energy to the flash; only the flash's later returns move. With shift=None the flash
        is absorbed, as if it escaped to infinity. The energy change at t_d (cross terms with
        the rest of the field) is reported as retime_energy.
        """
        s = s.copy()
        flash0 = self.flash(s) if retime is not None else None
        keys = ('t', 'x', 'L', 'd', 'energy', 'J', 'coulomb', 'exchange', 'own')
        out = {k: [] for k in keys}
        out['energy_parts'] = {k: [] for k in ('rest', 'kinetic', 'coulomb', 'transverse', 'deformation')}
        info = dict(contact_time=None, contact_point=None, released_energy=None,
                    unreleased_angular=None, min_separation=np.pi, retime_time=None, retime_energy=0.)
        maps = dict(t=[], transverse=[], deformation=[])
        if field_points is not None:
            Y = harmonics(self.ell, self.order, field_points)
            uY = harmonics(self.u_ell, self.u_order, field_points)
            next_map = 0.

        def record_map():
            nonlocal next_map
            if field_points is not None and s.t >= next_map - 1e-9:
                tr, de = self.field_density(s, Y, uY)
                maps['t'].append(s.t)
                maps['transverse'].append(tr)
                maps['deformation'].append(de)
                next_map = (np.floor(s.t/field_dt + 1e-9) + 1)*field_dt

        def record():
            e = self.energy_components(s)
            f = self.forces(s)
            for k, v in e.items():
                out['energy_parts'][k].append(v)
            vals = (s.t, s.x.copy(), s.L.copy(), self.separation(s.x), sum(e.values()),
                    self.angular_momentum(s), f['coulomb'], f['exchange'], f['own'])
            for k, v in zip(keys, vals):
                out[k].append(v)

        record()
        record_map()
        next_out = output_dt
        nsteps = int(np.ceil(end/dt - 1e-9))
        h = end/nsteps
        for _ in range(nsteps):
            if s.alive:
                prev = s.copy()
            self.step(s, h)
            if retime is not None and info['retime_time'] is None and s.t >= retime[0] - 1e-9:
                before = sum(self.energy_components(s).values())
                fq, fp = self.free_wave(*flash0, s.t)
                s.Q, s.P = s.Q - fq, s.P - fp
                if retime[1] is not None:
                    fq, fp = self.free_wave(*flash0, s.t + retime[1])
                    s.Q, s.P = s.Q + fq, s.P + fp
                info.update(retime_time=s.t, retime_energy=sum(self.energy_components(s).values()) - before)
                prev = s.copy()
            if s.alive:
                info['min_separation'] = min(info['min_separation'], self.separation(s.x))
            if self.contact(s):
                if refine_contact:
                    s = self._refine(prev, h)
                info.update(contact_time=s.t, contact_point=(s.x[0]+s.x[1]).tolist(),
                            contact_separation=self.separation(s.x), contact_speed=self.velocity(s).tolist())
                record()                                    # the ledger just before and just after surgery
                info['unreleased_angular'] = self.release(s).tolist()
                info['released_energy'] = s.release['energy']
                record()
                # continue on the original grid
                remaining = end - s.t
                if remaining > 1e-12:
                    m = int(np.ceil(remaining/h - 1e-9))
                    for _ in range(m):
                        self.step(s, remaining/m)
                        record_map()
                        if s.t >= next_out - 1e-9:
                            record()
                            next_out = (np.floor(s.t/output_dt+1e-9)+1)*output_dt
                break
            record_map()
            if s.t >= next_out - 1e-9:
                record()
                next_out += output_dt
        arrays = {k: np.asarray(v) for k, v in out.items() if k != 'energy_parts'}
        arrays['energy_parts'] = {k: np.asarray(v) for k, v in out['energy_parts'].items()}
        e = arrays['energy'] - info['retime_energy']*(info['retime_time'] is not None)*(
            arrays['t'] >= (info['retime_time'] or 0.) - 1e-9)
        J = arrays['J']
        if field_points is not None:
            arrays['maps'] = {k: np.asarray(v) for k, v in maps.items()}
        info.update(energy_drift=float(np.abs(e-e[0]).max()/e[0]),
                    angular_drift=float(np.abs(J-J[0]).max()),
                    final_state=s)
        return arrays, info

    def _refine(self, prev, h):
        """Secant iteration on the step length so the step ends at separation 2a."""
        target = 2*self.par.mouth_angle

        def trial(tau):
            s = prev.copy()
            self.step(s, tau)
            return s, self.separation(s.x) - target

        a, fa = 0., self.separation(prev.x) - target
        b = h
        s, fb = trial(b)
        for _ in range(30):
            if abs(fb) < 1e-12 or b == a:
                break
            c = b - fb*(b-a)/(fb-fa)
            c = min(max(c, 1e-6*h), h)
            a, fa = b, fb
            b = c
            s, fb = trial(b)
        return s


def circle_angle(x, P, T):
    """Signed angle along the great circle through P towards T (in (-pi, pi])."""
    x = np.asarray(x)
    return np.arctan2(x @ np.asarray(T), x @ np.asarray(P))


def front_crossings(phi, t, start=None, skip=.05, advance=0.):
    """Times when the birth flash's ray fronts cross each mouth, on the birth circle.

    The birth flash leaves each mouth's start point phi_k(0) as a ring. Its two rays on
    the circle sit at phi_k(0) +/- t (mod 2 pi), meet at the antipode of the start point
    at t = pi and come back. phi has shape (time, 2): the mouths' circle angles.
    Returns dicts with the time, the mouth hit, the source (own or partner), how many
    antipodal foci the ray has passed, and whether it meets the mouth head-on. A flash
    advanced by `advance` (the re-timed control) sits at phi_k(0) +/- (t + advance).
    """
    phi, t = np.asarray(phi), np.asarray(t)
    start = phi[0] if start is None else start
    vel = np.gradient(np.unwrap(phi, axis=0), t, axis=0)
    out = []
    for j in range(2):
        for k in range(2):
            for sign in (1, -1):
                g = np.angle(np.exp(1j*(start[k] + sign*(t + advance) - phi[:, j])))
                change = ((g[:-1] < 0) & (g[1:] >= 0)) | ((g[:-1] > 0) & (g[1:] <= 0))
                hit = change & (np.abs(g[:-1]-g[1:]) < np.pi) & (t[1:] > skip)
                for i in np.flatnonzero(hit):
                    tc = t[i] - g[i]*(t[i+1]-t[i])/(g[i+1]-g[i])
                    out.append(dict(time=float(tc), mouth=j, source='own' if j == k else 'partner',
                                    foci=int(np.floor((tc + advance)/np.pi + 1e-12)), direction=sign,
                                    head_on=bool(sign*vel[i, j] < 0)))
    return sorted(out, key=lambda e: e['time'])


def bare_timeline(width=.25, cutoff=64, operator='conformal', theta=None, t=None):
    """Field of an axisymmetric pulse released at rest at P, with no mouths.

    Returns the energy density along a meridian, e(theta, t) = (u_t^2 + u_theta^2 + m^2 u^2)/2,
    with m^2 = 1/4 for the conformal operator (omega = l + 1/2, exact refocusing) or 0 for
    the massless one (omega = sqrt(l(l+1)), which carries the transverse field above).
    """
    theta = np.linspace(0, np.pi, 361) if theta is None else np.asarray(theta)
    t = np.linspace(0, 4*np.pi, 801) if t is None else np.asarray(t)
    l = np.arange(cutoff+1)
    if operator == 'conformal':
        omega, m2 = l + .5, .25
    elif operator == 'massless':
        omega, m2 = np.sqrt(l*(l+1.)), 0.
    else:
        raise ValueError(operator)
    b = np.exp(-width**2*l*(l+1)/2)*(2*l+1)/(4*np.pi)
    mu = np.cos(theta)
    Pl = np.array([legendre.legval(mu, np.eye(cutoff+1)[k]) for k in range(cutoff+1)])           # (l, theta)
    dPl = np.array([legendre.legval(mu, legendre.legder(np.eye(cutoff+1)[k])) for k in range(cutoff+1)])
    dPl = -dPl*np.sin(theta)[None, :]                                                          # d/dtheta
    co, si = np.cos(np.outer(t, omega)), np.sin(np.outer(t, omega))
    u = (co*b) @ Pl
    ut = (-si*b*omega) @ Pl
    uth = (co*b) @ dPl
    energy = .5*(ut**2 + uth**2 + m2*u**2)
    return theta, t, energy, u, ut
