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
- **Complete history** (`run_history`), as in the archived v39 cycle:
  - The released pulse refocuses at the antipode of the contact point after pi.
  - There `recreate` starts the next pair, conserving energy and canonical angular
    momentum.
  - Each mouth carries its proper time since its own birth. The throat identifies
    equal proper times, so its clock offset is a functional of the whole velocity
    history since that birth (`clock_offsets`), never of where the mouths are now.
  - Offsets are zero at every birth, and clock memory does not survive the surgery.

**Neck channel** (`NeckChannel`, optional). The throat interior becomes a wave
channel, so the clock offset can act back on the motion.
- **Geometry.** #338's catenoid neck from the mouth angle alone: circumference
  2 pi f(s), f = sqrt(s^2 + b^2), b = sin^2 a, s in [-S, S], S = sin a cos a, so
  the length is L = sin 2a and the circumference matches the rims at both ends.
- **Field.** The neck carries the transverse field's dual scalar psi
  (B = psi_t, |E_T| = |grad psi|) in two parts:
  - the axisymmetric m = 0 part (weight 2 pi f);
  - the sideways m = 1 part, cos and sin of the rim angle (weight pi f, with the
    centrifugal term 1/f^2 that makes the narrow waist reflect most of it, as on
    the glued mesh of docs/mouth_kick_experiments.md).
  A mouth moving in a plane through the centre makes B odd under reflection in
  that plane, so the m = 0 part then vanishes and the sideways part carries the
  whole channel.
- **Ports.** Each end is tied by a port spring to the matching rim moment of the
  sphere-side psi around its mouth: the rim average for m = 0, and for m = 1 the
  cos and sin moments in the mouth's rim frame. As the stiffness g grows this
  enforces continuity of psi across the junction.
- **Rim frames and handle.** Each mouth carries a rim frame. Frames are transported
  with the mouth and turned by a spin s_k with moment of inertia M sin^2 a, so the
  port's twist is booked and canonical angular momentum stays exact. At birth e1
  points to the contact point. The neck joins rim angle psi at A to psi + pi at B
  (torus) or to pi - psi (Klein, orientation-reversing); the handle type enters
  only through the m = 1 part.
- **Clocks.** The neck is ultrastatic in its own time T. The throat joins equal
  proper times, so at external time t the neck slice runs from T = tau_A at mouth A
  to T = tau_B at mouth B. In that chart the neck has:
  - a lapse that interpolates the mouths' proper-time rates;
  - a tilt -A/L set by the history offset A = tau_A - tau_B, so signals cross in
    L + A one way and L - A the other.
- **Action.** The constraint dA/dt = 1/gamma_A - 1/gamma_B is imposed with a
  conjugate clock momentum mu, which shifts the rest energies by +mu and -mu. The
  lapse folds the neck and port energies into each mouth's rest energy.
  Altogether H = sum_k sqrt(m_k^2 + |p_k - q_k A_T|^2) + U + H_T, with
  m_k = M +/- mu + (mouth k's lapse share of the neck energy) + (port k)
  + s_k^2 / (2 M sin^2 a).
  It is autonomous: energy and canonical angular momentum stay conserved, and the
  offset's only reaction is mu.

**Scope.** Classical throughout: nothing is quantized. The neck shape and the
rest energy M (neck surface energy plus the flux energy inside the neck) are
held fixed.
"""
from dataclasses import dataclass, field

import numpy as np
from numpy.polynomial import legendre
from scipy import linalg, optimize, sparse, special

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


@dataclass(frozen=True)
class NeckChannel:
    """Throat interior as a wave channel; its geometry follows from the mouth angle."""
    nodes: int = 65              # nodes along the neck (linear elements, lumped masses)
    port: float = 100.           # port stiffness g (m = 0 row; m = 1 rows take g/2); 0 decouples the neck
    clock: str = 'history'       # offset the neck sees: 'history' (A), 'off' (0) or 'reversed' (-A)
    handle: str = 'torus'        # rim map through the neck: 'torus' (psi -> psi + pi) or 'klein' (pi - psi)
    max_tilt: float = .9         # chart stop when |A|/L reaches this


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


def birth_frames(x, point):
    """Rim frames at birth: e1 tangent at each mouth pointing to the contact point, e2 = x x e1."""
    e1 = point - (x @ point)[:, None]*x
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    return np.stack([e1, np.cross(x, e1)], 1)


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
    tau: np.ndarray = field(default_factory=lambda: np.zeros(2))    # mouth proper times since birth
    birth: float = 0.        # external time of the current birth
    cycle: int = 0           # births so far after the first
    eta: np.ndarray = None   # neck field (3, nodes): m = 0, m = 1 cos, m = 1 sin (neck channel only)
    pi: np.ndarray = None    # its momenta
    mu: float = 0.           # clock momentum conjugate to the offset tau_A - tau_B
    frames: np.ndarray = None    # (2, 2, 3) rim frames (e1, e2) of the two mouths
    spin: np.ndarray = None      # (2,) spin angular momenta about the mouth normals

    def copy(self):
        return PairState(self.t, self.x.copy(), self.L.copy(), self.Q.copy(), self.P.copy(), self.alive,
                         None if self.release is None else dict(self.release), self.tau.copy(), self.birth,
                         self.cycle, *(None if v is None else v.copy() for v in (self.eta, self.pi)), self.mu,
                         *(None if v is None else v.copy() for v in (self.frames, self.spin)))


class PairBirth:
    def __init__(self, parameters=Parameters(), neck=None):
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
        self.neck = neck
        self.copies = 2                 # transverse field copies: one per mouth, plus one for the ports
        if neck is not None:
            self._build_neck(neck)

    def _build_neck(self, nk):
        if nk.nodes < 3 or nk.port < 0 or nk.clock not in ('history', 'off', 'reversed') \
                or nk.handle not in ('torus', 'klein') or not 0 < nk.max_tilt < 1:
            raise ValueError('Neck needs >= 3 nodes, port >= 0, known clock and handle modes, 0 < max_tilt < 1')
        a = self.par.mouth_angle
        S, b = np.sin(a)*np.cos(a), np.sin(a)**2
        self.neck_length = L = 2*S
        sigma = np.linspace(0, L, nk.nodes)
        h = L/(nk.nodes - 1)
        mid = (sigma[:-1] + sigma[1:])/2

        def radius(z):
            return np.sqrt((z - S)**2 + b**2)

        def width(z):
            return 2*np.pi*radius(z)

        lumped = np.full(nk.nodes, h)
        lumped[[0, -1]] = h/2
        self.neck_nodes = sigma
        self._n_mass = width(sigma)*lumped                   # int w phi_i (lumped)
        self._n_stiff = width(mid)/h                         # w_e/h per element
        self._n_centre = self._n_mass/radius(sigma)**2       # m = 1 centrifugal term
        self._weight = np.array([1., .5, .5])                # int cos^2 = int sin^2 = pi against 2 pi
        tilt = np.zeros((nk.nodes, nk.nodes))                # int w phi_i phi_j'
        for e in range(nk.nodes - 1):
            tilt[e:e+2, e:e+2] += width(mid[e])*np.array([[-.5, .5], [-.5, .5]])
        self._n_tilt = tilt
        # Lapse shares: the neck time runs at mouth A's proper rate at sigma = 0, B's at sigma = L.
        self._n_share = np.array([1 - sigma/L, sigma/L])
        self._e_share = np.array([1 - mid/L, mid/L])
        # Rim moments of psi = sum (P/omega) Y around a rim at angle a (Funk-Hecke), with the footprint
        # smoothing: the average is P_l(cos a) psi_l(x); the cos and sin moments are kappa_l e.grad psi_l(x).
        rim = np.array([legendre.legval(np.cos(a), np.eye(self.par.cutoff+1)[l]) for l in range(self.par.cutoff+1)])
        self._rim0 = self.smoothing*rim[self.ell]/self.omega
        self._rim1 = self.smoothing*(-2*special.lpmv(1, self.ell, np.cos(a))/self.lam)/self.omega
        self._handle = np.array([[1., 1.], [-1., -1. if nk.handle == 'torus' else 1.]])
        self.inertia = self.par.mass*np.sin(a)**2
        self.clock_sign = dict(history=1., off=0., reversed=-1.)[nk.clock]
        self.copies = 3

    # ---- field pieces -------------------------------------------------------------------------
    def _harmonics(self, x):
        key = np.asarray(x).tobytes()
        if self._cache[0] != key:
            self._cache = (key, harmonics(self.ell, self.order, x).conj())                  # (2, n)
        return self._cache[1]

    def footprints(self, x):
        return self.smoothing*self._harmonics(x)

    def _pad(self, c):
        """Per-mouth copies (2, n) padded with an empty port copy when the neck is present."""
        return np.vstack([c, np.zeros((self.copies - len(c), self.n), complex)])

    def _derivatives(self, c):
        """D_i c for i = 0, 1, 2, shape (3, n)."""
        return (self._stacked @ c).reshape(3, self.n)

    def _rot(self, w, c):
        """Ambient vectors n x grad chi at each mouth for coefficients c, shape (2, 3)."""
        return (w.conj() @ self._derivatives(c).T).real

    def potentials(self, s, w=None):
        """Per source copy k: A[k, j], E_T[k, j] (ambient vectors) and B[k, j] at mouth j."""
        w = self.footprints(s.x) if w is None else w
        A = np.stack([self._rot(w, Q/self.omega) for Q in s.Q])
        E = np.stack([-self._rot(w, P/self.omega) for P in s.P])
        B = np.stack([-(w.conj() @ (self.omega*Q)).real for Q in s.Q])
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

    # ---- neck channel --------------------------------------------------------------------------
    def tilt(self, s):
        """Slope of the neck slice in neck time, dT/dsigma = -A/L as the neck sees the offset."""
        return -self.clock_sign*(s.tau[0] - s.tau[1])/self.neck_length

    def rim_targets(self, s):
        """Sphere-side rim moments of psi at each mouth: the m = 0 averages psi0 with their weights
        w0, the m = 1 moment vectors G (the filtered gradient, ambient) and the filtered field c1."""
        Ystar = self._harmonics(s.x)
        P = s.P.sum(0)
        w0 = self._rim0*Ystar
        c1 = self._rim1*P
        return (w0.conj() @ P).real, w0, np.cross(self._rot(Ystar, c1), s.x), c1

    def port_frames(self, s):
        """Ambient directions of the neck's m = 1 cos and sin rows at each end (handle map at B)."""
        return s.frames*self._handle[:, :, None]

    def neck_shares(self, s):
        """Neck energy split by lapse share: (mouth A's, mouth B's); they sum to the neck energy."""
        k = self.tilt(s)
        u = s.pi - k*(s.eta @ self._n_tilt.T)
        node = u**2/(2*(1 - k**2)*self._n_mass)
        node[1:] += self._n_centre*s.eta[1:]**2/2
        element = self._n_stiff*np.diff(s.eta, axis=1)**2/2
        return self._n_share @ (self._weight @ node) + self._e_share @ (self._weight @ element)

    def port_energies(self, s):
        psi0, _, G, _ = self.rim_targets(s)
        ends = s.eta[:, [0, -1]]
        vec = np.einsum('ck,kcj->kj', ends[1:], self.port_frames(s))
        return .5*self.neck.port*((ends[0] - psi0)**2 + .5*np.sum((vec - G)**2, 1))

    def masses(self, s):
        """Rest energies m_k: M, shifted by the clock momentum (+mu, -mu), by each mouth's share
        of the neck and port energy and by its spin energy."""
        m = np.full(2, self.par.mass)
        if self.neck is None or not s.alive:
            return m
        return (m + np.array([1., -1.])*s.mu + self.neck_shares(s) + self.port_energies(s)
                + s.spin**2/(2*self.inertia))

    def rates(self, s, m=None):
        """Proper-time rates d tau_k/dt = m_k/E_k."""
        m = self.masses(s) if m is None else m
        return 1/self.gamma(s.L, m)

    def _port_reaction(self, s, pull0, pull1, targets):
        """Reactions to the port pulls pull0 (2,) on the m = 0 ends and pull1 (2, 3) on the m = 1
        ends (each the time integral of rate x port force on the neck end): the change of the port
        copy's modes, of the mechanical angular momenta and of the spins. They conserve the
        canonical angular momentum exactly."""
        _, w0, G, c1 = targets
        Ystar = self._harmonics(s.x)
        arms = np.cross(s.x, pull1)
        dQ = -pull0 @ w0
        dL = pull0[:, None]*self._rot(w0, s.P.sum(0))
        for k in range(2):
            dQ = dQ + self._rim1*(arms[k] @ self._derivatives(Ystar[k]))
            dL[k] += self._rot(Ystar[k:k+1], arms[k] @ self._derivatives(c1))[0]
        return dQ, dL, np.einsum('kj,kj->k', s.x, np.cross(pull1, G))

    # ---- mechanics -----------------------------------------------------------------------------
    def gamma(self, L, m=None):
        m = self.par.mass if m is None else m
        return np.sqrt(1+np.sum(L*L, -1)/m**2)

    def velocity(self, s, m=None):
        m = self.masses(s) if m is None else m
        return np.cross(s.L, s.x)/(m*self.gamma(s.L, m))[:, None]

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
        x, vel, point, _ = birth_positions(self.par.mouth_angle, v, plane)
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
        s = PairState(0., x, np.cross(x, p), self._pad(Q), np.zeros((self.copies, self.n), complex))
        if self.neck is not None:
            self._empty_neck(s, point)
        return s

    def _empty_neck(self, s, point):
        """A newborn neck: clock momentum and spins zero, rim frames from the contact point, and
        the neck field at rest in the static continuation of the sphere field, the minimum of
        neck plus port energy for the current rim moments (zero for a field-free birth)."""
        n, g = self.neck.nodes, self.neck.port
        s.mu, s.spin, s.frames = 0., np.zeros(2), birth_frames(s.x, point)
        s.eta, s.pi = np.zeros((3, n)), np.zeros((3, n))
        if g > 0:
            psi0, _, G, _ = self.rim_targets(s)
            side = np.einsum('kj,kcj->ck', G, self.port_frames(s))
            for r, (a, b) in enumerate((psi0, *side)):
                rhs = np.zeros(n)
                rhs[[0, -1]] = g*a, g*b
                s.eta[r] = np.linalg.solve(self._neck_stiffness(np.ones(2), r > 0), rhs)

    def quasi_static(self, x, vel):
        """Magnetostatic transverse field of the currents q_k v_k at x_k, one copy per mouth."""
        return self.source(x, vel)/self.lam

    def free_wave(self, Q, P, tau):
        """Free transverse evolution of (Q, P) by tau (any sign)."""
        co, si = np.cos(self.omega*tau), np.sin(self.omega*tau)
        return Q*co + P*si/self.omega, P*co - self.omega*Q*si

    def flash(self, s0):
        """The birth flash of a default start: the free wave launched by -Q_qs at t = 0."""
        return -self._pad(self.quasi_static(s0.x, self.velocity(s0))), np.zeros((self.copies, self.n), complex)

    def free_wave_force(self, history, Q0, P0, t0=0.):
        """Force of a free transverse wave (Q0, P0 at time t0, one copy per source) on each
        mouth along a recorded trajectory: shape (time, 2 mouths, 2 sources, 3)."""
        out = np.zeros((len(history['t']), 2, 2, 3))
        mass = np.asarray(history['masses']) if 'masses' in history else self.par.mass
        g = np.sqrt(1 + np.sum(history['L']**2, -1)/mass**2)
        for i, (t, x, L) in enumerate(zip(history['t'], history['x'], history['L'])):
            if history['energy_parts']['rest'][i] == 0:
                continue
            Q, P = self.free_wave(Q0, P0, t - t0)
            w = self.footprints(x)
            v = np.cross(L, x)/(np.broadcast_to(mass, g.shape)[i]*g[i])[:, None]
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
        """Energy ledger. With the neck, rest + kinetic + neck + port + spin = sum_k E_k: the
        clock momentum shifts the two rest energies oppositely and so drops out of the sum."""
        out = dict(rest=0., kinetic=0., coulomb=0., transverse=self.transverse_energy(s),
                   deformation=self.deformation_energy(s))
        if self.neck is not None:
            out.update(neck=0., port=0., spin=0.)
        if s.alive:
            m = self.masses(s)
            g = self.gamma(s.L, m)
            out.update(rest=2*self.par.mass, kinetic=float(np.sum(m*(g-1))),
                       coulomb=float(self.coulomb(s.x)[0]))
            if self.neck is not None:
                out.update(neck=float(self.neck_shares(s).sum()), port=float(self.port_energies(s).sum()),
                           spin=float(np.sum(s.spin**2)/(2*self.inertia)))
        return out

    def field_angular(self, Q, P):
        return -np.array([np.vdot(P, g @ Q).real for g in self.generators])

    def angular_momentum(self, s):
        """Total canonical angular momentum: mouths x cross (p + qA), their spins and the field."""
        J = self.field_angular(s.Q.sum(0), s.P.sum(0)) + self.deformation_angular(s)
        if s.alive:
            A, _, _ = self.potentials(s)
            p = np.cross(s.L, s.x) + self.charges[:, None]*A.sum(0)
            J = J + np.cross(s.x, p).sum(0)
            if self.neck is not None:
                J = J + s.spin @ s.x
        return J

    def forces(self, s):
        """Instantaneous forces on each mouth, by source: Coulomb from the partner,
        transverse from the partner's copy, transverse from the mouth's own copy and, with
        the neck, transverse from the field the ports emitted (throat) and the port spring's
        pull (port, per unit proper time)."""
        keys = ('coulomb', 'exchange', 'own') + (('throat', 'port') if self.neck is not None else ())
        out = {k: np.zeros((2, 3)) for k in keys}
        if not s.alive:
            return out
        out['coulomb'] = self.coulomb(s.x)[1]
        if self.par.transverse:
            _, E, B = self.potentials(s)
            v = self.velocity(s)
            for j in range(2):
                for k in range(self.copies):
                    f = self.charges[j]*(E[k, j] + B[k, j]*np.cross(v[j], s.x[j]))
                    out['throat' if k == 2 else 'own' if k == j else 'exchange'][j] += f
        if self.neck is not None:
            targets = self.rim_targets(s)
            ends = s.eta[:, [0, -1]]
            pull0 = self.neck.port*(ends[0] - targets[0])
            pull1 = .5*self.neck.port*(np.einsum('ck,kcj->kj', ends[1:], self.port_frames(s)) - targets[2])
            _, dL, _ = self._port_reaction(s, pull0, pull1, targets)
            out['port'] = np.cross(dL, s.x)
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

    def _magnetic(self, s, tau, m):
        """Rotation of each mouth's momentum about its normal by the magnetic field, and by the
        gyroscopic force s (x cross v) of a spinning mouth: the sphere's curvature acts on the spin
        like a monopole, which keeps spin + orbital angular momentum conserved."""
        strength = np.zeros(2)
        if self.par.transverse:
            strength += self.charges*self.potentials(s)[2].sum(0)
        if s.spin is not None:
            strength -= s.spin
        if np.any(strength):
            s.L = rotate(s.L, s.x, -strength*tau/(m*self.gamma(s.L, m)))

    def _coulomb_kick(self, s, tau):
        s.L += tau*np.cross(s.x, self.coulomb(s.x)[1])

    def _particle_flow(self, s, dt):
        """Mouths in the frozen transverse field and Coulomb force; the field modes take
        the mouths' currents as a source (trapezoid along the drift). The rest energies
        are frozen over the flow; the proper times, and with them the offset, advance."""
        m = self.masses(s)
        self._coulomb_kick(s, dt/2)
        self._magnetic(s, dt/2, m)
        if self.par.transverse:
            s.P[:2] += dt/2*self.source(s.x, self.velocity(s, m))
        length = np.linalg.norm(s.L, axis=1)
        g = self.gamma(s.L, m)
        s.tau += dt/g                                   # proper time, midpoint rule over the step
        for j in range(2):
            if length[j] > 0:
                s.x[j] = rotate(s.x[j], s.L[j]/length[j], length[j]*dt/(m[j]*g[j]))
                if s.frames is not None:                # rim frames are carried along (parallel transport)
                    s.frames[j] = rotate(s.frames[j], s.L[j]/length[j], length[j]*dt/(m[j]*g[j]))
        s.x /= np.linalg.norm(s.x, axis=1, keepdims=True)
        if s.frames is not None:                        # spin turns them about the normal
            e1 = rotate(s.frames[:, 0], s.x, s.spin*dt/(self.inertia*g))
            e1 -= np.sum(e1*s.x, 1, keepdims=True)*s.x
            e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
            s.frames = np.stack([e1, np.cross(s.x, e1)], 1)
        if self.par.transverse:
            s.P[:2] += dt/2*self.source(s.x, self.velocity(s, m))
        self._magnetic(s, dt/2, m)
        self._coulomb_kick(s, dt/2)

    def _neck_stiffness(self, rates, sideways):
        """Lapse-weighted stiffness of one neck row, with the port springs on its end nodes."""
        g = self.neck.port
        kd = rates @ self._e_share*self._n_stiff
        K = np.diag(np.r_[kd, 0] + np.r_[0, kd]) - np.diag(kd, 1) - np.diag(kd, -1)
        if sideways:
            K += np.diag((rates @ self._n_share)*self._n_centre)
        K[0, 0] += rates[0]*g
        K[-1, -1] += rates[1]*g
        return K

    def neck_generator(self, k, rates, sideways=False, targets=((0., 0.),)):
        """Linear generator of the neck and port flow at fixed tilt k and proper-time rates:
        d/dt (eta, pi, e) = Z (eta, pi, e) for one row (m = 0, or an m = 1 row with sideways=True),
        with the port targets held fixed. e holds one unit source per target pair (end A, end B)."""
        n, g = self.neck.nodes, self.neck.port
        D = (rates @ self._n_share)/((1 - k**2)*self._n_mass)
        K = self._neck_stiffness(rates, sideways)
        B = self._n_tilt
        Z = np.zeros((2*n + len(targets),)*2)
        Z[:n, :n] = -k*D[:, None]*B
        Z[:n, n:2*n] = np.diag(D)
        Z[n:2*n, :n] = -(K + k**2*B.T @ (D[:, None]*B))
        Z[n:2*n, n:2*n] = k*B.T*D[None, :]
        for j, (a, b) in enumerate(targets):
            Z[n, 2*n+j] = rates[0]*g*a
            Z[2*n-1, 2*n+j] = rates[1]*g*b
        return Z

    def _tilt_force(self, k, rates, eta, pi):
        """d/dk of the lapse-weighted neck Hamiltonian at fixed (eta, pi), all rows."""
        D = (rates @ self._n_share)/((1 - k**2)*self._n_mass)
        Be = eta @ self._n_tilt.T
        u = pi - k*Be
        return float(self._weight @ np.sum(D*u*(k*u/(1 - k**2) - Be), 1))

    def _neck_flow(self, s, h):
        """Neck interior and ports for h, with the mouths' rates, the offset, the rim frames and
        the sphere field's rim moments frozen. The neck evolves exactly. Simpson's rule gives the
        integrated reactions: the ports pull on the field (port copy), on the mouths and on their
        spins, and the offset's tilt force moves the clock momentum, d mu/dt = -dH/dA."""
        k = self.tilt(s)
        if abs(k) >= self.neck.max_tilt:
            raise RuntimeError(f'Neck chart stop: |A|/L = {abs(k):.3f} at t = {s.t:.4f}')
        n, g = self.neck.nodes, self.neck.port
        rates = self.rates(s)
        targets = self.rim_targets(s)
        psi0, _, G, _ = targets
        F = self.port_frames(s)
        side = np.einsum('kj,kcj->ck', G, F)                          # m = 1 targets (row, end)
        flows = (linalg.expm(self.neck_generator(k, rates, False, [psi0])*h/2),
                 linalg.expm(self.neck_generator(k, rates, True, side)*h/2))
        rows = []
        for r, (Eh, e) in enumerate(((flows[0], [1.]), (flows[1], [1., 0.]), (flows[1], [0., 1.]))):
            z0 = np.r_[s.eta[r], s.pi[r], e]
            z1 = Eh @ z0
            rows.append((z0, z1, Eh @ z1))
        simpson = (1., 4., 1.)
        ends = np.array([sum(c*z[[0, n-1]] for c, z in zip(simpson, zs)) for zs in rows])*h/6     # (row, end)
        pull0 = rates*g*(ends[0] - h*psi0)
        pull1 = .5*rates[:, None]*g*(np.einsum('ck,kcj->kj', ends[1:], F) - h*G)
        dQ, dL, ds = self._port_reaction(s, pull0, pull1, targets)
        # Canonical momenta are fixed under this flow: the mechanical ones take the ports' pull
        # and, as in the field flow, -q dA from the port copy's change.
        s.Q[2] += dQ
        s.L += dL - self.charges[:, None]*np.cross(s.x, self._rot(self.footprints(s.x), dQ/self.omega))
        s.spin += ds
        eta = [np.array([zs[i][:n] for zs in rows]) for i in range(3)]
        pi = [np.array([zs[i][n:2*n] for zs in rows]) for i in range(3)]
        dk = sum(c*self._tilt_force(k, rates, eta[i], pi[i]) for i, c in enumerate(simpson))*h/6
        s.mu += self.clock_sign/self.neck_length*dk
        s.eta, s.pi = eta[2], pi[2]

    def step(self, s, dt):
        if s.alive:
            if self.neck is not None:
                self._neck_flow(s, dt/2)
            self._field_flow(s, dt/2)
            self._particle_flow(s, dt)
            self._field_flow(s, dt/2)
            if self.neck is not None:
                self._neck_flow(s, dt/2)
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
          which must vanish with the charges, and with the neck its field and port energy.
          The transverse field, including what the ports emitted, continues freely.
        - Pulse: u_t = alpha s + sum_i gamma_i D_i s and u = sum_i beta_i D_i s, with
          s = exp(-w^2 l(l+1)/2) Y*(c). The beta and gamma make the pulse carry the mouths'
          canonical angular momentum (gamma is needed for the spin about c, which an
          axisymmetric shape alone cannot hold), and alpha then fixes the energy. They are
          solved by least squares with the smallest beta and gamma.
        - Returns the angular momentum left unaccounted (zero to solver tolerance).
        """
        e = self.energy_components(s)
        released = sum(v for k, v in e.items() if k not in ('transverse', 'deformation'))
        J_before = self.angular_momentum(s)
        c = s.x[0] + s.x[1]
        c /= np.linalg.norm(c)
        s.alive = False
        if self.neck is not None:
            s.eta, s.pi, s.mu, s.spin = np.zeros_like(s.eta), np.zeros_like(s.pi), 0., np.zeros(2)
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
        s.release = dict(time=s.t, point=c, energy=released, p0=p0, a0=a0,
                         x=s.x.copy(), p=np.cross(s.L, s.x), tau=s.tau.copy())
        return J_before - self.angular_momentum(s)

    def recreate(self, s):
        """Inverse surgery: re-create the pair where the released pulse refocuses (v39 RECONNECT).

        At t_c + pi the conformal deformation pulse released at contact point c has refocused
        at -c. The rule removes the whole deformation field and puts a touching pair there:
        - Positions are the antipodal images of the contact positions, -x_k: still touching
          (separation 2a), midpoint -c.
        - Momenta start from the contact momenta. As ambient vectors these are tangent at -x_k
          and point apart. They are then corrected by least squares in the tangent planes, so
          the pair carries exactly the pulse's energy (rest + kinetic + Coulomb at 2a, plus
          the port energy of an empty neck) and canonical angular momentum. For a symmetric
          pair the uncorrected mirror already does.
        - Clocks restart: the new mouths are born together, so their offset is zero. Clock
          memory does not survive the topology change: the neck starts empty and mu at zero.
        The transverse field is untouched. Like the release, this is a stated rule, not
        derived dynamics. Returns the energy and angular-momentum residuals.
        """
        if s.alive or s.release is None:
            raise ValueError('Re-creation needs a released pair')
        E_target = self.deformation_energy(s)
        J_target = self.deformation_angular(s)
        x = -s.release['x']
        p0 = s.release['p']
        basis = []
        for k in range(2):
            e1 = np.cross(x[k], p0[k]) if np.linalg.norm(p0[k]) > 1e-12 else np.cross(x[k], [1., 0, 0])
            e1 /= np.linalg.norm(e1)
            basis.append((e1, np.cross(x[k], e1)))
        w = self.footprints(x)
        A = self._rot(w, s.Q.sum(0)/self.omega) if self.par.transverse else np.zeros((2, 3))
        U = self.coulomb(x)[0]
        born = s.copy()
        born.alive, born.x, born.tau = True, x.copy(), np.zeros(2)
        if self.neck is not None:
            self._empty_neck(born, -s.release['point'])
        m = self.masses(born)

        def momenta(b):
            return np.array([p0[k] + b[2*k]*basis[k][0] + b[2*k+1]*basis[k][1] for k in range(2)])

        def residual(b):
            p = momenta(b)
            J = np.cross(x, p + self.charges[:, None]*A).sum(0)
            E = np.sqrt(m**2 + np.sum(p*p, 1)).sum() + U
            return np.r_[J - J_target, (E - E_target)/E_target, 1e-8*b]

        b = optimize.least_squares(residual, np.zeros(4), xtol=1e-15, ftol=1e-15, gtol=1e-15).x
        if abs(residual(b)[3]) > 1e-9:
            raise RuntimeError(f'Re-creation cannot carry the pulse energy (relative residual {residual(b)[3]:.1e})')
        p = momenta(b)
        s.release = None
        s.alive = True
        s.x = x.copy()
        s.L = np.cross(x, p)
        s.tau = np.zeros(2)
        s.birth = s.t
        s.cycle += 1
        s.eta, s.pi, s.mu, s.spin, s.frames = born.eta, born.pi, born.mu, born.spin, born.frames
        if self.separation_rate(s) <= 0:
            raise RuntimeError('Re-created pair is not separating')
        return dict(energy=float(sum(self.energy_components(s).values()) - E_target - self.transverse_energy(s)),
                    angular=(np.cross(x, p + self.charges[:, None]*A).sum(0) - J_target).tolist(),
                    correction=float(np.linalg.norm(b)))

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
        out['energy_parts'] = {k: [] for k in self.energy_components(s)}
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


    def run_history(self, s, end, cycles=1, dt=.004, output_dt=.02, neck_length=.565, observer=None):
        """Complete history: birth, flight, contact, release, refocusing, re-birth, repeated.

        Each cycle runs from a birth with synchronised mouth clocks to rim contact. The
        released pulse refocuses at the antipode of the contact point after pi, where
        `recreate` starts the next cycle, up to `cycles` births in all. Besides the ledgers,
        each mouth's proper time since its birth is recorded.

        `clock_offsets` turns the proper times into the throat's external-time offset. The
        throat identifies equal proper times on the two mouths, so this offset depends on
        the whole velocity history since birth, not on where the mouths are.
        neck_length (#338's 2S) enters only the chronology margin d + L_neck - |offset|; with
        the neck channel its own length sin 2a is used.

        With the neck channel the offset also acts back. The record then adds the clock
        momentum mu, the spins, the neck end values (row, end) and the forces from the field the
        ports emitted ('throat'), from the port springs ('port') and, for scale, from the
        partner's field ('exchange'). The rest energies m_k are always recorded.

        observer(state), if given, sees the full state at every recorded time (for field maps
        and animations); it must not change it.
        """
        s = s.copy()
        neck_length = neck_length if self.neck is None else self.neck_length
        keys = ('t', 'x', 'L', 'd', 'tau', 'alive', 'cycle', 'energy', 'J', 'masses')
        if self.neck is not None:
            keys += ('mu', 'spin', 'ends', 'throat', 'port', 'exchange')
        out = {k: [] for k in keys}
        out['energy_parts'] = {k: [] for k in self.energy_components(s)}
        events = []

        def record():
            e = self.energy_components(s)
            for k, v in e.items():
                out['energy_parts'][k].append(v)
            vals = (s.t, s.x.copy(), s.L.copy(), self.separation(s.x), s.tau.copy(), s.alive, s.cycle,
                    sum(e.values()), self.angular_momentum(s), self.masses(s))
            if self.neck is not None:
                f = self.forces(s)
                vals += (s.mu, s.spin.copy(), s.eta[:, [0, -1]].copy(), f['throat'], f['port'], f['exchange'])
            for k, v in zip(keys, vals):
                out[k].append(v)
            if observer is not None:
                observer(s)

        record()
        next_out = output_dt
        rebirth = None
        while s.t < end - 1e-12:
            if s.alive:
                prev = s.copy()
                self.step(s, min(dt, end - s.t))
                if self.contact(s):
                    s = self._refine(prev, s.t - prev.t)
                    record()
                    residual = self.release(s)
                    events.append(dict(kind='contact', time=s.t, cycle=s.cycle, point=s.release['point'].tolist(),
                                       clocks=s.release['tau'].tolist(), energy=s.release['energy'],
                                       speed=np.linalg.norm(s.release['p'], axis=1).tolist(),
                                       unreleased_angular=residual.tolist()))
                    record()
                    rebirth = s.t + np.pi if s.cycle + 1 < cycles else None
                    next_out = (np.floor(s.t/output_dt + 1e-9) + 1)*output_dt
                    continue
            else:
                target = min(end, rebirth) if rebirth is not None else end
                n = max(1, int(np.ceil((target - s.t)/dt - 1e-9)))
                h = (target - s.t)/n
                for _ in range(n):
                    self.step(s, h)
                    if s.t >= next_out - 1e-9:
                        record()
                        next_out = (np.floor(s.t/output_dt + 1e-9) + 1)*output_dt
                if rebirth is not None and abs(s.t - rebirth) < 1e-9:
                    record()
                    check = self.recreate(s)
                    events.append(dict(kind='birth', time=s.t, cycle=s.cycle,
                                       point=((s.x[0]+s.x[1])/np.linalg.norm(s.x[0]+s.x[1])).tolist(),
                                       speed=np.linalg.norm(self.velocity(s), axis=1).tolist(), **check))
                    record()
                    rebirth = None
                    next_out = (np.floor(s.t/output_dt + 1e-9) + 1)*output_dt
                continue
            if s.t >= next_out - 1e-9:
                record()
                next_out += output_dt
        arrays = {k: np.asarray(v) for k, v in out.items() if k != 'energy_parts'}
        arrays['energy_parts'] = {k: np.asarray(v) for k, v in out['energy_parts'].items()}
        e, J = arrays['energy'], arrays['J']
        offsets = clock_offsets(arrays)
        margin = np.where(arrays['alive'], arrays['d'] + neck_length - np.abs(offsets['external']), np.nan)
        arrays.update(offsets, chronology_margin=margin)
        info = dict(events=events, energy_drift=float(np.abs(e-e[0]).max()/e[0]),
                    angular_drift=float(np.abs(J-J[0]).max()), final_state=s, neck_length=neck_length,
                    min_chronology_margin=float(np.nanmin(margin)),
                    max_clock_difference=float(np.nanmax(np.abs(offsets['clock']))),
                    max_external_offset=float(np.nanmax(np.abs(offsets['external']))))
        return arrays, info


def clock_offsets(history, mass=1.):
    """Throat clock offsets from recorded proper times (rest energies from history['masses']
    when recorded, else `mass`).

    - clock: tau_A - tau_B at equal external time, i.e. what one sees comparing the two
      mouth clocks through the throat.
    - external: the external-time shift t_B - t of the event on B that the throat
      identifies with A at time t, the one with tau_B(t_B) = tau_A(t). It is found
      within the same cycle by inverting tau_B(t), and is NaN where that event falls
      outside the recorded flight (after contact).
    - first_order: gamma_B (tau_A - tau_B), the small-offset approximation.
    Both offsets are zero at every birth and depend only on the velocity histories.
    """
    t, tau, alive, cycle = (np.asarray(history[k]) for k in ('t', 'tau', 'alive', 'cycle'))
    clock = np.where(alive, tau[:, 0] - tau[:, 1], np.nan)
    external = np.full(len(t), np.nan)
    for c in np.unique(cycle):
        sel = np.flatnonzero(alive & (cycle == c))
        if len(sel) < 2:
            continue
        tb, taub = t[sel], tau[sel, 1]
        keep = np.r_[True, np.diff(taub) > 0]
        tb, taub = tb[keep], taub[keep]
        inside = (tau[sel, 0] >= taub[0]) & (tau[sel, 0] <= taub[-1])
        external[sel[inside]] = np.interp(tau[sel, 0][inside], taub, tb) - t[sel][inside]
    L = np.asarray(history['L'])
    m_b = np.asarray(history['masses'])[:, 1] if 'masses' in history else mass
    gamma_b = np.sqrt(1 + np.sum(L[:, 1]**2, -1)/m_b**2)
    return dict(clock=clock, external=external, first_order=gamma_b*clock)


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
