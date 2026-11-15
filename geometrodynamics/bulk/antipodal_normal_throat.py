"""A real angular harmonic on a responding, normal-reversing throat.

Public freeze: 533f925. Two normal samples and a periodic phi lattice;
theta is a continuous coordinate retained as one exact Fourier harmonic.
F = a(phi,u) cos(n theta) + b(phi,u) sin(n theta).
The normalized basis makes L_theta = n sum(a p_b - b p_a).
This is a Hamiltonian toy model, not a GR or Maxwell throat solution.
"""

from dataclasses import dataclass

import numpy as np

PREREG = "533f9257fcb28f015740b9f08557959bc3c1a83b"
SEED = 2026111201
AMPLITUDES = (0.05, 0.1, 0.2, 0.4)
CELLS = {"shift_kept": (1, 1), "antipodal": (1, -1),
         "reflective": (-1, 1), "reflect_reversed": (-1, -1)}


def normal_parity(m, n, eta=1):
    """Allowed normal parity for a scalar or sign-twisted scalar."""
    if eta not in (-1, 1):
        raise ValueError("eta must be +1 or -1")
    return eta * (-1) ** (m + n)


def jacobian(kind):
    c, r = CELLS[kind]
    return np.diag([1., float(c), float(r)])


def identify(points, kind="antipodal"):
    """The actual coordinate identification, angles modulo 2 pi."""
    out = np.asarray(points, dtype=float).copy()
    out[..., 0] += np.pi
    c, r = CELLS[kind]
    if c == 1:
        out[..., 1] += np.pi
    else:
        out[..., 1] *= -1
    out[..., 2] *= r
    out[..., :2] %= 2*np.pi
    return out


@dataclass(frozen=True)
class Throat:
    """State [x, p, q, pq], x and p shaped (quadrature, face, phi).

    q is an outward-opening scalar. Its signed normal displacement is
    (q_plus, -q_minus), so it swaps without a sign under the antipodal deck.
    All public state methods accept one state or a batch of row states.
    """

    sites: int = 16
    n: int = 1
    kappa: float = 0.75
    mass_sq: float = 0.25
    barrier: float = 2.0
    mouth_mass: float = 1.0
    support: float = 4.0

    def __post_init__(self):
        if self.sites < 8 or self.sites % 2 or self.n < 1:
            raise ValueError("even sites >= 8 and positive n required")
        if min(self.kappa, self.mass_sq, self.barrier) < 0:
            raise ValueError("negative stiffness")
        if self.mouth_mass <= 0 or self.support <= 0:
            raise ValueError("positive mouth mass and support required")

    @property
    def field_size(self):
        return 4 * self.sites

    @property
    def dim(self):
        return 2 * self.field_size + 4

    @property
    def phi(self):
        return 2 * np.pi * np.arange(self.sites) / self.sites

    @property
    def profile(self):
        return 0.6 + 0.4 * np.cos(2 * self.phi)

    def unpack(self, z):
        z = np.atleast_2d(np.asarray(z, dtype=float))
        b, f = len(z), self.field_size
        return (z[:, :f].reshape(b, 2, 2, self.sites),
                z[:, f:2*f].reshape(b, 2, 2, self.sites),
                z[:, 2*f:2*f+2], z[:, 2*f+2:])

    def pack(self, x, p, q, pq):
        b = len(x)
        return np.concatenate((x.reshape(b, -1), p.reshape(b, -1), q, pq), axis=1)

    def potential(self, q):
        return self.barrier * np.exp(-2 * q[:, :, None] * self.profile)

    def apply_operator(self, x, q):
        dx = 2 * np.pi / self.sites
        lap = (2*x - np.roll(x, 1, axis=-1) - np.roll(x, -1, axis=-1)) / dx**2
        return (lap + (self.n**2 + self.mass_sq) * x
                + self.kappa * (x - x[:, :, ::-1])
                + self.potential(q)[:, None] * x)

    def loads(self, z):
        """Field contributions to outward-opening generalized mouth forces."""
        x, _, q, _ = self.unpack(z)
        return np.sum(self.profile * self.potential(q) * np.sum(x*x, axis=1), axis=-1)

    def rhs(self, z, respond=True):
        x, p, q, pq = self.unpack(z)
        dq = pq / self.mouth_mass if respond else np.zeros_like(q)
        dpq = self.loads(z) - self.support*q if respond else np.zeros_like(q)
        return self.pack(p, -self.apply_operator(x, q), dq, dpq)

    def energy(self, z):
        x, p, q, pq = self.unpack(z)
        return (0.5 * np.sum(p*p + x*self.apply_operator(x, q), axis=(1, 2, 3))
                + np.sum(pq*pq, axis=1)/(2*self.mouth_mass)
                + 0.5*self.support*np.sum(q*q, axis=1))

    def angular_momentum(self, z, per_face=False):
        x, p, _, _ = self.unpack(z)
        values = self.n * np.sum(x[:, 0]*p[:, 1] - x[:, 1]*p[:, 0], axis=-1)
        return values if per_face else values.sum(axis=1)

    def normal_angular_flux(self, z, local=False):
        """Into the + face: d L_plus / dt after summing over phi.

        The minus face receives the negative. This diagnoses transfer
        across the normal link, not circulation around theta.
        """
        x, _, _, _ = self.unpack(z)
        current = self.kappa*self.n*(x[:, 0, 0]*x[:, 1, 1]
                                     - x[:, 1, 0]*x[:, 0, 1])
        return current if local else current.sum(axis=-1)

    def deck(self, z, kind="antipodal", eta=1):
        if eta not in (-1, 1):
            raise ValueError("eta must be +1 or -1")
        x, p, q, pq = self.unpack(z)
        xx = np.roll(x, self.sites//2, axis=-1)
        pp = np.roll(p, self.sites//2, axis=-1)
        c, r = CELLS[kind]
        if c == 1:
            xx, pp = (-1)**self.n * xx, (-1)**self.n * pp
        else:
            sign = np.array([1., -1.])[None, :, None, None]
            xx, pp = xx*sign, pp*sign
        if r == -1:
            xx, pp = xx[:, :, ::-1], pp[:, :, ::-1]
            q, pq = q[:, ::-1], pq[:, ::-1]
        return self.pack(eta*xx, eta*pp, q, pq)

    def project(self, z, kind="antipodal", eta=1):
        return (np.atleast_2d(z) + self.deck(z, kind, eta))/2

    def symplectic(self):
        f = self.field_size
        w = np.zeros((self.dim, self.dim))
        w[:f, f:2*f], w[f:2*f, :f] = np.eye(f), -np.eye(f)
        w[-4:-2, -2:], w[-2:, -4:-2] = np.eye(2), -np.eye(2)
        return w

    def frequency_sq(self, m, epsilon):
        dx = 2*np.pi/self.sites
        return (4*np.sin(np.pi*m/self.sites)**2/dx**2 + self.n**2
                + self.mass_sq + self.barrier + self.kappa*(1-epsilon))

    def mode(self, m, epsilon, sine=False):
        if not 0 <= m <= self.sites//2 or epsilon not in (-1, 1):
            raise ValueError("unresolved mode or invalid normal parity")
        wave = np.sin(m*self.phi) if sine else np.cos(m*self.phi)
        if np.linalg.norm(wave) < 1e-12:
            raise ValueError("zero sine mode")
        g = np.array([1., epsilon])[:, None] * wave
        return g/np.linalg.norm(g)

    def travelling(self, amplitude, m=0, eta=1):
        epsilon = normal_parity(m, self.n, eta)
        x = np.zeros((1, 2, 2, self.sites))
        p = np.zeros_like(x)
        x[0, 0] = amplitude*self.mode(m, epsilon)
        p[0, 1] = np.sqrt(self.frequency_sq(m, epsilon))*x[0, 0]
        return self.pack(x, p, np.zeros((1, 2)), np.zeros((1, 2)))[0]

    def frozen_matrix(self):
        """Independently assembled full-grid operator, one quadrature."""
        n = self.sites
        lap = np.zeros((n, n))
        for j in range(n):
            lap[j, j] += 2
            lap[j, (j-1) % n] -= 1
            lap[j, (j+1) % n] -= 1
        lap /= (2*np.pi/n)**2
        base = lap + (self.n**2+self.mass_sq+self.barrier)*np.eye(n)
        return np.kron(np.eye(2), base) + self.kappa*np.kron(
            np.array([[1., -1.], [-1., 1.]]), np.eye(n))

    def susceptibility(self, m, omega):
        """Frozen two-face mode response; away from its poles."""
        lam = self.frequency_sq(m, 1)
        return np.linalg.inv(np.array([[lam+self.kappa-omega**2, -self.kappa],
                                       [-self.kappa, lam+self.kappa-omega**2]]))


def step(model, z, h, respond=True):
    """Velocity Verlet; each substep conserves total angular momentum."""
    z = np.atleast_2d(z).copy()
    x, p, q, pq = model.unpack(z)
    p -= 0.5*h*model.apply_operator(x, q)
    if respond:
        pq += 0.5*h*(model.loads(z)-model.support*q)
    x += h*p
    if respond:
        q += h*pq/model.mouth_mass
    p -= 0.5*h*model.apply_operator(x, q)
    if respond:
        pq += 0.5*h*(model.loads(z)-model.support*q)
    return z


def flow(model, start, duration=40., h=0.01, sample=10, respond=True):
    steps = int(round(duration/h))
    if sample < 1 or steps % sample or not np.isclose(steps*h, duration):
        raise ValueError("duration must contain an integral number of sample blocks")
    z = np.atleast_2d(start).copy()
    history = [z.copy()]
    for j in range(steps):
        z = step(model, z, h, respond)
        if (j+1) % sample == 0:
            history.append(z.copy())
    return np.stack(history)


def classify(checks, required):
    missing = sorted(set(required)-checks.keys())
    failed = sorted(k for k in required if k in checks and not checks[k])
    return {"numerical_validity": "INCOMPLETE" if missing else
            "UNRESOLVED_NUMERICALLY" if failed else "PASS",
            "failed_checks": failed, "missing_checks": missing,
            "electric_charge": "NOT_DERIVED", "spin_lift": "NOT_DERIVED",
            "time_extension": "NOT_DERIVED", "einstein_dynamics": "NOT_DERIVED",
            "support": "ASSUMED_OSCILLATOR", "orienting_wave": "INITIAL_DATA"}
