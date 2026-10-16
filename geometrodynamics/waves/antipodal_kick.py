"""Conservative S² scalar/attachment comparator; see docs/antipodal_kick_protocol.md.

This operator is -Delta+1/4 on the round unit sphere, NOT the handled
minimal scalar, a soap-film constitutive law, or a GR perturbation equation.
"""
from dataclasses import dataclass

import numpy as np
from scipy import sparse, special


def harmonics(ell, order, points):
    """Y_lm at unit Cartesian points, with scipy 1.10+ compatibility."""
    points = np.asarray(points)
    theta = np.arccos(np.clip(points[..., 2], -1, 1))
    phi = np.arctan2(points[..., 1], points[..., 0])
    if hasattr(special, "sph_harm_y"):
        return special.sph_harm_y(ell, order, theta[..., None], phi[..., None])
    return special.sph_harm(order, ell, phi[..., None], theta[..., None])


@dataclass(frozen=True)
class Parameters:
    cutoff: int = 32
    sigma: float = 0.20
    footprint: float = 0.12
    mass: float = 0.30
    shape_mass: float = 0.05
    shape_frequency: float = 5.0
    coupling: float = 1.0
    energy: float = 0.1


@dataclass
class State:
    a: np.ndarray
    p: np.ndarray
    x: np.ndarray
    angular: np.ndarray
    q: float = 0.0
    pq: float = 0.0

    def copy(self):
        return State(self.a.copy(), self.p.copy(), self.x.copy(),
                     self.angular.copy(), self.q, self.pq)


class AntipodalKick:
    """Hamiltonian spectral evolution with exact rotational momentum balance."""

    def __init__(self, parameters=Parameters()):
        self.par = parameters
        if (parameters.cutoff < 1 or parameters.mass <= 0
                or parameters.shape_mass <= 0 or parameters.shape_frequency <= 0
                or parameters.footprint <= 0 or parameters.sigma <= 0
                or parameters.coupling < 0 or parameters.energy <= 0):
            raise ValueError("Positive masses, widths, frequency and energy required")
        pairs = [(l, m) for l in range(parameters.cutoff + 1)
                 for m in range(-l, l + 1)]
        self.ell, self.order = np.array(pairs).T
        self.n = len(pairs)
        self.omega = self.ell + 0.5
        self.smoothing = np.exp(-parameters.footprint**2 * self.ell * (self.ell+1)/2)
        # D_i a are coefficients of (e_i x n).grad u = i L_i u.
        rows, cols, vals = [], [], []
        for j, (l, m) in enumerate(pairs):
            if m < l:
                rows.append(j+1)
                cols.append(j)
                vals.append(np.sqrt((l-m)*(l+m+1)))
        plus = sparse.csr_matrix((vals, (rows, cols)), shape=(self.n, self.n))
        minus = plus.T
        self.generators = (
            0.5j*(plus+minus), 0.5*(plus-minus),
            sparse.diags(1j*self.order, format="csr"),
        )
        self.conjugate_index = self.ell**2 + self.ell - self.order
        self.parity_m = (-1.0)**self.order

    def footprint(self, x):
        return self.smoothing * harmonics(self.ell, self.order, x).conj()

    def target_state(self, eta=0.7):
        """Smooth target displacement and directional velocity, energy-normalized."""
        f = (self.ell*(self.ell+1)
             * np.exp(-self.par.sigma**2*self.ell*(self.ell+1)/2)
             * harmonics(self.ell, self.order, [1., 0., 0.]).conj())
        v = -eta * (self.generators[2] @ f)
        scale = np.sqrt(self.par.energy / self.wave_energy(f, v))
        return f*scale, v*scale

    def initial(self, eta=0.7, dephased=False, seed=338):
        f, v = self.target_state(eta)
        parity = (-1.)**self.ell
        a, p = -parity*v/self.omega, parity*self.omega*f
        if dephased:
            phase = np.random.default_rng(seed).uniform(-np.pi, np.pi,
                                                        self.par.cutoff+1)[self.ell]
            a, p = (a*np.cos(phase)+p*np.sin(phase)/self.omega,
                    p*np.cos(phase)-self.omega*a*np.sin(phase))
        return State(a, p, np.array([1., 0., 0.]), np.zeros(3))

    def wave_energy(self, a, p):
        return float(0.5*(np.vdot(p, p).real + np.vdot(self.omega*a, self.omega*a).real))

    def interaction(self, s):
        w = self.footprint(s.x)
        mean = np.vdot(w, s.a).real
        mismatch = s.q - mean
        torque = self.par.coupling*mismatch*np.array(
            [np.vdot(w, d@s.a).real for d in self.generators])
        return w, mean, mismatch, torque

    def wave_angular(self, s):
        return -np.array([np.vdot(s.p, d@s.a).real for d in self.generators])

    def energy(self, s):
        _, _, e, _ = self.interaction(s)
        p = self.par
        return (self.wave_energy(s.a, s.p) + np.dot(s.angular, s.angular)/(2*p.mass)
                + s.pq**2/(2*p.shape_mass)
                + p.shape_mass*p.shape_frequency**2*s.q**2/2
                + p.coupling*e**2/2)

    def kick(self, s, dt):
        w, _, e, torque = self.interaction(s)
        s.p += dt*self.par.coupling*e*w
        s.pq -= dt*self.par.coupling*e
        s.angular += dt*torque

    def drift(self, s, dt):
        """Exact free field, shape oscillator and sphere-geodesic flows."""
        angle = self.omega*dt
        co, si = np.cos(angle), np.sin(angle)
        s.a, s.p = s.a*co+s.p*si/self.omega, s.p*co-self.omega*s.a*si
        mu, frequency = self.par.shape_mass, self.par.shape_frequency
        co, si = np.cos(frequency*dt), np.sin(frequency*dt)
        s.q, s.pq = s.q*co+s.pq*si/(mu*frequency), s.pq*co-mu*frequency*s.q*si
        length = np.linalg.norm(s.angular)
        if length:
            axis = s.angular/length
            angle = length*dt/self.par.mass
            co, si = np.cos(angle), np.sin(angle)
            s.x = (s.x*co + np.cross(axis, s.x)*si
                   + axis*np.dot(axis, s.x)*(1-co))

    def step(self, s, dt):
        self.kick(s, dt/2)
        self.drift(s, dt)
        self.kick(s, dt/2)

    def sample(self, a, points):
        return (harmonics(self.ell, self.order, points) @ a).real

    def evolve(self, state=None, dt=0.002, end=np.pi+1.2, output_dt=0.02,
               event_times=()):
        """Output coefficients, state, torque, and conservation diagnostics.

        Every integration step contributes to diagnostic maxima and torque
        integrals. The focus-window endpoints and optional event_times are
        exact integration boundaries. Displacement is relative to initial X.
        """
        if dt <= 0 or output_dt <= 0 or end <= 0:
            raise ValueError("Positive steps and horizon required")
        if any(not np.isfinite(t) or t < 0 or t > end for t in event_times):
            raise ValueError("Event times must lie in [0, end]")
        s = self.initial() if state is None else state.copy()
        initial_x = s.x.copy()
        times = np.unique(np.r_[np.arange(0, end, output_dt), end,
                                [t for t in [np.pi-.6, np.pi, np.pi+.6] if t <= end],
                                event_times])
        keys = ("t", "a", "p", "x", "angular", "q", "pq", "torque",
                "wave_angular", "energy", "shape_energy", "ubar")
        out = {k: [] for k in keys}
        initial_energy = self.energy(s)
        initial_j = s.angular+self.wave_angular(s)
        j_scale = max(initial_energy, np.linalg.norm(initial_j))
        stats = dict(max_energy_relative=0., max_angular_relative=0.,
                     max_sphere_error=0., max_tangency_error=0., max_reality_error=0.,
                     max_attachment_angular=0., torque_activity=0., focus_torque_activity=0.)

        def check():
            stats["max_energy_relative"] = max(stats["max_energy_relative"],
                    abs(self.energy(s)-initial_energy)/initial_energy)
            stats["max_angular_relative"] = max(stats["max_angular_relative"],
                    np.linalg.norm(s.angular+self.wave_angular(s)-initial_j)/j_scale)
            stats["max_sphere_error"] = max(stats["max_sphere_error"], abs(s.x@s.x-1))
            stats["max_tangency_error"] = max(stats["max_tangency_error"], abs(s.x@s.angular))
            stats["max_attachment_angular"] = max(stats["max_attachment_angular"],
                                                  float(np.linalg.norm(s.angular)))
            for a in [s.a, s.p]:
                stats["max_reality_error"] = max(stats["max_reality_error"],
                    float(np.max(np.abs(a[self.conjugate_index]-self.parity_m*a.conj()))))

        def record(t):
            _, mean, _, torque = self.interaction(s)
            vals = (t, s.a.copy(), s.p.copy(), s.x.copy(), s.angular.copy(), s.q,
                    s.pq, torque, self.wave_angular(s), self.energy(s),
                    s.pq**2/(2*self.par.shape_mass)
                    + self.par.shape_mass*self.par.shape_frequency**2*s.q**2/2, mean)
            for k, v in zip(keys, vals):
                out[k].append(v)

        t = 0.
        record(t)
        for target in times[1:]:
            steps = max(1, int(np.ceil((target-t)/dt-1e-10)))
            h = (target-t)/steps
            for _ in range(steps):
                old_torque = np.linalg.norm(self.interaction(s)[3])
                self.step(s, h)
                activity = h*(old_torque+np.linalg.norm(self.interaction(s)[3]))/2
                stats["torque_activity"] += activity
                if np.pi-.6 <= t+h/2 <= np.pi+.6:
                    stats["focus_torque_activity"] += activity
                t += h
                check()
            t = target
            record(t)
        arrays = {k: np.asarray(v) for k, v in out.items()}
        lo, hi = [np.argmin(abs(arrays["t"]-v)) for v in [np.pi-.6, np.pi+.6]]
        stats.update(initial_energy=initial_energy,
                     initial_wave_energy=self.wave_energy(out["a"][0], out["p"][0]),
                     initial_angular=initial_j.tolist(),
                     focus_impulse=(arrays["angular"][hi]-arrays["angular"][lo]).tolist(),
                     final_impulse=(s.angular-arrays["angular"][0]).tolist(),
                     max_displacement=float(np.max(np.arccos(np.clip(arrays["x"]@initial_x, -1, 1)))),
                     final_displacement=float(np.arccos(np.clip(s.x@initial_x, -1, 1))),
                     max_shape_energy=float(np.max(arrays["shape_energy"])),
                     focus_activity_fraction=stats["focus_torque_activity"]
                        /max(stats["torque_activity"], 1e-300))
        return arrays, stats
