"""A real wave on S² moving two antipodal mouth pairs, one throat loaded.

The model of ``moving_mouth_pairs`` (PR #340) with three changes, so that the
sphere-with-handles drawing can show what ``viz/wave_mouth_network`` shows on
a ring:

* pair kinetic energy is relativistic, ``sqrt(M² + |L|²) − M``;
* a LOADED mouth has the metastable support of ``viz/wave_mouth_network``,
  open at ``q = 0`` and collapsed at ``q = −1``, holding ``E_s``; the other
  mouths keep #340's harmonic shape well;
* every term of the energy is recorded.

    H = ½ Σ (|p_lm|² + ω_l² |a_lm|²)                      ω_l = l + ½
        + Σ_pairs [ sqrt(M² + |L_h|²) − M ]
        + Σ_mouths [ pq²/(2μ) + S(q) + κ/2 (q − <u>)² ]

``<u>`` is the field averaged over the mouth footprint. As in #340 the scalar
lives on the ROUND sphere: nothing propagates along the drawn necks, the
support is assumed, and this is not a GR solution.
"""
from dataclasses import dataclass, replace

import numpy as np

from .antipodal_kick import Parameters
from .moving_mouth_pairs import MovingMouthPairs, PairState, harmonics

COMPONENTS = ("wave", "kinetic", "shape", "coupling")


@dataclass(frozen=True)
class Support:
    stored: tuple = ((2.0, 0.0), (0.0, 0.0))    # E_s per [pair][mouth]
    barrier: float = 0.22


PARAMETERS = Parameters(cutoff=24, sigma=0.22, footprint=0.15, mass=6.0,
                        shape_mass=0.6, shape_frequency=(4 / 0.6) ** 0.5,
                        coupling=3.0, energy=0.6)
# The packet starts 0.3 rad from the antipode of the loaded mouth (+z), so it
# refocuses beside that mouth near t = pi and pushes it off centre.
SOURCE = (0.3 * np.cos(0.7), 0.3 * np.sin(0.7), -1.0)
AMPLITUDE = 1.0
SUBTHRESHOLD = 0.6
DURATION = 30.0
STEP = 0.005


class WaveMouthSphere(MovingMouthPairs):
    def __init__(self, parameters=PARAMETERS, support=Support()):
        super().__init__(parameters)
        self.stored = np.asarray(support.stored, float)
        self.barrier = support.barrier
        if self.stored.shape != (2, 2) or (self.stored < 0).any():
            raise ValueError("stored energy: one nonnegative number per mouth")
        if self.barrier <= 3 / 16:
            raise ValueError("barrier must exceed 3/16 for a metastable well")

    # ---- support ---------------------------------------------------------

    def support(self, q):
        u, h = 1.0 + q, self.barrier
        loaded = self.stored * (16 * h * u ** 2 * (1 - u) ** 2 + u ** 2 * (3 - 2 * u))
        p = self.par
        plain = 0.5 * p.shape_mass * p.shape_frequency ** 2 * q ** 2
        return np.where(self.stored > 0, loaded, plain)

    def support_force(self, q):
        u, h = 1.0 + q, self.barrier
        loaded = -self.stored * (32 * h * u * (1 - u) * (1 - 2 * u) + 6 * u * (1 - u))
        p = self.par
        return np.where(self.stored > 0, loaded,
                        -p.shape_mass * p.shape_frequency ** 2 * q)

    def barrier_top(self):
        u = np.linspace(0, 1, 2001)
        h = self.barrier
        well = 16 * h * u ** 2 * (1 - u) ** 2 + u ** 2 * (3 - 2 * u)
        peak = (well[1:-1] > well[:-2]) & (well[1:-1] > well[2:])
        return float(u[1:-1][peak][0] - 1.0)

    # ---- state -----------------------------------------------------------

    def start(self, amplitude=AMPLITUDE, source=SOURCE):
        """An outgoing ring wave released from ``source``; mouths at rest."""
        c = np.asarray(source, float)
        c = c / np.linalg.norm(c)
        ll = self.ell * (self.ell + 1)
        f = np.exp(-self.par.sigma ** 2 * ll / 2) * harmonics(
            self.ell, self.order, c).conj().ravel()
        f = f - f[0] * (self.ell == 0)            # no uniform offset
        v = np.zeros_like(f)
        if amplitude > 0:
            f = f * amplitude * np.sqrt(self.par.energy / self.wave_energy(f, v))
        else:
            f = 0 * f
        return PairState(f, v, np.array([[0., 0., 1.], [1., 0., 0.]]),
                         np.zeros((2, 3)), np.zeros((2, 2)), np.zeros((2, 2)),
                         np.zeros((2, 2, 3)))

    # ---- energy ----------------------------------------------------------

    def energy_components(self, s):
        p = self.par
        mismatch = self.interaction(s)[2]
        speed = np.sum(s.angular ** 2, axis=1)
        return dict(wave=self.wave_energy(s.a, s.p),
                    kinetic=np.sqrt(p.mass ** 2 + speed) - p.mass,
                    shape=s.pq ** 2 / (2 * p.shape_mass) + self.support(s.q),
                    coupling=p.coupling * mismatch ** 2 / 2)

    # ---- integrator: support in the kick, free motion in the drift --------

    def kick(self, s, dt):
        super().kick(s, dt)
        s.pq += dt * self.support_force(s.q)

    def drift(self, s, dt):
        co, si = np.cos(self.omega * dt), np.sin(self.omega * dt)
        s.a, s.p = s.a * co + s.p * si / self.omega, s.p * co - self.omega * s.a * si
        s.q = s.q + dt * s.pq / self.par.shape_mass
        for k, angular in enumerate(s.angular):
            length = np.linalg.norm(angular)
            if length:
                axis = angular / length
                angle = length * dt / np.sqrt(self.par.mass ** 2 + length ** 2)
                co, si = np.cos(angle), np.sin(angle)
                x = s.x[k]
                s.x[k] = x * co + np.cross(axis, x) * si + axis * (axis @ x) * (1 - co)

    def wave_energy_by_degree(self, a, p):
        """Field energy in each spherical-harmonic degree l."""
        e = 0.5 * (np.abs(p) ** 2 + (self.omega * np.abs(a)) ** 2)
        return np.bincount(self.ell, weights=e, minlength=self.par.cutoff + 1)

    def speed(self, angular):
        """Mouth speed on the unit sphere, in units of the wave speed."""
        length = np.linalg.norm(angular, axis=-1)
        return length / np.sqrt(self.par.mass ** 2 + length ** 2)

    def run(self, state=None, duration=DURATION, h=STEP, sample=20):
        s = self.start() if state is None else state.copy()
        steps = int(round(duration / h))
        if steps % sample:
            raise ValueError("duration must hold a whole number of samples")
        keys = ("t", "a", "p", "x", "angular", "q") + COMPONENTS
        out = {k: [] for k in keys}

        def record(t):
            c = self.energy_components(s)
            for key, val in zip(keys, (t, s.a, s.p, s.x, s.angular, s.q,
                                       *(c[n] for n in COMPONENTS))):
                out[key].append(np.array(val).copy())

        record(0.0)
        j0 = self.wave_angular(s) + s.angular.sum(axis=0)
        drift = 0.0
        for n in range(steps):
            self.step(s, h)
            if (n + 1) % sample == 0:
                record((n + 1) * h)
                drift = max(drift, float(np.linalg.norm(
                    self.wave_angular(s) + s.angular.sum(axis=0) - j0)))
        data = {k: np.asarray(v) for k, v in out.items()}
        data["total"] = (data["wave"] + data["kinetic"].sum(axis=1)
                         + data["shape"].sum(axis=(1, 2))
                         + data["coupling"].sum(axis=(1, 2)))
        data["angular_momentum_drift"] = drift
        return data


def scenario(amplitude=AMPLITUDE, stored=((2.0, 0.0), (0.0, 0.0)), duration=DURATION,
             h=STEP, sample=20, parameters=PARAMETERS):
    model = WaveMouthSphere(parameters, Support(stored=stored))
    return model, model.run(model.start(amplitude), duration, h, sample)


# ---- outward energy flux through a fixed circle on the sphere ---------------

CAP_CENTRE = (0.15, 0.29, 0.945)   # fixed; chosen so the cap holds mouth 1 only
CAP_RADIUS = 0.85                  # radians
CAP_WINDOW = 15.0                  # after this the mouths come too close to it


def _frame(centre):
    c = np.asarray(centre, float)
    c = c / np.linalg.norm(c)
    e1 = np.cross(c, [1.0, 0.0, 0.0] if abs(c[0]) < 0.9 else [0.0, 1.0, 0.0])
    e1 /= np.linalg.norm(e1)
    return c, e1, np.cross(c, e1)


class CapFlux:
    """Field energy in a fixed spherical cap and the flux through its edge.

    For ``u_tt = (Δ − ¼) u + f`` the density ``e = ½(u_t² + |∇u|² + ¼u²)``
    obeys ``e_t + div S = u_t f`` with ``S = −u_t ∇u``. The field is
    band-limited, so both quadratures below are exact: equispaced in azimuth,
    Gauss–Legendre in the cosine of the polar angle about the cap centre.
    ``∇u`` comes from the rotation generators, ``n × ∇u = (D_i a)_i``.
    """

    def __init__(self, model, centre=CAP_CENTRE, radius=CAP_RADIUS):
        self.model, self.radius = model, float(radius)
        c, e1, e2 = _frame(centre)
        self.centre = c
        band = model.par.cutoff
        n_phi = 4 * band + 8
        phi = 2 * np.pi * np.arange(n_phi) / n_phi
        ring = np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2
        edge = np.cos(radius) * c + np.sin(radius) * ring
        self.edge = edge
        self.normal = (edge * np.cos(radius) - c) / np.sin(radius)
        self.edge_weight = np.sin(radius) * 2 * np.pi / n_phi
        self.edge_y = harmonics(model.ell, model.order, edge)
        nodes, weights = np.polynomial.legendre.leggauss(band + 4)
        lo = np.cos(radius)
        mu = 0.5 * (1 - lo) * nodes + 0.5 * (1 + lo)
        w_mu = 0.5 * (1 - lo) * weights
        pts = (mu[:, None, None] * c
               + np.sqrt(1 - mu ** 2)[:, None, None] * ring[None, :, :])
        self.cap = pts.reshape(-1, 3)
        self.cap_weight = np.repeat(w_mu, n_phi) * 2 * np.pi / n_phi
        self.cap_y = harmonics(model.ell, model.order, self.cap)

    def _fields(self, y, points, a, p):
        g = np.stack([(y @ (d @ a)).real for d in self.model.generators], axis=1)
        return (y @ a).real, (y @ p).real, g

    def flux(self, a, p):
        """Outward energy flux through the edge, ``−∮ u_t (ν·∇u) dl``."""
        _, ut, g = self._fields(self.edge_y, self.edge, a, p)
        grad = np.cross(g, self.edge)
        return float(-self.edge_weight * np.sum(
            ut * np.einsum("ij,ij->i", self.normal, grad)))

    def energy(self, a, p):
        u, ut, g = self._fields(self.cap_y, self.cap, a, p)
        density = 0.5 * (ut ** 2 + np.sum(g ** 2, axis=1) + 0.25 * u ** 2)
        return float(self.cap_weight @ density)


def cap_flux_ledger(amplitude=AMPLITUDE, stored=((2.0, 0.0), (0.0, 0.0)),
                    duration=CAP_WINDOW, h=STEP, sample=20, centre=CAP_CENTRE,
                    radius=CAP_RADIUS, parameters=PARAMETERS):
    """Time-integrated flux out of a fixed cap around the loaded mouth.

    ``inside`` is the cap's field energy, plus mouth 1's throat and coupling
    energy, plus the work mouth 1 has done on its pair's rotor (the rotor is
    shared with the antipodal mouth, which is outside). The check is

        net_outward(t) = inside(0) − inside(t).
    """
    model = WaveMouthSphere(parameters, Support(stored=stored))
    cap = CapFlux(model, centre, radius)
    s = model.start(amplitude)
    steps = int(round(duration / h))
    mass = model.par.mass

    def rotor_power():
        torque = model.interaction(s)[3][0, 0]
        length = np.linalg.norm(s.angular[0])
        return float(torque @ s.angular[0] / np.sqrt(mass ** 2 + length ** 2))

    def inside(work):
        c = model.energy_components(s)
        return (cap.energy(s.a, s.p), float(c["shape"][0, 0] + c["coupling"][0, 0]),
                work)

    def distances():
        mouths = np.array([s.x[0], -s.x[0], s.x[1], -s.x[1]])
        return np.arccos(np.clip(mouths @ cap.centre, -1, 1))

    flux, power = cap.flux(s.a, s.p), rotor_power()
    net = work = gross_out = gross_in = 0.0
    rows = [(0.0, 0.0, *inside(0.0), *distances())]
    for n in range(steps):
        model.step(s, h)
        new_flux, new_power = cap.flux(s.a, s.p), rotor_power()
        inc = 0.5 * h * (flux + new_flux)
        net += inc
        gross_out += max(inc, 0.0)
        gross_in += max(-inc, 0.0)
        work += 0.5 * h * (power + new_power)
        flux, power = new_flux, new_power
        if (n + 1) % sample == 0:
            rows.append(((n + 1) * h, net, *inside(work), *distances()))
    r = np.array(rows)
    total_inside = r[:, 2] + r[:, 3] + r[:, 4]
    residual = r[:, 1] - (total_inside[0] - total_inside)
    return {
        "centre": cap.centre.tolist(), "radius": radius, "duration": duration,
        "step": h, "time": r[:, 0], "net_outward": r[:, 1],
        "cap_field_energy": r[:, 2], "mouth_energy": r[:, 3],
        "rotor_work": r[:, 4], "mouth_distance": r[:, 5:],
        "net_outward_final": float(net), "net_outward_minimum": float(r[:, 1].min()),
        "gross_outward": gross_out, "gross_inward": gross_in,
        "inside_initial": float(total_inside[0]),
        "inside_final": float(total_inside[-1]),
        "cap_field_energy_final": float(r[-1, 2]),
        "rotor_work_final": float(work),
        "balance_max_residual": float(np.abs(residual).max()),
        "loaded_mouth_margin": float((radius - r[:, 5]).min()),
        "other_mouth_margin": float((r[:, 6:] - radius).min()),
    }
