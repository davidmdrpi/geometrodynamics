"""Two independent antipodal rotors coupled reciprocally to a round-S² scalar.

Effective mouth mechanics only: no scalar propagation on the drawn necks.
See docs/moving_mouth_pairs_protocol.md for the Hamiltonian and admission gates.
"""
from dataclasses import dataclass

import numpy as np
from scipy import special

from .antipodal_kick import AntipodalKick, Parameters


def harmonics(ell, order, points):
    """Pole-stable colatitude for tiny solved displacements away from ±z.

    arccos(z) loses the small transverse displacement when z rounds to one.
    Keep #339's helper unchanged for byte-identical replay of its old archive.
    """
    points = np.asarray(points)
    theta = np.arctan2(np.hypot(points[..., 0], points[..., 1]), points[..., 2])
    phi = np.arctan2(points[..., 1], points[..., 0])
    if hasattr(special, 'sph_harm_y'):
        return special.sph_harm_y(ell, order, theta[..., None], phi[..., None])
    return special.sph_harm(order, ell, phi[..., None], theta[..., None])


@dataclass
class PairState:
    a: np.ndarray
    p: np.ndarray
    x: np.ndarray
    angular: np.ndarray
    q: np.ndarray
    pq: np.ndarray
    mouth_impulse: np.ndarray

    def copy(self):
        return PairState(*(getattr(self, key).copy() for key in self.__dataclass_fields__))


class MovingMouthPairs(AntipodalKick):
    """M is the TOTAL inertia per antipodal pair; each mouth has shape mass mu."""

    def __init__(self, parameters=Parameters(mass=.10)):
        super().__init__(parameters)

    def initial(self, eta=.7, empty=False):
        profiles = [self.ell*(self.ell+1)
                    * np.exp(-self.par.sigma**2*self.ell*(self.ell+1)/2)
                    * harmonics(self.ell, self.order, x).conj()
                    for x in [[1., 0., 0.], [0., 0., 1.]]]
        f = profiles[0]+.65*profiles[1]
        v = -eta*(self.generators[2]@profiles[0]+.65*self.generators[0]@profiles[1])
        scale = np.sqrt(self.par.energy/self.wave_energy(f, v))
        parity = (-1.)**self.ell
        a, p = -parity*v*scale/self.omega, parity*self.omega*f*scale
        if empty:
            a[:] = 0
            p[:] = 0
        return PairState(a, p, np.array([[0., 0., 1.], [1., 0., 0.]]),
                         np.zeros((2, 3)), np.zeros((2, 2)), np.zeros((2, 2)),
                         np.zeros((2, 2, 3)))

    def interaction(self, s):
        mouths = s.x[:, None, :]*np.array([1., -1.])[None, :, None]
        w = self.smoothing*harmonics(self.ell, self.order, mouths).conj()
        mean = np.einsum('hsn,n->hs', w.conj(), s.a).real
        mismatch = s.q-mean
        derivatives = np.stack([d@s.a for d in self.generators])
        gradient = np.einsum('hsn,in->hsi', w.conj(), derivatives).real
        torque = self.par.coupling*mismatch[:, :, None]*gradient
        return w, mean, mismatch, torque

    def energy_components(self, s):
        p = self.par
        mismatch = self.interaction(s)[2]
        return dict(wave=self.wave_energy(s.a, s.p),
                    kinetic=np.sum(s.angular**2, axis=1)/(2*p.mass),
                    shape=s.pq**2/(2*p.shape_mass)
                          +p.shape_mass*p.shape_frequency**2*s.q**2/2,
                    coupling=p.coupling*mismatch**2/2)

    def energy(self, s):
        return float(sum(np.sum(v) for v in self.energy_components(s).values()))

    def kick(self, s, dt):
        w, _, mismatch, torque = self.interaction(s)
        s.p += dt*self.par.coupling*np.einsum('hs,hsn->n', mismatch, w)
        s.pq -= dt*self.par.coupling*mismatch
        s.angular += dt*torque.sum(axis=1)
        s.mouth_impulse += dt*torque

    def drift(self, s, dt):
        co, si = np.cos(self.omega*dt), np.sin(self.omega*dt)
        s.a, s.p = s.a*co+s.p*si/self.omega, s.p*co-self.omega*s.a*si
        mu, freq = self.par.shape_mass, self.par.shape_frequency
        co, si = np.cos(freq*dt), np.sin(freq*dt)
        s.q, s.pq = s.q*co+s.pq*si/(mu*freq), s.pq*co-mu*freq*s.q*si
        for k, angular in enumerate(s.angular):
            length = np.linalg.norm(angular)
            if length:
                axis = angular/length
                co, si = np.cos(length*dt/self.par.mass), np.sin(length*dt/self.par.mass)
                x = s.x[k]
                s.x[k] = x*co+np.cross(axis, x)*si+axis*(axis@x)*(1-co)

    def evolve(self, state=None, dt=.002, end=np.pi+1.2, output_dt=.04):
        if not all(np.isfinite(v) and v > 0 for v in (dt, end, output_dt)):
            raise ValueError('Finite positive time steps and horizon required')
        s = self.initial() if state is None else state.copy()
        x0, l0, impulse0 = s.x.copy(), s.angular.copy(), s.mouth_impulse.copy()
        e0 = self.energy(s)
        j0 = self.wave_angular(s)+s.angular.sum(axis=0)
        scale = max(e0, np.linalg.norm(j0), 1e-12)
        events = [t for t in [.6, np.pi-.6, np.pi, np.pi+.6] if t <= end]
        times = np.unique(np.r_[np.arange(0, end, output_dt), events, end])
        keys = ['t', 'a', 'x', 'angular', 'q', 'pq', 'mouth_impulse', 'mouth_torque',
                'wave_angular', 'energy', 'wave_energy', 'kinetic_energy',
                'shape_energy', 'coupling_energy', 'separation']
        out = {k: [] for k in keys}
        stats = dict(max_energy_relative=0., max_angular_relative=0.,
                     max_sphere_error=0., max_tangency_error=0., max_reality_error=0.,
                     max_pair_angular=0., max_separation_change=0.,
                     minimum_cross_mouth_distance=np.pi/2)
        sep0 = np.arccos(np.clip(s.x[0]@s.x[1], -1, 1))

        def check():
            stats['max_energy_relative'] = max(stats['max_energy_relative'],
                    abs(self.energy(s)-e0)/max(e0, 1e-12))
            stats['max_angular_relative'] = max(stats['max_angular_relative'],
                    np.linalg.norm(self.wave_angular(s)+s.angular.sum(axis=0)-j0)/scale)
            stats['max_sphere_error'] = max(stats['max_sphere_error'],
                    float(np.max(abs(np.sum(s.x*s.x, axis=1)-1))))
            stats['max_tangency_error'] = max(stats['max_tangency_error'],
                    float(np.max(abs(np.sum(s.x*s.angular, axis=1)))))
            stats['max_pair_angular'] = max(stats['max_pair_angular'],
                    float(np.max(np.linalg.norm(s.angular, axis=1))))
            sep = np.arccos(np.clip(s.x[0]@s.x[1], -1, 1))
            stats['max_separation_change'] = max(stats['max_separation_change'], abs(sep-sep0))
            stats['minimum_cross_mouth_distance'] = min(stats['minimum_cross_mouth_distance'],
                                                        sep, np.pi-sep)
            for a in [s.a, s.p]:
                stats['max_reality_error'] = max(stats['max_reality_error'],
                    float(np.max(abs(a[self.conjugate_index]-self.parity_m*a.conj()))))

        def record(t):
            energies = self.energy_components(s)
            vals = [t, s.a, s.x, s.angular, s.q, s.pq, s.mouth_impulse,
                    self.interaction(s)[3], self.wave_angular(s), self.energy(s),
                    energies['wave'], energies['kinetic'], energies['shape'], energies['coupling'],
                    np.arccos(np.clip(s.x[0]@s.x[1], -1, 1))]
            for key, val in zip(keys, vals):
                out[key].append(np.array(val).copy())

        record(0.)
        t = 0.
        for target in times[1:]:
            steps = max(1, int(np.ceil((target-t)/dt-1e-10)))
            h = (target-t)/steps
            for _ in range(steps):
                self.step(s, h)
                check()
            t = target
            record(t)
        arrays = {k: np.asarray(v) for k, v in out.items()}
        stats.update(initial_energy=e0, initial_wave_energy=arrays['wave_energy'][0].item(),
                     initial_angular=j0.tolist(), final_pair_impulse=(s.angular-l0).tolist(),
                     final_mouth_impulse=(s.mouth_impulse-impulse0).tolist(),
                     final_displacement=np.arccos(np.clip(np.sum(s.x*x0, axis=1), -1, 1)).tolist())
        if end >= np.pi+.6:
            lo, hi, launch = [np.flatnonzero(times == t)[0] for t in [np.pi-.6, np.pi+.6, .6]]
            stats['focus_mouth_impulse'] = (arrays['mouth_impulse'][hi]-arrays['mouth_impulse'][lo]).tolist()
            stats['launch_mouth_impulse'] = (arrays['mouth_impulse'][launch]-impulse0).tolist()
        return arrays, stats
