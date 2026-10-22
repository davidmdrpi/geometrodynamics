"""Pair birth under the antipodal-only (BAM) interaction rule.

`pair_birth.PairBirth` lets every part of the field act on the mouths at every instant:
the Coulomb attraction, the exchange and self fields, and with `NeckChannel` the ports.
The BAM routing rule adopted in #347 restricts the interaction support instead:

    a front launched at (n, t) can next interact only at its antipodal focus (-n, t + pi).

This module rebuilds #345's pair history under that rule. Waves become **fronts**:
characteristic packets with an origin n, a launch time t, an energy E, an angular momentum
J about the sphere's centre (conserved in flight) and a complex phase. Mouths a front passes
on its way feel nothing.

- **Reflections.** If no open mouth's aperture contains the focus, the front passes through
  the caustic and continues as a front launched from (-n, t + pi). The 2-D Gouy phase
  multiplies its phase by -i, so after two foci it is back sign-flipped. That is what the
  field does exactly for the conformal operator (`pair_birth.bare_timeline`).
- **Mouth interactions.** If an open mouth's aperture (angular radius a) contains the focus at
  the focal time, the mouth captures the front there:
  - J's component across the mouth's normal goes into its orbital angular momentum, and the
    component along the normal into its spin. J is conserved exactly.
  - The rest of the front's energy, E less the mouth's kinetic and spin energy gain, enters
    the throat. It cannot be negative, because |J| <= E for a front.
- **Throat interactions.** The front crosses the neck in neck time L = sin 2a (#338):
  - It enters at sigma_in = tau_k, the capturing mouth's proper time since birth. It leaves
    the other mouth when tau_j = sigma_in + L, on the actual clock histories (#347's map).
  - It leaves from the aperture point mapped through the handle. In the sphere-oriented rim
    frames used here (e1 to the birth contact point, e2 = x cross e1), a Klein handle is a
    rotation, (u, v) -> (-u, -v), and a torus handle a reflection, (u, v) -> (-u, v). The
    point inside the aperture is mapped as a portal, which is an assumption: the glued
    handle's throat attaches only at the rim.
  - It leaves as an isotropic front (J = 0, no recoil) with the energy and phase it carried.
  - **Advanced exits.** If the exit lies in the destination's past (the offset exceeds L),
    the exit is an advanced edge. Nothing on its way interacts under the rule, so inserting
    it is consistent while its first focus is still ahead. Otherwise it is recorded as a
    causality boundary, with its energy booked as unresolved and never as absorbed.
  - Fronts still inside the throat at annihilation are released with the pair.
- **Static interaction.** The Coulomb attraction acts at every separation, so the strict rule
  removes it (static='none'). The flux energy is then part of the fixed rest energy, and the
  mouths coast on geodesics between captures. static='coulomb' keeps it, as a labelled
  control that breaks the rule for the static field only.
- **Emission.** Fronts are launched only by events, never continuously:
  - the birth flash: one front per mouth from its birth position, with the energy of #345's
    flash (the free wave of -Q_qs; cross terms split evenly) and its angular momentum, which
    is zero because the flash starts at rest;
  - throat exits;
  - the release pulse at rim contact (#345's rule). It refocuses at the antipode after pi,
    where the next pair is re-born, carrying the pulse's energy and angular momentum less
    those of the new flash.
- **Clocks.** Each mouth's proper time since its own birth, as in #345.

Everything is classical. The amplitude phases are bookkeeping for interference, which is not
modelled. The front-capture rules are stated rules, not derived from a field solution.
"""
from dataclasses import dataclass, field

import numpy as np
from scipy import optimize

from .pair_birth import PairBirth, Parameters, birth_frames, birth_positions, rotate


@dataclass(eq=False)
class Front:
    origin: np.ndarray       # where it was launched or last focused
    t: float                 # launch time
    energy: float
    J: np.ndarray            # angular momentum about the centre
    phase: complex = 1.
    kind: str = ''
    foci: int = 0            # caustics passed


@dataclass(eq=False)
class ThroatFront:
    source: int
    dest: int
    sigma_out: float         # the destination's proper time at exit
    local: np.ndarray        # aperture coordinates in the destination's frame (already mapped)
    energy: float
    phase: complex
    kind: str
    t_in: float


@dataclass
class AntipodalState:
    t: float
    x: np.ndarray
    L: np.ndarray            # mechanical angular momenta x cross p
    tau: np.ndarray
    frames: np.ndarray
    spin: np.ndarray
    alive: bool = True
    birth: float = 0.
    cycle: int = 0
    fronts: list = field(default_factory=list)
    throat: list = field(default_factory=list)
    pulse: dict = None
    unresolved: float = 0.

    def copy(self):
        return AntipodalState(self.t, self.x.copy(), self.L.copy(), self.tau.copy(), self.frames.copy(),
                              self.spin.copy(), self.alive, self.birth, self.cycle,
                              [Front(f.origin.copy(), f.t, f.energy, f.J.copy(), f.phase, f.kind, f.foci)
                               for f in self.fronts],
                              [ThroatFront(f.source, f.dest, f.sigma_out, f.local.copy(), f.energy, f.phase, f.kind,
                                           f.t_in) for f in self.throat],
                              None if self.pulse is None else dict(self.pulse), self.unresolved)


def angular_distance(x, y):
    return float(np.arccos(np.clip(np.dot(x, y), -1, 1)))


def local_coordinates(centre, frame, point):
    """Azimuthal-equidistant coordinates of `point` in the tangent frame at `centre`."""
    r = angular_distance(centre, point)
    tangent = point - np.dot(point, centre)*centre
    norm = np.linalg.norm(tangent)
    if norm < 1e-15:
        return np.zeros(2)
    return r*np.array([tangent @ frame[0], tangent @ frame[1]])/norm


def from_local(centre, frame, local):
    r = np.linalg.norm(local)
    if r == 0:
        return centre.copy()
    direction = (local[0]*frame[0] + local[1]*frame[1])/r
    p = np.cos(r)*centre + np.sin(r)*direction
    return p/np.linalg.norm(p)


class AntipodalPairBirth:
    def __init__(self, parameters=Parameters(), static='none', handle='klein', reflection_phase=-1j):
        if static not in ('none', 'coulomb') or handle not in ('klein', 'torus'):
            raise ValueError("static must be 'none' or 'coulomb', handle 'klein' or 'torus'")
        self.par = parameters
        self.field = PairBirth(parameters)       # flash, Coulomb energy and charges
        self.static = static
        self.handle = handle
        self.map = np.array([-1., -1.]) if handle == 'klein' else np.array([-1., 1.])
        self.reflection_phase = complex(reflection_phase)
        a = parameters.mouth_angle
        self.length = np.sin(2*a)
        self.inertia = parameters.mass*np.sin(a)**2

    # ---- births and fronts ----------------------------------------------------------------------
    def flash_fronts(self, x, vel, t):
        """#345's birth flash, the free wave of -Q_qs launched at rest, split into one front per
        mouth: energy (self part plus half the cross term) and angular momentum."""
        f = self.field
        Q = -f.quasi_static(x, vel)
        e = [.5*np.vdot(f.omega*Q[k], f.omega*Q[k]).real for k in range(2)]
        cross = np.vdot(f.omega*Q[0], f.omega*Q[1]).real
        zero = np.zeros(f.n, complex)
        J = [f.field_angular(Q[k], zero) for k in range(2)]           # zero: the flash starts at rest
        return [Front(x[k].copy(), t, e[k] + cross/2, J[k], 1., f'flash {"AB"[k]}') for k in range(2)]

    def initial(self, v=1/3, plane=((0., 0., 1.), (1., 0., 0.)), velocities=None):
        """Mouths in rim contact at P, kicked apart with speed v each (or `velocities`). As in
        #345 the mouths pay for the birth flash: it is emitted at the birth event, and the
        momenta are corrected so that pair plus flash carry the kicked pair's energy and
        angular momentum."""
        x, vel, point, _ = birth_positions(self.par.mouth_angle, v, plane)
        if velocities is not None:
            vel = np.asarray(velocities, float)
            vel = vel - np.sum(vel*x, -1, keepdims=True)*x
        speed = np.linalg.norm(vel, axis=1)
        if np.any(speed >= 1):
            raise ValueError('Mouth speeds must be below 1')
        M = self.par.mass
        p = M*vel/np.sqrt(1-speed**2)[:, None]
        U = self.field.coulomb(x)[0] if self.static == 'coulomb' else 0.
        E = np.sqrt(M**2 + np.sum(p*p, 1)).sum() + U
        p, frames, flash = self.rebirth_momenta(0., x, point, E, np.cross(x, p).sum(0), p)
        s = AntipodalState(0., x, np.cross(x, p), np.zeros(2), frames, np.zeros(2))
        s.fronts = flash
        return s

    # ---- mechanics ---------------------------------------------------------------------------
    def energies(self, s):
        M = self.par.mass
        E = np.sqrt(M**2 + np.sum(s.L**2, 1))
        return E

    def velocity(self, s):
        return np.cross(s.L, s.x)/self.energies(s)[:, None]

    def ledger(self, s):
        out = dict(rest=0., kinetic=0., spin=0., coulomb=0., fronts=float(sum(f.energy for f in s.fronts)),
                   throat=float(sum(f.energy for f in s.throat)), pulse=0. if s.pulse is None else s.pulse['energy'],
                   unresolved=s.unresolved)
        if s.alive:
            E = self.energies(s)
            out.update(rest=2*self.par.mass, kinetic=float(E.sum() - 2*self.par.mass),
                       spin=float(np.sum(s.spin**2)/(2*self.inertia)))
            if self.static == 'coulomb':
                out['coulomb'] = float(self.field.coulomb(s.x)[0])
        return out

    def angular_momentum(self, s):
        J = sum((f.J for f in s.fronts), np.zeros(3))
        if s.pulse is not None:
            J = J + s.pulse['J']
        if s.alive:
            J = J + s.L.sum(0) + s.spin @ s.x
        return J

    def separation(self, s):
        return angular_distance(s.x[0], s.x[1])

    def _kick(self, s, h):
        if self.static == 'coulomb':
            s.L += h*np.cross(s.x, self.field.coulomb(s.x)[1])

    def _gyro(self, s, h):
        if np.any(s.spin):
            s.L = rotate(s.L, s.x, s.spin*h/self.energies(s))

    def step(self, s, h):
        """Kick, gyroscopic turn, exact geodesic drift (frames carried, spin turning them),
        gyroscopic turn, kick. Without the static field and spin it is exact free motion."""
        if not s.alive:
            s.t += h
            return
        self._kick(s, h/2)
        self._gyro(s, h/2)
        E = self.energies(s)
        M = self.par.mass
        for k in range(2):
            length = np.linalg.norm(s.L[k])
            if length > 0:
                axis, angle = s.L[k]/length, length*h/E[k]
                s.x[k] = rotate(s.x[k], axis, angle)
                s.frames[k] = rotate(s.frames[k], axis, angle)
        s.x /= np.linalg.norm(s.x, axis=1, keepdims=True)
        e1 = rotate(s.frames[:, 0], s.x, s.spin*h*M/(self.inertia*E))
        e1 -= np.sum(e1*s.x, 1, keepdims=True)*s.x
        e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
        s.frames = np.stack([e1, np.cross(s.x, e1)], 1)
        s.tau += h*M/E
        self._gyro(s, h/2)
        self._kick(s, h/2)
        s.t += h

    # ---- interactions at foci ---------------------------------------------------------------
    def capture(self, s, front, k, focus):
        """Mouth k captures `front` at its focus: J into the mouth, the remaining energy into the
        throat. Returns the energy that entered the throat and the aperture coordinates."""
        before = self.energies(s)[k] + s.spin[k]**2/(2*self.inertia)
        x = s.x[k]
        s.L[k] += front.J - np.dot(front.J, x)*x
        s.spin[k] += np.dot(front.J, x)
        after = self.energies(s)[k] + s.spin[k]**2/(2*self.inertia)
        inside = front.energy - (after - before)
        if inside < -1e-12:
            raise RuntimeError('Captured front cannot pay for the mouth it pushes')
        return max(inside, 0.), local_coordinates(x, s.frames[k], focus)

    def rebirth_momenta(self, t, x, point, E_target, J_target, p0):
        """Least-squares momenta at a birth at time t: the pair and its flash carry exactly the
        given energy and angular momentum (cf. PairBirth.recreate)."""
        M = self.par.mass
        frames = birth_frames(x, point)
        basis = [(frames[k][0], frames[k][1]) for k in range(2)]
        U = self.field.coulomb(x)[0] if self.static == 'coulomb' else 0.

        def momenta(b):
            return np.array([p0[k] + b[2*k]*basis[k][0] + b[2*k+1]*basis[k][1] for k in range(2)])

        def residual(b):
            p = momenta(b)
            E = np.sqrt(M**2 + np.sum(p*p, 1))
            flash = self.flash_fronts(x, p/E[:, None], t)
            Ef = sum(f.energy for f in flash)
            Jf = sum((f.J for f in flash), np.zeros(3))
            return np.r_[np.cross(x, p).sum(0) + Jf - J_target, (E.sum() + U + Ef - E_target)/E_target, 1e-8*b]

        b = optimize.least_squares(residual, np.zeros(4), xtol=1e-15, ftol=1e-15, gtol=1e-15).x
        if abs(residual(b)[3]) > 1e-9:
            raise RuntimeError('Re-birth cannot carry the pulse energy')
        p = momenta(b)
        E = np.sqrt(M**2 + np.sum(p*p, 1))
        return p, frames, self.flash_fronts(x, p/E[:, None], t)

    # ---- the history --------------------------------------------------------------------------
    def run(self, s, end, cycles=1, dt=.004, output_dt=.02):
        """Event-driven history: free (or Coulomb) flight between events, and as events the foci
        (capture or reflection), throat exits, rim contact with release, and re-birth at the
        pulse's antipodal focus, up to `cycles` births. Returns records and an event list."""
        s = s.copy()
        keys = ('t', 'x', 'L', 'd', 'tau', 'alive', 'cycle', 'energy', 'J', 'spin', 'n_fronts', 'n_throat')
        out = {k: [] for k in keys}
        out['ledger'] = {k: [] for k in self.ledger(s)}
        events = []
        track = self._new_track(s)
        rebirth = None

        def record():
            e = self.ledger(s)
            for k, v in e.items():
                out['ledger'][k].append(v)
            vals = (s.t, s.x.copy(), s.L.copy(), self.separation(s), s.tau.copy(), s.alive, s.cycle, sum(e.values()),
                    self.angular_momentum(s), s.spin.copy(), len(s.fronts), len(s.throat))
            for k, v in zip(keys, vals):
                out[k].append(v)

        record()
        next_out = output_dt
        while s.t < end - 1e-12:
            focus_time = min((f.t + np.pi for f in s.fronts), default=np.inf)
            target = min(end, focus_time, rebirth if rebirth is not None else np.inf, next_out)
            h = target - s.t
            if s.alive:
                h = min(h, dt)
            prev = s.copy()
            self.step(s, h)
            if s.alive:
                event = self._continuous_event(prev, s, h)
                if event is not None:
                    kind, h_event = event
                    s = prev.copy()
                    self.step(s, h_event)
                    self._extend(track, s)
                    record()
                    if kind == 'contact':
                        rebirth = self._release(s, events, cycles)
                    else:
                        self._exit(s, s.throat.pop(kind), events)
                    record()
                    continue
                self._extend(track, s)
            if abs(s.t - focus_time) < 1e-9:
                record()
                for f in [f for f in s.fronts if abs(f.t + np.pi - s.t) < 1e-9]:
                    self._focus(s, f, track, events)
                record()
            if rebirth is not None and abs(s.t - rebirth) < 1e-9:
                self._rebirth(s, events)
                rebirth = None
                track = self._new_track(s)
                record()
            if s.t >= next_out - 1e-9:
                record()
                next_out = (np.floor(s.t/output_dt + 1e-9) + 1)*output_dt
        arrays = {k: np.asarray(v) for k, v in out.items() if k != 'ledger'}
        arrays['ledger'] = {k: np.asarray(v) for k, v in out['ledger'].items()}
        alive = arrays['alive']
        arrays['clock'] = np.where(alive, arrays['tau'][:, 0] - arrays['tau'][:, 1], np.nan)
        e, J = arrays['energy'], arrays['J']
        margin = np.where(alive, arrays['d'] + self.length - np.abs(arrays['clock']), np.nan)
        info = dict(events=events, energy_drift=float(np.abs(e - e[0]).max()/e[0]),
                    angular_drift=float(np.abs(J - J[0]).max()), final_state=s, length=self.length,
                    min_chronology_margin=float(np.nanmin(margin)))
        return arrays, info

    @staticmethod
    def _new_track(s):
        return dict(t=[s.t], tau=[s.tau.copy()], x=[s.x.copy()], frames=[s.frames.copy()])

    @staticmethod
    def _extend(track, s):
        for key, val in (('t', s.t), ('tau', s.tau.copy()), ('x', s.x.copy()), ('frames', s.frames.copy())):
            track[key].append(val)

    def _continuous_event(self, prev, s, h):
        """The earliest rim contact or throat exit inside the step just taken, if any."""
        a = self.par.mouth_angle
        candidates = []
        closing = np.dot(self.velocity(prev)[0] - self.velocity(prev)[1], prev.x[0] - prev.x[1]) < 0
        if self.separation(s) <= 2*a and closing:
            candidates.append(('contact', self._root(prev, h, lambda q: self.separation(q) - 2*a)))
        for i, tf in enumerate(s.throat):
            if s.tau[tf.dest] >= tf.sigma_out:
                candidates.append((i, self._root(prev, h, lambda q, tf=tf: q.tau[tf.dest] - tf.sigma_out)))
        return min(candidates, key=lambda c: c[1]) if candidates else None

    def _release(self, s, events, cycles):
        """Rim contact: the pair, its spin and the throat's contents leave as one pulse."""
        E = self.ledger(s)
        released = E['rest'] + E['kinetic'] + E['spin'] + E['coulomb'] + E['throat']
        c = (s.x[0] + s.x[1])/np.linalg.norm(s.x[0] + s.x[1])
        events.append(dict(kind='contact', time=s.t, cycle=s.cycle, point=c.tolist(), clocks=s.tau.tolist(),
                           energy=released, throat_released=E['throat'], speed=np.linalg.norm(self.velocity(s),
                                                                                              axis=1).tolist()))
        s.pulse = dict(point=c, time=s.t, energy=released, J=s.L.sum(0) + s.spin @ s.x, x=s.x.copy(),
                       p=np.cross(s.L, s.x))
        s.alive, s.throat, s.spin, s.L = False, [], np.zeros(2), np.zeros((2, 3))
        return s.t + np.pi if s.cycle + 1 < cycles else None

    def _exit(self, s, tf, events, t_exit=None, x_exit=None, frame_exit=None, advanced=False, t_in=None):
        t_exit = s.t if t_exit is None else t_exit
        x_exit = s.x[tf.dest] if x_exit is None else x_exit
        frame_exit = s.frames[tf.dest] if frame_exit is None else frame_exit
        f = Front(from_local(x_exit, frame_exit, tf.local), t_exit, tf.energy, np.zeros(3), tf.phase,
                  f'exit {"AB"[tf.dest]} ({tf.kind})')
        s.fronts.append(f)
        events.append(dict(kind='throat_exit', time=t_exit, t_in=tf.t_in, mouth=tf.dest, cycle=s.cycle,
                           energy=tf.energy, phase=[tf.phase.real, tf.phase.imag], advanced=advanced,
                           origin=f.origin.tolist(), front=tf.kind))

    def _focus(self, s, f, track, events):
        """A front at its antipodal focus: captured by an open mouth there, else reflected."""
        a = self.par.mouth_angle
        focus = -f.origin
        hits = [k for k in range(2) if s.alive and angular_distance(s.x[k], focus) <= a + 1e-12]
        if not hits:
            events.append(dict(kind='reflection', time=s.t, point=focus.tolist(), front=f.kind,
                               phase_in=[f.phase.real, f.phase.imag], foci=f.foci + 1))
            f.origin, f.t = focus, s.t
            f.phase *= self.reflection_phase
            f.foci += 1
            return
        k = min(hits, key=lambda k: angular_distance(s.x[k], focus))
        j = 1 - k
        s.fronts = [g for g in s.fronts if g is not f]
        energy, local = self.capture(s, f, k, focus)
        tf = ThroatFront(k, j, s.tau[k] + self.length, local*self.map, energy, f.phase, f.kind, s.t)
        event = dict(kind='capture', time=s.t, mouth=k, cycle=s.cycle, front=f.kind, energy=f.energy,
                     into_throat=energy, offset=float(s.tau[0] - s.tau[1]), J=f.J.tolist(),
                     phase=[f.phase.real, f.phase.imag], local=local.tolist())
        events.append(event)
        if s.tau[j] < tf.sigma_out:
            s.throat.append(tf)
            event['exit'] = 'scheduled'
            return
        # Advanced edge: the destination's clock passed sigma_out already. Nothing on the way
        # interacts, so the exit front is consistent while its first focus is still ahead.
        times = np.asarray(track['t'])
        t_exit = float(np.interp(tf.sigma_out, np.asarray(track['tau'])[:, j], times))
        i = int(np.clip(np.searchsorted(times, t_exit) - 1, 0, len(times) - 2))
        w = (t_exit - times[i])/(times[i+1] - times[i])
        x = (1-w)*track['x'][i][j] + w*track['x'][i+1][j]
        x /= np.linalg.norm(x)
        fr = (1-w)*track['frames'][i][j] + w*track['frames'][i+1][j]
        e1 = fr[0] - (fr[0] @ x)*x
        e1 /= np.linalg.norm(e1)
        if t_exit + np.pi > s.t:
            event['exit'] = 'advanced'
            self._exit(s, tf, events, t_exit, x, np.array([e1, np.cross(x, e1)]), advanced=True)
        else:
            event['exit'] = 'causality_boundary'
            s.unresolved += energy

    def _rebirth(self, s, events):
        pulse = s.pulse
        x = -pulse['x']
        p, frames, flash = self.rebirth_momenta(s.t, x, -pulse['point'], pulse['energy'], pulse['J'], pulse['p'])
        s.x, s.L, s.frames, s.spin = x.copy(), np.cross(x, p), frames, np.zeros(2)
        s.tau, s.alive, s.birth, s.cycle, s.pulse = np.zeros(2), True, s.t, s.cycle + 1, None
        s.fronts.extend(flash)
        events.append(dict(kind='birth', time=s.t, cycle=s.cycle, point=(-pulse['point']).tolist(),
                           speed=np.linalg.norm(self.velocity(s), axis=1).tolist(), flash=[f.energy for f in flash]))

    def _root(self, prev, h, f):
        """Earliest step length in (0, h] at which f(state) crosses zero (secant)."""
        lo, flo = 0., f(prev)
        hi = h
        q = prev.copy()
        self.step(q, hi)
        fhi = f(q)
        for _ in range(40):
            if abs(fhi) < 1e-13 or hi == lo:
                break
            c = hi - fhi*(hi - lo)/(fhi - flo)
            c = min(max(c, 1e-9*h), h)
            lo, flo = hi, fhi
            hi = c
            q = prev.copy()
            self.step(q, hi)
            fhi = f(q)
        return hi
