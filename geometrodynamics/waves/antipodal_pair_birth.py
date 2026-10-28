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
  - **Two-port throat** (throat=ThroatTwoPort(...)). Instead of crossing whole in time L, a
    captured front splits by the throat's Maxwell two-port (#349's operator for the matched
    bridge, `throat_two_port`), weighted by the front's spectrum over the field's frequencies:
    a transmitted part, its spectrum filtered by T(omega), leaves the other mouth after its group
    delay; a reflected part, filtered by R(omega), leaves the capturing mouth from the point it
    entered after its own delay, with its momentum reversed (each ray retraces its path). Delays
    are in the rims' proper time, as L is. The two-port's scattering phases are not used; the
    front keeps its phase, as a whole crossing does.
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
- **Signed exchange** (exchange='signed'; the rules above are exchange='rest'). It tests the
  Klein throat as a reversal:
  - The flash is emitted by the moving mouths, isotropically in each mouth's rest frame, so each
    half carries momentum E v. Its amplitude has the sign of its mouth's charge.
  - At capture the mouth takes c J, with c = sign(q) x sign(front) x Re(phase): a like-signed
    front pushes, an opposite-signed one pulls, and one in quadrature does neither. The throat
    holds the rest, (1 - c) J. A pull is capped where that front would carry more angular
    momentum than energy.
  - At exit the front leaves with its momentum mapped through the tube (`tube_map`). The Klein
    handle flips its sign; the torus does not. The exit mouth recoils by the angular momentum
    the map turned. The exit is capped in the same way.

Everything is classical. The amplitude phases are bookkeeping for interference, which is not
modelled. The front-capture rules are stated rules, not derived from a field solution.
"""
from dataclasses import dataclass, field

import numpy as np
from scipy import optimize

from .pair_birth import PairBirth, Parameters, birth_frames, birth_positions, rotate
from .throat_two_port import split


@dataclass(eq=False)
class Front:
    origin: np.ndarray       # where it was launched or last focused
    t: float                 # launch time
    energy: float
    J: np.ndarray            # angular momentum about the centre
    phase: complex = 1.
    kind: str = ''
    foci: int = 0            # caustics passed
    sign: float = 1.         # orientation-sensitive amplitude: compression +1, rarefaction -1
    spectrum: np.ndarray = None   # energy weights over the field's frequencies (two-port throat only)


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
    J: np.ndarray = None     # angular momentum the throat holds for it
    momentum: np.ndarray = None   # (a, b, s): momentum in the capturing mouth's frame and spin about its normal
    sign: float = 1.
    route: str = 'transit'   # 'transit' to the other mouth, or 'reflect' back out of the capturing one
    spectrum: np.ndarray = None


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
    unresolved_J: np.ndarray = field(default_factory=lambda: np.zeros(3))

    def copy(self):
        return AntipodalState(self.t, self.x.copy(), self.L.copy(), self.tau.copy(), self.frames.copy(),
                              self.spin.copy(), self.alive, self.birth, self.cycle,
                              [Front(f.origin.copy(), f.t, f.energy, f.J.copy(), f.phase, f.kind, f.foci, f.sign,
                                     f.spectrum) for f in self.fronts],
                              [ThroatFront(f.source, f.dest, f.sigma_out, f.local.copy(), f.energy, f.phase, f.kind,
                                           f.t_in, None if f.J is None else f.J.copy(),
                                           None if f.momentum is None else f.momentum.copy(), f.sign, f.route,
                                           f.spectrum)
                               for f in self.throat],
                              None if self.pulse is None else dict(self.pulse), self.unresolved,
                              self.unresolved_J.copy())


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
    def __init__(self, parameters=Parameters(), static='none', handle='klein', reflection_phase=-1j,
                 exchange='rest', flash_scale=1., neck_length=None, throat=None):
        if static not in ('none', 'coulomb') or handle not in ('klein', 'torus') or exchange not in ('rest', 'signed') \
                or flash_scale <= 0:
            raise ValueError("static 'none' or 'coulomb', handle 'klein' or 'torus', exchange 'rest' or 'signed', "
                             "flash_scale > 0")
        self.par = parameters
        self.exchange = exchange
        self.flash_scale = flash_scale
        self.field = PairBirth(parameters)       # flash, Coulomb energy and charges
        self.static = static
        self.handle = handle
        self.map = np.array([-1., -1.]) if handle == 'klein' else np.array([-1., 1.])
        self.reflection_phase = complex(reflection_phase)
        a = parameters.mouth_angle
        self.length = np.sin(2*a) if neck_length is None else float(neck_length)   # S3: s3_focusing
        self.inertia = parameters.mass*np.sin(a)**2
        self.throat = throat
        if throat is not None:
            if abs(throat.angle - a) > 1e-12:
                raise ValueError("The two-port's mouth angle must be the model's")
            if neck_length is None:
                self.length = throat.length
            # The field's distinct frequencies, sqrt(l(l + 1)), and each mode's place among them.
            self.spectral_omega, self._mode_frequency = np.unique(self.field.omega, return_inverse=True)
            self.response = throat.response(self.spectral_omega)

    # ---- births and fronts ----------------------------------------------------------------------
    def flash_fronts(self, x, vel, t):
        """#345's birth flash, the free wave of -Q_qs, split into one front per mouth: energy (self
        part plus half the cross term, times flash_scale) and angular momentum.
        - exchange='rest': launched at rest, as the field is, so it carries no angular momentum.
        - exchange='signed': each half is emitted by its moving mouth, isotropically in the mouth's
          rest frame, so it carries momentum E v and angular momentum E x cross v. Its amplitude
          carries the sign of its mouth's charge: compression from +q, rarefaction from -q."""
        f = self.field
        Q = -f.quasi_static(x, vel)
        e = [.5*np.vdot(f.omega*Q[k], f.omega*Q[k]).real for k in range(2)]
        cross = np.vdot(f.omega*Q[0], f.omega*Q[1]).real
        e = [self.flash_scale*(e[k] + cross/2) for k in range(2)]
        if self.exchange == 'signed':
            J = [e[k]*np.cross(x[k], vel[k]) for k in range(2)]
            signs = np.sign(f.charges)
        else:
            zero = np.zeros(f.n, complex)
            J = [self.flash_scale*f.field_angular(Q[k], zero) for k in range(2)]      # zero: it starts at rest
            signs = np.ones(2)
        spectra = [None, None]
        if self.throat is not None:
            for k in range(2):
                w = np.bincount(self._mode_frequency, .5*np.abs(f.omega*Q[k])**2, len(self.spectral_omega))
                spectra[k] = w/w.sum()
        return [Front(x[k].copy(), t, e[k], J[k], 1., f'flash {"AB"[k]}', 0, float(signs[k]), spectra[k])
                for k in range(2)]

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
        J = J + sum((f.J for f in s.throat if f.J is not None), np.zeros(3)) + s.unresolved_J
        if s.pulse is not None:
            J = J + s.pulse['J']
        if s.alive:
            J = J + s.L.sum(0) + s.spin @ s.x
        return J

    def separation(self, s):
        return angular_distance(s.x[0], s.x[1])

    def separation_rate(self, s):
        v = self.velocity(s)
        return float(-(v[0] @ s.x[1] + v[1] @ s.x[0])/max(np.sin(self.separation(s)), 1e-300))

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
    def coupling(self, front, k):
        """Share of the front's angular momentum mouth k takes at capture.
        - 'rest': all of it.
        - 'signed': c = sign(q_k) x sign(front) x Re(Gouy phase). A same-sign front pushes the
          mouth along its arrival, an opposite-sign front pulls it back, as #345's own returns
          and partner flash did. A front in quadrature (one empty focus behind it) does neither."""
        if self.exchange == 'rest':
            return 1.
        return float(np.sign(self.field.charges[k])*front.sign*np.real(front.phase))

    def _gain(self, s, k, J):
        """Kinetic and spin energy mouth k would gain by taking angular momentum J."""
        x, M = s.x[k], self.par.mass
        Lk = s.L[k] + J - np.dot(J, x)*x
        spin = s.spin[k] + np.dot(J, x)
        return (np.sqrt(M**2 + Lk @ Lk) - self.energies(s)[k]
                + (spin**2 - s.spin[k]**2)/(2*self.inertia))

    def _transfer(self, s, k, J, energy):
        """Give angular momentum J to mouth k (orbit across its normal, spin along it); returns
        the energy left over from `energy` after the mouth's kinetic and spin energy change."""
        left = energy - self._gain(s, k, J)
        if left < -1e-12:
            raise RuntimeError('A front cannot pay for the mouth it moves')
        x = s.x[k]
        s.L[k] += J - np.dot(J, x)*x
        s.spin[k] += np.dot(J, x)
        return max(left, 0.)

    def capture(self, s, front, k, focus):
        """Mouth k captures `front` at its focus. It takes c J (see `coupling`); the throat takes the
        rest of the angular momentum and of the energy. A pull is limited by the throat: the front
        it holds cannot carry more angular momentum than energy, so if c J would leave it more, c
        is raised to that bound (saturated). Returns (energy into the throat, aperture coordinates,
        angular momentum into the throat, c, saturated)."""
        J, E = front.J, front.energy

        def room(c):
            return E - self._gain(s, k, c*J) - abs(1 - c)*np.linalg.norm(J)

        c = self.coupling(front, k)
        saturated = bool(room(c) < 0)
        if saturated:
            c = optimize.brentq(room, c, 1., xtol=1e-15)
        inside = self._transfer(s, k, c*J, E)
        return inside, local_coordinates(s.x[k], s.frames[k], focus), (1 - c)*J, c, saturated

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
                if event is not None and event[0] == 'contact' and focus_time - prev.t <= event[1] + 1e-9:
                    # A front focusing as the rims touch is caught first, and leaves with the release.
                    s = prev.copy()
                    self.step(s, focus_time - prev.t)
                    event = None
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
        held = sum((f.J for f in s.throat if f.J is not None), np.zeros(3))
        s.pulse = dict(point=c, time=s.t, energy=released, J=s.L.sum(0) + s.spin @ s.x + held, x=s.x.copy(),
                       p=np.cross(s.L, s.x))
        s.alive, s.throat, s.spin, s.L = False, [], np.zeros(2), np.zeros((2, 3))
        return s.t + np.pi if s.cycle + 1 < cycles else None

    def tube_map(self, momentum, sign):
        """How a front leaves the tube, in sphere-oriented rim frames, from its momentum (a, b) in
        the capturing mouth's frame and its spin s about that mouth's outward normal. Following
        its rays through the tube (in radially, out radially at the mapped rim angle):
        - Klein: (a, b, s) unchanged, and the amplitude sign flips, because the surface normal
          reverses along any loop through a Klein handle (compression becomes rarefaction).
        - Torus: (a, b, s) -> (a, -b, -s) and the sign is kept.
        In ambient terms, for antipodal mouths with birth frames (e1 along the birth circle), both
        maps turn J by pi; neither reverses it as a whole. The Klein map reverses the spin and the
        orbital part from motion along e1, and keeps the part from motion along e2. The torus
        reverses both orbital parts and keeps the spin. So both reverse the orbital angular
        momentum of a front moving along the birth circle, and there they differ only in sign."""
        a, b, sp = momentum
        if self.handle == 'klein':
            return np.array([a, b, sp]), -sign
        return np.array([a, -b, -sp]), sign

    def _exit(self, s, tf, events, t_exit=None, x_exit=None, frame_exit=None, advanced=False, t_in=None):
        """A front leaves the throat at mouth j. With exchange='signed' it leaves with the tube-mapped
        momentum and sign, and mouth j recoils by the angular momentum the handle turned (the throat
        held tf.J). For an advanced exit the recoil is applied when the exit is resolved. Like a pull
        at capture, the exit is capped where the front would carry more angular momentum than energy."""
        j = tf.dest
        t_exit = s.t if t_exit is None else t_exit
        x_exit = s.x[j] if x_exit is None else x_exit
        frame_exit = s.frames[j] if frame_exit is None else frame_exit
        energy, J_out, sign, recoil, saturated = tf.energy, np.zeros(3), tf.sign, np.zeros(3), False
        rate = self.separation_rate(s)
        if self.exchange == 'signed' and tf.momentum is not None:
            if tf.route == 'reflect':
                (a, b, sp), sign = -tf.momentum, tf.sign     # each ray retraces its path; no handle crossed
            else:
                (a, b, sp), sign = self.tube_map(tf.momentum, tf.sign)
            pi = a*frame_exit[0] + b*frame_exit[1]
            J_map = np.cross(x_exit, pi) + sp*x_exit

            def room(w):
                return tf.energy - self._gain(s, j, tf.J - w*J_map) - w*np.linalg.norm(J_map)

            # The front leaves with the mapped momentum, scaled down only if, after paying for the
            # recoil, it would carry more angular momentum than energy (saturated).
            saturated = bool(room(1.) < 0)
            if saturated and room(0.) < 0:
                raise RuntimeError('The throat front cannot pay for the exit recoil')
            J_out = (optimize.brentq(room, 0., 1., xtol=1e-15) if saturated else 1.)*J_map
            recoil = tf.J - J_out
            energy = self._transfer(s, j, recoil, tf.energy)
        label = 'reflected' if tf.route == 'reflect' else 'exit'
        f = Front(from_local(x_exit, frame_exit, tf.local), t_exit, energy, J_out, tf.phase,
                  f'{label} {"AB"[j]} ({tf.kind})', 0, sign, tf.spectrum)
        s.fronts.append(f)
        events.append(dict(kind='throat_exit', time=t_exit, t_in=tf.t_in, mouth=j, cycle=s.cycle,
                           energy=energy, phase=[tf.phase.real, tf.phase.imag], advanced=advanced,
                           origin=f.origin.tolist(), front=tf.kind, sign=sign, recoil=recoil.tolist(),
                           saturated=saturated, drate=self.separation_rate(s) - rate, route=tf.route,
                           source=tf.source, J_out=J_out.tolist()))

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
        L_before = s.L[k].copy()
        rate = self.separation_rate(s)
        energy, local, J_in, c, saturated = self.capture(s, f, k, focus)
        x, frame = s.x[k], s.frames[k]
        pi = np.cross(J_in - np.dot(J_in, x)*x, x)
        momentum = np.array([pi @ frame[0], pi @ frame[1], np.dot(J_in, x)])
        event = dict(kind='capture', time=s.t, mouth=k, cycle=s.cycle, front=f.kind, energy=f.energy,
                     into_throat=energy, offset=float(s.tau[0] - s.tau[1]), J=f.J.tolist(),
                     phase=[f.phase.real, f.phase.imag], local=local.tolist(), sign=f.sign, coupling=c,
                     saturated=saturated,
                     kick=np.cross(s.L[k] - L_before, x).tolist(), drate=self.separation_rate(s) - rate)
        events.append(event)
        if self.throat is None:
            tf = ThroatFront(k, j, s.tau[k] + self.length, local*self.map, energy, f.phase, f.kind, s.t, J_in,
                             momentum, f.sign)
            event['exit'] = self._enter(s, tf, track, events)
            return
        if f.spectrum is None:
            raise ValueError('A front routed through the two-port needs a spectrum')
        parts = split(f.spectrum, self.response)
        event['transmitted'] = dict(fraction=parts['T'], delay=parts['tau_t'])
        event['reflected'] = dict(fraction=parts['R'], delay=parts['tau_r'])
        for share, dest, delay, where, route, spectrum in (
                (parts['T'], j, parts['tau_t'], local*self.map, 'transit', parts['transmitted']),
                (parts['R'], k, parts['tau_r'], local, 'reflect', parts['reflected'])):
            if share > 0:
                tf = ThroatFront(k, dest, s.tau[k] + delay, where, share*energy, f.phase, f.kind, s.t, share*J_in,
                                 share*momentum, f.sign, route, spectrum)
                event['transmitted' if route == 'transit' else 'reflected']['exit'] = self._enter(s, tf, track,
                                                                                                    events)
        event['exit'] = event['transmitted'].get('exit')

    def _enter(self, s, tf, track, events):
        """Schedule a throat front's exit; returns how it leaves."""
        j = tf.dest
        if s.tau[j] < tf.sigma_out:
            s.throat.append(tf)
            return 'scheduled'
        if tf.route == 'reflect':
            self._exit(s, tf, events)          # no delay: straight back out
            return 'immediate'
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
            self._exit(s, tf, events, t_exit, x, np.array([e1, np.cross(x, e1)]), advanced=True)
            return 'advanced'
        s.unresolved += tf.energy
        s.unresolved_J = s.unresolved_J + tf.J
        return 'causality_boundary'

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


# ---- antipodal time reversal ---------------------------------------------------------------------
def resonant_speeds(mouth_angle, count=4):
    """Birth speeds u_k = (pi - 2a)/((2k + 1) pi), k = 0, 1, ... A symmetric pair coasting at u_k
    reaches contact, at t_c = (pi - 2a)/u_k = (2k + 1) pi, exactly when each mouth's birth flash
    refocuses on its partner (a flash focuses at the antipode of its birth point at odd multiples
    of pi, and that is where the partner meets contact)."""
    return (np.pi - 2*mouth_angle)/((2*np.arange(count) + 1)*np.pi)


def closure_band(mouth_angle, k, length=None):
    """Birth speeds (u_low, u_k] at which a symmetric strict cycle closes on itself: the partner
    catches each flash at its (2k + 1)-th focus, (2k + 1) pi <= t_c <= (2k + 1) pi + a/u (inside
    the aperture, before contact), and the caught energy is still
    inside the throat when the rims touch, t_c < (2k + 1) pi + L/sqrt(1 - u^2), so it leaves with
    the release pulse. Only the upper end, u_k, makes emission and absorption exact images; there
    focus and contact coincide, and the front counts as caught just before contact."""
    from scipy.optimize import brentq
    a = mouth_angle
    L = np.sin(2*a) if length is None else length
    uk = (np.pi - 2*a)/((2*k + 1)*np.pi)
    throat = brentq(lambda u: (np.pi - 2*a)/u - (2*k + 1)*np.pi - L/np.sqrt(1 - u*u), 1e-6, uk*(1 - 1e-12),
                    xtol=1e-15)
    return max(throat, (np.pi - 3*a)/((2*k + 1)*np.pi)), uk     # and the flash must be caught at all


def two_port_closure_band(model, k):
    """The closure band of `closure_band` for a model with a two-port throat. Nearly all of a
    caught flash is reflected, so it stays in the throat only for the reflection delay tau_r of
    its spectrum, not for the neck's length; the band's lower edge solves
    t_c = (2k + 1) pi + tau_r(u)/sqrt(1 - u^2), with tau_r taken from the birth flash at speed u.
    The transmitted part, delayed at least the characteristic time, is then in the throat too."""
    from scipy.optimize import brentq
    a = model.par.mouth_angle
    uk = (np.pi - 2*a)/((2*k + 1)*np.pi)

    def kick_for(u):
        return brentq(lambda v: np.linalg.norm(model.velocity(model.initial(v)), axis=1)[0] - u, .02, .95,
                      xtol=1e-14)

    def tau_r(u):
        return split(model.initial(kick_for(u)).fronts[0].spectrum, model.response)['tau_r']

    def late(u):
        return (np.pi - 2*a)/u - (2*k + 1)*np.pi - tau_r(u)/np.sqrt(1 - u*u)

    edge = (np.pi - 3*a)/((2*k + 1)*np.pi)        # and the flash must be caught at all
    if late(edge) <= 0:
        return edge, uk
    return brentq(late, edge, uk*(1 - 1e-12), xtol=1e-12), uk


def _slerp_history(t, x, when):
    """Position at time `when` on a recorded strict history. Between events the strict flight is
    a constant-speed geodesic, so slerp between consecutive samples is exact there."""
    i = int(np.clip(np.searchsorted(t, when, side='right') - 1, 0, len(t) - 2))
    while i + 1 < len(t) - 1 and t[i + 1] <= t[i]:
        i += 1
    span = t[i + 1] - t[i]
    w = 0. if span <= 0 else (when - t[i])/span
    angle = angular_distance(x[i], x[i + 1])
    if angle < 1e-15:
        return x[i].copy()
    p = (np.sin((1 - w)*angle)*x[i] + np.sin(w*angle)*x[i + 1])/np.sin(angle)
    return p/np.linalg.norm(p)


def elliptic_mismatch(h, info):
    """How far the first flight of a strict history is from being invariant under the free
    involution (n, t, q) -> (-n, t_c - t, -q), with t_c the first rim contact: the antipodal map,
    a time reflection about mid-flight and charge conjugation, which exchanges the mouths.
    - trajectory: max over the flight of |x_A(t) + x_B(t_c - t)| (the image of A must be B);
    - flash: the image of a mouth emitting its birth flash at (x_k(0), 0) is its partner
      absorbing a wave at (-x_k(0), t_c). The flash is at that point only at odd multiples of pi,
      so the flash mismatch is the distance from t_c to the nearest one.
    Returns (trajectory mismatch, flash mismatch, t_c, captures of a birth flash in the flight)."""
    contact = next(e['time'] for e in info['events'] if e['kind'] == 'contact')
    t, x = np.asarray(h['t']), np.asarray(h['x'])
    keep = t <= contact + 1e-12
    worst = 0.
    for when in t[keep]:
        a = _slerp_history(t, x[:, 0], when)
        b = _slerp_history(t, x[:, 1], contact - when)
        worst = max(worst, float(np.linalg.norm(a + b)))
    n = max(1, 2*int(np.floor(contact/(2*np.pi))) + 1)
    flash = min(abs(contact - k*np.pi) for k in (n - 2, n, n + 2) if k > 0)
    captures = [e for e in info['events'] if e['kind'] == 'capture' and e['front'].startswith('flash')
                and e['time'] <= contact]
    return worst, float(flash), float(contact), captures
