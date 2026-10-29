"""Antipodal-only characteristic transport on supplied mouth histories.

The interaction support is a model rule: only a surface leg's antipodal focus
can hit a mouth. Intermediate wavefront/rim crossings have NO interaction.
This is a routing boundary problem, not the local scalar-port action in
legacy_history_throat.py, nor a solved recoil/GR field model. Fronts have arbitrary
continuous amplitudes. A route edge is not a quantum or a momentum impulse.
"""
from dataclasses import dataclass
import numpy as np
from scipy.interpolate import PchipInterpolator


class HistoryBoundary(ValueError):
    """A required event is outside the supplied worldtube history."""


def unit(x):
    x = np.asarray(x, dtype=float)
    if x.shape != (3,) or not np.all(np.isfinite(x)) or np.linalg.norm(x) == 0:
        raise ValueError('Expected a finite nonzero 3-vector')
    return x/np.linalg.norm(x)


@dataclass
class MouthHistory:
    """Prescribed meridional worldtube, with strictly increasing proper clock.

    q is unwrapped. sigma = orientation * s + anchor; orientation=-1 reverses
    the time map. Spatial parity is a separate two-dimensional rim map.
    chart_side selects outer (+1) or inner (-1) tangent chart convention; it
    does not change the worldline or the temporal orientation.
    """
    t: np.ndarray
    q: np.ndarray
    s: np.ndarray
    radius: float = .16
    orientation: int = 1
    anchor: float = 0.
    opened: bool = True
    chart_side: int = 1
    end_condition: str = 'dirichlet'

    def __post_init__(self):
        self.t, self.q, self.s = [np.asarray(x, dtype=float) for x in (self.t, self.q, self.s)]
        if (self.t.ndim != 1 or len(self.t) < 2 or self.q.shape != self.t.shape
                or self.s.shape != self.t.shape or not np.all(np.isfinite([self.t, self.q, self.s]))
                or np.any(np.diff(self.t) <= 0) or np.any(np.diff(self.s) <= 0)
                or not 0 < self.radius < np.pi/2 or self.orientation not in (-1, 1)
                or self.chart_side not in (-1, 1) or not np.isfinite(self.anchor)
                or self.end_condition not in ('dirichlet', 'neumann')):
            raise ValueError('Invalid mouth history')
        self._q = PchipInterpolator(self.t, self.q, extrapolate=False)
        self._s = PchipInterpolator(self.t, self.s, extrapolate=False)
        self._inverse = PchipInterpolator(self.s, self.t, extrapolate=False)

    def check(self, t):
        if not np.isfinite(t) or not self.t[0] <= t <= self.t[-1]:
            raise HistoryBoundary('Event outside supplied mouth history')

    def frame(self, t):
        self.check(t); q = float(self._q(t))
        return np.array([np.sin(q), 0., np.cos(q)]), np.array([
            [self.chart_side*np.cos(q), 0., -self.chart_side*np.sin(q)], [0., 1., 0.]])

    def sigma(self, t):
        self.check(t)
        return self.orientation*float(self._s(t))+self.anchor

    def external(self, sigma):
        s = (sigma-self.anchor)/self.orientation
        # Only tolerate floating point noise at an existing boundary.
        tol = 32*np.finfo(float).eps*max(1., abs(s), abs(self.s[-1]))
        if not self.s[0]-tol <= s <= self.s[-1]+tol:
            raise HistoryBoundary('Throat arrival outside supplied proper-clock history')
        return float(self._inverse(np.clip(s, self.s[0], self.s[-1])))

    def local(self, t, point):
        c, frame = self.frame(t); point = unit(point)
        cosine = np.clip(c@point, -1., 1.)
        tangent = point-cosine*c; norm = np.linalg.norm(tangent)
        angle = float(np.arctan2(norm, cosine))
        if norm < 1e-14:
            return np.zeros(2) if cosine > 0 else np.array([np.pi, 0.])
        return angle*(frame@tangent)/norm

    def point(self, t, local):
        c, frame = self.frame(t); local = np.asarray(local, dtype=float)
        angle = np.linalg.norm(local)
        if angle == 0:
            return c
        return unit(np.cos(angle)*c+np.sin(angle)*(local@frame)/angle)

    def export(self):
        return dict(t=self.t.tolist(), q=self.q.tolist(), s=self.s.tolist(),
                    radius=self.radius, orientation=self.orientation,
                    anchor=self.anchor, opened=self.opened, chart_side=self.chart_side,
                    end_condition=self.end_condition)


def inertial_mouth(q0, speed=0., horizon=20., **kwargs):
    """Explicit prescribed-motion control; no radiation force is inferred."""
    if not np.isfinite(speed) or abs(speed) >= 1 or horizon <= 0:
        raise ValueError('Invalid inertial worldline')
    t = np.linspace(0., horizon, 201)
    return MouthHistory(t, q0+speed*t, t*np.sqrt(1-speed*speed), **kwargs)


def phasor(z):
    """JSON representation of an arbitrary complex amplitude: [real, imag]."""
    return [float(np.real(z)), float(np.imag(z))]


class AntipodalTransport:
    def __init__(self, mouths, length=.565, handle='klein',
                 time_policy='transport_orientation', recurrence_tolerance=1e-9,
                 capture_phase='after_caustic'):
        self.mouths = tuple(mouths)
        if (len(self.mouths) != 2 or not np.isfinite(length) or length <= 0
                or handle not in ('klein', 'torus')
                or time_policy not in ('transport_orientation', 'channel_forward')
                or capture_phase not in ('after_caustic', 'before_caustic')
                or not np.isfinite(recurrence_tolerance) or recurrence_tolerance <= 0):
            raise ValueError('Invalid transport geometry or policy')
        self.length = float(length)
        self.handle = handle
        self.time_policy = time_policy
        self.capture_phase = capture_phase
        self.recurrence_tolerance = float(recurrence_tolerance)

    def rim_map(self, source, destination):
        """Derive the chart map from topology, not from a chart-dependent parity.

        In sphere-oriented charts the existing mesh has Klein map -I, and
        torus map diag(1,-1). Changing chart_side changes its matrix, not the
        orientability. Inside |local|<radius this is an ASSUMED DISC PORTAL
        extension: the glued surface itself contains only the rim circle.
        """
        cs = np.diag([self.mouths[source].chart_side, 1.])
        cd = np.diag([self.mouths[destination].chart_side, 1.])
        geometric = -np.eye(2) if self.handle == 'klein' else np.diag([1., -1.])
        return cd@geometric@cs

    @staticmethod
    def boundary_phase(mouth):
        return -1. if mouth.end_condition == 'dirichlet' else 1.

    def mouth_at_focus(self, t, focus):
        hits = []
        for i, mouth in enumerate(self.mouths):
            local = mouth.local(t, focus)
            if np.linalg.norm(local) <= mouth.radius+1e-12:
                hits.append((i, local))
        if len(hits) > 1:
            raise ValueError('Overlapping mouths at antipodal focus; no unique route')
        return hits[0] if hits else (None, None)

    def trace(self, emission, origin, amplitude=1., max_legs=24, channel='surface',
              mouth=0, local=(0., 0.), time_direction=1):
        """Trace a complex front with explicit orientation and recurrence checks.

        transport_orientation carries the time arrow through both attachments:
        dσ/dλ has sign ε_in * dt/dλ; dt_out/dλ has sign ε_out * dσ/dλ.
        channel_forward instead resets the arrow to + in every new channel.
        The latter is separate boundary data, not orientation transport.

        A finite history, a route cycle, and a leg budget remain unresolved.
        A geometric recurrence reports its homogeneous scalar gain condition;
        it is not a solved driven, broadband field boundary problem.
        """
        if (not np.isfinite(amplitude) or not np.isfinite(emission) or max_legs < 1
                or time_direction not in (-1, 1)):
            raise ValueError('Invalid front or route budget')
        if channel not in ('surface', 'throat') or mouth not in (0, 1):
            raise ValueError('Invalid initial channel')
        amplitude = complex(amplitude); gain = 1.+0j
        t = float(emission); origin = unit(origin); local = np.asarray(local, dtype=float)
        if local.shape != (2,) or not np.all(np.isfinite(local)):
            raise ValueError('Invalid aperture coordinate')
        if channel == 'throat' and np.linalg.norm(local) > self.mouths[mouth].radius:
            raise ValueError('Initial throat front outside aperture')
        for m in self.mouths:
            m.check(t)
        sigma = self.mouths[mouth].sigma(t) if channel == 'throat' else None
        direction = time_direction
        edges = []; seen = []; status = 'leg_budget'; continuation = None; closest = None
        seed = dict(t=t, origin=origin.tolist(), amplitude=phasor(amplitude),
                    channel=channel, time_direction=direction)
        for _ in range(max_legs+1):
            location = origin.copy() if channel == 'surface' else local.copy()
            clock = t if channel == 'surface' else sigma
            candidates = []
            for previous in seen:
                if (previous['channel'] != channel or previous['direction'] != direction
                        or (channel == 'throat' and previous['mouth'] != mouth)):
                    continue
                if channel == 'surface':
                    dist = float(np.arctan2(np.linalg.norm(np.cross(location, previous['location'])),
                                          np.clip(location@previous['location'], -1., 1.)))
                else:
                    dist = float(np.linalg.norm(location-previous['location']))
                dt = float(abs(clock-previous['clock']))
                candidates.append((max(dt, dist), dt, dist, previous))
            if candidates:
                metric, dt, dist, previous = min(candidates, key=lambda x: x[0])
                diagnostic = dict(from_edge=previous['edge'], to_edge=len(edges),
                    channel=channel, time_difference=dt, position_difference=dist,
                    max_difference=metric, tolerance=self.recurrence_tolerance)
                if closest is None or metric < closest['max_difference']:
                    closest = diagnostic
                if metric <= self.recurrence_tolerance:
                    loop_gain = gain/previous['gain']
                    residual = abs(amplitude-previous['amplitude'])
                    status = 'cycle_requires_boundary_solution'
                    continuation = dict(edge=previous['edge'], amplitude=phasor(amplitude),
                        loop_gain=phasor(loop_gain), seed_consistency_residual=float(residual),
                        homogeneous_solution='any_amplitude' if abs(loop_gain-1) < 1e-12 else 'zero_only',
                        closure=diagnostic,
                        interpretation='Gain condition a=G*a on a route closed within stated tolerance; not a driven field solution')
                    break
            if len(edges) >= max_legs:
                break
            seen.append(dict(channel=channel, direction=direction, mouth=mouth,
                location=location, clock=clock, edge=len(edges), gain=gain, amplitude=amplitude))
            if channel == 'surface':
                target = -origin; arrival = t+direction*np.pi
                edge = dict(channel='surface', t0=t, t1=arrival, origin=origin.tolist(),
                            target=target.tolist(), amplitude=phasor(amplitude),
                            distance=float(np.pi), time_direction=direction)
                try:
                    hit, rim = self.mouth_at_focus(arrival, target)
                except HistoryBoundary:
                    edge.update(event='history_boundary', mouth=None)
                    edges.append(edge); status = 'history_boundary'; break
                # Exact positive-frequency conformal S2 factor:
                # exp[-i(l+1/2)π] = -i (-1)^l. The antipodal map supplies
                # (-1)^l. Backward propagation is its inverse. Empty foci must
                # carry this phase. An open aperture replaces the caustic with
                # a boundary, whose capture-side phase is an explicit policy.
                before = hit is not None and self.mouths[hit].opened and self.capture_phase == 'before_caustic'
                factor = 1.+0j if before else -1j*direction
                amplitude *= factor; gain *= factor
                edge.update(caustic_factor=phasor(factor), focus_amplitude=phasor(amplitude))
                if hit is None or not self.mouths[hit].opened:
                    edge.update(event='antipodal_reflection', mouth=hit,
                                reason='empty_focus' if hit is None else 'closed_mouth')
                    if hit is not None:
                        factor = self.boundary_phase(self.mouths[hit])
                        amplitude *= factor; gain *= factor
                        edge.update(end_condition=self.mouths[hit].end_condition, boundary_factor=phasor(factor))
                    origin = target; t = arrival
                else:
                    edge.update(event='mouth_entry', mouth=hit, rim=rim.tolist())
                    mouth, local = hit, rim
                    t = arrival; sigma = self.mouths[mouth].sigma(t); channel = 'throat'
                    direction = (self.mouths[mouth].orientation*direction
                                 if self.time_policy == 'transport_orientation' else 1)
                edge['outgoing_amplitude'] = phasor(amplitude)
                edges.append(edge)
            else:
                dest = 1-mouth; outgoing_sigma = sigma+direction*self.length
                matrix = self.rim_map(mouth, dest); mapped = matrix@local
                edge = dict(channel='throat', source=mouth, destination=dest,
                            sigma0=float(sigma), sigma1=float(outgoing_sigma), t0=t,
                            t1=None, amplitude=phasor(amplitude), rim_in=local.tolist(),
                            rim_out=mapped.tolist(), distance=self.length, rim_map=matrix.tolist(),
                            time_direction=direction,
                            time_orientation=self.mouths[mouth].orientation*self.mouths[dest].orientation)
                try:
                    arrival = self.mouths[dest].external(outgoing_sigma)
                except HistoryBoundary:
                    edge.update(event='history_boundary'); edges.append(edge)
                    status = 'history_boundary'; break
                edge['t1'] = arrival
                if np.linalg.norm(mapped) > self.mouths[dest].radius+1e-12:
                    raise ValueError('Mapped front outside destination aperture')
                if self.mouths[dest].opened:
                    origin = self.mouths[dest].point(arrival, mapped)
                    edge.update(event='mouth_exit', target=origin.tolist())
                    t = arrival; channel = 'surface'
                    direction = (self.mouths[dest].orientation*direction
                                 if self.time_policy == 'transport_orientation' else 1)
                else:
                    factor = self.boundary_phase(self.mouths[dest])
                    amplitude *= factor; gain *= factor
                    edge.update(event='throat_reflection', end_condition=self.mouths[dest].end_condition,
                                boundary_factor=phasor(factor))
                    mouth, local, sigma, t = dest, mapped, outgoing_sigma, arrival
                edge['outgoing_amplitude'] = phasor(amplitude)
                edges.append(edge)
        return dict(seed=seed, edges=edges, status=status, continuation=continuation,
                    closest_return=closest, length=self.length, handle=self.handle,
                    time_policy=self.time_policy, capture_phase=self.capture_phase,
                    recurrence_tolerance=self.recurrence_tolerance,
                    chart_map_a_to_b=self.rim_map(0, 1).tolist(),
                    aperture_model='disc_portal_extension; mesh validates rim only',
                    caustic_convention='positive-frequency conformal S2: -i forward, +i backward per empty focus; open capture follows capture_phase',
                    mouths=[m.export() for m in self.mouths],
                    motion='prescribed; no recoil or momentum transfer solved',
                    interpretation='Characteristic transport under the imposed antipodal-only interaction rule')
