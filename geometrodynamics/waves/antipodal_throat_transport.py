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

    def __post_init__(self):
        self.t, self.q, self.s = [np.asarray(x, dtype=float) for x in (self.t, self.q, self.s)]
        if (self.t.ndim != 1 or len(self.t) < 2 or self.q.shape != self.t.shape
                or self.s.shape != self.t.shape or not np.all(np.isfinite([self.t, self.q, self.s]))
                or np.any(np.diff(self.t) <= 0) or np.any(np.diff(self.s) <= 0)
                or not 0 < self.radius < np.pi/2 or self.orientation not in (-1, 1)
                or self.chart_side not in (-1, 1) or not np.isfinite(self.anchor)):
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
                    anchor=self.anchor, opened=self.opened, chart_side=self.chart_side)


def inertial_mouth(q0, speed=0., horizon=20., **kwargs):
    """Explicit prescribed-motion control; no radiation force is inferred."""
    if not np.isfinite(speed) or abs(speed) >= 1 or horizon <= 0:
        raise ValueError('Invalid inertial worldline')
    t = np.linspace(0., horizon, 201)
    return MouthHistory(t, q0+speed*t, t*np.sqrt(1-speed*speed), **kwargs)


class AntipodalTransport:
    def __init__(self, mouths, length=.565, spatial_parity=-1, reflection_phase=-1):
        self.mouths = tuple(mouths)
        if (len(self.mouths) != 2 or not np.isfinite(length) or length <= 0
                or spatial_parity not in (-1, 1) or reflection_phase not in (-1, 1)):
            raise ValueError('Invalid transport geometry')
        self.length = float(length)
        self.spatial_parity = spatial_parity
        # This chosen caustic phase is boundary data, not derived from the old
        # S2 scalar spectrum. It must never be confused with the rim reflection.
        self.reflection_phase = reflection_phase

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
              mouth=0, local=(0., 0.)):
        """Trace one characteristic front, including backward exterior-time edges.

        A finite history, a route cycle, and a leg budget are explicit unresolved
        boundaries. None is absorption. No history is extrapolated, and no seed
        is repeatedly re-injected at a caustic or a throat exit.
        """
        if not np.isfinite(amplitude) or not np.isfinite(emission) or max_legs < 1:
            raise ValueError('Invalid front or route budget')
        if channel not in ('surface', 'throat') or mouth not in (0, 1):
            raise ValueError('Invalid initial channel')
        t = float(emission); origin = unit(origin); local = np.asarray(local, dtype=float)
        if local.shape != (2,) or not np.all(np.isfinite(local)):
            raise ValueError('Invalid aperture coordinate')
        if channel == 'throat' and np.linalg.norm(local) > self.mouths[mouth].radius:
            raise ValueError('Initial throat front outside aperture')
        for m in self.mouths:
            m.check(t)
        sigma = self.mouths[mouth].sigma(t) if channel == 'throat' else None
        edges = []; seen = {}; status = 'leg_budget'; continuation = None
        seed = dict(t=t, origin=origin.tolist(), amplitude=float(amplitude), channel=channel)
        for _ in range(max_legs):
            # Approximate geometric recurrence is reported as unresolved. It is
            # not treated as a solved self-consistent amplitude on a closed loop.
            location = origin if channel == 'surface' else np.r_[mouth, local]
            clock = t if channel == 'surface' else sigma
            key = (channel, *np.round(location, 10), round(float(clock), 10))
            if key in seen:
                status = 'cycle_requires_boundary_solution'
                continuation = dict(edge=seen[key], amplitude=float(amplitude))
                break
            seen[key] = len(edges)
            if channel == 'surface':
                target = -origin; arrival = t+np.pi
                edge = dict(channel='surface', t0=t, t1=arrival, origin=origin.tolist(),
                            target=target.tolist(), amplitude=float(amplitude), distance=float(np.pi))
                try:
                    hit, rim = self.mouth_at_focus(arrival, target)
                except HistoryBoundary:
                    edge.update(event='history_boundary', mouth=None)
                    edges.append(edge); status = 'history_boundary'; break
                if hit is None or not self.mouths[hit].opened:
                    edge.update(event='antipodal_reflection', mouth=hit,
                                reason='empty_focus' if hit is None else 'closed_mouth',
                                outgoing_amplitude=float(amplitude*self.reflection_phase))
                    amplitude *= self.reflection_phase
                    origin = target; t = arrival
                else:
                    edge.update(event='mouth_entry', mouth=hit, rim=rim.tolist())
                    mouth, local = hit, rim
                    t = arrival; sigma = self.mouths[mouth].sigma(t); channel = 'throat'
                edges.append(edge)
            else:
                dest = 1-mouth; outgoing_sigma = sigma+self.length
                # The spatial map acts on both aperture coordinates. A scalar
                # sign flip alone would not implement orientation reversal.
                mapped = local*np.array([1., self.spatial_parity])
                edge = dict(channel='throat', source=mouth, destination=dest,
                            sigma0=float(sigma), sigma1=float(outgoing_sigma), t0=t,
                            t1=None, amplitude=float(amplitude), rim_in=local.tolist(),
                            rim_out=mapped.tolist(), distance=self.length,
                            spatial_parity=self.spatial_parity,
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
                else:
                    edge.update(event='throat_reflection', outgoing_amplitude=float(amplitude*self.reflection_phase))
                    amplitude *= self.reflection_phase
                    mouth, local, sigma, t = dest, mapped, outgoing_sigma, arrival
                edges.append(edge)
        return dict(seed=seed, edges=edges, status=status, continuation=continuation,
                    length=self.length, spatial_parity=self.spatial_parity,
                    reflection_phase=self.reflection_phase,
                    mouths=[m.export() for m in self.mouths],
                    motion='prescribed; no recoil or momentum transfer solved',
                    interpretation='Characteristic transport under the imposed antipodal-only interaction rule')
