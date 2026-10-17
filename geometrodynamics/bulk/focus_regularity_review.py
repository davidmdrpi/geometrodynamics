"""Review controls: vacuum collapse and a smooth, supported static neck.

These are analytic controls, not a constitutive focus law or event evolution.
The original focus_regularity control and its archive remain unchanged.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.special import expit

from geometrodynamics.bulk.focus_regularity import RegularFocus


def vacuum_neck_control(radius=1., mouth_angle=.5):
    """Local Tangherlini vacuum development of the original spatial data.

    At the bifurcation waist, r(tau)^2=b^2-tau^2. Its singular endpoint
    has Kruskal (T,X)=(1,0). The open neck shields the entire approach
    from exterior initial data only when its initial mouth has X_m>1.
    Return log(X_m), avoiding overflow for very small mouth angles.
    """
    g = RegularFocus(radius, mouth_angle).geometry
    y = math.hypot(g.half_length, g.neck_radius)/g.neck_radius
    # y-1 = (S/b)^2/(y+1) avoids cancellation near a=pi/2.
    log_x = y + math.log(g.half_length/g.neck_radius) - math.log(y+1)
    return dict(neck_radius=g.neck_radius, half_length=g.half_length,
                waist_proper_time_to_singularity=g.neck_radius,
                log_kruskal_mouth_x=log_x,
                exterior_causally_excluded=log_x > 0,
                initial_waist_kretschmann=72/g.neck_radius**4)


def _step(z):
    """Flat C-infinity step and its first two derivatives w.r.t. z."""
    z = np.asarray(z, dtype=float)
    w = np.where(z >= 1, 1., 0.)
    wp, wpp = np.zeros_like(z), np.zeros_like(z)
    inside = (z > 0) & (z < 1)
    t = z[inside]
    # Below machine resolution the derivatives vanish numerically.
    with np.errstate(over='ignore', divide='ignore'):
        h = -1/t + 1/(1-t)
    active = abs(h) < 700
    v = expit(h)
    a, hp, hpp = t[active], np.zeros_like(t), np.zeros_like(t)
    hp[active] = 1/a**2 + 1/(1-a)**2
    hpp[active] = -2/a**3 + 2/(1-a)**3
    # expit(-h) retains the derivative tail when expit(h) rounds to one.
    prod = expit(h)*expit(-h)
    w[inside] = v
    wp[inside] = prod*hp
    wpp[inside] = prod*(hpp+(1-2*v)*hp*hp)
    return w, wp, wpp


@dataclass(frozen=True)
class SmoothRegularFocus(RegularFocus):
    """Identical C-infinity spatial collars, distinct C-infinity lapses.

    Collar [S-delta,S] lies inside each original mouth; outside it the
    radius is the round continuation R sin(a+(|s|-S)/R). The lapse bump
    is supported in |s|<S-delta and equals exp(epsilon) at the waist.
    delta/S is a disclosed control input, not a derived focus scale.
    Scalar density is generally nonzero in the collar, equally for both
    epsilon values; the entire effective stress is inverse-engineered.
    """
    collar_fraction: float = .1

    def __post_init__(self):
        super().__post_init__()
        if not math.isfinite(self.collar_fraction) or not 0 < self.collar_fraction <= .25:
            raise ValueError('require 0<collar_fraction<=.25')

    def fields(self, s):
        s = np.asarray(s, dtype=float)
        S, b, R = self.geometry.half_length, self.geometry.neck_radius, self.radius
        delta = self.collar_fraction*S
        L = S-delta
        if np.any(~np.isfinite(s)) or np.any(abs(s) > S+delta):
            raise ValueError('point must be in neck or specified exterior collar patch')
        x = abs(s)
        f = np.hypot(x, b)
        fp, fpp = x/f, b*b/f**3
        angle = self.mouth_angle+(x-S)/R
        g, gp, gpp = R*np.sin(angle), np.cos(angle), -np.sin(angle)/R
        w, wp, wpp = _step((x-L)/delta)
        wp, wpp = wp/delta, wpp/delta**2
        r = (1-w)*f+w*g
        rp = np.sign(s)*((1-w)*fp+w*gp+wp*(g-f))
        rpp = (1-w)*fpp+w*gpp+2*wp*(gp-fp)+wpp*(g-f)
        q, qp, qpp = np.zeros_like(s), np.zeros_like(s), np.zeros_like(s)
        inside = x < L
        u = s[inside]/L
        v = 1-u*u
        exponent = 1-1/v
        active = exponent > -700
        qa, qpa, qppa = np.zeros_like(u), np.zeros_like(u), np.zeros_like(u)
        ua, va = u[active], v[active]
        h1 = -2*ua/(L*va**2)
        h2 = -(2/va**2+8*ua**2/va**3)/L**2
        qa[active] = np.exp(exponent[active])
        qpa[active] = qa[active]*h1
        qppa[active] = qa[active]*(h1*h1+h2)
        q[inside], qp[inside], qpp[inside] = qa, qpa, qppa
        return r, rp, rpp, np.exp(self.epsilon*q), self.epsilon*qp, self.epsilon*qpp
