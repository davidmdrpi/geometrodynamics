"""Option C: prescribed smooth shell density; no topology-switching geometry.

R is a signed shell-center parameter once the packet carries it through zero.
psi, not R or any one isocontour, is the regular field in this analogy.
"""
from dataclasses import dataclass
import numpy as np

H, END = .6, 9.
X = np.linspace(-H, H, 801)
LEVEL = .5

@dataclass(frozen=True)
class Parameters:
    waist: float = .16
    amplitude: float = .18
    thickness: float = .025

class ContactNeck:
    def __init__(self, parameters=Parameters()):
        self.par = parameters
        if not (0 < parameters.thickness < parameters.waist < 1 and 0 <= parameters.amplitude <= .4):
            raise ValueError('Require 0 < thickness < waist < 1 and amplitude in [0,.4]')
        self.level_radius = parameters.thickness*np.sqrt(-np.log(1-np.sqrt(1-LEVEL)))
        self.band_events = self.crossings(parameters.thickness)
        self.field_events = self.crossings(self.level_radius)
        self.tc = self.band_events[0] if self.band_events else None
        self.xc = self.minimum(self.tc)[1] if self.tc is not None else None

    def rest(self, x):
        return self.par.waist*np.cosh(np.arccosh(1/self.par.waist)*np.asarray(x)/H)

    def wave(self, x, t):
        x = np.asarray(x); z = x-(-.85+.20*t)
        taper = np.maximum(0, 1-(x/H)**2)**2
        return -self.par.amplitude*np.exp(-(z/.18)**2)*np.cos(14*z)*taper

    def radius(self, x, t):
        """Signed center parameter, not an embedded surface when negative."""
        return self.rest(x)+self.wave(x, t)

    def gap(self, t):
        """Legacy radial-band diagnostic, not the field's topology criterion."""
        return 2*(float(np.min(self.radius(X, t)))-self.par.thickness)

    def minimum(self, t):
        """Grid bracket plus continuous golden-section refinement of the minimum."""
        values = self.radius(X, t); i = int(np.argmin(values))
        if i in (0, len(X)-1):
            return float(values[i]), float(X[i])
        a, b = float(X[i-1]), float(X[i+1]); ratio = (np.sqrt(5)-1)/2
        c, d = b-ratio*(b-a), a+ratio*(b-a)
        fc, fd = float(self.radius(c,t)), float(self.radius(d,t))
        for _ in range(28):
            if fc < fd:
                b, d, fd = d, c, fc; c = b-ratio*(b-a); fc = float(self.radius(c,t))
            else:
                a, c, fc = c, d, fd; d = a+ratio*(b-a); fd = float(self.radius(d,t))
        x = (a+b)/2
        return float(self.radius(x,t)), x

    def crossings(self, level):
        # All sign-changing crossings in the horizon. Equality at the analytic
        # threshold is a grazing event at 4.25, not a nonzero closure interval.
        if self.par.amplitude < self.par.waist-level:
            return []
        if abs(self.par.amplitude-(self.par.waist-level)) < 1e-14:
            return [4.25]
        f = lambda t: self.minimum(t)[0]-level
        result = []; lo = 0.; flo = f(lo)
        for hi in np.linspace(0, END, 361)[1:]:
            fhi = f(hi)
            if (flo > 0) != (fhi > 0):
                a, b, fa = lo, float(hi), flo
                for _ in range(36):
                    mid = (a+b)/2; fm = f(mid)
                    if (fa > 0) == (fm > 0): a, fa = mid, fm
                    else: b = mid
                result.append((a+b)/2)
            lo, flo = float(hi), fhi
        return result

    def density(self, x, rho, t):
        r = self.radius(x, t); e = self.par.thickness
        a = np.exp(-((np.asarray(rho)-r)/e)**2)
        b = np.exp(-((np.asarray(rho)+r)/e)**2)
        return a+b-a*b

    def radius_derivatives(self, x, t):
        x = np.asarray(x); z = x-(-.85+.2*t); w = .18; k = 14
        rate = np.arccosh(1/self.par.waist)/H
        f = np.exp(-(z/w)**2)*np.cos(k*z)
        fp = np.exp(-(z/w)**2)*(-2*z/w**2*np.cos(k*z)-k*np.sin(k*z))
        fpp = np.exp(-(z/w)**2)*((4*z*z/w**4-2/w**2-k*k)*np.cos(k*z)+4*k*z/w**2*np.sin(k*z))
        taper = (1-(x/H)**2)**2
        tp = -4*x/H**2*(1-(x/H)**2); tpp = -4/H**2+12*x*x/H**4
        rx = self.par.waist*rate*np.sinh(rate*x)-self.par.amplitude*(fp*taper+f*tp)
        rxx = rate**2*self.rest(x)-self.par.amplitude*(fpp*taper+2*fp*tp+f*tpp)
        return rx, rxx

    def hessian_components(self, x, rho, t):
        """Cartesian Hessian in axial/radial/tangential basis, including axis limit.

        Valid on the displayed neck |x|<=H. Returns xx, rr, tangent, xr.
        """
        r = self.radius(x, t); rx, rxx = self.radius_derivatives(x, t)
        rho = np.asarray(rho); e2 = self.par.thickness**2
        a = np.exp(-(rho-r)**2/e2); b = np.exp(-(rho+r)**2/e2); ab = a*b
        aa = (4*(rho-r)**2/e2**2-2/e2)*a
        bb = (4*(rho+r)**2/e2**2-2/e2)*b
        dr = -2*(rho-r)/e2*a-2*(rho+r)/e2*b+4*rho/e2*ab
        dR = 2*(rho-r)/e2*a-2*(rho+r)/e2*b+4*r/e2*ab
        drr = aa+bb-(16*rho*rho/e2**2-4/e2)*ab
        dRR = aa+bb-(16*r*r/e2**2-4/e2)*ab
        drR = -aa+bb-16*rho*r/e2**2*ab
        tangent = np.divide(dr, rho, out=np.array(drr, copy=True), where=np.abs(rho)>1e-10)
        return dRR*rx*rx+dR*rxx, drr, tangent, drR*rx

    def axis_peak(self, t):
        # If R crosses zero, continuity and positive endpoint R imply psi=1
        # somewhere on the axis, even if the sampling grid misses that root.
        rmin = self.minimum(t)[0]
        g = np.exp(-(max(rmin, 0)/self.par.thickness)**2)
        return 2*g-g*g

    def phase(self, t):
        if self.axis_peak(t) >= LEVEL:
            return 'Diffuse density bridges the axis', 'The field stays smooth; the white contour is a derived boundary.'
        if self.field_events and t > self.field_events[-1]:
            return 'The wave passes; the lumen reopens', 'No caps, permanent separation, or retraction rule is imposed.'
        return 'Opposing diffuse walls approach', 'The resting neck is fixed; the prescribed packet moves the density.'

    def record(self):
        e = self.par.thickness; b = self.par.waist
        samples = [(0., 0., 4.25), (-.129, .03, 3.53), (.1, .05, 4.), (0., -.04, 5.2)]
        return dict(status='KINEMATIC_NOT_DYNAMICAL', regularization='C_SMOOTH_FIELD_DERIVED_SURFACE',
                    parameters=self.par.__dict__, contour_level=LEVEL, contour_radius_threshold=self.level_radius,
                    contact_time=self.tc, contact_position=self.xc, band_crossings=self.band_events,
                    field_contour_crossings=self.field_events,
                    band_amplitude_threshold=b-e, field_amplitude_threshold=b-self.level_radius,
                    threshold_rule='A > threshold gives an interval; equality is grazing at t=4.25',
                    minimum_signed_R=b-self.par.amplitude,
                    weak_wave_minimum_band_gap=2*(b-.08-e),
                    samples=[dict(x=x,rho=rho,t=t,psi=float(self.density(x,rho,t))) for x,rho,t in samples])
