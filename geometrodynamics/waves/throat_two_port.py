"""The matched three-dimensional bridge as a Maxwell two-port, for routing fronts through it.

#349 (3f72d9a) derives the radiative Maxwell field inside the repository's bridge
(`physical_throat.VacuumThroat`: r_s^2 = 1 - b/r with b = sin^3 a, matched in radius and slope to
the round S3 at the rim) for a stated static lapse N = 1 - beta h, h = (1 - (s/S)^2)^3, with S the
half length. Each polarization and angular index l >= 1 obeys

    psi_tt = N d_s(N psi_s) - N^2 l(l + 1)/r^2 psi,

and the bridge is a two-port between the rims, with characteristic impedance omega at each rim.
This module computes that two-port independently of #349's code: it parametrizes the bridge as
r = b cosh^2 x (ds = 2r dx), propagates (u, N u_s) with closed-form cell exponentials, and
multiplies the cells pairwise. It adds what a router needs:
- the transmission and reflection probabilities T(omega) and R(omega) = 1 - T(omega);
- their group (Wigner) delays, d arg t/d omega and d arg r/d omega, in the rims' proper time
  (N = 1 at the rims). A transmitted front is never allowed out before the characteristic time,
  the integral of ds/N, and a reflected one never before it went in.

The port convention, a free line outside each rim, is #349's. It is least reliable at low
frequency, where the centrifugal barrier l(l + 1)/sin^2 chi already stands outside the rim.
Nothing here models the S3 exterior, a sourced emission or backreaction. The bridge is static
only if it is held open: unsupported, it contracts under Einstein's equation (#349, dc9c689).
"""
import numpy as np

from .physical_throat import VacuumThroat


class ThroatTwoPort:
    def __init__(self, mouth_angle=.3, lapse_depth=0., ell=1, cells=2**14):
        if not 0 < mouth_angle < np.pi/4 or not 0 <= lapse_depth < 1 or int(ell) != ell or ell < 1:
            raise ValueError('Need 0 < mouth_angle < pi/4, 0 <= lapse_depth < 1 and integer ell >= 1')
        self.angle, self.depth, self.ell = float(mouth_angle), float(lapse_depth), int(ell)
        self.b = float(VacuumThroat(mouth_angle).neck_radius())
        rim = np.sin(mouth_angle)
        self.X = float(np.arccosh(np.sqrt(rim/self.b)))
        edges = np.linspace(-self.X, self.X, int(cells) + 1)
        x = (edges[1:] + edges[:-1])/2
        self.dx = edges[1] - edges[0]
        self.half_length = self.b*(self.X + np.sinh(2*self.X)/2)
        self.r = self.b*np.cosh(x)**2
        s = self.b*(x + np.sinh(2*x)/2)
        self.lapse = 1 - self.depth*(1 - (s/self.half_length)**2)**3
        # The characteristic (light-cone) time between the rims.
        self.optical_length = float(np.sum(2*self.r*self.dx/self.lapse))

    @property
    def length(self):
        return 2*self.half_length

    def _amplitudes(self, omega):
        """(t, r) for a wave sent in at the left rim; outgoing-only data at the right rim are
        propagated back through the cells, and the left data split into incident and reflected."""
        L = self.ell*(self.ell + 1)
        alpha = 2*self.r*self.dx/self.lapse
        beta = 2*self.r*self.dx*(self.lapse*L/self.r**2 - omega**2/self.lapse)
        k = np.lib.scimath.sqrt(alpha*beta)
        c = np.cosh(k).real
        sk = np.where(np.abs(k) > 1e-12, np.sinh(k)/np.where(k == 0, 1, k), 1.).real
        cells = np.empty((len(k), 2, 2))
        cells[:, 0, 0], cells[:, 0, 1], cells[:, 1, 0], cells[:, 1, 1] = c, alpha*sk, beta*sk, c
        while len(cells) > 1:                       # pairwise product, right cell on the left
            if len(cells) % 2:
                cells = np.concatenate([cells, np.eye(2)[None]])
            cells = cells[1::2] @ cells[0::2]
        u, j = np.linalg.solve(cells[0], np.array([1., 1j*omega]))
        incident, reflected = (u + j/(1j*omega))/2, (u - j/(1j*omega))/2
        return 1/incident, reflected/incident

    def scattering(self, omega):
        """Transmission and reflection amplitudes t, r at frequency omega (> 0)."""
        if not omega > 0:
            raise ValueError('Positive frequency required')
        return self._amplitudes(float(omega))

    def response(self, omega, step=1e-5):
        """T, R and the group delays (transmitted, reflected) at each frequency, the delays held to
        causality: never below the characteristic time (transmitted) or zero (reflected)."""
        omega = np.atleast_1d(np.asarray(omega, float))
        T, tau_t, tau_r = [], [], []
        for w in omega:
            t, r = self._amplitudes(w)
            dw = step*max(w, 1.)
            tp, rp = self._amplitudes(w + dw)
            tm, rm = self._amplitudes(w - dw)
            T.append(abs(t)**2)
            tau_t.append(np.angle(tp/tm)/(2*dw))
            tau_r.append(np.angle(rp/rm)/(2*dw))
        T = np.array(T)
        return dict(omega=omega, T=T, R=1 - T, tau_t=np.maximum(tau_t, self.optical_length),
                    tau_r=np.maximum(tau_r, 0.), raw_tau_t=np.array(tau_t), raw_tau_r=np.array(tau_r))


def split(weights, response):
    """Split a front with spectral energy weights (summing to 1, on response['omega']) into its
    transmitted and reflected parts: their energy fractions, energy-weighted delays and filtered
    spectra."""
    w = np.asarray(weights, float)
    wt, wr = w*response['T'], w*response['R']
    T, R = float(wt.sum()), float(wr.sum())
    return dict(T=T, R=R,
                tau_t=float(wt @ response['tau_t']/T) if T > 0 else np.inf,
                tau_r=float(wr @ response['tau_r']/R) if R > 0 else 0.,
                transmitted=wt/T if T > 0 else w*0, reflected=wr/R if R > 0 else w*0)
