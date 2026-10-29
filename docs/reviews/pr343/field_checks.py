"""Independent checks of PR #343's option-C field (cb4b5ec).

    PYTHONPATH=scripts python docs/reviews/pr343/field_checks.py    # from a checkout of #343
"""
import numpy as np
from scipy import optimize
from neck_wave_contact import ContactNeck, Parameters, H, LEVEL

m = ContactNeck(); e = m.par.thickness
xs = np.linspace(-H, H, 24001)
def rmin(mm, t):
    i = np.argmin(mm.radius(xs, t)); lo, hi = xs[max(i-2, 0)], xs[min(i+2, len(xs)-1)]
    return optimize.minimize_scalar(lambda x: mm.radius(x, t), bounds=(lo, hi), method='bounded', options={'xatol': 1e-13}).fun
# 1. Event pairs from scratch: axis density 2g - g^2 = 1/2 with g = exp(-(R/eps)^2).
r_half = e*np.sqrt(-np.log(1-np.sqrt(1-LEVEL)))
for name, lev, rec in [('band R = eps', e, m.band_events), ('contour psi = 1/2', r_half, m.field_events)]:
    f = lambda t: rmin(m, t)-lev
    roots = [optimize.brentq(f, a, b, xtol=1e-14) for a, b in [(3.3, 4.25), (4.25, 5.2)]]
    print(f'{name}: mine {roots[0]:.10f} / {roots[1]:.10f}; PR {rec[0]:.10f} / {rec[1]:.10f}; '
          f'max diff {max(abs(np.array(roots)-rec)):.1e}; sum {roots[0]+roots[1]:.10f} (= 2 x 4.25: time-symmetric)')
# 2. psi is even in rho and in R: the formula cannot tell a wall at +R from one at -R.
G = lambda s: np.exp(-(s/e)**2); psi = lambda rho, R: G(rho-R)+G(rho+R)-G(rho-R)*G(rho+R)
rho, R = np.linspace(-.1, .1, 201)[:, None], np.linspace(-.1, .1, 201)[None, :]
print('psi(rho, R) = psi(-rho, R) = psi(rho, -R):', np.allclose(psi(rho, R), psi(-rho, R)), np.allclose(psi(rho, R), psi(rho, -R)))
# 3. Where the sampled Hessian maximum sits: contact region vs rims.
def hess_norm(xg, rg, t):
    xx, rr, tt, xr = m.hessian_components(xg, rg, t)
    return np.sqrt(xx**2+rr**2+tt**2+2*xr**2)
X, Rr = np.meshgrid(np.linspace(-H, H, 481), np.linspace(0, 1.08, 433))
best = max(((hess_norm(X, Rr, t).max(), t) for t in np.linspace(0, 9, 91)))
hn = hess_norm(X, Rr, best[1]); i = np.unravel_index(np.argmax(hn), hn.shape)
near = (abs(X-m.xc) < .15) & (Rr < .15)
cmax = max(hess_norm(X, Rr, t)[near].max() for t in np.linspace(3.3, 5.2, 39))
print(f'sampled Hessian max x eps^2 = {best[0]*e*e:.2f} at t = {best[1]:.1f}, x = {X[i]:+.3f}, rho = {Rr[i]:.3f} (near a rim); '
      f'max within 0.15 of the contact point over t in [3.3, 5.2]: {cmax*e*e:.2f}')
# 4. Large amplitude: signed R carries the walls through the axis.
for A in [.18, .1877, .19, .25, .30]:
    mm = ContactNeck(Parameters(amplitude=A)); xa = np.linspace(-H, H, 24001)
    psi_axis = mm.density(xa, 0., 4.25); R = mm.radius(xa, 4.25)
    bridged = psi_axis >= LEVEL
    runs = np.flatnonzero(np.diff(bridged.astype(int)))
    print(f'A = {A:.4f}: min R at t = 4.25 is {R.min():+.4f}; axis psi >= 1/2 on {len(runs)//2} separate interval(s); '
          f'axis psi at x = 0: {mm.density(0., 0., 4.25):.3f}; phase: "{mm.phase(4.25)[0]}"')
print(f'walls pass through each other (an open axis gap between two bridges) once A > b + r_half = {.16+r_half:.4f}')

# 5. A smooth positive part keeps the walls merged instead of crossing: R+ = (R + sqrt(R^2 + d^2))/2.
class Clamped(ContactNeck):
    d = e/4
    def radius(self, x, t):
        r = super().radius(x, t); return (r+np.sqrt(r*r+self.d**2))/2
for A in [.18, .30]:
    mm = Clamped(Parameters(amplitude=A)); xa = np.linspace(-H, H, 24001)
    bridged = mm.density(xa, 0., 4.25) >= LEVEL
    f = lambda t: rmin(mm, t)-r_half
    ev = [optimize.brentq(f, a, b, xtol=1e-12) for a, b in [(2.5, 4.25), (4.25, 6.)]]
    print(f'clamped (d = eps/4), A = {A:.2f}: axis psi >= 1/2 on {len(np.flatnonzero(np.diff(bridged.astype(int))))//2} interval(s) at t = 4.25; '
          f'contour closes/reopens at {ev[0]:.4f} / {ev[1]:.4f}')
