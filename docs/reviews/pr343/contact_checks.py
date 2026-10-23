"""Independent checks of PR #343's contact event: timing, curvature by surface, the jump, thresholds.

    PYTHONPATH=scripts python docs/reviews/pr343/contact_checks.py     # from a checkout of #343
"""
import numpy as np
from scipy import optimize
from neck_wave_contact import ContactNeck, Parameters, H, X, END

m = ContactNeck()
P = m.par
R = lambda x, t: m.radius(x, t)

def Rx(x, t, e=1e-6): return (R(x+e, t)-R(x-e, t))/(2*e)
def Rxx(x, t, e=1e-4): return (R(x+e, t)-2*R(x, t)+R(x-e, t))/e**2

# 1. Continuous minimum: refine min_x R on a fine grid then Brent, then root of gap in t.
def rmin(t):
    xs = np.linspace(-H, H, 20001); i = np.argmin(R(xs, t))
    lo, hi = xs[max(i-2, 0)], xs[min(i+2, len(xs)-1)]
    res = optimize.minimize_scalar(lambda x: R(x, t), bounds=(lo, hi), method='bounded', options={'xatol': 1e-13})
    return res.fun, res.x
tc = optimize.brentq(lambda t: rmin(t)[0]-P.thickness, 3.0, 3.531, xtol=1e-14)
rc, xc = rmin(tc)
print(f'PR tc = {m.tc:.12f}, xc = {m.xc:.6f} (grid)')
print(f'continuous tc = {tc:.12f}, xc = {xc:.9f}, diff = {m.tc-tc:.2e}')

# 2. Contact geometry: radius, slope, axial curvature, hoop curvature of mid/outer/inner surfaces.
k_ax = Rxx(xc, tc)/(1+Rx(xc, tc)**2)**1.5
print(f'at contact: mid radius R = {rc:.6f} (= eps), dR/dx = {Rx(xc,tc):.2e}, axial curvature R_xx = {Rxx(xc,tc):.4f}')
print(f'  mid surface hoop curvature 1/R = {1/rc:.2f}; outer 1/(R+eps) = {1/(rc+P.thickness):.2f}; inner 1/(R-eps) -> infinity')
print(f'  rest-neck waist hoop curvature 1/b = {1/P.waist:.2f}')

# 3. Approach: inner-surface radius and hoop curvature vs time-to-contact.
print('t_c - t     min(R-eps)    1/min(R-eps)   d/dt min R')
for s in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6]:
    r, _ = rmin(tc-s); dr = (rmin(tc-s+1e-7)[0]-rmin(tc-s-1e-7)[0])/2e-7
    print(f'{s:9.0e}   {r-P.thickness:.4e}   {1/(r-P.thickness):12.4e}   {dr:+.4f}')

# 4. Local model of the lumen near contact: R - eps ~ v (tc - t) + (R_xx/2)(x - xc)^2.
#    The inner boundary's zero set at tc is a single point; just after, an interval: Morse index 1 in (x,t).
v = -(rmin(tc+1e-6)[0]-rmin(tc-1e-6)[0])/2e-6
print(f'closing speed v = {v:.5f}; after tc the overlap half-length grows as sqrt(2 v (t-tc)/R_xx): '
      f'at t-tc = 1e-3 -> {np.sqrt(2*v*1e-3/Rxx(xc,tc)):.4f}')
# Interior extent of overlap R < eps after contact (what the PR's model would imply if not replaced).
for s in [1e-3, 1e-2, 5e-2]:
    xs = np.linspace(-H, H, 200001); inside = xs[R(xs, tc+s) < P.thickness]
    print(f'  t - tc = {s:.0e}: R < eps on x in [{inside.min():.4f}, {inside.max():.4f}], min R = {R(xs, tc+s).min():.4f}')

# 5. Discontinuities at reconnection.
pre = m.branches(m.tc-1e-9, 2001)[0]
post = m.branches(m.tc, 2001)
L, Rt = post
print(f'post-contact tips at x = {L[0][-1]:.6f}, {Rt[0][-1]:.6f}: separation {Rt[0][-1]-L[0][-1]:.6f} '
      f'(= 0.018 * 1.2 = {.018*1.2:.4f}) appears at t = tc')
xs = np.linspace(-H, H, 4001)
def interp(part, x):
    xx, rr = part[0], part[1]; o = np.argsort(xx); return np.interp(x, xx[o], rr[o], left=np.nan, right=np.nan)
rl, rr_ = interp(L, xs), interp(Rt, xs)
rpost = np.where(np.isfinite(rl), rl, rr_)
rpre = R(xs, m.tc)
ok = np.isfinite(rpost)
print(f'max |R_post - R_pre| over the surviving x = {np.nanmax(abs(rpost-rpre)[ok]):.4f}; '
      f'removed axial length {np.sum(~ok)*(xs[1]-xs[0]):.4f}')
# Cap-tip curvature from the drawn arrays: near the tip d = a r^2, so the tip curvature is 2a.
for x, r, _, _ in m.branches(m.tc, 20001):
    d = abs(x-x[-1]); k = (d > 0) & (d < 2e-4)
    a = np.polyfit(r[k]**2, d[k], 1)[0]
    print(f'cap tip at x = {x[-1]:+.4f}: curvature {2*a:.1f}, radius of curvature {1/(2*a):.4f} (eps = {P.thickness})')

# 6. Reopening under the prescribed wave alone, and the contact thresholds.
xs = np.linspace(-H, H, 40001)
ts = np.linspace(m.tc+1e-6, END, 9001)
vals = np.array([m.radius(xs, t).min()-P.thickness for t in ts])
print(f'prescribed lumen reopens at t = {ts[np.argmax(vals > 0)]:.3f}')
def bisect(has, lo, hi):
    for _ in range(30):
        mid = (lo+hi)/2
        if has(mid): hi = mid
        else: lo = mid
    return hi
A_star = bisect(lambda A: ContactNeck(Parameters(amplitude=A)).tc is not None, .08, .18)
b_star = bisect(lambda b: ContactNeck(Parameters(waist=b)).tc is None, .16, .30)
print(f'contact thresholds: A* = {A_star:.4f} at b = 0.16 (b - eps = {.16-P.thickness}); '
      f'b* = {b_star:.4f} at A = 0.18 (A + eps = {.18+P.thickness})')
w = ContactNeck(Parameters(amplitude=.08))
print(f'weak wave closest gap: {min(w.gap(t) for t in np.linspace(0, END, 3601)):.4f}')
