"""Review figure for PR #343: what is regular and what is singular at the wave-contact event.

    PYTHONPATH=scripts python docs/reviews/pr343/contact_figure.py docs/reviews/pr343/contact.png   # from a checkout of #343
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import optimize
from neck_wave_contact import ContactNeck, H, END

m = ContactNeck(); E = m.par.thickness
xs = np.linspace(-H, H, 40001)

def rmin(t):
    i = np.argmin(m.radius(xs, t)); lo, hi = xs[max(i-2, 0)], xs[min(i+2, len(xs)-1)]
    r = optimize.minimize_scalar(lambda x: m.radius(x, t), bounds=(lo, hi), method='bounded', options={'xatol': 1e-13})
    return r.fun, r.x
TC = optimize.brentq(lambda t: rmin(t)[0]-E, 3., 3.531, xtol=1e-14); XC = rmin(TC)[1]
def d2(x, t, e=1e-4): return (m.radius(x+e, t)-2*m.radius(x, t)+m.radius(x-e, t))/e**2

fig = plt.figure(figsize=(15, 9.2), constrained_layout=True)
gs = fig.add_gridspec(2, 3)
BAND, MID, INNER, CAP, MORSE = '#b98bc0', '#7a4f80', '#c0392b', '#138d75', '#1f5fa8'
zx = np.linspace(XC-.2, XC+.2, 2001)

def band(ax, t, title, alpha=1., zx=zx):
    r = m.radius(zx, t)
    for s in (-1, 1):
        ax.fill_between(zx, s*np.maximum(r-E, 0), s*(r+E), color=BAND, alpha=.55*alpha, lw=0)
        ax.plot(zx, s*r, color=MID, lw=1.2, alpha=alpha)
        ax.plot(zx, s*np.maximum(r-E, 0), color=INNER, lw=1.6, alpha=alpha)
    ax.axhline(0, color='k', lw=.6, ls=':')
    w = zx[-1]-zx[0]
    ax.set(xlim=(zx[0], zx[-1]), ylim=(-w/2, w/2), title=title, xlabel='x', ylabel='ρ (meridian section)')
    ax.set_aspect('equal')

ax = fig.add_subplot(gs[0, 0])
band(ax, TC-.3, f'(a) t = t$_c$ − 0.3: shell band [R−ε, R+ε], ε = {E}')
ax.text(.02, .04, 'red = facing (inner) surface R − ε\npurple line = mid-surface R', transform=ax.transAxes, fontsize=9, bbox=dict(fc='white', alpha=.85, ec='none'))

ax = fig.add_subplot(gs[0, 1])
band(ax, TC, '(b) t = t$_c$: the PR\'s "facing gap" = lumen diameter = 0')
r_c = rmin(TC)[0]
ax.annotate(f'inner surface reaches the axis at one point\n(x$_c$ = {XC:.4f}); in 3D its circumference → 0.\n'
            f'Mid-surface: R = ε, hoop κ = {1/r_c:.0f}, axial κ = {d2(XC, TC):.0f}\nOuter surface: R = 2ε, hoop κ = {1/(r_c+E):.0f}',
            xy=(XC, 0), xytext=(XC-.19, .105), fontsize=8.5, bbox=dict(fc='white', alpha=.85, ec='none'), arrowprops=dict(arrowstyle='->', color=INNER))

ax = fig.add_subplot(gs[0, 2])
zc = np.linspace(XC-.06, XC+.06, 2001)
band(ax, TC, '(c) PR rule at t$_c$ (zoom ×3): band before, caps after', alpha=.35, zx=zc)
for x, r, _, _ in m.branches(m.tc, 40001):
    k = (x > zc[0]) & (x < zc[-1])
    for s in (-1, 1): ax.plot(x[k], s*r[k], color=CAP, lw=2)
tips = [p[0][-1] for p in m.branches(m.tc, 4001)]
ax.annotate('', xy=(tips[0], -.035), xytext=(tips[1], -.035), arrowprops=dict(arrowstyle='<->', color='k'))
ax.text(np.mean(tips), -.04, f'tips appear {tips[1]-tips[0]:.4f} apart:\naxial length removed in one step', ha='center', va='top', fontsize=8.5, bbox=dict(fc='white', alpha=.85, ec='none'))
ax.text(zc[0]+.003, .057, 'green: caps drawn from t$_c$ on\ncap-tip curvature 164 and 140\n(radius of curvature 0.006–0.007 < ε)', fontsize=8.5, va='top', bbox=dict(fc='white', alpha=.85, ec='none'))

# (d) curvature approach, PR offset surface vs a smooth level-set (Morse) lumen
ax = fig.add_subplot(gs[1, 0])
s = np.logspace(-6, -0.6, 60)
rm = np.array([rmin(TC-q)[0] for q in s]); xm = np.array([rmin(TC-q)[1] for q in s])
rxx = np.array([d2(x, TC-q) for x, q in zip(xm, s)])
ax.loglog(s, 1/(rm-E), color=INNER, lw=2, label='PR inner surface R−ε: hoop curvature')
ax.loglog(s, 1/rm, color=MID, lw=2, label='mid-surface R: hoop curvature')
ax.loglog(s, rxx, color=MID, lw=1.2, ls='--', label='mid-surface: axial curvature')
ax.loglog(s, 1/(rm+E), color='0.5', lw=1.5, label='outer surface R+ε: hoop curvature')
rho = np.sqrt(rm**2-E**2)
ax.loglog(s, np.maximum(1/rho, rm*rxx/rho), color=MORSE, lw=2, label='level-set lumen ρ² = R²−ε²: max |κ|')
ax.axhline(1/E, color='k', lw=.7, ls=':'); ax.text(2e-6, 1/E*1.15, '1/ε', fontsize=9)
ax.text(1.5e-4, 2.5e4, 'slope −1', color=INNER, fontsize=9); ax.text(1.5e-4, 700, 'slope −½', color=MORSE, fontsize=9)
ax.set(xlabel='t$_c$ − t', ylabel='principal curvature', title='(d) approach to contact: what diverges')
ax.legend(fontsize=7.6, loc='upper right'); ax.grid(alpha=.3, which='both')

# (e) Morse alternative: continuous family through t_c with no imposed rule
ax = fig.add_subplot(gs[1, 1])
cols = plt.cm.viridis(np.linspace(0, .9, 5))
cols = ['#2c7fb8', '#000000', '#e6550d']
for c, q in zip(cols, [-.03, 0, .03]):
    r = m.radius(zx, TC+q); phi = r**2-E**2
    rho = np.where(phi > 0, np.sqrt(np.maximum(phi, 0)), np.nan)
    for sgn in (-1, 1): ax.plot(zx, sgn*rho, color=c, lw=1.8, label=f't − t$_c$ = {q:+.2f}' if sgn > 0 else None)
ax.axhline(0, color='k', lw=.6, ls=':')
ax.set(xlim=(XC-.05, XC+.05), ylim=(-.05, .05), xlabel='x', ylabel='ρ', title='(e) alternative lumen boundary ρ² = R² − ε², same R(x,t):\nneck → cone (at t$_c$ only) → two smooth caps')
ax.set_aspect('equal'); ax.legend(fontsize=8, loc='lower right')

# (f) the prescribed geometry would reopen; irreversibility is the rule
ax = fig.add_subplot(gs[1, 2])
ts = np.linspace(0, END, 1801)
for A, c, lab in [(.18, CAP, 'A = 0.18 (default)'), (.135, MORSE, 'A = b − ε = 0.135 (threshold)'), (.08, '0.5', 'A = 0.08 (weak)')]:
    mm = ContactNeck(type(m.par)(amplitude=A))
    ax.plot(ts, [mm.radius(xs, t).min()-E for t in ts], color=c, lw=2, label=lab)
ax.axhline(0, color='k', lw=.7)
ax.axvspan(m.tc, 4.9705, color=INNER, alpha=.12)
ax.text(5.05, .085, f'shaded: R < ε only for\nt ∈ [{m.tc:.3f}, 4.971] under the\nprescribed wave; the PR keeps\nthe caps separated after t$_c$', fontsize=8, bbox=dict(fc='white', alpha=.85, ec='none'))
ax.set(xlabel='t', ylabel='min$_x$ R − ε (half the facing gap)', title='(f) contact window; contact iff A > b − ε')
ax.legend(fontsize=8, loc='lower left'); ax.grid(alpha=.3)
fig.suptitle('PR #343 contact event: the mid-surface stays regular, the facing (inner) surface pinches at a point, and the caps appear by a jump', fontsize=12.5)
fig.savefig(sys.argv[1] if len(sys.argv) > 1 else 'contact.png', dpi=110)
print('t_c', TC, 'x_c', XC)
