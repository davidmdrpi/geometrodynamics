"""Diffuse-shell alternative for PR #343: the same prescribed R(x,t), drawn as a smooth density.

psi(x, y, t) = 1 - (1 - G(y - R))(1 - G(y + R)), G(s) = exp(-s^2/eps^2), on the meridian plane
(y = signed rho): a smooth union of the two wall profiles, saturating at 1 where they overlap.
psi is even in y, so it is a smooth function of rho^2 and hence of (x, y, z) in 3D. Nothing
is imposed at contact: the walls' density overlaps on the axis and the psi = 1/2 level set
changes topology at one Morse point while psi and all its derivatives stay bounded.

    PYTHONPATH=scripts python docs/reviews/pr343/diffuse_shell.py docs/reviews/pr343/diffuse.png   # from a checkout of #343
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from neck_wave_contact import ContactNeck

m = ContactNeck(); E = m.par.thickness; TC, XC = m.tc, m.xc
x = np.linspace(XC-.12, XC+.12, 721); y = np.linspace(-.12, .12, 721)
Xg, Yg = np.meshgrid(x, y)
G = lambda s: np.exp(-(s/E)**2)
def psi(t):
    R = m.radius(Xg, t)
    return 1-(1-G(Yg-R))*(1-G(Yg+R))
def hessian_max(p):
    dx, dy = x[1]-x[0], y[1]-y[0]
    pxx = np.gradient(np.gradient(p, dx, axis=1), dx, axis=1)
    pyy = np.gradient(np.gradient(p, dy, axis=0), dy, axis=0)
    pxy = np.gradient(np.gradient(p, dx, axis=1), dy, axis=0)
    return np.max(np.sqrt(pxx**2+pyy**2+2*pxy**2))
ts = np.linspace(TC-.6, TC+.6, 121)
axis = [psi(t)[len(y)//2].max() for t in ts]
hess = [hessian_max(psi(t)) for t in ts]
# time at which the psi = 1/2 level set first connects across the axis (its Morse instant)
t_half = ts[np.argmax(np.array(axis) >= .5)]
a = np.array(axis); t10, t90 = ts[np.argmax(a >= .1)], ts[np.argmax(a >= .9)]
print(f'axis density rises 0.1 -> 0.9 over t - t_c in [{t10-TC:+.3f}, {t90-TC:+.3f}], duration {t90-t10:.3f} (eps/v = {E/.12419:.3f})')
print(f'axis density at t_c: {psi(TC)[len(y)//2].max():.4f}; psi=1/2 level set closes the lumen at t ~ {t_half:.3f} '
      f'(t_c = {TC:.3f}); max Hessian norm over [t_c-.6, t_c+.6]: {max(hess):.0f} (1/eps^2 = {1/E**2:.0f})')
fig, axs = plt.subplots(1, 4, figsize=(16, 4.1), constrained_layout=True)
for ax, t, lab in zip(axs[:3], [TC-.3, t_half, TC+.3], ['t$_c$ − 0.3', f'{t_half:.2f} (ψ = ½ closes)', 't$_c$ + 0.3']):
    p = psi(t)
    im = ax.imshow(p, extent=(x[0], x[-1], y[0], y[-1]), origin='lower', cmap='magma', vmin=0, vmax=1)
    ax.contour(Xg, Yg, p, [.5], colors='w', linewidths=1.2)
    ax.set(title=f'shell density ψ, t = {lab}', xlabel='x', ylabel='ρ (signed)'); ax.set_aspect('equal')
fig.colorbar(im, ax=axs[2], shrink=.8)
ax = axs[3]
ax.plot(ts-TC, axis, lw=2, label='max ψ on the axis')
ax.plot(ts-TC, np.array(hess)*E**2, lw=2, label='max |∇²ψ| × ε²')
ax.axhline(.5, color='k', lw=.7, ls=':'); ax.axvline(0, color='0.5', lw=.7)
ax.set(xlabel='t − t$_c$', title='bounded and continuous through t$_c$'); ax.legend(fontsize=9); ax.grid(alpha=.3)
fig.suptitle('Alternative: the same prescribed walls as a smooth density of width ε. The lumen fills across a band ~ε wide in 0.22 time units (≈ ε/v), with no rule and no jump', fontsize=11.5)
fig.savefig(sys.argv[1] if len(sys.argv) > 1 else 'diffuse.png', dpi=105)
