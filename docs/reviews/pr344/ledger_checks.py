"""Independent checks of PR #344's sphere-side numbers: release work, cap refocusing, 2π return,
contact position and the causality of the release.

    PYTHONPATH=.:scripts python docs/reviews/pr344/ledger_checks.py     # from a checkout of #344
"""
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import minimize_scalar
from scipy.special import eval_legendre
from klein_release_model import KleinRelease

m = KleinRelease(); T = m.pulse_duration; w = m.omega
# 1. Work done by the source, integrated directly: W = ∫ F·p dt with the modes solved by an ODE
#    integrator (not Duhamel quadrature), starting from the incident state at t_c.
def rhs(t, y):
    a, p = y[:25], y[25:]; tau = (t-m.tc)/T
    F = m.scale*m.h*(30*tau**2*(1-tau)**2/T if 0 < tau < 1 else 0.)
    return np.r_[p, F-w*w*a]
a0, p0 = m.incident(m.tc)
sol = solve_ivp(rhs, (m.tc, m.tc+T), np.r_[a0, p0], rtol=1e-12, atol=1e-14, dense_output=True)
ts = np.linspace(m.tc, m.tc+T, 4001)
F = lambda t: m.scale*m.h*30*((t-m.tc)/T)**2*(1-(t-m.tc)/T)**2/T
W = np.trapezoid([F(t)@sol.sol(t)[25:] for t in ts], ts)
aT, pT = sol.sol(m.tc+T)[:25], sol.sol(m.tc+T)[25:]
ai, pi_ = m.incident(m.tc+T)
E = lambda a, p: .5*np.sum(p*p+w*w*a*a)
print(f'work ∫F·p dt = {W:.10f}; E(total) - E(incident) at pulse end = {E(aT, pT)-E(ai, pi_):.10f} (PR: 0.04)')
print(f'state vs PR Duhamel at pulse end: max |diff| = {max(np.abs(aT-m.state(m.tc+T)[0]).max(), np.abs(pT-m.state(m.tc+T)[1]).max()):.1e}')
# 2. Emitted cap energy with a separate quadrature (fine uniform grid in theta).
def cap(t):
    a, p = m.state(t, emitted_only=True); th = np.r_[np.linspace(0, .3, 3001), np.linspace(np.pi-.3, np.pi, 3001)]
    P = np.array([eval_legendre(l, np.cos(th)) for l in range(25)]).T*m.norm
    u = P@a; ut = P@p
    # d/dtheta via analytic derivative of the Legendre series in theta
    eps = 1e-6; up = (np.array([eval_legendre(l, np.cos(th+eps)) for l in range(25)]).T*m.norm)@a
    um = (np.array([eval_legendre(l, np.cos(th-eps)) for l in range(25)]).T*m.norm)@a
    dens = .5*(ut**2+((up-um)/(2*eps))**2+.25*u*u)*np.sin(th)*2*np.pi
    return np.trapezoid(dens[:3001], th[:3001])+np.trapezoid(dens[3001:], th[3001:])
c = m.tc+T/2
print(f'cap energy midflight {cap(c+np.pi/2):.6e}, antipode {cap(c+np.pi):.6e}, ratio {cap(c+np.pi)/cap(c+np.pi/2):.3f} (PR 47.18)')
# 3. 2π return: omega_l = l + 1/2, so every mode picks up cos(2π(l+1/2)) = -1.
a1, p1 = m.state(c+3); a2, p2 = m.state(c+3+2*np.pi)
print(f'2π return: max |a(t+2π) + a(t)| = {np.abs(a1+a2).max():.1e}')
# 4. Where contact happens, and when a unit-speed signal from it reaches each rim.
x = np.linspace(-.6, .6, 4001); i = np.argmin(m.radius(x, m.tc))
xc = minimize_scalar(lambda q: m.radius(q, m.tc), bounds=(x[i-2], x[i+2]), method='bounded', options={'xatol': 1e-12}).x
dA, dB = (xc+.6)/1.2*m.length, (.6-xc)/1.2*m.length
frac = lambda d: quad(lambda s: 30*s*s*(1-s)**2, 0, min(d/T, 1))[0]
print(f'contact x_c = {xc:+.4f}: intrinsic distance to rim A {dA:.4f}, rim B {dB:.4f}; release pulse T = {T}')
print(f'share of the release time profile emitted before a unit-speed signal from contact could reach A: {frac(dA):.3f}, B: {frac(dB):.3f}')
# 5. Geometry: #338 neck vs #343 profile used for contact.
a = .3; S, b = np.sin(a)*np.cos(a), np.sin(a)**2
print(f'#338 neck: half-length S = {S:.4f}, waist b = {b:.4f}, rim sin a = {np.sin(a):.4f}; in rim units h/a = {S/np.sin(a):.3f}, waist {b/np.sin(a):.4f}')
print(f'#343 profile used for contact: h/a = 0.6, waist 0.16; implied scale axial {m.length/1.2:.4f}, radial {np.sin(a):.4f} sphere units per neck unit')
print(f'contact needs gain > waist - r_half: #343 profile {.16-m.neck.level_radius:.4f}; #338 waist {b/np.sin(a)-m.neck.level_radius:.4f}')
