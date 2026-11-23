"""Full elliptic pulse followed by the Gaussian solution of #209's winding-force window."""
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit, ndtr, ndtri

from .elliptic_pointer_ring import CHI, SITES, WINDING

SIGMA, FORCE, TAU, END = 2., 0.5, 2., 20.
THETAS = np.arange(9)*np.pi/16


def split_weights(theta, eps, axis=np.pi/2, enabled=True):
    field = (WINDING+np.exp(2j*theta)*WINDING.conj())/np.sqrt(2)
    if enabled:
        # Same link geometry as #241, with an independently selectable orientation.
        diff = np.roll(np.eye(SITES), -1, axis=0)-np.eye(SITES)
        radius = 1+eps*np.cos(2*(CHI+np.pi/SITES)+axis)
        omega, vectors = np.linalg.eigh(diff.T @ (diff/radius[:, None]**2))
        duration = np.pi/(4*eps*(2-np.sqrt(2)))
        field = vectors @ (np.exp(-1j*duration*omega)*(vectors.T @ field))
    weights = abs(np.fft.fft(field)/np.sqrt(SITES))**2
    # The exact quadrupole preserves odd winding. For k=±1,±3 the
    # k-odd connection force sin(kπ/4)/sin(π/4) is exactly ±1.
    if np.sum(weights[::2]) > 1e-18:
        raise ValueError("even winding excited; two spatial packets no longer suffice")
    return float(sum(weights[[1,3]])), float(sum(weights[[3,5]]))


def packet(t):
    w = 2*np.pi/TAU
    if t <= TAU:
        d = FORCE*(t*t/2+(np.cos(w*t)-1)/w**2)
        u = FORCE*(t-np.sin(w*t)/w)
    else:
        d = FORCE*TAU*(t-TAU/2)
        u = FORCE*TAU
    width = SIGMA*np.sqrt(1+(t/(2*SIGMA**2))**2)
    return d, u, width, t/(4*SIGMA**4+t*t)


def transport(weight, ranks, scale):
    # Each initial position evolves with J/rho. No use of a random branch label.
    x0 = scale*SIGMA*ndtri(ranks)
    log_odds = np.log(np.clip(weight, 1e-300, 1))-np.log(np.clip(1-weight, 1e-300, 1))
    def rhs(t, x):
        d, u, width, dilation = packet(t)
        posterior = expit(log_odds+2*x*d/width**2)
        return dilation*x+(2*posterior-1)*(u-dilation*d)
    times = (2., 5., 10., END)
    sol = solve_ivp(rhs, (0, END), x0, t_eval=times, method="DOP853",
                    rtol=1e-10, atol=1e-12, max_step=0.25)
    if not sol.success:
        raise RuntimeError(sol.message)
    rank_error = 0.
    for t, x in zip(times, sol.y.T):
        d, _, width, _ = packet(t)
        cdf = weight*ndtr((x-d)/width)+(1-weight)*ndtr((x+d)/width)
        rank_error = max(rank_error,float(max(abs(cdf-ndtr(x0/SIGMA)))))
    d, _, width, _ = packet(END)
    right_mass = weight*ndtr(d/width)+(1-weight)*ndtr(-d/width)
    return dict(right_fraction=float(np.mean(sol.y[:,-1]>0)), right_mass=float(right_mass),
                max_rank_error=rank_error)


def run(eps, axis=np.pi/2, enabled=True, count=2048):
    ranks = (np.arange(count)+0.5)/count
    rows = []
    for theta in THETAS:
        weight, leakage = split_weights(theta,eps,axis,enabled)
        rows.append(dict(theta=float(theta), target=float(np.cos(theta)**2),
                         positive_force_weight=weight, leakage=leakage,
                         equilibrium=transport(weight,ranks,1.),
                         narrow=transport(weight,ranks,0.6)))
    return rows


def probe():
    d, _, width, _ = packet(END)
    # Balanced doublet with opposite freely accelerating windings at x=0.
    # rho_charge/rho_total=0; J_charge/rho_total = u - dilation*d.
    d2,u2,_,a2 = packet(TAU)
    # Orientation of the actual mouth, distinct from a displacement prior.
    axes = 2*np.pi*(np.arange(128)+0.5)/128
    phase_average = [float(np.mean([split_weights(th,.02,a)[0] for a in axes]))
                     for th in THETAS]
    return dict(random_mixer_axis_positive_force_weights=phase_average, count=2048, centroid_separation_in_widths=2*d/width,
                single_branch_wrong_side=float(ndtr(-d/width)),
                balanced_signed_current_at_center=dict(time=TAU, signed_density_over_norm_density=0.,
                                                       signed_current_over_norm_density=u2-a2*d2),
                primary=run(.02), halved_ellipticity=run(.01),
                mixer_off=run(.02,enabled=False), wrong_axis=run(.02,axis=0.))


if __name__ == "__main__":
    print(json.dumps(probe(), indent=2))
