# Direct check: one throat, fixed preparation, port-symmetric load -(k/2)(2v-1)^2 at the midpoint, Berger lambda = 2.
import numpy as np, sys
from scipy.optimize import minimize, minimize_scalar
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[4]))
from experiments.closure_ledger.antipodal_network.antipodal_two_time import arc_action, hopf, rank, unit
n = 8
def solve(x, kappa, seed):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=4); v = unit(v - v.dot(x) * x)
    t = np.arange(1, n)[:, None] * np.pi / n
    guess = unit(np.cos(t) * x + np.sin(t) * v + 0.05 * rng.normal(size=(n - 1, 4)))
    def nodes(p): return np.vstack((x, unit(p.reshape(-1, 4)), -x))
    def obj(p):
        N = nodes(p); return arc_action(N, 2.) - kappa / 2 * (2 * rank(N[n // 2]) - 1) ** 2
    r = minimize(obj, guess.ravel(), method="BFGS", options=dict(gtol=1e-10, maxiter=20000))
    N = nodes(r.x)
    fiber = np.abs(np.einsum('ij,ij->i', hopf(N[:-1]), N[1:])).max()     # largest vertical component per segment
    return r.fun, rank(N[n // 2]), fiber
def reduced(u, kappa):
    al = np.arccos(np.sqrt(u))
    f = lambda d: 2 * d * d - kappa / 2 * np.cos(2 * al + 2 * d) ** 2
    d = minimize_scalar(f, bounds=(-0.6, 0.6), method="bounded", options=dict(xatol=1e-12)).x
    return np.sin(al + d) ** 2
for kappa in (0.1, 0.3):
    for u in (0.2, 0.35, 0.5, 0.8):
        al = np.arccos(np.sqrt(u)); phase = np.exp(0.7j)
        z = np.array([np.cos(al), np.sin(al) * phase])      # a complex preparation, not just the real plane
        x = np.array([z[0].real, z[0].imag, z[1].real, z[1].imag])
        best = min((solve(x, kappa, s) for s in (1, 2, 3)), key=lambda r: r[0])
        sigma = u + best[1] - 1
        print('kappa %.1f u %.2f  defect %+.5f  reduced 1D %+.5f  first order %+.5f  max fibre component %.1e'
              % (kappa, u, sigma, u + reduced(u, kappa) - 1, -2 * kappa * u * (1 - u) * (2 * u - 1), best[2]))
