# c(theta): minimal action of one throat from x to e^{i theta} x, user's L = |zdot|^2 + lambda_s A^2.
# #367's arc_action is S = 1/2 int(|x_h|^2 + lam nu^2) with lam = 1 + lambda_s, so c = 2 S.
import numpy as np, sys
from scipy.optimize import minimize
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[4]))
from experiments.closure_ledger.antipodal_network.antipodal_two_time import arc_action, unit
n = 24
z = np.array([np.cos(0.6), np.sin(0.6) * np.exp(0.4j)])
x = np.array([z[0].real, z[0].imag, z[1].real, z[1].imag])
def target(theta):
    w = np.exp(1j * theta) * z
    return np.array([w[0].real, w[0].imag, w[1].real, w[1].imag])
def c_num(theta, lam, seeds=(1, 2, 3, 4)):
    end = target(theta); best = np.inf
    for s in seeds:
        rng = np.random.default_rng(s)
        # start from a mix of the vertical path and a random detour
        t = np.linspace(0, 1, n + 1)[1:-1, None]
        w = np.exp(1j * theta * t) * z[None, :]
        base = np.column_stack([w[:, 0].real, w[:, 0].imag, w[:, 1].real, w[:, 1].imag])
        v = unit(rng.normal(size=4))
        guess = unit(base + 0.8 * np.sin(np.pi * t) * v)
        obj = lambda p: arc_action(np.vstack((x, unit(p.reshape(-1, 4)), end)), lam)
        r = minimize(obj, guess.ravel(), method="BFGS", options=dict(gtol=1e-9, maxiter=40000))
        best = min(best, r.fun)
    return 2 * best
def c_paper(theta, ls):
    th = abs((theta + np.pi) % (2 * np.pi) - np.pi)         # fold to [0, pi]
    psi = np.pi - th
    near_pi = np.pi ** 2 - (1 + ls) * psi ** 2 / ls if np.pi - psi * (1 + 1 / ls) >= 0 else np.inf
    return min(near_pi, (1 + ls) * th ** 2, th * (2 * np.pi - th))
for ls in (2.0, 0.5):
    lam = 1 + ls
    print('lambda_s =', ls)
    for th in (np.pi, np.pi - 0.2, np.pi - 0.5, np.pi - 1.0, np.pi / 2, 0.5, 0.0):
        print('  theta %.3f  c numeric %.4f  c paper %.4f' % (th, c_num(th, lam), c_paper(th, ls)), flush=True)
    # shared constraint theta_A + theta_B = pi (mod 2 pi): total over the split
    grid = np.linspace(0, np.pi, 9)
    tot = [(th, c_num(th, lam) + c_num(np.pi - th, lam)) for th in grid]
    print('  theta_A + theta_B = pi: totals', ' '.join('%.2f:%.3f' % t for t in tot))
