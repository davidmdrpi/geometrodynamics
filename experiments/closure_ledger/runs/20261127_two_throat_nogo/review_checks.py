# Reproduces the five review findings on the two-throat note (PR #374 at 215aace).
# Action convention: S = 1/2 int L dt with L = |zdot|^2 + lambda_s A^2 (+ 2 lambda_c A_A A_B), loads added to S.
# This is #367's arc_action. Integral-of-L values are 2 S.
import numpy as np, sys
from scipy.optimize import minimize_scalar
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[4]))
from experiments.closure_ledger.antipodal_network.antipodal_two_time import arc_action, hopf, rank


def r4(z):
    return np.array([z[0].real, z[0].imag, z[1].real, z[1].imag])


def nodes(x, v, n):
    t = np.linspace(0, 1, n + 1)[:, None]
    return np.cos(np.pi * t) * x + np.sin(np.pi * t) * v


def joint_S(NA, NB, ls, lc):
    """Discrete joint action: two arc actions plus the cross term on the per-segment fibre velocities."""
    n = len(NA) - 1
    def fib(N):
        d = np.arccos(np.clip(np.sum(N[:-1] * N[1:], axis=-1), -1, 1))
        return np.sum(hopf(N[:-1]) * N[1:], axis=-1) / np.sinc(d / np.pi)
    cross = n * np.sum(2 * lc * fib(NA) * fib(NB)) / 2
    return arc_action(NA, 1 + ls) + arc_action(NB, 1 + ls) + cross


print("1. Non-horizontal counterexample: lambda_s = 2, lambda_c = 1.99, kappa = 0.3, x = (1,1)/sqrt2, first port")
ls, lc, kappa = 2.0, 1.99, 0.3
x = r4(np.array([1, 1]) / np.sqrt(2))
e1 = r4(np.array([1, 0]))
for n in (8, 64):
    NA, NB = nodes(x, hopf(e1), n), nodes(x, -hopf(e1), n)      # z = cos(pi t) x +/- i sin(pi t) e1
    S = joint_S(NA, NB, ls, lc)
    loads = sum(kappa / 2 * (2 * rank(N[n // 2]) - 1) ** 2 for N in (NA, NB))
    print("   n=%2d  midpoint intensities %.6f %.6f  S - loads = %.5f   (int L - loads = %.5f)"
          % (n, rank(NA[n // 2]), rank(NB[n // 2]), S - loads, 2 * S - loads))
print("   continuum: S = pi^2 + (pi^2/2)(lambda_s - lambda_c) = %.5f, minus loads 0.3 = %.5f"
      % (np.pi ** 2 + np.pi ** 2 / 2 * (ls - lc), np.pi ** 2 + np.pi ** 2 / 2 * (ls - lc) - 0.3))
f = lambda d: 2 * d * d - kappa / 2 * np.sin(2 * d) ** 2              # horizontal shift per throat, S units, u = 1/2
d = minimize_scalar(f, bounds=(-1.0, 1.0), method="bounded").x
print("   horizontal floor: best shift %.2e, joint S = %.5f (= pi^2 %.5f); in int-L units %.5f (= 2 pi^2)"
      % (d, np.pi ** 2 + 2 * f(d), np.pi ** 2, 2 * np.pi ** 2 + 2 * f(d)))
print("   local comparison per throat at small midpoint shift rho: horizontal 2 rho^2,"
      " uniform vertical tilt (pi^2/2) lambda rho^2, lambda = lambda_s -+ lambda_c by mode;"
      " 4/pi^2 = %.4f" % (4 / np.pi ** 2))

print("2. Shared holonomy with the cross term: lambda_s = 2, lambda_c = -1.5, theta_A + theta_B = pi")
ls, lc = 2.0, -1.5
th = np.pi / 2
print("   equal pure-fibre rotations through pi/2: int L = %.5f = %.3f pi^2" % (
    (1 + ls) * 2 * th ** 2 + 2 * lc * th ** 2, ((1 + ls) * 2 * th ** 2 + 2 * lc * th ** 2) / np.pi ** 2))
print("   frozen split (pi, 0): int L >= c(pi) = pi^2 = %.5f" % np.pi ** 2)

print("3. Central preparation u = 1/2 at stronger loads (S units): cost 2 d^2 - (kappa/2) sin^2(2d), curvature 4(1 - kappa)")
for kappa in (0.5, 0.99, 1.5, 2.0):
    best = min((minimize_scalar(lambda d: 2 * d * d - kappa / 2 * np.sin(2 * d) ** 2, bounds=b, method="bounded",
                                options=dict(xatol=1e-12)) for b in ((-1.0, 0.0), (0.0, 1.0))), key=lambda r: r.fun)
    v = np.sin(np.pi / 4 + best.x) ** 2
    print("   kappa %.2f  shift %+.4f  defect sigma = u + v - 1 = %+.4f (and %+.4f at the mirror minimum)"
          % (kappa, best.x, 0.5 + v - 1, -(0.5 + v - 1)))

print("5. Figure 2: handle exit versus cover antipode of the pulse source")
from geometrodynamics.viz import alice_strip as st
c = st.COLUMNS[st.WAVE_SOURCE]
p = st._on_sphere(c, True) / st.SHEET_RADIUS
print("   source . handle exit (bottom c) = %.5f;  source . refocus point (bottom -c) = %.5f"
      % (p @ (st._on_sphere(c, False) / st.SHEET_RADIUS), p @ (st._on_sphere(-c, False) / st.SHEET_RADIUS)))
