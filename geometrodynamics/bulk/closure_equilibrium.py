"""Classical equilibrium measures for the closed spin-frame itinerary.

The Hamiltonian/preparation model is specified in
``docs/closure_equilibrium_prereg.md``. It is an additional classical model,
not a coupling derived from the existing BAM action or a local detector.
The baseline potential is K ||Ad_G - I||_F^2 / 16; a round rotor in a
canonical ensemble then has density exp(-V/kT) with respect to sphere area.

The whole-sphere partition integral is computed before taking any closure
limit. The limiting density is positive phase coarea. A different potential
with exactly the same zero set gives a different measure. Residual
reparametrizations preserve the measure only when the stiffness is
transformed as well. Outcome-sector priors remain inputs throughout.

No projectors, Born weights, Fock space, or quantum target enters the
partition functions. The loop sign from J^2 is invisible to Ad_G.
"""

from functools import lru_cache
import math

import numpy as np
from scipy.special import roots_legendre


MODELS = ("frame", "normal", "normal_rescaled", "normal_covariant")


def _parameters(c, model, quartic, epsilon):
    if not np.isfinite(c) or not -1.0 < c < 1.0:
        raise ValueError("c = u.w must be finite and strictly between -1 and 1")
    if model not in MODELS:
        raise ValueError(f"unknown closure potential: {model}")
    if not np.isfinite(quartic) or quartic < 0:
        raise ValueError("quartic must be finite and nonnegative")
    if model != "frame" and quartic != 0:
        raise ValueError("the quartic deformation is defined only for the frame potential")
    if not np.isfinite(epsilon) or abs(epsilon) >= 0.5:
        raise ValueError("abs(epsilon) must be below 1/2, ensuring positive g everywhere")


def _unit(v):
    v = np.asarray(v, dtype=float)
    if v.shape != (3,) or not np.all(np.isfinite(v)) or np.linalg.norm(v) == 0:
        raise ValueError("a direction must be a finite nonzero three-vector")
    return v / np.linalg.norm(v)


def triangle_data(x, u, w):
    """N and D of the geodesic triangle, for one or more unit source directions.

    x has shape (..., 3). Analyzer directions are normalized; x is checked,
    not silently projected. The antipodal geodesic punctures are undefined.
    """
    u, w = _unit(u), _unit(w)
    c = float(u @ w)
    _parameters(c, "frame", 0, 0)
    x = np.asarray(x, dtype=float)
    if x.ndim == 0 or x.shape[-1] != 3 or not np.all(np.isfinite(x)):
        raise ValueError("x must contain finite unit three-vectors")
    if not np.allclose(np.linalg.norm(x, axis=-1), 1.0, rtol=0, atol=1e-10):
        raise ValueError("source directions x must have unit norm")
    N = x @ np.cross(u, w)
    D = 1.0 + c + x @ (u + w)
    if np.any(N * N + D * D < 1e-28):
        raise ValueError("geodesic holonomy is undefined at x = -u or -w")
    return N, D


def residual_energy(residual, stiffness=1.0):
    """V/K = a F^2/2; a positive physical stiffness may vary over the sphere."""
    F, a = np.asarray(residual, dtype=float), np.asarray(stiffness, dtype=float)
    if not np.all(np.isfinite(F)) or not np.all(np.isfinite(a)) or np.any(a <= 0):
        raise ValueError("residuals must be finite and stiffness finite and positive")
    return 0.5 * a * F * F


def coarea_weight(gradient_norm, stiffness=1.0):
    """Local low-temperature density 1/(sqrt(a) |grad F|) on a regular arc."""
    g, a = np.asarray(gradient_norm, dtype=float), np.asarray(stiffness, dtype=float)
    if (not np.all(np.isfinite(g)) or not np.all(np.isfinite(a))
            or np.any(g <= 0) or np.any(a <= 0)):
        raise ValueError("gradient norms and stiffness must be finite and positive")
    return 1.0 / (np.sqrt(a) * g)


def _energy(N, D, c, model, quartic, epsilon):
    if model == "frame":
        f2 = N * N / (N * N + D * D)
        return 0.5 * (f2 + quartic * f2 * f2)
    if model == "normal":
        return residual_energy(N)
    # x.(u+w) = D - (1+c); g is positive for abs(epsilon) < 1/2.
    g = 1.0 + epsilon * (D - 1.0 - c)
    a = 1.0 / (g * g) if model == "normal_covariant" else 1.0
    return residual_energy(g * N, a)


def rotor_potential(x, u, w, *, model="frame", quartic=0.0, epsilon=0.25):
    """Dimensionless V/K at physical source directions on S^2.

    frame: sin^2(theta)/2, optionally plus quartic*sin^4(theta)/2.
    normal: N^2/2, a distinct restoring model with the same closure locus.
    normal_rescaled: (g N)^2/2, a changed stiffness profile.
    normal_covariant: (g N)^2/(2g^2), the original normal model.
    """
    u, w = _unit(u), _unit(w)
    c = float(u @ w)
    _parameters(c, model, quartic, epsilon)
    N, D = triangle_data(x, u, w)
    return _energy(N, D, c, model, quartic, epsilon)


@lru_cache(maxsize=16)
def _quadrature(n_normal, n_azimuth, beta):
    """Positive hemisphere nodes, refined normally at the thermal width.

    z=sinh(t)/sqrt(max(beta,1)) covers the WHOLE interval [0,1]. The
    change of variable resolves the boundary layer without a thin-shell
    cutoff. All four potentials are even in z. Normalized Haar integration
    is int_0^1 dz times the azimuthal average. Returned arrays are read-only.
    """
    r, weights = roots_legendre(n_normal)
    scale = math.sqrt(max(beta, 1.0))
    edge = math.asinh(scale)
    t = (r + 1.0) * edge / 2.0
    z = np.sinh(t) / scale
    wz = weights * (edge / 2.0) * np.cosh(t) / scale
    phi = 2.0 * math.pi * (np.arange(n_azimuth) + 0.5) / n_azimuth
    cosphi = np.cos(phi)
    for arr in (z, wz, cosphi):
        arr.setflags(write=False)
    return z, wz, cosphi


def canonical_partition(c, beta, *, model="frame", quartic=0.0,
                        epsilon=0.25, n_normal=128, n_azimuth=1024):
    """Z=(1/4pi) int_S2 exp(-beta V/K) dOmega for c=u.w and beta=K/kT.

    Uses no closure restriction or limiting formula. In the bisector frame,
    N=sqrt(1-c^2) z and D=1+c+sqrt(2(1+c))*sqrt(1-z^2)*cos(phi).
    This coordinate reduction uses isotropy of the stated round-rotor model.
    Check numerical convergence in both n_normal and n_azimuth for new regimes.
    """
    _parameters(c, model, quartic, epsilon)
    if not np.isfinite(beta) or beta < 0:
        raise ValueError("beta must be finite and nonnegative")
    for n in (n_normal, n_azimuth):
        if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 8:
            raise ValueError("quadrature sizes must be integers of at least 8")
    if beta == 0:
        return 1.0
    z, wz, cosphi = _quadrature(n_normal, n_azimuth, float(beta))
    N = math.sqrt(1.0 - c * c) * z[:, None]
    D = (1.0 + c + math.sqrt(2.0 * (1.0 + c))
         * np.sqrt(1.0 - z[:, None] ** 2) * cosphi[None, :])
    value = _energy(N, D, c, model, quartic, epsilon)
    return float(wz @ np.mean(np.exp(-beta * value), axis=1))


def limiting_mass(c, *, model="frame", quartic=0.0, epsilon=0.25):
    """Unnormalized coarea integral M: sqrt(beta/2pi)*Z -> M/(4pi).

    The quartic frame term does not change the transverse Hessian. These
    closed forms are independent oracles for canonical_partition, never
    used by that integrator.
    """
    _parameters(c, model, quartic, epsilon)
    cross = math.sqrt(1.0 - c * c)
    if model == "frame":
        return 4.0 + 2.0 * (1.0 + c) * math.acos(-c) / cross
    mass = 2.0 * math.pi / cross
    if model == "normal_rescaled":
        mass /= math.sqrt(1.0 - 2.0 * epsilon * epsilon * (1.0 + c))
    return mass


def limiting_density(phi, c, *, model="frame", quartic=0.0, epsilon=0.25):
    """Density on the closure circle, azimuth measured from u+w.

    The phase density has the continuous value zero at the two punctures;
    this is a limiting density, not a definition of the holonomy there.
    """
    _parameters(c, model, quartic, epsilon)
    dot = math.sqrt(2.0 * (1.0 + c)) * np.cos(np.asarray(phi))
    cross = math.sqrt(1.0 - c * c)
    if model == "frame":
        return np.abs(1.0 + c + dot) / cross
    if model == "normal_rescaled":
        return 1.0 / ((1.0 + epsilon * dot) * cross)
    return np.ones_like(dot, dtype=float) / cross


def _joint(gamma, mass, prior_ratio):
    if not np.isfinite(gamma) or not 0 < gamma < math.pi:
        raise ValueError("gamma must be strictly between 0 and pi")
    if not np.isfinite(prior_ratio) or prior_ratio <= 0:
        raise ValueError("the like/unlike prior ratio must be finite and positive")
    c = math.cos(gamma)
    # Singlet itinerary: c_s = u.w = -s_A*s_B*cos(gamma).
    like, unlike = prior_ratio * mass(-c), mass(c)
    Z = 2.0 * (like + unlike)
    if not np.isfinite(Z) or Z <= 0:
        raise ValueError("partition normalization vanished; refine the quadrature")
    return {(1, 1): like / Z, (-1, -1): like / Z,
            (1, -1): unlike / Z, (-1, 1): unlike / Z}


def joint_probabilities(gamma, beta, *, prior_ratio=1.0, **kwargs):
    """History-sector probabilities in the assumed canonical preparation.

    A detector event law and the prior_ratio are NOT derived by this function.
    beta=None requests the analytic low-temperature limit explicitly.
    """
    if beta is None:
        mass = lambda c: limiting_mass(c, **kwargs)
    else:
        mass = lambda c: canonical_partition(c, beta, **kwargs)
    return _joint(gamma, mass, prior_ratio)


def correlation(probabilities):
    """Correlation of a normalized four-sector table."""
    return float(sum(sa * sb * value for (sa, sb), value in probabilities.items()))


def standard_chsh(beta=None, **kwargs):
    """CHSH at the standard coplanar angles; no global maximization claim."""
    e1 = correlation(joint_probabilities(math.pi / 4, beta, **kwargs))
    e3 = correlation(joint_probabilities(3 * math.pi / 4, beta, **kwargs))
    return abs(3.0 * e1 - e3)


def monte_carlo_joint(gamma, beta, *, n_samples=200000, seed=20260905,
                      model="frame", quartic=0.0, epsilon=0.25):
    """Independent full-sphere importance estimate, with ratio standard errors.

    Uses unreduced 3-vectors for each of the four sectors and shared random
    samples, without enforcing equal marginals. This is numerical integration,
    not a simulated equilibration process or detector trajectory. Priors are
    equal for this control. Standard errors include the shared normalization.
    """
    if not np.isfinite(gamma) or not 0 < gamma < math.pi:
        raise ValueError("gamma must be strictly between 0 and pi")
    if not np.isfinite(beta) or beta < 0:
        raise ValueError("beta must be finite and nonnegative")
    if not isinstance(n_samples, int) or isinstance(n_samples, bool) or n_samples < 2:
        raise ValueError("n_samples must be an integer of at least 2")
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n_samples, 3))
    x /= np.linalg.norm(x, axis=1)[:, None]
    a, b = np.array([0.0, 0.0, 1.0]), np.array([math.sin(gamma), 0.0, math.cos(gamma)])
    sectors = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    ws = np.column_stack([
        np.exp(-beta * rotor_potential(x, sa * a, -sb * b, model=model,
                                      quartic=quartic, epsilon=epsilon))
        for sa, sb in sectors
    ])
    total = ws.sum(axis=1)
    norm = float(total.mean())
    if norm == 0:
        raise ValueError("no thermal weight was sampled; lower beta or increase samples")
    p = ws.mean(axis=0) / norm
    influence = (ws - total[:, None] * p) / norm
    se = influence.std(axis=0, ddof=1) / math.sqrt(n_samples)
    signs = np.array([sa * sb for sa, sb in sectors])
    e_se = float((influence @ signs).std(ddof=1) / math.sqrt(n_samples))
    return {"probabilities": dict(zip(sectors, p)), "standard_errors": dict(zip(sectors, se)),
            "correlation": float(p @ signs), "correlation_standard_error": e_se,
            "n_samples": n_samples, "seed": seed}
