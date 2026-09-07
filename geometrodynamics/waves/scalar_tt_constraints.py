"""Conditional, linearized ESU constraint response for the smooth n=3 scalar.

Freeze e1706b1. Each time is an independent CMC slice with zero supporting
matter perturbations. These fields are not an Einstein evolution, and their
norm bounds do not bound the complete omitted scalar backreaction.
"""

from functools import lru_cache
import math

import numpy as np

from . import reciprocal_scalar_tt as rt


PUBLIC_PREREG = "e1706b166b24b69afceda85932e88a43b7fb61f9"
BASELINE = "67940cb32d26ffed0c24634348ed9896d3964d64"
SEED = 2026090712


def solve_hamiltonian(rho, degrees, radius=1., kappa=1.):
    """Invert on supplied scalar modes; never discard a sourced dipole."""
    rt.ReciprocalModel(radius=radius, kappa=kappa)
    rho, degrees = np.asarray(rho, dtype=float), np.asarray(degrees)
    if rho.shape != degrees.shape or np.any(degrees < 0) or np.any(degrees % 1):
        raise ValueError("one nonnegative integer degree per coefficient required")
    if not np.all(np.isfinite(rho)):
        raise ValueError("finite source coefficients required")
    kernel = degrees == 1
    if np.any(rho[kernel] != 0):
        raise ValueError("nonzero dipole source: Hamiltonian constraint incompatible")
    out = np.zeros_like(rho)
    out[~kernel] = -kappa*radius**2*rho[~kernel]/(4*(3-degrees[~kernel]*(degrees[~kernel]+2)))
    return out


class EvenHarmonicGrid:
    """All 84 real scalar harmonics of degrees 0,2,4,6, physically normalized."""

    def __init__(self, radial_order=8, angular_order=16, radius=1.):
        self.radius = float(radius)
        self.volume = 2*math.pi**2*self.radius**3
        self.points, self.weights = rt.sphere_quadrature(radial_order, angular_order, radius)
        columns = [np.ones((len(self.points), 1))/math.sqrt(self.volume)]
        degrees = [0]
        self.D = np.zeros((3, 84, 84))
        offset = 1
        for degree in (2, 4, 6):
            h = rt.harmonic_multiplet(degree)
            columns.append(rt.monomials(self.points, h.exponents) @ h.B/math.sqrt(self.volume))
            size = (degree+1)**2
            self.D[:, offset:offset+size, offset:offset+size] = h.D/self.radius
            degrees.extend([degree]*size)
            offset += size
        self.Y = np.concatenate(columns, axis=1)
        self.degrees = np.array(degrees)
        self.lambdas = self.degrees*(self.degrees+2)/self.radius**2
        self.frames = np.einsum("iab,pb->pia", rt.QUATERNION_DERIVATIVES, self.points)

    def project(self, values):
        return self.Y.T @ (self.weights*np.asarray(values))

    def evaluate(self, coefficients):
        return self.Y @ coefficients

    def gradient(self, coefficients):
        return np.stack([self.Y @ (d @ coefficients) for d in self.D], axis=-1)

    def laplacian_coefficients(self, coefficients):
        # Differentiated polynomial representations, independent of lambdas.
        return sum(d @ d @ coefficients for d in self.D)

    def longitudinal_coefficients(self, w):
        # 2 Hess(w) - (2/3) g Delta(w), in the invariant frame.
        lap = self.laplacian_coefficients(w)
        return np.array([[(self.D[i] @ self.D[j]+self.D[j] @ self.D[i]) @ w
                          - (2/3)*(i == j)*lap for j in range(3)] for i in range(3)])

    def divergence_coefficients(self, tensor_coefficients):
        # nabla^j A_ij, including both connection contractions explicitly.
        # Gamma_{ji}^k = epsilon_{jik}/a, Gamma_{jj}^k = 0 on round SU(2).
        eps = levi_civita()
        return np.array([sum(self.D[j] @ tensor_coefficients[i, j]
                             - sum(eps[j, i, k]*tensor_coefficients[k, j]/self.radius
                                   + eps[j, j, k]*tensor_coefficients[i, k]/self.radius
                                   for k in range(3)) for j in range(3)) for i in range(3)])

    def rms(self, coefficients):
        return float(np.linalg.norm(coefficients)/math.sqrt(self.volume))

    def momentum_charges(self, momentum):
        """Six ambient rotations and four gradient dipoles (unit frame)."""
        ambient = np.einsum("pi,pia->pa", momentum, self.frames)
        rotations = [self.weights @ (self.points[:, a]*ambient[:, b]
                                      - self.points[:, b]*ambient[:, a])
                     for a in range(4) for b in range(a+1, 4)]
        dipoles = np.einsum("p,pa->a", self.weights, ambient)/self.radius
        return np.r_[rotations, dipoles]

    def require_compatible_momentum(self, momentum, tolerance=1e-10):
        charges = self.momentum_charges(momentum)
        if np.max(np.abs(charges)) > tolerance:
            raise ValueError("momentum source has conformal Killing charge; no CMC solution")
        return charges


def levi_civita():
    eps = np.zeros((3, 3, 3))
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        eps[i, j, k], eps[i, k, j] = 1., -1.
    return eps


def round_sources(model, q, p, points):
    """Independent inherited improved stress on the round ESU, at a=1."""
    if model.radius != 1.:
        raise ValueError("inherited stress hardcodes radius=1")
    from .backreaction import stress_series
    jets, _ = rt.scalar_jets(model, q, p, -model.omega_scalar2*q, points)
    stress = stress_series(jets)[:, 0]
    return stress[:, 0, 0], -stress[:, 0, 1:]


class StandingWaveConstraints:
    """Scalar-only particular solutions for separately considered CMC slices."""

    def __init__(self, grid, amplitude=.2, kappa=1., axis=None):
        self.grid = grid
        self.model = rt.ReciprocalModel(radius=grid.radius, kappa=kappa)
        self.amplitude = float(amplitude)
        if not math.isfinite(self.amplitude):
            raise ValueError("finite amplitude required")
        axis = np.array([1., 2., 3.])/math.sqrt(14) if axis is None else np.asarray(axis)
        self.direction = self.model.multiplet.coherent(axis)
        self.mode = (rt.monomials(grid.points, self.model.multiplet.exponents)
                     @ self.model.multiplet.B @ self.direction/math.sqrt(grid.volume))
        self.h = grid.project(self.mode**2)

    def source_coefficients(self, time):
        omega = math.sqrt(self.model.omega_scalar2)
        A = self.amplitude*math.cos(omega*time)
        dA = -self.amplitude*omega*math.sin(omega*time)
        rho = (.5*self.amplitude**2*omega**2-self.grid.lambdas*A*A/12)*self.h
        J = -A*dA*self.h/6
        return rho, J, A, dA

    def response(self, time):
        rho, J, A, dA = self.source_coefficients(time)
        u = solve_hamiltonian(rho, self.grid.degrees, self.grid.radius, self.model.kappa)
        w = np.zeros_like(J)
        w[1:] = 3*self.model.kappa*J[1:]/(4*(3/self.grid.radius**2-self.grid.lambdas[1:]))
        return {"rho": rho, "J": J, "u": u, "w": w, "A": A, "dA": dA}

    def norms(self, time):
        r = self.response(time)
        g, kappa, a = self.grid, self.model.kappa, self.grid.radius
        u, rho, J, w = (r[key] for key in ("u", "rho", "J", "w"))
        j = math.sqrt(float(np.sum(g.lambdas*J*J)/g.volume))
        kl = math.sqrt(max(0., float((8/3)*np.sum(g.lambdas*(g.lambdas-3/a**2)*w*w)/g.volume)))
        return {"rho_rms": g.rms(rho), "rho_inhomogeneous_rms": g.rms(rho[1:]),
                "u_mean": float(u[0]/math.sqrt(g.volume)), "u_rms": g.rms(u),
                "u_inhomogeneous_rms": g.rms(u[1:]),
                "u_rms_upper": kappa*a*a*g.rms(rho)/12,
                "u_inhomogeneous_lower": kappa*a*a*g.rms(rho[1:])/180,
                "u_inhomogeneous_upper": kappa*a*a*g.rms(rho[1:])/20,
                "j_rms": j, "K_longitudinal_rms": kl,
                "K_longitudinal_lower": kappa*a*j/math.sqrt(30),
                "K_longitudinal_upper": kappa*a*math.sqrt(3/10)*j,
                "scalar_metric_rms": 4*math.sqrt(3)*g.rms(u),
                "scalar_metric_inhomogeneous_rms": 4*math.sqrt(3)*g.rms(u[1:]),
                "induced_TT_metric_frobenius": float(2*np.linalg.norm(self.induced_tt(time)))}

    def induced_tt(self, time):
        T, S = self.model.omega_tensor2, self.model.omega_scalar2
        cT, cS = math.cos(math.sqrt(T)*time), math.cos(2*math.sqrt(S)*time)
        B = (1-cT)/(2*T)+(cS-cT)/(2*(T-4*S))
        return self.model.source(self.amplitude*self.direction)*B/self.model.C

    def envelopes(self):
        """Continuous extrema, not sampled maxima. Valid for all times."""
        qperiod = math.pi/(2*math.sqrt(self.model.omega_scalar2))
        u1, u0 = self.response(0.)["u"], self.response(qperiod)["u"]
        out = {}
        for name, sl in (("u", slice(None)), ("u_inhomogeneous", slice(1, None))):
            z, d = u0[sl], (u1-u0)[sl]
            r = float(np.clip(-(z @ d)/(d @ d), 0., 1.)) if d @ d else 0.
            out[name] = {"minimum_rms": self.grid.rms(z+r*d),
                         "minimum_at_cos_squared": r,
                         "maximum_rms": max(self.grid.rms(u0[sl]), self.grid.rms(u1[sl]))}
        peak = self.norms(qperiod/2)
        out["K_longitudinal_maximum_rms"] = peak["K_longitudinal_rms"]
        out["j_maximum_rms"] = peak["j_rms"]
        T, S = self.model.omega_tensor2, self.model.omega_scalar2
        out["induced_TT_metric_all_time_upper"] = float(
            2*np.linalg.norm(self.model.source(self.amplitude*self.direction))/self.model.C
            *(1/T+1/abs(T-4*S)))
        # Independent exact l=2 certificate, derived after the freeze.
        # ||h_2||^2 = 27/(25 V), u_2 >= (11/30) kappa s^2 h_2 in norm.
        out["scalar_metric_inhomogeneous_all_time_lower"] = (
            66*self.model.kappa*self.amplitude**2/(25*self.grid.volume))
        return out


@lru_cache(None)
def homogeneous_tt_constraint_certificate():
    """Compute the old report's zeros from ADM curvature and connection."""
    import sympy as sp
    from geometrodynamics.bulk.tt_triangle_rotor import adm_quadratic_derivation
    adm = adm_quadratic_derivation()
    delta_g00 = sp.simplify(sp.sympify(adm["linear_scalar"])/2)
    a = sp.Symbol("a", positive=True)
    x, y, z, v, w = sp.symbols("x y z v w", real=True)
    rate = sp.Matrix([[x, z, v], [z, y, w], [v, w, -x-y]])
    momentum = [sp.simplify(-sum(sp.LeviCivita(j, i, k)*rate[k, j]/a
                                + sp.LeviCivita(j, j, k)*rate[i, k]/a
                                for j in range(3) for k in range(3))) for i in range(3)]
    return {"delta_G00": float(delta_g00),
            "delta_G0i_max_absolute": max(float(abs(v)) for v in momentum),
            "linear_spatial_curvature": adm["linear_scalar"],
            "momentum_connection_contraction": [str(v) for v in momentum],
            "exact": bool(adm["linear_vanishes"] and delta_g00 == 0
                          and all(v == 0 for v in momentum)),
            "scope": "linear homogeneous STF perturbation; background K=0; not the matter-sourced completion"}


@lru_cache(None)
def coherent_square_certificate():
    """Exact harmonic powers of the normalized coherent Y_3 squared.

    R=x0^2+(m.x)^2 is uniform on [0,1] under normalized S3 volume.
    Y^2=(4/V)[R^3+Re(x0+i m.x)^6]. Legendre modes have l=2k;
    the last term is a separate l=6 harmonic with phase average zero.
    """
    import sympy as sp
    R, x, z = sp.symbols("R x z", real=True)
    polynomials = [sp.legendre(k, 2*R-1).expand() for k in range(4)]
    coefficients = [sp.Rational(1, 4), sp.Rational(9, 20), sp.Rational(1, 4), sp.Rational(1, 20)]
    reconstruction = sp.expand(R**3-sum(c*p for c, p in zip(coefficients, polynomials)))
    eigen = [sp.expand(4*(R*(1-R)*sp.diff(p, R, 2)+(1-2*R)*sp.diff(p, R))
                       + 4*k*(k+1)*p) for k, p in enumerate(polynomials)]
    H6 = x**6-15*x**4*z**2+15*x**2*z**4-z**6
    harmonic = sp.diff(H6, x, 2)+sp.diff(H6, z, 2)
    powers = [sp.integrate((4*c*p)**2, (R, 0, 1)) for c, p in zip(coefficients, polynomials)]
    powers[-1] += 16*sp.integrate(R**6/2, (R, 0, 1))
    # Lower scalar metric norm > upper induced TT metric norm, after
    # canceling the common positive factor kappa s^2/V (s != 0).
    squared_gap = sp.Rational(66, 25)**2-sp.Rational(16*6, 49)
    return {"powers_times_volume": [str(p) for p in powers],
            "numeric_powers_times_volume": [float(p) for p in powers],
            "Legendre_reconstruction": str(reconstruction),
            "Legendre_eigen_residuals": [str(r) for r in eigen],
            "degree_six_Euclidean_laplacian": str(sp.expand(harmonic)),
            "squared_metric_bound_gap": str(squared_gap),
            "exact": bool(reconstruction == 0 and all(r == 0 for r in eigen)
                          and sp.expand(harmonic) == 0 and squared_gap > 0)}


VERDICT_FIELDS = ("conditional_constraints", "standing_wave_transverse_forcing",
                  "size_comparison", "complete_scalar_backreaction_bound", "triangle_map", "readout")


def verdict(checks):
    if not checks or not all(checks.values()):
        return {**dict.fromkeys(VERDICT_FIELDS, "UNRESOLVED"),
                "failed_checks": [k for k, v in checks.items() if not v]}
    return dict(zip(VERDICT_FIELDS, (
        "LINEAR_CMC_PARTICULAR_RESPONSE_BOUNDED_WITH_ZERO_SUPPORT_PERTURBATIONS",
        "ZERO_FOR_LEADING_STANDING_WAVE_ONLY", "CONTINUOUS_METRIC_BOUNDS_REPORTED",
        "NOT_ESTABLISHED_WITHOUT_SUPPORT_AND_EVOLUTION_CLOSURE", "NOT_DERIVED", "NOT_DERIVED")))
