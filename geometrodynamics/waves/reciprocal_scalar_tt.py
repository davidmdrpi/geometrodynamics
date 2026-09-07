"""Reciprocal conformal scalar / homogeneous TT action on the ESU.

Public freeze d8dc90d. This is a variational projection with all five tensor
components and a complete scalar harmonic multiplet. It does not solve the
omitted Einstein constraints or identify a two-boundary triangle history.
"""

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm, null_space

from geometrodynamics.bulk.tt_triangle_rotor import TensorModel


PUBLIC_PREREG = "d8dc90d6d66e14824c337f96ee93a512dc9ed84f"
SEED = 2026090711


def _stf_basis():
    out = [np.diag([1., -1., 0.]) / math.sqrt(2),
           np.diag([1., 1., -2.]) / math.sqrt(6)]
    for i, j in ((0, 1), (0, 2), (1, 2)):
        E = np.zeros((3, 3))
        E[i, j] = E[j, i] = 1 / math.sqrt(2)
        out.append(E)
    return np.array(out)


STF_BASIS = _stf_basis()


def tensor(b):
    return np.einsum("...a,aij->...ij", b, STF_BASIS)


def components(beta):
    return np.einsum("aij,...ij->...a", STF_BASIS, beta)


def powers(degree):
    return tuple((a, b, c, degree-a-b-c)
                 for a in range(degree, -1, -1)
                 for b in range(degree-a, -1, -1)
                 for c in range(degree-a-b, -1, -1))


@lru_cache(None)
def sphere_moment(exponents):
    """Exact monomial moment formula for normalized Haar on the unit S3."""
    if any(e % 2 for e in exponents):
        return 0.
    half = [e // 2 for e in exponents]
    num = math.prod(math.prod(range(1, 2*k, 2)) for k in half)
    den = math.prod(range(4, 4 + 2*sum(half), 2))
    return num / den


def monomials(points, exponents):
    return np.prod(np.asarray(points)[:, None, :] ** np.asarray(exponents)[None, :, :], axis=-1)


def cross_matrix(axis):
    x, y, z = axis
    return np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])


def quaternion_derivatives():
    """Matrices S_i with L_i(x)=x e_i=S_i x, in the inherited convention."""
    matrices = []
    for e in np.eye(3):
        S = np.zeros((4, 4))
        S[0, 1:] = -e
        S[1:, 0] = e
        S[1:, 1:] = -cross_matrix(e)
        matrices.append(S)
    return np.array(matrices)


QUATERNION_DERIVATIVES = quaternion_derivatives()


@lru_cache(None)
def harmonic_multiplet(degree):
    return HarmonicMultiplet(degree)


class HarmonicMultiplet:
    """Real harmonic polynomials; analytic moments fix the normalization.

    B maps modal coefficients into monomial coefficients. B is orthonormal
    for normalized unit-S3 Haar. Physical fields divide by sqrt(Vol(S3)),
    so their modal kinetic term is q_dot^2/2 at any radius.
    """

    def __init__(self, degree):
        if not isinstance(degree, (int, np.integer)) or not 1 <= degree <= 7:
            raise ValueError("degree must be an integer from 1 to 7")
        self.degree = int(degree)
        self.exponents = powers(degree)
        self.index = {e: i for i, e in enumerate(self.exponents)}
        self.moments = np.array([[sphere_moment(tuple(x+y for x, y in zip(a, b)))
                                 for b in self.exponents] for a in self.exponents])
        lower = powers(degree-2) if degree >= 2 else ()
        lower_index = {e: i for i, e in enumerate(lower)}
        lap = np.zeros((len(lower), len(self.exponents)))
        for j, e in enumerate(self.exponents):
            for axis, k in enumerate(e):
                if k >= 2:
                    target = list(e)
                    target[axis] -= 2
                    lap[lower_index[tuple(target)], j] += k*(k-1)
        H = null_space(lap) if len(lower) else np.eye(len(self.exponents))
        vals, vecs = np.linalg.eigh(H.T @ self.moments @ H)
        self.B = H @ (vecs / np.sqrt(vals)[None, :])
        self.laplacian = lap
        self.dimension = self.B.shape[1]
        self.polynomial_generators = np.array([self.polynomial_operator(S)
                                               for S in QUATERNION_DERIVATIVES])
        self.D = np.array([self.B.T @ self.moments @ M @ self.B
                           for M in self.polynomial_generators])
        rotations = []
        for axis in np.eye(3):
            S = np.zeros((4, 4))
            S[1:, 1:] = cross_matrix(axis)
            rotations.append(self.B.T @ self.moments @ self.polynomial_operator(S) @ self.B)
        self.rotation_generators = np.array(rotations)
        self.K = np.array([[-(Di @ Dj + Dj @ Di)/2 for Dj in self.D] for Di in self.D])

    def polynomial_operator(self, S):
        out = np.zeros((len(self.exponents), len(self.exponents)))
        for j, e in enumerate(self.exponents):
            for r, count in enumerate(e):
                if count:
                    for s in range(4):
                        if S[r, s]:
                            target = list(e)
                            target[r] -= 1
                            target[s] += 1
                            out[self.index[tuple(target)], j] += count*S[r, s]
        return out

    def coherent(self, axis):
        """Normalized Re[(x0 + i axis.x_vec)^n]; no mode projection is fitted."""
        axis = np.asarray(axis, dtype=float)
        if axis.shape != (3,) or not np.isfinite(axis).all() or abs(axis @ axis-1) > 1e-10:
            raise ValueError("axis must be a finite unit vector")
        v = np.r_[1.+0j, 1j*axis]
        poly = np.array([(math.factorial(self.degree) / math.prod(math.factorial(k) for k in e)
                          * np.prod(v ** np.array(e))).real for e in self.exponents])
        q = self.B.T @ self.moments @ poly
        return q / np.linalg.norm(q)

    def rotate(self, coefficients, axis, angle):
        axis = np.asarray(axis, dtype=float)
        axis = axis / np.linalg.norm(axis)
        K = np.einsum("i,ijk->jk", axis, self.rotation_generators)
        return expm(-angle*K) @ coefficients

    def algebra_checks(self):
        I = np.eye(self.dimension)
        return {
            "dimension": self.dimension,
            "expected_dimension": (self.degree+1)**2,
            "orthonormality": float(np.linalg.norm(self.B.T @ self.moments @ self.B-I)),
            "antisymmetry": float(max(np.linalg.norm(D+D.T) for D in self.D)),
            "casimir_scaled": float(np.linalg.norm(-sum(D @ D for D in self.D)
                                                   - self.degree*(self.degree+2)*I)
                                     / (self.degree*(self.degree+2)*np.linalg.norm(I))),
            "harmonic_residual_scaled": float(np.linalg.norm(self.laplacian @ self.B)
                                               / max(1., np.linalg.norm(self.B))),
            "generator_closure_scaled": float(max(
                np.linalg.norm(M @ self.B-self.B @ D) / max(1., np.linalg.norm(M @ self.B))
                for M, D in zip(self.polynomial_generators, self.D))),
        }

    def polynomial_jets(self, points, coefficients):
        coeff = self.B @ coefficients
        ex = np.array(self.exponents)
        value = monomials(points, ex) @ coeff
        grad = np.empty((len(points), 4))
        hess = np.empty((len(points), 4, 4))
        for i in range(4):
            exp_i = ex.copy()
            exp_i[:, i] = np.maximum(0, exp_i[:, i]-1)
            grad[:, i] = monomials(points, exp_i) @ (coeff*ex[:, i])
            for j in range(4):
                weights = ex[:, i] * (ex[:, j]-(i == j))
                exp_ij = exp_i.copy()
                exp_ij[:, j] = np.maximum(0, exp_ij[:, j]-1)
                hess[:, i, j] = monomials(points, exp_ij) @ (coeff*weights)
        return value, grad, hess


@dataclass(frozen=True)
class ReciprocalModel:
    degree: int = 3
    radius: float = 1.
    kappa: float = 1.

    def __post_init__(self):
        TensorModel(self.radius, self.kappa)
        harmonic_multiplet(self.degree)

    @property
    def multiplet(self):
        return harmonic_multiplet(self.degree)

    @property
    def C(self):
        return 2*math.pi**2*self.radius**3/self.kappa

    @property
    def volume(self):
        return 2*math.pi**2*self.radius**3

    @property
    def omega_tensor2(self):
        return 8/self.radius**2

    @property
    def omega_scalar2(self):
        return (self.degree+1)**2/self.radius**2

    @property
    def F(self):
        return np.einsum("aij,ijuv->auv", STF_BASIS, self.multiplet.K)/self.radius**2

    @property
    def size(self):
        return 2*(5+self.multiplet.dimension)

    def pack(self, b, P, q, p):
        return np.concatenate((b, P, q, p))

    def unpack(self, state):
        d = self.multiplet.dimension
        return state[..., :5], state[..., 5:10], state[..., 10:10+d], state[..., 10+d:]

    def source(self, q):
        return np.einsum("...i,aij,...j->...a", q, self.F, q)

    def energy_parts(self, state):
        b, P, q, p = self.unpack(np.asarray(state))
        return {
            "tensor": np.sum(P*P, axis=-1)/(2*self.C)
                      + self.C*self.omega_tensor2*np.sum(b*b, axis=-1)/2,
            "scalar": np.sum(p*p, axis=-1)/2 + self.omega_scalar2*np.sum(q*q, axis=-1)/2,
            "interaction": -np.sum(b*self.source(q), axis=-1),
        }

    def hamiltonian(self, state):
        return sum(self.energy_parts(state).values())

    def rhs(self, time, state, reciprocal=True):
        b, P, q, p = self.unpack(state)
        reaction = 2*np.einsum("a,aij,j->i", b, self.F, q) if reciprocal else np.zeros_like(q)
        return self.pack(P/self.C, -self.C*self.omega_tensor2*b+self.source(q),
                         p, -self.omega_scalar2*q+reaction)

    def integrate(self, state, times, rtol=1e-10, atol=1e-12, reciprocal=True):
        times = np.asarray(times)
        sol = solve_ivp(lambda t, y: self.rhs(t, y, reciprocal), (times[0], times[-1]),
                        state, t_eval=times, method="DOP853", rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.y.T

    def lagrangian(self, b, bdot, q, qdot):
        return (self.C*(bdot @ bdot-self.omega_tensor2*(b @ b))/2
                + (qdot @ qdot-self.omega_scalar2*(q @ q))/2+b @ self.source(q))

    def static_scalar_lagrangian(self, beta, q, qdot):
        """Exact scalar action on a STATIC homogeneous anisotropic metric.

        This is an independent first-variation control, not the full dynamic
        parent of the truncated ODE. The curvature term is retained.
        """
        eigenvalues = np.linalg.eigvalsh(beta)
        curvature = 2/self.radius**2 * (2*np.exp(-2*eigenvalues).sum()
                                        - np.exp(4*eigenvalues).sum())
        metric_inverse = expm(-2*beta)/self.radius**2
        gradient_energy = np.einsum("ij,ijuv,u,v->", metric_inverse, self.multiplet.K, q, q)
        return .5*(qdot @ qdot-gradient_energy-curvature*(q @ q)/6)


def primary_data(model=ReciprocalModel()):
    n, v = np.array([0., 0., 1.]), np.array([.4, 0., 0.])
    beta = .01*(np.outer(n, n)-np.eye(3)/3)
    velocity = .01*(np.outer(v, n)+np.outer(n, v))
    m = np.array([1., 2., 3.])/math.sqrt(14)
    q = .2*model.multiplet.coherent(m)
    return model.pack(components(beta), model.C*components(velocity), q, np.zeros_like(q))


def sphere_quadrature(radial_order=8, angular_order=16, radius=1.):
    """Hopf coordinates: u uniform on [0,1], two uniform circle angles."""
    z, w = np.polynomial.legendre.leggauss(radial_order)
    u, w = (z+1)/2, w/2
    angles = np.arange(angular_order)*2*math.pi/angular_order
    u, alpha, beta = np.meshgrid(u, angles, angles, indexing="ij")
    points = np.column_stack((np.sqrt(u).ravel()*np.cos(alpha).ravel(),
                              np.sqrt(u).ravel()*np.sin(alpha).ravel(),
                              np.sqrt(1-u).ravel()*np.cos(beta).ravel(),
                              np.sqrt(1-u).ravel()*np.sin(beta).ravel()))
    weights = np.repeat(w, angular_order**2)*2*math.pi**2*radius**3/angular_order**2
    return points, weights


def scalar_jets(model, q, p, qddot, points):
    """Pointwise jets in the repository's left-invariant orthonormal frame."""
    value, grad_e, hess_e = model.multiplet.polynomial_jets(points, q)
    dt, dtgrad_e, _ = model.multiplet.polynomial_jets(points, p)
    dtt = monomials(points, model.multiplet.exponents) @ model.multiplet.B @ qddot
    frames = np.einsum("iab,pb->pia", QUATERNION_DERIVATIVES, points)
    grad = np.einsum("pia,pa->pi", frames, grad_e)/model.radius
    dtgrad = np.einsum("pia,pa->pi", frames, dtgrad_e)/model.radius
    hess = (np.einsum("pia,pab,pjb->pij", frames, hess_e, frames)
            - model.degree*value[:, None, None]*np.eye(3))/model.radius**2
    norm = math.sqrt(model.volume)
    jets = {"phi": value[:, None]/norm, "dt": dt[:, None]/norm,
            "dtt": dtt[:, None]/norm, "grad": grad[:, None, :]/norm,
            "dtgrad": dtgrad[:, None, :]/norm, "hess": hess[:, None, :, :]/norm,
            "laplacian": np.trace(hess, axis1=1, axis2=2)[:, None]/norm}
    return jets, frames


def inherited_stress_diagnostics(model, state, radial_order=8, angular_order=16):
    """Independent improved-stress integral and omitted constraint sources.

    The inherited stress routine hardcodes the unit ESU. Values reported for
    the coupled history use that round-background functional to the retained
    leading order in phi^2; omitted metric corrections are not inferred zero.
    """
    if model.radius != 1.:
        raise ValueError("the inherited stress diagnostic requires radius=1")
    from .backreaction import shear_projection, stress_series
    points, weights = sphere_quadrature(radial_order, angular_order)
    b, P, q, p = model.unpack(state)
    qddot = model.unpack(model.rhs(0., state))[3]
    jets, frames = scalar_jets(model, q, p, qddot, points)
    stress = stress_series(jets)
    projected = shear_projection(stress, points, frames)
    integrated = np.einsum("p,ptij->ij", weights, projected)
    rho = stress[:, 0, 0, 0]
    flux = stress[:, 0, 0, 1:]
    mean = float(weights @ rho/model.volume)
    return {
        "integrated_tt": integrated,
        "action_tt": tensor(model.source(q)),
        "source_scaled_error": float(np.linalg.norm(integrated-tensor(model.source(q)))
                                      / max(1., np.linalg.norm(integrated))),
        "mean_energy_density": mean,
        "inhomogeneous_energy_rms": float(np.sqrt(weights @ (rho-mean)**2/model.volume)),
        "momentum_density_rms": float(np.sqrt(weights @ np.sum(flux*flux, axis=1)/model.volume)),
        "constraint_scope": "round-background matter sources at order phi^2; "
                            "homogeneous linear TT has delta G_00=delta G_0i=0",
        "points": len(points),
    }


VERDICT_FIELDS = ("reciprocal_dynamics", "scalar_modal_closure", "einstein_constraints",
                  "triangle_history_map", "source_local_readout", "probability_selection")


def verdict(checks):
    failed = sorted(k for k, v in checks.items() if not v) if checks else ["no checks supplied"]
    if failed:
        return {**dict.fromkeys(VERDICT_FIELDS, "UNRESOLVED"), "failed_checks": failed}
    return {"reciprocal_dynamics": "VARIATIONAL_TT_SCALAR_PROJECTION_VERIFIED",
            "scalar_modal_closure": "COMPLETE_MULTIPLETS_INVARIANT_IN_PROJECTED_SCALAR_EQUATION",
            "einstein_constraints": "OMITTED_METRIC_AND_SUPPORT_RESPONSE_REQUIRED",
            "triangle_history_map": "NOT_DERIVED",
            "source_local_readout": "NOT_DERIVED",
            "probability_selection": "NOT_DERIVED", "failed_checks": []}
