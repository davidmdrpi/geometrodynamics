"""Pointwise support by one real conformal scalar; freeze de55f3f.

The odd-sector exclusion is a global analytic argument, reproduced in the
derivation. Numerical witnesses are controls, not a classification by search.
The evolved response is about an EVEN homogeneous control, not BAM support.
"""

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.integrate import solve_ivp


BASELINE = "6631aa3057b16ee96c2a35f45fc43a848d65fca4"
PUBLIC_PREREG = "de55f3f3175adafcf2cb760e0aef5767ca3e5016"
SEED = 2026090915
REQUIRED_CHECKS = (
    "stress_identity", "global_obstruction", "homogeneous_background",
    "admissibility", "response_variation", "constraint_propagation",
    "nonfluid_control", "scope_and_order", "fail_closed",
)
VERDICT_FIELDS = (
    "odd_sector_exact_ESU_support", "homogeneous_even_control",
    "control_kinetic_regularity", "control_response", "control_constraints",
    "BAM_support_selection", "triangle_map", "Phi_selection", "causality_gate",
)


def stf(matrix):
    matrix = np.asarray(matrix)
    return matrix - np.trace(matrix, axis1=-2, axis2=-1)[..., None, None]*np.eye(3)/3


@dataclass(frozen=True)
class HomogeneousSupport:
    radius: float = 1.
    kappa: float = 1.
    phase: float = 0.

    def __post_init__(self):
        if not all(np.isfinite([self.radius, self.kappa, self.phase])):
            raise ValueError("finite parameters required")
        if self.radius <= 0 or self.kappa <= 0:
            raise ValueError("positive radius and kappa required")

    @property
    def amplitude(self):
        return math.sqrt(3/self.kappa)

    @property
    def cosmological_constant(self):
        return 3/(2*self.radius**2)

    @property
    def density(self):
        return 3/(2*self.kappa*self.radius**2)

    def jets(self, time):
        theta = np.asarray(time)/self.radius + self.phase
        phi = self.amplitude*np.cos(theta)
        velocity = -self.amplitude*np.sin(theta)/self.radius
        return phi, velocity, -phi/self.radius**2

    def kinetic(self, time):
        phi = self.jets(time)[0]
        F = 1/self.kappa-phi**2/6
        K = 1/(self.kappa*F)+3/(2*self.kappa)*(phi/(3*F))**2
        return F, K


def improved_stress(metric, inverse, einstein, phi, derivative, hessian):
    """Off-shell stress, including xi G phi^2; no wave substitution."""
    square = derivative @ inverse @ derivative
    box = np.sum(inverse*hessian)
    return (np.outer(derivative, derivative)-metric*square/2
            +(metric*2*(square+phi*box)
              -2*(np.outer(derivative, derivative)+phi*hessian)
              +einstein*phi**2)/6)


def curvature_from_jets(metric, dg, ddg):
    """Levi-Civita/Ricci directly from metric jets; dg[k,i,j]=partial_k g_ij.

    Independent of the linear scalar Einstein equations used for evolution.
    """
    inverse = np.linalg.inv(metric)
    dinverse = -np.einsum("ia,kab,bj->kij", inverse, dg, inverse)
    gamma = np.zeros((4, 4, 4))
    dgamma = np.zeros((4, 4, 4, 4))
    for u in range(4):
        for i in range(4):
            for j in range(4):
                for v in range(4):
                    first = dg[i, v, j]+dg[j, v, i]-dg[v, i, j]
                    gamma[u, i, j] += inverse[u, v]*first/2
                    for k in range(4):
                        second = ddg[k, i, v, j]+ddg[k, j, v, i]-ddg[k, v, i, j]
                        dgamma[k, u, i, j] += (dinverse[k, u, v]*first
                                              +inverse[u, v]*second)/2
    ricci = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            for k in range(4):
                ricci[i, j] += dgamma[k, k, i, j]-dgamma[j, k, i, k]
                for v in range(4):
                    ricci[i, j] += (gamma[k, i, j]*gamma[v, k, v]
                                    -gamma[v, i, k]*gamma[k, j, v])
    scalar = np.sum(inverse*ricci)
    return inverse, gamma, ricci, ricci-metric*scalar/2, scalar


def zonal(degree, angle):
    """Y_l(chi)=sin((l+1)chi)/sin chi; exact unit-S3 eigenfunction jets."""
    if degree < 0 or degree != int(degree):
        raise ValueError("nonnegative integer degree required")
    k = degree+1
    value = math.sin(k*angle)/math.sin(angle)
    first = k*math.cos(k*angle)/math.sin(angle)-value/math.tan(angle)
    second = -degree*(degree+2)*value-2*first/math.tan(angle)
    return value, first, second


def direct_geometry(model, time, degree, field_jets, lapse_jets, psi_jets,
                    epsilon=0., angle=.83, polar=1.07):
    """Exact local metric/field jets for exponentiated longitudinal metric.

    alpha, psi and chi are harmonic coefficients; each input is (f,dt f,dtt f).
    Exponentiation fixes a finite-epsilon continuation, with the required
    linear tangent. No Einstein or scalar equation is imposed on input jets.
    """
    a = model.radius
    Y, DY, DDY = zonal(degree, angle)
    alpha, ad, add = lapse_jets
    psi, pd, pdd = psi_jets
    chi, cd, cdd = field_jets
    logs = np.array([2*epsilon*alpha*Y, 2*math.log(a)-2*epsilon*psi*Y,
                     2*math.log(a*math.sin(angle))-2*epsilon*psi*Y,
                     2*math.log(a*math.sin(angle)*math.sin(polar))-2*epsilon*psi*Y])
    diag = np.exp(logs)*np.array([-1., 1., 1., 1.])
    dl = np.zeros((4, 4))
    ddl = np.zeros((4, 4, 4))
    dl[0] = 2*epsilon*Y*np.array([ad, -pd, -pd, -pd])
    dl[1] = 2*epsilon*DY*np.array([alpha, -psi, -psi, -psi])
    dl[1, 2:] += 2/math.tan(angle)
    dl[2, 3] = 2/math.tan(polar)
    ddl[0, 0] = 2*epsilon*Y*np.array([add, -pdd, -pdd, -pdd])
    ddl[0, 1] = ddl[1, 0] = 2*epsilon*DY*np.array([ad, -pd, -pd, -pd])
    ddl[1, 1] = 2*epsilon*DDY*np.array([alpha, -psi, -psi, -psi])
    ddl[1, 1, 2:] -= 2/math.sin(angle)**2
    ddl[2, 2, 3] = -2/math.sin(polar)**2
    g = np.diag(diag)
    dg = np.array([np.diag(diag*v) for v in dl])
    ddg = np.array([[np.diag(diag*(ddl[i, j]+dl[i]*dl[j]))
                     for j in range(4)] for i in range(4)])
    inverse, gamma, ricci, G, R = curvature_from_jets(g, dg, ddg)
    P, Pd, Pdd = model.jets(time)
    value = P+epsilon*chi*Y
    derivative = np.array([Pd+epsilon*cd*Y, epsilon*chi*DY, 0., 0.])
    dd = np.zeros((4, 4))
    dd[0, 0] = Pdd+epsilon*cdd*Y
    dd[0, 1] = dd[1, 0] = epsilon*cd*DY
    dd[1, 1] = epsilon*chi*DDY
    hessian = dd-np.einsum("kij,k->ij", gamma, derivative)
    stress = improved_stress(g, inverse, G, value, derivative, hessian)
    frame = 1/np.sqrt(np.abs(diag))
    physical = stress*np.outer(frame, frame)
    geometric = G*np.outer(frame, frame)
    wave = np.sum(inverse*hessian)-R*value/6
    return {"stress": physical, "einstein": geometric, "R": R,
            "trace": np.sum(inverse*stress), "wave": wave, "phi": value,
            "ricci": ricci, "metric": g}


def stress_response(model, time, degree, field_jets, lapse_jets, psi_jets):
    """All four physical scalar stress coefficients, off shell, one harmonic."""
    a = model.radius
    lam = degree*(degree+2)/a**2
    P, Pd, _ = model.jets(time)
    c, cd, cdd = field_jets
    al, ald, _ = lapse_jets
    ps, psd, psdd = psi_jets
    rho = (Pd*cd+P*c/a**2+P*lam*c/3-al*Pd**2
           -P**2*lam*ps/3+P**2*ps/a**2-P*Pd*psd)
    J = (P*cd-2*Pd*c-P**2*psd-P*Pd*al)/3
    p = (Pd*cd/3-P*cdd/3-al*Pd**2/3-2*al*P**2/(3*a**2)
         +P*Pd*ald/3+2*P*Pd*psd/3-2*P*lam*c/9
         +P**2*psdd/3-P**2*ps/(3*a**2)+P**2*lam*(ps-al)/9)
    Pi = -P*c/3+P**2*(ps-al)/6
    delta_R = -6*psdd+12*ps/a**2-4*lam*ps+2*lam*al
    KG = -cdd-(lam+1/a**2)*c-2*al*P/a**2+(ald+3*psd)*Pd-delta_R*P/6
    return dict(rho=rho, J=J, p=p, Pi=Pi, KG=KG, delta_R=delta_R)


def response_tensor(model, degree, response, angle=.83):
    """Reconstruct orthonormal response for the zonal harmonic."""
    Y, DY, DDY = zonal(degree, angle)
    hess = np.diag([DDY, DY/math.tan(angle), DY/math.tan(angle)])/model.radius**2
    result = np.zeros((4, 4))
    result[0, 0] = response["rho"]*Y
    result[0, 1] = result[1, 0] = -response["J"]*DY/model.radius
    result[1:, 1:] = np.eye(3)*response["p"]*Y+response["Pi"]*stf(hess)
    return result


def lapse_from_anisotropy(model, time, state):
    c, cd, ps, psd = state
    P, Pd, _ = model.jets(time)
    f = 1-model.kappa*P**2/6
    fd = -model.kappa*P*Pd/3
    factor = model.kappa*P/(3*f)
    derivative = model.kappa*(Pd/f-P*fd/f**2)/3
    return ps+factor*c, psd+factor*cd+derivative*c


def response_residuals(model, time, degree, state, accelerations=(0., 0.)):
    c, cd, ps, psd = state
    cdd, psdd = accelerations
    al, ald = lapse_from_anisotropy(model, time, state)
    r = stress_response(model, time, degree, (c, cd, cdd), (al, ald, 0.), (ps, psd, psdd))
    lam = degree*(degree+2)/model.radius**2
    return dict(
        hamiltonian=2*(3/model.radius**2-lam)*ps-model.kappa*r["rho"],
        momentum=-2*psd-model.kappa*r["J"],
        spatial=2*psdd-2*ps/model.radius**2+2*lam*(ps-al)/3-model.kappa*r["p"],
        anisotropic=ps-al-model.kappa*r["Pi"], KG=r["KG"],
        trace=-r["rho"]+3*r["p"]-model.jets(time)[0]*r["KG"],
        Pi=r["Pi"], rho=r["rho"], p=r["p"], alpha=al,
    )


def initial_response(model, degree, chi=1., velocity=.2, time=0.):
    if degree < 2 or degree != int(degree):
        raise ValueError("this scalar response gate covers integer degree >=2")
    def constraints(metric):
        r = response_residuals(model, time, degree, (chi, velocity, *metric))
        return np.array([r["hamiltonian"], r["momentum"]])
    zero = constraints((0., 0.))
    matrix = np.column_stack([constraints(e)-zero for e in np.eye(2)])
    metric = np.linalg.solve(matrix, -zero)
    return np.array([chi, velocity, *metric])


def response_rhs(model, degree, time, state):
    def equations(acceleration):
        r = response_residuals(model, time, degree, state, acceleration)
        return np.array([r["KG"], r["spatial"]])
    zero = equations((0., 0.))
    matrix = np.column_stack([equations(e)-zero for e in np.eye(2)])
    cdd, psdd = np.linalg.solve(matrix, -zero)
    return np.array([state[1], cdd, state[3], psdd])


def integrate_response(model, degree, times, initial=None, rtol=1e-10, atol=1e-12):
    times = np.asarray(times, dtype=float)
    if degree < 2 or degree != int(degree):
        raise ValueError("integer degree >=2 required")
    if times.ndim != 1 or len(times) < 2 or not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0):
        raise ValueError("finite increasing times required")
    y0 = initial_response(model, degree, time=times[0]) if initial is None else np.asarray(initial)
    sol = solve_ivp(lambda t, y: response_rhs(model, degree, t, y), (times[0], times[-1]),
                    y0, t_eval=times, method="DOP853", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y.T


@lru_cache(maxsize=1)
def exact_certificate():
    """Symbolic local identities. Global extension is proved in the document."""
    import sympy as s
    a, k, t = s.symbols("a k t", positive=True, real=True)
    P = s.sqrt(3/k)*s.cos(t/a)
    Pd, Pdd = s.diff(P, t), s.diff(P, t, 2)
    rho = (Pd**2+P**2/a**2)/2
    pressure = Pd**2/2+(-2*(Pd**2+P*Pdd)-P**2/a**2)/6
    F = 1/k-P**2/6
    K = 1/(k*F)+3/(2*k)*(P/(3*F))**2
    residuals = {
        "wave": Pdd+P/a**2,
        "hamiltonian": 3/a**2-3/(2*a**2)-k*rho,
        "spatial": -1/a**2+3/(2*a**2)-k*pressure,
        "radiation": rho-3*pressure,
        "kinetic_simplification": K-1/(1-k*P**2/6)**2,
    }
    # Independent symbolic isotropy derivation for generic first/second jets.
    v = s.Matrix(s.symbols("v0:3"))
    h00,h01,h02,h11,h12,h22,z = s.symbols("h00 h01 h02 h11 h12 h22 z", nonzero=True)
    H = s.Matrix([[h00,h01,h02],[h01,h11,h12],[h02,h12,h22]])
    tf = lambda M: M-s.trace(M)*s.eye(3)/3
    spatial = tf(v*v.T-s.Rational(1,6)*(2*v*v.T+2*z*H))
    reciprocal_hessian = 2*v*v.T/z**3-H/z**2
    identity = spatial-z**3*tf(reciprocal_hessian)/3
    residuals.update({f"isotropy_{i}{j}": identity[i,j] for i in range(3) for j in range(3)})
    # h=A+B.x, Hess h=-(B.x)g/a^2: phi=1/h has isotropic stress.
    A,bx = s.symbols("A bx")
    hh = A+bx
    grad_phi = -v/hh**2
    hess_phi = 2*v*v.T/hh**3+bx*s.eye(3)/(a**2*hh**2)
    reciprocal = tf(2*grad_phi*grad_phi.T/3-hess_phi/(3*hh))
    residuals.update({f"reciprocal_{i}{j}": reciprocal[i,j] for i in range(3) for j in range(3)})
    # div Hess h: grad lambda = 3 grad lambda + 2/a^2 grad h.
    dh = s.symbols("dh")
    dlambda = s.solve(s.Symbol("dl")-3*s.Symbol("dl")-2*dh/a**2, s.Symbol("dl"))[0]
    residuals["curvature_commutation"] = dlambda+dh/a**2
    # The response accelerations remain solvable at both types of turning point.
    pp = s.symbols("P", real=True)
    f = 1-k*pp**2/6
    acceleration_matrix = s.Matrix([[-1, pp], [k*pp/3, 2*f]])
    residuals["acceleration_determinant"] = acceleration_matrix.det()+2
    u, ell = s.symbols("u ell", nonnegative=True)
    determinant = 4*((1+u)/2)**2*(ell-3)-6*u
    residuals["initial_constraint_determinant"] = determinant-((ell-3)+(2*ell-12)*u+(ell-3)*u**2)
    reduced = {name: str(s.simplify(value)) for name,value in residuals.items()}
    return dict(residuals=reduced, all_zero=all(v=="0" for v in reduced.values()),
                F_lower_bound="1/(2*kappa)", K_range="[1,4]",
                acceleration_determinant="-2",
                initial_constraint_determinant_times_a2="(L-3)+(2L-12)sin^2(theta)+(L-3)sin^4(theta); L=l(l+2)>=8",
                global_argument="docs/scalar_esu_support.md#the-global-step",
                global_scope="one smooth real scalar; exact pointwise isotropy on complete round S3")


def failed_checks(checks):
    return [name for name in REQUIRED_CHECKS if checks.get(name) is not True]


def verdict(checks):
    failed = failed_checks(checks)
    if failed:
        return {**{name: "UNRESOLVED" for name in VERDICT_FIELDS}, "failed_checks": failed}
    return dict(
        odd_sector_exact_ESU_support="EXCLUDED_IN_STATED_SINGLE_REAL_SCALAR_CLASS",
        homogeneous_even_control="EXACT_BUT_OUTSIDE_ODD_SECTOR",
        control_kinetic_regularity="REGULAR_POSITIVE_KINETIC_COEFFICIENTS",
        control_response="GENERIC_SCALAR_ANISOTROPIC_STRESS",
        control_constraints="PROPAGATED_FOR_TESTED_L_GE_2_RESPONSES",
        BAM_support_selection="NOT_DERIVED", triangle_map="NOT_DERIVED",
        Phi_selection="NOT_DERIVED", causality_gate="OPEN",
    )
