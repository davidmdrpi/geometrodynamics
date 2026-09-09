"""Leading-order field-supported homogeneous TT rotor; freeze 11e625f.

The scalar is free at O(s), the metric projection is O(s^2). This is not
a full corrected Einstein-matter history or a selection of its initial data.
"""

from dataclasses import dataclass
from functools import cached_property, lru_cache
import math

import numpy as np

from geometrodynamics.bulk import tt_triangle_rotor as rotor
from . import reciprocal_scalar_tt as rt


BASELINE = "c08f46a56def92edc3ab1ad1146973ab90c40aa2"
PUBLIC_PREREG = "11e625fa877de2bb0669af06384347ede69d5777"
SEED = 2026090814


@lru_cache(None)
def scalar_pair():
    """Fixed polynomials, converted to the inherited basis without a fit."""
    h = rt.harmonic_multiplet(3)
    polynomials = np.zeros((2, len(h.exponents)))
    for row, terms in enumerate((
            {(3, 0, 0, 0): 1, (1, 0, 0, 2): -3},
            {(0, 3, 0, 0): 1, (0, 1, 2, 0): -3})):
        for exponent, value in terms.items():
            polynomials[row, h.index[exponent]] = math.sqrt(8)*value
    modes = polynomials @ h.moments @ h.B
    modes.flags.writeable = False
    return modes


@dataclass(frozen=True)
class DrivenRotor:
    amplitude: float = .02  # scalar amplitude s; tensor amplitude is signed
    radius: float = 1.
    kappa: float = 1.

    def __post_init__(self):
        rt.ReciprocalModel(radius=self.radius, kappa=self.kappa)
        if not np.isfinite(self.amplitude) or self.amplitude <= 0:
            raise ValueError("positive finite scalar amplitude required")

    @cached_property
    def model(self):
        return rt.ReciprocalModel(radius=self.radius, kappa=self.kappa)

    @property
    def omega_scalar(self):
        return math.sqrt(self.model.omega_scalar2)

    @property
    def omega(self):
        return math.sqrt(self.model.omega_tensor2)/2

    @property
    def k(self):
        return 6*self.amplitude**2/(self.model.C*self.radius**2)

    @property
    def A(self):
        return -2*self.k/self.model.omega_tensor2

    @property
    def period(self):
        return math.pi/self.omega  # RP2 director/tensor period, not scalar period

    def scalar(self, time, second_amplitude=1.):
        f, g = scalar_pair()
        c, s = np.cos(self.omega_scalar*time), np.sin(self.omega_scalar*time)
        return (self.amplitude*(c*f+second_amplitude*s*g),
                self.amplitude*self.omega_scalar*(-s*f+second_amplitude*c*g))

    def source(self, time, second_amplitude=1.):
        q, _ = self.scalar(time, second_amplitude)
        return rt.tensor(self.model.source(q))/self.model.C

    def orbit(self, time, speed_factor=1.):
        speed = self.omega*speed_factor
        n = np.array([np.cos(speed*time), np.sin(speed*time), 0.])
        v = speed*np.array([-n[1], n[0], 0.])
        beta = rotor.embedding(self.A, n)
        velocity = rotor.field_velocity(self.A, 0., n, v)
        acceleration = rotor.field_acceleration(self.A, 0., 0., n, v, -speed**2*n)
        return {"n": n, "v": v, "beta": beta, "velocity": velocity,
                "acceleration": acceleration,
                "required_source": acceleration+self.model.omega_tensor2*beta}

    def initial(self):
        z = self.orbit(0.)
        q, p = self.scalar(0.)
        return self.model.pack(rt.components(z["beta"]),
                               self.model.C*rt.components(z["velocity"]), q, p)

    def integrate(self, times, rtol=1e-12, atol=1e-14):
        # Inherited 42-component evolution, including all five tensor modes
        # and the complete free 16-mode scalar. Never restrict to the cone.
        return self.model.integrate(self.initial(), times, reciprocal=False,
                                    rtol=rtol, atol=atol)


def full_stress(field, momentum, points):
    """Unit-radius inherited improved stress for arbitrary supplied degrees."""
    from .backreaction import stress_series
    jets = None
    for degree in sorted(set(field) | set(momentum)):
        model = rt.ReciprocalModel(degree=degree)
        zero = np.zeros(model.multiplet.dimension)
        q, p = field.get(degree, zero), momentum.get(degree, zero)
        part, frames = rt.scalar_jets(model, q, p, -model.omega_scalar2*q, points)
        jets = part if jets is None else {key: jets[key]+part[key] for key in jets}
    if jets is None:
        raise ValueError("nonempty field or momentum degrees required")
    return stress_series(jets)[:, 0], frames


def constraint_charges(stress, points, weights, frames):
    """Four Hamiltonian dipoles; six rotations AND four gradient dipoles."""
    rho, current = stress[:, 0, 0], -stress[:, 0, 1:]
    ambient = np.einsum("pi,pia->pa", current, frames)
    rotations = [weights @ (points[:, a]*ambient[:, b]-points[:, b]*ambient[:, a])
                 for a in range(4) for b in range(a+1, 4)]
    return {"Hamiltonian": np.einsum("p,p,pa->a", weights, rho, points),
            "rotations": np.array(rotations),
            "gradient_dipoles": np.einsum("p,pa->a", weights, ambient)}


def compatible(charges, scale=1., tolerance=1e-9):
    if not np.isfinite(scale) or scale <= 0:
        return False
    for key, size in (("Hamiltonian", 4), ("rotations", 6), ("gradient_dipoles", 4)):
        value = np.asarray(charges.get(key, []))
        if value.shape != (size,) or not np.isfinite(value).all():
            return False
        if np.max(np.abs(value))/scale >= tolerance:
            return False
    return True


def spatial_checks(candidate, radial_order=8, angular_order=16):
    if candidate.radius != 1.:
        raise ValueError("inherited full stress fixes unit radius")
    points, weights = rt.sphere_quadrature(radial_order, angular_order)
    rows = []
    for phase in (0., math.pi/8, math.pi/4, math.pi/2, math.pi):
        t = phase/candidate.omega_scalar
        q, p = candidate.scalar(t)
        stress, frames = full_stress({3: q}, {3: p}, points)
        integrated = rotor.stf(np.einsum("p,pij->ij", weights, stress[:, 1:, 1:]))
        charges = constraint_charges(stress, points, weights, frames)
        rows.append({"time": t, "points": len(points),
                     "source_error": float(np.linalg.norm((integrated/candidate.model.C-candidate.source(t))/candidate.k)),
                     "charges_over_s2": {key: (value/candidate.amplitude**2).tolist() for key, value in charges.items()},
                     "compatible": compatible(charges, candidate.amplitude**2)})
    return rows


def negative_constraint_controls():
    points, weights = rt.sphere_quadrature(8, 16)
    rng = np.random.default_rng(SEED)
    h = rt.harmonic_multiplet(2)
    q = rng.normal(size=h.dimension); q /= np.linalg.norm(q)
    left = np.column_stack([D @ q for D in h.D])
    p = rng.normal(size=h.dimension)
    p -= left @ np.linalg.lstsq(left, p, rcond=None)[0]
    p /= np.linalg.norm(p)
    stress, frames = full_stress({2: q}, {2: p}, points)
    rotations = constraint_charges(stress, points, weights, frames)
    invariant = np.array([-p @ D @ q for D in h.D])
    q1 = rng.normal(size=4); q1 /= np.linalg.norm(q1)
    p2 = rng.normal(size=9); p2 /= np.linalg.norm(p2)
    stress, frames = full_stress({1: q1}, {2: p2}, points)
    gradients = constraint_charges(stress, points, weights, frames)
    return {"missed_rotations": {"three_invariant_charges": invariant.tolist(),
                "charges": {k: v.tolist() for k, v in rotations.items()},
                "rejected": not compatible(rotations)},
            "missed_gradients": {"charges": {k: v.tolist() for k, v in gradients.items()},
                "rejected": not compatible(gradients)}}


@lru_cache(None)
def exact_certificate():
    """Exact polynomials/moments, independent of numerical mode eigenvectors."""
    import sympy as sp
    x = sp.symbols("x0:4", real=True)
    c, d = sp.symbols("c d", real=True)
    f = sp.sqrt(8)*(x[0]**3-3*x[0]*x[3]**2)
    g = sp.sqrt(8)*(x[1]**3-3*x[1]*x[2]**2)

    def mean(expression):
        total = sp.S(0)
        for exponents, coefficient in sp.Poly(sp.expand(expression), *x).terms():
            if any(e % 2 for e in exponents):
                continue
            ks = [e//2 for e in exponents]
            moment = sp.prod(sp.factorial(2*k)/(4**k*sp.factorial(k)) for k in ks)/sp.factorial(sum(ks)+1)
            total += coefficient*moment
        return sp.simplify(total)

    def derivative(poly, matrix):
        return sp.expand(sum(int(matrix[i, j])*x[j]*sp.diff(poly, x[i])
                             for i in range(4) for j in range(4)))

    Df = [derivative(f, S) for S in rt.QUATERNION_DERIVATIVES]
    Dg = [derivative(g, S) for S in rt.QUATERNION_DERIVATIVES]
    Gf = sp.Matrix(3, 3, lambda i, j: mean(Df[i]*Df[j]))
    Gg = sp.Matrix(3, 3, lambda i, j: mean(Dg[i]*Dg[j]))
    cross = sp.Matrix(3, 3, lambda i, j: mean((Df[i]*Dg[j]+Df[j]*Dg[i])/2))
    wave = [sp.simplify(sum(derivative(derivative(h, S), S)
                            for S in rt.QUATERNION_DERIVATIVES)+15*h) for h in (f, g)]
    q, p = c*f+d*g, 4*(-d*f+c*g)
    gradq = [sp.diff(q, xx)-3*xx*q for xx in x]
    gradp = [sp.diff(p, xx)-3*xx*p for xx in x]
    j = [(-2*p*dq+q*dp)/3 for dq, dp in zip(gradq, gradp)]
    rho = p**2/2+sum(dq*dq for dq in gradq)/6+sp.Rational(11, 2)*q*q
    ham = [mean(rho*xx) for xx in x]
    rotations = [mean(x[a]*j[b]-x[b]*j[a]) for a in range(4) for b in range(a+1, 4)]
    gradients = [mean(jj) for jj in j]
    theta = sp.Symbol("theta", real=True)
    n = sp.Matrix([sp.cos(theta), sp.sin(theta), 0])
    beta = -sp.Rational(3, 2)*(n*n.T-sp.eye(3)/3)
    source = sp.diag(-2, -2, 4)
    equation = (2*beta.diff(theta, 2)+8*beta-source).applyfunc(sp.trigsimp)
    expected = sp.diag(3, 3, 9)
    residuals = [mean(f*f)-1, mean(g*g)-1, mean(f*g), *wave,
                 *(Gf-expected), *(Gg-expected), *cross, *ham, *rotations, *gradients, *equation]
    return {"norms": [str(mean(f*f)), str(mean(g*g))], "overlap": str(mean(f*g)),
            "gradient_f": [[str(v) for v in row] for row in Gf.tolist()],
            "gradient_g": [[str(v) for v in row] for row in Gg.tolist()],
            "symmetric_cross_gradient": [[str(v) for v in row] for row in cross.tolist()],
            "wave_residuals": [str(v) for v in wave],
            "Hamiltonian_dipoles": [str(v) for v in ham],
            "six_rotation_charges": [str(v) for v in rotations],
            "four_gradient_charges": [str(v) for v in gradients],
            "tensor_equation_residuals": [str(v) for v in equation],
            "A_times_C_over_s2": "-3/2", "Omega_times_a": "sqrt(2)",
            "all_exact_zero": all(v == 0 for v in residuals)}


@lru_cache(None)
def frequency_certificate():
    import sympy as sp
    z = sp.Symbol("z", real=True)
    a, b, c = sp.symbols("a b c", real=True)
    quadratic = a*sp.cos(z)**2+2*b*sp.cos(z)*sp.sin(z)+c*sp.sin(z)**2
    fourier = (a+c)/2+(a-c)*sp.cos(2*z)/2+b*sp.sin(2*z)
    residual = sp.trigsimp(sp.expand_trig(quadratic-fourier))
    # Frequency independence is the analytic argument in the freeze/write-up;
    # these algebraic identities check its coefficients and exceptional branch.
    slow = sp.Rational(2, 5)
    return {"quadratic_Fourier_residual": str(residual),
            "allowed_Omega_times_a": ["4", "sqrt(2)"],
            "constant_source_exception_coefficient": str(8-4*sp.sqrt(2)**2),
            "original_speed_oscillatory_coefficient": str(8-4*slow**2),
            "original_speed_excluded": bool(slow != 4 and 8-4*slow**2 != 0),
            "all_exact": bool(residual == 0 and 8-4*sp.sqrt(2)**2 == 0),
            "scope": "necessary frequencies for nonzero constant-A uniform rotation, not a classification of all trajectories"}


REQUIRED_CHECKS = ("exact_polynomial_certificate", "action_and_improved_stress",
                   "complete_constraint_compatibility", "omitted_charge_controls",
                   "all_tensor_equations", "unrestricted_tensor_evolution",
                   "necessity_controls", "amplitude_radius_and_smallness",
                   "frequency_restriction", "fail_closed_verdicts")
VERDICT_FIELDS = ("normal_balance", "constraint_compatibility", "homogeneous_tensor_evolution",
                  "full_Einstein_matter_evolution", "preparation_selection", "triangle_map", "Phi_selection")


def failed_checks(checks):
    return (["missing: "+name for name in REQUIRED_CHECKS if name not in checks]
            + [name for name, passed in checks.items() if not passed])


def verdict(checks):
    failures = failed_checks(checks)
    if failures:
        return {**dict.fromkeys(VERDICT_FIELDS, "UNRESOLVED"), "failed_checks": failures}
    return {**dict(zip(VERDICT_FIELDS, (
        "LEADING_ORDER_FIELD_SUPPORTED_ROTOR", "CHOSEN_CONSTRAINT_COMPATIBLE_PREPARATION",
        "VERIFIED_AT_ORDER_S2", "NOT_ESTABLISHED", "NOT_DERIVED", "NOT_DERIVED", "NOT_DERIVED"))),
        "failed_checks": []}


def verdict_controls():
    good = dict.fromkeys(REQUIRED_CHECKS, True)
    if verdict(good)["normal_balance"] != "LEADING_ORDER_FIELD_SUPPORTED_ROTOR":
        return False
    for name in REQUIRED_CHECKS:
        for remove in (False, True):
            checks = dict(good)
            if remove:
                checks.pop(name)
            else:
                checks[name] = False
            result = verdict(checks)
            expected = "missing: "+name if remove else name
            if result["failed_checks"] != [expected] or any(result[k] != "UNRESOLVED" for k in VERDICT_FIELDS):
                return False
    return all(verdict({})[k] == "UNRESOLVED" for k in VERDICT_FIELDS)
