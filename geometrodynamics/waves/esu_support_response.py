"""ESU perfect-fluid response and a preparation-scoped cubic scalar force.

Public freeze c07891d. The support's adiabatic derivative is a parameter,
not a derived BAM constitutive law. Scalar constraints propagate at order
s^2; the projected s^3 force is not a complete corrected scalar evolution.
"""

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.integrate import solve_ivp

from . import reciprocal_scalar_tt as rt
from . import scalar_tt_constraints as sc


PUBLIC_PREREG = "c07891d573ab719ba90e73f8f3309c8aad674350"
BASELINE = "080c1cc58cd7a01fc619ff5d933c5fdfb5b8cf33"
SEED = 2026090713
DEGREES = np.array([0, 2, 4, 6])
POWERS = np.array([1., 27/25, 1/5, 201/175])


@dataclass(frozen=True)
class SupportResponse:
    """Four harmonic-group amplitudes multiplying H_l=(Y^2)_l.

    Reduced state: (psi_l, psi_dot_l). Independent fluid state appends
    (rho_f,l, J_f,l). H_l has physical squared norm POWERS[l]/volume.
    """

    sound_speed_squared: float = 1/3
    amplitude: float = .02
    radius: float = 1.
    kappa: float = 1.

    def __post_init__(self):
        rt.ReciprocalModel(radius=self.radius, kappa=self.kappa)
        if not math.isfinite(self.sound_speed_squared) or not 0 <= self.sound_speed_squared <= 1:
            raise ValueError("this response family uses 0 <= c_s^2 <= 1")
        if not math.isfinite(self.amplitude) or self.amplitude <= 0:
            raise ValueError("positive finite amplitude required")

    @property
    def volume(self):
        return 2*math.pi**2*self.radius**3

    @property
    def enthalpy(self):
        return 2/(self.kappa*self.radius**2)

    @property
    def lambdas(self):
        return DEGREES*(DEGREES+2)/self.radius**2

    @property
    def L(self):
        return 3/self.radius**2-self.lambdas

    @property
    def omega2(self):
        return 16/self.radius**2

    @property
    def cubic_scale(self):
        return self.kappa*self.amplitude**3/(self.volume*self.radius**2)

    def sources(self, time):
        omega = math.sqrt(self.omega2)
        A = self.amplitude*math.cos(omega*time)
        dA = -self.amplitude*omega*math.sin(omega*time)
        rho = self.amplitude**2*self.omega2/2-self.lambdas*A*A/12
        J = np.full(4, -A*dA/6)
        dJ = np.full(4, -(dA*dA-self.omega2*A*A)/6)
        Pi = -A*A*(self.omega2-self.lambdas/12)/(2*self.L)
        dPi = -A*dA*(self.omega2-self.lambdas/12)/self.L
        # Constant potentials have no gradient/Hessian meaning.
        J[0] = dJ[0] = Pi[0] = dPi[0] = 0.
        return {"A": A, "dA": dA, "rho": rho, "J": J, "dJ": dJ, "Pi": Pi, "dPi": dPi}

    def initial(self, independent_fluid=False):
        psi = self.kappa*self.sources(0.)["rho"]/(2*self.L)
        state = np.r_[psi, np.zeros(4)]
        return np.r_[state, np.zeros(8)] if independent_fluid else state

    def fields(self, time, state):
        psi, dpsi = np.asarray(state[:4]), np.asarray(state[4:8])
        src = self.sources(time)
        if len(state) == 16:
            rho_f, J_f = np.asarray(state[8:12]), np.asarray(state[12:16])
        elif len(state) == 8:
            rho_f = 2*self.L*psi/self.kappa-src["rho"]
            J_f = -2*dpsi/self.kappa-src["J"]
            J_f[0] = 0.
        else:
            raise ValueError("expected reduced (8) or independent-fluid (16) state")
        alpha = psi-self.kappa*src["Pi"]
        dalpha = dpsi-self.kappa*src["dPi"]
        alpha[0] = dalpha[0] = 0.
        ddpsi = (psi/self.radius**2+self.kappa*(self.sound_speed_squared*rho_f+src["rho"]/3)/2
                 - self.kappa*self.lambdas*src["Pi"]/3)
        curvature = 4*self.L*psi-6*ddpsi+2*self.lambdas*alpha
        trace_curvature = self.kappa*(1-3*self.sound_speed_squared)*rho_f
        return {**src, "psi": psi, "dpsi": dpsi, "ddpsi": ddpsi,
                "alpha": alpha, "dalpha": dalpha, "rho_f": rho_f, "J_f": J_f,
                "curvature": curvature, "trace_curvature": trace_curvature}

    def reduced_rhs(self, time, state):
        src = self.sources(time)
        psi, dpsi = state[:4], state[4:8]
        ddpsi = ((1/self.radius**2+self.sound_speed_squared*self.L)*psi
                 + self.kappa*(1/3-self.sound_speed_squared)*src["rho"]/2
                 - self.kappa*self.lambdas*src["Pi"]/3)
        return np.r_[dpsi, ddpsi]

    def fluid_rhs(self, time, state):
        """Evolve fluid conservation and spatial Einstein; no constraint solve."""
        f = self.fields(time, state)
        drho = self.lambdas*f["J_f"]+3*self.enthalpy*f["dpsi"]
        dJ = -self.sound_speed_squared*f["rho_f"]-self.enthalpy*f["alpha"]
        dJ[0] = 0.
        return np.r_[f["dpsi"], f["ddpsi"], drho, dJ]

    def integrate(self, times, independent_fluid=False, rtol=1e-12, atol=1e-14):
        times = np.asarray(times)
        if times.ndim != 1 or len(times) < 2 or times[0] != 0 or np.any(np.diff(times) <= 0):
            raise ValueError("strictly increasing times starting at zero required")
        sol = solve_ivp(self.fluid_rhs if independent_fluid else self.reduced_rhs,
                        (0., float(times[-1])), self.initial(independent_fluid), t_eval=times,
                        method="DOP853", rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.y.T

    def constraints(self, time, state):
        f = self.fields(time, state)
        ham = 2*self.L*f["psi"]-self.kappa*(f["rho_f"]+f["rho"])
        mom = -2*f["dpsi"]-self.kappa*(f["J_f"]+f["J"])
        return ham, mom[1:]

    def coordinate_force(self, time, state):
        """Projection of the s^3 Newtonian-coordinate wave force onto Y."""
        f = self.fields(time, state)
        integrand = (-2*self.omega2*f["A"]*f["alpha"]
                     + f["dA"]*(f["dalpha"]+3*f["dpsi"])
                     - 30*f["A"]*f["psi"]/self.radius**2
                     + f["A"]*self.lambdas*(f["alpha"]-f["psi"])/2
                     - f["A"]*f["curvature"]/6)
        return float(POWERS @ integrand/self.volume)

    def initial_proper_force(self):
        """Initial free-fall fluid-clock force, not coordinate acceleration."""
        psi = self.initial()[:4]
        return float(-self.amplitude*(POWERS @ ((30/self.radius**2+self.lambdas/2)*psi))/self.volume)

    def induced_tt_force(self, time):
        T = 8/self.radius**2
        B = ((1-math.cos(math.sqrt(T)*time))/(2*T)
             +(math.cos(2*math.sqrt(self.omega2)*time)-math.cos(math.sqrt(T)*time))/(2*(T-4*self.omega2)))
        return 48*self.kappa*self.amplitude**3*B*math.cos(math.sqrt(self.omega2)*time)/(self.volume*self.radius**4)


    def continuous_potential_bounds(self, duration):
        """All-space, all-time upper bounds on [0,duration], including l=0.

        Bound the constant/cos(2 omega t) forcing convolution, then use the
        S3 addition theorem ||H_l||_infinity <= (l+1)sqrt(POWERS_l)/V.
        """
        if not math.isfinite(duration) or duration < 0:
            raise ValueError("finite nonnegative duration required")
        q = (1+self.sound_speed_squared*(3-DEGREES*(DEGREES+2)))/self.radius**2
        Cmax = np.cosh(np.sqrt(np.maximum(0., q))*duration)
        r_const = self.amplitude**2*(self.omega2/2-self.lambdas/24)
        r_cos = -self.amplitude**2*self.lambdas/24
        pi_const = -self.amplitude**2*(self.omega2-self.lambdas/12)/(4*self.L)
        pi_const[0] = 0.
        g_const = (self.kappa*(1/3-self.sound_speed_squared)*r_const/2
                   - self.kappa*self.lambdas*pi_const/3)
        g_cos = (self.kappa*(1/3-self.sound_speed_squared)*r_cos/2
                 - self.kappa*self.lambdas*pi_const/3)
        psi = (np.abs(self.initial()[:4])*Cmax
               + np.abs(g_const)*duration**2*Cmax/2
               + np.abs(g_cos)*(Cmax+1)/(q+4*self.omega2))
        alpha = psi+2*self.kappa*np.abs(pi_const)
        alpha[0] = 0.
        spatial_factor = (DEGREES+1)*np.sqrt(POWERS)/self.volume
        return {"psi": float(spatial_factor @ psi), "alpha": float(spatial_factor @ alpha)}


class SpatialFields:
    """Independent pointwise geometry and stress projections on S3."""

    def __init__(self, model, radial_order=8, angular_order=16, axis=None):
        self.model = model
        self.grid = sc.EvenHarmonicGrid(radial_order, angular_order, model.radius)
        self.wave = sc.StandingWaveConstraints(self.grid, amplitude=model.amplitude,
                                              kappa=model.kappa, axis=axis)
        self.H_coeff = np.array([np.where(self.grid.degrees == l, self.wave.h, 0.) for l in DEGREES])
        self.H = self.grid.Y @ self.H_coeff.T
        self.grad_H = np.stack([self.grid.gradient(c) for c in self.H_coeff], axis=1)
        self.lap_H = np.stack([self.grid.evaluate(self.grid.laplacian_coefficients(c)) for c in self.H_coeff], axis=1)

    def jets(self, time):
        src = self.model.sources(time)
        direction = self.wave.direction
        return rt.scalar_jets(self.wave.model, src["A"]*direction, src["dA"]*direction,
                              -self.model.omega2*src["A"]*direction, self.grid.points)[0]

    def stress_check(self, time):
        if self.model.radius != 1.:
            raise ValueError("inherited improved stress requires unit radius")
        from .backreaction import stress_series
        stress = stress_series(self.jets(time))[:, 0]
        rho, j = stress[:, 0, 0], -stress[:, 0, 1:]
        pressure = np.trace(stress[:, 1:, 1:], axis1=1, axis2=2)/3
        anis = stress[:, 1:, 1:]-pressure[:, None, None]*np.eye(3)
        coeff = np.einsum("pn,p,pij->ijn", self.grid.Y, self.grid.weights, anis)
        scalar_projection = np.zeros(84)
        for i in range(3):
            for jdx in range(3):
                operator = (self.grid.D[i] @ self.grid.D[jdx]+self.grid.D[jdx] @ self.grid.D[i])/2
                if i == jdx:
                    operator -= sum(d @ d for d in self.grid.D)/3
                scalar_projection += operator.T @ coeff[i, jdx]
        norm = (2/3)*self.grid.lambdas*(self.grid.lambdas-3)
        scalar_projection[1:] /= norm[1:]
        scalar_projection[0] = 0.
        src = self.model.sources(time)
        expected_pi = src["Pi"] @ self.H_coeff
        return {"rho": scaled_error(self.H @ src["rho"], rho),
                "pressure": scaled_error(self.H @ (src["rho"]/3), pressure),
                "momentum": scaled_error(np.einsum("pli,l->pi", self.grad_H, src["J"]), j),
                "anisotropic_scalar_projection": scaled_error(scalar_projection, expected_pi)}

    def force_fields(self, time, state):
        f, jets = self.model.fields(time, state), self.jets(time)
        phi, dt, dtt = (jets[k][:, 0] for k in ("phi", "dt", "dtt"))
        grad = jets["grad"][:, 0]
        lap = jets["laplacian"][:, 0]
        alpha, psi, dalpha, dpsi, curvature = (self.H @ f[k] for k in ("alpha", "psi", "dalpha", "dpsi", "curvature"))
        grad_alpha = np.einsum("pli,l->pi", self.grad_H, f["alpha"])
        grad_psi = np.einsum("pli,l->pi", self.grad_H, f["psi"])
        Ft = (2*alpha*dtt+(dalpha+3*dpsi)*dt+2*psi*lap
              + np.sum((grad_alpha-grad_psi)*grad, axis=1)-curvature*phi/6)
        # At t=0 the support is initially at rest with no pressure gradient.
        # This conversion includes the coordinate acceleration of the fluid.
        initial_proper = Ft-2*alpha*dtt-dalpha*dt-np.sum(grad_alpha*grad, axis=1)
        geometric_initial = 2*psi*lap-np.sum(grad_psi*grad, axis=1)-curvature*phi/6
        return {"coordinate": Ft, "initial_proper_conversion": initial_proper,
                "initial_proper_geometry": geometric_initial}

    def project_force(self, force):
        return float(self.grid.weights @ (self.wave.mode*force))

    def exact_coordinate_acceleration(self, time, state, epsilon):
        """Nonlinear exponential test metric with full ADM scalar curvature.

        This is a metric-variation control, not a nonlinear Einstein solution.
        """
        f, jets = self.model.fields(time, state), self.jets(time)
        alpha, psi, dalpha, dpsi, ddpsi = (self.H @ f[k] for k in ("alpha", "psi", "dalpha", "dpsi", "ddpsi"))
        ga, gp = (np.einsum("pli,l->pi", self.grad_H, f[k]) for k in ("alpha", "psi"))
        la, lp = (self.lap_H @ f[k] for k in ("alpha", "psi"))
        e = epsilon
        R = (np.exp(2*e*psi)*(6/self.model.radius**2+4*e*lp-2*e*la
                              -2*e*e*np.sum(gp*gp+ga*ga-gp*ga, axis=1))
             + np.exp(-2*e*alpha)*(12*e*e*dpsi*dpsi-6*e*ddpsi+6*e*e*dalpha*dpsi))
        return (e*(dalpha+3*dpsi)*jets["dt"][:, 0]
                + np.exp(2*e*(alpha+psi))*(jets["laplacian"][:, 0]
                     + e*np.sum((ga-gp)*jets["grad"][:, 0], axis=1))
                - np.exp(2*e*alpha)*R*jets["phi"][:, 0]/6)

    def exact_initial_clock_accelerations(self, epsilon):
        """Two derivatives of the same scalar in the exponential test metric.

        Normalize an initially resting geodesic, construct Gamma^mu_00 from
        metric derivatives, and contract its covariant scalar Hessian. This
        does not call the first-order force or clock-conversion formulas.
        Only the initial fluid clock is geodesic in the frozen preparation;
        this is not a later fluid history or a nonlinear Einstein solution.
        """
        state = self.model.initial()
        f, jets = self.model.fields(0., state), self.jets(0.)
        alpha, psi, dalpha = (self.H @ f[k] for k in ("alpha", "psi", "dalpha"))
        grad_alpha = np.einsum("pli,l->pi", self.grad_H, f["alpha"])
        # Spatial components use the background orthonormal frame. For
        # Gamma^mu_00, the frame's spatial commutators do not contribute.
        diagonal = np.column_stack((-np.exp(2*epsilon*alpha),
                                    np.repeat(np.exp(-2*epsilon*psi)[:, None], 3, axis=1)))
        d_g00 = 2*epsilon*diagonal[:, :1]*np.column_stack((dalpha, grad_alpha))
        connection_numerator = -d_g00
        connection_numerator[:, 0] += 2*d_g00[:, 0]
        gamma_00 = connection_numerator/(2*diagonal)
        U0_squared = -1/diagonal[:, 0]
        coordinate = self.exact_coordinate_acceleration(0., state, epsilon)
        d_phi = np.column_stack((jets["dt"][:, 0], jets["grad"][:, 0]))
        proper = U0_squared*(coordinate-np.sum(gamma_00*d_phi, axis=1))
        return {"coordinate": coordinate, "proper": proper}


def scaled_error(actual, expected):
    actual, expected = np.asarray(actual), np.asarray(expected)
    return float(np.linalg.norm(actual-expected)/max(1., np.linalg.norm(expected)))


@lru_cache(None)
def initial_coefficient_certificate():
    import sympy as sp
    powers = [sp.Rational(1), sp.Rational(27, 25), sp.Rational(1, 5), sp.Rational(201, 175)]
    proper, coordinate = [], []
    for l, power in zip((0, 2, 4, 6), powers):
        lam = sp.Integer(l*(l+2))
        psi = (8-lam/12)/(2*(3-lam))
        alpha = sp.Integer(0) if l == 0 else (24-lam/6)/(2*(3-lam))
        proper.append(sp.factor(-(30+lam/2)*psi*power))
        coordinate.append(sp.factor(((lam/2-32)*alpha-(30+lam/2)*psi)*power))
    return {"proper_by_degree": [str(x) for x in proper],
            "coordinate_by_degree": [str(x) for x in coordinate],
            "proper_total": str(sum(proper)), "coordinate_total": str(sum(coordinate)),
            "proper_expected": bool(sum(proper) == -sp.Rational(7976, 875)),
            "coordinate_expected": bool(sum(coordinate) == sp.Rational(55096, 875)),
            "scope": "initial zero-support-perturbation preparation; proper fluid clock and separately fixed coordinate clock"}


@lru_cache(None)
def metric_einstein_derivation():
    """First metric/connection variation in a round S3 coordinate chart.

    No support equation or constraint inverse is used to construct Ricci.
    Generic alpha(t,x) and psi(t,x) retain all time/spatial derivatives.
    """
    import sympy as sp
    t, chi, theta, az = sp.symbols("t chi theta azimuth", real=True)
    a = sp.Symbol("a", positive=True)
    coordinates = (t, chi, theta, az)
    alpha = sp.Function("alpha")(*coordinates)
    psi = sp.Function("psi")(*coordinates)
    phi = sp.Function("phi")(*coordinates)
    g = [-sp.Integer(1), a*a, a*a*sp.sin(chi)**2,
         a*a*sp.sin(chi)**2*sp.sin(theta)**2]
    inv = [1/x for x in g]
    h = [-2*alpha]+[-2*psi*x for x in g[1:]]
    def component(diagonal, i, j):
        return diagonal[i] if i == j else sp.Integer(0)
    def metric_derivative_sum(diagonal, d, b, c):
        return (sp.diff(component(diagonal, d, c), coordinates[b])
                + sp.diff(component(diagonal, d, b), coordinates[c])
                - sp.diff(component(diagonal, b, c), coordinates[d]))
    connection = sp.MutableDenseNDimArray.zeros(4, 4, 4)
    variation = sp.MutableDenseNDimArray.zeros(4, 4, 4)
    for d in range(4):
        for b in range(4):
            for c in range(4):
                base = metric_derivative_sum(g, d, b, c)
                connection[d, b, c] = sp.simplify(inv[d]*base/2)
                variation[d, b, c] = sp.simplify((inv[d]*metric_derivative_sum(h, d, b, c)
                                                  - inv[d]**2*h[d]*base)/2)
    ricci = sp.MutableDenseMatrix.zeros(4, 4)
    for i in range(4):
        for j in range(i, 4):
            value = sum(sp.diff(variation[c, i, j], coordinates[c])
                        - sp.diff(variation[c, i, c], coordinates[j])
                        for c in range(4))
            value += sum(variation[c, c, d]*connection[d, i, j]
                         + connection[c, c, d]*variation[d, i, j]
                         - variation[c, j, d]*connection[d, i, c]
                         - connection[c, j, d]*variation[d, i, c]
                         for c in range(4) for d in range(4))
            ricci[i, j] = ricci[j, i] = sp.simplify(value)
    ricci_bar = [sp.Integer(0)]+[2*x/a**2 for x in g[1:]]
    scalar = sp.simplify(sum(inv[i]*ricci[i, i]-inv[i]**2*h[i]*ricci_bar[i] for i in range(4)))
    def hessian(function, i, j):
        return (sp.diff(function, coordinates[i], coordinates[j])
                - sum(connection[k, i, j]*sp.diff(function, coordinates[k]) for k in range(4)))
    def laplacian(function):
        return sum(inv[i]*hessian(function, i, i) for i in range(1, 4))
    Q = psi-alpha
    geometric_scalar = 4*laplacian(psi)+12*psi/a**2-6*sp.diff(psi, t, 2)-2*laplacian(alpha)
    # Normal-frame energy variation includes the normalization of n^0.
    energy = ricci[0, 0]-h[0]*3/a**2+scalar/2-2*alpha*3/a**2
    energy_residual = sp.simplify(energy-2*(laplacian(psi)+3*psi/a**2))
    momentum = [sp.simplify(ricci[0, i]-2*sp.diff(psi, t, coordinates[i])) for i in range(1, 4)]
    spatial = []
    for i in range(1, 4):
        for j in range(1, 4):
            actual = inv[i]*ricci[i, j]
            if i == j:
                actual += -inv[i]**2*h[i]*ricci_bar[i]-scalar/2
            expected = inv[i]*hessian(Q, i, j)
            if i == j:
                expected += -laplacian(Q)+2*sp.diff(psi, t, 2)-2*psi/a**2
            spatial.append(sp.simplify(actual-expected))
    # Independent first variation of (sqrt(-g))^-1 d_mu(sqrt(-g)g^mu nu d_nu phi).
    density = a**3*sp.sin(chi)**2*sp.sin(theta)
    dvolume = alpha-3*psi
    dbox = sum(sp.diff(density*(-inv[i]**2*h[i]+dvolume*inv[i])*sp.diff(phi, coordinates[i]),
                       coordinates[i])/density for i in range(4))
    dbox -= dvolume*(-sp.diff(phi, t, 2)+laplacian(phi))
    expected_box = (2*alpha*sp.diff(phi, t, 2)+(sp.diff(alpha, t)+3*sp.diff(psi, t))*sp.diff(phi, t)
                    + 2*psi*laplacian(phi)
                    + sum(inv[i]*sp.diff(alpha-psi, coordinates[i])*sp.diff(phi, coordinates[i]) for i in range(1, 4)))
    curvature_residual = sp.simplify(scalar-geometric_scalar)
    wave_residual = sp.simplify(dbox-expected_box)
    cs = sp.Symbol("cs2", real=True)
    freq = cs*(sp.Symbol("ell", integer=True)*(sp.Symbol("ell", integer=True)+2)-3)-1
    residuals = [energy_residual, *momentum, *spatial, curvature_residual, wave_residual]
    return {"Hamiltonian_residual": str(energy_residual),
            "momentum_residuals": [str(x) for x in momentum],
            "spatial_Einstein_residuals": [str(x) for x in spatial],
            "curvature_residual": str(curvature_residual), "wave_operator_residual": str(wave_residual),
            "all_exact_zero": bool(all(x == 0 for x in residuals)),
            "fluid_only_frequency_times_a2": str(freq),
            "ell2_frequency_times_a2": str(5*cs-1),
            "ell0_frequency_times_a2": str(-3*cs-1),
            "ell2_threshold": str(sp.solve(5*cs-1, cs)[0]),
            "scope": "generic scalar metric variation; projected spatial scalar Einstein equations; no matter closure used in Ricci"}


REQUIRED_CHECKS = (
    "metric_Einstein_derivation", "improved_stress_and_scalar_anisotropy",
    "two_spatial_rules", "independent_fluid_evolution", "propagated_constraints",
    "ODE_refinement", "trace_curvature", "fluid_only_frequency",
    "nonlinear_metric_variation", "initial_proper_and_coordinate_coefficients",
    "clock_conversion", "amplitude_and_radius_scaling", "induced_TT_force",
    "small_metric_regime", "exact_preparation_obstruction")


def failed_checks(checks):
    return (["missing: "+name for name in REQUIRED_CHECKS if name not in checks]
            + [name for name, ok in checks.items() if not ok])


VERDICT_FIELDS = ("support_class", "scalar_constraints", "initial_proper_coefficient",
                  "cancellation", "complete_evolution", "BAM_support_selection", "readout")


def verdict(checks):
    failed = failed_checks(checks)
    if failed:
        return {**dict.fromkeys(VERDICT_FIELDS, "UNRESOLVED"), "failed_checks": failed}
    return dict(zip(VERDICT_FIELDS, (
        "PERFECT_FLUID_ADIABATIC_RESPONSE_ASSUMED_WITH_PARAMETER_CS2",
        "PROPAGATED_AT_ORDER_S2_IN_SCALAR_SECTOR",
        "MINUS_7976_OVER_875_TIMES_KAPPA_S3_OVER_V_A2",
        "IDENTICAL_CANCELLATION_EXCLUDED_FOR_FROZEN_PREPARATION",
        "NOT_ESTABLISHED", "NOT_DERIVED", "NOT_DERIVED")))
