# Pre-registration: coupled degree-1 multiplet–metric tensor response

Baseline: `fe3d423d322fd1ea8760dc4aee09dff7fdf2acb9` (#294).
Branch: `codex/coupled-multiplet-response`. Seed: `2026091014`.
Publish this document before implementation, finite variations, integration
or Floquet measurements. The formulas below are analytic predictions worked
out while scoping the question, not discoveries from a numerical scan.

## 1. Question and the first closed sector

#294 gives exact ESU support with four independent real conformal fields.
Does their metric response preserve #292's bare homogeneous tensor equation,
or replace it? Derive that response from the specified field action and
verify closure of this linear sector in the full Einstein–matter equations.

The first calculation keeps all five homogeneous symmetric trace-free
metric components and tests the claim that **zero scalar perturbations**
are a dynamically consistent tensor response. This is a claim to prove,
not a rigid-fluid assumption. General scalar/vector perturbations,
inhomogeneous tensor towers, degree-3/5 supporting backgrounds and the
response of degree-5 equal-stress families are outside this first PR.
Failure of sector closure must be reported as such; no omitted equation may
be discarded to retain the desired oscillator.

Use exactly #294's four-component action with xi=1/6, kappa>0 and a>0.
No additional signal component, perfect fluid, damping, fifth-dimensional
junction, imposed rotor or constitutive law is introduced. On unit S3
embedding coordinates x_I, sum x_I^2=1, the background is

    phi_I(t,x)=P(t) x_I,
    P(t)=sqrt(3/(4 kappa)) cos(theta), theta=2t/a+delta,
    Lambda=3/(2a^2), rho=3/(2 kappa a^2).

This equals phi_I=s Y_I cos(theta), Y_I=2x_I, s^2=3/(16 kappa).
Use proper time on the background: lapse=1, shift=0. A pure linear tensor
perturbation has no lapse or shift perturbation. Do not transfer #290's
scalar-sector coordinate/proper-clock correction into this calculation.

Write the spatial metric in the inherited invariant coframe as

    g_ab=a^2 (exp(2 epsilon beta(t)))_ab,
    beta=sum_A b_A E_A, tr(beta)=0, tr(E_A E_B)=delta_AB.

All results about response refer to the derivative at epsilon=0. A finite
metric continuation is a differentiation control, not a nonlinear solution.

## 2. P1: prove or refute closure of the tensor sector

The degree-1 identity is Hess(x_I)=-x_I g/a^2. For a transverse trace-free
spatial metric perturbation, delta(Box) phi_I and delta R are predicted to
vanish. Thus delta phi_I=delta phi_dot_I=0 initially is preserved by each
linearized Klein–Gordon equation. This statement must be checked for every
component, with odd parity retained, and at both field and velocity zeros.

All linear Hamiltonian, momentum and spatial-trace Einstein residuals are
predicted to vanish in this sector. Check them independently using full
metric curvature and the full summed improved stress, not by simply setting
scalar/vector coefficients to zero in the candidate tensor equation.
Arbitrary unsourced scalar perturbations are not asserted to vanish; their
absence is the initial-data restriction defining this tensor experiment.

## 3. P2: the action and full stress must give the same response

Define Q=P^2, f=1-kappa Q/6 and C=Vol(S3)/kappa. The predicted quadratic
homogeneous tensor action is

    L_2=C/2 [ f tr(beta_dot^2) - g tr(beta^2)/a^2 ],
    g=8+(2/3) kappa Q=8+(1/2) cos^2(theta).

Derive this from the full action, including the curvature coupling to the
nonzero supporting fields. Terms of order beta^2 phi_background^2 cannot be
omitted: the background is not small in the signal expansion. The inherited
bare TT/scalar cubic truncation is therefore not an adequate starting action.

Useful independent action inputs to verify are

    R_3=6/a^2 - (8/a^2) tr(beta^2) + O(beta^3),
    R_4=R_3+tr(beta_dot^2)+O(beta^3) for unit lapse/fixed volume,
    sum_I |grad phi_I|^2=(Q/a^2) tr(exp(-2 beta)).

The expected full mixed spatial trace-free stress variation is

    delta T^a_b|TF = -(2Q/a^2) beta
                       + (Q_dot/6) beta_dot
                       + (Q/6) delta G^a_b|TF,
    delta G^a_b|TF = beta_ddot+(8/a^2) beta.

This yields the same five-component equation as action variation:

    f beta_ddot + f_dot beta_dot + (g/a^2) beta = 0.

Tensor components are mixed or measured in the perturbed physical frame;
background pressure times the metric perturbation must not be mistaken for
anisotropic stress. Off-shell finite variations must test the sign and
normalization of the beta, beta_dot and beta_ddot terms separately, including
noncommuting symmetric matrices and independent spatial points.

At Q=0 instantaneously, the bare stiffness is recovered, but a zero crossing
is not an all-time decoupling. An auxiliary multiplier eta on the response
coefficients (eta=0 bare, eta=1 physical) is permitted only as a mathematical
control; eta!=1 is not another self-consistent ESU background.

## 4. P3: bare-frequency transfer and periodic linear evolution

In tau=t/a the predicted equation is

    f b'' + f' b' + g b=0,
    f=1-cos^2(2tau+delta)/8,
    f'=sin(4tau+2delta)/4,
    g=8+cos^2(2tau+delta)/2.

Its coefficient period is pi/2. At tau=delta=0, b=1, b'=0,
b''=-8 (the bare oscillator), the predicted residual is 3/2. This is a
precomputed non-transfer control, not a numerical instability criterion.
A nonzero bare cosine must fail the all-time physical equation. Reversing
the f' term or dropping the metric dependence of the supporting stress are
required negative controls at phases where the relevant coefficients do
not vanish.

Integrate the fundamental response for each of the two independent initial
conditions and all five STF components. The canonical variables are
(b,p=f b') after dividing momentum by C/a. The fundamental matrix must
preserve their symplectic form, equivalently f times the velocity Wronskian.
Check independent normal-form evolution of y=sqrt(f) b, whose coefficient is

    g/f - f''/(2f) + (f')^2/(4 f^2).

Do not treat the f' term as dissipative friction. It comes from a periodic
kinetic coefficient and has a Hamiltonian formulation. Compare coordinate
and physical time at a=.7,1,2; kappa=.4,1 changes field amplitude, not the
normalized tensor response.

Numerically measure monodromy over tau in [0,pi/2], its determinant, trace,
eigenvalues and convergence. Report trace-based classification only as
NUMERICALLY_ELLIPTIC, NUMERICALLY_HYPERBOLIC or UNRESOLVED near a band edge.
No rigorous Floquet stability proof or full-system stability claim follows
from a numerical period map. If giving a quasifrequency, state its modulo-4
ambiguity in tau units; compare the directly invariant trace instead of
choosing an unregistered frequency branch. Also integrate 20 coefficient
periods and compare with repeated monodromy action. No numerical Floquet
outcome, quasifrequency or amplification is presumed here.

## 5. Frozen checks and failure behavior

Use symbolic rational identities for the degree-1 Hessian/gradient sums,
curvature/action expansion, full stress coefficients and action/Einstien
agreement. Every advertised symbolic residual must be exactly zero.
For independent full-curvature checks use coordinate metric jets and
Christoffel/Ricci evaluation, with symmetric differences at epsilon=.002,
.001,.0005. At deterministic generic points/phase data require relative
response errors below 2e-4 on the finest step and second-order halving
ratios between 3.5 and 4.5 where the leading error is above 1e-10.
Record absolute errors and normalization max(1, norm(expected)). Include
pure beta, pure beta_dot, pure beta_ddot and generic mixed controls.

For DOP853 use (rtol,atol)=(1e-10,1e-12) and (1e-12,1e-14), maximum step
pi/100 in tau. Require period-map refinement, determinant-one, independent
normal-form, and 20-period comparisons below 1e-8 in normalized norm.
Classify elliptic only for |trace|<2-1e-7 and hyperbolic only for
|trace|>2+1e-7, with both solver tolerances agreeing. Otherwise UNRESOLVED.
These are numerical gates, not interval enclosures of the exact trace.

REQUIRED_CHECKS, fixed here before implementation:

| Gate | Meaning |
|---|---|
| `background_and_parity` | inherited exact ESU, normalization, four odd components |
| `scalar_sector_closure` | all four unsourced linear KG equations |
| `quadratic_action` | full action includes nonzero-background curvature terms |
| `independent_curvature_stress` | full geometric/stress variations, mixed-frame convention |
| `constraint_completion` | Hamiltonian, momentum and trace residuals |
| `bare_frequency_control` | bare solution fails; response-sign controls discriminate |
| `canonical_propagation` | proper-time scaling, symplectic map, normal-form comparison |
| `period_map_convergence` | two tolerances, determinant, 20 periods, finite trace |
| `scope_and_order` | linear tensor restriction and unestablished broader claims |
| `failure_paths` | every missing/false gate and malformed result fails closed |

All ten gates are required for each affirmative physical verdict in this
first coupled-response experiment. Check exactly this key set. Missing,
false or unknown keys, nonfinite period-map data, or an unclassified band
edge make the affected result UNRESOLVED; no inferred no-go or stability
claim replaces failed numerical evidence. Exercise every named gate as
missing and false, checking CLI nonzero exit and overwriting stale JSON and
Markdown. Computation exceptions must also overwrite stale reports.

Keep separate verdict fields: `linear_tensor_sector`,
`coupled_tensor_equation`, `bare_frequency_transfer`, `period_map_type`.
The first three require the full gate set; period-map type also requires
valid trace/convergence evidence. Fixed scope fields remain
`general_scalar_vector_response: NOT_DERIVED`,
`full_dynamical_stability: NOT_ESTABLISHED`,
`nonlinear_persistence: NOT_ESTABLISHED`,
`preparation_selection: NOT_DERIVED`, `Phi_selection: NOT_DERIVED`,
`causality_gate: OPEN`. A failed numerical witness cannot establish
nonexistence of a more general coupled solution.

Deliver an isolated module, exact and independent geometric checks, probe,
archive, targeted tests and derivation. Preserve this freeze. Any corrected
prediction or later analytic extension must be identified explicitly.
