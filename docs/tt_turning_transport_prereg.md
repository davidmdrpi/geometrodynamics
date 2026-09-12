# Pre-registration: TT turning point, instantaneous actions and future data

Computational baseline: `6099a0baba1a14536ebbaad50323468c3281836e` (#297).
Main at `2d9e62d1b83ae710655933bb08f6396e902f9ecd` is integrated without
changing that baseline tree. Branch: `codex/tt-turning-transport-prereg`.
Seed: `2026091218`.

Publish this freeze before implementation, symbolic verification, phase
scans, root finding or propagation for this round. This document contains
analytic predictions worked out while scoping and predictions reported in
the [independent review of #297](https://github.com/davidmdrpi/geometrodynamics/pull/297#issuecomment-5648068582).
They are predictions to verify, not discoveries of an unrun experiment.
Keep the freeze intact if a prediction fails; document corrections separately.

## 1. The question and the meaning of support independence

On #297's fully supported expanding FRW family, does the canonical normal
equation turn near A=3a for every matter phase? What happens to two explicitly
defined instantaneous oscillator actions? Can initial classical tensor data
be transported to a complete, normalized pair of future asymptotic data?

The prior is a unique simple zero of the normal potential within one percent
of A=3a, failure of adiabaticity there, and a finite constant tensor amplitude
in the asymptotic future accompanied by a second, decaying mode. The action
definitions below are not expected to approach a finite nonzero conserved
quantity through that departure. This is not a conjecture that every possible
classical invariant is absent.

Use exactly the homogeneous five-component TT block, four real conformal
degree-1 fields and action established in #294–#297. No general tensor tower,
rotor drive, alternative matter sector or perfect-fluid response is added.
“Nearly support-independent” means weak dependence on the phase/departure
parameters within this specified family. #294's higher-degree or equal-stress
Gram families have not inherited this closed TT operator; this freeze makes
no universal statement about their turning points.

## 2. P1: exact phase reduction and normalization budget

The physical coefficients, normalized coefficients and momentum are distinct:

    M=A^2/kappa-q^2/6, K=8A^2/kappa+2q^2/3,
    m=kappa M/a^2, k=kappa K/a^2,
    p=m beta', p_physical=Vol(S3_unit) (a^2/kappa) p.

The inherited coefficients() API returns (m,m',m'',k), not (M,M',M'',K).
Throughout this round y=sqrt(m) beta. The physical action associated with
the normalized pair differs by the fixed factor Vol(S3_unit) a^2/kappa.
Record that factor rather than comparing unlike action normalizations.

Restrict to the expanding branch above the ESU: initial departure d>0,
sigma=+1, q=a sqrt(3/(4 kappa)) cos(2eta+delta). Its future conformal endpoint
is eta_star=log((2+d)/d)/sqrt(2). Define

    x=eta_star-eta > 0,
    R=A/a=coth(x/sqrt(2)),
    alpha=delta+2eta_star (mod pi), theta=alpha-2x.

Then the predicted exact reduction of the normalized tensor operator is

    m=R^2-cos^2(theta)/8,
    k=8R^2+cos^2(theta)/2,
    (m beta')'+k beta=0,
    W=k/m-m''/(2m)+m'^2/(4m^2), y''+W y=0,

where primes still mean d/deta=-d/dx. The bare control is W0=9-R^2.
The phase alpha is the only remaining operator parameter: a and kappa set
normalization/clock units, and d fixes the time origin. This must be checked
against #297's API at matched times and phases before using the reduction.
Do not lose a derivative sign when switching from eta to x.

The entire open expanding interval has m>=7/8. The limit x=0 is a singular
conformal endpoint at infinite proper time, not a finite-time breakdown of
the supporting background. No integration may step across x=0.

## 3. P2: a turning-point conjecture with a genuine all-phase gate

Freeze the conjecture:

    For every alpha in [0,pi), W has exactly one zero on 1<R<infinity.
    That zero is simple and satisfies 2.97 <= R_turn <= 3.03.

This is a one-percent bound about the exact support-free marker R_turn=3.
The review's sampled values near 3.001–3.018 motivate the conjecture but do
not prove it for other phases. Do not tighten the frozen bound around the
observed extrema after scanning, and do not turn grid agreement into a proof.

First derive W explicitly in R, sin(theta), cos(theta). Use the positive
denominator m to seek analytic bounds, polynomial inequalities or certified
interval bounds. The global claim needs positivity below the strip,
negativity above it, and uniqueness/simplicity along theta'=2 and
R'=(R^2-1)/sqrt(2) within the strip. An unvalidated numerical interval package
or a solver residual is not a global certificate. If such a certificate is
not obtained, report UNRESOLVED_FOR_ALL_PHASES and retain the sampled result
as a separately labeled diagnostic. A counterexample falsifies the conjecture.

Root scans use alpha=j*pi/64, j=0..63, plus 16 seeded phases, a bracketed
solver and absolute root tolerance 1e-11. Include the bare operator as an
exact R=3 control. Report R_turn, x_turn, W residual and dW/deta. A derivative
bounded away from zero, rather than a sign-changing plot alone, supports
the simple-root interpretation. Check positive/negative regions independently.

The e-fold marker ln(R_turn) is measured from the asymptotic ESU radius a.
From an initial radius a(1+d), the corresponding expansion is
ln(R_turn/(1+d)). Neither convention specifies a universal proper departure
time without initial data.

## 4. P3: define the action before testing its survival

Two different instantaneous quantities are to be evaluated, in fixed
normalizations. With Omega=sqrt(k/m)>0 throughout the open interval, define

    H_beta=(p^2/m+k beta^2)/2,
    J_beta=H_beta/Omega
          =[p^2/(m Omega)+m Omega beta^2]/2.

For y=sqrt(m) beta, define only in the region W>0

    omega=sqrt(W), J_y=(y'^2+W y^2)/(2 sqrt(W)),
    epsilon_WKB=|omega'|/omega^2=|W'|/(2 W^(3/2)).

Do not replace W by |W| and call the continued expression an oscillatory
action. Do not call J_y undefined at W<=0 evidence that the original beta
equation is singular: m and k remain positive at the turning point.

Precomputed identity to verify:

    J_beta' = (m Omega)' [beta^2-p^2/(m Omega)^2]/2.

Thus the specific instantaneous action need not be conserved by a
time-dependent Hamiltonian even though the flow is symplectic. An action
derivative vanishing at an isolated oscillator phase is not conservation.
Test the coefficient identity for arbitrary canonical data and numerically
for 16 initial oscillator phases at R=1.05, normalized to J_beta=1:

    beta=sqrt(2/(m Omega)) cos(chi),
    p=sqrt(2m Omega) sin(chi), chi=2pi j/16.

At a simple zero of W, epsilon_WKB is predicted to diverge as distance to
the zero to the power -3/2. The review's description of small adiabaticity
“right up to” the marker does not apply to this standard diagnostic. Check
approach from W>0 at distances 1e-2,1e-3,1e-4,1e-5 in conformal time.

J_y generically diverges when y' at the zero is nonzero. It is not predicted
to diverge for every initial condition: include a solution initialized with
(y,y')=(1,0) at the zero and evolved backwards, for which J_y approaches zero.
The definition fails at the zero in either case. Do not infer the behavior
of the whole solution space from a generic example.

An exact quadratic invariant is always constructible from a chosen invertible
fundamental matrix S: for z=(beta,p), I=(S^-1 z)^T G(S^-1 z)/2 with fixed
positive G is constant. Verify this as a control against an overbroad
“no classical invariant exists” verdict. Its chosen S and G are extra basis
and normalization inputs; it is not automatically a preferred action or a
quantized count. No selection of Phi follows from either invariant test.

## 5. P4: the future keeps two coefficients, not only a frozen amplitude

The preliminary endpoint expansion predicts

    R=sqrt(2)/x+O(x), m=2/x^2+O(1), k/m=8+O(x^2),
    beta_xx-(2/x) beta_x+8 beta+...=0.

The Frobenius indices are 0 and 3. A real future basis is predicted to exist
with the following normalization in the inherited (beta,p=m beta') pair:

    b_c(x)=1+4x^2+O(x^4), coefficient of x^3 fixed to zero,
    b_d(x)=-x^3/6+O(x^5),
    p_c=-16/x+O(x), p_d=1+O(x^2).

The difference of indices is an integer. Explicitly check the resonant
recurrence and absence of a forced logarithm; do not assume these series are
valid because leading powers fit. If a logarithm is forced, preserve the
freeze and correct the basis construction and associated verdict.

The chosen sign and factor -1/6 predict the exact Wronskian

    m(b_c b_d'-b_c' b_d)=1.

Derive it and carry both columns of the future fundamental matrix. The
coefficient data are (C,D), beta=C b_c+D b_d. Then beta tends to C, but D
remains the independent decaying coefficient. Projection onto C alone is
rank one and cannot preserve a two-dimensional symplectic form.

For these particular instantaneous actions, further endpoint predictions are

    J_beta ~ sqrt(8) C^2/x^2       for C != 0,
    J_beta ~ D^2 x^2/(4 sqrt(8))  for C = 0, D != 0.

Verify coefficients, powers and exceptional cases. The growing normal
coordinate y~sqrt(2) C/x must not be mistaken for a growing physical tensor
amplitude beta. The turning point near 3a and the limit beta -> C occur at
different locations; this freeze does not identify them.

## 6. P5: input basis and the complete transport map

No preferred positive-frequency vacuum is inherited. Use a real canonical
input convention so that a future map is reproducible without asserting one.
Set eta_star=0, so eta=-x and the past reference ESU has phase alpha. Let
S_ESU(eta,0) be the #295 canonical fundamental matrix with S_ESU(0,0)=I.
Its elliptic Floquet block is a reference for past propagation, not a
quantization prescription. Check its monodromy, determinant and bounded
periodic/Floquet representation independently of the evolving background.

At the fixed matching point x_match=1 define the finite-cutoff input matrix

    B_in(L)=U_FRW(-1,-L) S_ESU(-L,0),

where U_FRW maps the normalized pair (beta,p). The prediction is that B_in(L)
converges as L->infinity, since the operator approaches the supported ESU
exponentially, while the reference flow remains bounded. The reference time
0 is a basis convention, not a claim that the physical endpoint is an ESU.

Construct B_out(-1) by the normalized Frobenius basis, propagated from small
positive x to the matching point. The proposed complete map is

    (C,D)^T = T z_in, T=B_out(-1)^-1 lim_L B_in(L).

Both input and output bases have unit canonical Wronskian, so det T=1 and
T^T J T=J, J=[[0,1],[-1,0]], are necessary checks. They are not a substitute
for cutoff, basis and full two-column convergence. Report all four matrix
entries and their errors, not only the first row or an amplitude histogram.

Use L=8,12,16,20 for input convergence; Frobenius orders 8,10,12 and starting
x=0.04,0.02,0.01 for output convergence. Match away from x=0 to avoid
subtracting divergent momenta or nearly parallel late-time solutions.
Check x_match=0.8,1.0,1.2 independently. Use alpha=j*pi/16, j=0..15, for the
full transport experiment, and the stated 64+16 phases for the root test.
These grids are fixed before results and do not constitute all-phase proofs.

Integrate with DOP853 at (rtol,atol)=(1e-10,1e-12),(1e-12,1e-14), with
max_step<=pi/100 in eta; reduce the step relative to distance from x=0 when
needed. Require normalized equation/series residuals <1e-8 on their stated
overlap, relative matrix disagreement <1e-8 for the final cutoff/order/tolerance
comparisons, and determinant/symplectic residuals <1e-8. Report raw absolute
errors too. If double precision cannot resolve a coefficient, use higher
precision or report it unresolved; never discard the decaying column.

Verify endpoint/action predictions both analytically and with independent
integration. Proper time is dt=A deta; repeat representative clock checks
at a=0.7,1,2 and kappa=0.4,1. Match d=0.05,0.15,0.30 with delta chosen to
hold alpha fixed, testing the reduction rather than treating time translates
as independent physics. The primary predictions are classical and in
conformal variables; no cosmic-time frequency convention is silently added.

## 7. Failure controls, verdicts and stopping rule

Required gate names and target dependencies:

| Gate | Targets |
|---|---|
| normalization | B, A, F, T |
| phase_reduction | B, A, F, T |
| turning_location | B |
| adiabatic_diagnostic | A |
| instantaneous_actions | A |
| frobenius_basis | F, T |
| past_basis | T |
| transport_convergence | T |
| symplectic_completion | T |
| negative_controls | B, A, F, T |
| scope | B, A, F, T |
| failure_paths | B, A, F, T |

B is the all-phase turning-bound verdict; A is the verdict for the two named
instantaneous actions and the WKB diagnostic; F is the asymptotic basis;
T is the complete real input/output map. Freeze and test the dependency map
itself. Every missing/false dependent gate or malformed/nonfinite evidence
must make its target UNRESOLVED. Evidence for unrelated targets remains
visible. A numerical scan without a global certificate cannot make B pass.

Required controls include the bare R=3 root; a constant oscillator with
conserved instantaneous action; the tuned y'=0 turning-point solution;
the pure decaying future solution; an exact pulled-back quadratic invariant;
an incorrect eta/x derivative sign; omitted input or output normalization;
and a rank-one frozen-amplitude projection falsely offered as a full map.
Exercise failure paths from affirmative, independently reproduced evidence,
including actual CLI failures, exceptions, missing certificates and stale
success artifacts. A string-ban test is not a substitute for the controls.

Use --output-dir, retaining --output as an optional alias. Archive formulas,
certificate status, numerical selections, raw maps and convergence evidence.
The next result may be a certified bound, a counterexample or an honestly
unresolved certificate; do not optimize its wording after seeing a residual.

Stop after these classical questions. No particle-production number,
preferred complex structure, quantization, count of histories, Phi selection,
rotor construction, general matter-sector result or operational causality
claim is licensed. Spatial harmonic labels, instantaneous actions, exact
basis-dependent invariants and asymptotic amplitude data remain distinct.
