# Round 12 freeze: odd multiplet support and preparation sensitivity

Baseline: `70fadb250273796bbbed2a674cf366b82a7d29ba` (#293), containing
main `68a9e86323d84c43427fd720979ff63497548ec7` and #292.
Branch: `codex/odd-multiplet-support`. Seed: `2026091012`.
This is a preregistration-only commit, to be published before this round's
implementation, amplitude scans or Gram-kernel calculations.

## 1. Prior state and question

#293 excludes exact round-ESU support by one smooth antipodally odd real
conformal scalar. Its post-freeze uniqueness argument classifies the
unrestricted single-field support as the homogeneous even family. It
explicitly leaves multiple fields outside the theorem's scope.

The user has already exhibited and checked an odd full-multiplet
construction, obtaining pointwise isotropic stress for k=1,3,5, constant
density through space and time, p=rho/3, and minimum kappa F values
0.875, 0.96875 and 0.986111.... These are **known candidate results supplied
before this freeze**, not discoveries to be credited to a blind search.
P1 and P2 below freeze them as explicit predictions for independent
re-derivation from the action. P3 has not been measured in this round;
its kernel, ranks and sensitivity outcomes are not presumed.

Question: does this expanded field content supply exact pointwise ESU
support, and which preparation perturbations preserve or break it?
Distinguish existence, constitutive response, preparation selection and
stability. A positive background is not a derived history measure.

## 2. Field content, normalization and admissibility

Use m independent real components of the same conformal action,

    I = integral sqrt(-g) [(R-2 Lambda)/(2 kappa)
              - (1/2) sum_I ((grad phi_I)^2 + R phi_I^2/6)],
    ds^2 = -dt^2 + a^2 dOmega_3^2,   a>0, kappa>0.

The sum is over **independent field components**, not over modes of a
single real field. The total stress is sum_I T[phi_I]. A single field
phi=sum_I phi_I instead has T[sum_I phi_I], including cross terms. These
two constructions must have distinct code paths and a discriminating
negative control. No ensemble or time average may be substituted for
the pointwise Einstein equations.

The extra component count is an assumption extending #293's action,
even though each term copies the existing scalar coupling. Every component
obeys the imposed odd antipodal condition; that condition is not derived.
The primary degree is k=1, with k=3,5 as higher-degree controls. A common
phase and equal coefficients are specified initial preparation, not an
equilibration mechanism. No self-interaction or phase-locking term is added.

Let N=(k+1)^2 and omega=(k+1)/a. Use real complete degree-k harmonics with
Haar normalization <Y_i Y_j>_S3=delta_ij. Thus sum_i Y_i^2=N; these are not
the repository's physically normalized Y_i/sqrt(V), V=2 pi^2 a^3.
Keep conversions explicit in code and reports.

## 3. Frozen predictions

### P1. The full multiplet gives pointwise stationary radiation stress

For m=N and phi_i=s Y_i cos(omega t+delta), re-derive

    sum_i Y_i^2 = N,
    sum_i D_a Y_i D_b Y_i = N k(k+2) g_ab/(3a^2),
    sum_i Y_i D_a Y_i = 0,
    sum_i T_i^TF = 0,   sum_i T_0a = 0,
    rho = N s^2 omega^2/2,   p = rho/3.

The first three are pointwise addition-theorem identities. Establish them
from the multiplet representation and independently check the **full
improved stress**, not just the gradient tensor. Require constancy over
space and time, including both scalar turning points. Individual members
must show nonzero anisotropy at a stated nonzero-amplitude phase, so the
cancellation is actually tested. Also collapse the component fields into
one coherent field and show that its pointwise support generally fails.

### P2. Einstein normalization and the full kinetic matrix

The predicted required relations are

    Lambda = 3/(2a^2),
    s^2 = 3/[kappa N (k+1)^2] = 3/(kappa N^2),
    sum_i phi_i^2 = 3 cos^2(omega t+delta)/(kappa N),
    f := kappa F = 1 - cos^2(omega t+delta)/(2N).

For physically normalized harmonics the coefficient squared is V s^2;
no hidden factor of N, V or a may be absorbed into fitted amplitudes.
Re-derive all Einstein components with the chosen Lambda. Reversing
Lambda or changing all amplitudes by 10% must fail the background equations.

For more than one scalar, the Einstein-frame kinetic coefficient is a
**matrix**, not #293's one-field scalar coefficient:

    K_IJ = delta_IJ/f + kappa phi_I phi_J/(6 f^2).

Predict transverse eigenvalues 1/f and the field-radial eigenvalue 1/f^2,
with the all-zero-field limit K=I. Check the entire matrix, its eigenvalues
and F>0 throughout a period. Larger minimum f is a larger margin from a
vanishing gravitational coefficient; it does not prove stronger dynamical
stability or constitute a general ordering of 'regularity'.

### P2a. A scoped minimality statement

Test/prove that at least four real components are necessary **within the
common-phase, single-degree, exact pointwise-support ansatz**. Four are
sufficient at k=1. The proof target is: write phi_I=f_I(x) cos(omega t).
At the simultaneous field-zero phase, constant nonzero density requires
sum_I f_I^2=r^2>0 constant. Isotropy and the eigenfunction equation then
require sum_I D_a f_I D_b f_I=k(k+2)r^2 g_ab/(3a^2), a rank-three form.
Since sum_I f_I D_a f_I=0, these three independent derivative vectors lie
in the (m-1)-dimensional tangent space perpendicular to f. Thus m>=4.

Do not label N=(k+1)^2 minimal for each k, or four minimal over arbitrary
relative phases, mixed degrees, interactions or other field actions.
At fixed k, lower-rank Gram matrices could furnish smaller constructions;
their existence is a separate question unless an explicit witness or
lower-bound proof is supplied.

## 4. P3: preparation sensitivity, with the full kernel retained

Run the requested diagonal amplitude test, and the larger basis-independent
test which determines what that diagonal test misses. Neither outcome is
frozen as positive or negative.

For common phase, write the component vector as

    phi(t,x)=C Y(x) cos(omega t),   G=C^T C >= 0.

G has N mode indices; C has m component rows and N mode columns. Internal
orthogonal rotations C -> O C leave G and stress exactly unchanged. Harmonic
basis changes transform G and the stress operator together. Report these
redundancies explicitly; a field-component count is not a harmonic-mode count.

The reference is G0=s^2 I_N. At fixed trace(G)=N s^2, perturb

    G(epsilon)=s^2 (I_N+epsilon H),
    H=H^T, tr(H)=0, ||H||_F/sqrt(N)=1.

Require positive semidefiniteness before constructing component fields
through a real square root or factorization. A simple sufficient interior
condition is |epsilon| ||H||_op < 1. Relative Frobenius Gram mismatch is
then exactly |epsilon|. No arbitrary normalization may be optimized to
reduce the measured stress after the preparation is chosen.

The improved stress is linear in G for this ansatz. Build two operators:

1. A_k: H -> spatial trace-free stress change only.
2. L_k: H -> **all** physical stress-component changes relative to the
   reference ESU source, including density, momentum and pressure through
   all space and time. A_k's kernel alone is not a family of ESU supports.

Use the output norm

    ||delta T||^2 = <delta rho^2 + 2 |delta j|^2
                       + ||delta T_spatial||_F^2>_(S3, one period) / rho0^2.

For A_k use the same normalization with only the trace-free spatial block.
These integrals define comparison norms; acceptance of a zero requires
the pointwise polynomial/Fourier coefficient identities, not an average
anisotropy that cancels across the sphere. Use the coefficient basis
{1, cos(2 omega t), sin(2 omega t)} and complete spatial polynomial content.

For k=1,3,5, report ranks, nullities, singular values, smallest positive
singular value and largest singular value under the fixed input/output
norms. Include the trace/amplitude direction separately. Distinguish
preparation redundancies from nonzero changes of G invisible to stress.

Certify kernel/rank claims with exact polynomial algebra. A rational harmonic
basis and its rational moment matrix can be used for rank and nullspace;
whitening for singular values may be numerical. A modular rank certificate
is only a rational-rank lower bound unless paired with a matching exact
kernel upper bound. A small singular value alone is UNRESOLVED, not proof
of an exact flat direction or uniqueness.

For an exact H in ker(L_k), linearity in G makes the whole PSD interval
G0+epsilon s^2 H an **exact** equal-stress family. Verify the factorized
fields' full stress independently. For nonkernel H, measure anisotropy
and full-source mismatch per unit relative Gram mismatch. Check the largest
and smallest-positive singular directions and 20 seeded random directions.

The diagonal control perturbs amplitudes by 1+epsilon d_i, sum_i d_i=0,
sqrt(mean(d_i^2))=1, then rescales the entire vector only to keep sum_i s_i^2
fixed at its prescribed value. Report this amplitude mismatch separately
from the induced Gram mismatch (whose first-order diagonal change is 2d_i).
Use epsilon=+/-0.01, +/-0.005, +/-0.0025, subject to the PSD gate. Include
explicit nonkernel negative controls and every certified kernel witness.

## 5. Meaning of sensitivity and handoff to round 13

If the fixed-trace kernel is trivial, the conclusion is isolated preparation
within the tested ansatz, with a measured finite susceptibility. Exact
isotropy requiring exact equality does not by itself establish infinite
fine-tuning, unstable evolution or improbability under an unspecified
preparation measure. If the kernel is nontrivial, characterize the exact
equal-stress family; this still does not provide a mechanism selecting it.

The amplitude experiments evolve free scalar phases on the fixed round ESU.
Perturbed sources which fail Einstein are *source-response controls*, not
self-consistent perturbed universes. This round does not equate those
controls with a coupled dynamical stability or attraction calculation.
Relative-phase robustness is not established by a common-phase Gram test.

Even a successful background does not automatically repair the perturbation
model behind #289--#292. Its field gradients can contribute anisotropic
stress when the metric and component fields respond. Reusing the bare
omega_T^2=8/a^2 as a coupled normal frequency requires a new derivation.
Likewise, using the same components for support and signal gives linear
cross terms; adding a distinct signal component is another field-content
choice. Record `TT_frequency_transfer: NOT_ESTABLISHED` and
`coupled_support_response: NOT_DERIVED` until those equations are derived.

Round 13 can reuse the harmonic and component-sum machinery, but must keep
'several degrees of one scalar' distinct from 'several independent scalar
components'. Its frequency restrictions must name which support response
is retained. O(s^4) persistence remains downstream of that response, not
merely downstream of finding a positive background.

## 6. Verification and separate verdicts

Construction controls: k=1,3,5; a=0.7,1,2; kappa=0.4,1; common-phase points
including 0, pi/4, pi/2, 3pi/4 and pi, plus a full-period scan. Use two
independent spatial rules sufficient for the polynomial degrees, e.g.
Hopf rules (8,16) and (12,24), and off-grid points. Derive the identities
symbolically and evaluate full stress through the inherited implementation
at a=1 and an explicitly radius-aware route for the other radii.

Require symbolic residuals to vanish exactly; full-stress and Einstein
errors <1e-10 normalized by rho0 or a^-2 as appropriate, with absolute
residuals and floors stated. Require direct-field reconstruction of each
Gram response <1e-9 normalized. Report measured mismatch scaling and
singular-value conditioning; do not substitute a convergence claim for an
exact rank/kernel certificate. Freeze the singular-value reporting cutoff
at 1e-10 times the largest singular value; it labels numerical candidates
for certification, never the final exact rank.

Keep separate outcomes for background existence, kinetic regularity,
scoped component-count bound, diagonal sensitivity, full preparation kernel,
coupled dynamical stability, field-content/preparation selection, tensor
frequency transfer, Phi selection and the causality gate. Failed numerical
witnesses yield UNRESOLVED, not nonexistence. Missing/failed named required
checks must invalidate the affected physical verdict and cause the CLI to
overwrite stale reports and exit nonzero. Exercise every such failure path.

No result in this freeze derives the multiplicity, common phase, chosen
configuration, preparation measure, triangle map or source-local apparatus.
A completed construction would settle existence in an explicit expanded
matter class; it would not settle those remaining selection/response issues.

Deliverables after this freeze: an isolated module, an independently checked
probe, targeted tests, exact kernel certificates, a derivation and archived
JSON/Markdown. Preserve this file and distinguish re-derived prior candidate
results from P3 outcomes first obtained after publication.

## 7. Review amendment, published before implementation

This appendix responds to the review of `a12cb2f`. Sections 1--6 above are
preserved verbatim. No P3 scan or implementation preceded this amendment.

P1 additionally requires **componentwise odd parity**, `Y_i(-x)+Y_i(x)=0`,
both as a polynomial identity and at quadrature and independent off-grid
points. An even-degree multiplet is a required discriminating control: its
isotropic stress must not qualify as odd support.

P2a uses two phases. Constant density at the simultaneous field-zero phase
fixes `sum f_I^2` constant in space. The subsequent spatial-isotropy step
uses a phase with nonzero cosine. At the zero-field phase the trace-free
stress vanishes without constraining the gradient tensor.

The required gate identifiers and their affected verdicts are now fixed.
In the table B=`background_existence`, K=`kinetic_regularity`,
M=`component_count_bound`, D=`diagonal_sensitivity`,
P=`full_preparation_kernel`. Every listed gate has at least one target;
unknown/missing targets or a changed gate set invalidate all five verdicts.
A failed check invalidates its listed targets even if other checks pass.

| Required gate identifier | Affected verdicts |
|---|---|
| `harmonic_parity` | B, M, D, P |
| `addition_identities` | B, K, M, D, P |
| `full_stress_background` | B, K, M, D, P |
| `einstein_normalization` | B, K, M, D, P |
| `kinetic_matrix` | K |
| `component_bound` | M |
| `independent_field_controls` | B, M, D, P |
| `exact_preparation_kernel` | P |
| `normalized_sensitivity` | D, P |
| `gram_reconstruction` | D, P |
| `scope_and_order` | B, K, M, D, P |
| `failure_paths` | B, K, M, D, P |

The tuple `REQUIRED_CHECKS` must contain exactly these twelve names, and
`CHECK_TARGETS` must be total with exactly these targets. Tests must compare
both against this published table and exercise each gate as missing and
false, checking every affected verdict, nonzero CLI status and replacement
of stale JSON/Markdown. Empty or unrelated-only check dictionaries cannot
pass. Missing/invalid kernel results cannot become an affirmative physical
classification through a gate set alone.

Fixed scope verdicts remain `coupled_dynamical_stability: NOT_ESTABLISHED`,
`field_content_preparation_selection: NOT_DERIVED`,
`TT_frequency_transfer: NOT_ESTABLISHED`,
`coupled_support_response: NOT_DERIVED`, `Phi_selection: NOT_DERIVED`, and
`causality_gate: OPEN`. No gate result promotes these statements.

The coherent-field control is a structural regression of #293, not a new
exclusion. P2 should explicitly compare its density and Lambda with #293's
homogeneous even control; their equality does not imply equal perturbation
response or select either preparation.
