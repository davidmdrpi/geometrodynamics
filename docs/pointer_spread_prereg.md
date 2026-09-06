# Preregistration: finite pointer momentum and instrument-modified histories

Follow-up to PR #284 at `9aa8b0a`. Publish this file before writing or running
the new experiment. Predictions below are analytic expectations, not fits to
a spread sweep. The question is whether a finite-width classical pointer
erases the early-record distinction in a specified coupled history model.

## Scope and new assumptions

This is a dynamical extension of the triangle-family model, not a derived
BAM field apparatus. The rotor direction is assumed locally accessible to
the pointer. The frame-to-field map missing from #284 remains missing. A
Gaussian reference preparation is a smooth classical choice, not a claim
that BAM supplies a calibrated experimental preparation. All physical scales
below are specified in dimensionless rotor units, not fitted to laboratory
data. No Born weights, quantum uncertainty relation or Hilbert tensor product
is to be inserted.

Unlike simply sampling P in `pointer_kick`, the experiment must integrate
the reciprocal source force, record before the future boundary, and recompute
the weights of the resulting complete paths. It must not condition only on
a triangle labelled by the perturbed endpoint while discarding the earlier
path's parallel transport.

## 1. Hamiltonian, times and reference preparation

Let x be on S2, p tangent to x, and (Q,P) be a canonical pointer. Use

    H(t) = |p|^2/2 + P^2/(2M) + h(t) P (m.x),
    m = (1,2,0)/sqrt(5), M = 1,
    h(t) = (2g/tau) sin^2(pi t/tau), 0 <= t <= tau,
    h(t) = 0 otherwise; g = 1, tau = 1.

The source has round inertia one. In embedded coordinates,

    xdot = p,
    pdot = -|p|^2 x - h(t) P [m - (m.x)x],
    Qdot = P/M + h(t)(m.x), Pdot = 0.

Read Y=Q(tau) at t=1. The future closure boundary is at t_B=2. After the
pulse the source moves freely to t_B, and the early record is not updated
or postselected using a later pointer value.

The reference source measure is uniform sphere area for x0 and a centered
isotropic tangent Gaussian for p0 with standard deviation sigma_s=0.1 in
each orthonormal component. The pointer reference is independent
P~Normal(0,sigma_P^2), Q0~Normal(0,0.15^2). Thus the main experiment has
nonzero source momentum width as well as nonzero pointer widths. Exactly
zero widths are controls, not descriptions of a realistic apparatus.

Freeze sigma_P = 0, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5. The primary
nonzero-width point is sigma_P=0.1, equal to the source's component width.
Fix a=e_z and future choices b0=e_x, b1=e_y. Use the same preparation,
pulse, readout axis, threshold and resolution at both settings.

## 2. Closure of the augmented path

Integrate the source's actual parallel transport C from x0 to x_B along
its Hamiltonian trajectory. With quaternion multiplication acting by
Ad_C on vectors,

    Cdot = (1/2) [0, x cross p] C, C(0)=1.

Append geodesic analyzer legs x_B -> u -> w -> x0, where u=s_A a and
w=-s_B b. The complete closed path has

    G_s = lift(w,x0) lift(u,w) lift(x_B,u) C,
    V_s = ||Ad_G_s-I||_F^2/16 = (1-G_s,0^2)/2,
    W_b = (1/4) sum_s exp(-beta V_s), beta=64.

Normalize quaternion expressions if needed for roundoff; do not discard
the path transport C. The common Pin loop minus sign is invisible to this
frame energy. Equal reference sector coefficients are still assumed.
The final leg returns to x0, so the source evolution is part of the closed
itinerary, not omitted from its holonomy. Geodesic antipodal punctures
remain undefined, of zero reference measure. Do not silently assign them
a preferred weight. The continuum parallel transport is approached by
refining the Hamiltonian integration.

The positive soft closure factor is an additional history-weight law. It
is not derived from the Hamiltonian's real-time evolution or from BAM.
We test its interaction with a specified local classical readout.

## 3. Two explicit preparation semantics

Let pi_s(dx0 dp0) be the source reference measure and nu_sigma(dP) the
pointer Gaussian. Let Z_b(P)=integral pi_s W_b(x0,p0,P).

**Primary: preserve the prepared pointer marginal.** Specify

    mu_b^fixed = nu_sigma(dP) pi_s(dx0 dp0) W_b / Z_b(P).

This defines the conditional source history law for each independently
prepared pointer momentum. Its P marginal is exactly the same Gaussian for
both future choices. Per-P normalization is an explicit choice, not a
consequence of Liouville measure. Q0 is independent and never enters W.

**Control: condition the entire joint reference.** Also calculate

    mu_b^joint = nu_sigma(dP) pi_s(dx0 dp0) W_b /
                 integral nu_sigma Z_b(P).

Its posterior P marginal is proportional to nu_sigma Z_b(P), and can be
different from the input Gaussian. Report posterior P moments, rather than
calling an unchanged reference distribution an unchanged actual preparation.

Both are precisely defined conditional history laws. Their realization by
an operational preparation and independent future-setting intervention in
BAM remains open. Preserving a marginal by definition does not derive a
laboratory mechanism for that law.

## 4. Readout, primary prediction and controls

Use B={|Y|>0.6}. Integrate the same Q0 noise analytically using Gaussian
tails around the calculated pointer shift. Do not give each setting a
different response kernel. Report both setting probabilities and their
signed contrast Delta=P(B|b1)-P(B|b0), early-record means/variances and
source/pointer diagnostics for both preparation semantics.

**P1: zero-width bridge.** With BOTH sigma_s=0 and sigma_P=0, the source is
stationary, C=1, Y=Q0+F(x0), and the construction must reduce to #284's
beta=64 calculation: tails about 0.1231586349 and 0.4897143360. Compare
pointwise loop energies with #283's triangle potential to 1e-11, and the
integrated probabilities with the independent whole-sphere quadrature to
3e-3 or better. This bridge does not assert that a moving zero-P source has
the same law as a stationary one.

**P2: primary nonzero spread.** Predict Delta>0.1 at sigma_P=0.1 for the
fixed-pointer-marginal law, with sigma_s=0.1. This is a robustness
hypothesis, not an assumed result. Freeze it before the sweep. Report a
failure if it does not hold. No monotonic trend, disappearance scale or
sign at the larger widths is preselected. Compute the joint-conditioned
control without requiring it to agree with the primary law.

**P3: remove the ensemble-update shortcut.** For each coupled trajectory
also evaluate its actual record using weights from the source path evolved
with the backreaction force removed (P=0 in the source equations). This
holds the old source-history posterior fixed while changing the readout.
Compare with the correctly reweighted law and report differences, even if
small or zero. Do not claim a nonzero difference by construction.

**P4: parity and preparation.** Simultaneous (x0,p0,P)->(-x0,-p0,-P), with
sector relabelling, preserves W and reverses the early shift. Means and
Gaussian-smoothed signs should remain zero and 1/2; this is not equality
of full laws across settings. Fixed-marginal posterior P moments must be
0 and sigma_P^2. A constant-response kernel must erase the distinction.
Check rotation covariance of the complete paths; with readout axis m=a,
the rotation taking b0 to b1 also leaves the instrument invariant, providing
a nondegenerate blind-readout control rather than assuming all projections
are informative.

**P5: finite force, not just additive noise.** Report source displacement
and momentum recoil relative to the zero-P trajectory. Verify nonzero
backreaction at nonzero P. Check sphere norm, tangent momentum and conserved
axial angular momentum (x cross p).m. At zero P, verify exact free-source
evolution; at nonzero P verify the equations against independent Hamilton
ODE integration, including the path quaternion and pointer record.
Check the work balance dH/dt=h'(t)P F and second-order convergence of the
chosen integration. No exact no-disturbance claim at finite spread.

## 5. Numerical design and decision rules

Use Strang splitting: exact free spherical geodesic drift plus its exact
parallel transport, alternating with exact canonical h(t_mid) P F kicks.
Use 64 pulse steps initially and exact free drift after the pulse. An
independent DOP853 ODE must check trajectory, record and transport, with
32/64/128-step convergence. Target component error below 1e-3 at 128 steps
on fixed moderate-momentum control states; failures must be investigated,
not renamed passes. Norm/tangency/axial-invariant residuals must be below
1e-9. Check quaternion products independently against the existing scalar
quaternion machinery and verify the full loop fixes x0.

Integrate P with 16-point standard-normal Gauss-Hermite quadrature. Integrate
the four source coordinates with scrambled Sobol points, 2^11 base points
plus explicit antipodal (x,p) partners. Use four fixed seeds 20260906 through
20260909. Reuse source samples across widths, settings and controls.
Report variation across scrambles; its standard error is a numerical
diagnostic, not a rigorous frequentist confidence interval.

At the primary width independently refine pulse steps to 128, Hermite nodes
to 32 and source base points to 2^12. The step refinement must move each
event probability by less than 1e-3; source and Hermite refinements must
move it by less than 1e-2. Require all four primary scrambles to satisfy
Delta>0.1, not just their mean. Also check Hermite refinement at sigma_P=0.5;
if wider-width results fail resolution, label them unresolved and increase
resolution under the same fixed model, without changing the primary width.
Report effective source sample sizes under the per-P weights.

For small sigma_P, the finite-beta, finite-time law should approach the
sigma_P=0 law continuously. Record the differences, without selecting a
power-law fit or declaring continuity of BAM from these numerical data.

Archive numerical values and all criterion statuses. The probe must exit
nonzero on any failed required criterion, and must not silently skip a failed
refinement or replace a failed primary prediction with a different width.
Implementation errors may be corrected with an explicit record. Preserve
this freeze and report any changes to the numerical design separately.

## 6. Interpretation and deliverables

A positive primary result would establish that smooth pointer momentum
spread does not automatically erase the record distinction in these stated
coupled classical history laws, including one with a fixed Gaussian pointer
marginal. It would close the exact-P=0 objection within this extension.
It would not establish an operational BAM signal without the missing field,
apparatus and intervention maps. A negative result for this one apparatus
would not prove universal dynamical non-readability. Either outcome is to be
reported at this scope; no program-wide verdict is preselected.

Deliver a reusable module, independent tests, an archived probe, a derivation
and assumption ledger, and a new PR. Update the status audit and the link
from #284. Do not modify mass fits or attempt composition in this experiment.
