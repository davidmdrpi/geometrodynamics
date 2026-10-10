# Half-clock symmetry and a held-out LRS 3/7 resonance

2026-10-10. Base: audit branch at `9015ad8`. This is an additive study;
the earlier breaking protocol, labels, sources and archives remain unchanged.

## Question and derivation

The earlier LRS 2/5 scan had a coherent harmonic-10 signal below its registered
resolution, contrary to a prior expectation of a detectable period-five chain.
Does an overlooked exact discrete symmetry account for the doubled harmonic?

In `r3_return_map.rhs`, the involution

    S(A,A',q,q',x,x') = (A,A',-q,-q',x,x')

commutes with the flow: the geometric equations depend on q squared and q q',
and the clock equation is odd in q. Let F take a downward q=0 crossing to the
next upward crossing. Then H=S F is a map on the original downward section,
and **P=H squared** wherever these crossings exist and are nondegenerate.
S preserves the canonical form (both clock position and momentum reverse).
This is a symmetry reduction; it does not assert that scalar sign is a gauge.

On a continued invariant circle, 2 omega_H=omega_P modulo 2 pi. The branch
must be measured, rather than assumed. Pre-freeze checks on the existing
circles a=.004 and .128, at phase .37, select

    rho_H = (rho_P + 1)/2.

The wrong branch misses by .0218 and .702 coordinate units respectively;
the right branch agrees within 2.1e-12; H squared and the full matrix P agree
within 7.2e-13. These are retrospective observations, not held-out evidence.

Thus, for coprime p,q, the resonance denominator of H is

    Q = 2q / gcd(p+q, 2q).

At P rotation 2/5, H rotation is **7/10**, and the first permitted angular
harmonic in the resonant normal form is 10. A fifth harmonic of P alone
does not survive the additional half-clock symmetry. At P rotation 3/7,
H rotation is **5/7**, and harmonic 7 is permitted. This removes the need
to invoke an extra integral to explain *why the fifth harmonic is absent*.

The statement concerns a normal-form angular coordinate. The archived
phase-sliced lambda uses an interpolated curve and a normalized Euclidean
symplectic-dual direction; these can introduce sidebands. Exact vanishing of
every nonmultiple of Q in that raw Fourier spectrum is not claimed.

An analytic resonant Hamiltonian can first have angle dependence at degree Q
in complex center coordinates (I^(Q/2)); its vector field can contain degree
Q-1 terms. This is order counting, not a computation of the coefficient or
a bound at the finite crossing amplitude. Ordinary normal-form suppression
is now a competing explanation for the tiny earlier signal. We do not predict
an absolute breaking amplitude, prove convergence of a normal form, identify
a new integral, or infer quantum action selection.

## Held-out prediction and failure

No map has been evaluated at the 3/7 crossing before this freeze. Read-only
inspection of the old continuation supplied the bracket
a=.18101933598375608 and .21526948230495083. The producer interpolates their
63-node circles linearly in omega to 6 pi/7, exactly as the earlier scan did.

**Prediction:** the raw obstruction spectrum has dominant harmonic **7**
among harmonics 1..21, its coefficient magnitude exceeds **10r**, and lambda
has **14 cyclic sign changes**. Here nu is the maximum, over the independent
checks below, of |lambda_alt-lambda_primary| plus the alternate residual,
and r=max(10nu,1e-11). Fourier amplitudes mean |FFT(lambda)|/84, without the
usual factor of two for real sine amplitude.

This is deliberately stronger than the symmetry theorem: symmetry permits
the coefficient but does not require it to be nonzero. If the numerical and
identity gates pass and any prediction condition fails, the registered label
is **SEVENTH_HARMONIC_PREDICTION_FAILED**, even if the signal is undetectably
small. Such a result leaves the exact square-root identity intact but defeats
the proposed detectable-leading-harmonic explanation at this crossing.

## Execution and gates

1. Retrospective identity checks at three existing circles (indices 0,20,24
   among successful archived circles), phases .37,1.21,2.43. Require H squared
   versus full-matrix P within 1e-9, correct circle lift within 1e-8, wrong
   lift at least 1e-3 away. Re-solve 12 of the old 2/5 scan points using H
   squared; these are descriptive reproduction checks, not a new detection.
2. At **84 uniform phases**, use the frozen `r3_breaking.scan_point`, with
   seven full-return nodes, P evaluated as two scalar half-clock integrations
   (DOP853 1e-13/1e-15), Newton tolerance 2e-13, maximum 12 iterations,
   centered Jacobian step 1e-7. Recompute the final residual at returned
   nodes, since the historical scanner's last history entry can precede
   its final update. Require every phase residual <=2e-12, maximum constraint
   <=1e-10 and phase-system condition number <=1e8.
3. At phases j=1,13,25,37,49,61,73, chord re-solve using the **full matrix
   equations** and direct full-period event detection with both DOP853
   1e-13/1e-15 and Radau 1e-12/1e-14. Six chord iterations maximum; require
   nu<=1e-10 and every check present. No averaging or missing-point FFT.
4. Numerical failure gives NUMERICALLY_UNRESOLVED. Valid numerics but failed
   identity checks give HALF_CLOCK_IDENTITY_FAILED. Otherwise apply the
   prediction's support/failure label above.

The coordinate norm is the maximum component norm on (A,p_A,x,p_x), in the
existing kappa=a_background=1 conventions, except the scanner's explicitly
Euclidean normalization of its tangent and dual. Lambda is a coordinate
obstruction, not an energy or invariant physical measure of breaking.

This is a multiple-shooting closure experiment, with seven independent
section nodes fine-tuned by Newton. It does not evolve generic initial data
for seven periods or stabilize the ESU hyperbolic direction. No corrective
kicks or controller are applied to trajectories. Integration stops at the
next clock crossing (half-map cap 2.2 conformal-time units); no long-horizon
stability inference is made.

Publish protocol, code and tests before running production. Inputs and
source hashes are bound in each JSON; output creation is exclusive. There
are no automatic retries with changed settings. Amendments or extra analyses
must be separately dated and described as post hoc. No old archive is edited.

## Scope and the next decision

If supported, prioritize the half-clock resonant normal form and its first
nonzero coefficient over fitting an unspecified extra integral. The earlier
2/5 result would then be less surprising than its full-map denominator
suggested. The circular #319 family needs a separate analysis of its
half-clock lift and D3 action; neither a new angular selection rule nor exact
integrability for that family is claimed here. This experiment does not yet
test evolving-background wave transport or any inhomogeneous perturbation.
