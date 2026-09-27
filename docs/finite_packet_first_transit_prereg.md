# Finite tensor packet: first antipodal transit (prospective freeze)

Date: 2026-09-27 UTC. Baseline: main at
`99c94fee6dc6f5158c1ccd8e791ccd51f52e3b5a` (#310).
Publish this file before implementing or measuring packet evolution.
Do not edit this freeze afterward; retain failures and use dated addenda.

## Question and scope

Can a finite-band, spatially localized TT metric disturbance deliver a
localized electric-Weyl signal near its first antipodal transit on the
breathing Einstein-frame four-scalar ESU of #310?
This first stage is linear in packet amplitude on the exact supported
background. It includes spatial reconstruction and an independently measured
linear instability budget. It does not implement nonlinear constraint
completion, a dynamical receiver, two-object recoil, or action quantization.
The next nonlinear experiment must retain scalar modes sourced at second
order. No result here certifies that later experiment.

Trapped neck spheres in #309 do not prove an event horizon or forbid all
neck routes. No neck is needed or evolved in this experiment. Do not rerun
#308. The four conformal scalars are retained assumptions, not vacuum GR.

## Explicit spatial carrier (derived algebraically before publication)

On unit S3, gamma = dchi^2 + sin(chi)^2 q_AB dx^A dx^B. Use the real
unit-normalized S2 harmonic Y_20, L=6, and its tracefree Hessian
Z_AB = D_A D_B Y + 3 q_AB Y. For each n=2,...,80 set

    A_n = C_(n-2)^3(cos chi)
    B_n = sin(chi)^2 [A_n' + 3 cot(chi) A_n]/6
    C_n = -sin(chi)^2 A_n/2
    D_n = sin(chi)^2 [B_n' + 2 cot(chi) B_n - A_n/2]/2
    H_chichi = A_n Y
    H_chiA = B_n D_A Y
    H_AB = C_n q_AB Y + D_n Z_AB.

The trace and divergence vanish; the radial eigen-equation is
A''+6 cot(chi) A'+[n(n+2)-8]A=0. Verify independently that the full tensor
rough Laplacian is -(n(n+2)-2)H, using coordinate covariant derivatives
for n=2,3,5 at nonsingular points. Normalize H in the S3 L2 norm. The
angular-integrated contraction of two such tensors is

    3 A_n A_m/2 + 12 B_n B_m/sin(chi)^2
                      + 12 D_n D_m/sin(chi)^4.

This constructs one fixed axisymmetric polarization, not all tensor states.
The family is smooth at the poles despite polar-coordinate singularities.
Its antipodal pullback parity is (-1)^n.

## Preparation and evolution

Use n=2..80; Gaussian spectral windows (center,width)=(12,3),(24,6),(40,10).
Unnormalized coefficients are exp(-(n-center)^2/(2 width^2)) times
A_n(0)/||H_n||. Normalize coefficients to sum c_n^2=1.
Initial conditions are h_n=epsilon c_n cos phase,
p_n=f(0)(n+1)epsilon c_n sin phase, phases 0 and pi/4,
epsilon=1e-6,1e-5,1e-4. Evolve using #310's supported tensor equation:

    h_n' = p_n/f
    p_n' = -[n(n+2) f + 2 R^2] h_n,
    R=sqrt(3)/2 cos(2 eta), f=1-R^2/6.

Record eta in [0,pi+.15] at 801 equally spaced samples, explicitly including
eta=0 and pi. Also sample pi+linspace(-.15,.15,301) to locate the arrival
peak. No phase, spectral window, aperture or criterion is tuned after a run.

Run both all-degree packets on the S3 cover and even-n packets. Even-n
metric data descend to RP3 and have paired antipodal lobes already at
eta=0. Their result is a recurrence of one quotient location; it cannot
establish a signal between two independently located antipodal objects.
The all-degree case is a cover-space control, not RP3-compatible data.

## Observables

Reconstruct h_ij and the electric Weyl tensor measured by background
comoving orthonormal observers:

    E_n = [n(n+2) h_n - h_n'']/(4 f).

Check this identity against a direct first-order coordinate Weyl calculation
on the conformally related static product metric for n=2,3. The background
Weyl tensor vanishes, so its first-order perturbation is gauge invariant.
E is a curvature observable, not an energy density or a completed receiver.

Also record the local mixed-index support-stress TT response
delta T^i_j(TT) = -R^2 h^i_j/f^2. This is a field-derived gravitational
coupling diagnostic; it supplies no extra mouth constitutive law.

Integrate h_ij h^ij and E_ij E^ij over the north cap chi<=.3, south cap
chi>=pi-.3, and equatorial belt |chi-pi/2|<=.15. Indices in these
contractions refer to the unit-sphere orthonormal basis; E includes its
physical f factor above. Archive full-sphere norms and fractions too.
Use cap-average versus belt-average squared-Weyl density (divide by the
respective volumes), never label it transferred energy.
Record south-cap tidal-power peak time within the fixed arrival window;
it is an envelope diagnostic and not a causal front. Finite-band initial
data have nonzero tails, which must be reported.

At pi, compare the full state (h_n,h_n'/(n+1)) with minus the antipodal
pullback of its initial state. Report relative norm error, normalized
overlap, and the analogous Weyl-profile mismatch. Compare arrival cap
fraction with the initial source-cap fraction. For even-n data use the
union of both caps for fractions, but still report each separately.

## Numerical validation and controls

G1: coordinate TT divergence, trace, rough-Laplacian and Weyl checks below
1e-9 normalized residual; normalized harmonic Gram matrix differs from I
by <1e-8 for all degrees. These checks must not reuse the radial eigenvalue
identity as a replacement for coordinate differentiation.

G2: primary DOP853 rtol=1e-11 atol=1e-13 in unit-amplitude states,
independent Radau rtol=1e-10 atol=1e-12, same times and initial conditions.
Maximum relative scaled-state difference <1e-7. Use unit-amplitude
solutions then scale: amplitude scaling is a linear consistency control,
not evidence about nonlinear dynamics.

G3: radial Gauss-Legendre quadrature on each separate integration interval,
256 and 512 nodes, differs by <1e-7 relative to full-sphere power;
global orthonormalization uses 512 nodes. There is no spectral-cutoff
convergence claim: the n<=80 preparation itself is fixed.

G4: artificial free conformal dispersion h''+(n+1)^2 h=0 reproduces
inverted antipodal phase-space refocusing to <1e-9. A bare static tensor
h''+n(n+2)h=0 is a diagnostic control. These are controls, not alternative
solutions of the retained supported action. Zero-amplitude h,E,stress
must vanish, and amplitude-scaled squared powers must scale as epsilon^2.

All gates must pass for a certified packet result. Failure => UNRESOLVED.
For each preparation separately define LOCALIZED_FIRST_TRANSIT if:
(a) initial source-cap metric fraction >=.5;
(b) relative full-state antipodal error <=.1;
(c) arrival Weyl target-cap fraction / initial Weyl source-cap fraction >=.9;
(d) target-cap / equatorial-belt mean squared-Weyl density >=10;
(e) south-cap peak time differs from pi by <=.05.
If a numerical gate passes but a physical criterion fails, report
NOT_ESTABLISHED with the failed criteria. Even-n results are named
PAIRED_RECURRENCE, never two-location delivery. No unconditional prediction
is made for these cutoffs; #310 makes approximate high-frequency
refocusing plausible, but localization and finite-band distortion are new.

## Instability budget

Independently integrate #310's scalar n=2 fundamental matrix with DOP853
and Radau; report spectral radius and maximum singular gain in the stated
(here alpha,alpha'/3,beta,beta'/3) coordinates, and the transfer norm into
Newtonian Phi/Psi through the first transit. Normalize the scalar zonal
Y_2=sin(3chi)/(sqrt(2) pi sin chi); max|Y_2|=3/(sqrt(2) pi).
Require DOP853/Radau scaled-map agreement <1e-7.

For the known homogeneous growing branch report gain exp(sqrt(2) eta)
and the seed bound to keep fractional scale-factor drift <=.01.
For n=2 report the seed bound to keep |Phi Y| and |Psi Y| <=.01 using
the induced transfer norm over the full first transit.
Report seeds 1e-8,1e-6,1e-4 and zero as a control. These are linear
sensitivity budgets in explicit coordinates; do not combine them into a
claimed nonlinear stability theorem or suppress the unstable modes by fiat.

## #310 replay follow-up in the same PR

Do not alter #310's frozen criteria, archived maps or historical verdicts.
Full replay must compare decisions obtained by scoring fresh measurements
with archived decisions exactly, as well as compare numerical evidence.
Cover the demonstrated jointly modified convergence-errors plus re-scored
result exploit. Extend the same protection to full extension replay.
Retain and explicitly authenticate legacy evidence when validator source
hashes change; do not rewrite historical source hashes to hide the change.

Archive preparation, spatial checks, primary/independent solutions, raw
quadratures and controls, derived decisions, source hashes, baseline and
this public freeze commit. Replay must rebuild decisions and reject joint
evidence/result tampering. Include a full fresh replay and regression tests.
