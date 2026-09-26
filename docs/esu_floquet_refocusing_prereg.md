# Prospective Floquet spectrum and antipodal refocusing of the four-scalar ESU

Date: 2026-09-26. Parent: `main` at `02c93c3`.
Publish this freeze before implementing or measuring any period map below.
Retain failed outcomes. Corrections are dated notes; nothing here is edited
after publication.

## 1. Question and why it is asked now

#308–#309 closed the throat as a classical channel. The Einstein-frame null
energy condition holds (f>0), the evolved neck is future-trapped, and it
pinches on a timescale near its own size, roughly .03. The remaining
geometric mechanism for antipodal interaction needs no throat: on R x S3
every conformally coupled free field obeys

    phi(eta+pi, x) = -phi(eta, -x),

because the modes have omega_n=n+1 and antipodal parity (-1)^n. A disturbance
reassembles, inverted, at the antipode after conformal time pi. This is
conjugate-point focusing with the point-caustic (Gouy) sign -1. It is exact
for free conformal fields on a fixed ESU. The model's four scalars are
conformally coupled, but their background breathes, supports the geometry,
and couples to metric perturbations.

For linear perturbations of every harmonic degree about that background:

1. Is each period map elliptic (bounded) or hyperbolic (growing)?
2. After one antipodal transit, eta=pi, does a perturbation reappear inverted
   at the antipode? Is that exact, asymptotic in degree, or absent?

This is a classical, linear question. A refocusing result would not be a
quantum, a measurement event, or a discreteness mechanism. Integer spectra
on S3 are kinematic consequences of compactness.

## 2. Inherited background and conventions

Use the action of #294–#300: four real conformal scalars, kappa=8piG=1,
ESU radius a=1, Lambda=3/2. Work in the Einstein frame, which is equivalent
(g_E=f g_J with unchanged fields):

    S = integral sqrt(-g_E)[R/2 - G_AB dphi^A.dphi^B/2 - U],
    f = 1-|phi|^2/6,  G = I/f + phi phi^T/(6 f^2),  U = Lambda/f^2.

Background, with eta the Jordan conformal time (also the Einstein conformal
time):

    g_E = f(eta)(-deta^2 + gamma_S3),   phi = R(eta) x,
    R = q cos(2 eta),  q = sqrt(3)/2,   f = 1 - R^2/6,
    H = f'/(2f) = -R R'/(6 f),  R'' = -4R,  R'^2 = 4(q^2-R^2).

Here x is the unit-S3 embedding. f and R^2 have period pi/2; R has period pi
with R(eta+pi/2)=-R(eta). The antipodal transit time equals the field period
pi. Einstein proper time over it is integral of sqrt(f), about .968 pi.

Perturbations are classified under the diagonal spatial/internal SO(4). With
a tangent/radial split of the internal vector at each point, the field
perturbation is `delta phi = alpha Y x + beta grad Y + v`, where v is a
transverse tangent vector field. Metric perturbations split into scalar,
vector and TT parts. By SO(4) covariance, every harmonic of a given type
and degree n gives the same reduced equation; only its Laplacian eigenvalue
enters. Write k=n(n+2).

## 3. Reduced equations (derived before this freeze; see section 6)

### T: tensor (n>=2)

For TT h_ij (g_ij=f(gamma_ij+h_ij)), no field perturbation is sourced at
linear order. The fixed gradient stress `G_AB d_i phi d_j phi=(R^2/f)gamma_ij`
has a mixed-index TT response under metric variation, which gives a
background-supported mass:

    h'' + (f'/f) h' + [n(n+2) + 2R^2/f] h = 0.

At n=2 this is algebraically identical to #297's
`(M beta')' + K beta = 0`, `M=f`, `K=8+2R^2/3`, because
`8/f + 2R^2/(3f) = 8 + 2R^2/f`. Canonical coordinates: (h, f h').

### V: vector (n>=2)

Field `v = w(eta) V`, with V a transverse vector harmonic,
`(nabla^2+2)V = -lambda_n V`, `lambda_n=(n+1)^2-4`. The metric vector
perturbation is g_0i=f S_i, S=s(eta)V (the F_i=0 gauge). Then

    w'' + (n+1)^2 w = R s' + 2 R' s                      (field)
    s = 2(R w' - R' w)/(2 R^2 + lambda_n f)              (0i constraint)
    (f s)' = -2 R w                                      (ij, unused check)

Without s the field equation is exactly the free conformal one,
omega=n+1. Evolve (w,w'), with s from the constraint and s' from the ij
equation. n=1 is pure gauge plus the exact global rotation w=R: with
lambda_1=0, s=(w/R)' satisfies the field equation identically.

### S: scalar (n>=2), Newtonian gauge

    g = f[-(1+2 Phi Y) deta^2 + (1-2 Psi Y) gamma],
    delta phi = alpha Y x + beta grad Y.

    Phi = Psi - 2 R beta/f                                       (ij traceless)
    Psi' = -H Phi + [(R beta' - R' beta)/f + R' alpha/f^2]/2      (0i)
    0 = -R' alpha'/f^2 + (R R'/f) Psi' + (k R/f) beta
        - R(12f-7) alpha/f^3 - 2(f k - 12 f + 9) Psi/f
        + 3(8f-7) Phi/f                                          (00)
    alpha'' = -(R R'/(3f)) alpha' + 3 R' Psi' + R' Phi' + 2 k beta
              - [k - 6 + (25 f - 14)/f^2] alpha - 6 R Psi - 8 R Phi
    beta''  = -(k + 2R^2/f) beta + 2 alpha/f.

Use (00) to solve Psi algebraically. Its coefficient is
`(-2k f^2 + 24 f^2 + 6 f - 21)/f^2`, which is nonzero for n>=2 and
f in [7/8,1]. Phi' comes from differentiating the traceless relation. The
state is (alpha, alpha', beta, beta'). The trace ij equation
`2Psi'' - R'alpha'/f^2 - (2RR'/(3f))Psi' - (RR'/(3f))Phi' + (kR/f)beta
 - R(6f-7)alpha/f^3 - 2(4f-3)Psi/f - (8f-9)Phi/f = 0`
is not used for evolution and serves as a consistency check.

n=0 and n=1 are not re-derived here. #296 already reports them: the
homogeneous map over pi is hyperbolic, with multipliers
85.019695223207 and 0.011761980531 (the Eddington mode). The constrained
degree-1 cover block is -I. They enter only as cited context.

## 4. Observables

For each sector X in {T,V,S} and n=2..80, compute the fundamental matrix
M_n(pi) over eta in [0,pi], starting at eta=0, in the coordinates above.
Also compute M_n(pi/2).

**Stability.** Floquet multipliers mu of M_n(pi) (d=2 for T,V; d=4 for S).
- ELLIPTIC if max|mu| <= 1+1e-7.
- HYPERBOLIC if max|mu| >= 1+1e-5.
- MARGINAL otherwise.

Also report det M and whether the eigenvalues occur in reciprocal pairs.

**Refocusing.** Let P_n be the antipodal parity of the physical carrier:
- the field delta phi for S and V: (-1)^(n+1) for S, (-1)^n for V;
- the metric perturbation for T: (-1)^n.

Define

    Rf_n = P_n M_n(pi),   F_n = -tr(Rf_n)/d,
    D_n  = || Rf_n + I ||_2   in (u, u'/(n+1)) coordinates at eta=0.

Perfect inverted antipodal refocusing is Rf_n=-I, i.e. F_n=1 and D_n=0.
F_n is invariant under symplectic change of basis and time origin; D_n is
a secondary, basis-dependent defect. Also record the refocusing phase
theta_n = arccos(-tr(Rf_n)/2) for T and V. Report the harmonic-amplitude
convention (-1)^n M_n(pi) for S alongside, because the factor x or grad in
delta phi flips the sign relative to Y.

Per sector, give a refocusing classification:
- EXACT: D_n < 1e-8 for all n.
- ASYMPTOTIC: F_n increases toward 1 over n in [20,80], with a positive
  fitted exponent p in 1-F_n ~ C n^-p, and F_n >= .99 for n in [40,80].
- PLATEAU: F_n converges to a limit below .999 (the fitted p is consistent
  with 0 within the fit uncertainty).
- ABSENT: none of the above.

Report the odd-sector (RP3-compatible) subsets separately: T n even,
V n odd, S n even. In these, the fields are componentwise odd and the metric
is antipodally invariant, as inherited from #294.

## 5. Predictions and controls, stated before measurement

Leading-order WKB with a time-averaged mass, applied to the canonical
amplitude sqrt(f) h for T:

    T: theta_n ~ pi m_T^2/(2(n+1)),
       m_T^2 = <2R^2/f> - <H^2> - 1,
    V: theta_n ~ pi m_V^2/(2(n+1)),
       m_V^2 = <2R^2/f>,
    <2R^2/f> = 12(sqrt(8/7) - 1) = 0.8285428...

Averages are over eta in [0,pi]; <H^2> is evaluated by quadrature. Both
predict ASYMPTOTIC with 1-F_n ~ theta_n^2/2, i.e. p=2. Gate: the fitted
n*theta_n over n in [40,80] must match the prediction within 10%.
Otherwise the prediction fails and is reported as failed.

S: no quantitative prediction. A heuristic, non-binding expectation: the
2k beta and 2alpha/f couplings split the frequencies to roughly
n+1 +/- 1/sqrt(f). That would leave a degree-independent phase error near
pi(<f^-1/2>-1), about .1 rad, i.e. a PLATEAU near F~.995. Gravity terms
could change this. This expectation is not a gate.

Stability: no prediction. Every free frequency n+1 sits on a parametric
resonance of the pi/2-periodic pump. Whether the O(1) mass terms move the
low-n modes out of their tongues is not known.

Controls (all must pass before any verdict is issued):

- **C1.** Free conformal test field, omega=n+1: Rf_n = -I to 1e-10, and
  F_n=1.
- **C2.** V with s set to zero reduces to C1.
- **C3.** T at n=2: tr M_2(pi/2) = -0.0963065402 to 1e-8, reproducing
  #294's published map. Disclosed: squaring that published trace already
  gives tr M_2(pi) = -1.99073, so the n=2 tensor refocusing value is not a
  new prediction.
- **C4.** A bare static-ESU tensor (R^2/f and f' removed, omega^2=n(n+2))
  gives the analytic theta_n = pi(sqrt(n(n+2)) - (n+1)).

## 6. What was computed before this freeze (disclosure)

Symbolic derivations only, in the Einstein frame. A coordinate
linearization engine was checked against the exact breathing background:
zero residual in all ten Einstein components and four field equations.
Vector (toroidal h(chi) d_phi) and scalar (zonal Y(chi), Newtonian gauge)
reductions were derived from it:

- In V, constraint plus field equation imply (f s)' = -2 R w exactly.
- In S, with the equations of section 3, the time derivative of the
  (00)-solved Psi equals the (0i) value, and the unused trace equation
  vanishes. Both were checked exactly at random rational points (residual
  0 to 160 digits).
- T was derived covariantly and matched to #297 at n=2 by algebra.
- The WKB averages in section 5 are analytic.

No period map, multiplier, fidelity or phase was computed for any sector or
degree, except the inference from #294's published n=2 trace stated in C3.

## 7. Numerics and evidence gates

- **G1. Re-derivation.** Re-run the symbolic reductions from the committed
  engine. Additionally verify all three sectors at explicit harmonics:
  - zonal Y = sin((n+1)chi)/sin(chi) for n=2,3,4;
  - toroidal h = Gegenbauer C^(2)_(n-1)(cos chi) for n=2,3.

  All residuals must vanish symbolically, or be below 1e-12 at random
  points. For T, verify the covariant TT stress response on an explicit
  n=2 left-invariant TT tensor.
- **G2. Unused equations.** Along numerical solutions, the unused V ij
  equation and the S trace equation must hold to 1e-8 relative.
- **G3. Integrators.**
  - DOP853 at rtol/atol 1e-12/1e-14 is primary.
  - DOP853 at 1e-10/1e-12 is a secondary check.
  - Classical RK4 with 2^16 and 2^17 fixed steps is independent.
  - Traces must agree to 1e-7 between primary and 2^17-step RK4, and the
    RK4 traces must converge with a ratio in [8,32].
- **G4. Structure.** |det M_n(pi) - 1| < 1e-8, reported for all sectors.
- **G5. Controls.** C1–C4.

If a gate fails for a sector, that sector's verdicts are UNRESOLVED; other
sectors stand. Do not tune tolerances, degrees or thresholds after
measurement. Do not change the start phase eta=0; phase independence of
traces is a reported check (eta0 = .3 for n=2..10).

Archive all matrices, multipliers, fidelities, gate results, source hashes
and this freeze commit. Replay must recompute the maps and reject tampered
evidence. Tests cover the gates, controls, parity conventions and
failed-output withdrawal.

## 8. Verdicts

Reported separately:

    T_STABILITY, V_STABILITY, S_STABILITY:
        ELLIPTIC_2_TO_80 | HYPERBOLIC_AT[list] | MARGINAL_AT[list] | UNRESOLVED
    T_REFOCUSING, V_REFOCUSING, S_REFOCUSING:
        EXACT | ASYMPTOTIC | PLATEAU | ABSENT | UNRESOLVED
    T_WKB_PREDICTION, V_WKB_PREDICTION: PASS | FAIL | UNRESOLVED

The same are also given restricted to the odd (RP3) subsets.

Not tested here:
- nonlinear evolution;
- the homogeneous instability's interplay with inhomogeneous modes (the
  Eddington mode grows by 85 per transit and is excluded only at linear
  order);
- handles or mouths;
- momentum transfer between objects;
- quantization.

A refocusing verdict says what a later two-object experiment could look
for, and when: a response at the antipode at eta=pi. It does not establish
that such an exchange occurs.
