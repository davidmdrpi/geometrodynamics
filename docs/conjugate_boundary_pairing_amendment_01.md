# Amendment 01: oracle degree, influence and a transmitting carrier

Date: **2026-11-08 (America/Los_Angeles)**. Applies before Stage 1.

This responds to the [review of PR #351 at 8b13f3d](https://github.com/davidmdrpi/geometrodynamics/pull/351#issuecomment-6734211935).
The [original protocol](conjugate_boundary_pairing_protocol.md), published at
836bb278, and its [baseline evidence](conjugate_boundary_pairing_baseline.md)
remain unchanged. This amendment adds prospective tests; no responding-mouth
trajectory or oracle experiment has been run. The spatial/time discretization
and exterior implementation still require the pre-run specification demanded
by Stage 1. This amendment does not close that requirement.

## 1. Source and theorem convention

Add David Turturean, [A proof of the Aaronson–Ambainis conjecture,
arXiv:2611.09773v1](https://arxiv.org/html/2611.09773v1), submitted November 6,
2026. This is a preprint with author-reported Lean kernel verification.
We have read its statement, not independently replayed the formal proof or
certified its correspondence with the printed conjecture. The
[formalization's scope statement](https://github.com/physical-superintelligence/aaronson-ambainis#verification)
distinguishes those tasks and excludes the query-complexity consequences from
its formal coverage. Retain the original Aaronson–Ambainis reference.

For uniform x in {-1,1}^n, use exactly the convention of Theorem 1.1:

    f_hat(S) = 2^-n sum_x f(x) product_{i in S} x_i,
    Var(f) = sum_{S != empty} f_hat(S)^2,
    Inf_j(f) = E[((f(x)-f(x^j))/2)^2]
             = sum_{S containing j} f_hat(S)^2.

If |f|<=1 and degree(f)<=d, d>=1, the stated bound is

    max_j Inf_j(f) >= c Var(f)^(5/2)/d^6,
    c = sqrt(15)/(500 pi^4).

For a [0,1] energy fraction p, use f=2p-1 throughout. Then
Var(f)=4 Var(p) and Inf_j(f)=E[(p(x)-p(x^j))^2]. This avoids mixing the
half-difference and full-difference influence constants. The polynomial-method
degree bound is a separate, older result; the AA bound does not establish it.

## 2. Define the first oracle comparison

Retain the original two-incidence phase grid as the physical response study.
Add a distinct small Boolean diagnostic with **n=6 and T=1**:

- Prepare six orthonormal incoming characteristic waveforms: three temporal
  modes at each of the two ports. Use the positive incoming-energy inner
  product, not the signed net boundary flux. On [0,3], Gram–Schmidt the ordered
  real seeds u(t)^m exp(-u(t)^2/2) cos(omega(t-0.6)), m=0,1,2,
  u(t)=(t-0.6)/0.15, separately at each port. Set them to zero outside the
  preparation interval. Freeze their normalization, sampling, endpoint
  treatment and convergence check before solving dynamics. These are incoming
  characteristic amplitudes; they are not prescribed Dirichlet field values.
- A single oracle layer is O_x=diag(x_1,...,x_6) on these prepared modes.
  Supply equal mode amplitudes A/sqrt(6), so incident energy is A^2 in this
  normalization for every bit string. This subexperiment has its own fixed
  energy normalization; it does not redefine A in the original phase grid.
- The sign layer acts once during preparation. There is no bit-controlled
  material parameter, source adaptation, stored input-dependent support
  state, or second oracle crossing. Subsequent propagation and nonlinear
  response are input-independent laws acting on the encoded field.
- Count T as calls to O_x, not throat crossings, solver time steps, or focal
  events. Multiple propagation/reflection events do not query a disconnected
  oracle. Any future return through O_x increments T and needs a new protocol.
- Evaluate **all 64 strings**, with equal weights, at A=1/64,1/32,1/16,1/8
  and both carrier frequencies. Include A=0 as a zero-source control, not a
  normalized probability. Use the original lapse and support parameter grid.

This makes the access model explicit: one simultaneous six-mode sign operation
resembles one quantum phase-oracle call. Loading six individual classical
bits costs six ordinary bit queries. Neither convention is silently exchanged
for the other. A six-bit enumeration establishes no asymptotic advantage.

For each run record n, T, incident modes, generated/resolved modes, spatial
unknowns, time steps, actual spectral bandwidth/tails, incident energy, initial
mechanical energy, support/pump work, observation time, and numerical/readout
precision. Report solver time and memory separately. The target additive
precision for the normalized readout is epsilon=1e-4; report the achieved
precision and the smaller tolerance actually needed for high-degree signals.
Bell angles remain a different input space, not Boolean oracle bits.

## 3. Readout and boundedness gate

Use W_A(x), the nonnegative energy collected at the right outgoing port during
[0,3], with left output, stored field energy and mechanical energy retained in
the full ledger. For the initial unpumped, zero-mechanical-energy comparator,
set E_ref(A)=A^2 and p_A=W_A/E_ref. Source work must be computed from incoming
characteristics, including all imposed boundary work. Positivity of the
complete comparator energy and energy balance must imply 0<=p_A<=1 for every
string. Boundedness must be checked over the full table; never clip it.

If pump or support work can supply additional energy, predeclare a finite
**input-independent** upper bound E_ref on all available energy, including
initial energy and supplied positive work, and prove/verify its adequacy.
Dividing each output by its own input-dependent work generally makes a ratio
with a different Fourier degree; that is not a harmless normalization.
An unbounded reservoir or an unavailable energy bound makes this test
INAPPLICABLE, rather than an AA counterexample. A bounded energy fraction
remains an observable, not an established detector acceptance probability.

## 4. Degree and amplitude-scaling gates

For T=1, input-independent linear propagation gives psi_1 linear in the six
bits and W_2 of degree at most 2. In the unpumped support comparator,

    psi=A psi_1 + A^3 psi_3 + O(A^5),
    q=A^2 q_2 + O(A^4),
    W_A=A^2 W_2 + A^4 W_4 + O(A^6).

Here q_2 is quadratic, psi_3 is at most cubic, and W_4 is at most quartic in
the bits. Reduce powers using x_i^2=1. Terms of degree above 2 are **permitted**
at order A^4; their coefficient need not be nonzero. Symmetry or dynamics can
cancel them. This is a conditional prediction to discriminate, not a fitted
target or a theorem that the nonlinear system has degree at most 4 exactly.

Define both quantities rather than calling both "Fourier weight":

    H_raw(A) = sqrt(sum_{|S|>2} W_hat_A(S)^2),
    Weight_raw(A) = H_raw(A)^2.

If a nonzero quartic excess term survives, H_raw scales as A^4 and
Weight_raw as A^8. With the fixed E_ref=A^2 normalization, H_p scales as A^2
and Weight_p as A^4. The excess norm relative to a nonzero leading raw
readout is O(A^2); its squared fraction is O(A^4). Thus the review's A^4
prediction applies to raw coefficients/norm, or normalized squared weight,
**not raw squared Fourier weight**.

The linear comparator's H_p must be <=1e-10. For the nonlinear calculation,
record the complete Walsh table, Parseval residual, H_raw, Weight_raw, H_p,
Weight_p, variance, and every influence at all three resolutions. Let eta be
the maximum full-table change in p between successive resolutions, taking
the larger of the two refinement changes and 1e-10. Require the changes to
decrease; eta is an empirical error indicator, not a rigorous error bound.
The normalized high-degree norm has at most the full L2 error, which is at
most the pointwise error. Require H_p>10 eta on at least three consecutive
nonzero amplitudes before fitting a power. Otherwise report UNRESOLVED_SCALING.

For resolved adjacent amplitudes, publish log2[H_raw(2A)/H_raw(A)]. A leading
quartic signature means both slopes on the smallest resolved triple lie in
[3.5,4.5] and remain there at the finest two resolutions. Publish all slopes;
do not discard a disagreeing larger amplitude. Below resolution report an
upper envelope, not exact absence. A signal resolved above degree 2 rejects
the registered one-query quantum-readout identification, conditional on the
access and readout definitions; it is not an experimental violation of QM.
Its physical explanation may be precisely the nonlinear classical response.

An analytic diagnostic control (not a BAM trajectory) is

    a_A(x)=A x_1 + A^3 x_2 x_3 x_4,
    W_A=a_A^2=A^2+2 A^4 x_1 x_2 x_3 x_4+A^6.

It has H_raw=2 A^4 and Weight_raw=4 A^8. For boundedness tests use the
input-independent denominator (A+A^3)^2. A linear amplitude supplies the
zero-excess control. These controls must be reproduced by the eventual
Walsh implementation before interpreting a numerical mouth experiment.

## 5. Influence gate and limits of the inference

First verify uniform input coverage, boundedness, access count and numerical
validity. Evaluate the stated AA lower bound at d=2T=2 only for a function
whose degree hypothesis is established. The fixed-background linear
comparator supplies that structural hypothesis. Numerically small Fourier
tails in the nonlinear run do not prove exact degree. Such a run receives
only a conditional numerical-consistency report, not a certified application
of the theorem. Alternatively d=n=6 always bounds the exact multilinear
degree of a six-bit table, but loses the proposed T-query interpretation.

Report max influence, variance, the bound, their difference, and propagated
uncertainty from the per-input readout errors. Declare a resolved consistency
pass only if its interval lies above zero; an overlapping interval is
UNRESOLVED_NUMERICALLY. A constant function is VACUOUS, not a success for the
mechanism. If a resolved discrepancy appears, audit degree, boundedness,
normalization, numerics, access model and theorem applicability before any
mathematical conclusion. Passing this weak necessary bound is not quantum
reconstruction. It is not a gate on CHSH or on the Born rule.

The preprint's consequence concerns approximate acceptance probabilities on
most uniformly distributed inputs and counts queries, not simulation time.
It does not turn an arbitrary classical nonlinear wave device into a quantum
algorithm. Structured promises remain important comparison cases; any BAM
claim of reproducing an oracle separation must expose the physical access
and resource cost rather than declaring a whole wave encounter one classical
bit query. No such complexity claim is made by this six-bit experiment.

## 6. Transmitting carrier and clocks

Add omega=60 alongside omega=20 to the original Stage-1 preparation grid.
Retain envelope width 0.15 and both lapse depths. Direct checks of the inherited
static solver at rtol=2e-12 give the following finite-port energy fractions:

| Lapse depth | Transmission at 20 | Transmission at 60 | Optical crossing time |
| --- | ---: | ---: | ---: |
| 0 | 0.0000567872153 | 0.8518771660 | 0.6621766691 |
| 0.4 | 0.0050303400 | 0.9999934036 | 0.8380207687 |

The clock is t at the rim; crossing time is integral ds/N for stretch=1.
These are carrier checks, not packet transmission or pairing measurements.
Report incident and outgoing spectra, integrated energy, spectral centroid
and bandwidth for each carrier and for any generated sidebands, with a
window/transform convention fixed before the run. Resolve the Gaussian tails
and oracle-mode bandwidth rather than treating the carrier as monochromatic.

Finite observation time bounds propagation opportunities only after geometry,
boundary access and a positive travel-time bound are specified. It does not
bound nonlinear polynomial degree by 2T. Repeated nonlinear mixing can
generate all six-bit Fourier levels without another oracle call. Also, on a
fixed finite Boolean cube even an infinite return resolvent reduces to a
multilinear function of degree at most n wherever defined; what is lost is a
degree bound controlled by a finite traversal count uniformly in n. The
original protocol's resolvent warning should be read in this sense.

The [electrovacuum completion](einstein_throat_evolution.md) has central
contraction time pi*M=0.04053978 at angle 0.3 and zero charge. The supported
background's [0,3] window spans central proper time 3 or 1.8 for lapse depths
0 or 0.4, respectively: roughly 74 or 44 such times, not 74 in both clocks.
These are different Lorentzian completions, not a common collapse prediction.
Retain EFFECTIVE_SUPPORT_ASSUMED; neither a positive oscillator stiffness nor
a long-lived numerical solution supplies the missing GR support dynamics.

The constructive discriminator is now sharper: can a physically matched
mouth response retain a common quadratic operational law as its geometry
responds, or do resolved higher-degree terms identify the regime in which
that proposed quantum correspondence breaks down? Either result guides the
next dynamical derivation. The static baseline and Stages 3–4 are unaffected.
