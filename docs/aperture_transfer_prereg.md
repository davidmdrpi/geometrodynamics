# Prospective finite-aperture capture and transfer experiment

Date: 2026-10-09 UTC. Parent: PR #322, commit
`b4aed67693a26e66882143a9507e9070b79ec59f`.
Publish this specification, implementation and producer before measuring the
validation schedule. No validation propagation pilot has been run. Pre-freeze
checks cover exact round refocusing, cap quadrature/tails, matrix identities,
and a small lmax=3 implementation test; they are not validation evidence.

## Question and scope

Does capture and a measurable shifted first-arrival return survive finite
aperture, nonround compact geometry and two temporal bandwidths, with an
explicit field/port energy ledger? A fixed-delay graph cannot answer this.

Resolve a real scalar field on ultrastatic R x Berger S3. This is the
three-dimensional compact wave layer, not the repository's four-spatial-
dimensional S4 bulk completion. Prescribed static geometry, scalar coupling
and distributed transducers are assumptions. No Einstein evolution, support
stress, gravitational mouth motion, topology surgery or excised boundary
matching is supplied. The ports are compact spatial profiles coupled to
one semi-infinite scalar lead each. They are NOT absorbing holes and need
not collect every wave crossing their support. Their coherent coupling,
including cancellations, is part of the tested action.

Family 344's Blaschke rigidity motivates testing imperfect focusing; its
sphere special case is classical. Family 365 motivates operator-valued
input/output measurements, but its smooth boundary-recovery theorem is not
an input here and does not identify this model from finite data. No result
from the new mathematical preprints is presumed verified or used as a lemma.

## Geometry, aperture and action

Use R=1/pi and c=1, so a round antipodal geodesic takes time one.
Write the Berger metric as R^2 times the round metric with Hopf fiber lengths
multiplied by b. Its scalar curvature is 2(4-b^2)/R^2. On a spin-l/2 SU(2)
representation the scalar Laplacian eigenvalues are

    R^-2 [l(l+2) + (b^-2 - 1) m^2],   m=-l,-l+2,...,l,

with multiplicity l+1 for each m. This follows by replacing the fiber part
of the round Casimir by b^-2 times that part. The wave operator is
-delta_g + xi*Scal_g, with primary xi=1/6. This conformal scalar assumption
makes round frequencies pi*(l+1) exactly; minimal xi=0 is an explicit
control, not silently equated to this theory. All tested operators are
nonnegative. The minimal constant mode is retained.

At identity A and minus identity B use radial compact profiles
p(theta)=(1-(theta/a)^2)^4 for round angular distance theta<a, zero outside,
normalized in physical L2(dV_b). Aperture angle a is a round-coordinate
footprint, not a fixed physical geodesic radius in every Berger metric.
Physical volume changes with b; normalization is declared and not fitted.
Port strength gamma=8 is fixed for all cases. It is not derived from GR.

Let w_A,w_B be these normalized functions. Impose at each lead endpoint
psi_j(0)=sqrt(gamma)*<w_j,phi>_L2. Vary the bulk action

    1/2 integral [phi_t^2 - |grad phi|^2 - xi*Scal*phi^2] dt dV

plus the two ordinary scalar lead actions. In orthonormal spectral coordinates,
with K the wave-frequency squares and W the two profile coefficient columns,

    q'' + Kq + W W^T q' = 2 W a_in,
    b_out = W^T q' - a_in,
    dE_bulk/dt = |a_in|^2 - |b_out|^2,
    E_bulk = (|q'|^2 + q^T K q)/2.

The port-observable subspace in each exact spectral block is compressed by
its 2x2 real Gram matrix. This is an exact reduction for these two ports and
zero initial field, not a fit of the bulk to a delay. The full retained
harmonic space has sum_(l<=L)(l+1)^2 modes; unexcited dark modes are omitted.
For displaced B=-exp(i delta sigma_axis), the cross Gram entries use Wigner
representation matrix elements. Both horizontal and fiber displacement are
tested. Character and norm identities are checked independently.

The source lead carries cos(2*pi*t)^4*cos(pi*w*t) on |t|<.25, zero otherwise.
A is the primary source, B has zero incoming field, and the bulk starts at
q=q'=0 at t=-.25. Reverse-source control independently drives B.
The outgoing B flux is the coherent captured energy. A reflection and the
remaining bulk energy are recorded; capture is not normalized by an inferred
or fitted 'successfully launched' fraction.

## Return and measurement conventions

Route the recorded B output through a unit-transmission handle with proper
transit h=.125, scalar sign -1 and prescribed clock offset Delta=1.5:

    t_return = t_B + h - Delta,    r(t_return) = -b_B(t_B).

Every recorded sample is relabelled once, including negative return times.
There is no FFT wrapping, interpolation, gain, lost shifted sample or second
copy of its energy. Delta=0 and a disconnected receiver are controls. This is
feed-forward first-arrival transport; the returned wave is NOT reinjected
into A. A new self-consistent feedback history with this measured transfer
operator remains NOT_ESTABLISHED. The clock preparation assumptions and work
budget of #322 are inherited; changing the spatial bulk does not solve them.

The first-arrival measurement window ends at t=1.75. Residual bulk energy is
retained in the ledger, not called a loss or assumed negligible. An extension
to 2.25 reports additional capture separately; it must reproduce the entire
original interval. Arrival mean and RMS spread are computed from all B flux
in the declared window, without adaptive peak picking.

Also report complex two-port scattering at omega=pi*(j+.37), j=1,...,18,
using exp(-i omega t). Check reciprocity and unitarity. The handle multiplies
the B output by -exp(i omega*(h-Delta)); its magnitude alone cannot identify
Delta, and its zero-frequency phase contains no time-shift information.
Only the given geometry/ports are being tested, not general inverse uniqueness.

## Frozen schedule and gates

35 cases, in the exact order in aperture_transfer_probe.schedule:

- b=(.9,1,1.1), a=(.4,.6), w=(8,12), each with (L=40,dt=1/512) and
  (L=56,dt=1/1024): 24 primary cases.
- At b=1.1,a=.4,w=12,L=56,dt=1/1024: disconnected B, reversed source,
  time refinement to 1/2048, mode refinement to 72, extension to 2.25,
  minimal coupling xi=0, round minimal coupling b=1,xi=0, and offsets
  delta=(.25,.5) separately in horizontal and fiber directions: 11 controls.

Numerical validity requires EVERY case and applicable comparison:

- cumulative discrete field-plus-lead energy error/input energy <1e-9;
- final archived state agrees with its energy to relative input scale <1e-9;
- independent reconstruction from recorded forces a_in-b_out agrees with
  every port sample <1e-8, every field energy/input <1e-8 and final state in
  unweighted (q,q') Euclidean norm <1e-8;
- cap truncation norm deficit <1e-4; no truncated-profile renormalization;
- premature B energy/input <1e-5 before the conservative characteristic
  bound -.25+min(1,b)*(pi-delta-2a)/pi;
- output-history L2 difference/source L2 <.03 for paired refinement,
  separate time/mode refinement and extension, and <1e-7 for port-reversed
  source after swapping output labels; interpolation onto coarse midpoints
  is only used for these comparisons, never for propagation or return;
- two-port frequency unitarity and reciprocity defects <1e-10.

After numerical validity, mechanism success requires ALL twelve primary
fine-grid cases to capture >1e-4 of input energy and return >1e-5 before
source onset. Each perturbed capture must retain >.1 of the matching round
capture. Disconnected capture and causal before-source return must be <1e-12.
At delta=.5 BOTH displacement axes must change the output history by >.01
in source-normalized L2. Minimal-coupling controls are reported separately;
their numerical validity is mandatory but conformal-model success cannot
be promoted to coupling-independent success.

All gates pass: FINITE_APERTURE_TRANSFER_SUPPORTED.
Numerically valid but any mechanism gate fails: FINITE_APERTURE_TRANSFER_FAILED.
Any numerical gate fails: NUMERICALLY_UNRESOLVED, with the failing case named.
No post-measurement threshold, geometry, strength, bandwidth or horizon changes.
No corrective kicks, timing fit or amplitude-dependent tuning is permitted.

## Evidence

Publish source hashes and a pre-measurement commit. Save every time grid,
incoming/outgoing waveform, field-energy history and final modal state in
lossless NPZ/base64. Pin all archives, metadata and results by a SHA-256
manifest. Replay reconstructs dynamics from archived port forces, recalculates
all diagnostics and decisions, and checks source and file fingerprints.

GR support, excised-mouth matching, gravitational recoil and closed-feedback
history remain NOT_ESTABLISHED regardless of the reduced verdict. Fine tuning
has not been eliminated: homogeneous Berger geometry, two selected aperture
profiles, conformal coupling, fixed lead strength and prescribed clock offset
are explicit inputs. Passing this finite parameter box proves no open-set or
arbitrary-geometry robustness theorem.

References: Berger spectra can be checked against Appendix A of Egidi,
Gittins, Habib and Peyerimhoff, J. Spectr. Theory 13 (2023), 1297-1343,
https://ems.press/content/serial-article-files/47035 . Motivating manuscripts:
https://github.com/openai/math/tree/main/preprints/The-metric-Blaschke-theorem-September-23-2026
and https://github.com/openai/math/tree/main/preprints/Determination-of-a-metric-and-a-unitary-connection-from-one-boundary-patch-October-5-2026 .
