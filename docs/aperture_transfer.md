# Finite-aperture capture on a resolved compact geometry

The 35 archived cases show substantial first-arrival capture and a shifted
return in the declared Berger-S3 scalar model. Primary fine-grid capture is
**9.86%–50.81% of incident lead energy**. A +/-10% Hopf-fiber deformation
retains **34.59%–69.36%** of the matching round-sphere capture. Nonround
geometry substantially reduces capture; it does not eliminate it in this box.

There is an important diagnostic qualification. The **unaltered frozen
verdict is NUMERICALLY_UNRESOLVED**. One reverse-source diagnostic mistakenly
counted reflection at the driven B port as premature receiving-port arrival.
A separately labelled post-measurement correction checks arrival at A and
reconstructs the same archived dynamics. All numerical and mechanism gates
then pass, giving
**FINITE_APERTURE_TRANSFER_SUPPORTED_AFTER_DIAGNOSTIC_CORRECTION**.
This is not represented as an unqualified prospective pass.

## What is now resolved, and what is still assumed

PR #322 supplied a bulk edge with a prescribed travel time. This experiment
instead propagates a scalar field on an ultrastatic Berger S3, using its exact
Laplacian spectrum and spatially compact aperture profiles. The fine cutoff
retains 63,365 scalar harmonics (l<=56). Only their exactly equivalent
port-observable subspaces need evolution; unexcited dark modes are omitted.
No travel-time distribution is fitted to the answer.

This is the three-spatial-dimensional compact wave layer, not the earlier
S4 spatial bulk completion or an evolved GR solution. Berger deformation is
homogeneous anisotropy, not arbitrary inhomogeneous metric perturbation.
Apertures are distributed coherent transducers coupled to scalar leads, not
holes excised from the geometry. The port action is an additional physical
assumption. Captured flux means energy emitted down the receiving lead; it
is not all geometric-optics flux through a surface or a local energy snapshot.

Primary coupling is the conformal scalar xi=1/6, explicitly responsible for
exact round refocusing. Two minimal-scalar controls are also retained. Neither
coupling is derived from the BAM matter sector here. Aperture angles .4 and
.6 are round-coordinate footprints and are L2-normalized in the physical
Berger volume. They are not equal physical-radius mouths across the sweep.
The scalar lead strength gamma=8 and clock offset Delta=1.5 are fixed inputs.
There are no corrective kicks, timing fits or post-outcome parameter changes.

## Model and measured operator

The radius is 1/pi, so round antipodal geodesic transit is one. If b multiplies
Hopf fiber lengths, the squared wave frequencies are

\[
\omega_{lm}^2=\pi^2\left[l(l+2)+(b^{-2}-1)m^2+2\xi(4-b^2)\right],
\quad m=-l,-l+2,\ldots,l.
\]

The sphere scalar curvature is 2 pi^2(4-b^2). A compact radial aperture has
profile (1-(theta/a)^2)^4 for theta<a and vanishes outside. Character
quadrature supplies its harmonic coefficients without renormalizing the
truncated profile. Rotation matrix elements give displaced receiver profiles.

Variation of the scalar bulk action plus the two leads with endpoint
constraint psi_j(0)=sqrt(gamma)<w_j,phi> gives

\[
\ddot q+Kq+WW^T\dot q=2W a_{\rm in},\qquad
b_{\rm out}=W^T\dot q-a_{\rm in}.
\]

Consequently, with E=(|qdot|^2+q^T Kq)/2,

\[
E(t)-E(t_0)=\int_{t_0}^t
\left(|a_{\rm in}|^2-|b_{\rm out}|^2\right)ds.
\]

The source is cos(2 pi t)^4 cos(pi w t) on |t|<.25, zero otherwise.
All initial bulk coordinates and velocities are zero. The denominator for
capture is the complete incident lead energy, including any prompt reflection.
Remaining bulk energy is explicitly retained; it is not called a loss.

The code also computes the complete complex 2x2 frequency scattering matrix
at 18 preselected nonresonant frequencies per primary fine case. With
exp(-i omega t), the measured first-arrival B output is routed as

\[
t_{\rm return}=t_B+0.125-1.5,\qquad r(t_{\rm return})=-b_B(t_B).
\]

Every sample and its flux is retained exactly once. The handle changes phase
by -exp[i omega(0.125-Delta)]; a magnitude-only measurement cannot identify
Delta, and zero frequency carries no time-shift information. That phase is
an imposed clock map, not recovered from geometry in this experiment.

**The returned packet is not reinjected into A.** This is measured bulk
capture followed by a lossless feed-forward handle. It does not establish a
new self-consistent closed-feedback history with the resolved bulk.

## Primary measurements

Values below use lmax=56, dt=1/1024 and the fixed window [-.25,1.75]. Times
are in units of round antipodal transit. Mean and RMS spread use the entire
receiving-port flux in that window, not selected peaks or individual rays.

| b | Aperture a | Carrier w | Capture (%) | Before-source return (%) | Arrival mean | RMS spread |
|---:|---:|---:|---:|---:|---:|---:|
| .9 | .4 | 8 | 25.923 | 25.740 | .94412 | .07682 |
| .9 | .4 | 12 | 14.652 | 14.594 | .94459 | .07530 |
| .9 | .6 | 8 | 33.725 | 33.458 | .91341 | .08694 |
| .9 | .6 | 12 | 10.169 | 10.110 | .91626 | .09097 |
| 1 | .4 | 8 | 50.809 | 49.990 | .97233 | .07042 |
| 1 | .4 | 12 | 42.361 | 42.075 | .97529 | .06893 |
| 1 | .6 | 8 | 50.107 | 48.905 | .95108 | .08521 |
| 1 | .6 | 12 | 18.008 | 17.483 | .96418 | .09357 |
| 1.1 | .4 | 8 | 29.763 | 27.983 | .99382 | .07592 |
| 1.1 | .4 | 12 | 17.018 | 16.522 | .98481 | .07031 |
| 1.1 | .6 | 8 | 34.754 | 33.380 | .96950 | .08754 |
| 1.1 | .6 | 12 | 9.865 | 9.231 | .97290 | .09149 |

![Waveforms and capture fractions](figures/aperture_transfer.png)

A larger normalized aperture is not always a better coherent receiver:
at w=12 it captures less than the smaller aperture. This is an overlap and
coupling result, not a statement that geometric collecting area reduces flux.
The field-energy remainder ranges from 11.55% to 48.51% in these primary
cases. Prompt A reflection plus B output plus the remainder closes the ledger.

For b=1.1,a=.4,w=12, horizontal/fiber displacements of .5 radians change
output histories by .50378/.54725 in source-normalized L2. These four
preselected off-antipodal measurements are a sensitivity scan; they do not
locate a global focal maximum or determine an inverse geometry uniquely.
The xi=0 controls capture 16.761% for b=1.1 and 41.957% for b=1, compared
with conformal values 17.018% and 42.361%. Only those two minimal-coupling
cases were tested; coupling-independent robustness is not established.

## Frozen verdict and diagnostic correction

The protocol and implementation were published at
`cc81c4c4dacb1322810dffd5b434e650cd4af08f` before any validation trajectory.
All 35 cases were run once. Original sources, provenance, raw histories,
result.json and manifest remain unchanged. There were no propagation pilots
in the validation parameter box. Pre-freeze implementation checks are listed
in the [protocol](aperture_transfer_prereg.md).

The sole frozen numerical rejection is reverse_source. Its source is B, but
the frozen diagnostic always interprets outgoing B as captured arrival. It
therefore reports 0.395708 of input as premature arrival: almost entirely
ordinary reflection at the source. The separately measured label-swapped
waveforms already agree to 5.90e-16 in source-normalized L2.

The new replay wrapper exchanges the port labels and flips the antisymmetric
modal coordinates of the archived reverse-source final state. This is the
exact port-exchange symmetry of the Gram representation, not a substitute
forward simulation. It reconstructs every field energy and outgoing sample
from the recorded forces and diagnoses arrival at physical A. No virtual
reversed-handle result is claimed. The corrected receiving-port capture
matches the forward case; its premature flux is below the frozen threshold.
All other diagnostics and all tolerances stay the same.

The original result remains NUMERICALLY_UNRESOLVED. The corrected support
label is explicitly post-measurement. A future independent protocol should
use source-relative receiver indexing from its initial freeze.

## Verification and evidence

| Diagnostic | Largest observed | Frozen limit |
|---|---:|---:|
| Cumulative field-plus-lead energy defect / input | 3.68e-14 | 1e-9 |
| Port reconstruction residual | 2.78e-15 | 1e-8 |
| Aperture omitted squared L2 norm | 3.62e-6 | 1e-4 |
| Premature receiver flux / input, with corrected receiver | <4.10e-19 | 1e-5 |
| Paired mode/time refinement, source-normalized L2 | .007866 | .03 |
| Separate time refinement | .0009902 | .03 |
| Separate mode refinement | 2.11e-9 | .03 |
| Reverse-source waveform comparison | 5.90e-16 | 1e-7 |
| Frequency unitarity / reciprocity defects | <2.04e-15 | 1e-10 |

The extra interval to 2.25 adds 2.85e-10 of input to B capture in its tested
case and reproduces the original waveform interval exactly. It changes the
remaining field energy as more radiation exits A. The positive-delay causal
return has zero before-source energy by construction; separate characteristic
arrival bounds test the resolved bulk's causality rather than crediting that
clock relabelling as independent propagation evidence.

The [35 lossless archives and manifest](../experiments/closure_ledger/runs/20261009_aperture_transfer/)
contain source/output histories, field energies and final modal states.
The pinned manifest SHA-256 is
`26622f2692e5b53daa26c7207cde8548fea08874e4895d431d00b0787fd454d4`.
Replay verifies hashes and source provenance, reconstructs dynamics from
recorded port forces, and reports both verdicts. It never reruns production.

```sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.aperture_transfer_replay
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_aperture_transfer.py
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.aperture_transfer_plot
```

All 21 focused tests pass, including the reverse-source regression, independent
dense stepping, analytic refocusing, frequency unitarity, and rejection of
altered dynamics, source/grid changes, missing/corrupt archives and schedules.
Production used Python 3.12.14, NumPy 2.5.3 and SciPy 1.18.1. Diagnostic
replay with NumPy 2.3.5 and SciPy 1.17.0 reproduces both verdicts and all
seven audited checks; this is a portability check, not a second simulation.

GR support, excised-mouth matching, gravitational recoil and closed-feedback
history remain NOT_ESTABLISHED. The next step is to couple this measured
bulk response back into a fully self-consistent history, retaining its stored
energy and finite aperture response; then replace distributed ports by actual
finite-boundary matching. Neither step is supplied by a Penrose inequality or
by an inverse-boundary theorem alone.
