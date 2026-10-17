# PR341 review follow-up

Source: [independent review of fd12a3f](https://github.com/davidmdrpi/geometrodynamics/pull/341#issuecomment-6239379678).
Reviewer scripts and figure are preserved at
[d35e81d](https://github.com/davidmdrpi/geometrodynamics/tree/d35e81d47b20dbf257c2e8574c20a64d2ddd6d11/docs/reviews/pr341).

## Accepted verification and corrections

The reviewer reports all 48 moving trajectories and four direct references
reproduce with identical numerical leaves, 17 passing tests, and byte-identical
gates and numerical report. Their independently written fixed-axis Strang
solver, extrapolated in timestep, agrees with the adaptive reference to about
3e-10 relatively. Their finite-amplitude moving runs approach that coefficient
with the predicted quadratic correction; those finite-amplitude values are
not themselves within 1e-9 of the limiting coefficient.

The worker completion order changed JSON key order, so equal values did not
produce equal archive bytes. Both checkpoint and report writers now use one
canonical sorted-key serializer. A nested-order permutation test covers the
reported failure. The existing archive was reserialized and its manifest
updated; no trajectories were rerun or numerical leaves changed. The protocol
and original results remain available in prior commits. The report no longer
writes /tmp/mouth_exchange.png. SVG bytes may still depend on Matplotlib version.

The remote protocol commit followed the stated local freeze. The unpublished
local 0206381 is not independently verifiable through GitHub history; the
existing replication/non-blind-preregistration qualification remains in place.

## What the review adds physically

The reviewer resolves pair 1's mixed y torque into two terms. At the orthogonal
preparation, e^x D_y U^z vanishes numerically, while e^z D_y U^x carries the
exchange. A wavefront from the other packet encounters unequal shape mismatches
at the two mouths: the source-side mouth has been excited, whereas the target
mouth has not yet experienced its own focus. Opposite gradients therefore need
not cancel. The strongest contribution is associated with the mid-transit
crossing near pi/2, consistent with our archived event-window accounting.

These term-resolved results and the start-offset scan are **reviewer results,
not independently rerun here**. Their reported delta=1.2 coefficients are
approximately 1.020e-3 and 1.269e-3, about 2.1 times the delta=0 coefficients.
Reciprocity persists; its magnitude is preparation-dependent. The scan does
not establish an asymptotic delta-independent value: the delta=.6 and 1.2
results still differ by several percent.

Omega=5 is the bare shape frequency, not a demonstrated frequency of the fully
coupled ringing. The mismatch dynamics include both coupling stiffness and
field response. Likewise, "not localized at the antipodal focus" describes
what was measured without suggesting that a quantum event was computed.

## Next experiment, not yet run or preregistered

Use theta in {60,75,90,105,120} degrees and explicitly defined start offsets,
including delta=0 and a full-source-focus preparation such as delta=1.2.
Rotate the corresponding packet center, mouth pair and directional preparation
consistently; preserve common joint normalization and document the resulting
initial energy. State precisely which pair and mouth is being compared.

There are **two geometric arrivals**, not a single universal peak at theta.
For source -n and pair mouths +/-z, with n dot z=cos(theta), distances are
theta to -z and pi-theta to +z. These coincide only at 90 degrees. Measure time
relative to the source focus (subtract delta). Compare per-mouth signed torque
and absolute torque with both arrival windows, alongside source and target
focus windows. The source-side ringing suggests a dominant contribution near
theta, but this must be tested, not imposed by selecting a peak after the run.
Shape phase and finite packet width can shift or cancel the net pair response.

Before running, freeze window widths, overlap accounting, timing tolerance,
resolution criteria, and a failure rule. Failure to track the predicted
mouth-specific arrivals, or persistence only under a special start convention,
would defeat the simple crossing explanation in that form. A focus-localized
interpretation must separately survive a quantitative inside/outside-focus
absolute-torque test. This follow-up does not promote either hypothesis to a
passed result or silently extend the frozen PR341 protocol.
