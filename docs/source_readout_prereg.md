# Preregistration: source readability, field grounding, and cross-round checks

Follow-up to the COMMENT review of PR #283 at `7d70728`.
This file is to be published on GitHub before the new implementations or
numerical probes are written. Predictions below follow from analytic work;
they are not fitted to a numerical result. Priority is the causality gate.

## A. Causality: the full distribution, and the effect of a readout

**A1 — parity is not distributional blindness.** For an antipodally even
source measure, an odd observable has zero mean and a symmetric law. This
does not imply that its law is the same for different settings. Freeze
`a = e_z`, future choices `b_0 = e_x`, `b_1 = e_y`, and the fixed source
readout axis `m = (1,2,0)/sqrt(5)`. Let `F(x) = m.x`. In the positive,
equal-sector, sharp phase-coarea ensemble, both means vanish and both
sign probabilities are 1/2, but the variances are exactly `1/10` and
`4/10`. Each orthogonal-setting circle has outcome-summed density
proportional to `1 + |sin phi| + |cos phi|`, so its two in-plane second
moments are 1/2. This is an explicit counterexample to the broad claim
"all odd observables are blind"; the existing zero-mean/sign checks stand.

Check independent circle integration, the legacy source-density function,
and a whole-sphere finite-temperature integration of #283. At `beta = 64`
the variance gap should exceed 0.15; separately double both integration
coordinates and require variance changes below `1e-4`.

**A2 — finite-resolution readout.** Use a source-local, setting-independent
response kernel `Y = F(x) + Z`, with `Z ~ Normal(0, 0.15^2)`. The event
`|Y| > 0.6` should have a probability difference above 0.2 between the two
frozen settings, in both the sharp and beta=64 ensembles. Evaluate its
Gaussian tail analytically at each source point, rather than sampling a
different kernel for each setting. Without noise the first setting has
zero tail and the second has the circle-density integral above
`|sin phi| > 0.3 sqrt(5)`; check the independent antiderivative. These are
conditional source-record statistics, not by themselves a BAM channel.

**A3 — an explicit conditional classical pointer.** Append a canonical
pointer `(Q,P)` with interaction `h(t) P F(x)`. Its impulse, of area `g`, is

    x' = x,  p' = p - g P dF,  Q' = Q + g F(x),  P' = P.

For `P=0`, the full source phase-space state is unchanged while the record
shifts. This is an exact Hamiltonian extension on a regular rotor chart,
not a one-way equation with the backreaction omitted. For a finite pulse
and free pointer Hamiltonian `P^2/(2M)`, `P=0` stays invariant under any
source Hamiltonian: the old source history and both its boundary values
survive, and `Q` records `int h(t) F(x(t)) dt`. The pointer's final position
is an output, not an independently fixed future boundary value. Energy
change from explicit time dependence is `h'(t) P F`, hence zero on this
preparation. A zero-width momentum preparation is an additional ideal
classical assumption; finite momentum spread produces a bounded recoil
`|delta p| <= |g P| |dF|`.

Check the impulse against independent Hamilton ODE integration (including
nonzero P), its symplectic Jacobian, and a finite-duration pulse coupled
to a nontrivial periodic source orbit. For P=0 require unchanged source
endpoint below `1e-9`; verify that a finite pointer-momentum spread changes
the source. This disproves a universal obstruction from "a measurement
adds a boundary condition" alone. It does **not** establish that BAM has
the required local field observable, coupling, preparation, or readable
record. Those must be audited separately and remain explicit inputs.

**A4 — the actual operational gate.** A setting-independent pointer kernel
is insufficient without a theory of interventions. The required test is
the full early-record law, using the instrument-modified ensemble and
independently selectable future settings, with no selection on a later
record. Any non-readability theorem must constrain the admissible kernels,
the modified source measure, or future-setting controllability. Demonstrate
that a constant-response kernel erases the information and the Gaussian
kernel above retains it. Do not infer a BAM signal from an unmodified
posterior alone; do not infer an obstruction from the absence of a map.

## B. Grounding and assumption budget

Audit actual field variables and couplings in `tangherlini/`, `waves/`,
`transaction/`, `shells/junction.py`, and the conservative Duffing-source
probe. Distinguish a scalar mouth charge or oscillator named q from the
quaternion frame coordinate q. Identify a field-to-triangle map if present;
otherwise state its absence within these implementations, not a theorem
excluding all BAM or nonspherical classical fields.

**B1 — a conditional local parent for the restoring energy.** The existing
massless tube has static Dirichlet-to-Neumann matrix
`(A/L)[[1,-1],[-1,1]]`. Extend it to a matrix-valued transported field,
with endpoint matrices I and R = Ad_G. Minimizing its positive gradient
energy yields `E_min = A ||R-I||_F^2/(2L)`, hence #283's energy with
`K = 8 A/L`. Check both the existing DtN routine and an independent
finite-difference interior solve. This is a realizability construction
requiring extra transported vector/frame channels and endpoint matching;
the present scalar monopole implementation does not supply those channels.
Do not claim the proper throat interval is already the analyzer itinerary,
or that static elimination proves causal or thermal dynamics for it.

**B2 — symmetry is not exact-potential uniqueness.** The Frobenius chordal
distance is bi-invariant, but that invariance permits nonlinear functions
of its square. The quadratic response and its scale are physical choices.
The choice of holonomy versus N is relocated into the restoring coupling;
it is not eliminated. The field inventory will determine whether that
choice is presently derived.

**B3 — round rotor is a stronger sufficient premise than Haar.** The
canonical momentum integral of `H = p_i G^{ij}(x) p_j/2 + V` gives density
`sqrt(det G) d^2x`. Round constant inertia gives sphere area, but does not
follow from a Haar position prior. A rotation-invariant non-Gaussian
momentum law also gives Haar; conversely position-dependent scalar inertia
gives a nonuniform density, and anisotropic tangent inertia with constant
determinant can retain Haar. Check a smooth positive variable-inertia
example and an explicitly integrated non-Gaussian radial momentum law.
No unchanged assumption count will be claimed. Equilibrium is not
equilibration. Record the weaker local-only provenance class of #283's
original freeze separately from this public preregistration.

## C. Cross-round consistency

For multiple non-collinear angles and all four singlet-sign sectors,
compare independent constructions of:

1. triangle quaternion holonomy, N/D phase, and the holonomy-trace scalar;
2. regular phase-gradient/coarea density and #283's transverse stiffness;
3. #281's positive component sum `M_0+M_pi` with #283's M(c), and its
   difference with the oriented integral;
4. normalized positive and oriented sector laws, accounting explicitly for
   the partner sign, factor-two phase convention, and sector priors;
5. source density, parity, and its probability normalization across rounds;
6. #282's normal-window result and #283's normal-energy result.

Use tolerances appropriate to independent numerical integrations; do not
equate finite-window, finite-temperature, and exact limiting values.
Report every discovered disagreement with its convention and scope. Tests
must compare actual values, not a narrative verdict flag. Frozen historical
outputs remain intact; corrections belong in current code and documentation.

## Deliverables and exclusions

Reusable source-readout and grounding controls, independent tests, a probe
with nonzero exit on failed criteria, and a derivation/audit document.
Update the status-of-record `docs/qft_emergence_audit.md`, the relevant
overbroad odd-observable claims, and PR #283. Keep the causality result
conditional unless the actual missing BAM field/intervention maps are
constructed. Composition is lower priority and is not to displace this
work. No claim of a completed quantum reconstruction or program-wide
non-readability theorem is preselected.
