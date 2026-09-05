# Pre-registration: classical frame stiffness and the closure measure

Frozen before `bulk/closure_equilibrium.py`, its tests, or its probe exist.
Base: `92a915bfaaabd02564accd7467b81cedb1ee8c16` (after PR #281).

## Question and scope

Does a specified classical energy penalizing frame mismatch reproduce the
phase-window measure of `closure_measurement.py`? Which parts of that
measure are geometric, which belong to the preparation or apparatus, and
can this mechanism select the holonomy-weighted branch?

This is a **conditional model**, not a claimed derivation of a detector
coupling from the existing BAM action. We keep the geodesic itinerary of
`closure_current.py`; its outcome signs are boundary sectors, not outputs
of a solved local detector. The frame direction is treated as physical.
Canonical equilibrium, a round constant-inertia rotor, equal sector
coefficients for the baseline, and the new restoring energy are explicit
assumptions. No equilibrium or local implementability theorem is claimed.

## Model fixed before calculation

Use the singlet-sign triangle `x -> u -> w -> x`, where `u = s_A a` and
`w = -s_B b`. Write

```
N = x.(u x w),   D = 1 + x.u + u.w + w.x,
theta = atan2(N, D),   G = cos(theta) + x sin(theta).
```

The full loop has an additional `-1` from `J^2`; the classical frame
`R = Ad_G` cannot distinguish it. The two punctures `x = -u, -w`, where
the geodesic holonomy is undefined, are excluded. All statements are for
fixed non-collinear analyzer directions (`0 < gamma < pi`); no uniform
collinear limit is claimed.

The proposed apparatus compares the returned oriented orthonormal frame
with its initial frame using an isotropic quadratic mismatch:

```
V_frame = (K/16) ||R - I||_F^2 = (K/2) sin(theta)^2,
H_s = |p|_round^2/(2I) + V_frame,
b = K/(k_B T),
Z_s(b) = (1/4pi) int_S2 exp[-b sin(theta_s)^2/2] dOmega.
```

Integrating the momenta of the round rotor gives the area measure. The
common momentum factor cancels between sectors only for identical rotor
kinetics. Outcome weights in this *history-ensemble model* are
`pi_s Z_s / sum(pi_t Z_t)`. Reading them as actual event frequencies remains
an assumption to be derived from a BAM apparatus later.

## Analytic predictions and falsifiers

**P1 — frame identity.** `||Ad_G-I||_F^2 = 8 sin(theta)^2`. Verify using
the existing quaternion geodesic rotations and `mouth_spin_frame`, rather
than building both sides from the same trigonometric expression. The cost
is unchanged by `G -> -G`, by simultaneous spatial rotations, and by the
initial frame. Falsifier: residual above `1e-11` on regular random loops.

**P2 — low-temperature measure.** Gaussian integration normal to the
closure circle gives, on every regular arc,

```
sqrt(b/(2pi)) Z_s(b) -> (1/4pi) int_Gamma |D_s|/|u x w| d sigma.
```

The two punctures carry no limiting mass: the excluded arclength integral
is `O(epsilon^2)`. For `c = u.w`, `delta = acos(c)`, the unnormalized
circle integral is

```
M(c) = 4 + 2(1+c)(pi-delta)/sqrt(1-c^2).
```

Thus equal sector coefficients recover the **positive** phase coarea
measure, with the singlet partner sign: `E(gamma) = -E_round5(gamma)` and
standard-angle `CHSH = 2.142283163...`, not `2 sqrt(2)`.

Compute `Z` on the whole sphere at `b = 16, 64, 256, 1024, 4096`, never
restricting the quadrature to the closure circle. At `gamma = pi/4, 1,
3pi/4`, require the final normalized joint law to agree with the limiting
one to `2e-3` and the scaled sector integrals to `1%`; refine both
quadrature directions independently to establish numerical error below
`1e-4` in the joint law. The finite-b trend is a result, not a frozen
monotonicity or convergence-order assertion. If thresholds are missed,
report the miss rather than changing them without a correction note.

**P3 — a physical countermodel on the same zero set.** Replace the
potential by `V_N = K N^2/2`. Now the same closure circle has limiting
density `1/|u x w|`, uniform in arclength. In fact `N_s^2` is identical in
all four sectors, so their partition functions are equal for every b and
`E = 0` for equal priors. Require equality to `1e-12`; independently check
against its one-dimensional Gaussian integral over the sphere's normal
coordinate. This is a different restoring energy, not a contradiction of
P2 or a change of coordinates in the same physical apparatus.

**P4 — reparametrization covariance.** For a smooth simple-zero residual
`F` and positive stiffness `a(x)`, `V = K a F^2/2` gives limiting density

```
d mu_Gamma proportional to d sigma / (sqrt(a) |grad F|).
```

Under `F -> g(x) F`, `a -> a/g(x)^2` with smooth positive g, the energy
and finite-b measure are exactly unchanged. Holding a fixed instead
changes the apparatus. Test both: agreement to `1e-12` for the covariant
case and a nonzero effect for fixed stiffness. Use the analytic control
`F=N`, `g=1+epsilon x.(u+w)`, `epsilon=1/4`; the limiting mass is
`2pi/[|u x w| sqrt(1-2 epsilon^2 (1+c))]`. This distinguishes a residual
coordinate from a physical response law.

**P5 — universality with limited premises.** The nonnegative deformation
`V = (K/2)(sin^2 theta + lambda sin^4 theta)`, `lambda >= 0`, changes
finite-b corrections but not the leading transverse Hessian or limiting
measure. Test `lambda=0,2` at `b=4096`, same P2 tolerances. This is a
restricted family of classical frame potentials, not all classical BAM
actions. No oriented cancellation can arise by assigning negative
probabilities to the thermal histories: their weights remain positive.
This positivity statement does **not** forbid a different positive model
from reproducing the integrated quantum law, or forbid coherent classical
field responses.

**P6 — unresolved inputs survive.** Vary the like/unlike sector prior
ratio through `0.5,1,2`; verify it changes E while the paired-sector
symmetry preserves half marginals. Report this rather than deriving the
equal prior from isotropy. Neither composition, source-readout causality,
nor a local detector event mechanism is supplied here.

## Independence and deliverables

Implement a small reusable module, tests, a probe with an explicit nonzero
exit on failed criteria, and a derivation document. Verify whole-sphere
partition integrals against the limiting one-dimensional integral and
against independent spherical Monte Carlo at a moderate b, with sampling
uncertainty reported. Tests should exercise the numerical and analytic
claims, not assert that a narrative verdict string is true.

Correct the existing statement “Haar conditioned on N=0 gives |D|” to
identify the **phase-window** prescription. Preserve its implemented
behavior and its historical provenance. Update the README with a concise
conditional result and a link. Do not relabel the program-wide probability
fork as closed.

Potential successful result: **the stated classical frame-equilibrium
model derives the positive phase-coarea law, and physical stiffness—not
the closure locus alone—determines the limiting measure**. A failed
oracle, unresolved quadrature, or a mismatch with this prediction must be
reported explicitly. No numerical prediction is selected because it is
closer to a quantum target.

## Integration and publication note (after implementation)

The text above was frozen in local commit `76ed50e` before the module,
tests, or probe were written. Shell GitHub push authentication was
unavailable, so this preregistration and the subsequent implementation
were uploaded afterwards through the authorized GitHub connection. This
is a local preregistration, not a claim of public timestamping before the
calculation; the original local hash is retained as provenance.

PR #282 landed while this work was being developed. Its independent
conditioning-variable correction is incorporated from main commit
`e96d48abf57fc97f135bb323fe0da116793d1077`. It supplies the correction
requested above; this work adds the conditional equilibrium mechanism.
No prediction, acceptance threshold, or model above changed on integration.
