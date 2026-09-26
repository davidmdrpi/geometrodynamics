# Next target: gravitating mouths interacting through the exterior bulk

Date: 2026-09-26 UTC. This is a design and entry-gate specification following
[the #308 surface audit](handle_surface_audit.md), not a claim that a new
two-body evolution or a quantitatively frozen experiment already exists.
Trapped mouths are permitted. Signals need not traverse the neck.

## 1. Construct two objects before claiming reciprocal exchange

Choose and document two distinct localized gravitating regions in the
closed spatial manifold, their enclosing surfaces S1,S2 and their exterior
connection. If the model describes two exterior mouths of one handle,
place the extraction surfaces separately in the two bulk collars and
explicitly describe the connecting neck. Do not use the two identified
faces at s=+/-L as two objects. If the target is two separate black holes,
construct two distinct trapped regions and check whether there is also a
common enclosing trapped surface. A sign change of the normal is not a
second degree of freedom.

Retain the existing four conformal scalars if the goal is to continue the
ESU-supported model. The present work is classical Einstein-scalar gravity,
not vacuum GR with the supporting matter already derived from geometry.
Neither an exact static background nor removal of its growing mode is
allowed during the interaction. Follow proper separation, individual areas,
bulk scale and the support fields. Trapping may be monitored in both
frames, with the Jordan frame identifying the physical probe metric.

Spherical #308 evolution cannot represent general directional recoil or
radiative tensor modes. Extend the spatial/metric ansatz to at least an
appropriate nonspherical sector before claiming that measurement. One
concrete symmetric encounter uses a reference S3 embedded in R4, with
mouth centers at +e1 and -e1 and both initial tangent directions along +e2.
Their reference great-circle paths meet at +e2; rotations in the (e3,e4)
plane preserve the preparation, leaving SO(2) axisymmetry. This defines
what "head-on" means here; the embedding labels are not imposed trajectories
for the gravitating evolution. Axisymmetry can test that sector, not
arbitrary three-component recoil. A second test particle on the old metric
is not this extension.

## 2. Complete and validate the constraints

Specify the free data, global topology, field bundle, time orientation,
extrinsic-curvature sign, gauge, boundary/excision or puncture treatment,
and both local momenta. Solve Hamiltonian and momentum constraints with
the full source included. An imposed symmetry may prepare equal-and-opposite
initial charges, but the code must still evolve both regions and evaluate
both observables separately. Allow a small asymmetric control so the
measurement cannot pass solely by constructing one result as minus the other.

For a puncture or inner-boundary construction, establish regularity and
constraint-compatible boundary conditions rather than declaring the old
Bowen-York momentum seed a completed solution. Determine whether the two
surfaces bound an exterior domain or are nonseparating on the chosen
manifold. No global charge-selection rule is imported from a different
topology.

## 3. Specify momentum and its balance before measuring it

Choose a finite-region charge definition and fix its observers, reference
subtraction (if any), normalization and comparison vector/frame. A surface
Hamiltonian or timelike-worldtube stress construction is a possible route;
[Brown and York (1992)](https://arxiv.org/abs/gr-qc/9209012) gives the
boundary-stress framework. Its full balance law must be derived for the
actual chosen boundary motion and scalar system before implementation.
A Brown-York energy formula by itself is not a recoil observable.

Name the slicing used by any conformal-Killing constraint-charge check.
With pi^ij=K^ij-K gamma^ij, a conformal-Killing comparison vector obeying
`D_(i X_j)=lambda gamma_ij` contributes `-2K lambda` to the divergence
identity. Use maximal slicing K=0 to remove that term, or retain it.
If the evolved vector is not conformal Killing, retain the full
`pi^ij D_(i X_j)` term. In particular #308's geodesic slicing does not
preserve K=0. Its radial `partial_s` charge is not the translational
observable needed for recoil.

Measure P1(t) and P2(t) separately on surfaces that remain geometrically
attached to the two regions. Track surface motion, changes in the comparison
frame and extraction-radius dependence. Horizon area growth is useful
alongside this measurement but is not linear momentum by itself. A closed
universe supplies no asymptotically flat ADM momentum frame.

Integrate the exterior contribution independently, including matter flux,
geometric/boundary stress, observer work and transport terms required by
the chosen formulation. Partition the geometry before writing a ledger:

- **Two separately enclosed objects:** their independently measured charges
  and the exterior contribution may form a balance of the schematic form
  `Delta P1 + Transport(Delta P2) + exterior contribution = residual`.
  The global constraints still restrict the allowed initial data; separate
  measurements do not mean freely specifiable independent momenta.
- **Two collars of one handle:** the two nonseparating extraction spheres
  together bound both a neck region N and a complementary bulk region B.
  Include the neck contribution as well:
  `Delta P1 + Transport(Delta P2) + exterior contribution + neck contribution = residual`.
  If the interior is excised, account for the inner-boundary flux/work in
  the retained-domain balance instead of evolving its content. Do not
  count both descriptions of the same contribution.

These are accounting requirements, not derived or implemented equations;
the chosen surface charges and volume partition must avoid double counting.
For the spatial constraint check on a retained one-handle slice, choose
normals outward from N and consistently restrict a global comparison X.
Then

    Q1_N + Q2_N = integral_N [j_j X^j + pi^ij D_(i X_j)] dV,
    Q1_B + Q2_B = integral_B [j_j X^j + pi^ij D_(i X_j)] dV,
    Q1_B = -Q1_N, Q2_B = -Q2_N.

This illustrates why equal-and-opposite boundary bookkeeping on a slice
cannot determine whether a disturbance propagated through the neck or
through the bulk. Do not force Delta P2=-Delta P1 at every time: retained
neck and bulk degrees of freedom can contribute while a disturbance
propagates. The spatial momentum-constraint integral
in the audit is an initial-data/constraint check, not a substitute for this
time-dependent ledger. General relativity has no unique local tensorial
gravitational momentum density; the boundary formulation must state what
its exterior term measures.

## 4. Separate propagation, backreaction and the controls

Prepare a localized perturbation near one region and solve the constraints
again. Compare with the unperturbed coupled evolution. Constraint solving
is elliptic, so the perturbed initial metric can already differ globally:
do not mistake that initial difference for superluminal propagation.
Track characteristic propagation in the evolved exterior geometry and
the subsequent change of the receiver's independently measured charge.
For a perturbation with genuinely localized causal support, a response
attributed to the exterior must not precede the earliest exterior null
arrival. Establish matching initial data on the relevant receiver domain
of dependence, or explicitly separate the global initial-data change;
otherwise a difference before that arrival is not a clean propagation test.
Use the physical metric and an explicit clock comparison, not coordinate
distance divided by a background speed.

For a retained handle, track the competing neck causal route too. The
measured local trapped spheres alone do not establish a global causal
barrier or prove that every response must use the bulk; distinctions
between trapped regions and their causal boundaries already arise in
[spherical examples](https://arxiv.org/abs/1009.0225). Establish the absence
of a relevant neck causal path over the measured interval, or quantify
and separate its contribution. With excision, verify that all relevant
characteristic fields at the inner boundary flow out of the retained
computational domain; merely labeling it "trapped" is insufficient.
Only after that check can an appropriately timed receiver response be
attributed to the exterior route. Traversability remains unnecessary.

Evolve for the measured bulk travel time, not merely the old t<=.01 neck
window. Retain the ESU support's instability and report any gauge, curvature,
positivity or constraint stop as a failure/inconclusive endpoint.

Use convergent spatial/time refinement, more than one extraction radius,
zero-perturbation and sign-reversed preparations, and an asymmetric
preparation. Report both raw charges and their matched-reference
differences with the same gauge, surface and frame prescriptions. The
existing audit gives a concrete warning: the final eta=.3 minus eta=0
cut-flux difference is +5.80584e-7, while the common collapse contribution
is about -0.003954. A reference difference can reveal preparation response,
but is not automatically a gauge-invariant physical momentum. Apply the
same subtraction to every term of the chosen balance. A fixed-background/test-particle calculation is a labelled
negative control for reciprocal backreaction. A pulse must be made from
retained dynamical fields/metric degrees of freedom, not a prescribed
mouth force. If gravitational radiation is the claimed carrier, resolve
its nonspherical degrees of freedom and distinguish its transfer from
scalar-matter flux.

Before any production run, freeze the actual initial-data family, solver,
gauge, duration, resolutions, extraction surfaces, numerical error budgets,
controls and failure criteria. Those choices depend on successfully
constructing the two-region geometry, and are not silently filled in by
reusing #308's reflection-symmetric one-neck data.

## 5. Keep the milestones separate

| Milestone | Evidence required | Current status |
|---|---|---|
| Future-trapped neck in #308 | Both future expansions negative, frame and convergence checks | Established at saved t=.01 |
| Two distinct tracked gravitating regions | Geometric identity, valid constraints and separate extraction surfaces | Not tested |
| Exterior-bulk reciprocal momentum response | Independent P1,P2, resolved carrier and convergent full balance | Not tested |
| Finite-duration discrete event/action selection | Mechanism and robust action selection without imposed detector or kick | Not tested |
| Emergent quantum statistics | Derived preparation/measurement probabilities and composition | Not established |

Smooth classical exchange can be a useful prerequisite without being a
quantum. Do not require an instantaneous momentum discontinuity as the
only possible later selection mechanism; equally, integrating a continuous
force over a chosen interval does not establish an action quantum.
