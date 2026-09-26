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
appropriate nonspherical sector before claiming that measurement. Axisymmetry
may suffice for a chosen head-on preparation; it does not test arbitrary
three-component recoil. A second test particle on the old metric is not
this extension.

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

Measure P1(t) and P2(t) separately on surfaces that remain geometrically
attached to the two regions. Track surface motion, changes in the comparison
frame and extraction-radius dependence. Horizon area growth is useful
alongside this measurement but is not linear momentum by itself. A closed
universe supplies no asymptotically flat ADM momentum frame.

Integrate the exterior contribution independently, including matter flux,
geometric/boundary stress, observer work and transport terms required by
the chosen formulation. The balance has the schematic structure

    Delta P1 + Transport(Delta P2) + exterior contribution = residual.

This is a design requirement, not yet a derived or implemented equation.
Do not force Delta P2=-Delta P1 at every time: the bulk can carry momentum
while a disturbance propagates. The spatial momentum-constraint integral
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
Evolve for the measured bulk travel time, not merely the old t<=.01 neck
window. Retain the ESU support's instability and report any gauge, curvature,
positivity or constraint stop as a failure/inconclusive endpoint.

Use convergent spatial/time refinement, more than one extraction radius,
zero-perturbation and sign-reversed preparations, and an asymmetric
preparation. A fixed-background/test-particle calculation is a labelled
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
