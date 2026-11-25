# Can an antipodal second-time boundary select a detector record?

A two-time classical model gives an exact conditional result. When motion along the Hopf fiber costs more than horizontal motion, **every globally minimizing history has complementary detector barriers**. No detector label enters the boundary condition. Two identical mechanical scattering probes then give zero doubles and eventual single records with Born frequencies under the invariant history measure. A common finite observation slice still has a positive measure of misses. The round metric, stationary-action principle alone, and a biased preparation each remove a different part of the result.

[Pre-run protocol](antipodal_two_time_protocol.md) · [Solver](../experiments/closure_ledger/antipodal_network/antipodal_two_time.py) · [Numerical archive](../experiments/closure_ledger/runs/20261125_antipodal_two_time/probe.json)

## The two-time variational problem

Let X=(z,w) lie on the unit three-sphere in two complex coordinates (four classical real coordinates). The fiber tangent at X is iX. With the real Euclidean inner product, write V=Vh+ν iX, where ν=〈iX,V〉. For λ>0 use the kinetic action

$$
 S_\lambda[X]=\frac12\int_0^1\left(|\dot X_h|^2+\lambda\nu^2\right)dt,
 \qquad X(0)=x,\quad X(1)=-x .
$$

This is a Berger metric on the mode sphere. Antipodal identification closes the two-time history in the quotient. Only the initial position and its second-time image are prescribed. The free initial velocity and midpoint are found by the boundary problem. Both ports sample the same component intensity at opposite visits, u=|X₁(0)|² and v=|X₁(1/2)|². The boundary condition includes neither the mixer angle nor detector records. The deformed metric is an effective input; the gravity sector would have to supply its sign and magnitude for a physical throat. [Background on the Berger metric](https://www.damtp.cam.ac.uk/user/md327/stephen_thesis.pdf)

The proof works over all absolutely continuous paths of finite action. Pointwise metric comparison, Cauchy–Schwarz on the unit time interval, and the antipodal round distance give

$$
 S_\lambda[X]\ \ge\ \tfrac12\min(1,\lambda)\int_0^1 |\dot X|^2dt
 \ \ge\ \tfrac12\min(1,\lambda)\pi^2 .
$$

Equality requires a constant-speed shortest round semicircle X(t)=x cos(πt)+y sin(πt), with unit real tangent y at x. Along it the vertical velocity is constant, ν=π〈ix,y〉. Equality in the metric comparison now classifies **all** minimizers.

| Vertical cost | Admissible midpoint y | Detector consequence |
|---|---|---|
| λ>1 | Complex-orthogonal to x. For x=(z,w), y=e^{iφ}(−conj(w),conj(z)). | v=1−u for every minimizing history. The arbitrary phase drops out. |
| λ=1 | Any real unit tangent y. | Antipodal closure leaves a continuum of inequivalent ranks. |
| 0<λ<1 | y=±ix. | v=u. At balance half the invariant ensemble eventually writes twice and half misses. |

Each listed family attains the bound, which completes the classification. In particular, at λ=2 the vertical path is still stationary: the fiber velocity is a constant-length Killing field, whose integral curves are geodesics by the Killing equation. Its action π² exceeds the horizontal minimum π²/2. The vertical solution at x=(1,0) misses both ports at balance. If a classical theory admits all stationary solutions, an additional history-selection principle is needed. This conclusion does not require counting the stationary paths or assigning weights to saddles.

As a numerical check, the optimizer varied every interior node of an eight-segment path, with fixed endpoints, three perturbed initial paths and no midpoint constraint. Each segment used the exact action of its round-arc interpolant. All nine solutions reached the analytic lower bound within 4×10⁻¹³. For λ=2, every fitted midpoint had |u+v−1| below 3×10⁻⁷. These finite fits check the calculation; the variational bound controls the continuum.

## How often does the round control fail?

Use the invariant law on ordered real orthonormal pairs (x,y). Conditional on x, the tangent y is uniform on a two-sphere. The component τ=〈ix,y〉 is therefore uniform on [−1,1]. The real orthogonality condition gives |〈x,y〉_complex|=|τ|. Writing xx* and yy* as rank-one projectors, their sum has eigenvalues 1±|τ|. Under joint U(2) rotation the first diagonal entry of the sum, u+v, is uniform between these eigenvalues (the first row of a Haar unitary has uniform squared first modulus). Hence the expected length of #365's failing interval is E|u+v−1|=1/4. With vertical motion it is E|2u−1|=1/2; with horizontal motion it vanishes on every history.

For the balanced control the two unit Bloch vectors for x and y have mutual angle δ=2 arccos |τ| (equivalently their dot product is 2|〈x,y〉|²−1). The threshold at 1/2 is a hemisphere. Two axes at angle δ viewed from a uniform direction have opposite signs with probability δ/π, by the area of their two spherical lunes. Since the two equal-energy probes need opposite rank signs for a single write, balanced misses and doubles each have probability (1−E[δ]/π)/2=1/2−1/π. These Bloch coordinates are only coordinates for two classical complex amplitudes. With 131,072 sampled pairs the mean rank defect was 0.25027, the double fraction 0.18153 and the miss fraction 0.18074 (each analytic probability 0.18169).

## Which measure produces the Born basin?

For λ>1 the complete minimizing histories are ordered complex orthonormal pairs (x,y), parametrized by U(2). Requiring invariance under the transitive U(2) action selects normalized Haar measure on this family. Its endpoint marginal is uniform round volume on S³. In coordinates x=(√u exp(iα),√(1−u) exp(iβ)), that volume is proportional to du dα dβ, so u is uniform. This uses a single setting-independent history measure. The boundary equations fix its support; imposing invariance fixes the weights. For a preparation weighted by 2u relative to invariant volume, all paths still minimize and still have complementary barriers, while the right basin has weight p² instead of p. The preparation of the physical apparatus therefore matters.

## When does the second detector stay dark?

Each visit prepares an identical unit-mass probe at q=−1, moving toward a compact barrier of height b equal to the sampled local intensity. Its potential is b(1−q²)² on [−1,1] and zero outside. A prescribed mixer gives incident energies p=cos²θ and 1−p. The right probe launches at t=0 and the left at t=1/2. A write occurs on crossing q=+1. Beyond that point the probe moves freely outward, giving a persistent transmitted-side record in this auxiliary channel. A subthreshold probe reflects and never writes. Both probes can receive an impulse. Darkness here means absence of a recorded crossing, and the reflected channel can still carry energy.

Newtonian energy conservation yields transmission iff E>b. For complementary heights (u,1−u), the right probe transmits when u<p and the left when u>p. At u=p both are on a separatrix or at rest and neither crosses. This tie has zero invariant volume, and every prepared trial is retained. Consequently the unconditional eventual right/left probabilities are p and 1−p with zero double probability. This obtains the classical switching test from barrier dynamics, using an explicit local intensity–barrier coupling and angle-dependent mixer input.

There is a finite-time cost. For E>b the travel time is

$$
 T(E,b)=\int_{-1}^1\frac{dq}{\sqrt{2[E-b(1-q^2)^2]}} .
$$

It increases continuously from √(2/E) at b=0 and diverges logarithmically as b approaches E from below for E>0. Thus, even after the free-flight time, a positive interval of otherwise transmitting histories misses any finite deadline. Before free-flight time, no history in that port has arrived. Integrating these intervals gives the unconditional readout. The archive includes both equal local elapsed times and simultaneous laboratory slices with the half-unit visit offset.

| Common observation time, balanced input | Right only | Left only | Double | Neither |
|---|---:|---:|---:|---:|
| 2 | 0 | 0 | 0 | 1 |
| 4 | 0.479845 | 0.457592 | 0 | 0.062563 |
| 8 | 0.499933 | 0.499864 | 0 | 0.000203 |
| Eventual | 0.5 | 0.5 | 0 | 0 |

The maximum miss fraction at t=8 on the nine-angle grid is 0.03347, at π/16. The asymmetry between mirrored angles on a common slice comes from the later launch of the left probe. An independent integration of Newton's equation checks the quadrature on both transmission and reflection trajectories.

## Does the correlation survive a load on the midpoint?

After the registered experiment, we added a point interaction −κv to the action at the second visit. This is a simple external load conjugate to the port intensity. It tests the stability of the boundary selection to a specified perturbation; a full reciprocal detector would also vary the apparatus degrees of freedom. The absolute minimizer exists by compactness of the target, coercivity of the kinetic action and continuity of the midpoint potential.

For λ>1 let x=(cos α,sin α) after choosing independent component phases, with 0<α<π/2. Follow the real horizontal semicircle from x to −x, moving at constant but different speeds on the two halves so the midpoint angle is α+π/2+δ. Its bare action is exactly π²/2+2δ², and its midpoint rank is v(δ)=sin²(α+δ). The bare minimizers all have v₀=1−u=sin²α. The change of loaded action relative to those histories is

$$
 2\delta^2-\kappa[\sin^2(\alpha+\delta)-\sin^2\alpha].
$$

The linear term is −κδ sin(2α). For either sign of any nonzero load, some sufficiently small δ gives strictly lower action. Because the bare action is bounded below by π²/2 on **all** histories, every global minimizer of the loaded problem must satisfy κ(v−v₀)>0. This breaks exact complementarity for every interior endpoint, even if the minimizing path leaves our trial family. Positive load opens a miss interval in power-split space, negative load a double interval.

At u=1/2 and κ=±0.1, taking δ=κ/4 certifies |u+v−1|≥0.0124895 for every loaded minimizer. The trial path improves the action by 0.00124896, and the minimum improves it by at least as much. This is a lower bound on a per-history interval of settings. It is not a click probability. A physical backreaction that preserves exclusivity would need a structural cancellation or symmetry that this localized load lacks.

## Physical scope and replay

This is a global BVP for a normalized mode coordinate with one-way detector coupling. It advances [#365](https://github.com/davidmdrpi/geometrodynamics/pull/365) by deriving complementary readouts from an action minimum and testing the records dynamically. The post-primary point load shows how a local coupling can break that geometric selection. The model supplies the fiber stiffness, absolute-minimum selection, invariant preparation, two visit times, mixer and barrier coupling. The outgoing probes retain their records without reinjection. No variation of the action feeds their stress back into the closed trajectory. A joint MTY spacetime and apparatus action would have to preserve the correlation under that backreaction and account for a common readout deadline. [#366's recirculating latches](https://github.com/davidmdrpi/geometrodynamics/pull/366) provide a separate durability test. The many-body and operator-algebra questions in [#232](https://github.com/davidmdrpi/geometrodynamics/pull/232) and [#236](https://github.com/davidmdrpi/geometrodynamics/pull/236) remain open on the classical route.

~~~sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.antipodal_network.antipodal_two_time
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p test_antipodal_two_time.py
~~~
