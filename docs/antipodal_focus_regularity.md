# Antipodal regularity: foundation, equations and a computable boundary test

11 October 2026. Foundational follow-up to #327. This revises the proposed
theory; it is not derived from unmodified vacuum GR. **No nonlinear antipodal
gravitational-wave event has been computed here.**

BAM assumes Einstein dynamics outside focal regions and an additional
regularity/matching postulate at the focus. An accepted event must continue
through a regular bulk region rather than terminate at a physical singularity.
The claim that the encounter is a quantum event remains a hypothesis.
The extra postulate must ultimately select a dynamical law. This document
shows why regularity alone does not do so.

## 1. A caustic is not a curvature singularity

For the exact ultrastatic spacetime

\[
 g=-dt^2+R^2d\Omega_n^2,\qquad J(\lambda)=R\sin(\lambda/R),
\]

the transverse Jacobi factor vanishes at \(\lambda=\pi R\), while
\(R_{abcd}R^{abcd}=2n(n-1)/R^4\) is constant. The spacetime is geodesically
complete. It is not vacuum: at \(\Lambda=0\),
\(\kappa\rho=n(n-1)/(2R^2)\) and
\(\kappa p=-(n-1)(n-2)/(2R^2)\).
Thus a ray caustic does not itself force a physical singularity in Einstein
theory with explicit support. Our control evaluates this distinction;
it does not evolve a wave packet. Caustic propagation effects also occur in
classical linear wave theory [4].

The revised postulate forbids physical singularities in accepted solutions.
It does not forbid a ray caustic, a horizon, or a future-trapped surface.
Both future null expansions must be measured to classify trapping.
A neck minimum is not enough; an ER bridge need not be traversable.
Bounded scalar invariants alone do not guarantee bounded parallel-propagated
curvature or geodesic completeness. A coordinate horizon is not a failure;
a solver crash is not a bounce.

## 2. AP-R0: an admissibility system, not yet a complete event law

Let \(\mathcal U_f\) be a finite region containing the encounter, bounded by
matching worldtubes \(\mathcal B_\pm\). Its location, scale and coupling must
be specified as part of the model, before a divergence occurs. In the intended
bulk \(D=5\); the observed closed \(S^3\) and the four-dimensional spatial
bulk are distinct. A reduced \(S^2\) visualization cannot supply this model.

Outside \(\mathcal U_f\), require

\[
 G_{ab}+\Lambda g_{ab}=\kappa_D T^{\rm support}_{ab}.
\]

Inside, specify either modified geometric equations or an effective sector:

\[
 G_{ab}+\Lambda g_{ab}=\kappa_D(T^{\rm support}_{ab}+T^f_{ab}),\qquad
 \nabla^a(T^{\rm support}_{ab}+T^f_{ab})=0.
\]

\(T^f\) needs a constitutive/evolution law. Defining it after choosing a
metric is an inverse construction, not a predictive matter theory.
More generally, for any sufficiently smooth candidate metric one can set
\(\kappa_D T^{\rm total}=G+\Lambda g\); the Einstein constraints then hold
by definition and conservation follows from the contracted Bianchi identity.
This does **not** establish boundary compatibility, a realizable matter action,
causal completeness or well-posed perturbation dynamics. It explains why an
unspecified effective stress cannot by itself select the focus response.
Gravitational radiation is evolved in \(g\); adding a separate exact local
gravitational-wave stress to the full Einstein equation would double-count it.
An averaged effective stress must declare its approximation and validity.

An accepted continuation must satisfy:

1. A nondegenerate Lorentz metric with a specified differentiability class:
   at least \(C^2\) for the target evolution, or explicitly formulated weak
   junctions with all distributional stress accounted for.
2. Finite curvature in regular parallel-propagated frames along causal probes
   and extension of those probes through the event. A numerical curvature
   bound alone is not a completeness proof.
3. In a handle description, a positive minimum areal radius \(r\ge b_*>0\),
   with lapse, shift, scale and possible evolution specified.
4. Constraint preservation, stress/flux balance and an explicit outgoing
   field/internal-state rule.
5. Existence, stability and uniqueness in the chosen class, or an explicit
   rule for multiple solutions. Regularity is not that rule.

These are **AP-R0 admissibility conditions**. A candidate must complete an
action or equivalent closed equations, for example

\[
 S=\frac{1}{2\kappa_D}\int_{\mathcal M\setminus\mathcal U_f}
       (R-2\Lambda)\,dV_g+S_{\rm GHY}+S_{\rm support}
       +S_f[g,\xi,\Psi;\lambda_i].
\]

\(\xi\) denotes proposed focus degrees of freedom, \(\Psi\) transported fields,
and \(\lambda_i\) every added parameter. An unspecified \(S_f\) is a **missing
equation**, not an implemented regularity law. Any candidate needs a constraint,
hyperbolicity and instability analysis. Limiting-curvature and black-bounce
models illustrate different additions to the theory, not interchangeable
prescriptions or guarantees that all regular cores are stable [5,6].

## 3. Inner/outer matching and Alice data

For timelike boundaries, specify \(\Phi:\mathcal B_-\to\mathcal B_+\), including
the angular map and any clock shift, and require \(h_-=\Phi^*h_+\).
With one consistently oriented transverse normal, thin-layer matching obeys

\[
 [K_{ij}]-h_{ij}[K]=-\kappa_D S_{ij},\qquad
 D_iS^i{}_j=-[T_{ab}n^a e^b{}_j].
\]

Normal reversal must be included when converting two outward-normal
conventions. Shell-free matching sets \(S_{ij}=0\). Null boundaries require
null-junction equations; the timelike formula cannot simply be carried
through a null boundary. A thick core solves its interior equations instead.

A possible lossless scalar identification, after pullback to one boundary,
uses outward canonical normal fluxes:

\[
 \Psi_+=U\Psi_-,\qquad \Pi_++U\Pi_-=0,\qquad U^TU=I.
\]

This is a candidate transmission condition, not a derived quantum operator.
Core storage, absorption and internal-state changes require their own terms.
An Alice lift is additional bundle data: complex conjugation, represented
as an orthogonal reflection on two real field components, reverses a specified
\(U(1)\) charge. Metric non-orientability alone does not force that gauge
action, spin statistics or a quantum amplitude.

On the local symmetric neck below, \((s,\Omega)\mapsto(-s,-\Omega)\) is a free,
time-preserving isometric involution: at \(s=0\), no point of \(S^3\) equals
its antipode. The angular antipodal map on \(S^3\) has degree \(+1\); radial
reflection supplies the spatial orientation reversal. This local quotient
does not establish a global extension over the specified exterior or select
the field-bundle lift.

An existing bulk handle whose observable access changes differs from creation
of new spatial topology. Smooth globally hyperbolic spacetimes have product
topology \(\mathbb R\times\Sigma\) [7]. A topology-creation model must specify
which hypothesis changes. If an MTY time shift creates closed causal histories,
the global consistency/uniqueness problem replaces a simple Cauchy evolution
claim. This PR's controls have zero mouth time shift.

## 4. Computable counterexample on the existing finite-mouth geometry

The repository's finite-mouth construction already assumes a round \(S^4_R\)
spatial bulk, mouth angle \(a\), a time-symmetric slice and a scalar-flat neck.
Those assumptions give

\[
 b=R\sin^2a,\quad S=R\sin a\cos a,\quad r(s)=\sqrt{s^2+b^2},\quad |s|\le S.
\]

Keep these data and introduce a regular lapse family:

\[
 ds_5^2=-N_\epsilon(s)^2dt^2+ds^2+r(s)^2d\Omega_3^2,\qquad
 N_\epsilon(s)=\exp[\epsilon(1-s^2/S^2)^3],\quad 0\le\epsilon\le0.1.
\]

Extend \(N=1\) outside the neck. At either mouth \(N=1,N'=N''=0\).
All members have the same spatial metric, \(K_{ij}=0\), neck/mouth radii,
boundary clocks and extrinsic matching data. The lapse obeys the uniform bound 1 <= N <= exp(0.1).
The old spatial join to the round exterior is Darmois \(C^1\) with a bounded
curvature jump, not globally \(C^2\); the new lapse adds no shell or jump.
This is a **piecewise-smooth junction control**. A fully \(C^2\) global candidate
needs a specified smoothing collar, its stress and recomputed transport.
The review supplement in §4b supplies that control without changing this
original producer or its archive.

Set \(A=r''/r,\ B=(1-r'^2)/r^2,\ C=N''/N,\ D=N'r'/(Nr)\).
The required orthonormal effective stress for \(\Lambda=0\) is

\[
 \kappa_5\rho=3(B-A)=0,\quad
 \kappa_5p_s=3(D-B),\quad
 \kappa_5p_\Omega=C+2D+2A-B,
\]

and

\[
 \mathcal K=4C^2+12D^2+12A^2+12B^2,\qquad
 \kappa_5(\rho+p_s)=3(D-A).
\]

All terms are bounded for \(b>0\). This static ansatz gives a positive
sum of squared independent curvature components, unlike a generic Lorentzian
invariant. The pressures differ across the family even though the scalar
constraint is identical. These are **inverse-engineered stresses**, not a
constitutive law for perturbations: there is no contradiction with uniqueness
for a fully specified Einstein-matter Cauchy problem.
The controls check regularity, initial constraints and matching. They do not
certify the dynamical well-posedness part of AP-R0 and are not yet admitted
as complete event models. Lapse is a gauge choice within a fixed development;
here the invariant stresses and boundary-clock transit change, so these are
different proposed developments with different unspecified support dynamics.

At the waist, \(\kappa_5(\rho+p_s)=-3/b^2<0\) for every member. The required
radial null-energy-condition violation remains explicit. An effective sector
must supply it, or modified geometric equations must replace that stress
interpretation. For \(k_\pm=N^{-1}\partial_t\pm\partial_s\),
\(\theta_\pm=\pm3r'/r\): both vanish at the waist and have opposite signs
elsewhere. These controls are not strictly future-trapped; that does not
exclude trapped mouths in the theory.

The radial null transit time in the common mouth clock is
\(\Delta t_\epsilon=\int_{-S}^S ds/N_\epsilon\). It changes with \(\epsilon\).
This is not a free clock rescaling: the boundary clocks stay fixed to \(N=1\),
and the effective stresses differ.

For a finite-frequency scalar control, put
\(\Psi=e^{-i\omega t}\phi(s)Y_\ell(\Omega)\) and \(p=Nr^3\phi'/b^2\).
The transfer between the actual mouths satisfies

\[
 \frac{d}{ds}\binom{\phi}{p}=
 \begin{pmatrix}
 0&b^2/(Nr^3)\\
 [Nr\,\ell(\ell+2)-\omega^2r^3/N]/b^2&0
 \end{pmatrix}\binom{\phi}{p}.
\]

The real generator has zero trace, so it preserves the directed Wronskian
flux and \(\det M=1\). Conservation does not select one transfer map.
This scalar diagnostic is **not a spin-2 wave calculation, a backreaction
solve, or an event/creation amplitude**.

For R=1, a=0.5 and epsilon=0 versus 0.1:

| Observable | epsilon=0 | epsilon=0.1 |
| --- | --- | --- |
| Neck radius | 0.229848847066 | 0.229848847066 |
| Null transit time | 0.841470984808 | 0.804399472467 |
| Largest sampled Kretschmann scalar (1001 neck samples) | 8598.88069 | 8644.83496 |
| Radial NEC contraction at the neck, kappa times stress | -56.7853877 | -56.7853877 |

The mouth-clock transit changes by 4.40556%.
All six fixed mode/frequency comparisons (ell=0,1; omega*S=0.5,2,5) give
different transfer matrices: the smallest maximum-entry difference is
0.0906232099, using the same dimensionless canonical state normalization.
The largest refinement difference is
3.04e-10; the largest determinant defect is
1.14e-11.
The full matrices and source hashes are saved in
[the control archive](../experiments/closure_ledger/runs/20261011_focus_regularity/result.json).
These numbers demonstrate nonuniqueness under the stated admissibility
conditions; they do not fit a BAM parameter or favor one epsilon.


## 4a. Review control: the vacuum development is singular

The original spatial neck with \(K_{ij}=0\) is the bifurcation slice of
the five-dimensional Schwarzschild–Tangherlini vacuum, locally

\[
 g=-Fdt^2+F^{-1}dr^2+r^2d\Omega_3^2,\qquad F=1-b^2/r^2.
\]

Here \(dr^2/F=ds^2\) on \(r^2=s^2+b^2\). The static coordinates break down
at the horizon; the vacuum development is understood in regular Kruskal
coordinates. This is **vacuum Einstein evolution of the local initial data**,
not vacuum matching to the same ultrastatic exterior. The latter is not
shell-free, as `finite_mouth.py` already records.

The timelike geodesic from the bifurcation waist has zero Killing energy:
\(\dot r^2=b^2/r^2-1\), hence

\[
 r(\tau)^2=b^2-\tau^2,\qquad
 \tau_{\rm end}=\int_0^b\frac{r\,dr}{\sqrt{b^2-r^2}}=b,
 \qquad \mathcal K=72b^4/r^8.
\]

Its finite-proper-time endpoint has divergent curvature. Initially the
waist has both null expansions zero; future interior spheres become
future-trapped, since both future radial null directions decrease \(r\).
Thus initial marginality does not mean the vacuum waist remains static.

To check whether exterior changes can affect this conclusion, fix the
Kruskal normalization

\[
 X^2-T^2=e^{2r/b}\frac{r-b}{r+b},\qquad
 X_m=e^{r_m/b}\sqrt{\frac{r_m-b}{r_m+b}},\quad r_m=R\sin a.
\]

On the initial slice \(T=0\), the neck occupies \(-X_m<X<X_m\).
The waist approaches the singular boundary at \((T,X)=(1,0)\).
Radial characteristics have slopes \(\pm1\), so its entire regular
approach is shielded from exterior initial data if **\(X_m>1\)**.
Strictly, the singular boundary is not a spacetime point in a domain of
dependence; the regular points approaching it are in that domain.
At \(R=1,a=0.5\), \(b=0.229848847066\) and \(X_m\simeq4.78>1\).
The supplement also checks \(a=1.3\), for which \(X_m<1\); that case
does not support the same exterior-independence conclusion.

Consequently, retaining these local data and vacuum dynamics until the
endpoint fails the BAM regularity postulate in the shielded case. Avoiding
that endpoint requires a change in the interior evolution or initial data
before it occurs. This calculation does **not** prove that every regular
alternative must have a static open neck or the NEC value below, nor that
every ray caustic collapses. A complete focus action is the next missing
input; its well-posedness, stability and quantum predictions would still
need separate demonstration.

## 4b. Review control: smooth collars preserve the ambiguity

Smoothing only \(r\) would not make the original spacetime \(C^\infty\):
the original lapse extended by one is \(C^2\), with a third-derivative jump.
The supplement therefore replaces **both** profiles and records new results
separately. Set \(\delta=cS\), \(L=S-\delta\), with disclosed control inputs
\(c=0.1,0.2\). On each collar \(L<|s|<S\), use

\[
 z=(|s|-L)/\delta,\quad
 w(z)=\frac{e^{-1/z}}{e^{-1/z}+e^{-1/(1-z)}},\quad
 r_c=(1-w)\sqrt{s^2+b^2}+wR\sin\bigl(a+(|s|-S)/R\bigr).
\]

Extend \(w=0\) for \(z\le0\), \(w=1\) for \(z\ge1\). Its derivatives
of all positive orders vanish at both endpoints. Keep the round continuation
outside the original mouths, and use the compactly supported lapse

\[
 \log N_\epsilon(s)=
 \begin{cases}
 \epsilon\exp\bigl(1-[1-(s/L)^2]^{-1}\bigr),& |s|<L,\\
 0,& |s|\ge L.
 \end{cases}
\]

For each fixed collar, both members have exactly the same spatial geometry,
original mouth locations/radii, boundary clocks and matching data. The waist
geometry is unchanged; all jets match the round ultrastatic exterior.
The local joined metric is \(C^\infty\) in regular charts. For positive
fixed collar width, radius and lapse stay positive, and the independent
curvature components are bounded. No limit \(\delta\to0\) is asserted.
This removes the differentiability qualification, **not** the missing
constitutive law, global Alice extension or dynamical AP-R0 checks.

The same warped-product stress formulas apply, but now
\(\kappa_5\rho=3(B-A)\) is generally nonzero in the collar. Its value is
identical across the two lapse choices. This support is explicitly included;
the smoothed slice is not claimed to be vacuum or everywhere scalar-flat.
Independent coordinate curvature and radial conservation tests cover the
collar and its exterior continuation. The supplement recomputes all six
scalar transfers for each pair, with refinement and determinant checks.
These are additional chosen regular geometries, not a derived collar scale.

| Collar fraction | Transit, epsilon=0 | Transit, epsilon=0.1 | Reduction |
| --- | --- | --- | --- |
| 0.1 | 0.841470984808 | 0.797579345456 | 5.21606% |
| 0.2 | 0.841470984808 | 0.802456194273 | 4.63650% |

All twelve paired mode comparisons differ. The largest refinement difference
is `5.94e-10`, and the largest determinant error is `1.76e-12`. A pilot at
`rtol=1e-11` failed the unchanged `1e-8` absolute refinement gate for the
high-gain `c=0.2, epsilon=0.1, ell=1, omega*S=0.5` transfer (`1.01e-7`).
Production uses `rtol=2e-13`, checked against `5e-14`; the failed pilot and
the precision change are disclosed in the supplement archive.

### Static-geometry admission test, with its scope fixed

For the specified static family (including these collars), at the waist
\(r=b,r'=0,r''=1/b\), hence exactly

\[
 \boxed{\kappa_5(\rho+p_s)=-3/b^2=-56.7853877\quad(R=1,a=0.5).}
\]

A candidate claiming to recover **this static geometry** must supply this
total radial null contraction (after any geometric modification is moved
to an effective source), as well as all other field equations. Merely
allowing "at least" this much NEC violation is not sufficient. For a more
general static minimum, the value is \(-3r''(0)/b\); for an evolving or
trapped throat, time-dependent terms enter. The boxed number is therefore
not a universal rejection threshold for all possible BAM focus laws.

The projected four-dimensional Weyl stress of #168 is not five-dimensional
matter support: that Tangherlini bulk is Ricci-flat and has zero physical
five-dimensional stress. Its regular exterior/horizon does not remove its
singular interior. Recovering the static open five-dimensional benchmark
requires genuine effective bulk support or modified bulk equations. We
correct the relevant #168 milestone in README/THESIS; the broader historical
dependency audit remains a separate follow-up.

## 5. What the cited mathematical papers contribute

We read the precise LaTeX statements at OpenAI/math commit
fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb. The repository says its manuscripts
are at different verification stages. This is an applicability audit, not
independent verification of the long proofs or external confirmation of BAM.

### Charged reduction in spatial dimensions at least four

Paper [1] explicitly assumes a neutral companion theorem. Its exterior is
oriented, complete with compact boundary and exactly one asymptotically
Euclidean end. Its matter constraint requires \(\mu_m\ge|J_m|\), its normalized
electric field is divergence-free, and the boundary satisfies its specified
weak trapping conditions. With the **full enclosing-area infimum**
\(A_{\min}\), \(X=(A_{\min}/\omega_{n-1})^{(n-2)/(n-1)}\) and
\(m=\sqrt{E_{\rm ADM}^2-|P_{\rm ADM}|^2}\), it states

\[
 m\ge|Q|,\qquad X\le m+\sqrt{m^2-Q^2}.
\]

The polynomial bound \(m\ge(X+Q^2/X)/2\) is the equivalent outer-branch
statement for \(X>|Q|\), not an unconditional substitute at small \(X\).
For a 5D bulk, \(n=4,\ \omega_3=2\pi^2,\ X=r_A^2\).
Its connected equality theorem adds MOTS, outermostness and
outer-area-minimizing conditions. It identifies an exterior and its regular
horizon attachment with Reissner-Nordström-Tangherlini, explicitly making no
claim behind the boundary, at extremality, or for disconnected equality.

This can constrain a proposed GR exterior after its hypotheses are checked.
It supplies neither a regular interior nor a radius evolution law.
A closed BAM bulk has no ADM end. A cavity mass is not automatically ADM mass;
a non-orientable region does not automatically meet orientability; an
effective support sector cannot silently be assumed to satisfy the energy
conditions. Passing an exterior inequality does not certify a regular focus.

### Enclosing area, charge and rigidity

Paper [2] separates three- and four-spatial-dimensional statements and the
additional rigidity assumptions. In \(n=4\) its neutral bound is
\(m\ge\frac12(A_{\min}/2\pi^2)^{2/3}\).
An arbitrary neck area or apparent-horizon area cannot simply replace the
required enclosing infimum. The paper's small-branch anti-de Sitter result
has different asymptotics and restrictions. These distinctions aid an exterior
audit; none is an inner/outer scattering law or nonsingular time evolution.

### Metric Blaschke theorem

Paper [3] states that a connected closed smooth Riemannian manifold with
\(\operatorname{inj}M=\operatorname{diam}M\) is, up to scale, a standard compact
rank-one symmetric space. This is a global metric condition, not a consequence
of spherical topology. It is also **not equivalent to point-antipode focusing
for every model in its conclusion**: for \(n\ge2\), the cut locus in round
\(\mathbb{RP}^n\) is \(\mathbb{RP}^{n-1}\). Standard complex projective space
likewise generally has a positive-dimensional cut locus. Cut locus, conjugate
locus and geodesic return period must be distinguished.

In a static spacetime \(-N^2dt^2+h\), null spatial paths follow the optical
metric \(N^{-2}h\), not generally \(h\). A dynamical spacetime need not admit
one fixed optical metric. A theorem about closed Riemannian metrics therefore
does not establish refocusing in an evolving Lorentzian bulk with handles.
Exact ray focusing also does not determine finite-bandwidth wave dispersion.

Its useful role is a rigidity benchmark and an assumption check: measure
arrival spread, focal displacement and finite-mouth capture under perturbation.
Adding a handle changes the geometry/topology; perfect Blaschke behavior must
not silently remain an independently free property.

## Dependency ledger

| Item | Status |
| --- | --- |
| Original neck under vacuum evolution | Singular waist endpoint after proper time \(b\); exterior-independence proved here only for \(X_m>1\) |
| Einstein dynamics outside the focus | Assumed classical law; support sector must be specified |
| No physical singularities in accepted events | Additional BAM postulate |
| Ideal round antipodal ray refocus | Background/control, with explicit support and optical metric |
| Positive-radius handle and Alice lift | Geometric and bundle inputs until derived |
| Regularity and conserved flux | Necessary checks, insufficient to select response, including the smooth-collar controls |
| Unique focus action or constitutive law | **Missing** |
| Well-posed nonlinear event and robust continuation | **Uncomputed** |
| Existing spectra/transport | Conditional on their background and boundary rules |
| Direct EGB mouth identification (#330) | Different interior: no finite neck; does not inherit neck transit \(2S\), waist \(b=R\sin^2a\), or neck-based `finite_mouth` spectra/DtN operators; transport requires new junction matching |
| QFT-on-background calculations | Use QFT ingredients; do not derive those ingredients |
| Quantum event, creation probabilities, Born rule, action quanta | **Not established here** |

## 6. Next executable gate and pre-stated failures

Before a caustic-event production run, publish one complete \(S_f\) or
equivalent closed system, all parameters and its region of operation.
Specify a covariant rule locating the matching worldtubes, the initial core
state, geometric/bundle identifications and every new length or time scale.
It must determine the metric/support response and field matching, rather
than impose a curvature inequality or choose a lapse after seeing results.
Calibrate proposed gates against historical controls before freezing them.

For a candidate claiming the static benchmark, first recover its regular
static solution, exact waist NEC contraction and junction ledger. A different
geometry must declare and test its own field-equation requirements. Then evolve
a finite, nonspherical gravitational-wave packet with the Einstein constraints
solved. Exact spherical symmetry would eliminate the intended radiative
degrees of freedom. Use the native bulk or label a reduction as a proxy.
Include the response of any existing mouth. A change in observable throat
access must not be relabelled as creation of full spatial topology.

For each candidate, fix these failures before interpreting a signal:

- No regular constraint-satisfying solution or no predictive continuation
  rejects that candidate in the stated regime.
- Resolved curvature divergence or geodesic termination violates the postulate.
  Failed numerical convergence is **unresolved**, never evidence of an event.
- Unbalanced flux/support work or constraint drift invalidates the calculation.
- Multiple admissible laws with different outputs and no selection principle
  leave the postulate underdetermined.
- Focusing, collapse thresholds, ordinary scattering or a classical topological
  sector alone do not establish a quantum event. Its discriminator must be
  stated separately, without inserting it into a quantized boundary rule,
  threshold or stochastic sampler.

This PR completes a foundation and computability audit. It does not claim
the missing focus law is derived. The measured ambiguity is a concrete target
for that law to eliminate.

Follow-up: the [first explicit action admission test](focus_action_admission.md)
now checks a closed ghost-scalar control, proves its original-mouth support
obstruction and reproduces the radial instability of its alternative regular
solution. The candidate is rejected under the stated gates; a viable focus
law remains missing. The [historical foundation claims audit](historical_foundation_claims_audit.md)
also reopens the pair-energy normalization and records the quantum-readout,
bundle and regular-event dependencies separately.

## Sources and reproduction

The source ledger pins the exact blobs read in the cited repository.

| Ref | Inspected source section | Git blob |
| --- | --- | --- |
| [1] | Charged reduction: build/sections/01-statements.tex | d10cafb1cfd94a9fe6601130542d0321f92e849d |
| [2] | Enclosing area: build/sections/00-introduction.tex | 83c0184a15b8f11968a755cb45046de44650d3df |
| [3] | Blaschke: build/sections/00-introduction.tex | 28d7f418fba62d70ea5d11ffe1f656097af2d635 |

- [1] [Charged reduction: precise statements](https://github.com/openai/math/blob/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb/preprints/A-Charged-Reduction-of-the-Spacetime-Penrose-Inequality-in-Spatial-Dimensions-at-Least-Four-October-5-2026/build/sections/01-statements.tex).
- [2] [Enclosing area, charge and rigidity: introduction](https://github.com/openai/math/blob/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb/preprints/Spacetime-Penrose-inequalities-enclosing-area-charge-and-rigidity-October-5-2026/build/sections/00-introduction.tex).
- [3] [Metric Blaschke theorem: statements](https://github.com/openai/math/blob/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb/preprints/The-metric-Blaschke-theorem-September-23-2026/build/sections/00-introduction.tex).
- [4] Harte and Drivas, [Caustics and wave propagation in curved spacetimes](https://arxiv.org/abs/1202.0540).
- [5] Chamseddine and Mukhanov, [Resolving Cosmological Singularities](https://arxiv.org/abs/1612.05860).
- [6] Simpson and Visser, [Black-bounce to traversable wormhole](https://arxiv.org/abs/1812.07114).
- [7] Bernal and Sánchez, [Smoothness of time functions and metric splitting](https://arxiv.org/abs/gr-qc/0401112).

This is an analytic counterexample with numerical verification, not a held-out
preregistered event experiment. The probe saves transfer matrices, source
hashes, refinement errors, constraints and curvature diagnostics. Units use
\(R=1\) and \(\kappa_5 T\) for stresses.

    OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.focus_regularity_probe
    python -m pytest tests/test_focus_regularity.py -q

All 11 focused tests pass under NumPy 2.3.5 and 2.2.6.
Tests compare the effective stress and Kretschmann scalar with a separate
coordinate-Christoffel/Riemann calculation, check radial stress conservation,
compare the static monopole with quadrature, refine transport integration,
check Wronskian conservation and verify scale covariance. Historical producers
and archives are unchanged.

Review supplement (original archive preserved):

    OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.focus_regularity_review_probe
    python -m pytest tests/test_focus_regularity.py tests/test_focus_regularity_review.py -q

The new archive lives in `runs/20261011_focus_regularity_review/result.json`.
The combined 22 tests additionally check vacuum proper-time quadrature,
geodesic integration, curvature in interior/exterior charts, the causal
countercontrol, smooth collar stress/matching and conservative transport.
