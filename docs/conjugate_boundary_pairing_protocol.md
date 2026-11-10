# Boundary pairing of a dynamically responding conjugate mouth pair

## Question and publication scope

Can a classical variational problem for one handle and its two responding
mouths supply an oriented boundary pairing, a positive physical readout and
a consistent rule for coherent joining? The target is a physical explanation
of quantum structure, not another fitted angular correlation.

Baseline: [`3293d5bad3c8e103c8514ad85afc0128050004e4`](https://github.com/davidmdrpi/geometrodynamics/commit/3293d5bad3c8e103c8514ad85afc0128050004e4),
the merge of PR #350. This PR prepares the next experiment. It contains this
prospective protocol and executable **baseline qualifications**, not a solved
dynamic GR mouth pair. No new source, support law, detector ensemble or
probability measure is silently promoted to a geometric result.

The static action, transport algebra and counterexamples below were known
before this protocol. Reproducing them is calibration, not discovery. The
protocol is committed before the new qualification code and numerical output.
The dynamic experiment is a later implementation against this protocol;
any changed hypothesis or threshold needs a dated amendment before its run.

## Prior PR findings that change the experiment

All source links in this table refer to the pinned baseline, not the current
heads of old PR branches. PR #349's branch continued after its merge; #350
imported the final sources and records. Those final files are the starting point.

| Prior finding | Consequence here |
| --- | --- |
| [#279](https://github.com/davidmdrpi/geometrodynamics/pull/279), [finite-mouth transport](https://github.com/davidmdrpi/geometrodynamics/blob/3293d5bad3c8e103c8514ad85afc0128050004e4/docs/finite_mouth_bundle_transport.md): the antipodal bulk handle is orientable; the brane and normal bundle twists, quotient and Pin lift are different data. Relative to the mouth Hopf fibre, transport is anti-linear. | Name the manifold, bundle, seam and representation separately. A two-port matrix is indexed by endpoints, not by spin components. Never identify those two-dimensional spaces merely because their sizes agree. |
| [#279's Born investigation](https://github.com/davidmdrpi/geometrodynamics/blob/3293d5bad3c8e103c8514ad85afc0128050004e4/docs/classical_born_rule.md): spin symmetry and Malus intensity do not determine event frequencies. | A positive flux or power remains a response observable until a record-producing detector and preparation ensemble are supplied. |
| [#280](https://github.com/davidmdrpi/geometrodynamics/pull/280): oriented closure gives the singlet-shaped law only with additional aggregation, sector and readout inputs. | Recover any relative sector coefficients from a single apparatus model, not separate normalizations chosen after inspecting an angle. |
| [#281](https://github.com/davidmdrpi/geometrodynamics/pull/281), including [review correction N11–N15](https://github.com/davidmdrpi/geometrodynamics/pull/281#issuecomment-5520069816): a holonomy trace is not the derived action; two classical quadratic readouts give CHSH `8 sqrt(2)/3` and `12 sqrt(2)/5`. | Vary an identified action. Quadratic homogeneity alone cannot select the readout, the Born rule or the quantum bound. |
| [#282](https://github.com/davidmdrpi/geometrodynamics/pull/282): phase and numerator windows have the same closure set but different conditional measures. | Specify preparation and any limiting procedure. Do not infer a measure from the set of solutions. |
| [#283](https://github.com/davidmdrpi/geometrodynamics/pull/283), [#284](https://github.com/davidmdrpi/geometrodynamics/pull/284), [#285](https://github.com/davidmdrpi/geometrodynamics/pull/285): equilibrium is conditional on a frame-restoring Hamiltonian; source readout and finite detector resolution change the causality question. | Record the support/preparation Hamiltonian and source-accessible observables. Detector marginals do not certify no signalling to the past. |
| [#286](https://github.com/davidmdrpi/geometrodynamics/pull/286): positive counting can yield both the quantum correlation and CHSH=4. | Oriented cancellation is not necessary for a fitted singlet curve. Positivity of separate probabilities is not strong positivity of a common history functional. |
| [#288](https://github.com/davidmdrpi/geometrodynamics/pull/288): independent composition selects no counting function; the common quaternion frame is an assumption; the on-closure product differential can have rank two. | Test coherent joining of interfaces with transported frames. Do not repeat a product-construction identity and call it physical selection. |
| [#341](https://github.com/davidmdrpi/geometrodynamics/pull/341): mixed scalar exchange is reciprocal but distributed in time. | Use complete encounter windows and independent mixed perturbations. A selected focus instant cannot represent the full interaction. |
| [#345](https://github.com/davidmdrpi/geometrodynamics/pull/345), [#349](https://github.com/davidmdrpi/geometrodynamics/pull/349), [two-port review](https://github.com/davidmdrpi/geometrodynamics/pull/349#issuecomment-6458618899): supported Maxwell scattering and loads come from one action; the waist load opens the throat; finite ports need exterior matching; packet routers discard complex scattering information. | Start with the full complex boundary map. Do not replace it with transmission fractions, mean delays or the conductor's opposite load sign. |
| [Einstein/support review](https://github.com/davidmdrpi/geometrodynamics/pull/349#issuecomment-6463315351), [#350](https://github.com/davidmdrpi/geometrodynamics/pull/350): initial constraints do not determine support dynamics; the same slice admits different radius evolutions. | An assumed positive mouth mass/stiffness is an effective comparator. It is not a derived GR support response. Keep that verdict separate from numerical validity. |

## Three pairings that must not be conflated

1. **Variational pairing:** the boundary term of the real action and its second
   variation. Its antisymmetrization gives a symplectic current; a Hessian is
   symmetric when its boundary conditions make it an on-shell action Hessian.
2. **Conjugate-mouth transport:** `J z = epsilon conjugate(z)`, with
   `epsilon=[[0,1],[-1,0]]` and `J^2=-I`, on the specified spin-frame bundle.
3. **Detector response:** a physical quadratic form in outgoing disturbances,
   followed, if available, by a rule for stable recorded events.

The symplectic form is not a positive norm. A real reciprocal port Hessian is
not automatically the spinor tensor epsilon. Inserting J in a coordinate
change can make epsilon appear algebraically without deriving a new physical
interaction. Such a result is `COORDINATE_PAIRING_ONLY`.

For normalized spinors with Bloch directions `s a` and `t b`, the known identity

    |z_s(a)^T epsilon w_t(b)|^2 = (1 - s t a.b)/2

gives the singlet probability after dividing by the four-channel sum, **if**
that squared pairing is the event weight with equal physical channel response.
This is `CONDITIONAL_SPINOR_IDENTITY`, never the experiment's success verdict.
A global minus sign of J cancels in that square. The protocol therefore needs
relative-path interference or another sign-sensitive observable to test a Pin
sector; the normalized two-mouth curve alone cannot distinguish it.

## Start from the inherited Maxwell action

Use `waves/maxwell_throat_response.py`, not the visualizer's neck. On the
specified supported bridge, one radiative master channel has

    S_EM = 1/2 int dt ds [lambda/N psi_t^2 - N/lambda psi_s^2
                         - lambda N ell(ell+1)/r_q^2 psi^2],
    r_q(s,t) = r_0(s) exp(q(t) h(s)),   h=[1-(s/S)^2]^3.

Here ell>=1. Static Wheeler flux is a separate sector. Two Maxwell
polarizations are integer-spin field channels; they are not a derived
spin-1/2 mouth state. The bridge in this reduction is the 3+1-dimensional
radial Maxwell comparator; the bulk spin-frame construction uses different
geometric data. An explicit intertwining/matching map is needed to connect
them. Duplicating the scalar equation and calling its copies a spinor is not
that map.

The spatial boundary variation is `-[j delta psi]_L^R`, with
`j=N psi_s/lambda`. Define outward fluxes `f=(-j_L,j_R)`. The outward
Dirichlet-to-Neumann map is `f=D u`; with fixed temporal endpoints,
`delta S_on=-f.delta u`. Thus the action's boundary Hessian is **-D**, not D,
under this convention. For two linearized solutions the spatial Wronskian is
`psi_1 j_2-j_1 psi_2`; both endpoint signs must be retained in joining.

At fixed lapse and intrinsic length the derived load is

    F_q = int ds lambda N ell(ell+1) h psi^2 / r_q^2 >= 0.

This is a waist-opening load, not an orbital force. A responding model needs

    S_total = S_EM + S_geometry+support + S_detector.

The latter two actions are not supplied by the static scattering calculation.
Stage 1 may use `S_support=int dt (M qdot^2/2-K q^2/2)` with explicit positive
M,K. It must say `EFFECTIVE_SUPPORT_ASSUMED`, retain its energy in the ledger,
and vary those parameters. It cannot close the GR or detector gates.

### A prospective order-counting prediction

For the reduced support comparator, expanding about psi=0, q=0 gives

    S_EM = S_EM,0 + int dt ds c(s) q psi^2 + O(q^2 psi^2),
    c(s)=lambda N ell(ell+1) h/r_0^2.

With zero initial mechanical perturbation and incident field amplitude A,
at fixed finite time away from singular response:

    psi=A psi_1+O(A^3),  q=A^2 q_2+O(A^4),
    W=A^2 W_2+O(A^4) for a quadratic field readout.

Consequently the field/mechanical mixed Hessian vanishes on the unpumped
background. A zero-field linear-response computation cannot demonstrate
selection by dynamically responding mouths. Around a nonzero solved pump,
the Jacobi equations have nonzero mixed blocks. Pump work and mechanical
channels must then be included; field-only scattering need not be unitary or
contractive. This is a conditional analytic prediction, not a numerical result.

## Staged experiment

### Stage 0: qualify the existing operators (this PR)

Run `python -m experiments.closure_ledger.boundary_pairing_baseline_probe`.
It produces a deterministic JSON report, exits nonzero on missing, nonfinite
or failed checks, and never declares the dynamic mechanism established.

- Recompute D and S independently on the supported bridge. Check the outward
  sign convention through `S=-(D-i omega I)^-1(D+i omega I)`.
- Split the same throat at an interior surface and at two interior surfaces.
  Eliminate shared Dirichlet data using flux continuity (Schur complement).
  Compare with the whole-throat map, for off-centre cuts as well as the waist.
  A Dirichlet pole or ill-conditioned elimination is unresolved, not clipped.
- Check quadratic load cross terms by polarization against inherited internal
  field integrals and against the independent scattering-phase derivative.
- Use right-port outgoing power for a positive Gram example. Its PSD property
  and vanishing third-order interference on a fixed background are structural
  controls, not evidence of a Born law or a common all-setting kernel.
- Keep a counterexample with four positive, normalized context-specific Gram
  matrices and CHSH=4. This prevents a false Tsirelson inference from per-setting
  positivity. The inconsistent common extension is not silently postulated.
- Reproduce the conditional epsilon identity including phases, antipodes,
  arbitrary common SU(2) rotations, unequal channel efficiencies and a
  sign-blindness control. No angle is used to calibrate a coupling.

### Stage 1: derive and solve a responding action

Before any trajectory, publish the complete chosen action, initial/boundary
data, support and source work, constraints, gauge conditions and observables.
Derive the field and mechanical equations from that same action. Establish
existence on the reported interval; stop before lapse, radius or coordinate
validity is lost. Do not prescribe q(t) and label it backreaction.

Start with one responding waist coordinate in the declared radial comparator.
Then introduce independent mouth displacement/shape coordinates only when
their inertias, support variations and field coupling are specified. One
waist mode is not two independently accelerating mouths. Opposite static flux
labels do not turn the radiative master field into a spinor.

Preparation grid, fixed here for the initial effective comparator: bridge
angle 0.3; lapse depths 0 and 0.4; neutral static flux; ell=1; characteristic
carrier omega=20 with a Gaussian temporal envelope of width 0.15 centred at
t=0.6, observed on [0,3]; A=0, 1/64, 1/32, 1/16, 1/8. Prepare both incidence
channels with phase differences 0, pi/2 and pi; include each one-sided source.
M=1,4 and K/M=1,4. Units are G=c=R=1 and the inherited master-field action
normalization. This grid tests sensitivity, not physical material constants.
Refine spatial/time resolution twice; finite-domain/exterior approximations
and incoming characteristic implementation must be frozen in an amendment
before execution. The Gaussian is a specified input, not a derived birth flash.

### Stage 2: obtain the causal tangent response and the boundary Hessian

Around each solved background, differentiate the coupled equations. For
boundary source variations u and interior/mechanical variations y, keep the
full block operator. Under a well-posed variational boundary prescription,

    H_eff = H_uu - H_uy H_yy^-1 H_yu.

Compute it independently by solving perturbations, by finite differences of
the full nonlinear solution, and by the second variation. Retain modes when
H_yy has a kernel; never hide a pseudoinverse or regularization.

A **retarded** response kernel is not generically the symmetric Hessian of an
ordinary one-copy on-shell action with fixed endpoints. State the boundary
prescription, compare retarded to adjoint advanced response where appropriate,
and distinguish the two objects. Retarded response alone fixes no fluctuation
ensemble or symmetric noise kernel. This prevents a response identity from
being reported as a probability-measure derivation.

Transport local frames across the physical seam. Classify the computed map as
real-linear, complex-linear or anti-linear relative to the inherited Hopf
complex structure before fitting an invariant tensor. If an actual two-spinor
bilinear B is obtained, test `U^T B U=B` and proportionality to epsilon, including
its magnitude and phase under independently changed endpoint frames. Do not
manufacture B by multiplying a scalar transmission by epsilon. The two-port
endpoint antisymmetric matrix is not the target spinor pairing.

### Stage 3: coherent joining, readout and distinguishability

Compute apparatus output disturbances v_alpha on one common solved background
with the same physical output channels and measure. If the action/detector
defines G, evaluate `D_alpha,beta=v_alpha^dagger G v_beta`. Derive positivity
on physical perturbations; do not replace negative eigenvalues with zero.
Record whether positivity comes from passive field flux, a stable reduced
sector or an extra assumption. Neither signed mechanical work nor net flux
is automatically positive.

Test direct versus joined solutions, including hidden mechanical output and
interface work. Change where the artificial interface is cut. Independent
device products and coherent cascades are different tests. Quantities true
by a Gram or product construction are calibration, not physical selection.

For interference use alternatives of the same linearized boundary problem.
Changing slits, sources or backreacted geometries changes the operator and is
not the same third-order-interference test. Report full nonlinear deviations
separately, with the amplitude scaling predicted above.

### Stage 4: conditional quantum consequences, only after physical matching

For Bell settings use common source preparation and one detector model, with
all reflected/transmitted/lost channels recorded. Test the four contexts
`a=z, a'=x, b=(z+x)/sqrt(2), b'=(z-x)/sqrt(2)` plus held-out directions.
Do not infer a shared history space from four separately positive response
matrices, or omit setting-dependent losses by coincidence normalization.

A common strongly positive decoherence functional reproducing the required
measurement marginals would support a Tsirelson inference. It is a target,
not an automatic consequence of passive scattering. Derive consistent
coarse-graining and joining before invoking that theorem. Quantum-measure
strong positivity is not by itself a reconstruction of all of QFT.

For event frequencies, supply an actual classical apparatus trajectory,
preparation ensemble, exclusive record rule and repeat-run frequencies. Any
common action per transition must be derived, not an inserted hbar or threshold.
Without these the result remains a response law. Test operational source
readout separately from detector no-signalling, including failed and lost runs.

## Registered gates and verdict fields

Stage-0 tolerances: transport/conditional identity residual 1e-12; DtN
reciprocity, whole-versus-cut maps and DtN-versus-S 2e-7 relative to
`max(1,norm(reference))`; load polarization and phase-force agreement 2e-6
on the same normalization; Gram/Hermiticity/I3 1e-12 relative; the PR-box
control must give CHSH=4 within 1e-12. Record the observed values, not just booleans.

Future numerical gates: full energy plus supplied work closes within 1e-5 of
nonzero prepared energy; tangent-versus-finite-difference output within 1e-4
relative for resolved nonzero outputs with a step-refinement plateau; cut
independence within 1e-4 with refinement; no claimed PSD eigenvalue below
`-1e-9 max(1,operator norm)` after numerical uncertainty. Report zero-signal
cases as unresolved for relative scaling, not as successful selection.

Every mandatory gate must be present, finite and true. Numerical failure or
missing physical input cannot produce a no-go theorem. Keep separate fields:

| Field | Permitted outcome examples |
| --- | --- |
| numerical validity | `PASS`, `UNRESOLVED_NUMERICALLY` |
| support dynamics | `EFFECTIVE_SUPPORT_ASSUMED`, `GR_SUPPORT_DERIVED`, `NOT_IMPLEMENTED` |
| spinor matching | `COORDINATE_PAIRING_ONLY`, `PHYSICAL_PAIRING_DERIVED`, `NOT_DERIVED` |
| positive response | `FIXED_BACKGROUND_FLUX_GRAM`, `PHYSICAL_COMMON_KERNEL_DERIVED`, `NOT_DERIVED` |
| event law | `CONDITIONAL_SPINOR_IDENTITY`, `DETECTOR_FREQUENCIES_DERIVED`, `NOT_DERIVED` |
| composition/causality | `QUALIFIED_STATIC_JOINING_ONLY`, `NOT_ESTABLISHED`, or a precisely scoped derived statement |

Stage 0's only success label is `BASELINE_QUALIFIED_DYNAMIC_EXPERIMENT_PENDING`.
Never aggregate it into a derived-Born, quantum-bound, or GR-completion label.

## Mathematical leads (not imported dynamical axioms)

- [Gudder, Hilbert Space Representations of Decoherence Functionals and Quantum Measures](https://arxiv.org/abs/1011.1694): representation once strong positivity is supplied.
- [Craig et al., A Bell Inequality Analog in Quantum Measure Theory](https://arxiv.org/abs/quant-ph/0605008): the common strongly positive joint functional and experimental marginals are load-bearing hypotheses.
- [Dowker and Wilkes, An Argument for Strong Positivity](https://arxiv.org/abs/2011.06120): compositional restrictions within the specified decoherence-functional framework; this does not derive tensor composition from BAM.
- [Aaronson and Ambainis, The Need for Structure in Quantum Speedups](https://arxiv.org/abs/0911.0996): a later operational test. A T-step affine-input linear transfer has amplitude degree at most T and quadratic readout degree at most 2T; a return resolvent hides repeated traversals. Define oracle access, boundedness and resource costs before applying an influence theorem.

The constructive target is an action-derived map that retains oriented cross
terms and composes through real interfaces. Even a partial result locating
which physical variable carries the antisymmetric pairing would advance this
program beyond another correct conditional angular identity.
