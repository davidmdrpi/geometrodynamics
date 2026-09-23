# Pre-registration: supported fields to the mouth mixing operation

Date: 2026-09-13.
Baseline: merged main `4cd86541d3b836b35561b0c4a3a54629d28851cd` (#300).
This document is to be committed and published before new implementation,
symbolic verification, numerical controls or response measurements for this
round. Source inspection and the analytic arguments below precede the freeze;
they are disclosed prior reasoning, not prospective discoveries.

## 1. Question and independent outcomes

Can the supported Einstein--four-scalar fields generate a physical mouth
deformation that coherently converts the throat's k=-1 and k=+1 modes,
without prescribing the missing mixing operator?

This connects the supported classical histories of #293--#300 to the
apparatus question of #238--#241. It does not identify those models by name
or notation alone. Five outcomes will be reported independently:

1. Whether an inherited bulk-to-mouth geometric and variational map exists.
2. Whether the reference lattice mixer survives a covariant geometric test.
3. Whether the complete quartet has an angular scalar-intensity source.
4. Whether a physical mixing response follows from the supported fields.
5. Which preparation, drive and apparatus inputs that response requires.

Noncommuting classical mode-conversion matrices are not canonical quantum
commutators. A successful conversion mechanism does not derive event
frequencies, a history measure, hbar, Phi, or the Born rule.

## 2. Prior state, including qualifications of the previous proposal

### P1. The exact quartet cancels its scalar-intensity modulation

Within #300's ansatz, with x.x=1,

    phi=B(q)x/A, B=q0 I+sum qi S_i,
    S_i S_j+S_j S_i=-2 delta_ij I,
    B^T B=Q I, Q=q.q,
    sum_alpha phi_alpha(x)^2=Q/A^2.

Consequently every nonconstant Fourier coefficient of the summed field
square on any closed curve x(chi) in S3 vanishes. In particular the m=2
coefficient vanishes. This is an exact prediction to verify algebraically,
with quadrature only as a check. It applies at arbitrary q and A>0, hence
at every time of #300's solution, without requiring a new integration.

Individual component squares can carry angular harmonics. Replacing the
sum by selected components changes the coupling and must be ledgered.
The prediction excludes only a source proportional to the summed scalar
intensity in this ansatz. It does not exclude derivative stresses, metric
anisotropy, localized perturbations, or different allowed preparations.

### P2. The previous ellipse is a specified lattice operator

The measurement capstone classifies the operator, occupation, and setting
orientation separately. Its `fiber_H` uses links

    R_(j+1/2)=1+epsilon cos(2 pi m (j+1/2)/N+varphi),
    t_(j+1/2)=R_(j+1/2)^(-2),

with a weighted graph Laplacian and the flat site inner product. At N=8
the reported first-order k=-1 to +1 matrix element for m=2 has magnitude
`(2-sqrt(2))*abs(epsilon)`, with a convention-dependent phase. This
coefficient is a reference-model prediction, not a result to refit.

The separate mouth-local potential model reports coefficient 1/sqrt(6).
It includes a radial envelope and different normalization. These two
coefficients are not interchangeable and need not agree.

Source inspection does not establish that the graph hopping prescription
is the wave operator obtained by varying a covariant mouth action. A
position-dependent material stiffness could realize that graph model;
calling the stiffness a metric is an additional identification to test.

### P3. Intrinsic-circle coordinate control, derived while drafting

For a static isolated circle with metric

    ds^2=R(chi)^2 dchi^2,
    R(chi)=R0[1+epsilon cos(m chi+varphi)], abs(epsilon)<1, m>=1,

arclength s=int R dchi has circumference 2 pi R0. A free scalar's action
has kinetic weight R and spatial stiffness 1/R:

    S=1/2 int dt dchi [R (partial_t u)^2-(partial_chi u)^2/R].

Thus its spatial generalized eigenproblem is

    -partial_chi(R^(-1) partial_chi u)=lambda R u.

It is isometric to the uniform circle and has spectrum k^2/R0^2. No
physical splitting of the +/-1 pair follows from this static intrinsic
metric alone. With kinetic matrix W and stiffness K, the degenerate
first-order calculation must use `delta K-lambda0 delta W`, not delta K
alone. Predict zero first-order splitting, including the off-diagonal
element in the unperturbed +/-1 eigenspace.

This is an analytic prior and a narrow control, not a verdict against an
embedded elliptic mouth. Extrinsic curvature, transverse modes, material
structure, lapse, connections, or boundary data can distinguish a real
mouth deformation. Those contributions must come from the actual model.
A time-dependent reparameterization also transforms lapse/shift and time
derivatives; changing R(t,chi) while holding those fixed is not this static
gauge control and must not inherit its conclusion automatically.

### P4. The existing classical support does not by itself define a mouth

#300 is a homogeneous 3+1-dimensional reduction on S3. It contains no
excised-mouth boundary, resolved throat, or detector material. Its bulk
action is the Einstein action plus four real conformal scalars:

    S_bulk=int sqrt(-g) [F R/2-Lambda/kappa
                        -1/2 sum_alpha (grad phi_alpha)^2],
    F=1/kappa-sum_alpha phi_alpha^2/6.

Boundary terms are needed when introducing boundaries. In particular,
the nonminimal F R coupling prevents blindly copying vacuum Israel
conditions and adding the scalar energy afterwards. Any proposed map must
derive the appropriate metric and scalar matching conditions, including
the relevant F and normal-derivative terms.

The expected repository outcome is an unspecified bulk-to-mouth map, not
a field-theoretic no-go. Finding a complete inherited map would overturn
that expectation. A missing map does not prevent completing P1--P3.

## 3. Fixed source inventory and notation discipline

Read implementations and their callers, not only status strings. Record
path, symbol, defining equation, provenance, applicable dimensions,
boundary conditions, supplied stress, and connection to the four fields.
The mandatory baseline inventory is:

| Source | What must be audited |
|---|---|
| `geometrodynamics/waves/nonlinear_supported_tt.py` | exact supported action, constraints, quaternion/coframe convention |
| `geometrodynamics/waves/coupled_multiplet_response.py` | full metric-dependent improved stress |
| `geometrodynamics/waves/reciprocal_scalar_tt.py` | variational source and reciprocal response |
| `geometrodynamics/waves/physical_throat.py` | scalar-flat spatial gluing; no supplied Lorentzian evolution |
| `geometrodynamics/waves/areal.py` | finite-mouth geometry, boundary channels and matching inputs |
| `geometrodynamics/waves/throat_operator.py` | self-adjoint boundary family versus selected physical boundary law |
| `geometrodynamics/waves/finite_throat.py` | resolved tube action, flux, and chosen geometry |
| `geometrodynamics/tangherlini/traversable_throat.py` | 5D supported benchmark and required stress |
| `geometrodynamics/shells/junction.py` | spherical Einstein thin-shell assumptions; applicability to m=2 |
| `geometrodynamics/transaction/network.py` | supplied ports, clock offsets and transmission data |
| `experiments/closure_ledger/throat_order_field_probe.py` | additional effective order field and unprovided metric coupling |
| `experiments/closure_ledger/throat_action_derivation_probe.py` | imposed action quantum versus a local field action |
| `experiments/closure_ledger/minimal_mixing_interaction_probe.py` | mouth-local m=2 potential, normalization and carrier assumptions |
| `experiments/closure_ledger/throat_apparatus_pointer_probe.py` | carrier preparation, quantum probability inputs, discrete charge |
| `experiments/closure_ledger/sigma_z_readout_capstone_probe.py` | graph kinetic measure, radius-to-hopping map and phase convention |

Follow direct references and search the repository for additional candidates;
archive the search commands and the exact file/commit inventory. Failure to
find a map is a statement about this audited baseline, not all possible BAM
extensions. A path that only mentions a throat is not a coupling.

Keep distinct: S3 harmonic degree, a local S2 angular multipole, the Hopf
fiber generator, the throat lattice winding k, and component-space rotations
of B. A map between any two must state coordinates, their periods, lift,
field transformation, antipodal action, and kinetic inner product.

The N=8 graph conserves a cyclic charge when wrap-around transitions are
included. It must not be certified by an integer winding commutator that
silently deletes those transitions. Continuum integer charge requires a
separate continuum definition. Likewise the 4D conformal coupling 1/6 and
the 5D coupling 3/16 cannot be interchanged without a derived reduction.

## 4. Geometry and response gates

### G. A complete bulk-to-mouth map

An affirmative map must specify the mouth as a physical surface or defect,
its embedding or construction, induced metric and relevant extrinsic data,
the fiber and time normalization, field restrictions, and a variational
interface law. It must either match a nonsingular geometry or supply the
surface action/stress it requires. A time-symmetric spatial constraint
solution is not a dynamical junction law.

If this map is missing, do not invent a surface stiffness, identify a
tensor component with ellipticity by fiat, or glue a 5D vacuum/supported
benchmark onto the 4D scalar solution by matching notation. Record the
missing equation or datum. Conditional operator controls still proceed.

### O. A covariant mode-conversion operator

Derive the full quadratic wave action in the candidate mouth geometry,
including kinetic weights, curvature coupling, connections, and interface
terms. Retain relevant transverse and other winding modes before any
two-mode restriction. Complex mode notation here denotes classical waves.

Report K, W, the physical normalization and the complete degenerate block.
For a static perturbation use delta K-lambda0 delta W. For a moving mouth
include the time-dependent basis/kinetic terms instead of importing that
static formula. An isolated off-diagonal matrix entry in an arbitrary
coordinate basis is insufficient: verify an invariant frequency splitting
or a conversion/scattering response relative to physically fixed ports.

Passing the old graph control alone yields only a lattice-control result.
If the covariant intrinsic-circle test cancels that apparent mixer, report
the cancellation at its exact scope. If an embedded-mouth action retains
a nonzero mixer, identify the term absent from the intrinsic-circle model.
Neither result is a universal no-go or a quantum probability derivation.

### R. A response generated by the supported fields

First examine the available bulk solution symbolically: the scalar-square
channel of P1 and the metric/stress channels must be reported separately.
Determine whether the inherited map supplies a forcing equation for a
physical m=2 mouth degree of freedom. Derive its source and reciprocal
bulk force from the same action. Reconstruct the full field/interface
residuals; a tensor projection alone cannot certify the response.

A prescribed static ellipse verifies operator sensitivity only. An
arbitrary initial anisotropic bulk metric is allowed preparation data,
but is not spontaneous generation of a mouth or a derived population law.
A nonzero source without a solved response is `SOURCE_TERM_ONLY`.

This freeze authorizes the source audit and analytic operator derivation.
If a complete inherited map permits a new coupled boundary-value or
evolution experiment, publish a prospective numerical addendum specifying
that equation, domain, data, grid, tolerances and residual gates before
solving it. Those quantities cannot honestly be fixed before the map is
known. No response sample may precede that addendum. It supplements this
freeze and cannot remove an original failed gate or replace a missing map.

### A. Preparation and apparatus accounting

Ledger independently the number and type of fields, initial amplitudes and
phases, mouth geometry, external drive, orientation, and coupling constants.
Classify the result as an inherited prepared-field solution, an added
apparatus/boundary interaction, or an externally prescribed deformation.
Keep operator existence, excitation/occupation, persistence, and orientation
control separate. Do not translate a classical amplitude into the older
carrier's quantum occupation nbar without deriving that normalization.

No preferred detector angle is required: orientation can be ordinary
preparation data. What must be derived is how the apparatus implements the
operation and how it responds back on the fields. Event frequencies would
additionally need a preparation measure and a physical record/basin map;
neither is supplied by a quadratic intensity identity.

## 5. Frozen controls and numerical reporting

No coupled mouth evolution is included before the section 4 addendum.
The following inexpensive checks are fully specified now. Use a distinct
random seed `2026091401`. Preserve all samples and failures.

1. **Quartet identity:** exact symbolic B^T B-Q I and Fourier cancellation;
   64 normal random q draws, A in {1,2}, and all six coordinate great
   circles x=e_i cos(chi)+e_j sin(chi). Uniform periodic quadrature with
   N in {32,64,128}; check harmonics m=1,2,3. Normalize by max(1,Q/A^2);
   numerical residual must be <1e-12. As a sensitivity control compute
   phi_0^2 for q=(1,0,0,0), A=1 on the (0,1) circle: its m=2 complex
   Fourier coefficient is 1/4, not zero. This changes the readout channel.
2. **Inherited graph:** reproduce `fiber_H` independently from its links,
   at N=8, m in {1,2,3}, varphi=j*pi/4 for j=0,...,7. Use the central
   epsilon derivative at h in {1e-2,5e-3,2.5e-3}. Compare with the exact
   derivative of the assembled graph; normalized error <1e-4 at finest h.
   Check the m=2 magnitude coefficient 2-sqrt(2), phase convention, and
   m=1,3 zero first-order controls; retain complex entries, not only norms.
3. **Covariant circle:** R0 in {1,2}, epsilon in {0,+/-0.05,+/-0.1},
   m in {1,2,3}, varphi in {0,pi/4,pi/2}. Derive the K,W cancellation
   exactly. Independently solve in arclength and verify the mapped modes
   u_k(chi)=exp(i k s(chi)/R0), k=+/-1,+/-2, by differentiating their
   flux R^(-1)u_k' and checking the differential equation pointwise at
   N in {32,64,128} chi nodes; normalized residual <1e-10. This is an
   independent formula check, not a finite-difference convergence claim.
   Check circumference and W-normalized overlaps by periodic quadrature
   (<1e-10 on the finest grid). Exact cancellation, not a tiny residual,
   is required for the zero-mixing verdict.
4. **Physical positive control:** a prescribed scalar potential
   V=v cos(m chi+varphi) on a uniform R0=1 circle, flat normalized Fourier
   modes, m=1,2,3 and the graph phase grid. At m=2 derive
   <+1|V|-1>=v exp(i varphi)/2; m=1,3 vanish. Use v in {0.05,0.1},
   the same quadrature grids and relative residual <1e-12. This is an
   added potential, not a prediction of the supported bulk action.
5. **Implementation-failure controls:** missing W, corrupted raw matrix
   elements, nonfinite data, a replaced symbolic zero certificate, and
   missing inventory entries must fail their own evidence gates. A
   passing control flag without raw evidence is never sufficient. A
   prescribed ellipse must never pass the inherited-map/field-response
   gate solely by making the operator test pass.

For second-order central differences, report consecutive error ratios and
require [3.5,4.5] when both errors exceed 1e-10. Report exclusions below
that floor. Use relative error ||a-b||/max(1,||a||,||b||); report absolute
errors too. Thresholds apply to each required case, not a mean. If exact
predictions or numerical gates fail, preserve the failure and diagnose it;
publish any changed grid or threshold prospectively as an extension.

## 6. Verdict contract and dependency isolation

There is no single all-pass headline. Every verdict includes its evidence
paths, model, assumptions and scope. The primary fields are:

| Field | Licensed outcomes |
|---|---|
| `bulk_mouth_map` | `INHERITED_MAP_DERIVED`, `BULK_MOUTH_MAP_UNSPECIFIED`, `SCOPED_MAP_OBSTRUCTION_PROVED`, `UNRESOLVED` |
| `operator_status` | separate `REFERENCE_LATTICE_REPRODUCED`, `INTRINSIC_CIRCLE_MIXING_CANCELS`, `PHYSICAL_MOUTH_MIXING_DERIVED`, or `UNRESOLVED` entries |
| `quartet_intensity` | `NO_ANGULAR_SOURCE_IN_EXACT_QUARTET`, `PREDICTION_REFUTED`, `UNRESOLVED` |
| `field_response` | `INHERITED_FIELD_RESPONSE_DERIVED`, `SOURCE_TERM_ONLY`, `CONDITIONAL_APPARATUS_RESPONSE`, `BLOCKED_BY_UNSPECIFIED_MAP`, `SCOPED_RESPONSE_OBSTRUCTION_PROVED`, `UNRESOLVED` |
| `preparation_status` | separate classifications for fields, geometry, drive, phase, orientation and population; absent mechanism = `NOT_DERIVED` |

Reference-lattice, intrinsic-circle and physical-mouth operator statements
can coexist: they concern different actions, not competing votes on one
number. Conditional apparatus response must name each added interaction.
`BLOCKED_BY_UNSPECIFIED_MAP` is not zero response; missing measurements are
null, never zero. A failed numerical construction licenses only unresolved
status, not an obstruction. A scoped obstruction needs an analytic proof.

Required evidence dependencies:

| Gate | Targets |
|---|---|
| `source_inventory` | map, physical operator, field response, preparation |
| `quartet_identity` | quartet intensity; any response claim using it |
| `graph_control` | reference lattice |
| `kinetic_measure_and_coordinate_control` | intrinsic circle, physical operator, field response |
| `physical_potential_control` | physical operator, field response |
| `map_and_interface_equations` | affirmative map, physical operator, field response |
| `reciprocal_action_and_full_residuals` | affirmative field response |
| `preparation_ledger` | physical operator, field response, preparation |
| `raw_evidence_and_failure_paths` | each target whose evidence is affected |
| `scope_and_provenance` | all targets |

Known missing downstream prerequisites leave the upstream quartet and
operator-control results intact. Unknown gate names fail all affirmative
claims. The eventual probe must actually exercise damaged evidence through
the verdict function, and exit nonzero on any required verification failure.
A correctly established repository gap is an allowed scoped audit result,
not a claim that the requested mixing response exists.

## 7. Deliverables and exclusions

This PR contains this freeze only. Implementation will retain the published
freeze SHA and baseline, archive source references and raw control evidence,
and add a results document separating prior predictions from discoveries.
Any numerical response addendum must be public before that response is run.

The results must answer what operation the audited field content actually
supplies, or identify the precise missing map/action that prevents an answer.
No CHSH optimization is needed. No Born state, quantum carrier occupation,
chosen event-frequency law, or history weight may be imported to certify a
classical result. The winding-to-spin map and closure-to-field map remain
separate until derived. No claim about generic inhomogeneous stability,
throat creation, universal non-readability, or operational causality follows
from the homogeneous persistence theorem or from these operator controls.

The intended advance is a defensible connection, or a localized obstruction,
between the classical field construction and the identified apparatus.
It is not another selection of Phi by agreement with the quantum answer.
