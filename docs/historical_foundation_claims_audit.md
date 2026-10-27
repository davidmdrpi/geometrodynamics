# Historical foundation claims: dependency audit

12 October 2026. Baseline `d42cb5e695f689dee254b51dbeabed3c21b6ad3c`.
This audits the event/regularity and quantum-foundation chain affected by
#328. It is **not** an exhaustive revalidation of the repository's mass,
flavor, QCD, thermodynamic or perturbative-QFT predictions. Historical
producers and their archives remain unchanged; the corrected current
interpretation is recorded here and in README/THESIS.

The audit distinguishes algebra within a supplied model, a chosen physical
identification, and an actual Einstein–focus evolution. A green legacy
`pass` flag establishes only the predicate its code evaluates.

## Findings by dependency

| Historical claim | What survives | Missing dependency / corrected status |
| --- | --- | --- |
| #55/#58 equilibrium implies electron rest energy and pair threshold | E=A/R+B R² has a minimum at (A/2B)^(1/3) | **Numerical normalization mismatch:** the implemented minimum is (3 alpha/4) m_e c², not m_e c². The threshold routine separately inserts the electron energy. |
| #166 antipodal focus triggers pair creation | Exact conformal finite-mode identity psi(chi,pi)=-psi(pi-chi,0) on the prescribed unit S³ | Nonlinear metric/core dynamics and the threshold law are not evolved. A finite smooth mode sum has finite pole limits; 1/sin²(chi) alone is not its physical energy density. |
| Hayward core is derived singularity avoidance | Chosen density/metric profiles give a regular static core and conditional stress relations | No focus action, formation calculation or stability/completeness proof follows from that static fit. SEC violation is not the null-convergence hypothesis of the Penrose trapped-surface theorem. |
| #168 global regular 5D embedding | Exact vacuum exterior and regular Killing horizon | Singular interior remains; corrected in #328. Projected 4D Weyl stress is not 5D support. |
| #63 C equals inner/outer swap | Opposite oriented integrals and the reflection involution | Two orientations of a cut do not by themselves supply two independent charged objects or an Alice bundle lift. Charge-conjugate field identification is extra data. |
| #66 scalar factorization derives a spacetime Dirac field | A ground-state Riccati factorization and partner spectra | A generic positive scalar operator admits such a factorization. A chosen spin bundle, Clifford principal symbol, covariant first-order equation and physical statistics do not follow from counting two mouths times two spin labels. |
| #170/#172/#196 Fermi statistics follows from metric topology | Conditional Pin/exchange algebra and Fermi-gas calculations with single occupancy | Pin-minus availability is not a unique Pin structure or a selected physical representation. Framing, history weights and particle identification remain inputs. A Fermi gas calculation cannot independently establish its imposed occupancy rule. |
| #200 cobordism completes dynamical pair creation / Pauli chain | Linking matrix diag(2,-2), homology, even framings and index-two handle construction | Topological cobordism is not a nondegenerate regular Lorentzian Einstein–focus solution. Avoiding indices 1 and 3 does not itself prove causal continuity or select pairs. A single +2 handle is also index two. Gauge charge assignment is additional. |
| #199 Bianchi identity selects guidance | A current/stress identity for a specified winding scalar and continuity checks in a reduced model | Bianchi enforces conservation on a solved Einstein–matter system; it does not derive that scalar from vacuum GR or determine a finite core's worldline. A nonsingular 5D core evolution is missing. |
| #198 Born rule derived | Continuity and equivariance under the specified phase-gradient flow; uniqueness within h(rho) | Selection of guidance, initial equilibrium, measurement records and general many-body statistics remain separate. Coarse-grained relaxation in one reduced run is not a universal attraction theorem. |
| Bell/CHSH from topology alone | Singlet projector calculation yields -cos(a-b) and 2 sqrt(2) | The code uses a tensor-product state, spin projectors and squared amplitudes. These are readout inputs, not an outcome of the uncomputed focus event. Binary outcomes and no-signaling alone do not enforce Tsirelson. |

## Reproduced normalization failure in #58

`pair_production_threshold_probe.test_T1_lowest_stable_configuration` sets

\[
 A=\alpha\hbar c/2,\quad R_*=\hbar/(m_ec),\quad
 B=A/(2R_*^3),\quad E(R_*)=3A/(2R_*)=(3\alpha/4)m_ec^2.
\]

Its test checks the **location and positive second derivative** of the
minimum; it never checks equality to m_e c². The separate
`test_T3_threshold_2mc2` sets `E_single=M_E_C2_MEV` and doubles it.
Both legacy tests pass while their connecting physical assertion fails:

| Quantity | Reproduced value |
| --- | --- |
| E_min / (m_e c²) | 0.005473014426975 |
| Twice the implemented minimum | 0.005593409251 MeV |
| Separately asserted pair energy | 1.021997899992 MeV |

This is not a new physical pair-production prediction: the toy energy
functional is incomplete as an event law. It is a specific failure of the
claimed equality in the historical derivation. Do not repair it by silently
renormalizing A or B. Any additional energy term needs its own derivation
and matching ledger. Conservation of a specified charge can exclude an
isolated charged product, but does not uniquely select exactly two products
among all charge-neutral channels.

## Finite bandwidth removes the apparent pole in #166

The actual field in the producer is a finite sum of
sin(n chi)/sin(chi). Each term approaches n at the north pole and
(-1)^(n+1)n at the south pole. The denominator cannot be taken to zero
while treating the Dirichlet numerator as a nonzero constant.
Our fixed 16-mode example has north-pole value 19.8948948916772;
near-pole error decreases quadratically, while the exact antipodal identity
holds to 2.14e-14. The geometric-optics caustic and finite-band wave solution
are different limits. Neither this check nor #166 solves backreaction.

## Readout countercontrol for the Bell claim

The production `bulk_probability` is literally the squared modulus of
`bulk_amplitude`; `history.closure.enumerate_bell_branches` likewise starts
from amplitude squared before applying history weights. Keep those same
singlet amplitudes and change only the normalized probability map:

\[
 P_p(s_a,s_b|a,b)=|A_{s_as_b}(a,b)|^p/\sum|A|^p.
\]

At standard CHSH settings, p=2 gives 2.82842712474619, whereas p=4 gives
3.77123616632825=8 sqrt(2)/3. Both have exactly half-valued marginals by
the pair symmetry (sampled error below 2.3e-16). Thus the available geometry,
binary outcomes and no-signaling leave a readout choice. The p=4 map is a
mathematical countercontrol, **not** a realizable quantum measurement or a
constructed GR mechanism. The usual Hilbert-space/Born assumptions recover
the quantum bound; this audit does not refute that theorem.

## Cobordism and guidance do not close the event gap

A single +2-framed 2-handle on S³ times an interval already gives an
index-two cobordism to RP³. Therefore the assertion that a single-geon
construction necessarily uses index 1 or 3 is false as a topological
statement. This does not prove that either the single or paired construction
admits the required regular causal metric, field bundle or conserved charge.
Changing Cauchy-surface topology is incompatible with a smooth globally
hyperbolic product evolution; the new focus postulate must say what changes.
An abstract history construction and chosen sum-over-histories weights are
not a solution of that problem.

The Bianchi identity is also not a particle equation of motion by itself.
For example, on a circle with fixed positive rho(x), both J=0 and a nonzero
constant J satisfy stationary continuity. They give different velocities
J/rho. A specified scalar action may select a current; linking that current
to the worldline of an independently tracked gravitating core still requires
the core solution. Equivariance then preserves an initially chosen ensemble;
it does not by itself select initial equilibrium or detector outcomes.

## Source coverage and reproduction

The archive binds source hashes. The audited predicates/functions are:

- #58: `test_T1_lowest_stable_configuration`, `test_T3_threshold_2mc2`;
  #55: `test_T4_self_consistent_equilibrium` uses the same anchored minimum.
- #166: `conformal_frequencies`, `reduced_field`, `physical_field` and the
  declared inherited threshold.
- #63: `test_T3_orientation_reversal`, `test_T4_curvature_flips`.
- #66: `_superpotential`, `test_T1_dirac_factorization`, discrete partner
  spectrum and two-mouth/spin interpretation.
- #200: `test_T2_construction_arithmetic`, declared Morse indices and the
  booleans in `test_T3_dowker_sorkin_conditions_closed`.
- #198/#199: displayed reduced field equations, `_velocity`,
  `kk_identity_check`, and the declared guidance/equilibrium identifications.
- Bell: `bulk_amplitude`, `bulk_probability`, and history branch enumeration.
- Core: `blackhole/interior.py`'s explicit Hayward metric ansatz.

Only the cheap threshold, finite-band and readout checks are rerun here;
the full historical simulation suite is not presented as independently
reproduced. The source audit concerns logical dependencies, not replacement
of all historical results. Unreviewed spectrum/mass/QFT claims retain their
previous evidence status, not a new audit endorsement.

    python -m experiments.closure_ledger.historical_foundation_audit_probe
    python -m pytest tests/test_historical_foundation_audit.py -q

The next focus-law gate is in [focus_action_admission.md](focus_action_admission.md).
The first explicit candidate fails; its rejection does not repair or erase
any of the quantum-event dependencies listed above.
