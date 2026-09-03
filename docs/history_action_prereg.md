# Pre-registration: does a classical BAM action select the oriented branch?

**Status: frozen before any implementation code.** Seventh round of the
finite-mouth chain, following `docs/closure_current_prereg.md` (`f954e3d`)
and the merge of PR #280 (`9303b22`), whose three underived inputs
— branch aggregation, sector coefficients, readout — this round attacks.
Success criteria may change afterwards only by an explicit correction note
appended to this file.

## 0. Scope correction recorded first

**C1 — the trunk is bifurcated; this round inherits only one strand.**
`geometrodynamics/bulk/closure_current.py` imports `mouth_spin_frame` and
`closure_measurement` and nothing else. The Bell/singlet strand therefore
descends from the identification of the Hopf bundle with `P_Spin(S²)`
(round 3) and **not** from the antipodal scalar boundary condition, the
`η` sign, or the quotient-versus-double-cover choice: the two throat
transports enter the reduced loop only as the global factor `J² = −1`,
which cancels in the normalisation. Accordingly **no quotient or `η`
assumption may be carried into any probability claim of this round**, and
none of this round's verdicts may be stated as conditional on them.

**C2 — "no local functional has the closure set as its critical set" is
false and is not pre-registered.** `N²` has that critical set trivially.
The question is not existence but *derivation*: whether BAM's own
classical structures single out a particular functional.

## 1. Results established before freezing (analytic oracle)

Let `θ(x; u, v) = Ω/2`, with `Ω = 2 atan2(N, D)`, `N = x·(u×v)`,
`D = 1 + x·u + u·v + v·x`. Let `G` be the derived SU(2) closure holonomy
of round 6, `G = cos(Ω/2) + sin(Ω/2) x`. Define the **holonomy-trace
functional**

```
S_H[γ] = −½ Tr G = −cos θ = −D / sqrt(D² + N²).
```

The following are established analytically and were confirmed numerically
before this file was committed. They are frozen as oracles; the module
must reproduce them, and a disagreement is a bug in the module, not a
result.

* **O1.** `−½ Tr G = −cos(Ω/2) = −D/sqrt(D²+N²)` identically (`3.3e-16`).
  In particular `S_H` is smooth wherever `D² + N² ≠ 0`, i.e. away from
  `x = −u, −v`, even though the `arg` chart is not.
* **O2.** `dS_H = sin θ · dθ`, so on `S²` every regular closed history is
  critical, and — because round 6 established
  `∇_{S²}θ|_{N=0} = (u×v)/D ≠ 0` — there are **no other** regular critical
  points. `Crit(S_H) = Γ_closure` away from `x = −u, −v`: the derived
  closure set is a **Morse–Bott critical manifold** of `S_H`.
* **O3.** On `Crit(S_H)` the transverse second derivative is
  `∇²_⊥ S_H = cos θ · |∇θ|² = sgn(D) · |u×v|²/D²`, so the index is `0`
  on the `θ = 0` (`D > 0`) component and `1` on the `θ = π` (`D < 0`)
  component.
* **O4.** Hence the one-dimensional transverse stationary-phase magnitude is
  `|∇²_⊥ S_H|^{−1/2} = |D| / |u×v|`, which is the **positive coarea density**
  of round 5 up to the common convention factor `2`.

O4 is the reason this functional and not another is the round's candidate:
its Morse–Bott saddle measure reproduces, with nothing tuned, the density
already derived by conditioning.

## 2. The falsifier that comes with the oracle

For a semiclassical aggregation `∫ e^{iκ S_H} …` over the two components,
with `S_H = −1` at `θ = 0` and `S_H = +1` at `θ = π`, and Maslov factors
`e^{±i(π/4)·sgn(κ ∇²_⊥ S_H)}` differing by one unit of index,

```
A_π / A_0  =  e^{2iκ} · e^{∓iπ/2}        (sign per orientation convention).
```

Three consequences are pre-registered:

* **F1.** `|A_π/A_0| = 1` for real `κ`. The magnitudes never distinguish the
  branches, so **the entire fork is a phase question**, and no refinement of
  the saddle *magnitude* can settle it. **[Superseded by note N5 — this is
  wrong as stated; only the local phase prefactor has unit modulus.]**

* **F2.** The ratio is real — i.e. the aggregation is either the positive
  count or the oriented sum, and not something complex that is neither — iff
  `4κ/π` is an **odd integer**. This statement is independent of the
  orientation convention, which shifts `κ` by `π/2`.
* **F3.** Given F2, the sign alternates with `κ` mod `π`. So even a real
  ratio does not choose the branch: `κ` mod `π` must itself come from
  somewhere. **Stationarity of `S_H` alone therefore cannot resolve the
  fork.** If BAM cannot derive `κ`, this round has located a *fourth*
  underived ingredient rather than removed one.

`κ` may not be fixed by requiring the answer. Choosing `κ` because it yields
`−1` is recorded as `NOT DERIVED`.

## 3. Five independent questions, five independent verdict fields

The round reports **five** verdicts. Collapsing them into one headline is
forbidden.

### A — Closed-history action

Does BAM's existing classical structure determine a single history action
whose angular part is `S_H` (or a functional with the same critical
manifold, index pattern and transverse measure), without inventing an
unrelated functional?

Permitted verdicts:
`HISTORY_ACTION_DERIVED_AND_SELECTS_ORIENTED_CURRENT`,
`HISTORY_ACTION_DERIVED_BUT_BRANCH_PHASE_UNDERDETERMINED`,
`HOLONOMY_TRACE_IS_A_STATIONARY_FUNCTIONAL_NOT_A_DERIVED_ACTION`.

**Rule.** Constructing any `F(cos θ)` after seeing which one works counts as
`NOT DERIVED`. To count as derived, the functional must follow from a
structure already in the repository, named by file and line, that was not
introduced for this purpose. A test asserts that the class of functionals
`F(cos θ)` with `F' ≠ 0` all share `Crit`, so agreement of critical sets is
explicitly *not* evidence of selection.

### B — Sector coefficients

Compute the **full symmetry group of the fixed-setting boundary problem**,
including detector boundary conditions, and its action on the four outcome
sectors. Discriminator: `G_phys` transitive on the four sectors **and**
preserving the induced measure ⟹ `r = 1`; two orbits ⟹ `r` free.

Pre-registered expectation, stated in advance so that finding it is not a
result: the obvious symmetries `++↔−−` and `+−↔−+` decompose the sectors as
`{++,−−} ⊔ {+−,−+}`, which leaves `r` free. Mouth exchange alone is
therefore **not** sufficient and testing only it does not discharge B.
The round must additionally test whether the four sectors are the fibres of
a free `Z₂ × Z₂` deck action by measure-preserving maps of the history
space, which would force `r = 1`; a negative answer there is a real result.

Permitted verdicts: `PHYSICAL_SYMMETRY_FORCES_EQUAL_SECTOR_MEASURE`,
`LIKE_UNLIKE_SECTOR_RATIO_REMAINS_FREE`.

### C — Classical readout

Identify an **actual** classical quantity arriving at a detector in BAM's
existing matter model — named by file and line — and derive the detector
response functional from its classical coupling. Abstract plausibility
arguments do not discharge C.

Control, fixed in advance: rescale the incoming field, `φ ↦ cφ`. If
`R[cφ] = c² R[φ]`, the response is an ordinary energy/intensity response and
does **not** justify frequencies linear in the oriented current. If the
object is a conserved topological/Poincaré-dual current with `R[cJ] = cR[J]`,
linearity has independent geometric meaning. Agreement with Born may not
choose between them.

Permitted verdicts: `CLASSICAL_DETECTOR_DERIVES_LINEAR_HISTORY_READOUT`,
`CLASSICAL_DETECTOR_RESPONDS_QUADRATICALLY`,
`NO_BAM_DETECTOR_COUPLING_CURRENTLY_DEFINES_THE_READOUT`.

### D — Compatibility with the existing radial action

`geometrodynamics/tangherlini/operator_audit.py:192` already defines
`closed_orbit_action(ω, ℓ) = 2∫√(ω²−V_ℓ) dr*`, described in its own
docstring as "the closure ledger's `∮p dq`". Ask whether one classical
phase-space construction yields both it and a mouth-holonomy term, ideally
`S_closed = ∮p_r dr* + S_frame + S_detector` from a single symplectic
structure. If the only available connection is to append `−cos(Ω/2)`
because it has the desired stationary set, the history action is
**independently postulated** and D fails. This is a compatibility test, not
permission to invent an action; its purpose is to stop the repository
growing two unrelated meanings of "action".

Permitted verdicts: `ONE_SYMPLECTIC_STRUCTURE_YIELDS_BOTH`,
`HISTORY_ACTION_INDEPENDENTLY_POSTULATED`.

### E — Mandatory causality gate: source readability

Round 5 established `ρ(x|a,b) ≠ ρ(x)`. If `x` is a physical classical source
variable and any ordinary local coupling can read an informative `F(x)`,
future settings are in principle observable at the source. The requirement is
strictly stronger than the detector no-signalling already shown:

```
P(O_S | a, b) = P(O_S)     for every source-local observable O_S available in BAM.
```

Enumerate the source-local observables and couplings that actually exist in
the repository and test each. A positive result requires a **dynamical**
reason for source inaccessibility; the word "hidden" does not discharge E.
If avoiding the signal requires declaring `x` gauge, that collides with
round 5's use of `x` as a physical source degree of freedom and must be
stated as an unresolved collision, not resolved by preference.

Permitted verdicts: `SOURCE_OBSERVABLES_OPERATIONALLY_NON_SIGNALLING`,
`SOURCE_READOUT_SIGNALS_FUTURE_SETTINGS`,
`X_PHYSICALITY_COLLIDES_WITH_NO_SIGNALLING`.

## 4. Overall bar, stated in advance

The headline
`CLASSICAL_BAM_DERIVES_THE_SINGLET_PROBABILITY_LAW`
may be printed **only** if all four of the following hold simultaneously:
one classical action selects the oriented branch (A), `r = 1` follows from
an actual symmetry and invariant measure (B), a BAM detector derives linear
event readout (C), and source observables remain operationally
non-signalling (E). Anything less retains the individual labels.

Expected verdict, stated in advance so that obtaining it is not a success:
the bar is **not** expected to be met, and F3 is the pre-registered reason.

## 5. Dependency ledger to be printed

Every input the round uses, classified `derived` / `chosen` / `imported` /
`open`, including at minimum: the geodesic-realignment ansatz, the itinerary,
the physicality of `x`, the sector prior, `κ`, the orientation convention,
and any detector coupling used in C.

## Correction notes (post-implementation)

Recorded here rather than applied silently. **No criterion changed**; each
note narrows or sharpens what may be claimed.

**N1 — E has a parity structure that was not anticipated.** The gate was
pre-registered as a single question. The measured answer splits: the
outcome-summed conditioned density is *exactly* antipodally even
(`ρ(−x|a,b) = ρ(x|a,b)`, residual `3.4e-20` on an antipodally paired grid), so
**every odd source-local observable is blind** — `E[x·m|a,b] = 0` and
`P(x·m>0|a,b) = 1/2` exactly, for every axis and setting pair. This is a
genuine partial protection and is recorded as a finding in its own right.
It does not change the verdict: every **even** observable signals
(`E[(x·m)²]` spread `0.0103`), non-coplanar settings give mutually singular
supports (total variation `1`, a one-shot signal), and by question C every
classical coupling in BAM is degree-2 homogeneous, hence even — so the
readable class is exactly the class BAM has.

**N2 — E's verdict label.** Both `SOURCE_READOUT_SIGNALS_FUTURE_SETTINGS` and
`X_PHYSICALITY_COLLIDES_WITH_NO_SIGNALLING` are defensible. The first is
reported, because it is what was measured; the collision is the
interpretation and is stated explicitly in `docs/history_action.md` rather
than used as the label.

**N3 — two implementation bugs found and fixed before the result was
recorded**, both of the same kind: a check that compared a quantity with
itself. (i) The first version of question D computed the outer leg as
`whole − inner`, making additivity true by construction, and the split point
lay outside the classically allowed region so the outer leg was zero; it now
integrates both legs independently over an interval split *inside* the allowed
region (`0.9924` and `0.7358`, additive to `2.0e-8`). (ii) The first version of
`fibre_action_is_weight_blind` recomputed the same base integrand twice; it now
measures that the Hopf action is vertical (`3.3e-16`) and re-runs round 6's
fibre-independence measurement (`2.5e-15`). Neither bug had reached a verdict.

**N4 — a consequence of C that was not pre-registered.** The quadratic readout
is not a *weaker* correlation but a **superquantum** one: `E = −2cos γ/(1+cos²γ)`
with exact `1/2` marginals and `S_max = 8√2/3 = 3.7712 > 2√2`. Recorded because
it inverts the natural reading of item 3: linearity is not a convenience that
happens to give the quantum answer, it is what keeps the model *at* the
Tsirelson bound rather than past it.

**N5 — F1 was wrong, and the corrected statement is stronger.** F1 treated the
local phase prefactor as the whole component amplitude. For a Morse–Bott
critical *manifold*, stationary phase carries an integral over each component,

```
A_j ~ e^{iκS_j} e^{iπσ_j/4} (2π/κ)^{1/2} M_j ,   M_j = ∫_{C_j} a/√|H_⊥| dσ,
```

so the true ratio is `A_π/A_0 = (M_π/M_0) e^{2iκ} e^{∓iπ/2}`. Since
`|H_⊥|^{−1/2} = |D|/|u×v|`, the masses `M_0` and `M_π` are the positive coarea
masses of the two closure arcs, and they are **generically unequal** — measured
`M_0 = 11.6708`, `M_π = 0.1695`, ratio `0.0145`. So `|A_π/A_0| ≠ 1`, and the
test `test_F1_the_saddle_magnitudes_never_separate_the_branches` was wrong and
is replaced.

What survives is cleaner than what was frozen. The two masses satisfy, to
residual `0.0`,

```
(M_0 − M_π)|u×v| = ∫_Γ D dσ      (the oriented sum)
(M_0 + M_π)|u×v| = ∫_Γ |D| dσ    (the positive count)
```

so stationary phase reproduces **both candidate aggregations exactly**, and the
sole undetermined element is their relative *phase*, `arg(A_π/A_0) = 2κ ∓ π/2`.
F2 and F3 are unaffected: the relative factor is real iff `4κ/π` is an odd
integer, and its sign alternates with `κ` mod `π`.

**N6 — C is downgraded to `NO_BAM_DETECTOR_COUPLING_CURRENTLY_DEFINES_THE_
READOUT`.** §3C required an actual detector coupling with its response
functional derived. That was not achieved: the five `source_audit.py` `T_kk`
are stress tensors, not detector models, and `mouth_flux` is a conserved flux,
not an event-counting coupling. The measured degree-2 homogeneity is real and
is retained as a supporting finding — *every existing matter/flux observable
examined is quadratic in field amplitude; none supplies a linear history-current
readout* — but it does not select a law. Two ordinary quadratic operations
disagree: squaring the integral gives `E = −2cos γ/(1+cos²γ)`,
`S_max = 8√2/3 = 3.7712`; squaring locally and then integrating gives, from the
exact `⟨D_s²⟩ = (1+c_s)(2+c_s)`, `E = −3cos γ/(2+cos²γ)`, `S_max = 3.3941`.
Both keep marginals at `1/2`. The earlier verdict
`CLASSICAL_DETECTOR_RESPONDS_QUADRATICALLY` claimed more than the code
establishes and is withdrawn.

**N7 — three overclaims in E, withdrawn.** (i) "Every even observable signals"
is false: constants and `x·x = 1` are even and blind, and both are now computed.
The supported claim is that *some* even functions separate the conditioned
ensembles. (ii) "BAM's couplings are degree 2 in `φ`, hence even in `x`" is a
non-sequitur — `x ↦ −x` on the mouth direction is not `φ ↦ −φ` on a scalar,
Maxwell or GL field, and no map from `x` to those field configurations exists in
the repository. Withdrawn. (iii) The observables used (`x·m`, `(x·m)²`, `|x·m|`)
are synthetic functions of `x`, not couplings BAM has, and nothing shows a BAM
apparatus could measure one without altering the global two-boundary solution —
in a two-boundary theory a source measurement is itself a further boundary
condition.

**N8 — E's label, and why it is outside the pre-registered set.** §3E offered
only `SOURCE_OBSERVABLES_OPERATIONALLY_NON_SIGNALLING`,
`SOURCE_READOUT_SIGNALS_FUTURE_SETTINGS` and
`X_PHYSICALITY_COLLIDES_WITH_NO_SIGNALLING`. That set has no term for what was
actually established, which is that the information is *present* while
readability is undemonstrated. The label used is
`SETTING_INFORMATION_IS_PRESENT_AT_SOURCE_READOUT_DYNAMICS_OPEN`. This is a
narrowing, not a weakening of the gate: it remains a serious causality hazard
requiring either a dynamical non-readability theorem or a reformulation in which
`x` is not operational.

**N9 — the A and D theorems are scoped.** "No single functional supplies both
the saddle magnitude and the branch phase" was too broad. What is established:
*no globally single-valued real class function of the holonomy alone is
simultaneously additive and has the closure manifold as its stationary set.* A
full classical action need not be a class function of `G` — it may contain
`∫p dq − H dt`, detector interactions, orientation-dependent potentials or
boundary terms. Likewise D establishes that the **repository's existing** radial
action plus the natural connection term does not generate `S_H`; it is not a
general no-go.

**N10 — the off-closure half of `Crit(S_H) = Γ` is now proved, not asserted.**
Checking `∇S_H = 0` *on* the closure set does not exclude critical points
elsewhere, and round 6's `∇θ|_{N=0} ≠ 0` concerns only `N = 0`. With `p = u+v`,
`q = u×v`, `A = 1+u·v`, `∇_{S²}θ = P_x(Dq − Np)/(D²+N²)`, so a regular critical
point needs `Dq − Np ∥ x`. As `p ⊥ q`, either `x` has a component along
`p̂ × q̂` — forcing `Dq − Np = 0`, hence `D = N = 0`, the excluded chart
singularity — or `x ∈ span(p,q)`, where the condition reduces to
`x_p = −|p|/A = −sec(γ/2)`, impossible for `0 < γ < π` since `|x_p| ≤ 1`.
Confirmed to `1.8e-15` against the closed form, with a global grid search
finding `min|∇θ| = 0.285` off closure.

**N11 — B's scope, and a narrowed label.** `sector_symmetry_group` enumerates
the `O(3)` **base isometries** fixing both setting axes, and
`fibre_action_is_weight_blind` closes the vertical `Spin(2)` extension. That is
not a classification of every discrete symmetry of the fixed-setting
boundary-value problem, so calling it "the full fixed-setting symmetry group"
was stronger than the code constructs. New `discrete_symmetry_extension` checks
the three further operations the model actually supplies — detector-label
exchange `A↔B`, history reversal of the derived loop (`D` is symmetric in
`u,w`, so every weight is invariant to `0.0`, and `N → −N` vanishes on `Γ`),
and the Pin deck (fibre-blind) — and all three preserve `s_A s_B`. The result
survives, but as an enumeration rather than a classification, so the verdict
label becomes `NO_IDENTIFIED_SYMMETRY_FORCES_EQUAL_SECTOR_MEASURE` in place of
`LIKE_UNLIKE_SECTOR_RATIO_REMAINS_FREE`.

**N12 — the Morse–Bott masses are conditional on an unrecorded measure.**
`morse_bott_component_masses` sets `a = 1` and integrates uniform arclength, so
`M_j = ∫|D|/|u×v| dσ` reproduces round 5 **conditional on the inherited
Haar/uniform source measure** — equivalently unit amplitude in this
finite-dimensional integral. That is a choice, and it is now in the dependency
ledger. New `amplitude_dependence` demonstrates it: `a = 1 + 0.5 x·u` takes
`M_0` from `11.6708` to `14.4736` and the oriented sum from `9.6780` to
`12.0975`. (An amplitude varying only along `u×v` is invisible, since
`x·(u×v) = 0` on `Γ`.) So "with nothing tuned" now applies to the *Hessian*,
not to the measure. Also recorded: this is stationary phase over the
two-parameter family of geodesic triangles parameterised by `x`, **not** over
BAM's full loop or history space, and that scope is stated wherever `S_H` is
called a history-action candidate.

**N13 — C's reconnaissance is narrowed in prose.** The code measures six
specific quantities — the five `source_audit.py` null stresses plus
`mouth_flux`. That supports *all audited stress and flux quantities are
degree-2 homogeneous; none supplies a linear detector readout*. It does **not**
support "every existing classical observable in the repository is quadratic in
field amplitude"; in particular the perfect-fluid row scales `u·k` rather than
a field amplitude, and unaudited field/configuration variables exist. The
narrow form is now used in `README.md`, `docs/history_action.md` and the audit.
The C verdict itself was already correct.

**N14 — the endpoint problem, closed by excision rather than downgraded.**
`S_H` is undefined at `x = −u, −v`, so the critical components are open arcs
ending at singular points of the phase; in local tangent coordinates there
`D ~ sin γ·ξ`, `N ~ sin γ·η`, so `S_H ~ −ξ/√(ξ²+η²)` is genuinely
direction-dependent, and textbook Morse–Bott stationary phase wants smoothness
near the critical manifold. New `excision_estimate` supplies the missing
argument: on `Γ`, `|D| ~ sin γ·s` in arclength from the puncture, so the mass
inside a disc of radius `ε` is `sin γ·ε²/|u×v|` — two-sided, the factor two
because the disc covers the arc on both sides, which an initial one-sided
version of this estimate got wrong by exactly `2`. It is `O(ε²)`, measured
within `0.3%` of the closed form at `ε = 0.2` and improving to `1e-5` at
`ε = 0.025`. The punctures carry no mass, so the quoted coefficients are the
genuine ones and need not be relabelled formal.

**N15 — exact closed form for the local-square maximum.** `S_max` for
`E = −3cos γ/(2+cos²γ)` is attained at `γ = π/4`, where `E = −3√2/5`, giving
exactly `12√2/5 = 3.394112549...` rather than the grid decimal.
