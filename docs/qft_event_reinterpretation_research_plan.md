# QFT-event reinterpretation via antipodal transactions — research plan

Opens a new dynamical thread orthogonal to the closure-ledger and
gauge-falsifier work. The closure-ledger thread (PRs #11–22, #23, #24)
established the *static* eigenvalue / Berry-phase content of BAM in
the lepton sector; the gauge-falsifier (PR #24) established that the
symmetric-gauge throat transport `T = iσ_y` is the standard spin-½
representation in a particularly clean section. This thread tests
whether BAM has *dynamical* content matching QFT event structure.

## The question

Can BAM's three composable dynamical elements

  - **throat worldlines** — each particle is a pair of mouths on `S³`
    connected by a throat with its own worldline trajectory;
  - **time dilation at mouth** — proper time at each mouth tracks the
    mouth's velocity in the bulk frame, generically different from
    coordinate time;
  - **antipodal closure** — the back mouth lives at the `S³` antipode
    of the front mouth and must close the transaction consistently,

reproduce QFT event structure for a simple local interaction?

The natural starting point is **Compton-like scattering**: a single
QED vertex pair `γ + e → γ + e` rewritten as a closed antipodal
transaction with mouth proper-time skew as the new dynamical variable.

## Why Compton first

  - **Two-vertex tree-level QED diagram** with a single Mandelstam
    invariant `s = (k + p)²` — minimal QFT structure.
  - **Low-energy limit (Thomson scattering)** has a clean analytic
    cross section `(r_e²/2)(1 + cos²θ)` that any reinterpretation must
    eventually reproduce.
  - **Crossing symmetry** (Compton `↔` pair annihilation under
    `s ↔ u`) is the natural arena for the antipodal map — antipodal
    closure on `S³` might *be* crossing symmetry in BAM language.
  - **Tree-level only** — no loop integrals to relate to BAM's bulk
    radial channel, which keeps the test purely kinematic at first.

Pair production, electron self-energy, vertex correction are larger
follow-on thread targets that need this kinematic foundation first.

## What "antipodal transaction" means concretely

A standard QED scattering event `(p, k) → (p', k')` is *local*: the
4-vertex sits at a single spacetime point `X`. BAM rewrites this as:

  - **Front-mouth vertex** at `X`: the front mouths of the electron's
    throat-pair and the photon's throat-pair (if photons admit a
    throat-pair description — see open questions) meet at `X`.
  - **Back-mouth vertex** at the `S³` antipode `X̃` of `X`: the back
    mouths meet at `X̃`. Antipodal in the spatial slice; same
    coordinate-time.
  - **Throat-worldline pinch**: between the in- and out-states, each
    throat traces out a worldsheet whose mid-line is the throat
    trajectory. The pinch at scattering deforms this worldsheet.
  - **Transaction closure**: the front and back vertices are linked
    by a Wheeler-Feynman handshake (retarded offer + advanced
    confirm + phase closure) on the closed `S³` slice.

The conservation law `p + k = p' + k'` holds at the front vertex; by
the antipodal map, an equivalent conservation `−p + −k = −p' + −k'`
holds at the back vertex automatically (it is the parity-flipped
front vertex). So conservation at the back vertex is not a new
constraint — it is forced by the closure on `S³`. This is the BAM
prediction: *no extra independent conservation laws beyond the standard
ones*. If the antipodal-closure picture imposed additional independent
constraints, those would either (a) over-determine standard Compton
kinematics (falsifying BAM) or (b) be new selection rules visible in
high-precision QED tests.

## Mouth proper-time skew — the new dynamical variable

The candidate BAM-specific dynamical content is the **proper-time
skew** Δτ between mouths. Several distinct quantities all carry this
label and need to be disambiguated:

  - **Inter-mouth skew Δτ_inter** between the two mouths of the same
    particle's throat-pair. If the front mouth moves at 3-velocity v
    in the bulk frame and the back mouth at antipodal point moves at
    `−v` (the antipodal map negates spatial vectors in `R⁴`), then
    both mouths have the same boost factor `γ_v` and **the inter-mouth
    skew vanishes identically**. This rules out one candidate
    interpretation immediately and is the first thing the probe
    should verify.
  - **Inter-particle skew Δτ_ext** between two distinct particles at
    a shared vertex. Standard relativistic kinematics: the electron's
    proper time differs from the photon's "proper time" (which is
    `0`) by the boost factor. This is not BAM-specific.
  - **Throat-pinch skew Δτ_throat** between the moment the front
    mouth experiences the scattering and the moment the back mouth
    "knows" about it through the throat. In a static throat geometry
    this is `0` (information travels instantly through the throat in
    the BAM picture); in a dynamic throat it could be non-zero. This
    is the most interesting candidate.
  - **Closure-cycle skew Δτ_cycle** integrated around a full
    antipodal cycle (front vertex → throat traverse → back vertex →
    throat back). For a closed transaction this must be a *natural*
    multiple of a closure quantum if BAM is to reproduce discrete
    QED amplitudes.

The first probe enumerates which of these are non-trivial and which
satisfy a clean closure quantum condition.

## Falsifiable predictions for the first probe

### P1. Antipodal closure does not over-constrain standard Compton kinematics

Given the in-state `(p_e, k_γ)` with `p_e + k_γ = P` (total
4-momentum), an out-state `(p_e', k_γ')` satisfying
`p_e' + k_γ' = P` exists on a 2-parameter family parametrised by the
scattering angles `(θ, φ)`. The BAM antipodal-closure constraint
(back vertex's 4-momenta are the spatial inversion of the front's)
must be satisfiable for *every* point on this family — otherwise the
closure imposes additional selection rules not present in standard
QED.

PASS: for every sampled `(θ, φ)`, the antipodal back-vertex
4-momenta are consistent (3-momentum conservation at the back vertex
holds identically by parity, energy conservation holds because energy
is invariant under spatial inversion).

FAIL: some `(θ, φ)` produces an inconsistent back vertex. This
would mean BAM forbids certain Compton scattering angles — a sharp
falsification.

### P2. Inter-mouth skew Δτ_inter vanishes identically

For the same-particle, both-mouths case, the bulk-frame velocity is
antipodal-flipped and the boost factor is identical, so
`Δτ_inter ≡ 0` exactly. The first probe verifies this for a sweep
of electron lab-frame velocities and confirms that the only
non-vanishing skews are inter-particle or throat-dynamic.

### P3. Throat-pinch skew Δτ_throat is the BAM-specific dynamical content

In the static-throat limit, `Δτ_throat = 0`. For Compton scattering,
the throat geometry deforms during the pinch — the two electron
mouths must "split" momentarily because the front mouth changes
4-velocity at the vertex while the back mouth cannot
instantaneously feel this (causality on `S³`). The probe should
compute `Δτ_throat` for a specific scattering geometry and identify
whether it equals `ℏ / E_CM` (would point to an "uncertainty-like"
reading) or a topological invariant (would point to a closure-quantum
reading).

### P4. Crossing symmetry as antipodal symmetry

The s-channel `(p_e + k_γ → p_e' + k_γ')` and the u-channel
`(p_e + k_γ' → p_e' + k_γ)` (with crossed photons) are related by
QED crossing symmetry. The antipodal map on the back vertex
`p ↔ −p` should map the s-channel front vertex to the u-channel
back vertex (or vice versa). If this map is exact, then crossing
symmetry in QED *is* antipodal closure in BAM. The probe samples
several `s, u` pairs and verifies the correspondence to numerical
precision.

## What the first probe is NOT

This is a kinematic feasibility test, not a cross-section derivation.

  - **No Klein-Nishina amplitude**. Reproducing the differential
    cross section from BAM's transaction amplitude algebra is a
    larger follow-on probe. The current probe checks only that BAM's
    kinematic skeleton is consistent with QED's tree-level kinematic
    skeleton.

  - **No loop corrections**. The bulk radial channel might encode
    one-loop self-energy / vertex corrections in some BAM
    reinterpretation. The current probe is at tree level only.

  - **No Lorentz-invariance test**. The antipodal map is defined in
    a specific bulk frame (the rest frame of `S³`). Whether the
    closure condition is Lorentz-covariant when transcribed to other
    bulk frames is a separate target.

## Stopping condition

The first probe closes when:

  (a) **P1–P4 all PASS.** Antipodal closure reproduces standard
      Compton kinematics without over- or under-constraining; the
      mouth proper-time skew zoo is disentangled and the
      BAM-specific element (throat-pinch skew) is identified. The
      thread then continues into the amplitude / cross-section
      probe.

  (b) **Any of P1–P4 FAILS.** Either the antipodal-closure picture
      is inconsistent with standard QED kinematics (over- or
      under-constrains), or one of the proper-time skews is
      ill-defined. The failure mode locates exactly where the BAM
      QFT-event reinterpretation breaks down, and the thread either
      pivots to a different scattering geometry or closes with a
      negative result.

## Cross-references

- `geometrodynamics/transaction/handshake.py` — existing
  Wheeler-Feynman handshake on `S³`. The amplitude algebra
  (offer × confirm × phase-closure) is the eventual substrate for
  the follow-on cross-section probe.
- `geometrodynamics/transaction/particles.py` — `Particle4` with
  `p4` (position on `S³`) and `vel4` (4-velocity); `MouthState`
  with throat modes. The probe's particle representation builds on
  these.
- `geometrodynamics/transaction/s3_geometry.py` — antipodal map
  `antipode4`, geodesic distance `geo4`, Green function. The
  probe uses `antipode4` as the spatial-inversion operator that
  defines the back vertex.
- `experiments/closure_ledger/compton_antipodal_kinematics_probe.py`
  — first probe in this thread.

## Cross-thread connections

The closure-ledger thread established that BAM's static eigenvalue
content is consistent with QM observables (Bell, spin-½ holonomy,
lepton ladder up to ℏ). This thread tests the *dynamical* content
that closure-ledger left untouched. A successful Compton probe would
extend BAM's claim from "particle = closed eigenstate of `S³`
geometry" to "particle interaction = closed transaction across the
`S³` antipode" — the natural composition.

The gauge-falsifier (PR #24) established that the symmetric-gauge
throat transport `T = iσ_y` is gauge-equivalent to the standard
Bloch construction, so the *transport map* at the throat is fixed.
The remaining freedom is in the throat's *worldline dynamics* — how
the throat trajectory moves through `S³` and how mouth proper times
accumulate. This is what the Compton probe targets.
