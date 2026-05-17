# Breit–Wheeler cross-process validation — research plan

Follow-on to PR #35 (Compton vertex resummation: exact closed-form
`F²(x, c) = 4·x³·(x² + 1 − x·sin²θ) / [(1+c²)·(1+x)²]` with
`x = ω'/ω`). PR #35 closed the Compton thread analytically. The
open follow-on, flagged at the bottom of the resummation plan, is:

> Cross-process generalisation. Does the same F work for pair
> production γγ → e⁺e⁻?

This probe answers that question.

## The question

The Compton BAM construction has two ingredients:

  - a **baseline** factor
    `f_BAM_baseline_C(x, c) = (1 + c²) · (1 + 1/x)² / 8`
    with `c = cos θ_lab`, `x = ω'/ω` — built from the Thomson
    angular factor `(1+c²)` and a lab-frame propagator-like factor
    `(1+1/x)²`.
  - a **vertex modification** `F²_C(x, c)` that, when multiplied
    into the baseline, exactly reproduces Klein–Nishina:
    `f_BAM_baseline · F²_C = f_KN_normalized` to machine precision.

Both ingredients were derived in the **Compton lab frame** with the
electron at rest. The question is whether the construction is
**process-general** — i.e. whether the same algebraic decomposition,
analytically continued to a different kinematic region, reproduces
the next QED tree process.

The minimal test is **Breit–Wheeler pair production** `γγ → e⁺ + e⁻`.
It is the cleanest crossed partner of Compton: the same two-vertex
tree diagram with photons in (rather than one in, one out) and a
lepton pair out (rather than one in, one out).

## Crossing relations

In Lorentz-invariant form, Compton `(γ(k) + e(p) → γ(k') + e(p'))`
uses

  - `s_C = (k + p)²`
  - `t_C = (k − k')²`
  - `u_C = (k − p')²`

For Breit–Wheeler `(γ(k₁) + γ(k₂) → e⁻(p₋) + e⁺(p₊))` we use

  - `s_BW = (k₁ + k₂)²`
  - `t_BW = (k₁ − p₋)²`
  - `u_BW = (k₁ − p₊)²`

The crossing that maps Compton to BW is

    k → k₁,    p → −p₊,    k' → −k₂,    p' → p₋

which yields

    s_C → u_BW,    u_C → t_BW,    t_C → s_BW.

Under this map, the Compton invariants used by the BAM construction
become

  - `x_C = (m² − u_C)/(s_C − m²)  →  x_⊗ = (m² − t_BW)/(u_BW − m²)`
  - `cos θ_lab,C = 1 + 2·t_C·m²/[(s_C − m²)·(m² − u_C)]
                →  c_⊗ = 1 + 2·s_BW·m²/[(u_BW − m²)·(m² − t_BW)]`

In the CM frame of BW with photon energy `E`, lepton velocity
`β = √(1 − m²/E²)`, and lepton CM angle `θ`:

    x_⊗ = −(1 − β·cosθ)/(1 + β·cosθ)               (negative)
    c_⊗ = (2β² − β²·cos²θ − 1)/(1 − β²·cos²θ)       (can be negative)
    sin²θ_⊗ = 4·β²·(1 − β²)·sin²θ/(1 − β²·cos²θ)²    (positive)

`x_⊗ < 0` is the *analytic continuation* of `ω'/ω` to the BW region.
`c_⊗` is the *crossed* cos-angle. These quantities are the BAM
construction's variables evaluated at BW kinematics.

## The standard BW prediction

The textbook unpolarised spin-summed squared amplitude for BW in
the CM frame is

    |M̄|²_BW / (8 e⁴)
        = 2·(1 + β²·c²)/(1 − β²·c²)
          + 4·(1 − β²)/(1 − β²·c²)
          − 4·(1 − β²)²/(1 − β²·c²)²

with `c = cosθ_CM`. The corresponding KN amplitude is

    |M̄|²_KN / (8 e⁴)
        = (s−m²)/(m²−u) + (m²−u)/(s−m²) + 2A + A²

with `A = 2t·m²/[(s−m²)(m²−u)]`. Standard QED crossing gives

    |M̄|²_BW(s, t, u) = −|M̄|²_KN evaluated at (s_C=u_BW, t_C=s_BW, u_C=t_BW)

— a **sign flip** from the single crossed fermion line. The
overall positivity of `|M̄|²_BW` is restored because the crossed
Compton expression evaluates to a negative number in the BW
physical region.

## Predictions

### P1. Invariant Compton amplitude crossed to BW kinematics reproduces |M̄|²_BW

Evaluate the Compton-side invariant `|M̄|²_KN(s, t, u)/(8e⁴)` at
the BW Mandelstam invariants `(s_BW, t_BW, u_BW)` and apply the
fermion-crossing sign flip. Compare to the direct textbook
`|M̄|²_BW(β, c)/(8e⁴)`. PASS at machine precision is the standard
QED prediction.

This is a consistency check on the Mandelstam algebra (not a BAM
test per se).

### P2. BAM baseline × F² crossed to BW equals the same |M̄|²_BW

The BAM Compton product `f_baseline_C · F²_C` equals
`f_KN_normalized = x·(x²+1−x·sin²θ)/2`. Express this in invariant
Mandelstam form and evaluate at BW kinematics. PASS if the BAM
crossed product matches the standard BW amplitude to the same
precision as P1.

If P2 PASSES: the BAM Compton construction commutes with crossing
— it is process-general at the amplitude level, and the same
closed-form F factor works for BW under analytic continuation.

If P2 FAILS: there is a Compton-specific algebraic feature in the
baseline or F that is not preserved by crossing. The construction
would then be a clever lab-frame fit rather than a process-general
amplitude.

### P3. BAM baseline and F² individually retain a sensible
       BW-side meaning

Even if the product P2 passes (as it should by construction if P1
passes), the *individual factors* `f_baseline_C(x_⊗, c_⊗)` and
`F²_C(x_⊗, c_⊗)` evaluated at BW kinematics may not have a
natural BW interpretation. Test:

  - Is `f_baseline_C` at BW kinematics positive and finite across
    the physical BW region (`0 < β < 1`, `−1 ≤ cosθ ≤ 1`)?
  - Does `F²_C` at BW kinematics retain the same monotonic
    structure as in the Compton region?
  - Is there a natural "BW baseline" that better matches the
    structural decomposition (e.g. a `(1 + β²c²)` Thomson-pair
    factor for the low-velocity BW limit)?

This is informative, not pass/fail.

### P4. BW total cross section recovery

Integrate the BAM-crossed differential cross section over solid
angle and compare to the standard Breit–Wheeler total
cross section

    σ_BW(s) = (π·r_e²/2)·(1 − β²)·[(3 − β⁴)·log((1+β)/(1−β))
                                   − 2β·(2 − β²)]

PASS at machine precision is the QED prediction. This is an
integrated end-to-end test of P2 across a wide `β` range.

### P5. Threshold and ultra-relativistic limits

  - **Threshold** (`β → 0`, `s → 4m²`): the BW total cross section
    vanishes linearly as `σ_BW ∝ β`. The BAM-crossed prediction
    must reproduce this linear suppression.
  - **Ultra-relativistic** (`β → 1`, `s ≫ 4m²`): the BW total
    cross section falls as `σ_BW ∼ (πr_e²·m²/s)·log(s/m²)`. The
    BAM-crossed prediction must reproduce this logarithmic falloff.

## Verdict structure

  - **PROCESS_GENERAL_UNDER_CROSSING**: P1, P2, P4, P5 all pass at
    machine precision. The BAM F is not a Compton-specific
    algebraic fit; it is the closed form of the invariant QED
    amplitude under analytic continuation. The same F (analytically
    continued) reproduces Breit–Wheeler pair production.

  - **CROSSING_AMBIGUOUS**: P1 passes (Mandelstam algebra works)
    but the BAM product P2 has a structural defect (e.g.
    f_baseline_C develops a pole or sign issue in the BW region
    that requires a compensating "regulator" in F²). The
    construction is process-general only after a frame-dependent
    fix-up.

  - **COMPTON_SPECIFIC**: P2 fails. The BAM construction has an
    intrinsic Compton lab-frame structure that cannot be analytically
    continued to BW. F would then be a Compton-only algebraic fit.

## What this leaves open

  - **Tree-level only.** Pair annihilation (`e⁺e⁻ → γγ`) and
    Bhabha / Møller scattering are obvious follow-ons but require
    new diagram structures (s-channel and t-channel together).
  - **Loop level.** One-loop corrections to the vertex factor are
    not addressed; the closed-form F is a tree-level reorganisation.
  - **BAM derivation of F from first principles.** Still open from
    PR #35.

## Cross-references

  - PR #35: `docs/compton_vertex_resummation_research_plan.md`,
    `experiments/closure_ledger/compton_vertex_resummation_probe.py`
    — the Compton thread culmination.
  - `experiments/closure_ledger/breit_wheeler_cross_process_probe.py`
    — this probe.
