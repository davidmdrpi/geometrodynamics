# Hopf-fibre helicity transport derivation of Q(x, c) — research plan

Follow-on to PR #38 (F² = K(x)²·Q(x, c) decomposition) and PR #39
(K(x) = 2x/(1+x) derived from equal-action splitting). PR #39 closed
the kinematic / caustic factor; this probe targets the **polarization
factor**

    Q(x, c) = x² + x·(1 − x)²/(1 + c²)

and tests whether it can be derived from a concrete BAM Hopf-fibre
helicity transport model.

## The two-channel structure

`Q` splits cleanly as a sum of two non-negative pieces (PR #38 T5):

    Q = |A_pres|² + |A_flip|²
    A_pres = x                                  (helicity-preserving)
    A_flip = √x · (1 − x) / √(1 + c²)           (helicity-flipping)

Equivalently, `A = (A_pres, A_flip)` is a 2-component "throat-pinch
polarization spinor" living in the Hopf-fibre helicity basis, and
`Q = ⟨A|A⟩`.

## BAM-native interpretation under test

### 1. Helicity-preserving amplitude `A_pres = x`

By the equal-action splitting of PR #39 applied to energy flux:

  - per-mouth photon amplitude through the throat = `√x`
    (the per-mouth pinch amplitude that reproduces `K(x) = 2x/(1+x)`
    after harmonic-mean normalisation)
  - two mouths contribute coherently → `A_pres = √x · √x = x`

This is the "diagonal" channel: helicity is preserved through both
pinches, both mouths contribute the energy-flux amplitude `√x`.

### 2. Helicity-flipping amplitude `A_flip = √x · (1−x) / √(1+c²)`

Same equal-action principle applied to spin/helicity:

  - per-pinch helicity-flip amplitude = `(1−x)/√(1+c²)`, where:
    * `(1 − x)` is the **recoil deficit** — the photon energy
      transferred to the electron, available to drive an angular
      momentum kick
    * `1/√(1+c²)` is the **inverse Thomson polarization sum
      normalisation** — the angular distribution of the helicity-flip
      amplitude must be normalised by the Hopf-fibre helicity sum
      `(1+c²)/2 = cos⁴(θ/2) + sin⁴(θ/2)`
  - one mouth flips, the other preserves → `A_flip = √x · (1−x)/√(1+c²)`

The per-pinch flip amplitude is **non-trivial only at finite recoil**
(`A_flip → 0` as `x → 1`), consistent with the absence of helicity
mixing in the Thomson limit.

### 3. Hopf-fibre helicity sum (Wigner-d¹)

The normalising factor `(1+c²)/2` is the sum

    Σ_{λ ∈ {+1, −1}} |d¹_{1, λ}(θ)|²  =  cos⁴(θ/2) + sin⁴(θ/2)
                                       =  (1 + cos²θ) / 2

over Wigner-d¹ matrix elements for spin-1 helicity transport through
angle θ. PR #38 T4 verified this identity to machine precision and
identified it with the Hopf-fibre photon helicity sum. This is the
**rigorously derived** ingredient; the Q derivation builds on it.

## Tests

  T1. **Wigner-d¹ Thomson sum** (recap PR #38 T4): verify
      `cos⁴(θ/2) + sin⁴(θ/2) = (1 + cos²θ)/2`, identifying the
      `(1+c²)` factor in Q's denominator as the Hopf-fibre
      helicity-sum normalisation.

  T2. **Polarization spinor ansatz**: define `A = (A_pres, A_flip)`
      with `A_pres = x`, `A_flip = √x·(1−x)/√(1+c²)`. Verify
      `⟨A|A⟩ = A_pres² + A_flip² = Q(x, c)` to machine precision
      across an `(x, c)` grid.

  T3. **`A_pres = x` from per-mouth energy amplitude**: in PR #39
      the harmonic-mean throat-rate derives from per-mouth amplitude
      `√x` (each mouth carries half the closure-quantum action with
      energy-flux factor `√x`). Two mouths in the preserving channel
      contribute `√x · √x = x`.

  T4. **`A_flip = √x·(1−x)/√(1+c²)` from equal-spin-action splitting**:
      analogously, the per-mouth flip amplitude `(1−x)/√(1+c²)` is
      the recoil-deficit `(1−x)` weighted by the inverse Thomson
      polarization-sum normalisation `1/√(1+c²)`. One mouth preserves
      (amplitude `√x`) and the other flips (amplitude `(1−x)/√(1+c²)`),
      giving `A_flip = √x · (1−x)/√(1+c²)`.

  T5. **Alternative flip-amplitude weightings rejected**: candidate
      per-mouth flip amplitudes
        - `(1−x)` (no Thomson normalisation)
        - `(1−x)·(1+c²)` (Thomson sum in numerator instead of denom)
        - `(1−x²)/√(1+c²)` (energy-rescaled deficit)
        - `√(1−x)/√(1+c²)` (square-root deficit)
      All fail to reproduce `A_flip² = x(1−x)²/(1+c²)` exactly.

  T6. **Helicity-channel orthogonality at Thomson limit (`x = 1`)**:
      verify `A_flip → 0` and `A_pres → 1` as `x → 1`, so Q → 1 (no
      vertex modification at Thomson).

  T7. **Complex-amplitude / quadrature reading**: the spinor
      `A = (A_pres, A_flip)` is equivalent to a complex Hopf-fibre
      amplitude `A_complex = A_pres + i·A_flip`. The imaginary unit
      represents a 90° Hopf-fibre rotation = helicity-flip. Verify
      `|A_complex|² = Q`.

  T8. **K²·Q reproduces F² (cross-check with PR #38 and PR #39)**:
      `F²(x, c) = K(x)² · Q(x, c)` to machine precision.

  T9. **Cross-process consistency**: under crossing `x → x_⊗ < 0`,
      both `A_pres(x_⊗)` and `A_flip(x_⊗)` continue analytically;
      the spinor ansatz preserves the Q decomposition in the
      BW/annihilation region.

## Verdict structure

  - **Q_DERIVED**: T1–T8 all pass at machine precision; T5 rejects
    alternative flip-amplitude weightings; T9 confirms cross-process
    consistency. The Q polarization factor is derived from
    Hopf-fibre helicity transport (Wigner-d¹) + equal-action
    splitting of per-mouth amplitudes for energy and spin channels,
    with the helicity-flipping channel normalised by the inverse
    Thomson polarization sum.

  - **Q_PARTIAL**: Layer 1 (Thomson Wigner-d sum) verified, but
    Layer 2 (the specific A_flip form) is not uniquely picked out
    by the BAM principle — alternative weightings match too. The
    derivation is informative but not unique.

  - **Q_DERIVATION_BROKEN**: T2 fails, meaning the proposed spinor
    decomposition doesn't reproduce Q algebraically — would indicate
    a structural error in the ansatz.

## What this leaves open

  - **First-principles BAM action**: like PR #39, the equal-action
    splitting (for both energy and spin) is the natural flux-continuity
    postulate, but a derivation from a specific BAM S³ throat action
    coupled to the Hopf bundle remains the deeper task.

  - **Helicity-resolved Compton comparison**: the standard QED
    helicity-resolved Compton amplitudes `|M(λ → λ')|²` have their
    own (different) algebraic split. The BAM decomposition is one
    consistent split; relating it to the QED helicity basis is a
    follow-on.

  - **Loop corrections**: still tree-level.

## Cross-references

  - PR #38: `throat_nucleation_caustic_derivation_probe.py` — F² =
    K(x)²·Q(x, c) decomposition; Thomson pol sum from Wigner-d¹.
  - PR #39: `two_mouth_flux_action_probe.py` — K(x) = 2x/(1+x) from
    equal-action splitting at the two mouths.
  - `geometrodynamics/hopf/spinor.py` — SU(2) spinor transport (the
    geometric infrastructure for Hopf-fibre helicity).
  - `geometrodynamics/hopf/connection.py` — Hopf connection
    `A_φ(χ) = ½cos(χ)`.
  - `experiments/closure_ledger/hopf_helicity_transport_probe.py` —
    this probe.
