# Compton vertex resummation probe — plan

Follow-on to PR #34 (O(ε²) closure at `(ν₀, ν₁, ν₂, ξ) = (9/4, −4,
7/4, −1/2)`). PR #34 surfaced two structural patterns:

  - `ν₀ = (−3/2)² = γ²(O(ε))` — recursive across orders
  - `ξ = −1/2 = −A_φ(0)` — Hopf-linked

The recursive pattern suggests the vertex modification factor
`F(ε, θ)` resums to a **closed form** at all orders in ε. This probe
derives F analytically, tests it numerically, and identifies its
BAM-interpretable factorisation.

## Exact analytic F

From the photon-structure construction:

    f_BAM_baseline_norm = (1+c²) · (1+1/x)² / 8
    f_KN_norm           = x² · (x + 1/x − sin²θ) / 2

The vertex factor F satisfies `f_BAM_modified = f_BAM_baseline · F²`,
so `F² = f_KN / f_BAM_baseline`. Direct algebra:

    F²(x, c) = 4 · x³ · (x² + 1 − x·sin²θ) / [(1 + c²) · (1 + x)²]

with `x = ω'/ω = 1/(1 + ε·(1 − cos θ))`. This is a **closed-form
rational function of x and c**. The resummation exists.

## Tests

### T1. Exact F² verification

Compute the closed-form `F²(ε, θ)` numerically and compare to
`f_KN/f_BAM_baseline` at a fine (ε, θ) grid. Verify agreement to
machine precision.

### T2. Small-ε Taylor expansion match

Expand `F(x, c)` to O(ε²) and compare to the (ν₀, ν₁, ν₂, ξ)
coefficients from PR #34. Verify that the closed form reproduces
the previously-derived perturbative coefficients exactly.

### T3. Recursive pattern test

Compute the O(ε^n) coefficients of `F` for n = 1, 2, 3, 4 and check
whether the constant piece at each order follows `(−3/2)ⁿ` (the
geometric resummation conjecture).

### T4. Natural factorisations

Test whether F factors as a product of BAM-interpretable pieces:
  - Kinematic Padé factor: `(2x/(1+x))^p`
  - Polarisation ratio: `((x²+1−x·sin²θ)/(1+c²))^q`
  - Throat-transport factor: `±x^r`

Find the best two-factor product representation.

### T5. Hopf-connection link

The ξ = −1/2 = −A_φ(0) pattern suggests F contains a factor tied
to the Hopf-connection charge. Test whether one factor of the
decomposition matches `(1 − A_φ(0)·something)` or a Hopf-fibre
holonomy.

### T6. KN reproduction with closed-form F

Plug the closed-form F into the BAM amplitude and verify exact
reproduction of f_KN at multiple ε values up to ε ~ 1.

## Predictions

**EXACT_RESUMMATION** (most likely): F² has the predicted closed
form; T1, T2, T6 all pass at machine precision. The structural
patterns from PR #34 are confirmed as consequences of the
underlying closed form.

**RECURSIVE_PATTERN_CONFIRMED**: T3 shows the (−3/2)ⁿ pattern
holds at each order n ≥ 1.

**HOPF_FACTOR_IDENTIFIED**: T5 finds a natural Hopf-linked factor.

## What this leaves open

  - **BAM derivation of F from first principles.** Even with the
    closed form identified, deriving it from a BAM Lagrangian /
    action remains the deeper analytic task.
  - **The closed-form F is `4·x³·(...) / ((1+c²)·(1+x)²)`.** The
    `(1+c²)` denominator is the polarisation sum factor itself —
    suggesting F should be derived AS a modification of the
    polarisation sum, not as a separate amplitude factor.
  - **Cross-process generalisation.** Does the same F work for
    pair production γγ → e⁺e⁻? Open follow-on.

## Cross-references

- PR #25–#34 in the QFT-event-reinterpretation thread.
- `experiments/closure_ledger/compton_vertex_resummation_probe.py`
  — this probe.
