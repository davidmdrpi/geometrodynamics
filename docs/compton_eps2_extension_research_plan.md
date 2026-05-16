# O(ε²) analytic extension probe — plan

Follow-on to PR #33 (d-scaling discrimination, falsifying candidate
C). The d-scaling test left 7 surviving candidates for γ = −3/2;
this probe attempts the second discrimination path: extending the
analytic derivation to O(ε²) and asking which surviving candidates
naturally predict the next-order coefficient.

## The O(ε²) residual after PR #31's vertex

The Taylor expansion at small ε:

    f_KN(ε, θ)        = (1+c²)/2  − ε·(1−c)(1+c²)        + ε²·(1−c)²·(4+3c²)/2  + O(ε³)
    f_BAM_baseline    = (1+c²)/2  + ε·(1−c)(1+c²)/2      + ε²·(1−c)²(1+c²)/8     + O(ε³)

After PR #31's vertex modification (`γ = −3/2`, ε-only):

    f_BAM_PR31_O(ε²) = −(1−c)²(1+c²)/4

The O(ε²) residual that an extended vertex must close:

    Δ_required_O(ε²)
       = f_KN_O(ε²) − f_BAM_PR31_O(ε²)
       = (1−c)²·(4+3c²)/2 + (1−c)²(1+c²)/4
       = (1−c)²·(9 + 7c²) / 4

This is the analytic target. The probe's central question: what BAM
vertex structure at O(ε²) produces exactly this angular form?

## Candidate O(ε²) vertex modifications

The natural extensions of the PR #30/31 vertex families:

### Family B' — extended angular modulation at O(ε²)

    V_B' = (ε·ε'*) · (1 + ε·μ₁ + ε²·μ₂)
    μ₁ = −3(1−c)/2     (PR #31 fixed)
    μ₂ = β'·sin²θ + γ'·(1−c) + δ'·(1−c)² + ζ'·cos θ + η'·cos θ·(1−c)

The O(ε²) contribution to f_BAM from μ₂ (via the polarisation sum
factor `(1+c²)`):

    2·μ₂·(1+c²)/2 = μ₂·(1+c²)

Setting equal to Δ_required = (1−c)²·(9+7c²)/4:

    μ₂·(1+c²) = (1−c)²·(9+7c²)/4

This requires μ₂ to be a polynomial in c whose product with (1+c²)
gives the polynomial `(1−c)²·(9+7c²)/4`. **The probe will derive
this exactly.**

### Family A' — ε·k coupling at O(ε)

    V_A = (ε·ε'*) + ε·α·(ε·k̂')(ε'*·k̂)

This contributes to the polarisation sum at O(ε²) through |V_A|²
expanded to quadratic order:

    Σ_pol |V_A|² = (1+c²) − 2εα·c·sin²θ + ε²·α²·sin⁴θ

The α² sin⁴θ contribution is at O(ε²) and provides a structure
NOT in the Family B (1+c²) polynomial basis. This could close the
gap that Family B alone cannot.

### Family C extended — per-channel kinematic weights

    M_x ∝ G_S3(ψ_x)^p

At O(ε²), this contributes corrections to the propagator factor.
The probe will analytically expand and check.

## Discrimination targets from PR #33

7 surviving candidates for γ = −3/2. Their natural O(ε²) extensions:

  - **A (doubled electron Casimir)**: extends to C_2² for j=1/2 →
    predict O(ε²) ∝ (3/4)² = 9/16, doubled = 9/8.
  - **B (photon Casimir − Hopf)**: extends to (C_2(j=1))² − Hopf² →
    (2)² − (1/2)² = 4 − 1/4 = 15/4.
  - **D (Hopf × photon multiplicity)**: extends to (Hopf·2j+1)² →
    (1/2·3)² = 9/4.
  - **F (closure + Hopf)**: extends to (n_closure + 1/2)² = (3/2)² = 9/4.
  - **G (Pauli trace)**: extends to Σ_i Tr(σ_i⁴) — but σ_i² = I so
    Tr(σ_i⁴) = Tr(I) = 2, sum over i gives 6. Same as O(ε).
    Suggests Pauli gives no O(ε²) correction → would FAIL closure.
  - **H (two-mouth × spin-½ Casimir)**: extends to n_mouth · C_2² =
    2·9/16 = 9/8.

Different candidates predict different O(ε²) coefficients. The probe
solves the matching equation for the O(ε²) vertex coefficient
analytically and compares to each candidate.

## Tests

### T1. O(ε²) analytic expansion verification

Numerically extract O(ε²) Taylor coefficients of f_KN and f_BAM_PR31
at small ε and compare to closed-form predictions.

### T2. Δ_required_O(ε²) angular decomposition

Show that Δ_required = (1−c)²·(9+7c²)/4 has specific coefficient
structure in {1, c, c², c³, c⁴} basis: (9, −18, 16, −14, 7)/4.

### T3. Family B' matching equation

Solve for μ₂ such that μ₂·(1+c²) = Δ_required·4 in polynomial
basis. Identify whether a polynomial μ₂ exists (rank-check the
system) and report.

### T4. Family A contribution at O(ε²)

Include the α² sin⁴θ term and re-solve. Check whether adding this
structure makes the system solvable.

### T5. Comparison to candidate O(ε²) predictions

For each surviving PR #32/33 candidate, compute its natural O(ε²)
extension and compare to the structural coefficients of μ₂. Identify
which candidates (if any) predict the right next-order coefficient.

## Expected outcome

**CLOSED_WITH_FAMILY_A** if adding Family A (α² sin⁴θ term) closes
O(ε²) with clean (α, μ₂) values. This would identify ε·k coupling
as the natural BAM extension at O(ε²) — even though it was useless
at O(ε) (per PR #30/31).

**PARTIAL_CLOSURE** if any combination reduces O(ε²) residual
significantly but doesn't close exactly.

**FAMILY_B_INSUFFICIENT** if the matching equation for Family B'
alone has no polynomial solution (predicted, per the rational
function with (1+c²) denominator).

**CANDIDATE_DISCRIMINATION** if a specific candidate from the
surviving 7 uniquely predicts the right O(ε²) structure.

## What this leaves open even with closure

  - **O(ε³) and beyond.** The pattern might continue, requiring
    additional vertex structure at each order.
  - **Cross-process consistency.** The same vertex coupling should
    work for pair production γγ → e⁺e⁻ — testable as the next
    discrimination probe.
  - **BAM-Lagrangian derivation.** Even with a clean coefficient
    structure identified, deriving it from a first-principles BAM
    action remains open.

## Cross-references

- PR #25–#33 in the QFT-event-reinterpretation thread.
- PR #31 — analytic vertex derivation, γ = −3/2 at O(ε).
- PR #32 — 8 candidate derivations.
- PR #33 — d-scaling falsifies candidate C; 7 candidates survive.
- `experiments/closure_ledger/compton_eps2_extension_probe.py` —
  this probe.
