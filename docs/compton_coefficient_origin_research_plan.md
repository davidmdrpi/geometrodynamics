# Coefficient-origin probe for γ = −3/2 — plan

Follow-on to PR #31 (analytic vertex derivation establishing
`(α, β, γ) = (0, 0, −3/2)` as the exact O(ε) closure for Compton).
The clean rational `γ = −3/2` is a strong signal of a natural
BAM-derivable coefficient; this probe asks **which BAM ingredient
actually produces −3/2**.

## The number to derive

The vertex coupling at O(ω/m):

    V = (ε · ε'*) · [1 + γ · (ω/m) · (1 − cos θ)],    γ = −3/2

The dimensionful piece `(ω/m)·(1−cos θ) = (k·k')/m²` is the standard
QED kinematic invariant (Mandelstam-t style). The dimensionless
`γ = −3/2` is the new structural target. Sign and magnitude both
matter.

## Candidate BAM derivations

The probe tests a curated list of natural BAM ingredients that
plausibly produce `|γ| = 3/2`:

### (A) Spin-1/2 quadratic Casimir, doubled

    γ_A = −2 · C_2(SU(2), j=1/2)
        = −2 · (1/2)·(3/2)
        = −3/2  ✓

Source: the electron is spin-1/2 in BAM (Hopf-fibre holonomy at the
poles). The Casimir `j(j+1) = 3/4`, doubled to `3/2` by the
double-mouth structure (two throat traversals per scattering
event).

### (B) Photon spin-1 Casimir minus 1/2

    γ_B = −(C_2(SU(2), j=1) − 1/2)
        = −(1·2 − 1/2)
        = −3/2  ✓

Source: the photon polarisation sum gives spin-1 structure; the
−1/2 subtraction is the Hopf-connection factor at χ = 0.

### (C) S² embedding dimension over polarisation count

    γ_C = −dim(R³ ambient of S²) / dim(transverse photon polarisations)
        = −3 / 2
        = −3/2  ✓

Source: the photon's polarisation lives in the tangent bundle of
S²; the ambient R³ has dim 3, the transverse subspace has dim 2.

### (D) Hopf-connection charge times spin-1 multiplicity

    γ_D = −A_φ(χ=0) · (2j_γ + 1)|_{j=1}
        = −(1/2) · 3
        = −3/2  ✓

Source: Hopf charge ½ at the BAM lock (χ=0) coupling to the
photon's (2j+1) = 3 spin states.

### (E) Throat-pinch O(ε) coefficient from antipodal traverse

The kinematic probe (PR #25) identified Δτ_throat = (π·R_S3/c)·(ε)
as the natural BAM-specific dynamical quantity. If `γ = −Δτ_throat/Δτ_natural`
where `Δτ_natural = (π·R_S3/c)/( something giving 3/2)`...

    γ_E = ? — speculative; tested for completeness.

### (F) Sum of Hopf-fibre windings

    γ_F = −(1 + 1/2) = −3/2  ✓

Source: 1 from antipodal closure on S³ (great-circle holonomy) +
1/2 from Hopf-fibre half-charge.

### (G) Trace of Pauli matrices squared

    γ_G = −(1/2) · Tr(σ_i σ_i)/2  for i = x, y, z
        = −(1/2) · (3·2)/2 = −3/2  ✓

Source: trace of Pauli generators of SU(2). The (1/2) comes from
the natural normalisation; the 3 from the three generators; the
factor 2 is the Pauli matrix dimension.

## Discriminating criteria

Multiple candidates evaluate to −3/2. The probe ranks them by:

  1. **Single-ingredient parsimony.** A derivation from one
     BAM-native object (e.g. just the Hopf charge) is more
     compelling than one combining several.
  2. **Sign consistency.** The sign should fall out of the
     derivation, not be added by hand.
  3. **Generalisability.** A derivation that predicts coefficients
     for *other* vertex factors (sin²θ, ε·k contractions) and
     correctly gives the zeros for those is the strongest. (PR #31:
     the analytic prediction is α = 0, β = 0 for these.)
  4. **Connection to existing BAM derivations.** Coefficients that
     already appear elsewhere in BAM (e.g. the Hopf holonomy
     `π·cos(χ)`, the closure-quantum `2π`, etc.) get higher rank.

## Tests

### T1. Numerical evaluation of each candidate

Compute each γ_X and verify it equals −3/2 to numerical precision.

### T2. Predict α and β simultaneously

For each derivation, check whether the same ingredient
**predicts α = 0 and β = 0 in addition to γ = −3/2**. A
derivation that only gets γ right but doesn't constrain α and β is
weaker.

### T3. Naturalness ranking

Rank candidates by the criteria above. Identify the cleanest
single-ingredient natural derivation.

### T4. Connection to existing BAM thread

Cross-reference each candidate's ingredients with the closure-ledger
results (PRs #11-#22). Coefficients already derived in those
threads get a connection score.

### T5. Verification by next-order prediction

The leading-O(ε) coefficient `γ = −3/2` is known. If the candidate
derivation extends naturally to O(ε²), use it to predict the next
coefficient and verify against the residual structure from PR #31
T5 (where the analytic optimum had residual ∝ ε^1.88).

## Expected outcome

**Best case (CLEAN_SINGLE_INGREDIENT)**: One candidate is
significantly more natural than others, predicts α = β = 0
simultaneously, and its ingredient appears elsewhere in BAM. The
−3/2 has a clear single source.

**Likely case (MULTIPLE_PLAUSIBLE)**: Several candidates equal
−3/2 but no clear discriminator. The probe identifies the leading
candidates and the discriminating experiments/derivations needed.

**Pessimistic case (NUMEROLOGY)**: All candidates are post-hoc
constructions giving 3/2 by coincidence. The probe rejects them
all and identifies what kind of analytic input would distinguish.

## Cross-references

- PR #25–#31 in the QFT-event-reinterpretation thread.
- `geometrodynamics/hopf/connection.py` — Hopf connection
  `A_φ = ½cos(χ)`.
- `geometrodynamics/embedding/transport.py` — `T = iσ_y`,
  `T² = −I`.
- `experiments/closure_ledger/compton_coefficient_origin_probe.py`
  — this probe.
