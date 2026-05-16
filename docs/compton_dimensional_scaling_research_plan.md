# Tangherlini dimensional-scaling discrimination probe — plan

Follow-on to PR #32, which identified 8 plausible BAM derivations
of `γ = −3/2` with no unique winner. PR #32 listed three concrete
discrimination paths; this probe executes the first one: **does
`γ` scale with bulk dimension d, or is it d-independent?**

## The discrimination

The 8 candidates from PR #32 split into two structural classes:

  - **d-independent** (predict `γ = −3/2` in any d):
    A (electron Casimir × 2 mouths), B (photon Casimir − Hopf),
    D (Hopf charge × photon multiplicity), F (closure winding +
    Hopf), G (Pauli trace), H (two-mouth × spin-½ Casimir).

  - **d-dependent** (predict `γ` rescaling with d):
    C (embedding-dim / pol-count: predicts `γ_d = −d_spatial/(d_spatial−1)`),
    E (natural units, speculative).

If γ depends on d, the d-independent candidates are falsified.
If γ is constant in d, candidate C is falsified.

## Setup — BAM in d_spatial dimensions

Generalising the photon-structure probe (PR #28):

  - Transverse photon polarisations: `d_spatial − 1` modes.
  - Polarisation completeness: `Σ_λ ε^λ_i ε^λ_j = δ_ij − k̂_i·k̂_j`.
  - Sum over polarisations:
    `Σ_{λ,λ'} |ε^λ(k)·ε^{λ'}(k')|² = (d_spatial − 2) + cos²θ`.
  - Initial polarisation average: divide by `(d_spatial − 1)`.

The BAM amplitude at O(ε):

    f_BAM_d_O(ε) = (1 − c)·((d−2) + c²) / (d−1)

The d-dim Klein-Nishina analog (natural generalisation, derived
in the probe code):

    f_KN_d(ε, θ)/f_KN_d(ε, 0)
      = (x²/(d−1)) · ((d−2)(x + 1/x) − sin²θ + (3−d))

    f_KN_d_O(ε) = −2·(1 − c)·((d−2) + c²) / (d−1)

Δ_required_d at O(ε):

    Δ_required_d = −3·(1 − c)·((d−2) + c²) / (d−1)

## The matching equation

For a Family B vertex modification `(ε·ε'*)·(1 + ε·μ)` with
`μ = β·sin²θ + γ·(1−c)`:

    f_BAM_d_modified_O(ε)
      = (1−c)·((d−2)+c²)/(d−1)  +  2μ·((d−2)+c²)/(d−1)

Setting equal to f_KN_d_O(ε):

    (1−c) + 2μ = −2·(1−c)
    μ = −3·(1−c) / 2
    →  β = 0,  γ = −3/2  (independent of d)

The factor `((d−2)+c²)/(d−1)` **cancels exactly**. The result is
that the matching vertex coupling is universally `γ = −3/2`,
regardless of `d`.

## Predictions

### T1. Analytic d-independence

Verify analytically (by re-deriving in the probe) that the matching
equation reduces to `μ = −3(1−c)/2` independently of d.

### T2. Numerical d-scaling

Compute `γ_BAM_d` numerically for `d ∈ {3, 4, 5, 6}` by extracting
the O(ε) coefficient from the d-dim BAM amplitude. Verify that all
give γ = −3/2 to numerical precision.

### T3. Candidate predictions table

For each of the 8 PR #32 candidates, tabulate the predicted γ(d) for
d ∈ {3, 4, 5, 6}. Mark which candidates are consistent with the
numerical observation `γ_d = −3/2`.

### T4. Falsified candidates

The probe falsifies any candidate predicting `γ_d ≠ −3/2`. The
specific candidates falsified are reported.

### T5. Surviving candidates

The probe identifies the candidates surviving the d-scaling test.
This is the discrimination output — a narrower set of plausible
BAM derivations.

## Expected outcome

**FALSIFIES_C** (most likely): γ is d-independent, candidate C
(embedding/pol) is falsified, Casimir-based candidates survive.

**FALSIFIES_CASIMIRS** (unexpected): γ varies with d, the Casimir-
based candidates are falsified, candidate C is consistent.

**INCONSISTENT** (probe construction issue): numerical γ(d) doesn't
match either prediction class — would indicate the d-dim BAM analog
isn't well-defined.

## What this leaves open even after the d-scaling discrimination

Even if the d-independence is confirmed and candidate C is
falsified, **6 candidates remain**. Further discrimination requires:

  - O(ε²) analytic extension (different candidates predict different
    next-order coefficients).
  - Cross-process test (pair production γγ → e⁺e⁻).
  - Polarised cross sections (different candidates may give
    different polarisation-dependent corrections).

## Cross-references

- PR #25–#32 in the QFT-event-reinterpretation thread.
- PR #31 — analytic vertex derivation establishing γ = −3/2 in d=4.
- PR #32 — coefficient-origin probe enumerating 8 candidate
  derivations.
- `geometrodynamics/tangherlini/` — Tangherlini 5D bulk machinery.
- `experiments/closure_ledger/compton_dimensional_scaling_probe.py`
  — this probe.
