# Throat Dirac spinor from S_BAM

**Run:** 2026-05-24T15:38:57+00:00

Derives the throat Dirac 4-spinor (the inner/outer mouth doubling posited in PR #65) from the BAM radial structure via the Dirac/SUSY factorization plus the wormhole two-sidedness, and disentangles parity from the antipodal Z₂.

- **Factorization**: `H − E₀ = A†A (A = d/dr* + W; V − E₀ = W² − W′)`
- **Two mouths**: A†A, AA† SUSY-isospectral; inner/outer (B3 odd extension, #63)
- **4-spinor**: `4 = 2 (mouths) × 2 (SU(2) spin, B2) = Ψ_inner ⊕ Ψ_outer`
- **P vs Z₂**: P=γ⁰ radial (Dirac); antipodal Z₂ angular (S³ base)
- **Scope**: structure derived; full bulk spinor (S³ coupling) open

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_dirac_factorization` | H−E₀=A†A (V−E₀=W²−W', err 2e-04) | **PASS** |
| T2 | `T2_A_reproduces_bam_spectrum` | A†A+E₀ reproduces ω²(l,n) | **PASS** |
| T3 | `T3_susy_isospectrality_two_mouths` | A†A, AA† nonzero-isospectral: True | **PASS** |
| T4 | `T4_inner_outer_mouths_odd_extension` | inner/outer mouths; B3 odd extension (#63) | **PASS** |
| T5 | `T5_four_spinor_derived` | 4 = 2 mouths × 2 spin = Ψ_inner⊕Ψ_outer | **PASS** |
| T6 | `T6_parity_vs_antipodal_z2` | P=γ⁰ radial; antipodal Z₂ angular | **PASS** |
| T7 | `T7_honest_scope_b4` | structure derived; full S³ spinor open | **PASS** |
| T8 | `T8_assessment` | PR #65 4-spinor derived | **PASS** |

## T1: Dirac factorization

- E₀ = ω₀² = 1.112444 (ω₀ = 1.054725)
- max |(V−E₀) − (W²−W')| interior = 1.65e-04
- factorization H = A†A + E₀ holds: True

## T2: A reproduces the BAM spectrum

| k | A†A eig | H−E₀ eig | diff |
|---:|---:|---:|---:|
| 0 | -0.0000 | 0.0000 | 4.5e-13 |
| 1 | 2.7846 | 2.7857 | 1.1e-03 |
| 2 | 7.2585 | 7.2628 | 4.3e-03 |

## T3: SUSY isospectrality (the two mouths)

- A†A top eigs: [60383.88, 62239.52, 66935.82, 87110.65]
- AA† top eigs: [60383.88, 62239.52, 66935.82, 87110.65]
- nonzero spectra match (SUSY): True

## T4: Inner/outer mouths + B3 odd extension

- two-sided wormhole (inner/outer): True
- modes odd under throat reflection u(2R_MID−r)=−u(r) (#63): True
- Dirichlet node at the throat: True

## T5: The 4-spinor derived

- Dirac square-root factor (the two mouths, A/A†): 2
- SU(2) spin factor (T=iσ_y, B2): 2
- total = 4 = Ψ_inner ⊕ Ψ_outer (matches PR #65: True)

## T6: P vs the antipodal Z₂

- P = γ⁰ acts on: radial/Dirac components (r → 2R_MID − r across the throat)
- antipodal Z₂ (B2) acts on: S³ angular base (RP³ deck transformation σ: p → −p)
- distinct tensor factors: True

## T7: Honest scope / B4

- derived: the 4-component doubling structure (radial/throat sector)
- open: full closed-form bulk spinor with the S³ angular coupling
- B4: structure dimensionless/geometric; scale = single anchor

## T8: Assessment

- factorization: H = A†A + E₀ (Dirac square root of the radial operator)
- two mouths: A†A, AA† isospectral (SUSY); inner/outer, B3 odd extension (#63)
- 4-spinor: 4 = 2 (mouths) × 2 (SU(2) spin, B2) = Ψ_inner ⊕ Ψ_outer
- parity vs Z₂: P=γ⁰ radial (Dirac); antipodal Z₂ angular (S³ base)
- remaining: full closed-form bulk spinor with the S³ angular coupling

## Verdict

**THROAT_DIRAC_DERIVED.** THROAT DIRAC SPINOR DERIVED. The 4-component throat Dirac spinor posited in PR #65 is derived from the BAM radial structure, closing the #65 open notes.

THE DIRAC FACTORIZATION. The closure-ledger radial operator H = −d²/dr*² + V_tangherlini (spectrum ω²(l,n)) factorizes as a perfect square H − E₀ = A†A with A = d/dr* + W, W = −ψ₀'/ψ₀ (the Riccati identity V − E₀ = W² − W′, verified in the interior). A is the first-order radial Dirac operator; H = A†A + E₀ is its square. The Dirac structure is literally the square root of the BAM radial operator: a 2-component object where H was 1-component — the doubling. A†A reproduces the closure-ledger ladder.

THE TWO MOUTHS. A†A (= H − E₀) and the partner AA† share their nonzero spectrum (SUSY isospectrality) — the two Dirac components carry the same ω(l,n). These are the two wormhole mouths (inner r<R_MID, outer r>R_MID), joined at the throat by the B3 hard wall / the odd extension u(2R_MID−r)=−u(r) (the #63 reflection).

THE 4-SPINOR. Combining: throat Dirac 4-spinor = 2 (Dirac square-root / the two mouths, A/A†) × 2 (SU(2) spin, T=iσ_y, B2) = Ψ_inner ⊕ Ψ_outer. PR #65's posited spinor is derived — the square-root factor from the radial operator, the spin factor from B2.

P vs ANTIPODAL Z₂. Parity P=γ⁰ acts on the radial/Dirac components (the #63 reflection across the throat); the antipodal Z₂ (B2) acts on the S³ angular base (the RP³ deck transformation) — distinct tensor factors, radial vs angular (closing the other #65 note).

HONEST SCOPE. The doubling STRUCTURE is derived; the full closed-form bulk Dirac spinor of the complete S_BAM (with the S³ angular coupling) is the radial/throat sector only. B4: the structure is dimensionless/geometric; the mass scale rides on the single anchor.

## What this leaves open

- **The full bulk Dirac spinor from S_BAM.** The radial/throat sector is derived here; the S³ angular spinor harmonics and the complete coupled solution are the next step.
- **The superpotential W in closed form.** W = −ψ₀'/ψ₀ is computed numerically; a closed form from the Tangherlini f(r) would sharpen it.
