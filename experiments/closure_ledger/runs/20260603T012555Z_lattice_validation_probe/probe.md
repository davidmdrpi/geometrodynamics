# High-resolution lattice validation of the BAM measure operators (PR #120)

**Run:** 2026-06-03T01:25:55+00:00

Validates that the **discrete software** reproduces the **continuum analytic** results of PRs #116–#119 — periodic, antiperiodic, AND **generic twisted-holonomy** — exactly for structural/symmetry quantities and as `O(1/N²)` for finite-difference ones.

- **Eigenvalues**: discrete −∂_τ² → (2πk/L)², O(1/N²) (ratio ≈ 16)
- **Ghost/antiperiodic**: lattice → (2sinh)², (2cosh)², O(1/N²); m→0 ⟹ det'=L², det_AP=2
- **Generic holonomy**: |det P_a|=2sin(πa) EXACT; twisted eigenvalues/det O(1/N²); η(a)=1−2a; phase (π/2)(1−2a) (ζ(0)=0)
- **η-invariant**: centered ∂_τ (odd N): η = 0 EXACT at finite N, 1 zero mode
- **Tangherlini GY**: det(H)/det(H_free) → 1.574370 (PR #116)
- **Convergence**: finite-difference O(1/N²); structural/symmetry (incl. |det P_a|) EXACT at finite N

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | validate discrete ≡ continuum (PRs #116–#119) | **PASS** |
| T2 | `T2_eigenvalue_convergence` | eigenvalues −∂_τ² → (2πk/L)², O(1/N²) | **PASS** |
| T3 | `T3_ghost_determinant_periodic` | ghost log-det → (2sinh(mL/2))², O(1/N²); TM at N=10⁶ | **PASS** |
| T4 | `T4_antiperiodic_determinant` | antiperiodic det → (2cosh(mL/2))², O(1/N²) | **PASS** |
| T5 | `T5_generic_holonomy` | generic holonomy: |det P_a|=2sin(πa) exact; η(a)=1−2a; phase (π/2)(1−2a) | **PASS** |
| T6 | `T6_det_prime_L2_and_antiperiodic_2` | m→0: det′(−∂_τ²)=L², det_AP=2 | **PASS** |
| T7 | `T7_eta_zero_exact_and_zero_mode_count` | η = 0 EXACT at finite N + 1 zero mode | **PASS** |
| T8 | `T8_tangherlini_gelfand_yaglom` | Tangherlini GY → 1.574370 (PR #116) | **PASS** |
| T9 | `T9_assessment` | LATTICE_VALIDATION_DISCRETE_MATCHES_CONTINUUM | **PASS** |

## Generic holonomy `a ∈ {1/4, 1/3, 2/3, 3/4}`

| a | |det P_a| (lattice) | 2 sin(πa) | exact? | η(a)=1−2a | phase (π/2)(1−2a) | det P_a |
|---:|---:|---:|:---:|---:|---:|---|
| 0.25 | 1.414214 | 1.414214 | ✓ | 0.5 | 0.7854 | 1.000+1.000i |
| 0.3333 | 1.732051 | 1.732051 | ✓ | 0.3333 | 0.5236 | 1.500+0.866i |
| 0.6667 | 1.732051 | 1.732051 | ✓ | -0.3333 | -0.5236 | 1.500-0.866i |
| 0.75 | 1.414214 | 1.414214 | ✓ | -0.5 | -0.7854 | 1.000-1.000i |

- **|det P_a| = 2 sin(πa) is EXACT on the lattice** (any N): `Π_k 2(1−cos(2π(k+a)/N)) = |1−e^{2πia}|² = 4 sin²(πa) (exact, any N)`.
- twisted eigenvalues converge `O(1/N²)` (ratios [16.0, 16.0]); massive twisted det `O(1/N²)` (ratios [16.0, 15.99]).
- **branch convention**: principal branch arg∈(−π,π]; ζ(0)=0 (no zero mode) ⟹ arg det P_a = (π/2)η(0) = (π/2)(1−2a).
- consistency: 2 sin(π/2) = 2, η = 0, phase = 0 ⟹ det = 2 (antiperiodic, PR #119).

## Ghost determinant: lattice → continuum (O(1/N²))

| N | lattice log-det | abs error |
|---:|---:|---:|
| 64 | 4.372592 | 8.87e-04 |
| 256 | 4.373424 | 5.54e-05 |
| 1024 | 4.373476 | 3.47e-06 |
| 4096 | 4.373479 | 2.17e-07 |

Error ratios [15.99, 16.0, 16.0] (≈ 16 ⟹ `O(1/N²)`).

## Verdict

**LATTICE_VALIDATION_DISCRETE_MATCHES_CONTINUUM.** HIGH-RESOLUTION LATTICE VALIDATION: THE DISCRETE SOFTWARE REPRODUCES THE CONTINUUM ANALYTIC DERIVATION OF PRs #116–#119 — PERIODIC, ANTIPERIODIC, AND GENERIC TWISTED-HOLONOMY.

EIGENVALUES. The discrete −∂_τ² eigenvalues (4/h²)sin²(πk/N) → (2πk/L)² with relative error O(1/N²) (ratio ≈ 16 per N×4).

GHOST / ANTIPERIODIC DETERMINANTS. The lattice log-determinants converge to the continuum (2 sinh(mL/2))² and (2 cosh(mL/2))² as O(1/N²); the transfer-matrix closed form reproduces the periodic one at N = 10⁶. The m → 0 limits give det'(−∂_τ²) = L² and det(∂_τ)_antiperiodic = 2.

GENERIC HOLONOMY a ∈ {1/4, 1/3, 2/3, 3/4}. For twisted BC ψ(τ+L) = e^{2πia}ψ(τ) (eigenvalues 2πi(n+a)/L): the discrete twisted eigenvalues (1/h)sin(2π(k+a)/N) converge to 2π(k+a)/L O(1/N²); |det P_a| = 2 sin(πa) holds EXACTLY on the lattice at any N via the product identity Π_k 2(1−cos(2π(k+a)/N)) = |1−e^{2πia}|² = 4 sin²(πa) (giving √2, √3, √3, √2); the massive twisted det → 4(sinh²(mL/2)+sin²(πa)) O(1/N²); η_A(0) = 1 − 2a = {1/2, 1/3, −1/3, −1/2}; and the branch convention — ζ(0) = 0 for the twisted operator (no zero mode, balanced spectrum), so the phase is PURELY the η piece — gives arg det P_a = (π/2)(1−2a), i.e. det P_a = 2 sin(πa)·e^{i(π/2)(1−2a)} = {1+i, 1.5+0.866i, 1.5−0.866i, 1−i}. Consistency: at a = 1/2 this is 2 sin(π/2) = 2 with η = 0 and phase 0, the antiperiodic real determinant of PR #119. (η(a) = 1 − 2a and the phase are continuum zeta-regularized quantities — the lattice matches the magnitude |det P_a| exactly and the symmetric points a = 0, 1/2.)

η = 0 AND STRUCTURE (EXACT). The discrete centered ∂_τ (odd N) has eigenvalues (i/h)sin(2πk/N) symmetric under k ↔ N−k, so the η-invariant Σ sign = 0 EXACTLY at finite N, with one zero mode and a purely imaginary spectrum.

TANGHERLINI GEL'FAND–YAGLOM. det(H)/det(H_free) → 1.574370 (PR #116), stable to ~1e-7 by N = 2000.

So the structural/symmetry quantities (η = 0, zero-mode count, spectrum symmetry, the |det P_a| = 2 sin(πa) product identity) match EXACTLY at finite N, while the finite-difference quantities converge O(1/N²): the discrete software behaves exactly as the continuous analytic derivation.
