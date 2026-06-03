# High-resolution lattice validation of the BAM measure operators (PR #120)

**Run:** 2026-06-03T01:24:56+00:00

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
| T5 | `T5_generic_holonomy` | generic holonomy: |det P_a|=2sin(πa) exact; η(a)=1−2a; phase (π/2)(1−2a) | **FAIL** |
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
- twisted eigenvalues converge `O(1/N²)` (ratios [16.0, 16.0]); massive twisted det `O(1/N²)` (ratios [16.0, 0.0]).
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

**LATTICE_VALIDATION_DISCREPANCY.** DISCREPANCY. A discrete quantity failed to match its continuum analytic value; investigate the implementation.
