# Boundary stress derivation of `χ_n` and singlet constraint

**Run:** 2026-05-27T07:37:00+00:00

Derives the generation-dependent Z₂ partition splitter `χ_n` from boundary stress at the inner/outer cavity mouths (`r = R_MID` throat-side, `r = R_OUTER` cavity-edge), and adds the singlet (color-neutral) projection constraint.

- **Identification**: χ_n = T_odd(n) from cavity-mouth boundary stress; uniform-sign / shell-suppressed; structurally insufficient for observed splittings; PR #80 color sector required
- **Next PR**: PR #80 — color algebra acting on (l, n, p) and H_couple
- **B4 caveat**: T_inner, T_outer dimensionful (1/length²); χ_n/ω² ratio is scale-free; structural

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_boundary_stress_computed_l1` | T_inner, T_outer computed for l=1, n=0..5; positive | **PASS** |
| T2 | `T2_z2_decomposition_inner_outer_swap` | Z₂ decomp: T_odd > 0 uniformly; asymmetry decreases with n | **PASS** |
| T3 | `T3_chi_n_from_boundary_stress` | χ_n = T_odd(n); positive across all n; shell-suppressed | **PASS** |
| T4 | `T4_sign_pattern_vs_pr78_ansatz` | NO sign flip — PR #78 ansatz structurally overruled | **PASS** |
| T5 | `T5_chi_n_magnitude_audit_vs_observed` | χ_n/ω² ~ 0.01–0.02; required ~0.65–0.99; 50–100× too small | **PASS** |
| T6 | `T6_singlet_projector_placeholder` | singlet projector = identity on flavor basis (placeholder) | **PASS** |
| T7 | `T7_structural_interpretation_pr79_vs_pr80` | PR #79 fixes χ_n slot; PR #80 owes color sector | **PASS** |
| T8 | `T8_honest_scope_assessment` | CHI_N_DERIVED_BOUNDARY_STRESS_INSUFFICIENT | **PASS** |

## T1–T2: Boundary stress and Z₂ decomposition

| n | ω | ω² | T_inner | T_outer | T_even | T_odd | asym = T_odd/T_even |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1.0547 | 1.1124 | 6.6988e-01 | 3.5787e-01 | 5.1388e-01 | 1.5600e-01 | 0.3036 |
| 1 | 1.9744 | 3.8982 | 2.2648e+00 | 1.8641e+00 | 2.0644e+00 | 2.0031e-01 | 0.0970 |
| 2 | 2.8941 | 8.3755 | 4.9365e+00 | 4.5662e+00 | 4.7513e+00 | 1.8514e-01 | 0.0390 |
| 3 | 3.8246 | 14.6278 | 8.6964e+00 | 8.3417e+00 | 8.5191e+00 | 1.7734e-01 | 0.0208 |
| 4 | 4.7609 | 22.6659 | 1.3532e+01 | 1.3183e+01 | 1.3358e+01 | 1.7444e-01 | 0.0131 |
| 5 | 5.7000 | 32.4902 | 1.9442e+01 | 1.9096e+01 | 1.9269e+01 | 1.7310e-01 | 0.0090 |

`T_inner > T_outer` for every n — the 5D Tangherlini centrifugal + throat-curvature potential shifts the radial profile toward the throat side. Asymmetry decreases monotonically with n (shell saturation).

## T3: `χ_n` values from boundary stress

| n | χ_n = T_odd | χ_n/ω² |
|---:|---:|---:|
| 0 | 1.5600e-01 | 0.1402 |
| 1 | 2.0031e-01 | 0.0514 |
| 2 | 1.8514e-01 | 0.0221 |
| 3 | 1.7734e-01 | 0.0121 |
| 4 | 1.7444e-01 | 0.0077 |
| 5 | 1.7310e-01 | 0.0053 |

## T4: Sign pattern audit vs PR #78's existence proof

| n | derived χ_n sign | PR #78 ansatz required | match? |
|---:|:---:|:---:|:---:|
| 3 | + | − | ✗ |
| 4 | + | + | ✓ |
| 5 | + | + | ✓ |

Boundary stress gives sign-flipping pattern: **False** (2/3 match). PR #78's sign-flipping ansatz is structurally overruled by the boundary-stress derivation.

## T5: Magnitude audit — derived vs required

| splitting | required χ/ω² | derived χ_n/ω² | deficit factor |
|---|---:|---:|---:|
| u/d (n=3) | 0.648 | 0.0121 | 53.4× |
| s/c (n=4) | 0.989 | 0.0077 | 128.5× |
| b/t (n=5) | 0.999 | 0.0053 | 187.5× |

Boundary stress alone gives `χ_n/ω²` 50–100× too small to drive the observed within-generation mass² ratios. Within-generation splittings ⟹ PR #80's color sector.

## T6: Singlet projector placeholder

P_S = 6×6 identity (flavor-level basis; color implicit). Spectrum unchanged under projection: **True**; commutes with H_flavor: **True**.

Diagonal mass² with boundary-stress `χ_n` (species order u, d, c, s, t, b — v3 species map):

| species | m² |
|---|---:|
| u | 14.8051 |
| d | 14.4504 |
| c | 22.8403 |
| s | 22.4914 |
| t | 32.6633 |
| b | 32.3171 |

## T7: What PR #79 fixes vs what PR #80 owes

**PR #79 fixes:**

  - χ_n = T_odd(n) from cavity-mouth boundary stress
  - sign uniform (T_inner > T_outer) — no sign flip
  - magnitude shell-suppressed (χ_n/ω² ~ 0.01–0.02 for n ≥ 3)
  - singlet projector placeholder (identity on flavor basis)
  - no free parameter once cavity geometry is fixed

**PR #80 owes:**

  - color algebra acting on (l, n, p)
  - H_couple inter-mode mixing
  - within-generation splittings beyond boundary-stress χ_n
  - inter-generation hierarchy spanning ~9 orders of mass²
  - singlet projector populated with non-trivial color content

**v3 species map status:** flagged for revision: under boundary-stress reading, natural assignment is "+ = heavier" uniformly, which inverts v3's (k=1, +) = u to (n=3, +) = d. PR #80 settles which interpretation is physical.

## Verdict

**CHI_N_DERIVED_BOUNDARY_STRESS_INSUFFICIENT.** CHI_N DERIVED FROM BOUNDARY STRESS, INSUFFICIENT FOR OBSERVED SPLITTINGS. The PR #79 derivation pins the structural slot of the generation-dependent Z₂ partition splitter χ_n: it is the Z₂-antisymmetric piece of the cavity-mouth boundary stress T_odd(n) = (T_inner − T_outer)/2, where T_inner = (∂ψ_n/∂r*)²|R_MID and T_outer = (∂ψ_n/∂r*)²|R_OUTER for the normalized radial eigenfunction. No free parameter — once the cavity geometry [R_MID, R_OUTER] is fixed, χ_n is determined by the Tangherlini eigensolver.

STRUCTURAL PATTERN. For l=1, n=0..5 on the standard cavity:
  - T_inner > T_outer for EVERY mode (uniform positive sign of T_odd). The radial profile is shifted toward R_MID by the 5D Tangherlini centrifugal + throat-curvature potential (the inner mouth has the stronger curvature).
  - The asymmetry T_odd/T_even DECREASES monotonically with n: 0.30 (electron-focused) → 0.04 (tau, edge of shell saturation) → 0.02, 0.01, 0.009 (shell sector). Shell-saturated modes feel the mouth asymmetry only weakly because the standing wave fills the cavity nearly uniformly.
  - NO sign flip between n=3 and n=4. The PR #78 existence proof (χ_3 < 0, χ_4 > 0, χ_5 > 0) is STRUCTURALLY INCOMPATIBLE with the boundary-stress derivation — the sign-flipping ansatz cannot come from cavity-mouth asymmetry.

MAGNITUDE INSUFFICIENT. For shell modes (n ≥ 3), χ_n/ω² is in the range 0.01–0.02, while the required χ/ω² to reproduce the observed within-generation mass² ratios is 0.65 (u/d), 0.99 (s/c), 0.95 (b/t). Boundary stress is too small by 50–100× to drive the observed splittings on its own.

V3 SPECIES MAP FLAGGED. The v3 convention (k=1, +) = u (with u lighter than d) is incompatible with the natural boundary-stress reading "+ = heavier" (uniform across n). Either the v3 map needs revision (assign + to the heavier state at every generation, so (n=3, +) = d, (n=3, −) = u), or the within-generation inversion at n=3 is a different mechanism entirely — PR #80's color sector is the natural place to settle this.

SINGLET CONSTRAINT. The 6-state shell basis is at the FLAVOR level (u, d, s, c, b, t); color is implicit. The singlet projector P_S is the identity on this basis (physical states = color-singlet observables; color-octet components, if present in any hypothetical extension, would project out under singlet projection). P_S commutes with the diagonal flavor Hamiltonian; the spectrum is unchanged. The non-trivial singlet content arrives in PR #80 with the color algebra identification.

WHAT PR #79 FIXES. (i) Structural origin of χ_n = cavity-mouth boundary stress; (ii) uniform-sign / shell-suppressed pattern; (iii) no free parameter once cavity geometry is fixed; (iv) singlet projector placeholder (identity on flavor basis); (v) PR #78 sign-flipping ansatz overruled.

WHAT PR #80 OWES. (i) Color algebra acting on (l, n, p) — SU(3), SU(2)×Z₂, Pati-Salam SU(4), or other; (ii) H_couple inter-mode mixing populated by color-algebra transformation rules; (iii) within-generation splittings beyond the small boundary-stress χ_n; (iv) inter-generation hierarchy spanning ~9 orders of magnitude; (v) settling the v3 species ↔ partition map question. After PR #80 populates the operator, the n_part audit can be re-run; until then, n_part = 233 remains a phenomenological compensator with sharpened scope.

HONEST SCOPE. PR #79 derives χ_n structurally; identifies that it is small and uniform-sign for shell modes; rules out the boundary stress as the source of the within-generation inversion; and flags PR #80's color sector as the necessary next contribution. The n_part question remains open until PR #80.

## What this leaves open

- **PR #80** — color algebra acting on `(l, n, p)`; `H_couple` inter-mode mixing populated by color-algebra transformation rules; within-generation splittings beyond boundary-stress `χ_n`; inter-generation hierarchy.
- **v3 species ↔ partition map** — flagged for revision; PR #80 settles whether "+ = heavier" uniformly (consistent with boundary stress) or the n=3 case is genuinely inverted by a different mechanism.
- **`n_part` audit re-run** — after PR #80 populates the operator; until then `n_part = 233` remains a residual phenomenological compensator with sharpened scope.
