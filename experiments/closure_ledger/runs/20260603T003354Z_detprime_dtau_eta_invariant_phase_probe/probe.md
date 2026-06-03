# The phase / η-invariant framework for det′(∂_τ) (PR #119)

**Run:** 2026-06-03T00:33:54+00:00

Builds the full mathematical framework for the phase of the first-order determinant `det'(∂_τ)` (PR #118 only asserted `η = 0`). The phase is governed by the **Singer/APS formula** `det'(A) = |det'|·exp[±i(π/2)(ζ_{|A|}(0) − η_A(0))]`, splitting into a local (scaling) and a topological (η-invariant) piece. Both BAM closure sectors sit at **η = 0** (real determinant).

- **Phase formula**: `det'(A) = |det'|·exp[±i(π/2)(ζ_{|A|}(0) − η_A(0))]`
- **Split**: phase = (π/2)ζ(0) [local/scaling] − (π/2)η(0) [topological/asymmetry]
- **η of flux**: η_A(0) = 1 − 2a (Hopf holonomy a = kχ/2π)
- **BAM sectors**: orientable a=0 ⟹ η=0; Möbius a=1/2 ⟹ η=0; both ⟹ real det'
- **Concrete**: det'_periodic = L; det_antiperiodic = 2
- **Justifies**: PR #118 det'(∂_τ) = +L (no anomalous phase)
- **Open**: genuine η-phase for intermediate holonomy; interplay with matter phase

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | build the phase/η framework (PR #118 asserted η = 0) | **PASS** |
| T2 | `T2_setup_and_modulus` | ∂_τ anti-s.a., A = −i∂_τ s.a.; |det′| = L unambiguous | **PASS** |
| T3 | `T3_phase_formula` | det′(A) = |det′|·exp[±i(π/2)(ζ(0) − η(0))]; local + topological | **PASS** |
| T4 | `T4_eta_invariant_with_flux` | η(0) = 1 − 2a; reduced η ≡ 0 for periodic & antiperiodic | **PASS** |
| T5 | `T5_concrete_determinants` | det′(∂_τ)_periodic = L; det(∂_τ)_antiperiodic = 2 | **PASS** |
| T6 | `T6_bam_hopf_holonomy_connection` | Hopf holonomy a = kχ/2π ⟹ η = 1 − 2a; a=0, 1/2 ⟹ η = 0 | **PASS** |
| T7 | `T7_scope` | rigorously justifies PR #118; intermediate-holonomy η-phase open | **PASS** |
| T8 | `T8_assessment` | DETPRIME_DTAU_PHASE_ETA_INVARIANT_FRAMEWORK_BAM_SECTORS_ETA_ZERO | **PASS** |

## The η-invariant vs the flux (Hopf holonomy)

| flux a | η(0) = 1 − 2a | symmetric? |
|---:|---:|:---:|
| 0.0 | 1.0 |  |
| 0.25 | 0.5 |  |
| 0.5 | 0.0 | ✓ (η = 0) |
| 0.75 | -0.5 |  |

`a = 0` (orientable, periodic) and `a = 1/2` (Möbius, antiperiodic) both give `η = 0` — the symmetric spectra where `det'(∂_τ)` is real. Generic `a` gives a genuine η-phase.

## Concrete determinants (closed forms)

| sector | closed form | m → 0 | det |
|---|---|---|---|
| periodic (orientable) | `2 sinh(mL/2)` | → 0 (zero mode = CKV) | `det′(∂_τ) = L = 6.28319` |
| antiperiodic (Möbius) | `2 cosh(mL/2)` | → 2 (no zero mode) | `det(∂_τ) = 2.0` |

The periodic residue reproduces PR #118's `det'(∂_τ) = L` (real, positive) — now *derived*, not asserted.

## Verdict

**DETPRIME_DTAU_PHASE_ETA_INVARIANT_FRAMEWORK_BAM_SECTORS_ETA_ZERO.** THE PHASE OF det'(∂_τ) IS GOVERNED BY THE η-INVARIANT: det'(A) = |det'|·exp[±i(π/2)(ζ(0) − η(0))], AND BOTH BAM CLOSURE SECTORS SIT AT η = 0 (REAL DETERMINANT). PR #118 asserted η(−i∂_τ) = 0 from the n → −n symmetry; this probe builds the full framework.

SETUP AND MODULUS. P = ∂_τ is anti-self-adjoint (eigenvalues 2πin/L, imaginary); A = −i∂_τ is self-adjoint (μ_n = 2πn/L). The modulus |det'(∂_τ)| = det'(P†P)^{1/2} = L is unambiguous; the phase is the subtle part.

THE PHASE FORMULA (Singer / APS). Choosing a branch for the negative eigenvalues, ζ_A(s) = ζ_+(s) + e^{∓iπs}ζ_−(s), gives det'(A) = |det'(A)|·exp[±i(π/2)(ζ_{|A|}(0) − η_A(0))]. The phase splits into a LOCAL (heat-kernel / scaling) piece (π/2)ζ_{|A|}(0) and an INTRINSIC (spectral-asymmetry) piece −(π/2)η_A(0). For periodic A: ζ_{|A|}(0) = 2ζ_R(0) = −1 and η(0) = 0, so the phase is the −π/2 scaling rotation, removed by the symmetric branch ⟹ det'(∂_τ) real.

THE η-INVARIANT AND THE FLUX. Threading a U(1) holonomy a ∈ [0,1) (eigenvalues 2π(n+a)/L), the Hurwitz zeta ζ_H(0,a) = ½ − a gives η_A(0) = ζ_H(0,a) − ζ_H(0,1−a) = 1 − 2a — the η-invariant tracks the holonomy linearly. The reduced η of the symmetric sectors vanishes identically: periodic with the zero mode (the CKV) removed, η' ≡ 0; antiperiodic {2π(n+½)/L}, η ≡ 0. (The naive value 1 at a = 0 is the spectral-flow jump as the zero mode crosses.)

CONCRETE DETERMINANTS. With an IR mass, det(∂_τ + m)_periodic = 2 sinh(mL/2) → 0 (the zero mode), whose residue gives det'(∂_τ)_periodic = L (real, positive — exactly PR #118's value, now derived); det(∂_τ + m)_antiperiodic = 2 cosh(mL/2) → 2, so det(∂_τ)_antiperiodic = 2 (L-independent, no zero mode).

THE BAM CONNECTION. The Hopf holonomy ∮A = e^{ikχ} is the flux a = kχ/2π threading the closure loop, and η(0) = 1 − 2a is the spectral-asymmetry phase it induces. The non-orientable (Möbius) Z₂ half-twist is the antiperiodic shift a = 1/2 ⟹ η = 0; the orientable sector is a = 0 ⟹ η = 0. Both physical BAM sectors are η = 0 (real determinant), which is precisely why PR #118's det'(∂_τ) = +L carries no anomalous phase. A generic intermediate holonomy (0 < a < 1/2) would give a genuine η-phase exp[−i(π/2)(1−2a)].

SCOPE. ESTABLISHED: the phase formula, the local/topological split, the periodic (det' = L) and antiperiodic (det = 2) sectors, and the rigorous justification of PR #118 (both BAM sectors η = 0 ⟹ real det'). OPEN: the genuine η-phase for intermediate Hopf holonomy and its interplay with the matter determinant phase in the full measure.

## What this establishes (and leaves open)

- **Establishes:** the phase formula `det′(A) = |det′|·exp[±i(π/2)(ζ(0) − η(0))]`, the local/topological split, the periodic (`det′ = L`) and antiperiodic (`det = 2`) sectors, and the rigorous justification of PR #118 — both BAM sectors (orientable `a=0`, Möbius `a=1/2`) have `η = 0`, so `det′(∂_τ)` is real with no anomalous phase.
- **Open:** the genuine η-phase `exp[−i(π/2)(1−2a)]` for intermediate Hopf holonomy, and its interplay with the matter determinant phase in the full measure.
