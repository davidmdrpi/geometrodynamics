# Geometric throat arc synthesis (PR #131)

**Run:** 2026-06-04T20:49:47+00:00

Capstone of the geometric throat arc — PRs #116 and #127–#130. Re-verifies a keystone from each arc member together, consolidates the unified picture (the antipodal identification of the 5D Tangherlini horizon), and lays out the epistemic ledger.

- **Horizon (#116/#127)**: f(rs)=0, K=72rs⁴/r⁸ regular, T_H=1/2πrs, k₅=D_bulk=5 (#116/#127)
- **Regular & antipodal (#128)**: EF/Kruskal regular, proper=√(2rs ε), F(rs)=4e⁻², (U,V)→(−U,−V)=C-swap (#128)
- **Antipodal BC (#129)**: Y_l(−x)=(−1)^l ⟹ even-l Neumann/odd-l Dirichlet, unitary mirror (#129)
- **Spectrum (#130)**: antipodal real undamped (stable matter) vs absorbing complex ringdown (#130)
- **Unified object**: one antipodal primitive, five faces (#58/#63/#128/#129/#130)
- **Open**: AdS scale k=κ₅²/Λ₅ (#112); nucleation rate (#58/#88); global brane solution; horizon QNM tower (#130)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | synthesise the geometric throat arc (#116, #127–#130) | **PASS** |
| T2 | `T2_throat_is_5d_tangherlini_horizon` | throat = 5D Tangherlini horizon: f(rs)=0, K=72rs⁴/r⁸, T_H=1/2πrs | **PASS** |
| T3 | `T3_horizon_regular_and_antipodal` | horizon-regular & antipodal: EF det finite, F(rs)=4e⁻², (U,V)→(−U,−V) | **PASS** |
| T4 | `T4_antipodal_l_parity_bc` | antipodal BC: Y_l(−x)=(−1)^l ⟹ even-N/odd-D, unitary mirror | **PASS** |
| T5 | `T5_spectral_consequence` | antipodal real undamped vs absorbing complex ringdown | **PASS** |
| T6 | `T6_one_primitive_five_faces` | one antipodal primitive, five faces (#58/#63/#128/#129/#130) | **PASS** |
| T7 | `T7_epistemic_ledger` | epistemic ledger: derived / postulated / open | **PASS** |
| T8 | `T8_assessment` | GEOMETRIC_THROAT_ARC_SYNTHESIS_ANTIPODAL_HORIZON_UNIFIED | **PASS** |

## The arc keystones (re-verified together)

| PR | keystone | value |
|---|---|---|
| #116/#127 | Kretschmann at throat (regular) | K = 72.0 |
| #116/#127 | Hawking temperature | T_H = 0.159155 = 1/(2π rs) |
| #128 | EF det at throat (nondegenerate) | -0.29904 |
| #128 | Kruskal factor at throat | F(rs) = 0.54134 = 4 e⁻² |
| #128 | proper distance to throat | √(2 rs ε) = 0.2 |
| #129 | antipodal BC (l=0 fundamental) | real ω = 1.186+0.000i |
| #130 | absorbing BC (l=0 fundamental) | complex ω = 1.893-1.159i |

## One primitive, five faces

| face | PR | form |
|---|---|---|
| C inner/outer swap | #63 | c₁ → −c₁ |
| throat ↔ antithroat | #58 | nucleation channel |
| Kruskal antipode | #128 | (U,V,Ω) → (−U,−V,Ω̄) |
| unitary-mirror BC | #129 | l-parity Neumann/Dirichlet |
| stable-matter selector | #130 | real undamped spectrum |

All five are the same geometric object — the antipodal identification of the 5D Tangherlini horizon. **"Bulk Antipodal Mechanics" is the mechanics of this one identification.**

## Verdict

**GEOMETRIC_THROAT_ARC_SYNTHESIS_ANTIPODAL_HORIZON_UNIFIED.** THE GEOMETRIC THROAT ARC UNIFIES INTO ONE PICTURE: THE ANTIPODAL IDENTIFICATION OF THE 5D TANGHERLINI HORIZON. PRs #116 and #127–#130 lifted the BAM throat from a radial cavity operator to its full five-dimensional geometry and read off the physics; this capstone re-verifies the arc's keystones together and consolidates them.

THE THROAT IS THE 5D TANGHERLINI HORIZON (#116/#127). f(rs) = 0; the Kretschmann scalar K = 72 rs⁴/r⁸ is finite on the whole cavity (the only true singularity is at r = 0); T_H = 1/(2π rs) carries the closure quantum; and the cavity potential's coefficients are the D=5 reductions (centrifugal l(l+D−3), curvature D−2), so k₅ = D_bulk = 5. The throat's parent is a genuine, curvature-regular D=5 vacuum.

IT IS HORIZON-REGULAR AND ANTIPODAL (#128). Eddington–Finkelstein and Kruskal coordinates remove the coordinate singularity (det g = −r⁶ sin⁴χ sin²θ finite at the throat, proper distance √(2 rs ε) = the ε healing length, F_Kruskal(rs) = 4 e⁻²), and the maximal extension's antipodal map (U,V) → (−U,−V) preserves UV — it is the throat ↔ antithroat C-swap.

THE ANTIPODAL IDENTIFICATION FIXES THE BC (#129). The S³ harmonics carry Y_l(−x) = (−1)^l Y_l(x), so the antipodal identification forces the wave BC by l-parity — even-l Neumann, odd-l Dirichlet — a real BC with zero throat flux: a unitary mirror.

THE SPECTRUM FOLLOWS (#130). The antipodal BC gives a real, undamped spectrum (sharp, stable matter); an absorbing horizon gives complex frequencies (Im ω < 0, damped ringdown). Stable matter requires the unitary antipodal throat.

ONE PRIMITIVE, FIVE FACES. The antipodal identification of the 5D Tangherlini horizon appears across the program as the C = inner/outer swap (#63), the throat ↔ antithroat nucleation channel (#58), the antipodal Kruskal map (#128), the l-parity unitary-mirror BC (#129), and the selector of the real, stable matter spectrum (#130). "Bulk Antipodal Mechanics" is the mechanics of this one identification.

THE EPISTEMIC LEDGER. DERIVED: the throat's parent bulk is a genuine D=5 Tangherlini vacuum (Ricci-flat, curvature-regular cavity, S³ horizon, T_H = 1/2πrs, k₅ = D_bulk); the coordinate singularity is removable; the antipodal identification fixes the l-parity BC and the unitary mirror; the antipodal spectrum is real (stable matter), the absorbing one complex. POSTULATED: the antipodal identification itself — BAM's defining axiom — shown self-consistent (unitary, stable-matter-supporting), not forced. OPEN: the exact AdS scale k = κ₅²/Λ₅ (#112); the dynamical throat ↔ antithroat nucleation rate (#58/#88); the global brane-localised solution; the idealised r* → −∞ horizon QNM tower and GW coupling (#130).

SCOPE. A synthesis/consistency capstone: it re-verifies the arc's keystones together and organises them; it does not add new derivations, remove any open item, or strengthen the antipodal postulate from "self-consistent" to "forced".
