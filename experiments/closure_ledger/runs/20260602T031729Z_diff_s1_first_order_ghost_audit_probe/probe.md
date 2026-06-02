# First-order Diff(S¹) Faddeev–Popov ghost audit (PR #118)

**Run:** 2026-06-02T03:17:29+00:00

Audits the PR #117 ghost determinant rigorously: distinguishes `P = ∂_τ`, `P†P = −∂_τ²`, `det'(P)`, `det'(P†P)`; handles the phase/η-invariant of `det'(∂_τ)`; treats CKV zero-mode division and zero-mode norms; and fixes the ghost L-power. **Result: the FP ghost contributes `L` (first order); `L²` only under an explicit second-order ghost convention.**

- **Four objects**: P=∂_τ (1st, 1 zero mode=CKV); P†P=−∂_τ² (2nd, 1 zero mode); det'(P)=L; det'(P†P)=L²
- **η-invariant**: η(−i∂_τ) = 0 (symmetric spectrum) ⟹ det'(∂_τ)=+L, no phase; antiperiodic Möbius: η=0, no CKV
- **Convention**: first-order ghost = det'(P) = L; second-order = det'(P†P) = L² (over-counts by L)
- **Zero modes**: norm √L each; 2 zero modes (CKV+modulus) ⟹ L; ghost NET first-order L/L=1
- **Measure**: Z = Σ_sectors ∫ (dL/L) det^{−1/2}_matter e^{−S}; net dL·L^{−1−d/2}
- **L-power fixed**: ghost L¹ in det'(P); net L⁰ after norms; L² only with a second-order convention

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | audit goal: distinguish the 4 objects, fix the L-power | **PASS** |
| T2 | `T2_four_objects_and_spectra` | P=∂_τ (1st, eigvals 2πin/L), P†P=−∂_τ² (2nd); 1 zero mode each | **PASS** |
| T3 | `T3_two_determinants` | det'(P†P)=L²; |det'(P)|=det'(P†P)^{1/2}=L (computed) | **PASS** |
| T4 | `T4_eta_invariant_phase` | η(−i∂_τ)=0 ⟹ det'(∂_τ)=+L, no phase; antiperiodic: no CKV | **PASS** |
| T5 | `T5_first_vs_second_order_convention` | first-order det'(P)=L vs second-order det'(P†P)=L² (×L) | **PASS** |
| T6 | `T6_ckv_division_and_zero_mode_norms` | zero-mode norms √L·√L=L ⟹ ghost net L/L=1 (1st), L²/L=L (2nd) | **PASS** |
| T7 | `T7_measure_table_and_net_L_power` | measure Z=Σ∫(dL/L)det^{−1/2}_matter e^{−S}; net dL·L^{−1−d/2} | **PASS** |
| T8 | `T8_assessment` | BAM_LOOP_MEASURE_GHOST_FIRST_ORDER_DETPRIME_P_EQUALS_L_NET_dL_OVER_L | **PASS** |

## The two determinants

| L | det′(P†P) | det′(P) | det′(P†P)^½ | match |
|---:|---:|---:|---:|:---:|
| 6.28319 | 39.47842 | 6.28319 | 6.28319 | ✓ |
| 1.0 | 1.0 | 1.0 | 1.0 | ✓ |
| 3.32408 | 11.04951 | 3.32408 | 3.32408 | ✓ |
| 5.0 | 25.0 | 5.0 | 5.0 | ✓ |

`det'(P†P) = L²` (second order); `det'(P) = det'(P†P)^{1/2} = L` (first order). The η-invariant of `−i∂_τ` vanishes (symmetric spectrum), so `det'(∂_τ) = +L` with no anomalous phase.

## The measure table (net L-power)

| factor | source | L-power |
|---|---|---|
| moduli-space measure dL/L | modulus dL × CKV automorphism 1/L | `L^{-1} dL` |
| FP ghost det'(P) (first order) | nonzero ghost modes | `L^{+1}` |
| ghost zero-mode norms (÷) | 1 CKV + 1 modulus, √L each | `L^{-1}` |
| ghost NET | det'(P)/norms | `L^{0}` |
| matter det^{-1/2} | Tangherlini/heat-kernel Weyl | `L^{-d/2}` |
| NET MEASURE (first order) | product | `dL · L^{-1-d/2}` |

**Measure:** `Z = Σ_sectors ∫ (dL/L) · det^{−1/2}_matter · e^{−S}`. Net L-power (first order): `dL · L^{-1-d/2}`. Second-order convention: `dL · L^{-d/2}  (one spurious power of L)`.

## Verdict

**BAM_LOOP_MEASURE_GHOST_FIRST_ORDER_DETPRIME_P_EQUALS_L_NET_dL_OVER_L.** THE Diff(S¹) FP GHOST IS FIRST-ORDER: det'(P) = L; THE L² POWER APPEARS ONLY UNDER AN EXPLICIT SECOND-ORDER GHOST CONVENTION. This audit follows PR #117 and distinguishes the four objects the review named.

THE FOUR OBJECTS. P = ∂_τ is first order (eigenvalues 2πin/L, one zero mode = the CKV); P†P = −∂_τ² is second order (eigenvalues (2πn/L)², one zero mode). Their zeta-regularized determinants are det'(P†P) = L² and |det'(P)| = det'(P†P)^{1/2} = L (verified to machine precision).

PHASE / η-INVARIANT. For A = −i∂_τ (self-adjoint, eigenvalues 2πn/L) the spectrum is symmetric under n → −n, so the η-invariant vanishes identically, η(0) = 0. Hence det'(∂_τ) carries no anomalous phase and equals +L (real, positive). In the antiperiodic / Möbius sector the eigenvalues are 2πi(n+½)/L — still symmetric (η = 0) but with NO zero mode, so there is no CKV in the non-orientable sector (a clean tie to the odd-k structure).

FIRST- VS SECOND-ORDER CONVENTION. The physical FP for one real reparametrization symmetry is the first-order bc system, Δ_FP = ∫ Db Dc e^{−∮ b ∂_τ c} = det'(P) = L. A second-order ghost convention (action b(P†P)c, i.e. a doubled/complex ghost) would instead give det'(P†P) = L². The two differ by exactly one power of L; the minimal Faddeev–Popov construction is first-order ⟹ L¹.

CKV DIVISION + ZERO-MODE NORMS. The constant mode on a circle of circumference L has norm √L. The ghost system has two zero modes — the c-ghost (the CKV, rigid rotation) and the b-ghost (dual of the one Teichmüller modulus) — with norms √L·√L = L. The ghost NET factor det'(P)/norms is therefore L/L = 1 in the first-order convention, but L²/L = L in the second-order one (a spurious extra power).

THE MEASURE. The moduli-space measure (one modulus, one CKV automorphism) is the standard worldline dL/L. With the first-order ghost NET = 1 and the matter determinant det^{−1/2}_matter ~ L^{−d/2} (the heat-kernel Weyl term), the BAM loop measure is Z = Σ_sectors ∫ (dL/L) det^{−1/2}_matter e^{−S}, net L-power dL·L^{−1−d/2}. The second-order convention would multiply by an extra L (dL·L^{−d/2}), confirming that the first-order det'(P) = L is the correct fixing. The absolute normalization (the κ₅²/Λ₅ anchor) and the multi-loop measure remain open.

## What this fixes (and leaves open)

- **Fixed:** the ghost L-power. The FP ghost is first-order, `det'(P) = det'(P†P)^{1/2} = L` (`η = 0`, real positive); after the CKV + modulus zero-mode norms (`√L` each) the ghost NET factor is `1`, giving the measure `Z = Σ ∫ (dL/L) det^{−1/2}_matter e^{−S}` (net `dL·L^{−1−d/2}`). `L²` arises only under an explicit second-order ghost convention, which over-counts by one power of `L`.
- **Open:** the absolute normalization of `Z` (the `κ₅²/Λ₅` bulk anchor, PR #112) and the multi-loop measure.
