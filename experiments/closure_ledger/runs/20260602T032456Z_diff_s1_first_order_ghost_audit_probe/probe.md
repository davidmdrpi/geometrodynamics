# First-order Diff(S¹) Faddeev–Popov ghost audit (PR #118)

**Run:** 2026-06-02T03:24:56+00:00

Audits the PR #117 ghost determinant rigorously: distinguishes `P = ∂_τ`, `P†P = −∂_τ²`, `det'(P)`, `det'(P†P)`; handles the phase/η-invariant of `det'(∂_τ)`; treats CKV zero-mode division and zero-mode norms; and fixes the ghost L-power. **Result: the FP ghost contributes `L` (first order); `L²` only under an explicit second-order ghost convention.**

- **Four objects**: P=∂_τ (1st, 1 zero mode=CKV); P†P=−∂_τ² (2nd, 1 zero mode); det'(P)=L; det'(P†P)=L²
- **η-invariant**: η(−i∂_τ) = 0 (symmetric spectrum) ⟹ det'(∂_τ)=+L, no phase; antiperiodic Möbius: η=0, no CKV
- **Convention**: first-order ghost = det'(P) = L; second-order = det'(P†P) = L² (over-counts by L)
- **No double-count**: det'(P) (primed, nonzero modes) excludes the CKV; CKV norm only in Vol(CKG) — divided once (proven by SVD: 1 zero singular value)
- **Measure**: Z = Σ_sectors ∫ (dL/L) det^{−1/2}_matter e^{−S}; net dL·L^{−1−d/2}; 1/L = single CKV factor
- **L-power fixed**: ghost L¹ in det'(P); L² only with a second-order convention

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | audit goal: distinguish the 4 objects, fix the L-power | **PASS** |
| T2 | `T2_four_objects_and_spectra` | P=∂_τ (1st, eigvals 2πin/L), P†P=−∂_τ² (2nd); 1 zero mode each | **PASS** |
| T3 | `T3_two_determinants` | det'(P†P)=L²; |det'(P)|=det'(P†P)^{1/2}=L (computed) | **PASS** |
| T4 | `T4_eta_invariant_phase` | η(−i∂_τ)=0 ⟹ det'(∂_τ)=+L, no phase; antiperiodic: no CKV | **PASS** |
| T5 | `T5_first_vs_second_order_convention` | first-order det'(P)=L vs second-order det'(P†P)=L² (×L) | **PASS** |
| T6 | `T6_no_double_counting_ckv_division` | no double-count: SVD 1 zero SV; CKV norm only in Vol(CKG) (once) | **PASS** |
| T7 | `T7_measure_table_single_counted` | single-counted measure Z=Σ∫(dL/L)det^{−1/2}_matter e^{−S}; net dL·L^{−1−d/2} | **PASS** |
| T8 | `T8_assessment` | BAM_LOOP_MEASURE_GHOST_FIRST_ORDER_DETPRIME_P_EQUALS_L_NET_dL_OVER_L | **PASS** |

## The two determinants

| L | det′(P†P) | det′(P) | det′(P†P)^½ | match |
|---:|---:|---:|---:|:---:|
| 6.28319 | 39.47842 | 6.28319 | 6.28319 | ✓ |
| 1.0 | 1.0 | 1.0 | 1.0 | ✓ |
| 3.32408 | 11.04951 | 3.32408 | 3.32408 | ✓ |
| 5.0 | 25.0 | 5.0 | 5.0 | ✓ |

`det'(P†P) = L²` (second order); `det'(P) = det'(P†P)^{1/2} = L` (first order). The η-invariant of `−i∂_τ` vanishes (symmetric spectrum), so `det'(∂_τ) = +L` with no anomalous phase.

## No double-counting of the CKV (SVD partition)

- ghost space = `ker(P)`[CKV] ⊕ `ker(P†)`[modulus] ⊕ (nonzero modes); SVD of `∂_τ`: **1 zero singular value** (right-null = CKV, left-null = modulus), **200 nonzero** in `det'(P)`.
- `det'(P)` (primed, nonzero modes) **excludes both kernels** ⟹ the CKV norm `‖φ‖` enters **only** `Vol(CKG)`, the modulus norm `‖ψ‖` **only** `dL` — each divided **once**.
- The first draft's extra `÷‖φ_CKV‖` (the `√L·√L` norms) alongside `1/Vol(CKG)` double-counted the single CKV; removed here.

## The corrected measure table (single-counted)

| factor | zero mode / sector (once) | L-power |
|---|---|---|
| modulus measure dL | ker(P†) = Teichmüller (once) | `dL` |
| CKV division 1/Vol(CKG) = 1/L | ker(P) = CKV (once; Vol from ‖φ‖) | `L^{-1}` |
| det'(P) = L (ghost nonzero) | nonzero modes (excl. kernels) | `folds into matter Tr e^{-LH}` |
| matter det^{-1/2} | Tr e^{-LH} (incl. ghost det) | `L^{-d/2}` |
| NET MEASURE | product (each zero mode once) | `dL · L^{-1-d/2}` |

**Measure:** `Z = Σ_sectors ∫ (dL/L) · det^{−1/2}_matter · e^{−S}`. The single `1/L` is the CKV factor (`1/L = 1/Vol(CKG) (single; = PR #74 1/(2π) at L=2π)`); `det'(P) = L` folds into the matter heat kernel. Net L-power: `dL · L^{-1-d/2}`.

## Verdict

**BAM_LOOP_MEASURE_GHOST_FIRST_ORDER_DETPRIME_P_EQUALS_L_NET_dL_OVER_L.** THE Diff(S¹) FP GHOST IS FIRST-ORDER: det'(P) = L; THE L² POWER APPEARS ONLY UNDER AN EXPLICIT SECOND-ORDER GHOST CONVENTION. This audit follows PR #117 and distinguishes the four objects the review named.

THE FOUR OBJECTS. P = ∂_τ is first order (eigenvalues 2πin/L, one zero mode = the CKV); P†P = −∂_τ² is second order (eigenvalues (2πn/L)², one zero mode). Their zeta-regularized determinants are det'(P†P) = L² and |det'(P)| = det'(P†P)^{1/2} = L (verified to machine precision).

PHASE / η-INVARIANT. For A = −i∂_τ (self-adjoint, eigenvalues 2πn/L) the spectrum is symmetric under n → −n, so the η-invariant vanishes identically, η(0) = 0. Hence det'(∂_τ) carries no anomalous phase and equals +L (real, positive). In the antiperiodic / Möbius sector the eigenvalues are 2πi(n+½)/L — still symmetric (η = 0) but with NO zero mode, so there is no CKV in the non-orientable sector (a clean tie to the odd-k structure).

FIRST- VS SECOND-ORDER CONVENTION. The physical FP for one real reparametrization symmetry is the first-order bc system, Δ_FP = ∫ Db Dc e^{−∮ b ∂_τ c} = det'(P) = L. A second-order ghost convention (action b(P†P)c, i.e. a doubled/complex ghost) would instead give det'(P†P) = L². The two differ by exactly one power of L; the minimal Faddeev–Popov construction is first-order ⟹ L¹.

CKV DIVISION — NO DOUBLE-COUNTING. The ghost field space splits orthogonally as ker(P)[CKV] ⊕ ker(P†)[modulus] ⊕ (nonzero modes), and the FP determinant det'(P) = det'(P†P)^{1/2} is the PRIMED determinant over the nonzero modes only. The SVD of ∂_τ shows exactly one zero singular value — its right-null vector is the CKV, its left-null the modulus — and N−1 nonzero singular values in det'(P). So the CKV norm ‖φ‖ enters ONLY the gauge-orbit volume Vol(CKG), the modulus norm ‖ψ‖ ONLY the dL measure, and det'(P) (a primed det over nonzero modes) excludes both kernels — each zero mode is divided exactly once. The first draft of PR #118 divided additionally by the zero-mode norms (√L·√L) on top of the CKV automorphism 1/Vol(CKG); since the CKV norm is already inside Vol(CKG), that double-counted the single CKV. Corrected here.

THE MEASURE (single-counted). Each zero mode appears once: the modulus (ker P†) as dL, the CKV (ker P) as 1/Vol(CKG) = 1/L (Vol computed from ‖φ‖), and det'(P) = L (ghost nonzero modes) folds into the matter heat kernel Tr e^{−LH} per the standard FP procedure — it is not an independent floating factor. With the matter determinant det^{−1/2}_matter ~ L^{−d/2}, the BAM loop measure is Z = Σ_sectors ∫ (dL/L) det^{−1/2}_matter e^{−S}, net L-power dL·L^{−1−d/2}, with the single 1/L the CKV factor (= PR #74's 1/(2π) at the closure loop L = 2π). The absolute normalization (the κ₅²/Λ₅ anchor) and the multi-loop measure remain open.

## What this fixes (and leaves open)

- **Fixed:** the ghost L-power and the CKV counting. The FP ghost is first-order, `det'(P) = det'(P†P)^{1/2} = L` (`η = 0`, real positive); `L²` arises only under an explicit second-order ghost convention. The CKV is divided **exactly once** — `det'(P)` (primed, nonzero modes) excludes the CKV (SVD: 1 zero singular value), so the CKV norm enters only `Vol(CKG)`; the first draft's extra `√L·√L` norm division was the double-count, now removed. Measure `Z = Σ ∫ (dL/L) det^{−1/2}_matter e^{−S}` (net `dL·L^{−1−d/2}`), the single `1/L` being the CKV factor `= 1/Vol(CKG)`.
- **Open:** the absolute normalization of `Z` (the `κ₅²/Λ₅` bulk anchor, PR #112) and the multi-loop measure.
