# Faddeev–Popov / ghost determinant for the Diff(S¹) quotient (PR #117, revised)

**Run:** 2026-06-02T02:42:10+00:00

Supplies the gauge sector of the S_BAM measure (PR #115 flagged it; PR #116 did the matter determinant). **Revised per review:** the ghost determinant is `det'(P†P)^{1/2} = det'(P) = L` — the square root of `det'(P†P) = L²`; the first version wrongly used `L²`. The one-loop measure's ghost L-power is fixed to `L¹` accordingly.

- **Which determinant**: Δ_FP = det'(P) = det'(P†P)^{1/2} = L  (NOT det'(P†P) = L²)
- **L-power fix**: ghost L-power L¹ (was wrongly L²); measure = ∫ (dL/L) det^{−1/2}_matter e^{−S}
- **CKV factor**: CKV volume L ⟹ 1/L; for L = 2π, 1/(2π) = PR #74 (unchanged)
- **Anomaly**: anomaly-free (no 1D conformal anomaly); only the discrete Z₂ (odd-k, PR #115) is nontrivial
- **Open**: absolute Z normalization (κ₅²/Λ₅); multi-loop; closed form

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_recap` | gauge sector of the measure (revised per review) | **PASS** |
| T2 | `T2_gauge_structure` | Diff(S¹) → constant einbein ⟹ 1 modulus L + 1 CKV | **PASS** |
| T3 | `T3_ghost_operator_and_kernel` | P = d/dτ; P†P = −d²/dτ²; ker(P) = constants = 1 CKV | **PASS** |
| T4 | `T4_which_determinant_is_the_ghost_det` | Δ_FP = det'(P) = det'(P†P)^{1/2} = L (NOT det'(P†P) = L²) | **PASS** |
| T5 | `T5_corrected_one_loop_measure` | measure Z = Σ ∫ (dL/L) det^{−1/2}_matter e^{−S}; ghost L¹ | **PASS** |
| T6 | `T6_ckv_is_closure_quantum_factor` | CKV 1/L; for L = 2π, 1/(2π) = PR #74 (unchanged) | **PASS** |
| T7 | `T7_anomaly_free_1d_worldline` | anomaly-free (1D; vs string c = −26) | **PASS** |
| T8 | `T8_assessment` | DIFF_S1_FP_GHOST_DET_IS_SQRT_DETPTP_EQUALS_L_ANOMALY_FREE | **PASS** |

## T4: which determinant is the ghost determinant

| candidate | value | role |
|---|---|---|
| `det'(P)` | L | the bc-ghost integral (FP determinant) |
| `det'(P†P)` | L² | intermediate Laplacian determinant (the **square** — not the ghost det) |
| `det'(P†P)^{1/2}` | L | = `det'(P)`; **THE Faddeev–Popov ghost determinant** |

| L | det′(P†P) | det′(P) | det′(P†P)^½ | match |
|---:|---:|---:|---:|:---:|
| 6.28319 | 39.47842 | 6.28319 | 6.28319 | ✓ |
| 1.0 | 1.0 | 1.0 | 1.0 | ✓ |
| 3.32408 | 11.04951 | 3.32408 | 3.32408 | ✓ |
| 5.0 | 25.0 | 5.0 | 5.0 | ✓ |

The Faddeev–Popov determinant is the `bc`-ghost path integral `∫ Db Dc e^{−∮ b(Pc)} = det'(P) = det'(P†P)^{1/2} = L`. The `L²` is `det'(P†P)` (the intermediate Laplacian determinant); the ghost determinant is its **square root**, `L`. The first version displayed `L²` — off by one square root.

## The corrected one-loop measure

```
Z = Σ_sectors  ∫ (dL/L)  ·  det^{−1/2}_matter  ·  e^{−S}
              └ Δ_FP = det′(P†P)^{1/2} = L is the einbein→proper-length
                Jacobian (⟹ modulus measure dL); CKV vol L ⟹ 1/L ┘
```

Ghost L-power = **L¹** (the proper-length Jacobian), not the spurious **L²** of the first version. The `1/L` is the CKV factor — `1/(2π)` at the closure loop `L = 2π` (PR #74).

## Verdict

**DIFF_S1_FP_GHOST_DET_IS_SQRT_DETPTP_EQUALS_L_ANOMALY_FREE.** THE Diff(S¹) FADDEEV–POPOV GHOST DETERMINANT IS det'(P†P)^{1/2} = L (NOT L²); THE ONE-LOOP MEASURE'S GHOST L-POWER IS FIXED ACCORDINGLY. (Revised per review.) PR #115 named the reparametrization Faddeev–Popov determinant as a piece of the measure; PR #116 regularized the matter determinant. This probe supplies the gauge sector — and the review correctly flagged that the first version conflated the ghost determinant with det'(P†P).

WHICH DETERMINANT. The Faddeev–Popov determinant is the bc-ghost path integral Δ_FP = ∫ Db Dc e^{−∮ b(Pc)} = det'(P), with P = d/dτ the operator taking the vector ghost c to the einbein variation. Because the b- and c-ghost spaces have equal dimension here, det'(P) = det'(P†P)^{1/2}. The intermediate Laplacian determinant is det'(P†P) = det'(−d²/dτ²) = L² (ζ(0) = −1, ζ'(0) = −2 ln L), so the GHOST determinant is its square root: Δ_FP = det'(P) = det'(P†P)^{1/2} = L. Of the three candidates {det'(P), det'(P†P), det'(P†P)^{1/2}}, it is det'(P†P)^{1/2} (= det'(P) in 1D by the ±n mode pairing) = L. The first version's "det'(−d²/dτ²) = L²" was the ghost determinant off by one square root.

THE CORRECTED MEASURE (L-POWER FIXED). Δ_FP = L is the Jacobian of the einbein → proper-length gauge fixing — it makes the modulus measure the proper circumference dL. Dividing by the conformal Killing volume Vol(U(1)) = L (the rigid rotation) gives the symmetry factor 1/L. So the gauge sector yields the standard worldline proper-time measure Z = Σ_sectors ∫ (dL/L) det^{−1/2}_matter e^{−S}: the ghost L-power is L¹ (the dL Jacobian), and there is NO separate det'_ghost = L² factor (the v1 error).

PR #74 UNCHANGED. The 1/(2π) is the CKV (c-ghost zero-mode) factor 1/Vol(U(1)) = 1/L at the closure loop L = 2π — independent of the determinant power, so the correction leaves it intact: PR #74's 1/(2π) is still the ghost zero-mode of the Diff(S¹) quotient.

ANOMALY-FREE. The 1D worldline carries no Weyl/conformal anomaly (no traceless symmetric 2-tensor in 1D), unlike the 2D string (bc-ghost c = −26 ⟹ D = 26); the only nontrivial anomaly is the discrete Z₂ orientation (odd-k, PR #115). With PR #116's matter determinant, the one-loop measure is finite/computable; the absolute normalization (κ₅²/Λ₅) and the multi-loop measure remain open.

## What this completes (and does not)

- **Completes:** the Diff(S¹) gauge sector — the FP ghost determinant `Δ_FP = det'(P†P)^{1/2} = L`, the corrected measure `Z = Σ ∫ (dL/L) det^{−1/2}_matter e^{−S}` (ghost L-power `L¹`), anomaly-free, with the CKV `1/L = 1/(2π)` (PR #74) intact. With PR #116's matter determinant the one-loop measure is finite/computable.
- **Open:** the absolute normalization of `Z` (the `κ₅²/Λ₅` bulk anchor, PR #112), the multi-loop / interacting measure, and a closed-form expression.
