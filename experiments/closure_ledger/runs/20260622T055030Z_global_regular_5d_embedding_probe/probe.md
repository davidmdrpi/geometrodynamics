# The global regular 5D embedding of the BAM throat (PR #168)

**Run:** 2026-06-22T05:50:30+00:00

Supplies the 5D derivation PR #167 flagged, as a **global regular** exact embedding (not Campbell–Magaard local existence): the BAM throat is the equatorial totally-geodesic slice of the 5D Schwarzschild–Tangherlini bulk. *(QFT on the classical throat, not quantum gravity.)*

- **Embedding**: equatorial χ=π/2 totally-geodesic slice of 5D Tangherlini (μ=r_s²)
- **Check 1**: induced metric = f=1−(r_s/r)²; K_μν=0 (totally geodesic)
- **Check 2**: projected bulk Weyl E_μν = −G⁴_μν (tidal fluid)
- **Check 3**: bulk Ricci-flat (ordinary 5D vacuum)
- **Regularity gate**: K₅=72μ²/ρ⁸ finite ρ≥r_s; singularity behind the regular horizon; compact χ — PASSES
- **Closes**: PR #167 positively; f=0 throat is the REGULAR 5D Killing horizon

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal_global_not_local` | goal: build the GLOBAL REGULAR embedding (not local existence) | **PASS** |
| T2 | `T2_construction_selects_bam_metric` | construction selects BAM's M=0 metric (μ=r_s²); the gate has teeth | **PASS** |
| T3 | `T3_check1_induced_metric_and_K` | check 1: induced metric = f; K_μν=0 (totally geodesic) | **PASS** |
| T4 | `T4_check2_projected_weyl_equals_tidal` | check 2: projected Weyl E_μν = −G⁴_μν (tidal fluid) | **PASS** |
| T5 | `T5_check3_bulk_field_equations` | check 3: bulk Ricci-flat (ordinary 5D vacuum) | **PASS** |
| T6 | `T6_regularity_gate` | regularity gate: K₅=72μ²/ρ⁸ finite ρ≥r_s — PASSES | **PASS** |
| T7 | `T7_what_closes_and_residue` | closes #167; f=0 throat is the REGULAR 5D Killing horizon | **PASS** |
| T8 | `T8_assessment` | GLOBAL_REGULAR_5D_EMBEDDING_EXISTS | **PASS** |

## The checks + the gate

| quantity | value |
|---|---:|
| induced metric = f=1−(r_s/r)² | yes; K_μν = 0.0 |
| projected Weyl error max\|E+G⁴\| | 1.1e-08 |
| bulk Ricci max\|R⁵_MN\| | 2.7e-07 |
| 5D Kretschmann at throat (72/r_s⁴) | 72.0 (finite) |
| regularity gate | PASSES |

## Verdict

**GLOBAL_REGULAR_5D_EMBEDDING_EXISTS_BAM_THROAT_IS_EQUATORIAL_TANGHERLINI_SLICE.** THE GAP CLOSES. The 5D derivation PR #167 flagged is supplied by an explicit, global, regular embedding — not local existence.

THE EMBEDDING. The BAM throat is the EQUATORIAL (χ = π/2) totally-geodesic slice of the 5D Schwarzschild–Tangherlini bulk ds²₅ = −F dt² + dρ²/F + ρ²dΩ₃², F = 1 − μ/ρ², with μ = r_s². The equator is a Z₂ fixed-point set, hence tension-free and matter-free; the construction works only for the pure-tidal (M = 0) form — exactly BAM's.

THE THREE CHECKS. (1) The induced metric is exactly f = 1−(r_s/r)² with K_μν = 0 (max |K| = 0e+00). (2) The projected bulk Weyl equals the 4D tidal stress, E_μν = −G⁴_μν (max |E + G⁴| = 1e-08) — the bulk-Weyl mechanism realised by an actual solution. (3) The bulk is Ricci-flat (max |R⁵_MN| = 3e-07) — an ordinary 5D vacuum.

THE REGULARITY GATE. K₅ = 72 μ²/ρ⁸ is finite throughout the exterior ρ ≥ r_s (= 72 at the throat); the only singularity ρ = 0 is behind the regular 5D Killing horizon; the extra dimension χ is compact and regular. The embedding is GLOBAL and REGULAR — the gate passes.

WHAT CLOSES. PR #167's bulk-Weyl reading is realised by an explicit regular 5D vacuum: no exotic brane matter, no brane gauge field, and the f = 0 throat is identified as the REGULAR 5D Killing horizon (the singularity ρ = 0 safely behind it). Honest residue: the throat sits at a (regular) horizon, the brane is the tension-free totally-geodesic slice (μ = r_s² fixed), and it is the exterior embedding.
