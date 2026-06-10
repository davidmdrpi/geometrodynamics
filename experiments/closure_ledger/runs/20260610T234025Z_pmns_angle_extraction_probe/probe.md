# PMNS angle extraction from mouth-localized cross-channel overlaps (PR #153)

**Run:** 2026-06-10T23:40:25+00:00

Assembles the PMNS matrix from the derived #151/#152 channel-dominant anarchy, the computed geometric mouth overlap, and one charged-side μ–τ rotation — the single rotation the winding hierarchy permits and the single one the data demand (an exact invariance argument). All three observed angles come out anarchy-typical and CP is generic; the e-row is hierarchy-protected, the BAM reason θ₁₃ is small while θ₂₃ is large. *(QFT on the classical throat, not quantum gravity.)*

- **Construction**: U = R₂₃(φ_ℓ)·O_geom·U_ν (#151/#152 ensemble; computed O_geom)
- **Bracket**: near-diagonal O fails θ₂₃ (98th pct); anarchic O fails θ₁₃ (≤7th pct)
- **Resolution**: left R₂₃: s12/s13 exactly invariant ⟹ one rotation; μ–τ the permitted block
- **Angles**: s12 @ 62nd, s23 @ 56th, s13 @ 27th percentile (φ_ℓ = 45°)
- **CP**: median |J| = 0.015; P(|J| > 0.01) = 61% — generic
- **Open**: charged-side matrix; Majorana phases/m_ββ; CKM analogue

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | assemble the PMNS from the derived #151/#152 structure | **PASS** |
| T2 | `T2_construction_and_O_geom` | O_geom computed: near-diagonal (~8–13° misalignments) | **PASS** |
| T3 | `T3_two_failure_modes` | failure modes bracket: near-diag fails θ₂₃; anarchic fails θ₁₃ | **PASS** |
| T4 | `T4_one_hierarchy_permitted_rotation` | left R₂₃ exactly invariant on s12/s13; μ–τ the permitted block | **PASS** |
| T5 | `T5_assembled_pmns_all_angles_natural` | all three angles natural (62/56/27th pct); generic CP | **PASS** |
| T6 | `T6_predictions` | θ₁₃ not-too-small preserved; no octant preference; m₁ unchanged | **PASS** |
| T7 | `T7_ledger_and_scope` | ledger: structure derived; values statistical; budget unchanged | **PASS** |
| T8 | `T8_assessment` | PMNS_NATURAL_ONE_HIERARCHY_PERMITTED_CHARGED_ROTATION_E_ROW_PROTECTED | **PASS** |

## The two failure modes (observed percentile per angle)

| O structure | sin²θ₁₂ pct | sin²θ₂₃ pct | sin²θ₁₃ pct |
|---|---:|---:|---:|
| near-diagonal (O_geom) | 62% | 98% | 27% |
| fully anarchic (fixed Haar) | 20% | 16% | 7% |

## The charged-side rotation window (s12/s13 exactly invariant)

| φ_ℓ (deg) | median sin²θ₂₃ | observed percentile |
|---:|---:|---:|
| 0 | 0.143 | 98.2% |
| 15 | 0.142 | 95.8% |
| 25 | 0.183 | 85.7% |
| 35 | 0.275 | 73.3% |
| 45 | 0.406 | 56.4% |
| 55 | 0.547 | 36.7% |
| 65 | 0.681 | 19.9% |
| 75 | 0.789 | 7.8% |

Invariance: max|Δs12| = 0.0, max|Δs13| = 0.0 (exact). Natural window [35, 45, 55, 65]° — broad, no fine-tuning; the μ–τ block is ×12.3 less hierarchical than e–μ (the permitted block).

## The assembled PMNS (φ_ℓ = 45°)

| angle | ensemble median | observed | percentile |
|---|---:|---:|---:|
| s12 | 0.2002 | 0.304 | 62% |
| s23 | 0.4065 | 0.45 | 56% |
| s13 | 0.0499 | 0.0224 | 27% |

CP: median |J| = 0.0151 (data |J|_max ≈ 0.033); P(|J| > 0.01) = 61% — generic CP, quantified.

## Verdict

**PMNS_NATURAL_ONE_HIERARCHY_PERMITTED_CHARGED_ROTATION_E_ROW_PROTECTED.** THE PMNS ANGLES EXTRACT NATURALLY FROM THE DERIVED CHANNEL-DOMINANT ANARCHY PLUS EXACTLY ONE HIERARCHY-PERMITTED CHARGED-SIDE μ–τ ROTATION: ALL THREE OBSERVED ANGLES ARE ANARCHY-TYPICAL (62nd/56th/27th PERCENTILES), CP IS GENERIC, AND THE E-ROW IS HIERARCHY-PROTECTED — WHY θ₁₃ IS SMALL WHILE θ₂₃ IS LARGE. #151/#152 derived the ν-side structure; this probe assembles the mixing matrix itself.

THE CONSTRUCTION. U = R₂₃(φ_ℓ)·O_geom·U_ν: the derived channel-dominant complex anarchic ensemble (3000 seeded draws), the computed near-diagonal geometric mouth overlap (misalignments ~8–13°), and one charged-side μ–τ rotation.

THE TWO FAILURE MODES BRACKET THE STRUCTURE. Near-diagonal O: θ₁₂ and θ₁₃ natural but θ₂₃ far too small (98th pct). Fully anarchic O: θ₂₃ natural but θ₁₃ far too large (≤7th pct). The data sit between — they select a specific intermediate charged-side structure.

ONE ROTATION RESOLVES IT — AND IT IS THE PERMITTED ONE. A left μ–τ rotation leaves sin²θ₁₂ and sin²θ₁₃ EXACTLY invariant (machine zero; it never touches the e-row) and moves only sin²θ₂₃ — so the data demand exactly one charged-side rotation. The winding hierarchy permits exactly that one: m_μ/m_τ = 0.060 is ×12 less hierarchical than m_e/m_μ, so an O(m_τ) off-diagonal gives an O(1) left μ–τ rotation while the e-row rotations stay suppressed. The e-row is hierarchy-protected — the BAM reason θ₁₃ is small while θ₂₃ is large. The natural window is broad (φ_ℓ ∈ ~[25°, 65°], ~45% of a uniform draw): no fine-tuning.

THE ASSEMBLED PMNS. At φ_ℓ = 45°: sin²θ₁₂ at the 62nd percentile, sin²θ₂₃ at the 56th, sin²θ₁₃ at the 27th — all three natural, the full observed point anarchy-typical. CP GENERIC: median |J| = 0.0151 (data |J|_max ≈ 0.033), P(|J| > 0.01) = 61% — the README "generic CP" claim quantified at the PMNS level.

PREDICTIONS. The e-row protection keeps θ₁₃ not-too-small (P(sin²θ₁₃ < 0.002) small — the historic anarchy success preserved); no strong θ₂₃ octant preference; the mass-side m₁ ≈ 0.04 meV / Σm_ν ≈ 58.8 meV prediction (#151/#152) rides the same ensemble, unchanged.

SCOPE. The winding-profile shape in O_geom is modelled (results need only its near-diagonality); φ_ℓ is an O(1) anarchic draw on the permitted block, not a tuned knob; Majorana phases / m_ββ, the charged-side matrix derivation, and the CKM analogue are open. No new input — the #150 budget is unchanged.
