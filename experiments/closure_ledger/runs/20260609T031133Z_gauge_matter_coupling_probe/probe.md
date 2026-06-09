# Gauge–matter coupling from the antipodal throat boundary (PR #141)

**Run:** 2026-06-09T03:11:33+00:00

Joins the gauge sector (the U(1)_Hopf photon, #42–#44) to the matter sector (the antipodal cavity modes, #129–#140) at the throat. The C-swap (#63) is inversion × charge conjugation (the throat is the C-surface), the gauge–matter vertex inherits the antipodal Z₂ selection rule, and U(1) charge is conserved by the unitary mirror — only α stays input. *(QFT on the fixed classical throat geometry, not quantum gravity.)*

- **Minimal coupling**: D_μ = ∂_μ − i c₁ A_μ; vertex c₁ ∫ A_μ j^μ
- **C-swap**: inversion (Y_l → (−1)^l) × charge conjugation (c₁ → −c₁); throat = C-surface (#63/#64)
- **Gauge vertex**: triple overlap ∫ Y_{l_γ}Y_{l₁}Y_{l₂}, Σl even (antipodal Z₂, #137/#140)
- **Charge conservation**: unitary mirror (#129) + C-swap flip ⟹ Σc₁ = 0 (#58)
- **Coupling strength**: α (the EM coupling, #105) — universal input residual
- **Open**: α (#105/#108); EM normalisation; higher gauge vertices; bulk scale (#133); flavor (#134)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | gauge–matter coupling from the antipodal throat boundary | **PASS** |
| T2 | `T2_minimal_gauge_coupling` | minimal coupling D_μ = ∂_μ − i c₁ A_μ; vertex c₁ ∫ A_μ j^μ | **PASS** |
| T3 | `T3_cswap_is_inversion_times_charge_conjugation` | C-swap = inversion (Y_l→(−1)^l) × charge conjugation (c₁→−c₁) | **PASS** |
| T4 | `T4_gauge_vertex_antipodal_z2_selection` | gauge vertex = Σl-even triple overlap (antipodal Z₂, #137/#140) | **PASS** |
| T5 | `T5_u1_charge_conserved_at_throat` | U(1) charge conserved at the throat: Σc₁ = 0 (#129/#58) | **PASS** |
| T6 | `T6_coupling_strength_is_alpha` | coupling strength = α (#105) — input residual | **PASS** |
| T7 | `T7_ledger_and_scope` | ledger/scope: structure derived; α / normalisation open | **PASS** |
| T8 | `T8_assessment` | GAUGE_MATTER_COUPLING_MINIMAL_Z2_SELECTED_CHARGE_CONSERVED_ALPHA_INPUT | **PASS** |

## The gauge–matter vertex obeys the antipodal Z₂ selection rule

| vertex (photon · matter·matter) | Σl | even? | ∫YYY | allowed? |
|---|---:|:---:|---:|:---:|
| γ(1) · matter(1,0) | 2 | ✓ | 0.25 | ✓ |
| γ(1) · matter(1,2) | 4 | ✓ | 0.041667 | ✓ |
| γ(0) · matter(1,1) | 2 | ✓ | 0.25 | ✓ |
| γ(1) · matter(1,1) | 3 | ✗ | 0.0 | ✗ |

The photon-matter-matter angular overlap is the cubic-vertex triple integral with one gauge leg; antipodal invariance forces `Σl = l_γ + l₁ + l₂` even — the same Z₂ as the matter vertices (#137/#140).

## Verdict

**GAUGE_MATTER_COUPLING_MINIMAL_Z2_SELECTED_CHARGE_CONSERVED_ALPHA_INPUT.** THE U(1)_HOPF GAUGE FIELD COUPLES MINIMALLY TO THE ANTIPODAL MATTER AT THE THROAT — THE C-SURFACE — WITH THE ANTIPODAL Z₂ SELECTION RULE AND CONSERVED CHARGE; ONLY α IS INPUT. The matter arc (#129–#140) built the matter sector on the antipodal throat and the gauge arc (#42–#44) the photon exchange kernel; this probe joins them at the throat boundary.

MINIMAL COUPLING. Matter of Hopf charge c₁ (|c₁| = 1, #58/#74) couples through D_μ = ∂_μ − i c₁ A_μ, giving the gauge–matter vertex c₁ ∫ A_μ j^μ with the matter current j^μ — the standard gauge-invariant coupling of the #42–#44 photon to the matter current.

THE C-SWAP = INVERSION × CHARGE CONJUGATION. The antipodal map A : x → −x (the C-swap #63) acts on both sectors at once: the matter harmonics carry Y_l → (−1)^l Y_l (#129/#140) and the Hopf charge flips c₁ → −c₁ (#63). One operation, two effects — a spatial inversion and a charge conjugation — so the throat is the particle ↔ antiparticle (C) surface (#63/#64), which is why the gauge field can couple to the matter there.

THE GAUGE VERTEX INHERITS THE ANTIPODAL Z₂ SELECTION RULE. The vertex couples a photon (l_γ) to two matter legs (l₁, l₂); its angular part ∫ Y_{l_γ}Y_{l₁}Y_{l₂} is the cubic-vertex triple overlap (#137) with one photon leg. S_BAM's antipodal invariance (the #140 Ward identity) forces Σl = l_γ + l₁ + l₂ even — the same antipodal Z₂ selection rule, now with the gauge leg (verified exactly: odd-Σl forbidden).

U(1) CHARGE IS CONSERVED AT THE THROAT. The antipodal mirror (#129) is unitary — zero net matter flux — and conserves the gauge charge flux; with the C-swap flip (c₁ → −c₁), outgoing charge re-emerges as the conjugate on the antipodal sheet, so Σc₁ = 0 (#58) and the throat balances particle against antiparticle. Charge conservation at the throat is the gauge face of the unitary mirror.

THE COUPLING STRENGTH IS α (INPUT). The minimal-coupling structure — the covariant derivative, the Σl-even triple-overlap vertex, the charge conservation, the C-surface — is derived from the antipodal geometry. The coupling strength is the EM coupling α (the "137 problem", #105), a universal residual not fixed by the geometry: structure derived, magnitude input.

SCOPE. Derives the gauge–matter coupling STRUCTURE at the antipodal throat. It does NOT fix α (#105) or the EM normalisation; higher gauge vertices and the running of α are not addressed. The α (#105/#108), bulk-scale (#133), and flavor (#134) residuals stand.
