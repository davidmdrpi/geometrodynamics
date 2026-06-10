# Electric and magnetic form factors: the EM gauge-arc capstone (PR #147)

**Run:** 2026-06-10T03:45:58+00:00

Assembles the relativistic vertex decomposition Γ^μ = γ^μF₁ + iσ^{μν}q_νF₂/2m on the antipodal cavity and capstones the EM gauge arc (#141–#146, with the #61/#62 magnetic keystones). The Gordon split is exact Dirac algebra; the Ward identity pins F₁ (exact charge) and leaves F₂ free (dressed moment) — one identity explains the arc's asymmetry; g = 2 (tree) and a = α/2π (loop) re-verify together; the radii are geometric. Only the value α(μ₀) is input. *(QFT on the classical throat, not quantum gravity.)*

- **Decomposition**: Γ^μ = γ^μF₁(q²) + iσ^{μν}q_νF₂(q²)/2m (Gordon-exact)
- **Ward asymmetry**: q_μσ^{μν}q_ν = 0 + ūq̸u = 0 ⟹ F₁ pinned, F₂ free
- **Keystones**: g_tree = 2 (#61); F₂(0) = α/2π = 0.00116141 (#62); vs a_e 0.00115965218
- **Sachs**: G_E(0) = c₁ exact; G_M(0) = c₁ + α/2π; r_M = r_E geometric (minimal model)
- **Arc**: #141→#147 one primitive: the unitary antipodal throat + Hopf charge
- **Open**: α² term; r_E−r_M splitting; recoil; normalisation (#133); flavor (#134)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | the EM-arc capstone: assemble the F₁/F₂ decomposition | **PASS** |
| T2 | `T2_gordon_decomposition` | Gordon identity exact (~1e-15): E/M split = Dirac algebra | **PASS** |
| T3 | `T3_ward_pins_f1_only` | Ward pins F₁ only: q_μσ^{μν}q_ν = 0, ūq̸u = 0 ⟹ charge exact, moment free | **PASS** |
| T4 | `T4_tree_keystone_g_equals_2` | tree keystone (σ·D)² = D² − σ·B ⟹ g = 2 (#61) | **PASS** |
| T5 | `T5_loop_keystone_f2_alpha_over_2pi` | loop keystone: simplex = 1 ⟹ a = α/2π vs a_e +0.15% (#62) | **PASS** |
| T6 | `T6_sachs_assembly_on_the_cavity` | Sachs: G_E(0) = c₁ exact, G_M(0) = c₁ + α/2π; r_M = r_E geometric | **PASS** |
| T7 | `T7_capstone_ledger` | arc ledger: #141→#147 one primitive; only α(μ₀) input | **PASS** |
| T8 | `T8_assessment` | EM_FORM_FACTORS_F1_WARD_PINNED_F2_ALPHA_OVER_2PI_RADII_GEOMETRIC_CAPSTONE | **PASS** |

## The exact identities (Dirac algebra)

| identity | max error |
|---|---:|
| Clifford {γ^μ,γ^ν} = 2η^μν | 0.0 |
| Gordon ū′γ^μu = ū′[(p+p′)^μ + iσq]u/2m | 1.33e-15 |
| Ward kills F₂: q_μσ^{μν}q_ν | 0.0 |
| on-shell ū′q̸u | 6.08e-16 |
| Pauli (σ·D)² = D² − σ·B | 8.75e-07 |

## The g ledger (tree + one loop vs experiment)

| quantity | value |
|---|---:|
| Schwinger simplex integral (expect 1) | 0.9999998 |
| a = F₂(0) = α/2π | 0.00116141 |
| a_e measured | 0.00115965218 |
| relative difference (α² and beyond) | 0.0015 |
| g (tree + α/2π) | 2.002322819 |
| g measured | 2.002319304 |

## The Sachs assembly (the arc's asymmetry, explicit)

| quantity | value |
|---|---:|
| G_E(0) (Ward-pinned, exact) | 1.0 |
| G_M(0) = c₁ + α/2π (dressed) | 1.0011614097 |
| g/2 = G_M(0)/G_E(0) | 1.00116141 |
| r_E (tortoise units, geometric, #146) | 0.2649 |
| r_M (minimal model) | 0.2649 |
| r_M/r_E − 1 | 0.0 |

| q | G_E(q)/G_E(0) | G_M(q)/G_M(0) |
|---:|---:|---:|
| 0.5 | 0.991257 | 0.991257 |
| 1.0 | 0.965393 | 0.965393 |
| 2.0 | 0.867278 | 0.867278 |
| 4.0 | 0.550363 | 0.550363 |

## Verdict

**EM_FORM_FACTORS_F1_WARD_PINNED_F2_ALPHA_OVER_2PI_RADII_GEOMETRIC_CAPSTONE.** THE ELECTRIC/MAGNETIC FORM-FACTOR DECOMPOSITION IS ASSEMBLED ON THE ANTIPODAL CAVITY AND THE EM GAUGE ARC IS CAPSTONED: THE GORDON SPLIT IS EXACT DIRAC ALGEBRA, THE WARD IDENTITY PINS F₁ (EXACT CHARGE) AND LEAVES F₂ FREE (DRESSED MOMENT), THE KEYSTONES g = 2 AND a = α/2π RE-VERIFY TOGETHER, AND THE RADII ARE GEOMETRIC — ONLY α(μ₀) IS INPUT.

THE GORDON DECOMPOSITION. ū′γ^μu = ū′[(p+p′)^μ + iσ^{μν}q_ν]u/2m verified with explicit Dirac spinors to ~1e-15: the EM current splits exactly into convection (electric) + spin (magnetic) parts. The F₁/F₂ decomposition is the Dirac algebra of the #141 minimal coupling, not an ansatz.

WHY THE CHARGE IS EXACT AND THE MOMENT IS NOT — ONE IDENTITY. The Ward contraction kills the F₂ term twice over: q_μσ^{μν}q_ν = 0 identically (exact) and ū′q̸u = 0 on shell (~1e-16). So the #142/#145 Ward identity constrains F₁ only: F₁(0) = c₁ exact and coupling-independent (#145/#146), while F₂ is gauge-free and dresses at every loop — the arc's asymmetry is Dirac algebra.

THE KEYSTONES, RE-VERIFIED TOGETHER (the #131 convention). TREE: (σ·D)² = D² − σ·B (~1e-6, finite differences) — the SU(2) anticommutator factor of 2 IS g_s = 2, F₂(0) = 0 at tree level (#61). LOOP: the Schwinger simplex integral = 1 (numerically 0.999998) ⟹ F₂(0) = α/2π = 0.00116141 vs the measured a_e = 0.00115965 (+0.15%, the α² Sommerfield term and beyond); g = 2(1 + α/2π) = 2.0023228 vs 2.0023193 (#62).

THE SACHS ASSEMBLY. G_E(q) is the #146 dressed-density transform with the geometric radius r_E = 0.265 (tortoise units); the magnetization rides the same charged-mode profile (the spin is on the charged line; φ is spinless), so r_M = r_E and G_M/G_M(0) = G_E/G_E(0) — form-factor scaling in the minimal model. At q = 0: G_E(0) = c₁ = 1 EXACT, G_M(0) = 1 + α/2π DRESSED, g/2 = G_M(0)/c₁.

THE ARC, ONE PRIMITIVE. #141 coupling, #142 Ward, #143 α ledger, #144 Π/running, #145 exact universal charge, #146 G_E/geometric radius, #147 the decomposition — every face derives from the unitary antipodal throat carrying the integer Hopf charge; the single EM input is the value α(μ₀) (#143).

SCOPE. Capstone assembly: the F₂(q) shape is modelled (flat ⟹ scaling), recoil/O(q²/m²) mixing unresolved; the α² Sommerfield term, the r_E − r_M splitting, the absolute normalisation (#133), and the flavor residuals (#134) stand.
