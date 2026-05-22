# B4 dimensional-bridge audit + Maslov closure-ledger probe

**Run:** 2026-05-22T20:35:47+00:00

Targets the last scaffold barrier B4 (the dimensional bridge ℏ = m_e·R_MID·c) and formalizes the Maslov index of the closure-ledger Bohr–Sommerfeld quantization. The Maslov machinery makes the ledger scale-free; scale-freeness is why B4 is irreducible.

## Maslov dictionary

- Bohr–Sommerfeld: `∮ p dr* = 2π(n + μ/4)`
- soft turning point: β = 0.25 (reflection phase π/2)
- hard wall (Dirichlet): β = 0.5 (reflection phase π)
- Tangherlini cavity: μ = 4 → radial ledger integer `n+1 (= μ/4 with μ=4)`

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_maslov_dictionary` | Maslov dict {soft:¼, hard:½} ✓ box & oscillator | **PASS** |
| T2 | `T2_tangherlini_doubly_dirichlet_mu4` | μ→4.00 (doubly-Dirichlet, →4) | **PASS** |
| T3 | `T3_throat_pi_is_three_things` | throat π = Dirichlet = closure-half = dwell | **PASS** |
| T4 | `T4_ledger_integers_dimensionless_quanta` | (3,6,109) = sums of dimensionless quanta | **PASS** |
| T5 | `T5_scale_invariance_core_audit` | scale-free (max dev ω·R 1e-11, S 5e-11) | **PASS** |
| T6 | `T6_dimensionless_residuals_closed` | R_OUTER, ε=7π/(100·5⁴), 8π, 7π/100 closed | **PASS** |
| T7 | `T7_1054_does_not_close_b4` | 1.054 dimensionless → orthogonal to B4 | **PASS** |
| T8 | `T8_b4_irreducibility` | 1 dimensionful input; ℏ=m_e·R_MID·c consistent | **PASS** |
| T9 | `T9_scaffold_final_assessment` | B1,B2,B3,B5 closed; B4 irreducible | **PASS** |

## T1: Maslov dictionary from analytic cavities

| cavity | n | ∫k dx/π | β_sum measured | β_sum predicted | match |
|---|---:|---:|---:|---:|:---:|
| square_well_hard+hard | 0 | 1.0000 | 1.0000 | 1.0000 | True |
| square_well_hard+hard | 1 | 2.0000 | 1.0000 | 1.0000 | True |
| square_well_hard+hard | 2 | 3.0000 | 1.0000 | 1.0000 | True |
| square_well_hard+hard | 3 | 4.0000 | 1.0000 | 1.0000 | True |
| harmonic_oscillator_soft+soft | 0 | 0.5000 | 0.5000 | 0.5000 | True |
| harmonic_oscillator_soft+soft | 1 | 1.5000 | 0.5000 | 0.5000 | True |
| harmonic_oscillator_soft+soft | 2 | 2.5000 | 0.5000 | 0.5000 | True |
| harmonic_oscillator_soft+soft | 3 | 3.5000 | 0.5000 | 0.5000 | True |

## T2: Tangherlini cavity is doubly-Dirichlet (μ=4)

| l | n | ω | WKB action | action/π | n+1 | μ implied | dev |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0 | 1.0547 | 2.7703 | 0.8818 | 1 | 3.527 | 0.118 |
| 1 | 1 | 1.9744 | 6.2655 | 1.9944 | 2 | 3.977 | 0.006 |
| 1 | 2 | 2.8941 | 9.4230 | 2.9994 | 3 | 3.998 | 0.001 |
| 1 | 3 | 3.8247 | 12.5659 | 3.9998 | 4 | 3.999 | 0.000 |
| 1 | 4 | 4.7609 | 15.7077 | 4.9999 | 5 | 4.000 | 0.000 |
| 1 | 5 | 5.7002 | 18.8494 | 6.0000 | 6 | 4.000 | 0.000 |

High-n implied Maslov index μ → 4.000 (target 4, doubly-Dirichlet).

## T3: The throat π is three things at once

- Dirichlet hard-wall reflection phase (B3): 3.141593
- closure half-quantum (2π / antipodal Z₂): 3.141593
- throat dwell phase ω·τ (τ=π/ω): [3.141592653589793, 3.141592653589793, 3.141592653589793, 3.141592653589793]
- throat Maslov contribution μ = 2 (β=½)

## T4: Ledger integers are sums of dimensionless quanta

| species | antipodal k | Hopf/throat | β-uplift | radial (n+1) | N_total |
|---|---:|---:|---:|---:|---:|
| electron | 1 | 1 | 0 | 1 | 3 |
| muon | 3 | 1 | 0 | 2 | 6 |
| tau | 5 | 1 | 100 | 3 | 109 |

## T5: Scale invariance — the core B4 audit

| λ | R_MID·λ | ω₀ (absolute) | ω₀·R_MID | WKB action S | ω₁/ω₀ |
|---:|---:|---:|---:|---:|---:|
| 1.00 | 1.0000 | 1.054725 | 1.054725 | 2.770585 | 1.871922 |
| 2.00 | 2.0000 | 0.527362 | 1.054725 | 2.770585 | 1.871922 |
| 0.50 | 0.5000 | 2.109449 | 1.054725 | 2.770585 | 1.871922 |
| 3.70 | 3.7000 | 0.285061 | 1.054725 | 2.770585 | 1.871922 |

- max |Δ(ω·R_MID)| = 1.05e-11
- max |ΔS| = 5.44e-11
- max |Δ(mass ratio)| = 1.14e-11
- **scale-free: True** — the machinery cannot produce a dimensionful scale.

## T6: Dimensionless residuals closed

| quantity | value | closed form | dimensionless |
|---|---:|---|:---:|
| R_OUTER | 1.2623 | cross-species fixed point γ_geom=ΣV_max=22.5 | True |
| epsilon_inner_cutoff | 0.000351858 | 7π/(100·5⁴) = 7π/62500 | True |
| transport | 25.1327 | 8π = 4·(2π) | True |
| resistance | 0.219911 | 7π/100 | True |

## T7: The 1.054 factor does not close B4

- ω(1,0) canonical = 1.054727, γ-lock = 1.053694 (dimensionless)
- closed form found in factor_1054_search: False (clean negative)
- closing 1.054 would close B4: False — a dimensionless number cannot supply the MeV scale.

## T8: B4 irreducibility

- bridge: `ℏ = m_e·R_MID·c`
- dim(m_e·R_MID·c) = (1, 2, -1) = dim(ℏ) = (1, 2, -1) (mass, length, time exponents)
- dimensionally consistent: True
- dimensionful inputs: 1 (m_e only)
- derived & dimensionless: closure quantum 2π, winding integers k,n, Maslov index μ, action ratios S/2π, ω·R_MID, mass ratios m_μ/m_e, m_τ/m_e, R_OUTER, ε, transport, resistance, 1.054

## T9: Scaffold final assessment

- **B1_closure_quantum**: CLOSED (winding θ-term)
- **B2_antipodal_Z2**: CLOSED (RP³ + spin structure)
- **B3_hard_wall**: CLOSED (T²=−I → Dirichlet; Maslov π reflection)
- **B5_reduction**: CLOSED (C×S³ master integral)
- **B4_dimensional_bridge**: IRREDUCIBLE (single m_e anchor, by scale-freeness)

## Verdict

**B4_IRREDUCIBLE.** B4 IRREDUCIBLE; MASLOV CLOSURE-LEDGER FORMALIZED. Two results in one probe:

PART A — the Maslov closure-ledger. Bohr–Sommerfeld with a Maslov correction at each boundary, ∮ p dr* = 2π(n + μ/4), with per-boundary deficits {soft turning point: ¼ (phase π/2); hard wall (Dirichlet): ½ (phase π)} — verified against the square well (hard+hard → n+1) and the harmonic oscillator (soft+soft → n+½). The Tangherlini throat cavity is bounded by a hard wall at each end (throat = Dirichlet from B3/T²=−I; outer wall), so μ=4 and ∮/2π → n+1: the closure-ledger Layer-2 "+1" per radial mode IS the Maslov index of the doubly-Dirichlet throat cavity. Its throat half (μ=2, reflection phase π) is the B3 hard wall — and that π is simultaneously the closure half-quantum (2π split by the antipodal Z₂) and the throat dwell phase ω·τ=π (τ=π/ω, the K-channel impedance). The per-species ledger integers (3, 6, 109) are sums of dimensionless quanta (antipodal k + Hopf/throat 1 + β-uplift + radial Maslov n+1).

PART B — the B4 audit. Every quantity the closure-ledger + Maslov machinery produces is dimensionless: winding integers, the Maslov index μ, action ratios S/2π, ω·R_MID, mass ratios, and the geometric residuals (R_OUTER, ε = 7π/(100·5⁴), transport = 8π, resistance = 7π/100, the 1.054 factor). Scale invariance is explicit: rescaling R_MID → λ·R_MID leaves every dimensionless output unchanged to machine precision and shifts only the absolute scale. A scale-free theory cannot produce a dimensionful scale, so exactly ONE dimensionful anchor is mathematically required; m_e (equivalently R_MID via ℏ = m_e·R_MID·c, dimensionally consistent) is the minimal choice. B4 is therefore IRREDUCIBLE by dimensional necessity — a structural feature, not an unsolved gap (cf. SI fixing c, ℏ, e by convention). The 1.054 factor is dimensionless and orthogonal to B4: even a closed form would not supply the MeV scale.

The scaffold is closed: B1, B2, B3, B5 derived; B4 shown to be the single mandatory dimensionful unit. What remains genuinely open is m_e itself — which, by this audit, cannot come from scale-free geometry.

## What this leaves open

- **m_e itself.** The single anchor is not derived — by this audit it *cannot* be, from scale-free geometry. Deriving m_e would require new dimensionful physics outside the closure-ledger scaffold.
- **First-principles internal action.** The Maslov indices are read from the WKB/boundary structure; an explicit covariant action whose stationary phase reproduces them is the standing follow-on (shared with the master-integral residual).
