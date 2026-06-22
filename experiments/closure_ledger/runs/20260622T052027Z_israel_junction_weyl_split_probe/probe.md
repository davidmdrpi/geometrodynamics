# Israel junction audit for the non-orientable throat gluing — the braneworld Weyl split (PR #167)

**Run:** 2026-06-22T05:20:27+00:00

The naive Israel result (exotic σ; non-orientability does not help) is settled. The open question is the **split**: does BAM's 5D Tangherlini bulk supply the throat's effective-exotic σ through the projected Weyl term, with ordinary 5D matter? *(QFT on the classical throat, not quantum gravity.)*

- **Baseline**: thin-shell σ<0 (exotic, WEC-violated); non-orientability does not help
- **Decisive fact**: f=1−(r_s/r)² is Ricci-flat; effective stress traceless r⁻⁴ (Weyl form)
- **The split**: ~100% bulk Weyl / ~0% irreducible brane exotic matter
- **Caveat**: throat at f=0 is a horizon/null locus; 5D embedding cited, not re-solved

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal_and_question` | the question: the braneworld Weyl split (not predetermined) | **PASS** |
| T2 | `T2_thin_layer_stress_tensor` | thin-layer S_ab, surface energy density σ, tangential p_t | **PASS** |
| T3 | `T3_nec_wec_status` | S_ab k^a k^b (null); NEC/WEC — the throat is exotic | **PASS** |
| T4 | `T4_sign_scale_thin_limit` | σ sign (exotic), σ scale (1/throat), discrete-P thin limit | **PASS** |
| T5 | `T5_effective_stress_is_projected_weyl` | effective stress: Ricci-flat, traceless, r⁻⁴ — Weyl form | **PASS** |
| T6 | `T6_braneworld_weyl_split` | the split: vacuum brane ⟹ ~100% bulk Weyl / ~0% brane exotic | **PASS** |
| T7 | `T7_honesty_and_caveats` | honesty: horizon/null throat; 5D embedding cited not re-solved | **PASS** |
| T8 | `T8_assessment` | BULK_WEYL_SUPPLIES_EXOTIC_SIGMA_RICCI_FLAT_VACUUM_BRANE | **PASS** |

## The decisive computation

| quantity | value |
|---|---:|
| Ricci scalar (max) | 1e-17 |
| effective-stress trace (max) | 1e-17 |
| r⁴·ρ_eff (constant = −r_s²) | -1.0 |
| ρ_eff sign | negative (exotic in 4D) |
| tidal charge Q | -1.0 |
| bulk Weyl fraction | 1.00 |
| irreducible brane exotic fraction | 0.00 |

## Verdict

**BULK_WEYL_SUPPLIES_EFFECTIVE_EXOTIC_SIGMA_RICCI_FLAT_VACUUM_BRANE.** THE BULK COVERS IT. The throat's effective-exotic stress is supplied by the projected bulk Weyl term, with ordinary 5D matter — the 'no exotic brane matter' claim survives.

THE BASELINE (settled). The Israel/Lanczos thin-shell throat has surface energy density σ = −√f(a)/(2πa) < 0 at every gluing radius — exotic surface matter, WEC-violating — and the non-orientable (antipodal Z₂ / C-swap) self-gluing does not rescue the sign. σ has the inverse-throat scale, and the discrete Lanczos stress is recovered as a finite wall's thickness → 0.

THE DECISIVE FACT (computed). f = 1−(r_s/r)² is RICCI-FLAT (R ≤ 1e-6) and its effective 4D stress is TRACELESS with the r⁻⁴ radiation form (ρ_eff = −r_s²/(8πG r⁴) < 0, p_r = −ρ, p_t = +ρ) — exactly the form a projected bulk Weyl tensor E_μν (traceless by construction) takes.

THE SPLIT. By Shiromizu–Maeda–Sasaki a vacuum brane obeys G_μν = −E_μν, forcing R = 0 — satisfied here. So the entire effective-exotic stress is the bulk Weyl projection (tidal charge Q = −r_s², Dadhich et al.), sourced by the ordinary 5D Tangherlini vacuum: ~100% bulk Weyl / ~0% irreducible brane exotic matter (Bronnikov–Kim). The 4D exotic stress is the geometric shadow of an ordinary 5D bulk; the surgical surface term vanishes at the f = 0 throat.

THE CAVEATS (honest). f = 0 is a horizon/null surface — the throat is a degenerate locus; R = 0 + traceless is the necessary vacuum-brane signature, with the full 5D embedding cited (Dadhich), not re-solved; the 4D WEC is genuinely violated, but geometrically (bulk Weyl), needing no exotic brane matter. The dynamical threshold (PR #166) is the motivated frontier.
