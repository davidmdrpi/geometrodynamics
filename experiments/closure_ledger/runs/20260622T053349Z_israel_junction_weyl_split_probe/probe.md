# Israel junction audit for the non-orientable throat gluing — the braneworld Weyl split (PR #167)

**Run:** 2026-06-22T05:33:49+00:00

The naive Israel result (exotic σ; non-orientability does not help) is settled. The open question is the **split**: does BAM's 5D Tangherlini bulk supply the throat's effective-exotic σ through the projected Weyl term, with ordinary 5D matter? *(QFT on the classical throat, not quantum gravity.)*

- **Baseline**: thin-shell σ<0 (exotic, WEC-violated); non-orientability does not help
- **Decisive fact**: f=1−(r_s/r)² is Ricci-flat; effective stress traceless r⁻⁴, ρ<0 (tidal-charge/Weyl form)
- **The split**: 100% bulk-Weyl FORM; bulk-Weyl attribution consistent (necessary conditions met), 5D derivation pending
- **Caveat**: consistent-with, not proven; f=0 is a horizon/null locus (not evaded); needs no fundamental brane gauge field

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal_and_question` | the question: the braneworld Weyl split (not predetermined) | **PASS** |
| T2 | `T2_thin_layer_stress_tensor` | thin-layer S_ab, surface energy density σ, tangential p_t | **PASS** |
| T3 | `T3_nec_wec_status` | S_ab k^a k^b (null); NEC/WEC — the throat is exotic | **PASS** |
| T4 | `T4_sign_scale_thin_limit` | σ sign (exotic), σ scale (1/throat), discrete-P thin limit | **PASS** |
| T5 | `T5_effective_stress_is_projected_weyl` | effective stress: Ricci-flat, traceless, r⁻⁴ — Weyl form | **PASS** |
| T6 | `T6_braneworld_weyl_split` | the split: 100% Weyl FORM; attribution consistent, 5D pending | **PASS** |
| T7 | `T7_honesty_and_caveats` | honesty: horizon/null throat; 5D embedding cited not re-solved | **PASS** |
| T8 | `T8_assessment` | tidal-charge/bulk-Weyl form; necessary conditions met, 5D pending | **PASS** |

## The decisive computation

| quantity | value |
|---|---:|
| Ricci scalar (max) | 1e-17 |
| effective-stress trace (max) | 1e-17 |
| r⁴·ρ_eff (constant = −r_s²) | -1.0 |
| ρ_eff sign | negative (exotic in 4D; tidal charge Q = -1.0) |
| bulk-Weyl FORM fraction | 1.00 (computed) |
| Ricci-flat (necessary) | True |
| excludes real brane Maxwell (sign) | True |
| bulk-Weyl attribution proven | False (5D derivation pending (cited: Dadhich/Bronnikov–Kim)) |

## Verdict

**THROAT_IS_TIDAL_CHARGE_BULK_WEYL_FORM_NECESSARY_CONDITIONS_MET_5D_PENDING.** CONSISTENT-WITH, NOT PROVEN — the program's narrowest, most closable gap to date. BAM's throat is the tidal-charge metric whose 4D exotic stress is CONSISTENT WITH being entirely a bulk-Weyl projection of an ordinary 5D vacuum. This does NOT prove it is, and does NOT evade the f = 0 horizon.

THE BASELINE (settled). The Israel/Lanczos thin-shell throat has surface energy density σ = −√f(a)/(2πa) < 0 at every gluing radius — exotic surface matter, WEC-violating — and the non-orientable (antipodal Z₂ / C-swap) self-gluing does not rescue the sign. σ has the inverse-throat scale, and the discrete Lanczos stress is recovered as a finite wall's thickness → 0.

THE DECISIVE FACT (computed). f = 1−(r_s/r)² is RICCI-FLAT (R ≤ 1e-6) and its effective 4D stress is TRACELESS with the r⁻⁴ radiation form (ρ_eff = −r_s²/(8πG r⁴) < 0, p_r = −ρ, p_t = +ρ) — the tidal-charge / bulk-Weyl form. This is the form a projected bulk Weyl tensor takes; it is ALSO the form a real on-brane Maxwell field (Reissner–Nordström) takes, and only the SIGN distinguishes them.

THE SPLIT (honest). The effective stress is 100% of the bulk-Weyl FORM. The on-brane-exotic-free ATTRIBUTION is conditional: by Shiromizu–Maeda–Sasaki a vacuum brane obeys G_μν = −E_μν, forcing R = 0 (met); and the NEGATIVE tidal charge (ρ_eff < 0) excludes a real brane Maxwell source (which gives ρ > 0), so the reading also requires that BAM has no fundamental brane gauge field (a model assumption). The NECESSARY conditions are met; the SUFFICIENT step — the explicit 5D embedding whose Weyl projection sources exactly this E_μν — is PENDING (the Dadhich/Bronnikov–Kim construction is cited, not re-solved for BAM's bulk).

THE HONEST LINE. Throat stress is the tidal-charge / bulk-Weyl form; on-brane exotic matter is avoidable IF the 5D embedding sources E_μν and BAM has no fundamental brane gauge field — necessary conditions met, 5D derivation pending. And f = 0 is a horizon (null): the surgical surface term vanishes there, which relocates σ, it does not evade the horizon. This is the strongest result the audits have reached — narrow, specific, closable — and still a 'consistent-with'.
