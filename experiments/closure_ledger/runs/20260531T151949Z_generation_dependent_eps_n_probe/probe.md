# Generation-dependent ε_n and the neutrino hierarchy spread (PR #113)

**Run:** 2026-05-31T15:19:49+00:00

Tests PR #91's proposed fix for the spread PR #112 left open: a generation-dependent healing length `ε_n` driven by the overtone boundary stress `χ_n` (decreasing with n). **Result: the direction is right, the magnitude overshoots.** `ε_n ∝ 1/χ_n` reproduces **normal ordering** untuned, but the natural `χ_n` variation (×8) is amplified by the steep bounce (`m_ν ∝ ε^{4.8}`, PR #112) into orders of magnitude in mass — `m_ν3/m_ν2 ≈ 162` vs observed 5.85 (×28). The spread stays a **residual**.

- **Direction**: normal ordering, untuned (ε_n ∝ 1/χ_n) — derived (sharpens PR #91)
- **Magnitude**: χ_n overshoots ×28; required ε_n gentle (1, 1.18, 1.57), not principled
- **Cause**: the steep bounce (m_ν ∝ ε^{4.8}, PR #112) amplifies the ×8 χ_n into ~10⁴ in mass
- **Status**: ε_n accommodates the spread but does not derive it — stays a residual
- **Plausibly**: the spread belongs to the mixing/anarchy sector (PR #92)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_setup` | uniform ε ⟹ ×2.7 spread; PR #91 proposed χ_n-driven ε_n | **PASS** |
| T2 | `T2_mechanism_and_direction` | ε_n ∝ 1/χ_n ⟹ normal ordering (direction right, untuned) | **PASS** |
| T3 | `T3_required_spread_is_gentle` | required ε_n gentle: (1, 1.18, 1.57) to hit observed m_2, m_3 | **PASS** |
| T4 | `T4_chi_driven_overshoots` | χ_n-driven overshoots: m_3/m_2 = 162 vs 5.85 (×28) | **PASS** |
| T5 | `T5_steep_bounce_amplifies_chi` | steep bounce (×4.8) amplifies ×8 χ_n; power p 0.15→0.31 (≠1) | **PASS** |
| T6 | `T6_accommodates_not_predicts` | no clean ε_n(χ_n) law — accommodates (fit), not predicts | **PASS** |
| T7 | `T7_honest_scope` | direction derived; spread residual (mixing/anarchy, PR #92) | **PASS** |
| T8 | `T8_assessment` | HIERARCHY_SPREAD_DIRECTION_FROM_EPSILON_N_MAGNITUDE_OVERSHOOTS_STAYS_RESIDUAL | **PASS** |

## Required (gentle) vs χ_n-driven (overshoots)

| gen | χ_n | required ε_n | required ratio | χ-driven ε_n ratio | χ-driven m_ν (meV) |
|---|---:|---:|---:|---:|---:|
| 1 | 0.304 | 0.011 | 1.0 | 1.0 | 2.1 |
| 2 | 0.097 | 0.013 | 1.18 | 3.13 | 1037.8 |
| 3 | 0.039 | 0.0172 | 1.57 | 7.79 | 167650.2 |

The required `ε_n` rises only ×1.57 over three generations; the principled `ε_n ∝ 1/χ_n` rises ×7.79, giving `m_ν3/m_ν2 = 161.5` vs observed `5.82` — a ×27.8 overshoot.

## Verdict

**HIERARCHY_SPREAD_DIRECTION_FROM_EPSILON_N_MAGNITUDE_OVERSHOOTS_STAYS_RESIDUAL.** A GENERATION-DEPENDENT ε_n GETS THE HIERARCHY DIRECTION RIGHT BUT OVERSHOOTS ITS MAGNITUDE — THE SPREAD STAYS A RESIDUAL. PR #112 derived the neutrino-mass SMALLNESS from a sub-throat, bulk-geometric healing length ε ~ R_c³, but a uniform ε gives a uniform bounce action, hence m_ν ∝ m_D — only a ×2.7 generation spread, far short of the observed hierarchy. PR #91 proposed the fix: the generations are cavity radial overtones n, and the overtone boundary stress χ_n (PR #79) decreases with n (0.304, 0.097, 0.039), so higher-overtone necks are more compliant — a generation-dependent ε_n. This probe makes that quantitative.

THE DIRECTION IS RIGHT. Compliance is the inverse of stiffness, so the natural law is ε_n ∝ 1/χ_n: with χ_n decreasing, ε_n increases with n, the tortoise length L*(ε_n) decreases, the action S(n) decreases, and the suppression weakens — so m_ν increases with n, and with m_D also increasing the result is NORMAL ORDERING, untuned. The qualitative hierarchy is reproduced.

THE MAGNITUDE OVERSHOOTS. But the observed hierarchy needs only a GENTLE ε_n variation. Anchoring gen 1 at ε_1 = R_c³ (the PR #112 value, m_ν1 ≈ 2.08 meV) and demanding the observed m_2 = 8.65, m_3 = 50.34 meV requires ε_n = (0.011, 0.013, 0.017), ratios (1, 1.18, 1.57). The χ_n-driven law badly overshoots: ε_n ∝ 1/χ_n gives ε ratios (1, 3.13, 7.79), hence m_ν = (2.1, 1038, 167650) meV and m_ν3/m_ν2 ≈ 162 against the observed 5.85 — a ×28 overshoot on the spread ratio, and orders of magnitude in absolute mass.

WHY: THE STEEP BOUNCE. The culprit is the bounce steepness established in PR #112, ∂ln m_ν/∂ln ε ≈ 4.8: the factor-~8 variation in χ_n is amplified into ~4 orders of magnitude in mass. The power in ε_n ∝ χ_n^{−p} that would reproduce the data is not the principled p = 1 but an inconsistent fractional p ≈ 0.15 (gen 1→2) to 0.31 (gen 2→3) — no single clean law fits both ratios. So a generation-dependent ε_n can ACCOMMODATE the spread (by fitting a gentle profile) but cannot PREDICT it from χ_n.

THE HONEST VERDICT. The hierarchy DIRECTION is derived (overtone compliance ⟹ normal ordering, untuned — a sharpening of PR #91), but its MAGNITUDE is not: the natural χ_n driver overshoots, and the same bounce steepness that made ε's absolute value a residual (PR #112) now blocks the natural overtone variation from setting the spread. The neutrino hierarchy spread stays a residual — ε_n accommodates it but does not derive it — and it plausibly belongs to the mixing / anarchy sector (PR #92) rather than a generation-dependent healing length.

## What this establishes (and does not)

- **Derived:** the hierarchy DIRECTION — overtone compliance (`ε_n ∝ 1/χ_n`) gives normal ordering, untuned (sharpens PR #91).
- **Not derived:** the hierarchy MAGNITUDE — the natural `χ_n` driver overshoots by ×28; the required `ε_n` is gentle (×1.57) and not a principled function of `χ_n`. The same bounce steepness that made `ε`'s value a residual (PR #112) blocks the spread, which plausibly belongs to the mixing/anarchy sector (PR #92).
