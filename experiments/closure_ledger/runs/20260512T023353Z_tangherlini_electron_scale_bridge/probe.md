# Tangherlini eigenvalue → electron-scale bridge probe

**Run:** 2026-05-12T02:33:53+00:00

Tests whether the canonical Compton identification `ℏ · ω(1, 0) = m_e c²` (which sets `R_MID = ℏ/(m_e c) = λ_C_reduced`) corresponds to a self-consistent value of R_OUTER, and how that compares with the locked γ_lepton = 22.5 pinhole geometry.

**Reference:** λ_C_reduced = ℏ / (m_e c) = 3.8616e-11 cm. Under the Compton bridge, geometric `R_MID = 1` corresponds to this length.

## (a) Compton bridge: R_OUTER such that ω(1, 0) = 1

- R_OUTER (geometric) = `1.449005`
- ω(l=1, n=0) = `1.000000` (target: 1 exactly)
- Σ V_max[l=0..5] at this R_OUTER = `23.6308`
- **Predicted R_MID** = ω · λ_C_reduced = `3.8616e-11` cm
- **Predicted R_OUTER (cm)** = `5.5955e-11` cm

## (b) γ_lepton lock: R_OUTER such that Σ V_max[0..5] = 22.5

- R_OUTER (geometric) = `1.262266`
- Σ V_max[l=0..5] = `22.5000` (target: 22.5)
- ω(l=1, n=0) at this R_OUTER = `1.053694` (deviation from Compton bridge: +5.369 %)
- Implied R_MID = ω · λ_C_reduced = `4.0689e-11` cm

## (c) Canonical baseline (R_OUTER = 1.26)

- ω(l=1, n=0) = `1.054727`
- Σ V_max[l=0..5] = `22.4527`
- Implied R_MID = `4.0729e-11` cm

## Structural tension between the two conditions

| condition | R_OUTER (geom) | ω(1, 0) | Σ V_max |
|---|---:|---:|---:|
| ω(1, 0) = 1 (Compton bridge) | 1.449005 | 1.000000 | 23.6308 |
| Σ V_max = 22.5 (γ_lepton lock) | 1.262266 | 1.053694 | 22.5000 |

**Relative difference:** `14.79 %`.

Two natural conditions on R_OUTER give DIFFERENT values: the Compton bridge (ω = 1 exactly) wants R_OUTER ≈ 1.4490, while the γ_lepton lock (Σ V_max[0..5] = 22.5) wants R_OUTER ≈ 1.2623. The difference is 14.79 %. Both cannot hold simultaneously: the current locked geometry chose the γ-lock, leaving a 5 % deviation from the Compton bridge.

## Verdict

**Substantial structural tension** (14.79 %). The two natural conditions on R_OUTER point to genuinely different geometries. The framework cannot satisfy both the Compton bridge AND the γ_lepton lock with a single R_OUTER under the canonical Tangherlini metric — one must be relaxed or refined.

**Bridge prediction at the Compton geometry:** under ω(1, 0) = 1, the throat radius is R_MID = λ_C_reduced ≈ 3.8616e-11 cm, exactly the reduced electron Compton wavelength. The outer cavity sits at R_OUTER = 1.4490 · λ_C_reduced ≈ 5.5955e-11 cm. **At this geometry, ℏ is fixed by m_e and c alone:** ℏ = m_e R_MID c, with R_MID predicted from ω(1, 0) = 1.

**Bridge prediction at the γ-locked geometry:** under Σ V_max = 22.5, R_OUTER = 1.2623 and ω(1, 0) = 1.0537. The implied R_MID ≈ 4.0689e-11 cm — about 5 % larger than λ_C_reduced. ℏ in physical units still requires m_e as anchor (the lock value is consistent with the lepton mass ladder, but doesn't predict ℏ independently).

## What's next

The probe identifies the **structural tension** between two natural R_OUTER conditions: the Compton bridge (ω = 1) and the γ_lepton lock (Σ V_max = 22.5). The next concrete sub-probes are:

- **Test the locked surrogate's mass ratios at the Compton bridge geometry.** Re-run `_build_generation_block` with R_OUTER set to the Compton-bridge value (where ω = 1 exactly) and see whether the lepton mass ratios m_μ/m_e and m_τ/m_e still come out at sub-percent. If yes: the Compton bridge is consistent with the locked spectrum and is the natural self-consistent geometry. If no: the γ-lock is the physical geometry and the 5 % Compton deviation is real.
- **Identify the structural meaning of R_OUTER ≈ 1.45.** The Compton-bridge R_OUTER is not obviously a clean geometric ratio (not √2, π/2, golden ratio, etc.). Either it's an emergent number with no closed form, or it's expressible as some combination of (k_5, π, …) in a way the current scan misses.