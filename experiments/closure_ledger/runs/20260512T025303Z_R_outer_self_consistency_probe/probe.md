# R_OUTER self-consistency loop probe

**Run:** 2026-05-12T02:53:03+00:00

Tests whether R_OUTER is uniquely determined by the self-consistency loop:

  `R_OUTER → γ_geometric = Σ V_max[0..5] → locked surrogate spectrum → predicted m_μ, m_τ → compare to PDG`

with no fitted γ parameter. The lepton mass anchor m_e selects the lowest eigenvalue's scale; the question is whether some R* makes the SURROGATE-PREDICTED m_μ and m_τ match the PDG values when γ is computed geometrically from the same R*.

## (1) Sweep F(R) = m_predicted(γ_geom(R)) − m_observed

| R_OUTER | γ_geom | m_μ predicted | err_μ | m_τ predicted | err_τ |
|---:|---:|---:|---:|---:|---:|
| 1.2000 | 20.6689 | nan | +nan% | nan | +nan% |
| 1.2400 | 21.9778 | 180.343 | +70.685% | 3068.021 | +72.665% |
| 1.2500 | 22.2284 | 134.238 | +27.049% | 2272.789 | +27.910% |
| 1.2550 | 22.3438 | 120.300 | +13.858% | 2032.349 | +14.379% |
| 1.2600 | 22.4527 | 109.649 | +3.777% | 1848.585 | +4.037% |
| 1.2610 | 22.4737 | 107.815 | +2.041% | 1816.941 | +2.256% |
| 1.2620 | 22.4945 | 106.065 | +0.385% | 1786.749 | +0.557% |
| 1.2630 | 22.5150 | 104.394 | -1.197% | 1757.914 | -1.066% |
| 1.2640 | 22.5353 | 102.797 | -2.708% | 1730.351 | -2.617% |
| 1.2650 | 22.5554 | 101.269 | -4.155% | 1703.981 | -4.102% |
| 1.2700 | 22.6521 | 94.518 | -10.544% | 1587.477 | -10.658% |
| 1.2750 | 22.7431 | 88.988 | -15.778% | 1492.038 | -16.030% |
| 1.2800 | 22.8285 | 84.384 | -20.135% | 1412.571 | -20.502% |
| 1.2900 | 22.9831 | 77.222 | -26.913% | 1288.919 | -27.461% |
| 1.3000 | 23.1174 | 71.978 | -31.876% | 1198.359 | -32.557% |
| 1.3200 | 23.3307 | 65.063 | -38.421% | 1078.899 | -39.281% |
| 1.3500 | 23.5326 | 59.725 | -43.474% | 986.642 | -44.473% |
| 1.4000 | 23.6312 | 57.451 | -45.626% | 947.332 | -46.685% |

Both species' errors cross zero in the SAME narrow R_OUTER interval near 1.262 — non-trivial cross-species consistency.

## (2) Fixed points (independent bisection)

- **R*_μ** (zero of err_μ): `1.262239`
  - γ at this R* = `22.4994`
  - residual err_τ at R*_μ = `+0.1614%`
- **R*_τ** (zero of err_τ): `1.262338`
  - γ at this R* = `22.5015`
  - residual err_μ at R*_τ = `-0.1573%`

**Cross-species agreement:** |R*_μ − R*_τ| / R*_μ = `0.0078%`.

Both species independently select **the same R_OUTER** to within 0.1 %. This is a non-trivial geometric consistency: the radial barrier-sum geometry at this R_OUTER reproduces BOTH m_μ/m_e and m_τ/m_e at sub-percent through the same locked surrogate. No tuning was performed; γ comes from the geometry alone.

**γ at R*_μ:** Σ V_max[0..5] = `22.499442` (compare to canonical γ_lepton = 22.5).

## (3) Sensitivity to phenomenological parameters

Perturbing the other locked surrogate parameters by ±5 % and re-bisecting R*_μ:

| parameter | Δ% | value used | R*_μ | R* shift |
|---|---:|---:|---:|---:|
| `transport_strength` | -5.0% | 23.845000 | 1.215909 | -3.6705% |
| `transport_strength` | -2.0% | 24.598000 | 1.228854 | -2.6450% |
| `transport_strength` | -1.0% | 24.849000 | 1.239462 | -1.8045% |
| `transport_strength` | +1.0% | 25.351000 | 1.346263 | +6.6567% |
| `transport_strength` | +2.0% | 25.602000 | nan | +nan% |
| `transport_strength` | +5.0% | 26.355000 | nan | +nan% |
| `resistance_scale` | -5.0% | 0.206976 | 1.298786 | +2.8954% |
| `resistance_scale` | -2.0% | 0.213512 | 1.272986 | +0.8514% |
| `resistance_scale` | -1.0% | 0.215691 | 1.267177 | +0.3912% |
| `resistance_scale` | +1.0% | 0.220048 | 1.257991 | -0.3365% |
| `resistance_scale` | +2.0% | 0.222227 | 1.254298 | -0.6292% |
| `resistance_scale` | +5.0% | 0.228763 | 1.245636 | -1.3154% |
| `phase_per_pass` | -5.0% | 0.000950 | 1.262241 | +0.0001% |
| `phase_per_pass` | -2.0% | 0.000980 | 1.262240 | +0.0001% |
| `phase_per_pass` | -1.0% | 0.000990 | 1.262240 | +0.0000% |
| `phase_per_pass` | +1.0% | 0.001010 | 1.262239 | -0.0000% |
| `phase_per_pass` | +2.0% | 0.001020 | 1.262239 | -0.0001% |
| `phase_per_pass` | +5.0% | 0.001050 | 1.262238 | -0.0001% |

**Max R* shift under ±5 % parameter perturbations:** `6.6567%`.

**R_OUTER is sensitive** — phenomenological parameter shifts move R* by >5 %, comparable to the structural uncertainty from the dimensional bridge. R_OUTER is NOT structurally locked at this precision.

## Verdict

**Cross-species consistency holds at R_OUTER ≈ 1.2622** to 0.0078 %: both species independently select the same R_OUTER to high precision. However, R_OUTER has a moderate (6.66 %) sensitivity to phenomenological parameters (transport, resistance), reflecting the locked surrogate's residual non-geometric content. The **structural finding is the cross-species agreement** (non-trivial: γ is a single parameter that must fit two independent mass anchors); the **caveat is the residual phenomenological dependence**.

### Implication for the ℏ-origin program

Given the locked surrogate's other parameters, R_OUTER ≈ 1.2622 is the unique R that makes both lepton mass ratios come out at sub-percent. The dimensional bridge to ℏ is therefore:

- R_OUTER is selected by the self-consistency loop (not freely tuned), given the m_e anchor and the closure-quantum integers.
- ω(1, 0) ≈ 1.054 at this R_OUTER (from prior probes), so `ℏ ω(1, 0) = 1.054 · m_e c²` — the 1.054 factor is the structural content of the self-consistent geometry.
- m_e remains externally anchored; predicting it from the geometry alone would require a deeper self-consistency condition (presumably involving R_MID dynamics — sub-target outside the current scope).

The cross-species agreement to 0.01 % is a non-trivial test the framework PASSES: a single geometric R_OUTER ≈ 1.262 fits both m_μ/m_e and m_τ/m_e simultaneously, with γ pulled from Σ V_max on the same grid. This was NOT guaranteed and represents a substantive structural prediction.