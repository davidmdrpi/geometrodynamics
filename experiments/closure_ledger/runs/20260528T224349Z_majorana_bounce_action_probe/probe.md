# Majorana instanton-bounce action probe (PR #88)

**Run:** 2026-05-28T22:43:49+00:00

PR #87 reframed the seesaw scale as a bounce action `S = ln(m_D/m_ν) ≈ 15–18` but left `S` open. This probe builds a reduced Euclidean bounce for the `ΔL=2` throat↔antithroat flip along the **non-orientable tortoise path** and tests whether the throat tension + boundary compliance produce `S ≈ 15–18`.

- **Identification**: Majorana bounce = non-orientable tortoise logarithm: rigid throat ⟹ massless ν; S ∝ ln(1/ε), naturally O(10)/gen-stable; EM-throat tension under-produces S by ~40×
- **Mechanism**: reduced Euclidean bounce S = √(2 μ E_c)·L*(ε) along the odd tortoise path
- **Structural wins**: rigid throat = massless ν; S is a tortoise log
- **Open**: ΔL=2 (B−L) tension ratio t ≈ 6–12 + compliance ε (not yet derived)
- **B4 caveat**: tension/inertia from electron-throat calibration; S not produced by EM tension; closing needs a stiffer B−L sector

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_required_bounce_action` | required S = ln(m_D/m_ν) ≈ 14.7–17.6, gen-stable | **PASS** |
| T2 | `T2_non_orientable_tortoise_path` | bounce coord = odd tortoise path; r* log-diverges (slope rs/2) | **PASS** |
| T3 | `T3_rigid_throat_gives_massless_neutrino` | rigid throat ε→0 ⟹ S→∞ ⟹ massless ν (compliance = mass) | **PASS** |
| T4 | `T4_reduced_bounce_em_calibrated` | EM-throat bounce S ≲ 1 even near-rigid (~40× short) | **PASS** |
| T5 | `T5_required_delta_L_tension_ratio` | matching S∈[15,18] needs ΔL=2 tension ratio t ≈ 6–12 | **PASS** |
| T6 | `T6_generation_structure` | uniform S ⟹ m_ν∝m_D; ×18 spread needs gen-compliance/mixing | **PASS** |
| T7 | `T7_honest_scope` | structure established; open input → O(10) tension ratio | **PASS** |
| T8 | `T8_assessment` | MAJORANA_BOUNCE_IS_TORTOISE_LOG_S_OPEN_AS_TENSION_RATIO | **PASS** |

## T3: Rigid throat ⟹ massless neutrino

| compliance ε | tortoise length L* | bounce S | m_ν/m_D = e^{−S} |
|---:|---:|---:|---:|
| 1e-02 | 1.820 | 0.110 | 8.957e-01 |
| 1e-04 | 4.130 | 0.250 | 7.789e-01 |
| 1e-08 | 8.736 | 0.528 | 5.896e-01 |
| 1e-16 | 17.946 | 1.085 | 3.378e-01 |
| 1e-32 | 36.367 | 2.200 | 1.109e-01 |

As `ε → 0` (a rigid, incompressible throat) the tortoise length and the bounce diverge, so `m_ν → 0`: a perfectly rigid throat gives an **exactly massless** Majorana neutrino. Boundary compliance `ε > 0` is the mass-generating parameter.

## T4: Reduced bounce (EM-throat-calibrated)

- `σ = 2.6526e-02`, `ρ = 2.3873e-01`, `E_c = 5.4870e-03`, `μ = 0.3333`, `√(2μE_c) = 0.0605`

| compliance ε | tortoise length L* | bounce S |
|---:|---:|---:|
| 1e-01 | 0.601 | 0.0364 |
| 1e-02 | 1.820 | 0.1101 |
| 1e-04 | 4.130 | 0.2498 |
| 1e-08 | 8.736 | 0.5283 |
| 1e-16 | 17.946 | 1.0854 |

Even an absurdly rigid `ε = 1e-16` gives `S ≈ 1.085` — **~40× short** of the required ~16. The EM-throat tension under-produces the action.

## T5: Required ΔL=2 tension ratio

| compliance ε | tortoise L* | required tension ratio t | barrier boost t³ |
|---:|---:|---:|---:|
| 1e-02 | 1.820 | 12.05 | 1752× |
| 1e-03 | 2.978 | 9.42 | 837× |
| 1e-06 | 6.433 | 6.41 | 264× |

Matching `S ≈ 16` at a sane near-rigid compliance needs the `ΔL=2` (B−L) throat tension `t ≈ 6.4–12.1×` the EM-throat tension. Pure compliance (`t=1`) would need an absurd `ε ~ 6e-231`.

## Verdict

**MAJORANA_BOUNCE_IS_TORTOISE_LOG_S_OPEN_AS_TENSION_RATIO.** MAJORANA BOUNCE IS THE NON-ORIENTABLE TORTOISE LOGARITHM; THE OPEN INPUT IS NOW A ΔL=2 TENSION RATIO. PR #87 reframed the seesaw scale as a bounce action S = ln(m_D/m_ν) ≈ 15–18 but left S open. This probe builds a reduced Euclidean bounce for the throat↔antithroat flip and sharpens the verdict.

THE NON-ORIENTABLE TORTOISE PATH. The ΔL=2 Majorana flip is the throat↔antithroat odd extension (c₁→−c₁, PR #63 inner/outer swap; u_full = [−u[::-1], u]) — the orientation-reversing identification across the throat. The 5D tortoise coordinate r* = r + (rs/2)ln|(r−rs)/(r+rs)| diverges LOGARITHMICALLY at the throat (slope rs/2): the throat is infinitely deep in tortoise distance.

RIGID THROAT ⟹ MASSLESS NEUTRINO. The reduced bounce S = √(2 μ E_c)·L*(ε), with the PR #58 pinch barrier E_c and brane inertia μ, runs over the tortoise path length L*(ε) to within a compliance cutoff ε of the throat. As ε→0 (a rigid, incompressible throat) L*→∞ and S→∞, so m_ν = m_D·e^{−S}→0: a perfectly rigid throat gives an EXACTLY MASSLESS Majorana neutrino. A boundary compliance ε>0 is precisely what gives the neutrino mass — the smallness of m_ν is the near-rigidity of the throat.

S IS A TORTOISE LOGARITHM. Because L*(ε) = −(rs/2)ln ε + const, the bounce action is S ∝ ln(1/ε): naturally O(10) and slowly (logarithmically) varying — exactly the O(15), generation-stable form PR #87 required. The huge keV→TeV gap is carried by a logarithm, not a hierarchy.

THE MAGNITUDE FALLS SHORT (HONEST). With the throat tension fixed by the ELECTRON throat (PR #58/#87: σ=1/12π, ρ=3/4π), the under-barrier momentum √(2 μ E_c) ≈ 0.060, so even a near-rigid compliance gives S ≈ 0.04–1 — ~40× SHORT of the required ~16. Pure compliance cannot close the gap (it would need an absurd ε ~ 1e-230). Matching S∈[15,18] at a sane near-rigid compliance (ε ~ 1e-3–1e-6) requires the ΔL=2 (lepton-number-violating / B−L) throat tension to be a factor t ≈ 6–12 stiffer than the EM-throat tension (barrier ~10²–10³× larger, since E_c ∝ σ³).

PROGRESSIVE LOCALISATION. The open input has been sharpened three times: a mysterious ~TeV mass (PR #86) → an O(15) instanton action S (PR #87) → an O(10) dimensionless B−L/EM tension ratio (this probe). Each step ties the unknown more tightly to BAM geometry.

HONEST SCOPE. ESTABLISHED (BAM-native): the bounce coordinate is the non-orientable tortoise path; the rigid throat gives a massless neutrino (compliance is the mass-generating parameter); the bounce action is a tortoise logarithm, hence naturally O(10) and coarsely generation-stable. NOT established: that S≈15–18 follows from the EM-throat tension (it under-produces by ~40×); matching needs a ΔL=2 tension ~6–12× stiffer (a dimensionless ratio, not yet derived), the compliance ε (a sub-throat length), and the detailed m_ν spectrum (the ×18 m_ν/m_D spread needs gen-dependent compliance or the mixing sector).

## What this leaves open

- **The ΔL=2 (B−L) tension ratio** `t ≈ 6–12` — the open input is now this dimensionless number, not a `~TeV` mass. Progressive localisation: ~TeV mass (PR #86) → `O(15)` action `S` (PR #87) → `O(10)` tension ratio (this probe).
- **The boundary compliance `ε`** — a sub-throat length, not yet derived from the bulk.
- **The detailed `m_ν` spectrum** — the `×18` spread in `m_ν/m_D` needs a generation-dependent compliance or the mixing sector.
