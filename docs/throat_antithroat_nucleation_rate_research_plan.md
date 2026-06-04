# Throat–antithroat dynamical nucleation rate on the horizon-regular 5D background (PR #132)

The geometric throat arc (#127–#131) built the horizon-regular 5D Tangherlini
background and identified the antipodal map `(U,V) → (−U,−V)` as the
throat ↔ antithroat C-swap, but left the **dynamical nucleation rate** open
(the synthesis #131 listed it as the lead open item). The Majorana bounce arc
(#87–#90) had computed the bounce action `S` that controls `m_ν = m_D e^{−S}`
on the EM/tortoise picture. This PR puts the nucleation **on the horizon-regular
background** and shows what that geometry contributes that the earlier arc could
only posit.

## The nucleation rate

The throat ↔ antithroat transition (the `ΔL = 2` Majorana / pair-production
channel, PR #58) is the **region I ↔ III crossing** of the maximal Kruskal
extension (PR #128), mediated by the **odd** (`c₁ → −c₁`, the C-swap PR #63)
Euclidean instanton. Its rate has the standard bounce form

```
Γ ~ [det'(H)/det(H_free)]^{−1/2} · e^{−S},
```

with the one-loop prefactor the Tangherlini fluctuation determinant (PR #116)
and `S` the reduced Euclidean bounce action on the odd tortoise path.

## What the horizon-regular background contributes

### 1. The smooth Euclidean cigar (Gibbons–Hawking) — β = 2π rs

Wick-rotating `t → −iτ`, the near-horizon Euclidean metric in the proper radius
`ρ = √(2 rs (r − rs))` is `ds²_E ≈ dρ² + κ²ρ² dτ²` with `κ = f'(rs)/2 = 1/rs`
(PR #128). This is the flat plane in polar coordinates `(ρ, κτ)` and is smooth —
**no conical defect** — iff `κτ` has period `2π`, i.e. `β = 2π/κ = 2π rs`. The
deficit angle `2π − κβ` vanishes exactly there (and is `±1.26` at `0.8×`/`1.2×`).
So the Euclidean throat closes off smoothly, the **nucleation temperature is the
Hawking temperature** `T_nuc = 1/β = 1/(2π rs) = T_H`, and the period is the
**closure quantum 2π** (PR #127). (The near-horizon `f ≈ κ²ρ²` is verified:
ratio → 1 as `dr → 0`.)

### 2. The bounce action is the horizon tortoise divergence — S ∝ ln(1/ε)

The odd bounce path runs in to the throat; its tortoise length to the `ε`
healing length is `L*(ε) = |r*(R_OUTER) − r*(rs+ε)|` with the D=5 tortoise
`r*(r) = r + (rs/2) ln|(r−rs)/(r+rs)|`. As `ε → 0` this **diverges
logarithmically**, `L*(ε) → (rs/2) ln(1/ε) + const` — asymptotic slope `rs/2 =
0.5`, verified to 4 digits:

| ε | L*(ε) |
|---:|---:|
| 1e-2 | 1.820 |
| 1e-3 | 2.979 |
| 1e-4 | 4.130 |
| 1e-5 | 5.281 |
| 1e-6 | 6.433 |

The reduced action `S = (tension)·√(2μ E_c)·L*(ε)` therefore grows as
`ln(1/ε)`: the exact-horizon limit `ε → 0` costs infinite tortoise length ⟹
`S → ∞`, `Γ → 0`, `m_ν → 0` — the **"rigid throat ⟹ massless neutrino"** of #88,
now read off directly as the horizon tortoise divergence, regulated by the
finite `ε` healing length (PR #112).

### 3. The prefactor closes the arc — the #116 determinant

The one-loop nucleation prefactor is the Tangherlini fluctuation determinant of
PR #116, `det(H)/det(H_free) = 1.574370`. The geometric arc closes on itself:
**#116 is the bounce prefactor, #127/#128 the regular stage, #58/#87–#90 the
bounce.**

## The rate, with the inherited residuals

With the `ΔL = 2` / B−L tension window `t ∈ [2π, k₅√(2π)] ≈ [6.28, 12.53]`
(PR #89) and `ε ~ R_c³` (PR #112), the #88–#90 chain gives `S ≈ 15–18` and
`m_ν = m_D e^{−S} ~ few meV` — the observed scale, retrodicted to order of
magnitude. This PR places that chain on the regular background and supplies the
smoothness condition, the `ln(1/ε)` origin, and the prefactor; it does **not**
remove the inherited residuals.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | throat↔antithroat nucleation rate on the regular 5D background |
| T2 | smooth cigar | deficit `2π − κβ = 0` at `β = 2π/κ = 2π rs`; `T_nuc = T_H` |
| T3 | antipodal instanton | region I↔III (#128), `c₁→−c₁` (#63), ΔL=2 (#58) |
| T4 | log-ε tortoise | `L*(ε) = (rs/2) ln(1/ε) + const`, slope `rs/2` verified |
| T5 | action & rate | `t ∈ [2π, k₅√(2π)]` (#89), `ε ~ R_c³` ⟹ `S ≈ 15–18`, `m_ν ~ meV` |
| T6 | prefactor | the #116 Tangherlini determinant `1.574370` |
| T7 | scope | new vs inherited residuals |
| T8 | assessment | `THROAT_ANTITHROAT_NUCLEATION_REGULAR_CIGAR_LOG_EPSILON_BOUNCE` |

## Established and open

  - **Established (BAM-native, new here):** the throat ↔ antithroat nucleation
    on the horizon-regular background is the odd antipodal instanton on a smooth
    Euclidean cigar — `β = 2π rs = 2π/κ` (the closure quantum, Gibbons–Hawking
    smoothness, `T_nuc = T_H`); the bounce action `S ∝ ln(1/ε)` is the horizon
    tortoise divergence (slope `rs/2`, the "rigid throat ⟹ massless ν" read off
    geometrically); the rate prefactor is the #116 Tangherlini determinant
    (the arc closes on itself).

  - **Does not / open (inherited):** the exact `ε` value (PR #112), the absolute
    gravitational scale `κ₅²/Λ₅` (PR #112), and hence the precise `S` and `m_ν`
    (PRs #88–#90); the rate is order-of-magnitude (meV), not pinned.

## Cross-references

  - `docs/five_d_tangherlini_throat_horizon_lift_research_plan.md` — #128, the
    regular background + antipodal map.
  - `docs/majorana_bounce_action_research_plan.md` /
    `docs/b_minus_l_tension_ratio_research_plan.md` /
    `docs/boundary_compliance_bulk_geometry_research_plan.md` — #88/#89/#90, the
    bounce chain placed here.
  - `docs/tangherlini_fluctuation_determinant_research_plan.md` — #116, the
    prefactor.
  - `docs/seesaw_scale_nucleation_compliance_research_plan.md` — #87, `M_R =
    m_D e^{S}`, `m_ν = m_D e^{−S}`.

## Run

```
python -m experiments.closure_ledger.throat_antithroat_nucleation_rate_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_throat_antithroat_nucleation_rate_probe/`.
Expected verdict:
`THROAT_ANTITHROAT_NUCLEATION_REGULAR_CIGAR_LOG_EPSILON_BOUNCE`, 8/8 PASS.
