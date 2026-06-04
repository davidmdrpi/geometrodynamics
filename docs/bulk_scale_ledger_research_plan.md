# Bulk scale ledger for κ₅²/Λ₅ and ΔR normalization (PR #133)

The absolute bulk scale — the 5D gravitational coupling `κ₅²` and cosmological
constant `Λ₅` — has surfaced as an **open residual at every step** of the
program: the RS tuning (#57) fixed only the dimensionless `√6` and left
`k = √(|Λ₅|/6)` open; the ε healing length (#112) left its absolute
normalization to `κ₅²/Λ₅`; the bulk lift (#127) left the exact AdS scale `k`
open; the nucleation rate (#132) left the absolute scale open. This PR is the
consolidating **ledger**: it counts the bulk dimensionful content, separates the
scale **modulus** (the unit) from the genuine **residual**, and shows the
recurring "`κ₅²/Λ₅`" mystery reduces to **one bounded dimensionless number**.

## The bulk dimensionful content (D=5, ℏ = c = 1)

| parameter | dimension | meaning |
|---|---|---|
| `κ₅²` | `[L³]` | the 5D gravitational coupling (G₅) |
| `Λ₅` | `[L⁻²]` | the 5D cosmological constant ⟺ `k = √(|Λ₅|/6) [L⁻¹]`, `L_AdS = 1/k` |
| `λ_crit` | `[L⁻⁴]` | the 4D brane tension `= 6k/κ₅²` |
| `R_MID`, `ΔR` | `[L]` | throat radius, bulk separation `ΔR = R_OUTER − R_INNER = 0.52 R_MID` |

## Three categories, not one mystery

### 1. ΔR is the scale **modulus** — the unit, not a residual (#52/#53)

The B4 scale-modulus theorem (#52) proved BAM cannot derive an absolute unit
from scale-free topology: exactly **one** external dimensionful anchor is
required. `ΔR = R_OUTER − R_INNER = 0.52 R_MID` is that anchor — a proper,
cosmologically-invariant length (#53) — and it sets the unit (`R_MID = 1`). The
model-geometry ratios `ΔR/R_MID = 0.52`, `R_OUTER/R_MID = 1.26` are fixed. So
ΔR is **units**, not a residual.

### 2. √6 is the one **fixed** dimensionless tuning (#57)

The Randall–Sundrum flatness condition is `λ_crit κ₅²/√|Λ₅| = √6 ≈ 2.449` — one
condition among `(λ, Λ₅, κ₅)`, with the `6` the AdS₅ curvature coefficient
`Λ₅ = −6k²`. Dimensionless and derived.

### 3. The **open** bulk number: the AdS scale `k·R_MID`, bounded ≲ 0.1

Once the unit (ΔR) and the tuning (√6) are fixed, the only remaining
dimensionless bulk freedom is the AdS scale in throat units,

```
k · R_MID = R_MID / L_AdS  =  (κ₅²/Λ₅ in throat units).
```

This is **the** recurring residual. It is **not pinned**, but it is **bounded**:
the cavity correction to the pure-Tangherlini background is `(k r)²` (#127), so

| `k·R_MID` | cavity correction | `L_AdS/R_MID` |
|---:|---:|---:|
| 0.05 | 0.40% | 20 |
| 0.10 | 1.59% | 10 |
| 0.20 | 6.35% | 5 |

`k·R_MID ≲ 0.1` keeps it `≲ 1.6%`, so `R_MID ≲ L_AdS/10`: the throat sits **deep
in the near-flat region of the AdS bulk**, which is exactly why the
pure-Tangherlini cavity (#116/#127) is a good approximation. The 5D Newton
constant in throat units, `κ₅²/ΔR³`, sets the gravity strength — the
dimensionful anchor `G` (#105/#106).

## The consolidated ledger

```
{κ₅², Λ₅}  ⟶  { G  (gravity strength κ₅²/ΔR³ = the dimensionful anchor) }
              + { √6 (RS flatness tuning, FIXED, #57) }
              + { k·R_MID (AdS scale, OPEN but bounded ≲ 0.1, #127/#112) }
            with ΔR the unit (scale modulus, #52/#53).
```

So the "`κ₅²/Λ₅` mystery" is **one bounded dimensionless number** (`k·R_MID ≲
0.1`), not a multi-parameter freedom. The ledger does not pin it; it **bounds**
it and separates it cleanly from the unit and the tuning.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | bulk scale ledger for κ₅²/Λ₅ and ΔR (#57/#112/#127/#132) |
| T2 | dimensionful content | κ₅²[L³], Λ₅[L⁻²], k[L⁻¹], λ_crit[L⁻⁴], R_MID/ΔR[L] |
| T3 | ΔR = scale modulus | the unit (#52/#53), not a residual; ΔR = 0.52 R_MID |
| T4 | √6 fixed tuning | `λ_crit κ₅²/√|Λ₅| = √6 ≈ 2.449` (#57) |
| T5 | open ratio bounded | `k·R_MID ≲ 0.1` by the cavity correction `(k r)²` (#127) |
| T6 | consolidated ledger | `{κ₅²,Λ₅} → {G} + {√6} + {k·rs bounded} + {ΔR unit}` |
| T7 | scope | bounds the residual, does not pin it (= the #112 residual) |
| T8 | assessment | `BULK_SCALE_LEDGER_ONE_BOUNDED_ADS_RATIO_DELTAR_IS_SCALE_MODULUS` |

## Established and open

  - **Established (BAM-native):** the bulk dimensionful content `{κ₅², Λ₅}`
    reduces, once `ΔR` sets the unit (scale modulus, #52/#53) and `√6` fixes the
    RS tuning (#57), to **one** open but **bounded** dimensionless number — the
    AdS scale `k·R_MID ≲ 0.1` (#127), the recurring `κ₅²/Λ₅` residual. The
    ledger bounds and isolates it.

  - **Does not / open:** the ledger does **not** pin `k·R_MID` (still open, = the
    #112 residual), nor the absolute `G` normalization (the dimensionful anchor
    itself); it adds **no** new free parameter — it is the same absolute-scale
    residual, now shown to be singular and bounded.

## Cross-references

  - `docs/brane_tension_tuning_research_plan.md` — #57, the `√6` RS tuning,
    `Λ₅ = −6k²`, `λ_crit = 6k/κ₅²`.
  - `docs/delta_r_scale_modulus_research_plan.md` /
    `docs/maslov_dimensional_bridge_research_plan.md` — #53/#52, the B4
    scale-modulus theorem (one anchor; ΔR the proper invariant length).
  - `docs/epsilon_bulk_compliance_research_plan.md` — #112, the ε absolute
    normalization left to `κ₅²/Λ₅`.
  - `docs/five_d_tangherlini_bulk_lift_research_plan.md` — #127, the cavity
    correction `(k r)²` that bounds `k·R_MID`.
  - `docs/program_synthesis_research_plan.md` — #104, the input budget (G the
    one dimensionful anchor).

## Run

```
python -m experiments.closure_ledger.bulk_scale_ledger_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_bulk_scale_ledger_probe/`.
Expected verdict:
`BULK_SCALE_LEDGER_ONE_BOUNDED_ADS_RATIO_DELTAR_IS_SCALE_MODULUS`, 8/8 PASS.
