# Horizon-regular coordinate lift for the 5D Tangherlini throat (PR #128)

PR #127 lifted the BAM throat to its explicit D=5 Schwarzschild–Tangherlini
bulk and showed the throat `r = rs = R_MID` is a **coordinate** (horizon)
singularity, not a curvature one — the Kretschmann scalar `K = 72 rs⁴/r⁸` is
finite there. But in Schwarzschild-type coordinates the metric still
degenerates at the throat (`g_rr = 1/f → ∞`). This PR constructs the
**horizon-regular coordinates** — Eddington–Finkelstein and Kruskal–Szekeres —
that remove the coordinate singularity, make the throat crossing manifestly
smooth, and exhibit the maximally-extended geometry whose **antipodal
identification is the geometric home of BAM's throat ↔ antithroat C-swap**.

All curvature is computed by the same self-contained numerical GR routine as
PR #127 (metric → Christoffel → Riemann → Ricci/Kretschmann), now applied to
the non-diagonal EF metric.

## The coordinate singularity is removable

In `(t, r)` coordinates `ds² = −f dt² + f⁻¹dr² + r²dΩ₃²`, `f = 1 − (rs/r)²`,
has `g_rr = 1/f → ∞` as `r → rs` (≈ 500 at `r = 1.001`). But `K = 72 rs⁴/r⁸`
is finite, so the divergence is a coordinate artifact, removable by a change of
chart.

## Eddington–Finkelstein: regular across the throat

With the D=5 tortoise coordinate `r*(r) = r + (rs/2) ln|(r−rs)/(r+rs)|` (from
`tangherlini/radial.py`) and the ingoing null time `v = t + r*`,

```
ds² = −f dv² + 2 dv dr + r² dΩ₃².
```

At the throat `g_vv = −f → 0` but `g_vr = 1`, so the `(v, r)` block has
determinant `−1` and the full metric determinant is `det g = −r⁶ sin⁴χ sin²θ`
— **finite and nonzero**. The Kretschmann scalar computed *in EF coordinates*
is still `72 rs⁴/r⁸` (verified by the numerical GR routine: `K = 72` at the
throat), confirming EF describes the same regular geometry with a nondegenerate
metric.

## Tortoise vs proper distance

| Δr = r − rs | tortoise `r*` | proper `∫dr/√f` | `√(2 rs Δr)` |
|---:|---:|---:|---:|
| 0.1 | −0.42 | 0.462 | 0.447 |
| 0.01 | −1.64 | 0.141 | 0.141 |
| 0.001 | −2.80 | 0.044 | 0.045 |

The tortoise (optical) distance to the throat is **infinite** (`r* → −∞`), but
the **proper** radial distance `∫dr/√f ≈ √(2 rs (r−rs))` is **finite** —
exactly the ε healing length `√(2 rs ε)` (PR #112). The throat is reachable in
finite proper distance even though it is infinitely far in the tortoise
coordinate.

## Kruskal–Szekeres: the maximal extension

The surface gravity is `κ = f'(rs)/2 = 1/rs`, so `κ·rs = 1`. With `u = t − r*`,
`v = t + r*` and `U = −(1/κ) e^{−κu}`, `V = (1/κ) e^{κv}`,

```
ds² = −F(r) dU dV + r² dΩ₃²,   F(r) = −f · e^{−2κr*} = (r+rs)²/r² · e^{−2r/rs}.
```

`F` is finite and nonzero at the throat — `F(rs) = 4 e⁻² ≈ 0.541` — **precisely
because `κ·rs = 1`** makes the `(r − rs)^{−κ rs}` factor cancel the simple zero
of `f`. The product `UV = −(1/κ²) e^{2κr*} → 0` at the throat: the **bifurcate
Killing horizon** `U = V = 0`. The full extension has four regions (exterior I,
interior II, antipodal exterior III, white hole IV).

## The antipodal identification = BAM's C-swap

The antipodal map `(U, V, Ω) → (−U, −V, Ω_antipodal)` is an **isometry** of the
extension that preserves `UV` (hence `r`) and maps region I ↔ region III. This
is the geometric realisation of BAM's throat ↔ antithroat identification — the
`C = inner/outer swap` (PR #63) with `c₁ → −c₁` (PR #58). The
maximally-extended 5D Tangherlini horizon with its antipodal gluing is the
geometric stage of **"Bulk Antipodal Mechanics"** itself: the program's
defining antipodal structure *is* the antipodal identification of the throat's
Kruskal horizon.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | construct horizon-regular coordinates (PR #127 coord. singularity) |
| T2 | singularity removable | `g_rr = 1/f → ∞` but `K = 72 rs⁴/r⁸` finite ⟹ coordinate |
| T3 | Eddington–Finkelstein | `det g = −r⁶ sin⁴χ sin²θ` finite/nonzero; `K` finite in EF |
| T4 | tortoise vs proper | `r* → −∞` (infinite); proper `√(2 rs Δr)` finite (= ε) |
| T5 | surface gravity & Kruskal | `κ = 1/rs` (`κ·rs = 1`); `F(rs) = 4 e⁻²`; `T_H = 1/(2π rs)` |
| T6 | maximal extension | `UV → 0` bifurcate Killing horizon; four regions |
| T7 | antipodal = C-swap | `(U,V) → (−U,−V)` isometry, I ↔ III = throat ↔ antithroat |
| T8 | assessment | `FIVE_D_TANGHERLINI_THROAT_HORIZON_REGULAR_ANTIPODAL_BIFURCATION` |

## Established and open

  - **Established (BAM-native):** the 5D Tangherlini throat's coordinate
    singularity is removable — Eddington–Finkelstein and Kruskal coordinates
    extend smoothly across it (nondegenerate metric, `K = 72 rs⁴/r⁸` finite,
    finite proper distance = the ε healing length), with `κ = 1/rs` and the
    Kruskal factor `F(rs) = 4 e⁻²`; the maximal extension's antipodal
    identification `(U,V) → (−U,−V)` is the geometric home of BAM's
    throat ↔ antithroat C-swap.

  - **Does not / open:** this is the classical, maximally-extended vacuum
    geometry — it does **not** compute the dynamical throat ↔ antithroat
    *nucleation amplitude* (the bounce rate, PR #58/#88); the lift provides that
    process's kinematic *stage* (the smooth crossing and the antipodal gluing),
    not its rate. The exact AdS scale `k` (= `κ₅²/Λ₅`, PR #112) and the global
    brane-localised solution (PR #127) remain open.

## Cross-references

  - `docs/five_d_tangherlini_bulk_lift_research_plan.md` — PR #127, the bulk
    lift whose coordinate horizon is regularised here.
  - `docs/charge_conjugation_swap_research_plan.md` — PR #63, `C = inner/outer
    swap` (the antipodal map realised geometrically here).
  - `docs/epsilon_bulk_compliance_research_plan.md` — PR #112, the ε healing
    length (= the finite proper distance to the throat).

## Run

```
python -m experiments.closure_ledger.five_d_tangherlini_throat_horizon_lift_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_five_d_tangherlini_throat_horizon_lift_probe/`.
Expected verdict:
`FIVE_D_TANGHERLINI_THROAT_HORIZON_REGULAR_ANTIPODAL_BIFURCATION`, 8/8 PASS.
