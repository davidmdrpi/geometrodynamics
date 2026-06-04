# Geometric throat arc synthesis (PR #131)

This is the **capstone** of the geometric throat arc — PRs #116 and #127–#130 —
which lifted the BAM throat from a radial cavity operator to its full
five-dimensional geometry and read off the physics. This PR introduces **no new
physics**; it (a) re-verifies, in one place, a keystone invariant from each arc
member (a cross-arc consistency check), (b) consolidates the unified picture —
the antipodal identification of the 5D Tangherlini horizon is the geometric
object **"Bulk Antipodal Mechanics"** is named for — and (c) lays out the honest
epistemic ledger.

## The arc

| PR | result |
|---|---|
| **#116** | Tangherlini fluctuation determinant — the matter cavity operator `H = −d²/dr*² + V_l`, `V_l = f[l(l+2)/r² + 3rs²/r⁴]`, `f = 1 − (rs/r)²` |
| **#127** | 5D Tangherlini bulk lift — the cavity is the boundary of a genuine D=5 vacuum: Ricci-flat (`Λ=0`), `K = 72 rs⁴/r⁸` regular on the cavity, `S³` horizon = Hopf base, `T_H = 1/(2π rs)`, `k₅ = D_bulk = 5`; reconciled with AdS₅/RS |
| **#128** | Horizon-regular lift — EF/Kruskal charts remove the coordinate singularity (`det g` finite, proper distance `√(2 rs ε)` = ε healing length, `F_Kruskal(rs) = 4 e⁻²`); the antipodal map `(U,V) → (−U,−V)` = the throat ↔ antithroat C-swap |
| **#129** | Null throat BC — the antipodal identification fixes the wave BC by l-parity (`Y_l(−x) = (−1)^l Y_l` ⟹ even-l Neumann, odd-l Dirichlet), a unitary mirror (zero throat flux) |
| **#130** | Antipodal vs absorbing spectrum — antipodal BC ⟹ real, undamped spectrum (stable matter); absorbing horizon ⟹ complex ringdown |

## The keystones, re-verified together

| PR | keystone | value |
|---|---|---|
| #116/#127 | Kretschmann at throat (regular) | `K = 72` |
| #116/#127 | Hawking temperature | `T_H = 0.159155 = 1/(2π rs)` |
| #128 | EF det at throat (nondegenerate) | `−0.299` |
| #128 | Kruskal factor at throat | `F(rs) = 0.541 = 4 e⁻²` |
| #128 | proper distance to throat | `√(2 rs ε) = 0.2` |
| #129 | antipodal BC (l=0 fundamental) | real `ω = 1.186` |
| #130 | absorbing BC (l=0 fundamental) | complex `ω = 1.893 − 1.159i` |

All keystones are mutually consistent and reproduce together in one run.

## One primitive, five faces

The antipodal identification of the 5D Tangherlini horizon is a single
geometric primitive that appears across the program as:

| face | PR | form |
|---|---|---|
| C = inner/outer swap | #63 | `c₁ → −c₁` |
| throat ↔ antithroat nucleation | #58 | nucleation channel |
| antipodal Kruskal map | #128 | `(U,V,Ω) → (−U,−V,Ω̄)` |
| l-parity unitary-mirror BC | #129 | even-l Neumann / odd-l Dirichlet |
| stable-matter selector | #130 | real undamped spectrum |

**"Bulk Antipodal Mechanics" is, literally, the mechanics of this one antipodal
identification on the bulk Tangherlini horizon.**

## The epistemic ledger

  - **Derived (within the arc):** the throat's parent bulk is a genuine D=5
    Tangherlini vacuum (Ricci-flat, curvature-regular cavity, `S³` horizon,
    `T_H = 1/2πrs`, `k₅ = D_bulk`); the coordinate singularity is removable
    (EF/Kruskal); the antipodal identification fixes the l-parity BC and makes
    the throat a unitary mirror; the antipodal spectrum is real (stable matter),
    the absorbing one complex (ringdown).
  - **Postulated (BAM's defining axiom):** the antipodal identification itself —
    that the throat is glued antipodally rather than absorbing. The arc shows
    this postulate is **self-consistent** (unitary, stable-matter-supporting),
    not that it is **forced**.
  - **Open:** the exact AdS scale `k = κ₅²/Λ₅` (PR #112); the dynamical
    throat ↔ antithroat nucleation rate (the bounce action, PRs #58/#88); the
    global brane-localised solution; the idealised `r* → −∞` horizon QNM tower
    and gravitational-radiation coupling (PR #130).

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | synthesise the geometric throat arc (#116, #127–#130) |
| T2 | throat = horizon | `f(rs)=0`, `K=72rs⁴/r⁸`, `T_H=1/2πrs`, `k₅=D_bulk=5` |
| T3 | regular & antipodal | EF det finite, `F(rs)=4e⁻²`, `(U,V)→(−U,−V)` = C-swap |
| T4 | antipodal BC | `Y_l(−x)=(−1)^l` ⟹ even-N/odd-D, unitary mirror |
| T5 | spectral consequence | antipodal real undamped vs absorbing complex ringdown |
| T6 | one primitive, five faces | #58/#63/#128/#129/#130 |
| T7 | epistemic ledger | derived / postulated / open |
| T8 | assessment | `GEOMETRIC_THROAT_ARC_SYNTHESIS_ANTIPODAL_HORIZON_UNIFIED` |

## Established and open

  - **Established (BAM-native):** the arc's keystones are mutually consistent
    and unify into one picture — the antipodal identification of the 5D
    Tangherlini horizon, derived as a genuine curvature-regular D=5 vacuum
    throat, postulated as glued antipodally, yielding a unitary mirror and a
    real, stable matter spectrum.

  - **Does not / open:** a synthesis/consistency capstone — it re-verifies and
    organises the arc; it does **not** add new derivations, remove any open
    item, or strengthen the antipodal postulate from "self-consistent" to
    "forced". The exact AdS scale, the nucleation rate, the global brane
    solution, and the idealised horizon QNMs remain open.

## Cross-references

  - `docs/tangherlini_fluctuation_determinant_research_plan.md` — #116
  - `docs/five_d_tangherlini_bulk_lift_research_plan.md` — #127
  - `docs/five_d_tangherlini_throat_horizon_lift_research_plan.md` — #128
  - `docs/null_throat_boundary_conditions_research_plan.md` — #129
  - `docs/antipodal_vs_absorbing_qnm_research_plan.md` — #130
  - `docs/charge_conjugation_swap_research_plan.md` — #63 (C-swap, one of the five faces)

## Run

```
python -m experiments.closure_ledger.geometric_throat_arc_synthesis_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_geometric_throat_arc_synthesis_probe/`.
Expected verdict:
`GEOMETRIC_THROAT_ARC_SYNTHESIS_ANTIPODAL_HORIZON_UNIFIED`, 8/8 PASS.
