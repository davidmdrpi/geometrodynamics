# The 5D Tangherlini bulk lift (PR #127)

PR #116 regularized the **Tangherlini fluctuation determinant** — the
one-loop matter background of the BAM measure — by running the radial cavity
operator

```
H = −d²/dr*² + V_tangherlini,   V = f(r)[l(l+2)/r² + 3 rs²/r⁴],  f(r) = 1 − (rs/r)².
```

That cavity is a *reduced, radial* object. This PR **lifts** it to its
explicit 5D parent — the Schwarzschild–Tangherlini (D=5) bulk metric — and
verifies that the BAM throat is the boundary trace of a **genuine 5D vacuum
geometry**, not a 4D ansatz dressed up. It then reconciles that
asymptotically-flat bulk with the AdS₅ / Randall–Sundrum brane bulk of PR #57.

All curvature in the probe is computed by a **self-contained numerical GR
routine** (metric → Christoffel → Riemann → Ricci/Kretschmann by high-order
finite differences); no symbolic-algebra dependency.

## The 5D Tangherlini metric

```
ds² = −f(r) dt² + f(r)⁻¹ dr² + r² dΩ₃²,   f(r) = 1 − (rs/r)^{D−3} = 1 − (rs/r)²   (D = 5),
```

with `dΩ₃²` the round unit 3-sphere and `rs = R_MID = 1`. The verified facts:

| property | value | meaning |
|---|---|---|
| metric power `D−3` | `2` | the **squared** Tangherlini power (D=5; matches spin-½ `T²=−I`, PR #73) |
| horizon | `f(rs)=0` at `r=rs=R_MID` | the **throat is the 5D horizon** |
| horizon topology | round `S³` | exactly BAM's **Hopf base** `S¹→S³→S²` (#59–#66) |
| Ricci | `R_μν = 0`, `Λ = 0` | a **genuine vacuum Einstein solution**, asymptotically flat |
| Kretschmann | `K = 72 rs⁴/r⁸` | finite on the cavity; **singularity only at `r=0`** |
| surface gravity | `κ = f'(rs)/2 = 1/rs` | — |
| Hawking temperature | `T_H = κ/2π = 1/(2π rs)` | the **closure quantum 2π** is the thermal period |

### Ricci-flat vacuum

`R_μν = 0` is confirmed numerically across the cavity (the residual is
finite-difference noise, `~10⁻⁴` near the horizon, `~10⁻⁷` away from it). The
throat's parent bulk is a real asymptotically-flat vacuum — **distinct** from
the AdS₅ RS bulk (`Λ₅ = −6k² < 0`). This is the honest tension the lift must
resolve (see below).

### Cavity curvature-regular

`K = 72 rs⁴/r⁸` (numeric matches analytic to `10⁻³`) is finite on the whole
cavity `[R_MID, R_OUTER]` — `72` at the throat down to `≈ 11.3` at
`R_OUTER = 1.26`. The **only** true curvature singularity is at `r = 0`,
behind the throat; `r = rs` is a *coordinate* (horizon) singularity, not a
curvature one. The ε healing length (PR #112) keeps the matter modes off the
horizon, in the regular exterior.

## The cavity potential descends from D=5

Both structural coefficients of the PR #116 potential are exactly the D=5
reductions of this metric:

  - **centrifugal `l(l+2)`** = the S³ scalar-Laplacian Casimir `l(l + D − 3)`
    at `D = 5` (the offset `D − 3 = 2` is also the metric power in `f`);
  - **curvature term `3 rs²/r⁴`** = `(D−2)/(2r)·f'(r)` at `D = 5`, with
    `f' = 2 rs²/r³` ⟹ coefficient `D − 2 = 3`.

So the program's "5" — `k₅ = D_bulk` (PR #73) — is realised here as the
**genuine bulk dimension of the metric**, not a fitted label. The throat is
the boundary of a real D=5 vacuum.

## Reconciliation with the AdS₅ / RS bulk

PR #57 tunes the BAM brane against an **AdS₅** bulk (`Λ₅ = −6k²`, the `√6`
Randall–Sundrum fine-tuning `λ_crit = √(6|Λ₅|)/κ₅²`). That bulk and the
Ricci-flat Tangherlini bulk are reconciled by the
**Schwarzschild–Tangherlini–AdS₅** metric

```
f(r) = 1 − rs²/r² + k² r²,
```

verified here to be **Einstein** with `R_μν = −4k² g_μν`, `Λ₅ = −6k²`. It
**interpolates**:

  - **near the throat** the AdS term `k²r² → 0` and `f → 1 − rs²/r²`, the pure
    Tangherlini neck (PR #116);
  - **far away** `f → k²r²`, the AdS₅ / RS asymptote (PR #57).

On the BAM cavity the AdS correction `k²r²` is `O(10⁻²)` for `k·rs ≲ 0.1`, so
the **pure-Tangherlini cavity of PR #116 is the correct near-throat limit**,
good to ~1%. The exact AdS scale `k` is the unpinned bulk ratio `κ₅²/Λ₅`
(PR #112) — the known open residual, organized here, not removed.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | lift the PR #116 cavity to its 5D bulk; reconcile with AdS₅/RS |
| T2 | 5D metric & horizon | `f = 1 − (rs/r)²`; throat = `S³` horizon (Hopf base) |
| T3 | Ricci-flat vacuum | `R_μν = 0`, `Λ = 0`, asymptotically flat |
| T4 | Kretschmann regularity | `K = 72 rs⁴/r⁸` regular on cavity; sing. only at `r=0` |
| T5 | potential descends | `l(l+2)` = S³ Casimir (`D−3=2`), `3rs²/r⁴` coeff `D−2=3` ⟹ `k₅=5` |
| T6 | Hawking period | `κ = 1/rs`, `T_H = 1/(2π rs)` carries the closure quantum |
| T7 | AdS₅/RS reconciliation | `f = 1 − rs²/r² + k²r²` Einstein (`Λ₅=−6k²`); interpolates; `k` open |
| T8 | assessment | `FIVE_D_TANGHERLINI_BULK_LIFT_RICCI_FLAT_VACUUM_CAVITY_REGULAR` |

## Established and open

  - **Established (BAM-native):** the BAM throat lifts to a genuine D=5
    Tangherlini vacuum bulk — Ricci-flat (`Λ=0`), Kretschmann `K = 72 rs⁴/r⁸`
    regular on the cavity (singularity only at `r=0`), `S³` horizon = the Hopf
    base, `T_H = 1/(2π rs)` carrying the closure quantum; the cavity
    potential's coefficients are the D=5 reductions (`k₅ = D_bulk = 5`); and it
    reconciles with the AdS₅/RS bulk via the Schwarzschild–Tangherlini–AdS₅
    interpolation.

  - **Does not / open:** the exact AdS scale `k` (= the unpinned `κ₅²/Λ₅`,
    PR #112) is not fixed; the full **global brane-localised**
    black-hole-on-brane solution (the junction-matched interpolating metric,
    not just the local neck and the asymptote) is not constructed; bulk
    backreaction / quantum gravity beyond the classical metric is not
    addressed.

## Cross-references

  - `docs/tangherlini_fluctuation_determinant_research_plan.md` — PR #116, the
    cavity operator lifted here.
  - `docs/k5_origin_research_plan.md` — PR #73, `k₅ = D_bulk = dim(S³)+2 = 5`.
  - `docs/brane_tension_tuning_research_plan.md` — PR #57, the AdS₅ / RS bulk
    and the `√6` tuning reconciled here.
  - `docs/epsilon_bulk_compliance_research_plan.md` — PR #112, the `κ₅²/Λ₅`
    bulk-ratio residual (the exact AdS scale `k`).

## Run

```
python -m experiments.closure_ledger.five_d_tangherlini_bulk_lift_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_five_d_tangherlini_bulk_lift_probe/`.
Expected verdict:
`FIVE_D_TANGHERLINI_BULK_LIFT_RICCI_FLAT_VACUUM_CAVITY_REGULAR`, 8/8 PASS.
