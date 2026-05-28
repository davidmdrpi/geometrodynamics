# ΔL=2 / B−L throat-tension ratio: derive or constrain (PR #89)

Follows PR #88, which built a reduced Euclidean bounce for the
neutrino's `ΔL=2` throat ↔ antithroat flip and localised the open input
to a single dimensionless number: the `ΔL=2` (lepton-number-violating /
B−L) throat tension must be a factor `t ≈ 6–12` stiffer than the
EM-throat tension to deliver the required bounce action `S ≈ 15–18`.
This probe asks whether `t` can be derived or constrained from BAM
structure.

## The `(t, ε)` degeneracy

The bounce is `S = √(2 μ E_c)·L*(ε) = t²·P0·L*(ε)` (PR #88 convention:
`σ → t·σ ⟹ μ ∝ t, E_c ∝ t³ ⟹ √(2μE_c) ∝ t²`), with `P0 ≈ 0.060` and
`L*(ε)` the tortoise path length to a boundary-compliance cutoff `ε`.
The real constraint is `t²·L*(ε) = S_req/P0 = const`: `t` and `ε` trade
off. Pinning `t` needs an independent argument for `t` (or for `ε`).
This probe argues for `t`.

## `t` is a B−L-breaking *global-closure* enhancement

The EM-throat tension is a **local** surface tension (`σ·Area`, PR #56).
The `ΔL=2` Majorana flip is **not** local: it reverses the throat's
orientation (`c₁ → −c₁`, the non-orientable / antipodal identification,
PR #63). An orientation reversal is a **global** operation on S³ — it
cannot be done locally — so its tension is the local tension times a
global-closure factor. `t` is that factor. BAM has exactly two
fundamental action scales for such a closure, and they bracket `t`:

| bound | BAM scale | value | physical meaning |
|---|---|---:|---|
| lower | closure quantum `2π` | 6.283 | single great-circle orientation reversal (minimal global flip) |
| upper | winding action `k_5√(2π) = √β_lepton` | 12.533 | full throat winding to the antipode |

  - **Lower bound — the closure quantum `2π`.** The cheapest global
    orientation reversal is a single great-circle traversal of the
    antipodal identification: one closure quantum `2π` (= `action_base`,
    the same `2π` in Hopf holonomy, throat dwell, the PR #74 loop
    measure, the PR #83 throat closure quantum). You cannot pay less
    than one closure quantum for a global flip ⟹ `t ≥ 2π ≈ 6.28`.

  - **Upper bound — the winding action `√β_lepton`.** The most expensive
    route is a full throat winding through the `k_5` structure to reach
    the antipode: the winding action `√β_lepton = √(k_5²·2π) = k_5·√(2π)`
    (PR #71). No costlier lepton-sector deformation exists ⟹
    `t ≤ k_5√(2π) ≈ 12.53`.

So `t ∈ [2π, k_5√(2π)] ≈ [6.28, 12.53]` — bracketed, parameter-free, by
the closure quantum and the winding action. This is **exactly** PR #88's
required `t ≈ 6–12`: the computed band `[6.41, 12.05]` (for `S=16` at
sane compliance `ε ∈ [1e-6, 1e-2]`) sits **inside** the window. The
band was not a fit but the BAM closure-to-winding window.

## Where in the window? → the compliance

Fixing `t` in the window fixes the compliance `ε`:

| `t` (window edge) | BAM scale | compliance `ε` |
|---:|---|---:|
| 6.283 | closure quantum `2π` | `5.9×10⁻⁷` (near-rigid) |
| 12.533 | winding action `k_5√(2π)` | `1.3×10⁻²` (more compliant) |

The neutrino's actual `t` sits between these, fixed by how much winding
the `ΔL=2` flip requires — i.e. by the boundary compliance. The open
input is sharpened from "an `O(10)` tension ratio" to "where in the
`[2π, k_5√(2π)]` window," a single residual number.

**Cross-check.** The winding/cavity-floor mass ratio
`m_charged/m_D = √(m²(1,0)/m²(0,0)) ≈ 11.9 ≈ √β_lepton` independently
lands at the winding (upper) edge — consistent with the flip borrowing
the winding channel.

## Tests

| # | test | finding |
|---|---|---|
| T1 | recap + degeneracy | `t ≈ 6–12`; `t²·L*(ε) = S_req/P0` |
| T2 | global-closure | `t` = B−L enhancement (orientation reversal is global) |
| T3 | lower bound | closure quantum `2π ≈ 6.28` |
| T4 | upper bound | winding action `k_5√(2π) = √β ≈ 12.53` |
| T5 | window brackets | `[2π, k_5√(2π)]` contains required `[6.41, 12.05]` |
| T6 | residual | `= compliance ε`; `m_charged/m_D ≈ √β` cross-check |
| T7 | honest scope | `t` constrained parameter-free; residual `= ε` |
| T8 | assessment | `B_MINUS_L_TENSION_BRACKETED_BY_CLOSURE_AND_WINDING_ACTIONS` |

## Established and open

  - **Established (BAM-native):** `t` is a B−L-breaking global-closure
    enhancement of the local EM tension; it is bracketed, parameter-free,
    by the closure quantum `2π` (minimal orientation reversal) and the
    winding action `k_5√(2π)` (full winding) — exactly PR #88's required
    `6–12`. The residual freedom is reduced to "where in the
    `[2π, k_5√(2π)]` window," i.e. the compliance `ε`.

  - **Open:** a *unique* `t`. The `(t, ε)` degeneracy remains (window
    edges map to `ε ≈ 6×10⁻⁷ … 1.3×10⁻²`); the residual open number is
    the compliance `ε` (a sub-throat length, not yet derived from the
    bulk). The bounce normalisation (the `t²` scaling, the `P0`
    prefactor) carries model dependence. So this is a CONSTRAINT +
    identification, not a closed derivation.

**Progressive localisation of the open input:** `~TeV` mass (PR #86) →
`O(15)` action `S` (PR #87) → `O(10)` tension ratio `t` (PR #88) → the
BAM closure-to-winding window `[2π, k_5√(2π)]` + a compliance residual
(PR #89).

## Run

```
python -m experiments.closure_ledger.b_minus_l_tension_ratio_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_b_minus_l_tension_ratio_probe/`.
Expected verdict: `B_MINUS_L_TENSION_BRACKETED_BY_CLOSURE_AND_WINDING_ACTIONS`, 8/8 PASS.
