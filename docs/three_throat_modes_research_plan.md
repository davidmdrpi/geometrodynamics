# Three throat modes from `k_5`: `#gen = (k_5+1)/2 = 3`

Closes the PR #70 follow-on "why exactly 3 throat modes." The clean
structural answer ties directly to the same `k_5 = 5` primitive used in
PR #71's `β_lepton = k_5²·(2π)`. The charged-lepton sector is the odd
subset of the `γ-lock` `l`-range `[0, k_5]`; the count is

```
#generations  =  (k_5 + 1) / 2  =  3   (for k_5 = 5) .
```

The cavity geometry independently shows the throat-localized ladder
saturates at ~3 modes (PR #68's saturating crossover) — consistent.

## The structural form

Combining the three established results:

  - **#67** (odd-k fermionic selection): charged leptons sit at **odd**
    closure depths `k`.
  - **γ-lock from `hbar_origin_status`**: the lepton sector's `l`-range
    is `[0, k_5]` (the cross-species fixed point uses
    `Σ V_max[0..k_5] = 22.5` with `k_5 = 5`).
  - **#71** (β_lepton structural form): the same `k_5` is the
    topological charge / Tangherlini bulk dimension.

The lepton depths are therefore the **odd integers in `[0, k_5]`**:

```
{1, 3, 5}   ←   odd integers in [0, 5]   →   (e, μ, τ) .
```

The count is `(k_5+1)/2`. For `k_5 = 5`: 3 charged-lepton generations.

## Dependence on `k_5`

| `k_5` | odd `l` in `[0, k_5]` | `#generations = (k_5+1)/2` |
|---:|---|---:|
| 3 | `{1, 3}` | 2 |
| **5** | **`{1, 3, 5}`** | **3** ← BAM |
| 7 | `{1, 3, 5, 7}` | 4 |
| 9 | `{1, 3, 5, 7, 9}` | 5 |

The 3-generation count is **locked to `k_5 = 5`**. A different bulk
dimension would predict a different number of charged-lepton
generations.

## Independent cavity-geometry confirmation (#68)

The radial overtone ladder at `l=1` on the cavity `[R_MID, R_OUTER]`:
the throat-localization metrics (⟨r⟩−R_MID, throat fraction,
participation ratio) **saturate** as `n` grows, with the lepton modes
(n=0,1,2 = e,μ,τ) showing the monotonic rise and n≳3 reaching the
shell-asymptote 2/3. The cavity geometry supports ~3 throat-localized
modes before shell saturation (PR #68's framing: a saturating crossover,
right after the third generation). This is independent confirmation that
the algebraic `(k_5+1)/2 = 3` matches the cavity-geometry mode count.

## Unification with PR #71

Both BAM lepton-sector structural results come from the same primitive:

```
β_lepton           =  k_5² · (2π)           =  50π    (PR #71)
#generations       =  (k_5 + 1) / 2          =  3      (this PR)
```

`k_5 = 5` (the Tangherlini bulk dimension / γ-lock cap) sets both the
within-lepton mass-scaling coupling and the generation count. The two
quadratic/linear faces of the same topological charge.

## Honest scope

  - **Is:** a clean structural derivation of `#generations = 3` from
    `k_5 = 5` via the odd-k selection (#67) and the γ-lock `l`-cap
    (`hbar_origin_status`). Unified with #71's `β_lepton = k_5²·(2π)`
    through the same `k_5` primitive.
  - **Is not:** a first-principles derivation of `k_5 = 5` itself.
    `k_5` is the Tangherlini bulk dimension, tied to the γ-lock
    cross-species fixed point at `R_OUTER ≈ 1.262` (established but not
    first-principles derived in `hbar_origin_status`).
  - **Honest note on the cavity geometry:** the throat-localized count
    is a **saturating crossover** (PR #68), not a razor-sharp cutoff
    — 3 modes under the saturation criterion, 4 under a tighter
    threshold. The *algebraic* `(k_5+1)/2 = 3` is exact; the cavity
    geometry is independently consistent.

## B4 accounting

`k_5` is a dimensionless integer (topological charge). The count
`(k_5+1)/2` is structural/topological — independent of the dimensionful
anchor.

## Tests

  T1. **Odd-k fermionic selection (recap #67).**
  T2. **Topological-charge cap `k ≤ k_5`.** lepton sector is the odd
      subset of `[0, k_5]`.
  T3. **Count `(k_5+1)/2 = 3`.** odd integers in `[0, 5]` = `{1, 3, 5}` =
      `(e, μ, τ)`.
  T4. **`k_5`-dependence.** `(k_5+1)/2`: 2/3/4/5 for `k_5 = 3/5/7/9`.
  T5. **Cavity-geometry consistency (#68).** independent throat-mode
      count from the saturating crossover is ~3 (consistent, soft).
  T6. **Unification with #71.** same `k_5` primitive:
      `β_lepton = k_5²·(2π)` and `#gen = (k_5+1)/2`.
  T7. **Falsification / B4.** a 4th charged lepton would falsify (would
      require `k_5 ≥ 7`); observation = 3. `k_5` dimensionless integer.
  T8. **Assessment.**

## Verdict structure

  - **THREE_GENERATIONS_FROM_K5** (expected): the 3 charged-lepton
    generations follow from `(k_5+1)/2 = 3` for `k_5 = 5` — the same
    Tangherlini topological-charge primitive that gives `β_lepton =
    k_5²·(2π) = 50π` (PR #71). The lepton depths `(1, 3, 5)` are the odd
    integers in `[0, k_5]`, combining the odd-k fermionic selection
    (#67) and the γ-lock `l`-cap (`hbar_origin_status`). The cavity
    geometry independently shows the throat-localized ladder saturates
    at ~3 modes (#68's crossover) — consistent. Closes the PR #70
    follow-on, leaving `k_5 = 5` itself (Tangherlini bulk dimension /
    γ-lock fixed point) as the established input.

  - **NO_K5_DERIVATION**: the structural form does not match.

## What this leaves open

  - **First-principles `k_5 = 5`.** Tied to the γ-lock cross-species
    fixed point at `R_OUTER ≈ 1.262`, `Σ V_max[0..5] = 22.5` (per
    `hbar_origin_status`); a deeper geometric derivation of why the
    `l`-range caps at 5 specifically is future work.

## Cross-references

  - `docs/hbar_origin_status.md` — the γ-lock cross-species fixed point
    (R_OUTER, `Σ V_max[0..5] = 22.5`).
  - `docs/even_k_absence_research_plan.md` — odd-k fermionic selection
    (#67).
  - `docs/throat_to_shell_transition_research_plan.md` — the cavity
    saturating crossover (#68).
  - `docs/beta_lepton_derivation_research_plan.md` — `β_lepton =
    k_5²·(2π)` (#71); same `k_5` primitive.
  - `experiments/closure_ledger/three_throat_modes_probe.py` — this
    probe.
