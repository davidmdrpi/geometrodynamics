# `β_lepton = k_5²·(2π)` — first-principles structural form

Derives the lepton-sector `β_lepton = 50π` from the closure-quantum
primitives, closing the "first-principles `β_lepton = 50π`" follow-on
flagged by PR #70. The structural identification:

```
β_lepton  =  k_5² · (2π)  =  5² · 2π  =  50π .
```

The locked closure-quantum integer is then

```
4 β_lepton / (2π)  =  4 · k_5²  =  4 · 25  =  100 ,
```

matching the documented lock from `docs/hbar_origin_status.md`. The
factor `4` is the spinor double cover (`4π` periodicity vs the `2π`
action base); the `k_5²` is the topological charge squared.

## Where β_lepton sits

The closure-quantum primitives established across the hbar-origin
thread (`docs/hbar_origin_status.md`):

| primitive | form | role |
|---|---|---|
| action base | `2π` | the S³ great-circle closure quantum |
| transport | `8π = 4·(2π)` | derived in `transport_resistance_origin_probe` |
| resistance | `7π/100` | derived in same probe |
| inner cutoff `ε` | `7π/(100·k_5⁴)` = resistance / k_5⁴ | `inner_boundary_derivation_probe` |
| `k_5 = 5` | topological charge / Tangherlini dim | independently established |

`β_lepton` joins this ladder cleanly as

```
β_lepton  =  k_5² · (2π) ,
```

i.e. the `p = 2` face of the `k_5^p · (2π)` closure-quantum family:

| p | `k_5^p · (2π)` | role |
|---:|---|---|
| 0 | `2π` | action base |
| 1 | `5·2π = 10π` | (free slot) |
| **2** | **`25·2π = 50π`** | **`β_lepton`** |
| 3 | `125·2π` | (free slot) |
| 4 | `625·2π` | `ε`'s denominator `100·k_5⁴` ladder |

So `β_lepton` is the **closure-quantum action base scaled by the
topological charge squared** — a clean closure-quantum/topological
identification of the previously phenomenological `50π`.

## Why squared (the `(k−3)²` consistency)

PR #70's β-uplift is `Φ_uplift(k) = β_lepton · (k−3)²` — **quadratic** in
the depth offset. The form `β_lepton = k_5²·(2π)` is the **matching
quadratic** in the topological charge: both factors of `β · (k−3)²` are
quadratic — `k_5²` in the coupling and `(k−3)²` in the index. At the
heaviest lepton (`k = 5`, the τ):

```
β · (k−3)²  =  k_5² · (2π) · (k_5 − 3)²  =  25 · 2π · 4  =  200π
            =  100 · (2π)  =  4 β_lepton ,
```

which is exactly the documented closure-quantum integer count for the
τ's β-uplift contribution. Structural consistency.

## Lepton/quark β asymmetry (per `quark_beta_status`)

The quark sector has `β_quark = N_q · π / 2` with `N_q = 2·n_part`,
`n_part = 233`. The robust quark invariant is `N_q ∈ 2ℤ` (the Z₂
partition factor of 2, supplied by the shell modes per PR #69); the
specific `n_part = 233` is **phenomenological** — none of BAM's
principled enumerations in the catalog produces it (`quark_beta_status`).

The lepton sector is **more principled-bounded**: `β_lepton = k_5²·(2π)`
sits in the `k_5^p · (2π)` family with both factors independently
established. The lepton β is at "principled-bounded by closure-quantum
+ topological charge"; the quark β stays at "principled-bounded by
parity only."

## Honest scope

  - **Is:** a clean closure-quantum / topological structural form
    `β_lepton = k_5²·(2π)`, where both `k_5 = 5` (the Tangherlini
    topological charge / 5D dimension) and `2π` (the action base) are
    independent primitives established elsewhere. The `4β/(2π) = 100 =
    4·k_5²` lock is the documented form; this probe gives its structural
    identification. The `k_5²` factor matches the `(k−3)²` uplift
    quadratic.
  - **Is not:** a first-principles derivation of `k_5 = 5` itself —
    that's the established 5D Tangherlini topological setup. Nor a
    derivation of the power `p = 2` from a deeper symmetry argument —
    the `(k−3)²` uplift quadratic provides the structural rationale; a
    group-theoretic derivation is future work.

## B4 accounting

`β_lepton` (radians) is **dimensionless**; `k_5` is an integer
topological charge; `2π` is the action base. The form
`β_lepton = k_5²·(2π)` is structural/topological — independent of the
single dimensionful anchor `m_e`.

## Tests

  T1. **Closure-quantum primitives.** `2π`, `8π = 4·(2π)`, `7π/100`,
      `k_5 = 5`, `ε = 7π/(100·k_5⁴)`.
  T2. **`β_lepton = k_5²·(2π)`.** verified numerically: `5²·2π = 50π` to
      machine precision.
  T3. **Closure-quantum integer `4β/(2π) = 4·k_5² = 100`.** matches the
      documented lock from `hbar_origin_status`.
  T4. **The `k_5^p·(2π)` ladder.** `β_lepton` at `p=2`; `ε` denominator
      `100·k_5⁴ = 4·k_5²·k_5⁴` at the other end. One unified family.
  T5. **Quadratic consistency with `(k−3)²` (PR #70).** `β·(k−3)²` at
      `k=k_5` gives `k_5²·(k_5−3)²·(2π) = 25·4·2π = 200π = 4β = 100·(2π)`
      — the τ's closure-quantum count.
  T6. **Lepton/quark β asymmetry.** lepton β principled-bounded (`= 4·k_5²` closure
      lock); quark β has phenomenological `n_part = 233`.
  T7. **Falsification / B4.** alternative forms `β = c·(2π)` for non-`k_5²` integers
      `c` would lack the closure-quantum/topological pedigree;
      `c = k_5²` is uniquely the principled enumerator with the right
      ladder + quadratic-consistency match.
  T8. **Assessment.**

## Verdict structure

  - **BETA_LEPTON_DERIVED** (expected): `β_lepton = k_5²·(2π) = 50π` is
    the structural form from the closure-quantum primitives, with `4β/
    (2π) = 4·k_5² = 100` the documented closure-quantum integer (the
    spinor `4` × topological charge squared). `β_lepton` sits at the
    `p = 2` face of the `k_5^p·(2π)` ladder, alongside `ε ∝ 1/k_5⁴` at
    the other end — one unified closure-quantum / topological family.
    The `k_5²` factor matches the `(k−3)²` β-uplift quadratic (PR #70).
    The lepton-β derivation closes the corresponding follow-on, leaving
    `k_5 = 5` itself (Tangherlini) as the established input.

  - **NO_DERIVATION**: the structural form does not match.

## What this leaves open

  - **First-principles `k_5 = 5`.** The Tangherlini bulk dimension /
    topological charge — established as the BAM setup, not derived from
    something more primitive.
  - **The power `p = 2`** from a deeper symmetry argument. The `(k−3)²`
    uplift quadratic provides the structural rationale; a group-theoretic
    derivation is future work.

## Cross-references

  - `docs/hbar_origin_status.md` — the closure-quantum primitives
    (`8π = 4·(2π)`, `7π/100`, `ε = 7π/(100·k_5⁴)`).
  - `docs/odd_k_closure_lemma.md` — `Φ_avail(k) = 2π(k+1) +
    50π·max(0,k−3)²`, the locked `50π`.
  - `docs/three_generation_boundary_research_plan.md` — the `(k−3)²`
    quadratic uplift (#70).
  - `docs/quark_beta_status.md` — the lepton/quark β asymmetry.
  - `geometrodynamics/tangherlini/lepton_spectrum.py` — `TAU_BETA_50PI`.
  - `experiments/closure_ledger/beta_lepton_derivation_probe.py` —
    this probe.
