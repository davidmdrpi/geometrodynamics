# The BAM sector-phase ledger: continuous η-phases vs discrete Z₂ topology (PR #121)

**Run:** 2026-06-03T05:28:49+00:00

Converts the validated `det'(∂_τ)` η-invariant machinery (PRs #117–#120) into a BAM **sector-phase ledger**: the loop-measure phase factorizes as a **continuous η-phase** (from the U(1) holonomy) times a **discrete Z₂ sign** (from the Möbius/odd-k orientation), with **no double-counting**.

- **Continuous**: η-phase e^{i(π/2)(1−2a)} (U(1) holonomy a; θ ∈ (−π/2,+π/2), never −1)
- **Discrete**: Z₂ sign (−1)^k (Möbius/odd-k orientation, w₁)
- **No double-count**: different groups (U(1) vs Z₂); different geometry (connection vs orientation); η-phase never reaches −1
- **Factorization**: det_full = |det P_a| · e^{i(π/2)(1−2a)} · (−1)^k, each factor once
- **Consistency**: at a = 1/2: η-phase = +1 ⟹ antiperiodic det real, Möbius sign purely (−1)^k

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | convert η-machinery → sector-phase ledger | **PASS** |
| T2 | `T2_two_independent_structures` | two structures: U(1) holonomy (connection) vs Z₂ (w₁) | **PASS** |
| T3 | `T3_continuous_eta_phase` | continuous η-phase θ(a)=(π/2)(1−2a) ∈ (−π/2,+π/2), never −1 | **PASS** |
| T4 | `T4_discrete_z2_sign` | discrete Z₂ sign (−1)^k: +1 torus, −1 Möbius | **PASS** |
| T5 | `T5_sector_phase_ledger` | the sector-phase ledger table | **PASS** |
| T6 | `T6_no_double_counting_proof` | no double-count: groups, geometry, no collision on −1 | **PASS** |
| T7 | `T7_factorized_measure_phase` | factorization det_full = |det|·η-phase·Z₂, each once | **PASS** |
| T8 | `T8_assessment` | BAM_SECTOR_PHASE_LEDGER_..._NO_DOUBLE_COUNT | **PASS** |

## The sector-phase ledger

| contribution | source | object | group | value | type |
|---|---|---|---|---|---|
| magnitude | twisted spectrum | `|det P_a| = 2 sin(πa)` | ℝ₊ | `—` | continuous |
| η-phase | U(1) holonomy a (Hopf/Wilson) | `spectral asymmetry η(0) = 1−2a` | U(1) | `e^{i(π/2)(1−2a)}` | CONTINUOUS |
| Z₂ sign | Möbius / odd-k orientation | `w₁ / e^{ikπ}` | Z₂ | `(−1)^k` | DISCRETE |
| (scaling ζ(0) phase) | local heat-kernel | `ζ(0)` | — | `absorbed (ζ(0)=0 for twisted; branch for periodic)` | removed |

## Continuous η-phase (confined to the right half-circle)

| a | θ = (π/2)(1−2a) | e^{iθ} | Re |
|---:|---:|---|---:|
| 0.0 | 1.5708 | +0.0000+1.0000i | 0.0 |
| 0.25 | 0.7854 | +0.7071+0.7071i | 0.7071 |
| 0.3333 | 0.5236 | +0.8660+0.5000i | 0.866 |
| 0.5 | 0.0 | +1.0000+0.0000i | 1.0 |
| 0.6667 | -0.5236 | +0.8660-0.5000i | 0.866 |
| 0.75 | -0.7854 | +0.7071-0.7071i | 0.7071 |
| 1.0 | -1.5708 | +0.0000-1.0000i | 0.0 |

`θ(a) ∈ (−π/2, +π/2)` for `a ∈ (0,1)` ⟹ the η-phase lies in the **open right half-circle** (`Re > 0`), reaching `+1` only at `a = 1/2` and **never `−1`** — so the Möbius `−1` is inaccessible to it.

## Factorized measure phase

```
det_full = |det P_a| · e^{i(π/2)(1−2a)} · (−1)^k
         = [magnitude] · [continuous η-phase] · [discrete Z₂ sign]
```

| a | k | \|det\| | η-phase | Z₂ | det_full |
|---:|---:|---:|---|---:|---|
| 0.25 | 1 | 1.414 | +0.707+0.707i | -1 | -1.000-1.000i |
| 0.25 | 2 | 1.414 | +0.707+0.707i | 1 | +1.000+1.000i |
| 0.333 | 3 | 1.732 | +0.866+0.500i | -1 | -1.500-0.866i |
| 0.5 | 1 | 2.0 | +1.000+0.000i | -1 | -2.000+0.000i |
| 0.5 | 2 | 2.0 | +1.000+0.000i | 1 | +2.000+0.000i |

At `a = 1/2` the η-phase is `+1` (real det), so the antiperiodic determinant's Möbius sign is **entirely** the discrete `(−1)^k` — the cleanest demonstration of no double-counting.

## Verdict

**BAM_SECTOR_PHASE_LEDGER_CONTINUOUS_ETA_TIMES_DISCRETE_Z2_NO_DOUBLE_COUNT.** THE BAM LOOP-MEASURE PHASE FACTORIZES INTO A CONTINUOUS η-PHASE (U(1) HOLONOMY) TIMES A DISCRETE Z₂ SIGN (MÖBIUS / ODD-K ORIENTATION), WITH NO DOUBLE-COUNTING. PRs #117–#120 built and lattice-validated the det'(∂_τ) η-invariant machinery; this ledger separates its two — and only two — phase/sign sources.

THE TWO STRUCTURES. A closure loop carries two independent twists: a U(1) holonomy a ∈ [0,1) (the connection / Hopf-Wilson line ∮A = e^{2πia}, continuous, π₁ of the fibre) and an orientation class (orientable vs Möbius, the winding parity k mod 2 = w₁ ∈ H¹(loop;Z₂), discrete). These are orthogonal topological data — a connection and an orientability.

THE CONTINUOUS η-PHASE. The twisted determinant is det P_a = 2 sin(πa)·e^{iθ(a)} with θ(a) = (π/2)η_A(0) = (π/2)(1−2a). Since the twisted operator has no zero mode, ζ(0) = 0 and the phase is PURELY the η-invariant piece. As a sweeps (0,1), θ(a) sweeps (−π/2,+π/2): the η-phase is confined to the OPEN RIGHT HALF-CIRCLE (Re > 0), reaching +1 only at a = 1/2 and NEVER reaching −1.

THE DISCRETE Z₂ SIGN. The Möbius / non-orientable closure contributes the orientation sign of the odd-k lemma, e^{ikπ} = (−1)^k: +1 for even k (orientable / torus cover), −1 for odd k (non-orientable / Möbius half-twist) — a discrete element of {±1}, the w₁ holonomy, independent of the connection a.

NO DOUBLE-COUNTING. The two are not the same contribution counted twice, for three independent reasons: (a) different groups — U(1) (continuous) vs Z₂ (discrete); (b) different geometry — the connection/holonomy (η, spectral asymmetry) vs the orientation/w₁ (orientability); (c) no collision on the nontrivial element — θ(a) ∈ (−π/2,+π/2) confines the η-phase to the open right half-circle, so it can NEVER equal −1; the Möbius −1 is inaccessible to the continuous phase and is purely the discrete Z₂. In particular, at the antiperiodic point a = 1/2 the η-phase is exactly +1 (det P_{1/2} = 2, real), so all of the Möbius character there is carried by the separate (−1)^k.

THE FACTORIZED PHASE. det_full = |det P_a| · e^{i(π/2)(1−2a)} · (−1)^k = magnitude × continuous η-phase × discrete Z₂ sign, each factor counted exactly once. The continuous and discrete parts multiply; they do not overlap.
