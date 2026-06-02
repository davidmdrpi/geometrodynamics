# First-order Diff(S¹) Faddeev–Popov ghost audit (PR #118)

PR #117 (revised) fixed the ghost determinant to `det'(P) = det'(P†P)^{1/2}
= L`. This follow-up **audits** that result rigorously: it distinguishes the
four objects the reviewer named — `P = ∂_τ`, `P†P = −∂_τ²`, `det'(P)`, and
`det'(P†P)` — handles the **phase / η-invariant** of the first-order
determinant `det'(∂_τ)`, treats the **CKV zero-mode division and the
zero-mode norms** explicitly, and ends with a **measure table** fixing the
net `L`-power. The headline: the FP ghost contributes `L` (first order); it
contributes `L²` **only** if one explicitly adopts a second-order ghost
convention (which over-counts by one power of `L`).

## The four objects

| object | order | spectrum | zero modes |
|---|---|---|---|
| `P = ∂_τ` | first | `2πi n / L`, `n ∈ ℤ` | 1 (`n=0` = CKV) |
| `P†P = −∂_τ²` | second | `(2πn/L)²` | 1 |
| `det'(P)` | — | first-order determinant (n ≠ 0) | — |
| `det'(P†P)` | — | second-order determinant (n ≠ 0) | — |

Zeta-regularized (`ζ_R(0) = −1/2`, `ζ_R'(0) = −½ ln 2π`):

```
det'(P†P) = L²        (ζ(0) = −1, ζ'(0) = −2 ln L),
|det'(P)| = L         (ζ_{|P|}(0) = −1, ζ_{|P|}'(0) = −ln L) = det'(P†P)^{1/2}.
```

verified to machine precision for `L = 2π, 1, 3.32, 5`.

| L | det′(P†P) | det′(P) | det′(P†P)^½ |
|---:|---:|---:|---:|
| 2π | 39.478 | 6.283 | 6.283 |
| 1 | 1.000 | 1.000 | 1.000 |
| 3.324 | 11.050 | 3.324 | 3.324 |
| 5 | 25.000 | 5.000 | 5.000 |

## Phase / η-invariant of `det'(∂_τ)`

The first-order determinant can carry a phase set by the η-invariant. For
`A = −i∂_τ` (self-adjoint, eigenvalues `2πn/L`), the spectrum is symmetric
under `n → −n`, so

```
η_A(s) = Σ_{n≠0} sign(2πn/L) |2πn/L|^{−s} ≡ 0   ⟹   η_A(0) = 0.
```

With `η = 0` there is **no anomalous phase**: `det'(∂_τ) = +L` (real,
positive). In the antiperiodic / Möbius sector the eigenvalues are
`2πi(n+½)/L` — still symmetric (`η = 0`) but with **no zero mode**, hence
**no CKV** in the non-orientable sector (a clean tie to the odd-k structure).

## First- vs second-order ghost convention

The physical FP for one real reparametrization symmetry is the first-order
`bc` system:

```
Δ_FP = ∫ Db Dc e^{−∮ b ∂_τ c} = det'(P) = L      (first order).
```

A second-order ghost convention (action `b(P†P)c`, i.e. a doubled/complex
ghost) would instead give `det'(P†P) = L²`. The two differ by exactly one
power of `L`; the minimal Faddeev–Popov construction is first-order ⟹ `L¹`.

## CKV zero-mode division + zero-mode norms

The constant mode on a circle of circumference `L` has norm `‖1‖ = √L`. The
ghost system has two zero modes — the c-ghost (the CKV, rigid rotation) and
the b-ghost (dual of the one Teichmüller modulus) — norms `√L·√L = L`. The
ghost **NET** factor `det'(P)/norms` is therefore

```
first-order:   L / L  = 1,
second-order:  L² / L = L      (spurious extra power).
```

## The measure table (net L-power)

| factor | source | L-power |
|---|---|---|
| moduli-space measure `dL/L` | modulus `dL` × CKV automorphism `1/L` | `L^{−1} dL` |
| FP ghost det `det'(P)` (first order) | nonzero ghost modes | `L^{+1}` |
| ghost zero-mode norms (÷) | 1 CKV + 1 modulus, `√L` each | `L^{−1}` |
| **ghost NET** | `det'(P)/norms` | `L^{0}` |
| matter `det^{−1/2}` | Tangherlini / heat-kernel Weyl | `L^{−d/2}` |
| **NET MEASURE (first order)** | product | `dL · L^{−1−d/2}` |

So the BAM loop measure is

```
Z = Σ_sectors ∫ (dL/L) · det^{−1/2}_matter · e^{−S},
```

with the ghost `L`-power fixed at `L¹` in `det'(P)` (net `L⁰` after the
zero-mode norms). The second-order convention would multiply by an extra
`L` (giving `dL·L^{−d/2}`), confirming the first-order fixing is the correct
one.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | audit; distinguish the 4 objects; fix the L-power |
| T2 | four objects | `P=∂_τ` (1st, 1 zero mode = CKV), `P†P=−∂_τ²` (2nd, 1 zero mode) |
| T3 | determinants | `det'(P†P)=L²`; `det'(P)=det'(P†P)^{1/2}=L` (computed) |
| T4 | η-invariant | `η(−i∂_τ)=0` ⟹ `det'(∂_τ)=+L`, no phase; antiperiodic: no CKV |
| T5 | convention | first-order `det'(P)=L` vs second-order `det'(P†P)=L²` (×L) |
| T6 | CKV + norms | `√L·√L=L` ⟹ ghost net `L/L=1` (1st), `L²/L=L` (2nd) |
| T7 | measure table | `Z=Σ∫(dL/L)det^{−1/2}_matter e^{−S}`; net `dL·L^{−1−d/2}` |
| T8 | assessment | `BAM_LOOP_MEASURE_GHOST_FIRST_ORDER_DETPRIME_P_EQUALS_L_NET_dL_OVER_L` |

## Established and open

  - **Established (BAM-native):** the Diff(S¹) FP ghost is first-order —
    `det'(P) = det'(∂_τ) = det'(P†P)^{1/2} = L` (`η = 0`, real positive;
    `det'(P†P) = L²` is the second-order square). After CKV + modulus
    zero-mode norms (`√L` each) the ghost NET is `1`, and the measure is
    `Z = Σ ∫ (dL/L) det^{−1/2}_matter e^{−S}` (net `dL·L^{−1−d/2}`). The `L²`
    power appears only under an explicit second-order ghost convention.

  - **Open:** the absolute normalization of `Z` (the `κ₅²/Λ₅` bulk anchor,
    PR #112) and the multi-loop measure.

## Cross-references

  - `docs/diff_s1_ghost_determinant_research_plan.md` — PR #117, the ghost
    determinant this audit verifies in detail.
  - `docs/tangherlini_fluctuation_determinant_research_plan.md` — PR #116,
    the matter determinant `det^{−1/2}_matter` in the measure.
  - `docs/s_bam_loop_measure_research_plan.md` — PR #74, the `1/(2π)` factor
    (the CKV automorphism `1/L` at `L = 2π`).
  - `docs/odd_k_closure_lemma.md` — the non-orientable sector (antiperiodic
    ghosts, no CKV).

## Run

```
python -m experiments.closure_ledger.diff_s1_first_order_ghost_audit_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_diff_s1_first_order_ghost_audit_probe/`.
Expected verdict:
`BAM_LOOP_MEASURE_GHOST_FIRST_ORDER_DETPRIME_P_EQUALS_L_NET_dL_OVER_L`, 8/8 PASS.
