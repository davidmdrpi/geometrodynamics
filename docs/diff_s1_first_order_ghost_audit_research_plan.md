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

## CKV zero-mode division — no double-counting (proof)

The ghost field space splits **orthogonally**:

```
c-space = ker(P) ⊕ (row space),     b-space = ker(P†) ⊕ (column space),
ker(P)  = CKVs (rigid rotation),     ker(P†) = moduli (Teichmüller).
```

The Faddeev–Popov determinant is the **primed** determinant `det'(P) =
det'(P†P)^{1/2}` over the **nonzero modes only** — it excludes both kernels.
Numerically (SVD of `∂_τ` on `S¹`, odd `N`): exactly **one zero singular
value** (its right-null vector = the CKV, its left-null = the modulus) and
**`N−1` nonzero** singular values in `det'(P)`. Therefore:

  - the CKV norm `‖φ‖` enters **only** the gauge-orbit volume `Vol(CKG)`; it
    is **not** in `det'(P)`;
  - the modulus norm `‖ψ‖` enters **only** the modulus measure `dL`;
  - each zero mode is divided **exactly once**.

The first draft of this PR divided additionally by the zero-mode norms
(`√L·√L`) **on top of** the CKV automorphism `1/Vol(CKG)` — since the CKV
norm is already inside `Vol(CKG)`, that double-counted the single CKV.
Corrected here: the CKV appears once, in `1/Vol(CKG) = 1/L`.

## The corrected measure table (single-counted)

| factor | zero mode / sector (once) | L-power |
|---|---|---|
| modulus measure `dL` | `ker(P†)` = Teichmüller | `dL` |
| CKV division `1/Vol(CKG) = 1/L` | `ker(P)` = CKV (Vol from `‖φ‖`) | `L^{−1}` |
| `det'(P) = L` (ghost nonzero) | nonzero modes (excl. kernels) | folds into matter `Tr e^{−LH}` |
| matter `det^{−1/2}` | `Tr e^{−LH}` (incl. ghost det) | `L^{−d/2}` |
| **NET MEASURE** | product (each zero mode once) | `dL · L^{−1−d/2}` |

So the BAM loop measure is

```
Z = Σ_sectors ∫ (dL/L) · det^{−1/2}_matter · e^{−S},
```

with the **single** `1/L` the CKV factor (`= 1/Vol(CKG)` = PR #74's `1/(2π)`
at `L = 2π`). The `det'(P) = L` (ghost nonzero modes) folds into the matter
heat kernel `Tr e^{−LH}` per the standard FP procedure — it is **not** an
independent floating factor, and the CKV norm is **not** divided twice.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | audit; distinguish the 4 objects; fix the L-power |
| T2 | four objects | `P=∂_τ` (1st, 1 zero mode = CKV), `P†P=−∂_τ²` (2nd, 1 zero mode) |
| T3 | determinants | `det'(P†P)=L²`; `det'(P)=det'(P†P)^{1/2}=L` (computed) |
| T4 | η-invariant | `η(−i∂_τ)=0` ⟹ `det'(∂_τ)=+L`, no phase; antiperiodic: no CKV |
| T5 | convention | first-order `det'(P)=L` vs second-order `det'(P†P)=L²` (×L) |
| T6 | no double-count | SVD: 1 zero singular value; `det'(P)` excludes the CKV; CKV norm only in `Vol(CKG)` (once) |
| T7 | measure table | single-counted `Z=Σ∫(dL/L)det^{−1/2}_matter e^{−S}`; net `dL·L^{−1−d/2}` |
| T8 | assessment | `BAM_LOOP_MEASURE_GHOST_FIRST_ORDER_DETPRIME_P_EQUALS_L_NET_dL_OVER_L` |

## Established and open

  - **Established (BAM-native):** the Diff(S¹) FP ghost is first-order —
    `det'(P) = det'(∂_τ) = det'(P†P)^{1/2} = L` (`η = 0`, real positive;
    `det'(P†P) = L²` is the second-order square). The CKV is divided **exactly
    once**: `det'(P)` (primed, nonzero modes) excludes the CKV (SVD: 1 zero
    singular value), so the CKV norm enters only `Vol(CKG)`. The measure is
    `Z = Σ ∫ (dL/L) det^{−1/2}_matter e^{−S}` (net `dL·L^{−1−d/2}`), the
    single `1/L` being the CKV factor `= 1/Vol(CKG)`. The `L²` power appears
    only under an explicit second-order ghost convention.

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
