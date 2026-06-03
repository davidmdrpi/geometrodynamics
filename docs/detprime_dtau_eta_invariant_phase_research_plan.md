# The phase / η-invariant framework for det′(∂_τ) (PR #119)

PR #118 used the fact that the η-invariant of `−i∂_τ` on `S¹` vanishes (so
the ghost determinant `det'(∂_τ)` is real, `= +L`) but only **asserted** it
from the `n → −n` symmetry. This PR builds the full **mathematical
framework** for the phase of the first-order determinant `det'(∂_τ)` via the
η-invariant.

## Setup and modulus

`P = ∂_τ` on the circle of circumference `L` is **anti**-self-adjoint
(`P† = −P`); its eigenvalues `2πin/L` (`n ∈ ℤ`) are pure imaginary. The
associated **self-adjoint** operator is `A = −i∂_τ`, eigenvalues `μ_n =
2πn/L`. The **modulus** is unambiguous,

```
|det'(∂_τ)| = det'(P†P)^{1/2} = L      (PR #116/#117/#118).
```

The **phase** is the subtle part this framework fixes.

## The phase formula (Singer / APS)

Choosing a branch for the negative eigenvalues, `ζ_A(s) = ζ_+(s) +
e^{∓iπs} ζ_−(s)` with `ζ_±(s) = Σ_{±μ_n>0} |μ_n|^{−s}`, and `ln det'(A) =
−ζ_A'(0)` gives

```
det'(A) = |det'(A)| · exp[ ± i(π/2)(ζ_{|A|}(0) − η_A(0)) ],
```

so the phase splits cleanly:

```
arg det'(A) = (π/2)·ζ_{|A|}(0)   [local — heat-kernel / scaling]
            − (π/2)·η_A(0)        [topological — spectral asymmetry].
```

`ζ_{|A|}(0) = ζ_+(0) + ζ_−(0)` is the **local** coefficient; `η_A(0) =
ζ_+(0) − ζ_−(0)` is the **η-invariant**, the intrinsic spectral asymmetry.
For periodic `A`: `ζ_{|A|}(0) = 2ζ_R(0) = −1`, `η(0) = 0`, so the phase is
the `−π/2` scaling rotation — removed by the symmetric branch, leaving
`det'(∂_τ)` real.

## The η-invariant with a flux (the Hopf holonomy)

Thread a `U(1)` holonomy `a ∈ [0,1)` through the loop (eigenvalues
`2π(n+a)/L`) — physically the Hopf/Wilson holonomy `∮A = e^{ikχ}`, `a =
kχ/2π`. With the Hurwitz zeta `ζ_H(0,a) = ½ − a`,

```
η_A(0) = ζ_H(0,a) − ζ_H(0,1−a) = 1 − 2a       (0 < a < 1).
```

| flux a | η(0) = 1 − 2a | symmetric? |
|---:|---:|:---:|
| 0 | 1 (→ reduced 0) | ✓ |
| 1/4 | 0.5 | |
| 1/2 | 0 | ✓ |
| 3/4 | −0.5 | |

The reduced η of the **symmetric sectors** vanishes identically: periodic
(zero mode = the CKV, removed) `η' ≡ 0`; antiperiodic `η ≡ 0`. The naive
formula's value `1` at `a = 0` is the **spectral-flow jump** as the zero
mode crosses.

## Concrete determinants (closed forms)

With an IR mass `m`,

```
det(∂_τ + m)_periodic     = 2 sinh(mL/2),
det(∂_τ + m)_antiperiodic = 2 cosh(mL/2).
```

- **Periodic (orientable):** `2 sinh(mL/2) → 0` as `m → 0` (the zero mode);
  the residue gives `det'(∂_τ) = L` (real, positive) — exactly PR #118's
  value, now *derived*.
- **Antiperiodic (Möbius):** `2 cosh(mL/2) → 2`; `det(∂_τ) = 2`
  (L-independent, no zero mode / CKV).

## The BAM connection

The Hopf holonomy `∮A = e^{ikχ}` is the `U(1)` flux `a = kχ/2π` threading
the closure loop, and `η_A(0) = 1 − 2a` is the spectral-asymmetry phase it
induces. The non-orientable (Möbius) `Z₂` half-twist is precisely the
antiperiodic shift `a = 1/2 ⟹ η = 0`; the orientable sector is `a = 0 ⟹ η =
0`. So both physical BAM sectors are `η = 0` (real determinant) — which is
why PR #118's `det'(∂_τ) = +L` carries no anomalous phase. A generic
intermediate holonomy (`0 < a < 1/2`) would give a genuine η-phase
`exp[−i(π/2)(1−2a)]`.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | build the phase/η framework (PR #118 asserted η = 0) |
| T2 | setup + modulus | `∂_τ` anti-s.a., `A = −i∂_τ` s.a.; `|det′| = L` unambiguous |
| T3 | phase formula | `det′(A) = |det′|·exp[±i(π/2)(ζ(0) − η(0))]`; local + topological |
| T4 | η-invariant | `η(0) = 1 − 2a`; reduced `η ≡ 0` for periodic & antiperiodic |
| T5 | concrete dets | `det′(∂_τ)_periodic = L`; `det(∂_τ)_antiperiodic = 2` |
| T6 | BAM connection | Hopf holonomy `a = kχ/2π ⟹ η = 1 − 2a`; `a=0, 1/2 ⟹ η = 0` |
| T7 | scope | rigorously justifies PR #118; intermediate-holonomy η-phase open |
| T8 | assessment | `DETPRIME_DTAU_PHASE_ETA_INVARIANT_FRAMEWORK_BAM_SECTORS_ETA_ZERO` |

## Established and open

  - **Established (BAM-native):** the phase formula `det'(A) =
    |det'|·exp[±i(π/2)(ζ_{|A|}(0) − η_A(0))]`; the local/topological split;
    `η_A(0) = 1 − 2a` tracking the Hopf holonomy; both BAM closure sectors
    (orientable `a=0`, Möbius `a=1/2`) at `η = 0` ⟹ `det'(∂_τ)` real
    (`= +L` periodic, `= 2` antiperiodic) — rigorously justifying PR #118.

  - **Open:** the genuine η-phase `exp[−i(π/2)(1−2a)]` for intermediate Hopf
    holonomy, and its interplay with the matter determinant phase in the
    full measure.

## Cross-references

  - `docs/diff_s1_first_order_ghost_audit_research_plan.md` — PR #118, which
    asserted `η = 0` (here derived) and fixed the ghost determinant.
  - `docs/diff_s1_ghost_determinant_research_plan.md` — PR #117, the ghost
    determinant `det'(P) = L`.
  - `docs/s_bam_loop_measure_research_plan.md` — PR #74, the closure quantum
    / Hopf holonomy.
  - `docs/odd_k_closure_lemma.md` — the non-orientable (antiperiodic, `a =
    1/2`) Möbius sector.

## Run

```
python -m experiments.closure_ledger.detprime_dtau_eta_invariant_phase_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_detprime_dtau_eta_invariant_phase_probe/`.
Expected verdict:
`DETPRIME_DTAU_PHASE_ETA_INVARIANT_FRAMEWORK_BAM_SECTORS_ETA_ZERO`, 8/8 PASS.
