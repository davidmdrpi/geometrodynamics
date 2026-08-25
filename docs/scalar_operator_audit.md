# Correcting the radial scalar operator, and pricing the correction

**Module:** `geometrodynamics/tangherlini/operator_audit.py`
**Operators:** `geometrodynamics/tangherlini/radial.py`
**Probe:** `python -m experiments.closure_ledger.scalar_operator_audit_probe` (8/8)
**Tests:** `tests/test_scalar_operator_audit.py` (32)
**Renderer:** none — this round is an operator migration and an audit, not a new picture

---

PR #270 found, while doing something else, that `tangherlini.radial.V_tangherlini`
was not the master potential of a minimally coupled massless scalar. It reported
the discrepancy and changed nothing, because a silent one-line replacement would
have turned fifty green regression tests into an unlabelled mixture of genuine
regressions and expected scientific changes.

This round makes the correction, and prices it.

## 0 · the correction

For `ds² = −A dt² + A⁻¹dr² + r²dΩ_n²` a minimally coupled massless scalar
`Φ = e^{−iωt}Y_ℓ(Ω)R(r)` obeys

    (1/r^n) ∂_r(r^n A R') + (ω²/A − ℓ(ℓ+n−1)/r²) R = 0 ,

and the unique first-derivative-free Schrödinger form, with `dr* = dr/A` and
`ψ = r^{n/2}R`, carries

    V_scalar = A[ ℓ(ℓ+n−1)/r² + n(n−2)A/(4r²) + n A'/(2r) ] .

At `n = 3` that is `A[(ℓ(ℓ+2) + 3/4)/r² + (9/4)r_h²/r⁴]`, while the repository
carried `A[ℓ(ℓ+2)/r² + 3r_h²/r⁴]`, so

    V_scalar − V_legacy = 3A²/(4r²) .

**Three independent confirmations, none of them a citation:**

| check | result |
|--|--|
| the gap equals `3A²/(4r²)` in closed form | `2.4e-15` |
| the gap carries no `ℓ` | `2.4e-15` |
| the flat limit reproduces Bessel | `2.2e-16` |
| bitwise agreement with `dynamics.master_potential` | exact |

The flat limit is what settles *which* operator is the scalar one, with no
appeal to authority: at `r_h → 0` the regular solution is
`φ = J_{ℓ+1}(ωr)/r`, so `ψ = r^{1/2}J_{ℓ+1}(ωr)`, and Bessel's equation gives
`V → ((ℓ+1)² − ¼)/r² = (ℓ(ℓ+2) + 3/4)/r²` — including the `3/4`. The legacy
operator **fails** that same test, which is what makes it a discriminator rather
than a formality. The general-`n` form was verified symbolically at `n = 2 … 6`.

**This is a bug, not a convention.** The old generic name implied the canonical
scalar operator and the implementation was short of it by an `ℓ`-independent
term.

### what the module now exposes

| name | role |
|--|--|
| `V_tangherlini_legacy` | **frozen.** The pre-#271 operator, kept only to reproduce archived runs, and documented as *not* the scalar master potential |
| `V_scalar_tangherlini` | the corrected operator, general `n` |
| `V_tangherlini` | the canonical name, now delegating to the corrected operator |

## 1 · the eigenvalues barely move, and move less as `ℓ` rises

| `ℓ` | legacy `ω_{ℓ,0}` | corrected | shift | min overlap |
|--|--|--|--|--|
| 0 | `1.00065891` | `1.00198000` | `+0.1320%` | `0.999998` |
| **1** | **`1.05472694`** | **`1.05582653`** | **`+0.1043%`** | `0.999998` |
| 2 | `1.13156946` | `1.13239953` | `+0.0734%` | `0.999998` |
| 3 | `1.21908274` | `1.21966785` | `+0.0480%` | `0.999999` |
| 4 | `1.30869618` | `1.30909388` | `+0.0304%` | `0.999999` |
| 5 | `1.39597349` | `1.39624205` | `+0.0192%` | `0.999999` |

The monotone fall with `ℓ` is not a coincidence. **An eigenvalue averages the
potential against a bound state**, so an `ℓ`-independent shift matters least
where the centrifugal term already dominates. This is the reassuring half of the
audit: the spectrum is stable at the `10⁻³` level and the eigenfunctions are
essentially unmoved.

## 2 · the barrier sums are not protected, and that is where the meaning moves

A barrier height **reads** the potential directly instead of averaging it, so
`V_max` and its sums shift by an order of magnitude more than `ω` does:

| channels | legacy sum | residual vs `22.5` | corrected sum | residual | `R_OUTER` legacy | corrected |
|--|--|--|--|--|--|--|
| `ℓ = 1..5` | `22.00824` | `−2.19%` | `22.33119` | **`−0.75%`** | `1.28737` | `1.26788` |
| `ℓ = 0..5` | `22.45268` | `−0.21%` | `22.83642` | **`+1.50%`** | `1.26227` | `1.24614` |

> **The two `γ` statements in the tree move in opposite directions.**

- The **canonical** claim — README's *"Pinhole γ ≈ Σ V_max[1..5] … −2.2 % off the
  locked γ = 22.5"* — **improves threefold**, to `−0.75%`, with nothing tuned.
- The **`ℓ = 0..5`** claim — that adding the `ℓ = 0` 5D channel closes the
  pinhole gap — **breaks**: it goes from a `−0.21%` undershoot to a `+1.50%`
  overshoot, and the sum closest to `22.5` **swaps** from `ℓ = 0..5` to
  `ℓ = 1..5`.

That claim is **withdrawn, not replaced.** The corrected `ℓ = 1..5` sum does land
nearer `22.5` than the old one did, but that is an observation, not a derivation
of why `22.5`.

Where it is load-bearing: `docs/lepton_axioms.md`
(*"`hard_pinhole_gamma ≈ 22.5` matches `Σ_{l=0..5} V_max(l)` ≈ 22.453"*) and
`docs/hbar_origin_status.md` (*"γ_lepton lock `Σ V_max[0..5] = 22.5` → R_OUTER ≈
1.2623, ω = 1.054"*). Both now carry a marker pointing here.

## 3 · what survives exactly

`ΔV` carries no `ℓ`. So the cross-`ℓ` perturbation operator is unchanged:

    max | (V_{ℓ₂} − V_{ℓ₁})_corrected − (V_{ℓ₂} − V_{ℓ₁})_legacy |  =  3.6e-15

over every pair to `ℓ = 5` — **algebraically exact**. Its *matrix elements* are
not, because they are taken between eigenfunctions that drift.

> **Structure invariant, numbers shifted.** That distinction is the whole
> partition of the audit.

Ratios behave similarly: `α_q(ℓ,0)` is a ratio of throat derivatives normalised
to `ℓ = 1`, so the common part of the shift divides out and only the differential
survives. Measured rather than assumed — a ratio is not automatically safe.

Closed-orbit WKB actions sit in between: they integrate `√(ω² − V)` against the
potential directly, without the bound state's averaging, so they move more than
the eigenvalues and less than the barrier sums. Each action is evaluated at *its
own* operator's ground frequency, so the comparison is between two
self-consistent orbits rather than one orbit against the other's potential.

## 4 · the one narrow downstream re-derivation

Rather than recompute the historical tree — which would only propagate an
unresolved branch choice — exactly three geometries were passed through the
**locked** lepton Hamiltonian with **nothing retuned**:

- **A** — keep `R_OUTER = 1.26` fixed and let `γ` float
- **B** — enforce `Σ_{1..5} V_max = 22.5`
- **C** — enforce `Σ_{0..5} V_max = 22.5`

| case | `R_OUTER` | `γ` | `m_μ` err | `m_τ` err |
|--|--|--|--|--|
| baseline, locked | — | `22.50000` | `−0.04%` | `+0.12%` |
| legacy `R=1.26`, `γ[0..5]` | `1.26000` | `22.45268` | `+3.78%` | `+4.04%` |
| legacy `R=1.26`, `γ[1..5]` | `1.26000` | `22.00824` | `+63.78%` | `+65.58%` |
| **A** corrected `R=1.26`, `γ[1..5]` | `1.26000` | `22.33119` | **`+15.16%`** | `+15.71%` |
| **A** corrected `R=1.26`, `γ[0..5]` | `1.26000` | `22.83642` | **`−20.52%`** | `−20.89%` |
| **B** corrected root `[1..5]=22.5` | `1.26788` | `22.50000` | **`−0.04%`** | `+0.12%` |
| **C** corrected root `[0..5]=22.5` | `1.24614` | `22.50000` | **`−0.04%`** | `+0.12%` |

The baseline row reproduces the locked spectrum and the legacy `R=1.26` row
reproduces the documented *"recovering the muon mass within 3.8%"*, so the
harness is the real chain and not a re-implementation of it.

**B and C are bit-identical.** `compute_knotted_lepton_spectrum` discards
`r_outer` outright (`del l, n_points, rs, r_outer`); the locked block consumes
the geometry **only** through the scalar `hard_pinhole_gamma`. So once `γ` is
enforced, the channel-set choice leaves **no trace in any observable**, and this
comparison cannot decide it.

> **`γ = 22.5` is the selector; `R_OUTER` is downstream of it.**

Which inverts the anticipated outcome. It is not that `γ` was diagnostic and
`R_OUTER` the real selector — it is *fixing* `R_OUTER` and letting `γ` float
that breaks the ladder, by `+15%` or `−21%` depending on channel set.

**And the correction weakens the "geometry supplies `γ`" story**, even while
improving the `1..5` residual in isolation. Under the legacy operator the
canonical `R = 1.26` geometry produced `γ[0..5] = 22.453`, `0.21%` from the lock,
and the masses landed within `3.8%` — that near-coincidence *was* the claim.
Under the corrected operator no channel set at `R = 1.26` lands near `22.5`, and
both damage the ladder at the `15–21%` level.

The reason is sensitivity: `d ln m_μ / d ln γ = −16.6`, so a **sub-percent
geometric residual is not a small residual** in this chain. A `0.5%` error in `γ`
is an `8%` error in the muon.

**What this does not settle.** The channel-set question is not decidable by the
lepton observables — it has to be settled by what `γ` *means* geometrically,
because the masses only ever saw the scalar.

## 5 · the ledger

The question is **not** whether the old tests stay green. It is which published
claims are algebraically untouched, which keep their meaning with different
digits, and which no longer say what they said.

| claim | verdict |
|--|--|
| cross-`ℓ` perturbation operator `V_{ℓ+2} − V_ℓ` | **exactly invariant** |
| Hopf, Pin⁻, odd-`k` ladder, antipodal parity | **exactly invariant** (no dependence; not re-run) |
| `α_q(ℓ,0)` throat-flux ratios | numerically shifted |
| `ω_{ℓn}` and eigenfunctions | numerically shifted |
| cross-`ℓ` transport matrix elements | numerically shifted |
| closed-orbit / WKB radial actions | numerically shifted |
| the `1.054` factor, `ω(1,0)` at the γ-locked geometry | numerically shifted |
| pinhole `γ = Σ V_max[1..5]` vs the locked `22.5` | numerically shifted (**improves**) |
| `R_OUTER` as the `γ = 22.5` fixed point | numerically shifted |
| *"adding the `ℓ = 0` 5D channel closes the γ discrepancy"* | **interpretation changed** |
| any generation or mass result whose chain runs through `γ`, `R_OUTER`, or a radial eigenvalue | **interpretation changed** |

The `1.054` entry deserves its own line: `ω(1,0)` becomes `1.0558`, so the quoted
`1.054` becomes `1.056`. That exceeds the `0.04%` Compton-bridge tolerance the
README states, so it needs **re-quoting, not re-deriving** — the mechanism is
untouched, the number is not.

### what is deliberately not re-run

The Hopf fibration, the Pin⁻ structure on the exchange `ℝP²`, the odd-`k` winding
ladder and antipodal parity have **no dependence** on the radial scalar operator.
Proximity is not dependence, and re-running them would manufacture the appearance
of an audit without adding one.

## 6 · what the test suite did *not* protect

Flipping the operator broke exactly **two** tests out of `1582`, and both were
PR #270's own bookkeeping — the ones asserting that #270 had changed nothing.
Every other test passed unchanged.

That is an audit finding in its own right, and not a comforting one. The γ sums,
the `R_OUTER` fixed point and the `1.054` factor are **not regression-locked
anywhere in the suite**; they live in prose, in probe narratives and in docs. A
silent replacement would have sailed through CI. The correction was caught by a
derivation, not by a test, and nothing in the tree would have caught it
otherwise.

## Scope

- **The operator correction is complete and exact.** The audit's impact numbers
  are measured on the canonical construction at `r_h = 1`, `R_OUTER = 1.26`.
- **Downstream mass and generation chains are classified, not re-derived.** Every
  chain running through `γ`, `R_OUTER` or a radial eigenvalue must be re-derived
  from the corrected operator before its number is quoted again. Doing that
  silently here would repeat the original error in the opposite direction.
- **The γ narrative is reopened, not resolved.** No replacement claim is offered.
