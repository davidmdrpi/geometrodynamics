# The global regular 5D embedding of the BAM throat (PR #168)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The closable step, closed

PR #167 left a precise gap. The throat's 4D exotic stress is of the
tidal-charge / bulk-Weyl form, and the no-on-brane-exotic reading meets its
**necessary** conditions (Ricci-scalar zero; a negative tidal charge that
excludes a real on-brane Maxwell source). But the **sufficient** step — an
explicit **global regular** 5D bulk whose projected Weyl tensor sources the
brane stress — was pending. This probe builds it, as a global, regular,
*exact* embedding (not a Campbell–Magaard local-existence series).

## The construction

The BAM brane is the **equatorial (χ = π/2) slice** of the 5D
Schwarzschild–Tangherlini bulk

```
ds²₅ = −F(ρ) dt² + dρ²/F(ρ) + ρ²[dχ² + sin²χ (dθ² + sin²θ dφ²)],
F(ρ) = 1 − μ/ρ²,   μ = r_s².
```

The equator χ = π/2 is a fixed-point set of the reflection χ → π−χ, hence
**totally geodesic** (`K_μν = 0`) — a tension-free, matter-free Z₂ brane.
Its induced 4D metric is exactly the BAM throat `f(r) = 1 − (r_s/r)²`
(with `r = ρ`). The extra dimension is the **compact** polar angle χ. The
construction works *only* for the pure-tidal `M = 0` form — a Schwarzschild
`1/r` term has no counterpart in the 5D Tangherlini `1/ρ²` potential — so
the gate has teeth: it selects exactly BAM's metric and ties the bulk mass
to the throat scale (`μ = r_s²`).

## The three printed checks

| check | result |
|---|---|
| 1. induced metric = `f`, and `K_μν = 0` | exact; `K_μν = 0` (totally geodesic) |
| 2. projected bulk Weyl `E_μν = −G⁴_μν` | max `\|E + G⁴\| ≈ 1e-8` (the tidal fluid) |
| 3. 5D field equations (bulk vacuum) | Ricci-flat, max `\|R⁵_MN\| ≈ 3e-7` |

Check 2 is the **bulk-Weyl mechanism, made explicit**: the brane's
effective exotic stress (`ρ_eff < 0`) *is* the projected Weyl tensor of the
ordinary 5D vacuum, verified by an actual solution rather than asserted.

## The regularity gate

The 5D Kretschmann — a coordinate **invariant** — is `K₅ = 72 μ²/ρ⁸`. The
closed form is validated against the numerical curvature to `~1e-6` at
`ρ ≥ 1.5 r_s` (finite-differencing is unreliable nearer the horizon, where
`g_ρρ = 1/F → ∞` is a **coordinate**, not curvature, breakdown). On the
exterior domain `ρ ∈ [r_s, ∞)` the invariant is **finite everywhere**,
rising to its maximum `72/r_s⁴` at the throat `ρ = r_s` — a regular value.
The only curvature singularity, `ρ = 0`, lies **behind** the regular 5D
Killing horizon `ρ = r_s`, off the brane's exterior domain; the extra
dimension `χ ∈ [0, π]` is compact and regular at its poles. The embedding
is **global and regular — the gate passes**.

## What this closes

PR #167's bulk-Weyl reading is no longer a "consistent-with": it is
**realised** by an explicit global regular 5D vacuum.

- No exotic brane matter; the brane is matter-free and carries no
  fundamental gauge field (the two conditions #167 flagged are both met by
  an actual bulk).
- The `f = 0` throat is identified as the **regular 5D Killing horizon**
  `ρ = r_s` — an improvement on #167's "horizon/null" caveat: it is regular
  (`K₅ = 72/r_s⁴` finite), not singular, with the true singularity `ρ = 0`
  safely behind it.

## Honest residue

- The throat still sits at a horizon — now known **regular**.
- The brane is the tension-free totally-geodesic equatorial slice (a clean
  special case; `μ = r_s²` fixes the bulk mass).
- It is the **exterior** embedding `ρ ≥ r_s`.

## Reproduce

```bash
python -m experiments.closure_ledger.global_regular_5d_embedding_probe
# Verdict: GLOBAL_REGULAR_5D_EMBEDDING_EXISTS_BAM_THROAT_IS_EQUATORIAL_TANGHERLINI_SLICE
```
