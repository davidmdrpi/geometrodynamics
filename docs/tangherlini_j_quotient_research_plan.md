# Tangherlini J-quotient consistency and brane non-orientability (PR #169)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The topological root of the non-orientable throat

PR #167 found the throat gluing is **non-orientable**; PR #168 embedded the
throat as the equatorial slice of the 5D Tangherlini bulk. This probe says
**why** — a clean dimension-parity statement about the antipodal (J)
quotient:

```
bulk angular sphere   S³ / antipodal = RP³   ORIENTABLE      (det +1)
brane angular slice    S² / antipodal = RP²   NON-ORIENTABLE  (det −1)
```

The antipodal involution `J: x ↦ −x` on `Sⁿ` has orientation determinant
`(−1)^{n+1}`: orientation-**preserving** for **odd** `n` (orientable
quotient), orientation-**reversing** for **even** `n` (non-orientable
quotient). The bulk's angular sphere is `S³` (odd ⟹ `RP³` orientable); the
throat mouth is the brane's angular `S²` (even ⟹ `RP²` non-orientable).

## The consistency

The **same** free isometric involution `J = −I` acts oppositely on the two
spheres because they sit one dimension apart, on opposite sides of the
parity:

| `n` | orientation det `(−1)^{n+1}` | `RPⁿ` orientable? |
|---|---:|---|
| 1 | +1 | yes |
| 2 | **−1** (brane mouth) | **no** |
| 3 | **+1** (bulk) | **yes** |
| 4 | −1 | no |
| 5 | +1 | yes |

- `J` is a **free** involution (`−x = x ⟹ x = 0 ∉ Sⁿ`) → the quotient is a
  smooth manifold (no orbifold points).
- `J` is an **isometry** of the round sphere (`JᵀJ = I`) → the 5D
  Tangherlini metric (round angular part) **descends** to the quotient.

So the J-quotient is consistently orientable on the bulk and
non-orientable on the brane mouth — the non-orientable throat is **forced**
by the single-dimension drop from bulk to mouth, not an extra assumption.

## The Tangherlini realization (ties to #167 / #168)

In the #168 coordinates the `S³` antipodal map is
`(χ, θ, φ) ↦ (π−χ, π−θ, φ+π)` (= `x ↦ −x`). It **fixes** the equatorial
`χ = π/2` brane and **restricts** there to the `S²` antipodal map
`(θ, φ) ↦ (π−θ, φ+π)`. The round-angular Tangherlini metric is `J`-invariant
and descends: the bulk becomes the **orientable `RP³`-angular** Tangherlini,
the throat mouth becomes the **non-orientable `RP²`**. The #167
non-orientable throat is exactly this `RP²` cross-cap inside the orientable
`RP³` bulk of #168.

## Consistency with spin / Pin (a remark, not a new derivation)

`RP³ ≅ SO(3)` is orientable, parallelizable, and admits a spin structure —
a consistent arena for bulk fields. `RP²` is non-orientable and admits
**no** spin structure, only a **Pin** structure: the half-twist carrying
the spin-½ / fermionic character. This is the same orientability grading
already in BAM as the C-swap (`C = iσ_y`, `T² = −1`; PR #63) and the
even-`k` absence (`k mod 2` = orientability / spin-statistics; PR #67). The
throat's `RP²` non-orientability is the topological seat of its Pin/spin-½
structure.

## Reproduce

```bash
python -m experiments.closure_ledger.tangherlini_j_quotient_probe
# Verdict: J_QUOTIENT_CONSISTENT_BULK_RP3_ORIENTABLE_BRANE_RP2_NON_ORIENTABLE
```
