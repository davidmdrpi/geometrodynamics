# Odd-k / generation-sector survival under a deformed bulk geometry (PR #183)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The generation sector is topologically protected

PR #174 derived the odd-k charged-lepton ladder {1, 3, 5} (= 3 generations)
from the **non-orientable** bulk: the throat closure `T = iσ_y` has `T² = −I`
(the Pin⁻ structure on RP²), so `T^k` is off-diagonal for odd k (the
orientation-reversing, fermionic sector) and diagonal for even k (orientable,
bosonic). PR #181/#182 showed the order-field winding survives continuous
evolution and changes only at an amplitude zero. This probe closes that
structure one level up — at the **bulk geometry**: does the odd-k generation
sector survive when the bulk is deformed?

It does, and for the same topological reason. The grading is set by two
**metric-independent** invariants, and a smooth deformation acts on them by
orientation-preserving conjugation — leaving them exactly invariant. The
sector can change only at a genuine **topology change**.

## The two invariants (both metric-independent)

| invariant | value | meaning |
|---|---:|---|
| deck det (brane `S²/antipodal`) | `−1` | RP² **non-orientable** |
| deck det (bulk `S³/antipodal`) | `+1` | RP³ orientable |
| `½ tr T²` | `−1` | Pin⁻ (`T² = −I`) |

The antipodal deck map is `−I` in **any** linear frame, so `det = (−1)^dim` —
purely a function of the dimension parity, metric-independent. The spin
closure `T² = −I` is the Pin⁻ structure forced by `w₁² = w₂` on RP²
(Stiefel–Whitney classes, topological). The grading:

```
tr(T^k) = 2 cos(kπ/2) = 0   for odd k  (off-diagonal, fermion, non-orientable)
                        = ±2  for even k (diagonal, boson, orientable)
```

## Survival under smooth deformation

A smooth metric/frame deformation acts on the holonomy by
orientation-preserving **conjugation** `T → U T U†` and on the deck map by a
**GL⁺ frame change** `A → F A F⁻¹` — neither of which can change a determinant
sign or a trace. Measured across **1000 random** such deformations:

| quantity | stays at | to |
|---|---:|---:|
| `½ tr T²` | `−1` | `10⁻¹⁵` |
| deck det (brane) | `−1` | `10⁻¹⁵` |
| deck det (bulk) | `+1` | `10⁻¹⁵` |

Named deformations — a **Berger squash** of the S³ Hopf fiber, a
**tidal-charge** radial stretch — leave the bulk orientability `+1` exactly.
The odd-k grading rides every smooth bulk deformation untouched.

## The generation count survives

The realized rungs are the odd `k ≤ k₅ = D_bulk = 5`:

```
{1, 3, 5} = (k₅ + 1)/2 = 3 generations     (= LEPTON_BASELINE_DEPTHS)
```

Both inputs are **topological**: `D_bulk = 5` is the bulk dimension (an
integer, unchanged by any metric deformation) and the odd-parity selection is
the orientability grading (preserved by deformation). So the 3-generation
count is not an artifact of the round metric — it is fixed by the bulk
dimension and the orientability class, and survives every smooth deformation.

## Changes only at a topology change

The grading changes **only** at a genuine topology change — never by a smooth
metric deformation. Metric deformations act by conjugation and exactly
preserve `½ tr T²`. The only path that flips the sector is a **non-metric**
deformation of the spin closure itself,

```
T(θ) = exp(iθσ_y) ,   ½ tr T² = cos 2θ ,
```

driving `T²` from `−I` (`θ = π/2`, non-orientable / fermionic) to `+I`
(`θ = 0`, orientable / bosonic). Its orientability invariant crosses **zero**
at `θ = π/4` — a **degenerate spin structure**, the topology-change event.
This is the exact bulk-level analog of the #182 phase slip: the invariant is
locally constant and jumps only where it passes through the singular
configuration (`½ tr T² = 0` here, `|q| = 0` there). Smooth bulk deformations
never move `θ`, so they can never reach the crossing.

## Robustness

Over **500 random** smooth bulk deformations, the odd-k grading is preserved
in every case — `tr(T^k) ≈ 0` for `k ∈ {1, 3, 5}` (the fermionic rungs stay
fermionic) and `½ tr T² = −1` — **500/500**. The generation sector is rigid
against arbitrary smooth deformation of the bulk geometry.

## Unity with #181/#182, and scope

The generation sector is to the **bulk geometry** what the order-field winding
is to the **soliton** (#181/#182): a topological charge robust to smooth
deformation, changing only at a singular / topology-change event
(`½ tr T² = 0` here; `|q| = 0` there). So the #174 round-metric derivation of
the odd-k {1,3,5} ladder is not special — the 3-generation structure is
topologically **protected** against any smooth deformation of the bulk.

**Scope.** The invariance is **exact** (topological: the deck determinant and
the Stiefel–Whitney / spin-closure class are metric-independent). The
deformations are within the orientability/spin-preserving class (smooth
metric/frame changes); a genuine topology change (a different antipodal
quotient, or the degenerate `θ = π/4` spin structure) lies **outside** the
smooth family by construction. The odd-k ladder and the count are read from
the #174 closure algebra — this probe establishes their **robustness**, not a
re-derivation. The result is purely topological; weak-field is not even
invoked.

## Reproduce

```bash
python -m experiments.closure_ledger.odd_k_generation_survival_probe
# Verdict: ODD_K_GENERATION_SECTOR_SURVIVES_BULK_DEFORMATION_TOPOLOGICALLY_PROTECTED_CHANGES_ONLY_AT_A_TOPOLOGY_CHANGE
```
