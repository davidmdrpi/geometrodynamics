# Spin 0 against spin 2 on the same S²

> **Framing.** A spin-2 **field on a fixed** `S²` — the tensor analogue of the
> scalar wave, *not* a solution of the 4-D linearised Einstein equations.
> General relativity in 2+1 dimensions has no propagating tensor modes at all.

## The correction

`warped_sphere.py` displays a **scalar** `u(t, θ, φ)` extrinsically, as a radial
height perturbation. A metric perturbation is not that kind of object: `h_ab` is
symmetric and trace-free — spin 2 — and it does not push a surface outward at
all. It *shears* it. Everything below is the list of consequences, each
measured rather than drawn.

| | scalar `u` | tensor `h_ab` |
|---|---|---|
| displayed as | radial height | tidal ellipses |
| local effect | area changes — it breathes | area preserved to `O(h²)` |
| at a pole | free to peak there | **must vanish** |
| at the focus | peaks **on** the antipode | a **ring** around it |
| multipoles | all `ℓ ≥ 0` | **`ℓ ≥ 2` only** |
| smallest source | a point | a **ring** |

## The field and its equation

An axisymmetric spin-`s` field on the unit sphere obeys

```
∂²_t h = ∂²_d h + cot d ∂_d h − (s²/sin²d) h
```

whose eigenvalue on `ₛY_ℓ0` is `−ℓ(ℓ+1)`. So the tensor shares the scalar's
dispersion `ω_ℓ = √(ℓ(ℓ+1))` and therefore its `t = π` refocus; what the
centrifugal term does instead is force `h → 0` at both poles.

Writing `h = sin²d · q` removes that term exactly,

```
∂²_t q = (1/sin⁵d) ∂_d( sin⁵d ∂_d q ) − 6 q
```

and this is integrated in conservative form on a cell-centred grid, so the poles
carry no flux and never appear in a denominator.

## The solver, against exact modes

| `q` | `ℓ` | `ω` | shape error | invariant drift |
|---|---:|---:|---:|---:|
| `1` | 2 | 2.4495 | `3.5e-08` | `3.8e-11` |
| `cos d` | 3 | 3.4641 | `5.6e-05` | `7.0e-13` |
| `7cos²d − 1` | 4 | 4.4721 | `1.3e-04` | `1.2e-13` |

after 10 periods.

## A spin-2 field cannot sit on a pole

`h = sin²d · q` vanishes at both poles for **every** `q`. Two consequences, one
at each end of the run:

| quantity | value |
|---|---:|
| tensor peak at distance | 2.9439 |
| — a ring of radius | 0.1977 |
| amplitude *on* the antipode | `2.2e-06` |
| **scalar** peak at distance | 3.1416 |

At the focus the tensor is a **ring** around the antipode, exactly where the
scalar piles up. And at the source end the same fact says the smallest source a
spin-2 field admits is a ring: **there is no point source of tidal shear.**

## Pure shear, not breathing

| quantity | value |
|---|---:|
| trace | `0.0e+00` |
| area ratio | 1.0000000117 |
| first-order area change | `1.17e-08` |
| `ε²h²/2` prediction | `1.17e-08` |

The ellipse has semi-axes `1 ± εh`, so the first-order area change vanishes
identically and what remains is exactly the second-order term. A radial height
perturbation changes area at first order — that is the difference you can see.

## Spin weight 2

| rotation of the frame | strain |
|---|---:|
| 180° | identical, `1.1e-15` |
| 90° | inverted, 2.00× amplitude |

The strain is `h₊cos2β + h_ˣsin2β`; the factor of two in the angle *is* the
spin weight, visible without any decomposition.

## The caustic: a quarter turn, not a flip

The obvious guess — that passing the antipodal focus swaps the stretch and
compression axes — is **not** what happens. Correlating the outbound waveform
on a fixed ring against all four candidates:

| the outbound front is the inbound one… | correlation |
|---|---:|
| unchanged | +0.3542 |
| **inverted** (a polarisation flip) | −0.3542 |
| phase +90° | −0.8165 |
| **phase −90°** | **+0.8165** |

The outbound front is the **Hilbert transform** of the inbound one: a 90° phase
shift, which turns a peak into a peak-trough pair rather than inverting it.
This is the Gouy shift, Maslov index 1 — and it belongs to the *wave equation*,
not to the spin: the scalar has it too.

## The round trip, where the axes do swap

| | correlation |
|---|---:|
| `h(2π)` with `−h(0)` | **+0.9974** |
| `h(2π)` with `h(0)` | −0.9974 |
| `h(4π)` with `h(0)` | +0.9877 |
| inversion residual | 0.0064 |

Two focal passages — the antipode at `t = π`, home at `t = 2π` — compose to a
half turn, and the field returns as minus itself. **That** is where the
polarisation axes swap. It is not exact because `ω_ℓ = √(ℓ(ℓ+1))` only
approaches `ℓ + ½`; the residual is that dispersion.

## Honest scope

* **Not linearised GR on a spacetime.** 2+1 gravity has no propagating tensor
  modes; this is the spin-2 analogue of the scalar wave on a fixed `S²`.
* **What is faithful**: two polarisations, spin weight 2, `ℓ ≥ 2` only,
  area-preserving shear, and the behaviour at a focus.
* **The Gouy shift is not a spin-2 effect.** It is shared with the scalar, and
  the probe says so rather than claiming it as a difference.
* **One polarisation at a time.** An axisymmetric source drives even parity
  (`h₊`) or odd (`h_ˣ`), never both; the module provides each.
* **The ellipse sizes are a display gain.** Their *shapes* are the solved
  field; the gain is fixed from the run's own peak strain.

## Reproduce

```bash
python -m experiments.closure_ledger.spin2_tidal_probe
# Verdict: A_TENSOR_WAVE_CANNOT_SIT_ON_ITS_OWN_FOCUS  (7/7)

python scripts/geometrodynamics_v43_tidal_sphere.py              # animate
python scripts/geometrodynamics_v43_tidal_sphere.py --still sheet.png
```

```python
from geometrodynamics.viz.spin2_tidal import TidalField
f = TidalField()
f.advance_to(3.14159)
f.peak_ring()          # {'distance': ..., 'amplitude': ..., 'axis': 'radial'}
f.ellipse(2.9, 0.0)    # the curve a ring of test particles becomes
```

## Where this goes next

1. **A ring source of tidal shear.** The spin-2 field already insists on a
   ring; `radial_caustic.py` says a ring is what folds. Those are the same
   statement approached from two sides.
2. **Both polarisations at once.** A non-axisymmetric source drives `h₊` and
   `h_ˣ` together, and the ellipse pattern rotates rather than pulsing.
3. **The throat.** What a spin-2 field does at a catenoidal neck — where the
   frame is not degenerate but the curvature is — is the natural next surface.
