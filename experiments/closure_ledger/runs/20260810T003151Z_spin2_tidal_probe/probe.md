# Spin 0 against spin 2 on the same S²

_2026-08-10T00:31:51+00:00_

> A spin-2 **field on a fixed** `S²` — the tensor analogue of the scalar wave, not linearised GR on a spacetime.

## The solver, against exact modes

| `q` | `ℓ` | `ω` | shape error | invariant drift |
|---|---:|---:|---:|---:|
| `1` | 2 | 2.4495 | `3.5e-08` | `3.8e-11` |
| `cos d` | 3 | 3.4641 | `5.6e-05` | `7.0e-13` |
| `7cos²d − 1` | 4 | 4.4721 | `1.3e-04` | `1.2e-13` |

After 10 periods.

## A spin-2 field cannot sit on a pole

| quantity | value |
|---|---:|
| tensor peak at distance | 2.9439 |
| — i.e. a ring of radius | 0.1977 |
| amplitude *on* the antipode | `2.2e-06` |
| **scalar** peak at distance | 3.1416 |
| lowest multipole | ℓ = 2 |

## Pure shear, not breathing

| quantity | value |
|---|---:|
| trace | `0.0e+00` |
| area ratio | 1.0000000117 |
| first-order area change | `1.17e-08` |
| `ε²h²/2` prediction | `1.17e-08` |

## Spin weight 2

| rotation of the frame | strain |
|---|---:|
| 180° | identical, `1.1e-15` |
| 90° | inverted, 2.00× amplitude |

## The caustic: a quarter turn, not a flip

| the outbound front is the inbound one… | correlation |
|---|---:|
| unchanged | +0.3542 |
| inverted | -0.3542 |
| phase +90 | -0.8165 |
| phase -90 | +0.8165 |

Best match: **phase_-90** — the Gouy shift, Maslov index 1.

## The round trip, where the axes do swap

| | correlation |
|---|---:|
| `h(2π)` with `−h(0)` | +0.9974 |
| `h(2π)` with `h(0)` | -0.9974 |
| `h(4π)` with `h(0)` | +0.9877 |
| inversion residual | 0.0064 |

## Verdict

**7/7 checks passed.**

**A_TENSOR_WAVE_CANNOT_SIT_ON_ITS_OWN_FOCUS.** A TENSOR WAVE CANNOT SIT ON ITS OWN FOCUS. Replacing the scalar with a spin-2 field on the same sphere changes the picture structurally, not cosmetically, and the changes are measured rather than drawn.

THE SOLVER IS EXACT WHERE IT CAN BE CHECKED. Three closed-form modes — q = 1 at ℓ = 2, q = cos d at ℓ = 3, q = 7cos²d − 1 at ℓ = 4 — keep their shape to 1.3e-04 after 10 periods at ω = √6, √12, √20, with the invariant conserved to round-off.

NO AMPLITUDE AT A POLE. h = sin²d·q vanishes at both poles for every q, so a spin-2 field cannot be nonzero where the frame degenerates. At the focus it is therefore a RING of radius 0.198 about the antipode, with 2.2e-06 on the antipode itself — while the scalar, on the same clock, piles up exactly there at d = 3.1416. The same fact at the other end says the smallest source a spin-2 field admits is a ring: there is no such thing as a point source of tidal shear.

PURE SHEAR, NOT BREATHING. The tensor is symmetric and trace-free (0.0e+00), so the ellipse it makes has semi-axes 1 ± εh and its first-order area change vanishes identically: measured 1.17e-08 against the second-order prediction 1.17e-08. A radial height perturbation changes area at first order — that is the difference you can see.

SPIN WEIGHT 2, DIRECTLY. The strain pattern is identical after a 180° rotation of the frame (1.1e-15) and inverted after 90° (2.00× the amplitude). That factor of two in the angle is the spin weight, visible without any decomposition.

THE CAUSTIC IS A QUARTER TURN, NOT A FLIP. The obvious guess — that passing the antipodal focus swaps the stretch and compression axes — is not what happens. Correlating the outbound waveform against all four candidates gives +0.816 for a −90° phase shift against -0.354 for an inversion: the outbound front is the HILBERT TRANSFORM of the inbound one. That is the Gouy shift, Maslov index 1, and it is a property of passing through a focus rather than of the spin — the scalar has it too.

THE ROUND TRIP IS WHERE THE AXES SWAP. Two focal passages, the antipode and then home, compose to a half turn: at t = 2π the field is minus what it started as, correlation 0.9974, and after two round trips it is itself again (0.9877). The residual 0.0064 is the dispersion left over because ω_ℓ = √(ℓ(ℓ+1)) only approaches ℓ + ½.

SCOPE. A spin-2 field on a FIXED S², not a solution of the 4-D linearised Einstein equations — 2+1 gravity has no propagating tensor modes at all. What is faithful is the polarisation structure: two components, spin weight 2, ℓ ≥ 2 only, area-preserving shear, and the behaviour at a focus.