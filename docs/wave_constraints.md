# What a drawn wave has to obey

> An audit of the height-field representation used by the circle-slice work,
> against constraints the physics it stands in for actually imposes.

Three objections prompted this. They do not all land the same way, and being
precise about which is which is the whole value of the exercise.

## 1. "The antipode moves before the wavefront arrives"

**This one does not land on the solve.** The amplitude at the antipode:

| t | amplitude at the antipode | ahead of the front |
|---:|---:|---:|
| 0.00 | `5.1e-133` | `1.1e-07` |
| 0.80 | `9.9e-73` | `8.1e-08` |
| 1.40 | `6.2e-40` | `9.1e-08` |
| 2.00 | `6.9e-17` | `1.3e-07` |
| 2.80 | `1.5e-01` | — |
| 3.00 | `9.3e-01` | — |

It climbs 130 orders of magnitude in step with `t → π`. Signal ahead of the
front sits at `~1e-07`, the leapfrog's own noise floor. **Nothing outruns the
front.**

## 2. But a constant never leaves — and that *is* a defect

On a closed surface `d²/dt² ∫u dA = ∫Δu dA = 0`, so the mean field is at worst
linear in time. `throat_wavefront._outgoing_velocity` already kills the linear
term. What it cannot kill is the **constant**:

| quantity | value |
|---|---:|
| monopole at launch | `0.008056` |
| predicted `w²/4` | `0.008100` |
| drift over a full return | `8.1e-04` |
| time-averaged displacement, far side | `0.008810` |
| quietest instant, `max\|u\|` | `0.0936` |

`ℓ = 0` has `ω = 0`, so this offset never propagates, radiates or decays.

It is **not** an early response — ahead of the front the higher-`ℓ` modes cancel
it exactly, which is why the far side really is dark. It is a **permanent**
one. Every point of the sphere carries a time-averaged displacement of `w²/4`,
and the quietest instant of a whole run still leaves the surface visibly
deformed.

**So the objection lands one step over from where it was aimed: nothing moves
early, and nothing ever stops moving.**

### And the mode is forbidden anyway

Electromagnetism has no monopole radiation; gravity has none at `ℓ = 0` or
`ℓ = 1`. The offending mode is precisely the one real radiation cannot have —
not a blemish on the analogy but outside it.

The checkable form is the spectrum, not `∫h dA`: `h` is a tensor component in a
frame, not a scalar, so its plain sphere-average has no reason to vanish and
measuring *that* would be measuring the wrong thing. What matters is that the
spin-2 operator's spectrum starts at `ℓ = 2`, so it admits no zero-frequency
mode and can hold no DC at all. Measured as the time-averaged value over four
returns: **`2e-06` for the spin-2 field against `8.3e-03` for the scalar.**

### The fix has to stay inside the pulse

| source | monopole | far side before arrival |
|---|---:|---:|
| gaussian | `+8.1e-03` | `1.5e-20` |
| gaussian, monopole-free | `-4.3e-07` | `1.4e-05` |
| compact, monopole-free | `-3.1e-07` | `3.8e-22` |

Subtracting the mean is the obvious repair and is worse than the disease: it
leaves `−w²/4` at the antipode, *exactly* the offset it removed, smeared over the
whole far side.

Cancelling with a wider Gaussian fixes the mean and costs four orders of
far-side quiet — the corrector's tail is fatter than the pulse it corrects, and
that is analytic, not numerical (it is resolution-independent).

Built from **compactly supported** bumps, `(1 − x²)⁴`, both conditions hold at
once and the far side ends up *quieter than the original*: exactly zero rather
than merely small. A Gaussian can never say that.

## 3. "Height should increase as the ring compresses"

**The physics is right and the picture hides it.** `A²·(circumference)` is
conserved along the converging front to `8.4%` — that is `A ∝ 1/√(sin d)`, exact
ring bookkeeping.

Two things hide it:

| quantity | value |
|---|---:|
| growth from the equator | `5.97×` |
| flat fraction of the trip | `0.72` |
| growth happens after | `t = 3.11` (of `π = 3.14`) |
| focal height / **launch** height | `0.935` |

- `1/√(sin d)` is flat across the middle and only diverges within a pulse width
  of the antipode, so `72%` of the trip shows nothing.
- **On a compact surface the launch is itself a focus.** The source is a point
  and so is the antipode, so the focal height does not merely fail to exceed the
  launch height — it comes in *below* it, at `0.935×`. The `5.97×` growth is
  relative to the *equator*.

An expectation of unbounded growth belongs to an open geometry, where the wave
starts on a finite front and converges to a point. It does not belong here.

## 4. "A deforming surface is not the right representation"

**This is right, and the sharpest reason is measurable.** The energy density is
`u̇² + |∇u|²`, so a *constant* offset:

| | |
|---|---:|
| displaces every point of the surface | yes |
| energy it carries | `0.000` |
| energy of the pulse, for comparison | `0.497` |

The most global feature of the drawn shape is invisible to the physics. The eye
goes to height; the energy is in the gradient. No choice of source repairs that
— it is a property of drawing a field as a displacement.

Together with the two earlier rounds this is now a list of four things the
height representation cannot do:

1. it cannot wind (`docs/circle_slice_bulk.md`, `docs/seam_scale.md`)
2. it cannot self-intersect (`docs/ring_and_fold.md`)
3. it holds a DC offset the geometry cannot shed — fixable in the source
4. its most salient feature carries none of the energy — not fixable

## Scope

This audits one representation; it does not propose its replacement. What it
establishes is which objections are about the solve (**none**), which are about
the initial data (the monopole — fixable, and fixed here), and which are about
drawing a field as a displacement at all.

## Reproduce

```bash
python -m experiments.closure_ledger.wave_constraints_probe
# Verdict: NOTHING_MOVES_EARLY_AND_NOTHING_EVER_SETTLES  (8/8)

python -m pytest tests/test_viz_wave_constraints.py -q       # 16 passed
```
