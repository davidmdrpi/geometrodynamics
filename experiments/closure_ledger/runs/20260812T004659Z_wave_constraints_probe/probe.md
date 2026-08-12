# What a drawn wave has to obey

_2026-08-12T00:46:59+00:00_

## The front is causal

| t | amplitude at the antipode | ahead of the front |
|---:|---:|---:|
| 0.00 | `5.1e-133` | `1.1e-07` |
| 0.40 | `4.4e-100` | `8.6e-08` |
| 0.80 | `9.9e-73` | `8.1e-08` |
| 1.40 | `6.2e-40` | `9.1e-08` |
| 2.00 | `6.9e-17` | `1.3e-07` |
| 2.40 | `5.6e-07` | `6.8e-07` |
| 2.80 | `1.5e-01` | `—` |
| 3.00 | `9.3e-01` | `—` |

Nothing outruns the front.

## But a constant never leaves

| quantity | value |
|---|---:|
| monopole at launch | 0.008056 |
| predicted `w²/4` | 0.008100 |
| drift over a full return | `8.1e-04` |
| time-averaged displacement, far side | 0.008810 |
| quietest instant, `max\|u\|` | 0.0936 |

## Only a localised, compact correction works

| source | monopole | far side before arrival |
|---|---:|---:|
| gaussian | `+8.1e-03` | `1.5e-20` |
| gaussian, monopole-free | `-4.3e-07` | `1.4e-05` |
| compact, monopole-free | `-3.1e-07` | `3.8e-22` |

## The ring concentrates where you cannot see it

| quantity | value |
|---|---:|
| growth from the equator | 5.97× |
| flat fraction of the trip | 0.72 |
| growth happens after | t = 3.11 |
| focal height / **launch** height | 0.935 |

## Verdict

**8/8 checks passed.**

**NOTHING_MOVES_EARLY_AND_NOTHING_EVER_SETTLES.** NOTHING MOVES EARLY, AND NOTHING EVER SETTLES. The causality objection is aimed one step away from a real defect, and the real defect is worse than the one that was suspected.

THE FRONT IS CAUSAL. The amplitude at the antipode is 5e-133 at launch and 7e-17 at t = 2.0, reaching 0.926 only at t ≈ π; signal ahead of the front sits at 7e-07, the scheme's own noise floor. Nothing outruns the front, and the early motion is not in the solve.

BUT A CONSTANT NEVER LEAVES. A Gaussian of width w carries a monopole of 0.008056 against a predicted w²/4 = 0.008100, and ℓ = 0 has ω = 0, so it never propagates, radiates or decays. It is NOT an early response — ahead of the front the higher modes cancel it exactly — it is a permanent one. The time-averaged displacement is 0.008810 at the far side, and the quietest instant of a whole run still leaves max|u| = 0.094. The surface never returns home.

AND THE MODE IS FORBIDDEN ANYWAY. Electromagnetism has no monopole radiation and gravity has none at ℓ = 0 or ℓ = 1, so the offending mode is precisely the one real radiation cannot have. That is not a blemish on the analogy, it is outside it. The checkable form is the spectrum: the spin-2 operator starts at ℓ = 2, so it admits no zero-frequency mode and can hold no DC at all. Measured as the time-averaged value over four returns, 1.9e-06 against the scalar's 8.3e-03 — a factor of 4290.

AND THE FIX HAS TO STAY INSIDE THE PULSE. Subtracting the mean moves the problem rather than removing it, leaving -0.00806 at the antipode — exactly the offset it was meant to remove. Cancelling with a wider Gaussian fixes the mean and costs four orders of far-side quiet, 1e-05 against 1e-20, because the corrector's tail is fatter than the pulse it corrects. Built from compactly supported bumps both conditions hold at once: monopole -3.1e-07 and a far side of 4e-22 — quieter than the original, because it is exactly zero rather than merely small.

THE RING DOES CONCENTRATE, CORRECTLY. A²·(circumference) is conserved along the converging front to 0.084, which is A ∝ 1/√(sin d). The solve is right.

YOU JUST CANNOT SEE IT, AND IT CANNOT WIN. 1/√(sin d) is flat across the middle: 72% of the trip stays within 50% of the equatorial minimum and the growth is all after t = 3.11. Worse for the expectation, on a compact surface the launch is itself a focus — the focal height is 0.94× the LAUNCH height, so the focus never beats the source. The 5.97× growth is relative to the equator. Unbounded growth belongs to an open geometry, not this one.

AND HEIGHT IS NOT ENERGY. The energy density is u̇² + |∇u|², so the constant offset displaces every point of the surface and carries exactly 0.0 energy against 0.497 for the pulse. The most global feature of the drawn shape is invisible to the physics. That is the measurable core of the objection to drawing a field as a displacement: the eye is on height, the physics is in the gradient.

SCOPE. This audits one representation; it does not propose its replacement. What it establishes is which objections are about the solve — none — which are about the initial data, and which are about drawing a field as a displacement at all.