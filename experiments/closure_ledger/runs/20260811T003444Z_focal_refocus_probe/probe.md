# The antipodal refocus of a trace-free deformation

_2026-08-11T00:34:44+00:00_

> A faithful **representation** of a spin-2 field in the ℝ³ embedding of a fixed `S²` — not backreaction.

## The focus is a node

| quantity | value |
|---|---:|
| peak density at the focus | 0.096453 |
| amplification | 2.2790× |
| focal time | 2.9897 |
| ring radius | 0.1662 |
| density **on** the antipode | `2.9e-08` of peak |
| conserved invariant drift | `1.7e-13` |

`∫ρ_E dA` is the *kinetic* half and oscillates against the gradient half (swing `0.99`), so the solver's invariant is the conservation check, not that.

## The ring radius is the pulse width

| pulse width | ring radius | ratio | amplification |
|---:|---:|---:|---:|
| 0.24 | 0.2238 | 0.933 | 2.320× |
| 0.18 | 0.1662 | 0.924 | 2.279× |
| 0.12 | 0.1139 | 0.949 | 2.238× |
| 0.09 | 0.0877 | 0.974 | 2.224× |
| 0.06 | 0.0589 | 0.982 | 2.154× |

Mean ratio 0.952, spread 0.058.

## ...and the amplification is *not* a spin-2 effect

| pulse width | tensor | scalar |
|---:|---:|---:|
| 0.24 | 2.320× | 2.067× |
| 0.18 | 2.269× | 2.059× |
| 0.12 | 2.252× | 2.034× |
| 0.09 | 2.188× | 1.951× |
| 0.06 | 2.176× | 1.863× |

Both are `O(1)` and neither runs away as the pulse narrows. What belongs to the spin is the **node** and the **ring**.

## The patch reports the eigenvector, in the limit of a point

| patch radius | aspect ratio | alignment |
|---:|---:|---:|
| 0.24 | 1.1824 | 0.9591 |
| 0.18 | 1.2889 | 0.9874 |
| 0.12 | 1.3944 | 0.9981 |
| 0.08 | 1.4460 | 0.9997 |
| 0.05 | 1.4710 | 1.0000 |

At mid-latitude the alignment is 1.000 at any of these sizes; only near the focus does the patch straddle a sign change.

## Where the area law runs out

| quantity | mid-latitude | focal ring |
|---|---:|---:|
| worst area change | 1.93% | 25.95% |
| max aspect ratio | 1.0467 | 1.3974 |

Fitted residual exponent `ε^2.00` — second order, at display gain `7.87`.

## Verdict

**7/7 checks passed.**

**THE_STRAINS_REFOCUS_ON_A_RING.** THE STRAINS REFOCUS ON A RING, NEVER ON THE POINT. Every principal-strain history on the sphere reconverges at the antipode, and what they build there is an annulus with a hole in the middle.

THE FOCUS IS A NODE. Because h = sin²d·q vanishes at both poles for every q, so does ḣ, and the effective density ∝ ḣ_ab ḣ^ab measured ON the antipode is 2.9e-08 of the peak. There is nowhere for a spin-2 field to sit at its own focus — the same fact that forbids a spin-2 point SOURCE, seen at the other end of the trip.

IT PILES INTO A RING WHOSE RADIUS IS THE PULSE. Across pulse widths 0.24 down to 0.06 the focal ring radius stays at 0.952 × the width, spread 0.058. The hole is not a numerical floor; it is the wave's own length scale.

AND THE AMPLIFICATION IS NOT A SPIN-2 EFFECT. The peak density grows by 2.18–2.32× — finite, and tempting to read as the spin protecting itself from a singularity. It is not. A scalar pulse refocused the same way amplifies by 1.86–2.07×, and neither runs away as the pulse is narrowed, because launch and focus are geometrically the same situation on a sphere. What belongs to the spin is the node and the ring, not the factor.

THE PATCH IS THE PICTURE. A disc of labelled particles carried by the displacement reports the local stretch axis to 1.000 where the field is smooth, and to 1.0000 on the focal ring once it is small enough not to straddle the sign change there — the convergence being the check that the disagreement is the patch's size rather than the construction. At a gain where the first-order statement is the whole statement it changes shape without changing size, holding its area to 1.9e-07.

THE AREA LAW FAILS FIRST, AND HARDEST, AT THE FOCUS. Pushed to the display gain ε = 7.87, the same patch moves its area by 1.9% at mid-latitude and 26.0% on the focal ring — 13× more — while distorting 8.5× harder. The residual scales as ε^2.00, so this is the second-order term carrying the LOCAL GRADIENT of the field: away from the focus that scale is the wavelength, on the focal ring it is the pulse width. The linear description runs out exactly where the wave reconverges.

SCOPE. No singularity forms here and none can: a linear field on a fixed round background, with no backreaction and no bulk crossing rule. What this establishes is the geometry such a rule would have to act on — an annulus of finite radius with a node at its centre — and the amplitude at which the linear description stops being trustworthy. It does not supply the rule.