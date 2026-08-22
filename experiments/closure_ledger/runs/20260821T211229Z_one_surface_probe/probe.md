# One-surface probe — one field, one curve, one question

**Question.** read as ONE scalar field on ONE surface rather than two overlaid curves, where do two oppositely-oriented refocusing contributions best span the bulk gap?

**Answer.** not at 'the antipode' in general -- at half a wavelength for a travelling mode, and at the antipode for every ODD zonal mode, which is where the two models disagree; coincident foci cancel exactly in both

**8/8 checks pass.**

| id | check | result |
|----|-------|--------|
| T1 | *** the repair costs nothing *** | PASS |
| T2 | *** coincidence cancels exactly *** | PASS |
| T3 | the optimum is half a wavelength | PASS |
| T4 | *** the antipode is parity-dependent *** | PASS |
| T5 | *** the zonal optimum IS the antipode *** | PASS |
| T6 | a pulse saturates where a mode diverges | PASS |
| T7 | the chord shrinks with frequency | PASS |
| T8 | fixed energy is not fixed amplitude | PASS |

## The repair costs nothing

| | |
|--|--|
| v66 separation vs one-surface field | `2.0` ulp of `R_mid` |
| the opposed field at `α = 0` | `0.0e+00` — exactly zero |
| its span threshold | `inf` |
| the like pair's, unaffected | `0.4908` |

> v66's 'like-signed' separation is the one-surface OPPOSED deformation, and v66's 'opposed' separation is the one-surface LIKE-signed one

## The monochromatic law, and the parity

| mode `m` | `α*` = half a wavelength | at the antipode | verdict |
|--|--|--|--|
| `1` | `1.0000π` | `2.0000` | **maximal** |
| `2` | `0.5000π` | `0.0000` | **cancels** |
| `3` | `0.3333π` | `2.0000` | **maximal** |
| `4` | `0.2500π` | `0.0000` | **cancels** |
| `5` | `0.2000π` | `2.0000` | **maximal** |
| `6` | `0.1667π` | `0.0000` | **cancels** |
| `7` | — | `2.0000` | **maximal** |
| `8` | `0.1250π` | `0.0000` | **cancels** |

> sin(m pi/2) is +-1 for odd m and 0 for even m; on S^3 the same parity appears as Z_n(pi) = (-1)^n, so it is a property of the spectrum and not of the plane-wave reduction

## Where the two models part company

| zonal order `n` | `ω = n+1` | measured `α*` | peak strength | at the antipode | at half a wavelength |
|--|--|--|--|--|--|
| `1` | `2` | `1.0000π` | `2.0000` | `2.00e+00` | `1.4142` |
| `2` | `3` | `0.5002π` | `1.3333` | `1.88e-12` | `1.1547` |
| `3` | `4` | `1.0000π` | `2.0000` | `2.00e+00` | `1.1217` |
| `4` | `5` | `0.2902π` | `1.2500` | `1.13e-12` | `1.1097` |
| `5` | `6` | `1.0000π` | `2.0000` | `2.00e+00` | `1.1039` |
| `6` | `7` | `0.7941π` | `1.2330` | `8.08e-13` | `1.1006` |
| `8` | `9` | `0.1596π` | `1.2266` | `6.28e-13` | `1.0971` |
| `10` | `11` | `0.1303π` | `1.2234` | `2.06e-12` | `1.0954` |

> a zonal harmonic is centred, so the antipode is where two centres coincide with opposite sign; a plane wave has no centre, so only its wavelength sets a scale

> this programme's kernel is n = 1, which is odd, so for it the antipode is both optimal and saturating

## A pulse is not a mode

Pulse width `0.18` rad. The threshold plateaus at `0.2163` (spread `4.2e-04`).

| offset `α/π` | in pulse widths | monochromatic `A_req` | pulse threshold | ratio |
|--|--|--|--|--|
| `0.02` | `0.35` | `4.1387` | `0.6833` | `6.06×` |
| `0.04` | `0.70` | `2.0704` | `0.3655` | `5.66×` |
| `0.06` | `1.05` | `1.3814` | `0.2725` | `5.07×` |
| `0.10` | `1.75` | `0.8310` | `0.2217` | `3.75×` |
| `0.20` | `3.49` | `0.4207` | `0.2156` | `1.95×` |
| `0.40` | `6.98` | `0.2212` | `0.2161` | `1.02×` |
| `0.60` | `10.47` | `0.1607` | `0.2164` | `0.74×` |
| `0.80` | `13.96` | `0.1367` | `0.2165` | `0.63×` |
| `1.00` | `17.45` | `0.1300` | `0.2166` | `0.60×` |

> two localized pulses cancel only while they overlap; a mode fills the whole circle and cancels everywhere

## The chord, and what it costs

| mode `m` | separation | span / `A` | bulk chord | `L/D` | span at fixed energy | spans the gap? |
|--|--|--|--|--|--|--|
| `1` | `1.0000π` | `2.0000` | `2.0000` | `3.846` | `4.0000` | yes |
| `2` | `0.5000π` | `2.0000` | `1.4612` | `2.810` | `2.0000` | yes |
| `3` | `0.3333π` | `2.0000` | `1.0967` | `2.109` | `1.3333` | yes |
| `4` | `0.2500π` | `2.0000` | `0.9037` | `1.738` | `1.0000` | yes |
| `5` | `0.2000π` | `2.0000` | `0.7915` | `1.522` | — | — |
| `6` | `0.1667π` | `2.0000` | `0.7213` | `1.387` | `0.6667` | yes |
| `8` | `0.1250π` | `2.0000` | `0.6421` | `1.235` | `0.5000` | **no** |
| `16` | `0.0625π` | `2.0000` | `0.5534` | `1.064` | `0.2500` | **no** |

The chord falls from `2.0000` to `0.5534`, with the limit the purely radial gap `0.52`.

> low frequency buys a large deformation on a long connection; high frequency buys a short connection with a small deformation

At fixed energy the highest mode that still spans the gap is `m = 6`.

**Not claimed:** no favourable frequency, because that needs an energy normalisation and a packet focusing law that this model does not contain
