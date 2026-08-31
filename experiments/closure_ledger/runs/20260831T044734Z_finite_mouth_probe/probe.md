# The finite-mouth scalar-flat handle

**6/6 pre-registered predictions reproduced.**

Every number below was frozen in `docs/finite_mouth_prereg.md` **before this module existed**. P1 and P4 are falsification attempts rather than confirmations.

> ## The construction

> One assumption: the closed `S³` universe is the totally geodesic equator of a round `S⁴_R` spatial bulk. It adds no scale and no shape parameter. Demanding `⁴R = 0` in the throat forces `f = √(s²+b²)`, and Darmois matching to the exterior then fixes **both** constants:

> ```
> b = R sin²a ,   S = R sin a cos a ,   L = R sin 2a
> ```

> There is no tube area, neck radius or throat length left to choose — which is exactly the freedom that carried the answer through PR #263–#265.

---

## P1 — the geometry is parameter-free

| `R` | `a` | `b` | `S` | `L` | matching residual |
|--|--|--|--|--|--|
| `1` | `0.3` | `0.087332193` | `0.282321237` | `0.564642473` | `1.1e-16` |
| `1` | `0.8` | `0.514599761` | `0.499786802` | `0.999573603` | `0.0e+00` |
| `2.5` | `0.3` | `0.218330481` | `0.705803092` | `1.411606183` | `1.1e-16` |
| `1` | `0.05` | `0.002497917` | `0.049916708` | `0.099833417` | `0.0e+00` |

**Falsifier.** Holding the areal radius correct and moving `b` by ±5% leaves a slope error of `4.70e-03`: there is no second matching pair.

> Darmois matching is two conditions (areal radius and its normal derivative) on two constants (b, S), so it has a unique solution. Holding the radius right and moving b by 5% leaves a slope error of order 1e-2: there is no second matching pair, and therefore no tube area, neck radius or length left to choose. PR #263-#265 spent three rounds discovering that a chosen area was carrying the answer.

## P2 — the quasi-local mass does not jump

| `R` | `a` | `μ` inside | `μ_ext(a)` | jump |
|--|--|--|--|--|
| `1` | `0.3` | `0.007626912` | `0.007626912` | `1.6e-17` |
| `1` | `0.8` | `0.264812914` | `0.264812914` | `5.6e-17` |
| `2.5` | `0.3` | `0.047668199` | `0.047668199` | `1.6e-16` |
| `1` | `0.05` | `0.000006240` | `0.000006240` | `2.1e-19` |

> mu = f^2(1-f'^2) equals b^2 everywhere inside and R^2 sin^4 chi outside; they agree at chi = a because b = R sin^2 a. The seam is smooth exactly when the quasi-local mass parameter is continuous -- the 5D lift of PR #265's Hawking-mass matching.

## P3 — no Israel shell

| `R` | `a` | `[h]` | `[K]` | `p_n` inside | `p_n` outside |
|--|--|--|--|--|--|
| `1` | `0.3` | `5.6e-17` | `1.3e-15` | `-3.000000` | `-3.000000` |
| `1` | `0.8` | `0.0e+00` | `0.0e+00` | `-3.000000` | `-3.000000` |
| `2.5` | `0.3` | `0.0e+00` | `2.2e-16` | `-0.480000` | `-0.480000` |
| `1` | `0.05` | `0.0e+00` | `0.0e+00` | `-3.000000` | `-3.000000` |

> [h_ab] = [K_ab] = 0 gives S_ab = 0. The normal pressure agreeing at -3/(8 pi G_5 R^2) on both sides is the Gauss-Codazzi constraint that must hold when no shell is present, and it is computed from the two sides separately.

> f'' = b^2/f^3 inside and -R sin chi outside do NOT agree, so the geometry is C^1 and not C^2. That is a finite STEP in bulk stress, not a delta function: the Israel layer depends on [K_ab], which vanishes. Tangential pressure and density jump, which is allowed because they involve second normal derivatives.

## P4 — the neck NEC price, attacked

Predicted `8πG₅(ρ+p_s)|₀ = −3/b² = -393.343997824` for **every** smooth lapse with `N(0) > 0`. Seven hostile lapses:

| lapse | `(ρ+p_s)` at the neck | deviation |
|--|--|--|
| `N = 1 (ultrastatic)` | `-393.343997824` | `5.7e-14` |
| `N = 1 + 0.7 s (asymmetric)` | `-393.343997824` | `5.7e-14` |
| `N = 1 + 3 s^2` | `-393.343997824` | `5.7e-14` |
| `N = 1 - 2 s^2 + 5 s^3` | `-393.343997824` | `5.7e-14` |
| `N = exp(4 s)` | `-393.343997824` | `5.7e-14` |
| `N = 2 + cos(9 s)` | `-393.343997824` | `5.7e-14` |
| `N = 0.05 + 8 s^2 (nearly null)` | `-393.343997824` | `5.7e-14` |

**None evades it.** Worst deviation `5.7e-14`.

> The lapse enters p_s only through 3 f'N'/(fN), and f'(0) = 0 is what MAKES s = 0 a neck. So the term vanishes there whatever N' does. This is stronger than the proposal stated: it needs no reflection symmetry, which is why an asymmetric and an oscillating lapse give the identical value.

> **N(0) = 0 -- the Tangherlini horizon branch, which is vacuum and non-traversable. Smooth AND traversable implies radial NEC violation at the neck.**

## P5 — the admittance, against an independent solver

| `ℓ` | closed-form diagonal | off-diagonal | relative error |
|--|--|--|--|
| 0 | `0.078793811` | `-0.078793811` | `2.6e-09` |
| 1 | `1.800864275` | `-0.003597716` | `1.4e-07` |
| 2 | `3.524730783` | `-0.000123268` | `2.8e-07` |
| 3 | `5.248599166` | `-0.000003754` | `4.8e-07` |
| 5 | `8.696335933` | `-0.000000003` | `1.0e-06` |

Second-order convergence of the independent solve at `ℓ = 2`:

| steps | relative error | ratio |
|--|--|--|
| `1000` | `4.47e-06` | — |
| `2000` | `1.12e-06` | `4.00` |
| `4000` | `2.79e-07` | `4.00` |

Monopole: `G = 0.078793811` = `π²R²sin⁴a/cos a`, and `2π²/I₃` agrees to `1.4e-17` with `I₃ = 250.517249267`.

Row sums vanish to `1.8e-16`. Row sums vanish exactly: a constant is in the kernel, so there is no static monopole shunt through the handle.

> The BVP solve never uses the sinh/cosh reduction, so this is a regression and not an identity. A shooting basis was tried first and rejected: its two solutions span e^{+-kx} over a rapidity 2kX ~ 23 at l = 5, so the far boundary condition costs ten digits and the error grew under refinement instead of shrinking. The conservative tridiagonal form converges at second order at every l tested.

## P6 — one spatial geometry, two lapses

| branch | `N(0)` | stress | causal character |
|--|--|--|--|
| Tangherlini | `0.0e+00` | vacuum (radial NEC `5.6e-07`) | horizon, non-traversable |
| ultrastatic | `1.0` | anisotropic, NEC-violating | traversable |

`⁴R` vanishes on the shared profile to `2.3e-13`, and `N_vac²` reproduces the Tangherlini `F(r)` to `3.3e-16`.

> N_vac^2 = |s|^2/(s^2+b^2) is exactly the Tangherlini F(r) = 1 - b^2/r^2 with r^2 = s^2+b^2, on the SAME spatial profile that the ultrastatic branch uses. The repository's vacuum throat and transaction throat are one spatial geometry with two lapses; the entire physical fork is the number N(0).

---

## What this does not settle

The ANEC cost of the finite ultrastatic throat is `-53.352064/(8πG₅)` at `R = 1, a = 0.3`, diverging as `−3π/(2Ra²)` for small mouths — so the point-mouth limit is singularly expensive. Whether any classical BAM degree of freedom can supply that stress, and therefore whether `N(0) > 0` is available at all, is untouched here. If none can, the geometry collapses onto the Tangherlini branch and the MTY transaction mechanism is unavailable.

The discrete BAM identification is deliberately absent: `Φ_spatial`, `(−1)^ℓ`, `η_orientation`, `η_wrap` and `U_spin` remain five separate objects, none folded into the metric.
