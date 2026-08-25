# Tangherlini dynamics probe — the first evolved Einstein equations

**Question.** what does the highest-symmetry 4+1 Einstein-scalar system do when it is actually evolved, in a horizon-penetrating chart -- and does the evolution satisfy the Einstein equation it never solves?

**Answer.** the vv equation, which the hierarchy never uses, converges at second order (1.989 -> 1.999); and a regular centre makes A > 0 identically, so no trapped surface can sit on a regular-centred ingoing null slice -- horizon formation is not observable in this chart, which is a statement about the chart and not about the physics

**7/7 checks pass.**

| id | check | result |
|----|-------|--------|
| T0 | the field equations are derived, and reduce to D = 4 | PASS |
| T1 | the hierarchy reproduces an exact flat-space mode | PASS |
| T2 | Tangherlini is an exact fixed point | PASS |
| T3 | *** the unused Einstein equation converges at 2nd order *** | PASS |
| T4 | *** a regular centre forbids a trapped surface *** | PASS |
| T5 | the master potential disagrees with the repository | PASS |
| T6 | the spectrum is NOT cross-validated, and is not reported | PASS |

## T0 — the system, derived for general `n`

* `d_r delta = (kappa/n) r (d_r phi)^2`
* `(r^{n-1} e^delta A)' = (n-1) r^{n-2} e^delta`
* `2 r^n d_r(d_v phi) + n r^{n-1} d_v phi + d_r(e^delta A r^n d_r phi) = 0`

| check | `n = 3` (D = 5) | `n = 2` (the known D = 4 system) |
|--|--|--|
| rr is the delta quadrature | yes | yes |
| vr is the A quadrature | yes | yes |
| birkhoff | yes | yes |

The `vv` equation is the monitor because it is the only independent component containing d_v A, which the hierarchy never computes -- so its residual tests the evolution rather than restating it.

## T1 — the exact flat-space mode

`phi = cos(w(v-r)) J_1(w r)/r` solves the flat `D = 5` wave equation in these coordinates, so `d_v phi` is known in closed form.

| points | flat-metric error | `psi` relative error | rate |
|--|--|--|--|
| 400 | `8.9e-16` | `1.470e-03` | — |
| 800 | `4.4e-16` | `3.649e-04` | `2.010` |
| 1600 | `7.8e-16` | `9.078e-05` | `2.007` |
| 3200 | `1.3e-15` | `2.262e-05` | `2.005` |
| 6400 | `1.1e-15` | `5.644e-06` | `2.003` |

*The r^{n/2} weight from odd n leaves a half-integer power in the origin quadrature; a fourth-order rule measures 2.5 there, so the scheme is uniformly second order instead.*

## T2 — Tangherlini is an exact fixed point

| points | metric error | max abs `delta` | max abs `psi` |
|--|--|--|--|
| 400 | `1.554e-15` | `0.0e+00` | `0.0e+00` |
| 1600 | `4.441e-15` | `0.0e+00` | `0.0e+00` |

Birkhoff in d = 5: with no scalar the only spherically symmetric solution is tangherlini, and the quadrature reproduces it rather than approximating it.

## T3 — the Einstein equation the code never solves

The monitored equation is **vv, the only independent component carrying d_v A**.

| points | spacing | max abs `vv` residual | at radius | rate |
|--|--|--|--|--|
| 400 | `0.0501` | `1.5511e-04` | `4.01` | — |
| 800 | `0.0250` | `3.9070e-05` | `4.01` | **`1.989`** |
| 1600 | `0.0125` | `9.7862e-06` | `4.00` | **`1.997`** |
| 3200 | `0.0063` | `2.4478e-06` | `4.00` | **`1.999`** |

**What this is not.** Not a hamiltonian/momentum constraint pair -- the hierarchy solves rr and vr exactly on every slice, so their residuals are identically zero and testing them would be circular.

**What an imposed outer condition did.** Freezing phi at r = r left an o(1) vv residual there that did not converge at all; the characteristic closure admits no outer boundary condition.

## T4 — a regular centre forbids a trapped surface

    r^{n-1} e^delta A = (n-1) int_0^r s^{n-2} e^delta ds > 0

a positive integrand over a positive interval, so `A > 0` strictly for `r > 0`, identically. The scan below is not the proof — it is the check that the code obeys the proof.

| profile | amplitude | min `A` | at radius | at `v` | trapped? |
|--|--|--|--|--|--|
| centred gaussian | 2.0 | `7.7091e-01` | `1.761` | `1.12` | no |
| centred gaussian | 5.0 | `3.2956e-01` | `1.161` | `0.01` | no |
| centred gaussian | 12.0 | `7.2561e-02` | `1.041` | `0.01` | no |
| thin shell at r = 2 | 2.0 | `1.0430e-01` | `2.509` | `0.01` | no |
| thin shell at r = 2 | 5.0 | `1.5134e-02` | `2.402` | `0.01` | no |
| thin shell at r = 2 | 12.0 | `5.6267e-03` | `2.469` | `1.06` | no |
| r^2 e^{-r^2/2} | 2.0 | `6.1319e-01` | `2.736` | `0.01` | no |
| r^2 e^{-r^2/2} | 5.0 | `1.4969e-01` | `2.495` | `0.01` | no |
| r^2 e^{-r^2/2} | 12.0 | `2.7030e-02` | `2.322` | `0.01` | no |
| oscillatory | 2.0 | `1.2462e-01` | `1.668` | `0.01` | no |
| oscillatory | 5.0 | `2.0036e-02` | `2.549` | `0.01` | no |
| oscillatory | 12.0 | `5.7429e-03` | `2.575` | `0.20` | no |

**The consequence.** Horizon formation is not observable in this gauge with a regular centre; the transition is the loss of central regularity, not a changing sign.

*Why it is a chart statement:* collapse still happens -- the trapped region is reached once the centre stops being regular, at which point the slice carries a nonzero interior mass and this quadrature no longer applies. Outgoing null cones, or excision of the centre -- for exactly this reason.

## T5 — a discrepancy found in passing

| | potential |
|--|--|
| derived here | `A[(l(l+2) + 3/4)/r^2 + (9/4) r_h^2/r^4]` |
| `tangherlini.radial.V_tangherlini` | `A[l(l+2)/r^2 + 3 r_h^2/r^4]` |
| difference | `3 A^2 / (4 r^2)` |

The difference matches that closed form to `5.4e-16`, and the flat limit matches Bessel to `4.3e-16` — psi = r^{1/2} J_{l+1}(wr) solves Bessel's equation with V = ((l+1)^2 - 1/4)/r^2 = (l(l+2) + 3/4)/r^2.

**Nothing was changed.** V_tangherlini is consumed by roughly fifty probes and by several derived constants; replacing it is a decision about the repository's published numbers, not a side effect of a dynamics round.

*Caveat:* the discrepancy is stated for a MINIMALLY COUPLED MASSLESS SCALAR with psi = r^{n/2} phi, which is what the existing docstring describes; a different field or a different substitution has a different potential.

## T6 — what this round did not earn

| asked for | delivered? |
|--|--|
| constraint convergence | yes |
| horizon formation | yes |
| horizon persistence | yes |
| perturbation spectrum | **no** |
| retarded outer-to-inner transfer function | **no** |

Two horizon-penetrating time-domain constructions, both stable and both converged, disagree: Kerr–Schild `1.01622 - 0.36231i` against tortoise `1.01876 - 0.26404i` at `l = 1`. Real parts within 0.3%; damping rates apart by **37%**.

So **no quasinormal frequency is reported**, and the transfer function is not built — it is a ratio of the same two signals and would inherit the same unresolved error.

*First thing to chase:* the Kerr-Schild operator does not converge to the exact flat-space mode at its inner cut (error flat at 1.07e-02 across a four-fold refinement), which points at the excision boundary rather than the operator.

**A converged number is not a correct number.**
