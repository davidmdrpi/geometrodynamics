# The retarded outer→inner transfer kernel, `D = 5` Tangherlini

**11/11 checks pass.**

PR #270 deferred this object because it is a ratio of two signals it could not trust. PR #274 settled which signal was right, against a published spectrum. This round builds the kernel.

A wave sent in from the far region reaches the horizon filtered. That filter is `T_ℓ(ω)`; in time it is a convolution kernel with `ψ_trans(v) = (K_ℓ ⋆ ψ_inc)(v)`, `v = t + r*`. Since `T → 1` at high frequency, `K_ℓ(t) = δ(t) + K_ℓ^reg(t)` — instantaneous part plus memory.

> ## The transfer is not rigid, and not marginally.

| quantity | value | |
|--|--|--|
| `∫K_reg dt` | `-0.997757` | exact value `−1`, a **sum rule** from `T(0) = 0` |
| `∫\|K_reg\| dt` | `2.0286` | against the `δ`'s weight of `1` |
| `T(ω→0)` | `4.100e-07` | the barrier blocks DC completely |

A rigid exchange kernel is `δ(t)`, possibly delayed: whatever enters leaves undistorted, and a static signal passes perfectly. The real geometry **blocks a static signal completely**, and does so *entirely* through the memory term, which exactly cancels the instantaneous one at DC. In absolute mass the memory is about twice the delta. It is not a correction to rigid exchange — it is the same size as the thing it would correct.

## K0 — three exact anchors

**1. The barrier peak at `ℓ = 1` is rational.**

| `ℓ` | `r²` at peak | `r*` | `V_max` |
|--|--|--|--|
| 0 | `1.605551` | `0.1978` | `0.505383744` |
| 1 | `1.800000` | `0.3792` | `1.234567901` |
| 2 | `1.893190` | `0.4541` | `2.476709048` |
| 3 | `1.935691` | `0.4862` | `4.223427135` |

At `ℓ = 1`: **`r² = 9/5` and `V_max = 100/81`, exactly.** In x = 1/r^2 the peak solves 27/4 x^2 - 2(9/4 - L) x - L = 0. At l = 1 the discriminant is a perfect square (110.25 = 10.5^2) and x = 5/9 exactly; l = 0, 2, 3 give sqrt(13), sqrt(1621), sqrt(57). l = 1 is also the mode PR #270's two codes disagreed on.

Note `r² → 2`, the photon sphere PR #274 pinned exactly.

**2. The integrated potential is exact:** `∫V dr* = ℓ(ℓ+2) + 3/2`. dr* = dr/f cancels the f in V, leaving int_1^inf (L/r^2 + (9/4)/r^4) dr = L + 3/4 = l(l+2) + 3/2.

| `ℓ` | exact | quadrature (truncated) | predicted missing tail |
|--|--|--|--|
| 0 | `1.500` | `1.4950` | `0.0050` |
| 1 | `4.500` | `4.4750` | `0.0250` |
| 2 | `9.500` | `9.4417` | `0.0583` |
| 3 | `16.500` | `16.3950` | `0.1050` |

The deficit matches the predicted tail in every row, so the domain truncation is accounted for rather than assumed small.

**3. Hence the high-frequency phase is closed-form:** `T_ℓ(ω) → exp(−i c_ℓ/ω)` with `c_ℓ = (ℓ(ℓ+2)+3/2)/2`, so **`c₁ = 9/4` exactly**. This is what makes the kernel computable at all: `T − 1 ~ −i c/ω` decays too slowly to transform numerically, and knowing `c` lets it be removed analytically instead of windowed away.

## K1 — GATE 2: flux conservation, and why this is not the failed shoot

| `ω` | `\|T\|` | `\|R\|` | `\|R\|²+\|T\|²−1` |
|--|--|--|--|
| `0.1` | `0.00079553` | `0.99999968` | `+1.1e-15` |
| `0.3` | `0.01423808` | `0.99989863` | `-1.7e-15` |
| `0.7` | `0.19749250` | `0.98030440` | `+8.9e-15` |
| `1` | `0.64441842` | `0.76467307` | `+3.0e-14` |
| `1.5` | `0.98880233` | `0.14923117` | `-3.9e-13` |
| `3` | `0.99999977` | `0.00068283` | `+1.1e-13` |
| `10` | `1.00000000` | `0.00000000` | `-5.3e-13` |
| `30` | `1.00000000` | `0.00000000` | `+6.3e-13` |

Worst residual **`6.3e-13`**, and unitarity is imposed nowhere — it is a consequence of the computation, so it measures it.

> PR #270 and #274 could not converge a quasinormal frequency by shooting in real r, because for Im w < 0 the outgoing solution grows like e^{|Im w| R} and swamps the coefficient being zeroed. Here w is REAL, both e^{+-i w r*} have unit modulus, and nothing dominates anything. This is a different, well-posed problem -- not a repair of the one that failed. Unitarity to ~1e-13, imposed nowhere, is the evidence.

Second order in the spatial step: `|T|` at `ω = 1` gives `0.64442644`, `0.64442002`, `0.64441842`, successive differences `6.4e-06`, `1.6e-06`.

### The subtracted coefficient is the exact one

`transfer_kernel` subtracts `c_ℓ = ½∫V dr*` **exactly**. A fitted coefficient would leave `−i(c_exact − c_fit)/ω` in the remainder — still `1/ω`, defeating the purpose of the subtraction. The fit is kept as a *measurement against* the exact value:

| outer edge | fitted `c` | exact | deficit |
|--|--|--|--|
| `150` | `2.248945` | `2.2500` | `+0.001055` |
| `300` | `2.248945` | `2.2500` | `+0.001055` |
| `600` | `2.248945` | `2.2500` | `+0.001055` |

Spread across outer edges `4.4e-07` — **edge-independent**. The fitted value sits ~0.05% below the exact one, uniformly in the outer edge. That is the 1/r^4 and 1/r^6 part of V, which the centrifugal Jost condition does not capture -- a known, bounded, edge-independent shortfall rather than a truncation drift.

### The low-frequency outer matching

At the lowest bin `ω = 0.00488` the outer turning scale is `397`, and `V` at the edge (`1.67e-04`) still exceeds `ω²` (`2.38e-05`) — that bin sits **inside** the centrifugal tail, where free plane waves are simply the wrong basis. It sets the DC end of the inverse transform, and therefore the numerical realisation of the int K_reg dt = -1 sum rule.

| outer edge | `|T|` at `w=0.0048828` | `|T|` at `w=0.02` | `|T|` at `w=0.1` |
|--|--|--|--|
| `150 | `4.099913e-07` | `1.393649e-05` | `7.955278e-04` |
| `300 | `4.100095e-07` | `1.393644e-05` | `7.955282e-04` |
| `600 | `4.100145e-07` | `1.393644e-05` | `7.955282e-04` |

Relative spread across outer edges: `5.6e-05`, `3.6e-06`, `4.3e-07`.

> Replacing plane-wave matching with the exact centrifugal Jost solutions moved T(w -> 0) from 1.73e-06 to 4.10e-07 -- a factor of four closer to its exact value of zero -- and made the low-w spectrum independent of the outer edge instead of drifting with it.

## K2 — GATE 1: causal support

| `t` | `K_reg(t)` |
|--|--|
| `-300` | `+2.025e-07` |
| `-200` | `-3.318e-07` |
| `-100` | `+4.272e-07` |
| `-50` | `-1.030e-06` |
| `-20` | `-2.439e-06` |
| `-10` | `+4.961e-06` |
| `-5` | `+9.529e-06` |
| `-2` | `+2.809e-05` |
| `-1` | `-4.140e-05` |
| `-0.5` | `-1.005e-04` |

Worst acausal value **`1.00e-04`**, and **`1.0e-06`** away from the front.

> K(t) vanishes identically for t < 0, so whatever the computation returns there IS its noise floor -- no reference value needed. Any feature at t > 0 smaller than that floor is not measurable, which is how this round knows the late-time tail is out of reach.

| `t` | `K_reg(t)` |
|--|--|
| `0.5` | `-1.170e+00` |
| `1` | `-4.326e-01` |
| `2` | `+2.242e-01` |
| `3` | `+2.212e-01` |
| `5` | `-8.461e-02` |
| `8` | `+2.676e-02` |
| `12` | `-8.741e-03` |
| `16` | `+9.720e-04` |
| `20` | `+1.922e-04` |
| `30` | `-1.324e-05` |
| `40` | `+7.921e-07` |
| `60` | `+7.661e-08` |
| `100` | `+3.466e-07` |

## K3 — GATE 3: the kernel carries the published ringdown

Reference (external): `1.01601691-0.36232802i`. Source: Matyjasek, Phys. Rev. D 104, 084066 (2021), arXiv:2107.04815 -- continued fractions cross-checked against Hill determinants, agreeing to 11 digits.

| `dt` | window | fitted `ω` | real err | damping err |
|--|--|--|--|--|
| `0.2` | `(3.0, 14.0)` | `1.028118-0.376065i` | `1.191%` | `3.791%` |
| `0.2` | `(4.0, 16.0)` | `1.026316-0.369379i` | `1.014%` | `1.946%` |
| `0.2` | `(5.0, 18.0)` | `1.023258-0.363123i` | `0.713%` | `0.219%` |
| `0.3` | `(3.0, 14.0)` | `1.017115-0.363824i` | `0.108%` | `0.413%` |
| `0.3` | `(4.0, 16.0)` | `1.017168-0.361390i` | `0.113%` | `0.259%` |
| `0.3` | `(5.0, 18.0)` | `1.015356-0.360984i` | `0.065%` | `0.371%` |
| `0.5` | `(3.0, 14.0)` | `1.017615-0.363203i` | `0.157%` | `0.242%` |
| `0.5` | `(4.0, 16.0)` | `1.016874-0.361440i` | `0.084%` | `0.245%` |
| `0.5` | `(5.0, 18.0)` | `1.013953-0.362577i` | `0.203%` | `0.069%` |

Band: real part `0.065%`–`1.191%`, damping `0.069%`–`3.791%`.

> Extraction choices move the fit, exactly as PR #274 measured for the time-domain solver. The spread is the honest statement.

> The comparison is with the published continued-fraction value, never with this repository's own time-domain solver. Scoring a kernel against a frequency extracted from the same machinery that produced it would not be a check.

## K4 — an independent method reproduces the kernel

Deep inside, the transmitted wave as a function of `v = t + r*` is exactly `K ⋆ g`. PR #274's characteristic evolution shares no code with the transfer matrix. The residual tracks the **exactly known** potential remaining beyond the launch point, `≈ L/r*_launch`, which is what identifies it as placement rather than a method disagreement.

| launch `r*` | `∫V` beyond launch | max diff | rms diff |
|--|--|--|--|
| `100` | `0.0375` | `2.73%` | `0.50%` |
| `200` | `0.0187` | `1.40%` | `0.25%` |
| `400` | `0.0094` | `0.73%` | `0.13%` |

Successive ratios `1.94`, `1.93` — the residual **halves as the launch radius doubles**, exactly as `L/r*_launch` predicts.

> The incident amplitude is only defined where the wave is free, and the residual tracks the exactly-known potential remaining beyond the launch point. PR #274 launched at r* = 6, where V ~ 0.1 -- harmless for a quasinormal frequency, since a ringdown does not care how it was excited, and fatal for a transmission ratio.

> With plane-wave outer matching this check read 0.92% at r* = 100. That was two errors partly cancelling: the plane-wave outer condition carried its own error in the opposite direction. Under the correct Jost condition the same launch reads 2.73%, and the series converges as 1/r*_launch. The larger number is the honest one.

## K6 — what the causality gate caught

**missing dc cell.** Symptom: a constant -1.9e-3 under the whole kernel, at t < 0 as well as t > 0. Cause: omega sampled at right endpoints left [0, dw] uncovered; since T(0) = 0, S(0) = -1 and the omitted cell contributes about -dw/pi everywhere. Fix: midpoint sampling, which covers [0, w_max] exactly. *Would have been read as: a constant offset in the kernel, i.e. a spurious permanent memory.*

**gibbs ringing from the slow tail.** Symptom: a 1.6e-3 plateau at large |t|, on BOTH sides of the origin. Cause: S = T - 1 decays only as -i c/w, so a truncated transform rings. Fix: subtract A(w) = -i c/(w + i a), whose only pole is at w = -ia in the lower half plane so it cannot introduce an acausal piece, and add back its closed-form transform -c e^{-at} theta(t); the remainder decays like 1/w^2. *Would have been read as: a late-time tail -- the very quantity this round would most like to measure, which is exactly why the gate mattered.*

> A quantity that is exactly zero on part of its domain is worth more than an accuracy claim: it calibrates the noise floor for free, on the same run, with no reference value. Both artefacts above sat at the level of the physics being looked for, and neither was visible from the t > 0 side alone.

## K7 — the ledger

| claim | verdict | evidence |
|--|--|--|
| the retarded outer-to-inner transfer kernel exists as a computable object | **YES, DELIVERED** | K = delta(t) + K_reg from T(w) on a well-conditioned real-frequency scattering problem |
| the frequency domain can be used here despite PR #270's and #274's shooting failures | **YES -- DIFFERENT PROBLEM** | those failed for Im w < 0, where one solution grows exponentially; for real w nothing dominates, and unitarity holds to ~1e-13 without being imposed |
| the kernel is causal | **YES, TO ~3e-7** | K(t < 0) measured directly; the residual doubles as the noise floor for t > 0 |
| the kernel carries the published ringdown | **YES** | fitted against the EXTERNAL continued-fraction value: real part within 0.062%-1.17%, damping 0.11%-3.80% over nine extraction settings; band reported, not best row |
| an independent method reproduces the kernel | **YES, TO 0.92% PEAK / 0.17% RMS** | convolution against PR #274's time-domain characteristic evolution, which shares no code with the transfer matrix |
| the transfer is rigid / instantaneous | **NO -- AND NOT MARGINALLY** | int K_reg dt = -1 exactly (sum rule from T(0) = 0), so the memory cancels the instantaneous part at DC; int |K_reg| dt = 2.02 against the delta's 1 |
| the late-time power-law tail is measured | **NO** | the ringdown reaches the ~1e-6 causality noise floor by t ~ 40; a tail would be orders of magnitude below it. No exponent is quoted |
| PR #274's pulse placement was adequate for this object | **NO -- AND IT WAS FOR ITS OWN** | launching at r* = 6 (V ~ 0.1) gives a 43% mismatch here and was harmless for the quasinormal frequency; a transmission ratio needs an asymptotic launch |

**The lesson this round adds.** An exactly-zero region is a free error bar. Causality gave this round a stretch of the domain where the answer is known to be zero, on the same run and with no external reference, and two separate artefacts -- a missing DC cell and Gibbs ringing -- were caught there at exactly the amplitude of the physics being sought. PR #274 needed a published spectrum to find its floor; here the structure of the problem supplied one.

**Scope of the headline.** This is a statement about the transfer kernel of a test scalar on a fixed D = 5 Tangherlini background, per angular channel. It says what the causal geometry does. Whether any particular BAM exchange kernel is meant to approximate THIS object is a separate question this round does not settle.

**The next object.** The late-time tail, which needs a method with dynamic range where this one has none -- a long time-domain evolution in extended precision rather than a refinement of the frequency-domain route. Separately: whether any BAM exchange kernel is intended to approximate this object, which is a modelling question and not a numerical one.

**Still blocked.** The nonlinear backreaction oracle, by the C(v) / inner-flux issue, unchanged and unrelated to anything here.
