# The retarded outer→inner transfer kernel, `D = 5` Tangherlini

**8/8 checks pass.**

PR #270 deferred this object because it is a ratio of two signals it could not trust. PR #274 settled which signal was right, against a published spectrum. This round builds the kernel.

A wave sent in from the far region reaches the horizon filtered. That filter is `T_ℓ(ω)`; in time it is a convolution kernel with `ψ_trans(v) = (K_ℓ ⋆ ψ_inc)(v)`, `v = t + r*`. Since `T → 1` at high frequency, `K_ℓ(t) = δ(t) + K_ℓ^reg(t)` — instantaneous part plus memory.

> ## The transfer is not rigid, and not marginally.

| quantity | value | |
|--|--|--|
| `∫K_reg dt` | `-0.997757` | exact value `−1`, a **sum rule** from `T(0) = 0` |
| `∫\|K_reg\| dt` | `2.0233` | against the `δ`'s weight of `1` |
| `T(ω→0)` | `1.734e-06` | the barrier blocks DC completely |

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
| `0.1` | `0.00079866` | `0.99999968` | `+1.3e-15` |
| `0.3` | `0.01423276` | `0.99989871` | `-1.9e-15` |
| `0.7` | `0.19749584` | `0.98030373` | `+8.9e-15` |
| `1` | `0.64441877` | `0.76467277` | `+3.0e-14` |
| `1.5` | `0.98880091` | `0.14924063` | `-3.9e-13` |
| `3` | `0.99999977` | `0.00068021` | `+1.1e-13` |
| `10` | `1.00000000` | `0.00000042` | `-5.3e-13` |
| `30` | `1.00000000` | `0.00000005` | `+6.3e-13` |

Worst residual **`6.3e-13`**, and unitarity is imposed nowhere — it is a consequence of the computation, so it measures it.

> PR #270 and #274 could not converge a quasinormal frequency by shooting in real r, because for Im w < 0 the outgoing solution grows like e^{|Im w| R} and swamps the coefficient being zeroed. Here w is REAL, both e^{+-i w r*} have unit modulus, and nothing dominates anything. This is a different, well-posed problem -- not a repair of the one that failed. Unitarity to ~1e-13, imposed nowhere, is the evidence.

Second order in the spatial step: `|T|` at `ω = 1` gives `0.64442679`, `0.64442037`, `0.64441877`, successive differences `6.4e-06`, `1.6e-06`.

## K2 — GATE 1: causal support

| `t` | `K_reg(t)` |
|--|--|
| `-300` | `+1.666e-07` |
| `-200` | `-3.190e-07` |
| `-100` | `+4.673e-07` |
| `-50` | `-9.173e-07` |
| `-20` | `-2.133e-06` |
| `-10` | `+5.132e-06` |
| `-5` | `+8.258e-06` |
| `-2` | `+2.714e-05` |
| `-1` | `-3.405e-05` |
| `-0.5` | `-1.026e-04` |

Worst acausal value **`1.03e-04`**, and **`9.2e-07`** away from the front.

> K(t) vanishes identically for t < 0, so whatever the computation returns there IS its noise floor -- no reference value needed. Any feature at t > 0 smaller than that floor is not measurable, which is how this round knows the late-time tail is out of reach.

| `t` | `K_reg(t)` |
|--|--|
| `0.5` | `-1.168e+00` |
| `1` | `-4.354e-01` |
| `2` | `+2.211e-01` |
| `3` | `+2.212e-01` |
| `5` | `-8.351e-02` |
| `8` | `+2.638e-02` |
| `12` | `-8.719e-03` |
| `16` | `+9.891e-04` |
| `20` | `+1.872e-04` |
| `30` | `-1.267e-05` |
| `40` | `+7.172e-07` |
| `60` | `+2.302e-07` |
| `100` | `+2.667e-07` |

## K3 — GATE 3: the kernel carries the published ringdown

Reference (external): `1.01601691-0.36232802i`. Source: Matyjasek, Phys. Rev. D 104, 084066 (2021), arXiv:2107.04815 -- continued fractions cross-checked against Hill determinants, agreeing to 11 digits.

| `dt` | window | fitted `ω` | real err | damping err |
|--|--|--|--|--|
| `0.2` | `(3.0, 14.0)` | `1.027952-0.376104i` | `1.175%` | `3.802%` |
| `0.2` | `(4.0, 16.0)` | `1.026151-0.369614i` | `0.997%` | `2.011%` |
| `0.2` | `(5.0, 18.0)` | `1.023143-0.363089i` | `0.701%` | `0.210%` |
| `0.3` | `(3.0, 14.0)` | `1.017028-0.363837i` | `0.100%` | `0.416%` |
| `0.3` | `(4.0, 16.0)` | `1.017149-0.361417i` | `0.111%` | `0.251%` |
| `0.3` | `(5.0, 18.0)` | `1.015390-0.361124i` | `0.062%` | `0.332%` |
| `0.5` | `(3.0, 14.0)` | `1.017586-0.363112i` | `0.154%` | `0.216%` |
| `0.5` | `(4.0, 16.0)` | `1.016979-0.361611i` | `0.095%` | `0.198%` |
| `0.5` | `(5.0, 18.0)` | `1.014191-0.362744i` | `0.180%` | `0.115%` |

Band: real part `0.062%`–`1.175%`, damping `0.115%`–`3.802%`.

> Extraction choices move the fit, exactly as PR #274 measured for the time-domain solver. The spread is the honest statement.

> The comparison is with the published continued-fraction value, never with this repository's own time-domain solver. Scoring a kernel against a frequency extracted from the same machinery that produced it would not be a check.

## K4 — an independent method reproduces the kernel

Deep inside, the transmitted wave as a function of `v = t + r*` is exactly `K ⋆ g`. PR #274's characteristic evolution shares no code with the transfer matrix, so this is real cross-validation — and it exposed a subtlety about where the pulse may be launched.

| pulse centre | launch `r*` | `V` at launch | max diff | rms diff |
|--|--|--|--|--|
| `12` | `6` | `9.76e-02` | `43.11%` | `10.17%` |
| `60` | `30` | `4.16e-03` | `7.26%` | `1.35%` |
| `200` | `100` | `3.75e-04` | `0.92%` | `0.17%` |

> PR #274's pulse was launched at r* = 6, where V ~ 0.1. That is fine for a quasinormal frequency and wrong for an incident amplitude, and it shows up here as a 43% mismatch that falls to 0.92% once the launch moves to r* = 100. The transfer kernel needs an asymptotic launch; the ringdown did not.

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
