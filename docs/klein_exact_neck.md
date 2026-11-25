# The exact-neck law for Klein splitting: held-out result

This tests `docs/klein_exact_neck_prereg.md`, frozen at `834ac70`.

* **Archive.** `experiments/closure_ledger/runs/20261124_klein_exact_neck/`
  (`probe.json`, `probe.md`).
* **Code.** `experiments/closure_ledger/klein_exact_neck_probe.py`, committed
  at `39adb40` before the run, and
  `experiments/closure_ledger/klein_exact_neck_predictions.py`, which is part
  of the freeze.
* **Solver.** PR #355's `geometrodynamics/bulk/alice_handle_ports.py`,
  unchanged from `07d58de`.

## Verdict

**`CONFIRMED`.** The reproduction check and all three gates pass at both
resolutions.

* **Reproduction.** All 17 in-sample splittings of #355 and #362 are
  reproduced to 5.1e-8 relative, against a gate of 1e-6.
* **G1, accuracy.** At every gated cell, `|r − 1| ≤ 0.2 a⁴ + 1e-3`, with
  `r = Δ_measured/E`. The worst cell uses 0.38 of its bound.
* **G2, the neck is the error.** At every gated cell,
  `|r − 1| ≤ 0.2 |r_H − 1|`. The worst ratio is 0.068: the exact neck removes
  at least 93% of `H`'s discrepancy.
* **G3, resolution.** The worst `|Δ₂₀₀₀/Δ₄₀₀₀ − 1|` is 3.0e-4, against 1e-3.

| a | m | predicted Δ = E | measured Δ, N=4000 | r − 1, N=2000 | r − 1, N=4000 | bound | r_H − 1, N=4000 | \|Δ₂₀₀₀/Δ₄₀₀₀ − 1\| |
|---|---|---|---|---|---|---|---|---|
| 0.25 | 1 | 2.1249e-3 | 2.1258e-3 | +0.00041 | +0.00043 | 1.78e-3 | 0.034 | 2.8e-5 |
| 0.25 | 2 | 3.0666e-6 | 3.0671e-6 | −0.00014 | +0.00017 | 1.78e-3 | 0.016 | 3.0e-4 |
| 0.35 | 1 | 8.3890e-3 | 8.4019e-3 | +0.00153 | +0.00154 | 4.00e-3 | 0.073 | 7.8e-6 |
| 0.35 | 2 | 4.7691e-5 | 4.7742e-5 | +0.00098 | +0.00107 | 4.00e-3 | 0.035 | 8.2e-5 |
| 0.35 | 3 | 2.3234e-7 | 2.3246e-7 | +0.00020 | +0.00049 | 4.00e-3 | 0.023 | 2.9e-4 |
| 0.45 | 1 | 2.3761e-2 | 2.3834e-2 | +0.00310 | +0.00310 | 9.20e-3 | 0.136 | 3.1e-6 |
| 0.45 | 2 | 3.8202e-4 | 3.8323e-4 | +0.00313 | +0.00316 | 9.20e-3 | 0.063 | 3.2e-5 |
| 0.45 | 3 | 5.2604e-6 | 5.2694e-6 | +0.00161 | +0.00172 | 9.20e-3 | 0.041 | 1.1e-4 |
| 0.45 | 4 | 6.7961e-8 | 6.8029e-8 | +0.00074 | +0.00101 | 9.20e-3 | 0.030 | 2.6e-4 |
| 0.55 | 1 | 5.5167e-2 | 5.5336e-2 | +0.00305 | +0.00305 | 1.93e-2 | 0.244 | 1.5e-6 |
| 0.55 | 2 | 2.0630e-3 | 2.0779e-3 | +0.00724 | +0.00725 | 1.93e-2 | 0.106 | 1.5e-5 |
| 0.55 | 3 | 6.6050e-5 | 6.6334e-5 | +0.00426 | +0.00431 | 1.93e-2 | 0.067 | 5.1e-5 |
| 0.55 | 4 | 1.9820e-6 | 1.9874e-6 | +0.00261 | +0.00273 | 1.93e-2 | 0.049 | 1.2e-4 |
| 0.55 | 5 | 5.7432e-8 | 5.7537e-8 | +0.00159 | +0.00183 | 1.93e-2 | 0.039 | 2.4e-4 |

Over these cells the splitting runs from 5.5e-2 down to 5.8e-8, six
decades. The law predicts it to at most 0.73%, and to at most 0.16% for
`a ≤ 0.35`. Nothing is fitted. On the same cells `H` is off by up to 24%.

## What it means

The pre-registration put the remaining error of `H` in the WKB step through
the neck. On fresh cells, replacing that step with the exact neck solution
removes 93–99% of the discrepancy. The splitting of the Klein doublets on
#355's handle is then fixed by three things:
* the exact band state at the sphere's frequency;
* the exact neck solution that vanishes at the far rim;
* Herring's formula, with the twist as a Bloch phase.

The reversal time `π/(2E)` follows from the geometry to the same accuracy.
Per cell it is in `probe.json` as `reversal_time`. It runs from 28 at
(0.55, 1) to 2.7e7 at (0.55, 5).

The empirical envelope held out of sample. The pre-registration recorded
`|r − 1| ≤ 0.11 a⁴` in-sample; held out, `max |r − 1|/a⁴` at N = 4000 is
0.111. At N = 4000 `r − 1` is positive at every gated cell.

## Diagnostics, not gated

**a = 0.65 and 0.7, beyond the in-sample range.** `|r − 1|` stays within
the formula `0.2 a⁴ + 1e-3`, using at most 0.37 of it. The worst is 1.68% at
(0.7, 2). `E` still removes at least 90% of `H`'s error, which there reaches
73%. At `m = 1` the sign of `r − 1` reverses: −0.29% at `a = 0.65` and
−0.93% at `a = 0.7`.

| a | m = 1 | m = 2 | m = 3 | m = 4 | m = 5 |
|---|---|---|---|---|---|
| 0.65, r − 1 | −0.0029 | +0.0135 | +0.0092 | +0.0062 | +0.0043 |
| 0.7, r − 1 | −0.0093 | +0.0168 | +0.0126 | +0.0088 | +0.0063 |

**Excited pairs, n = 1, 2.** `r` is the measured splitting divided by `E`,
using the band state `P_{m+n}^m`.
* For the 34 of 48 pairs with predicted `Δ/ω < 1e-2`, `r` is 0.983–1.019.
* For the other 14 it is 0.65–1.016, falling as the splitting grows. The
  values are 0.914 at (0.55, 1, n = 2), 0.740 at (0.65, 1, n = 2) and 0.651
  at (0.7, 1, n = 2), where the predicted `Δ/ω` is 0.11, 0.23 and 0.32.
  There Herring's leading order does not apply.
* #362's WKB form reached 2.8 on excited pairs and could not be evaluated
  past `ω = m/sin a`. The exact neck has neither problem.

**The rim kink.** The pre-registration noted, without using it, that
in-sample `r − 1` is close in size to the slope mismatch at the rim,
`|D′(−L)/(−m cos a) − 1|`. On the held-out cells with `m ≥ 2` and `r − 1`
well above the resolution, the ratio of the two is 0.92–1.04. At `m = 1` it
is looser, from 1.2 to 2.1. At `m = 1` and `a ≥ 0.65`, `r − 1` changes sign
but the kink does not. This was looked at after the fact. It points to the
next correction, a `C¹` single-well state, but it is not a result.

## Checks

* **Sector purity.** The archive's `purity` is the worst cosine/sine mixing
  over the whole spectrum, as in #362; it reaches 0.49. Checked directly at
  (0.25, 1), (0.45, 1) and (0.7, 1) with N = 2000:
  * the six lowest modes have cosine fractions of exactly 0 or 1 to 14
    digits;
  * the most mixed mode is at rank 3972–3990 of 4000, with `ω` near 1300.
* **Resolution.** All 24 signed `Δ₂₀₀₀/Δ₄₀₀₀ − 1` are negative, from −6.5e-7
  to −3.0e-4.

## Procedure notes

* **One coarse look before the run.** While smoke-testing the probe, after
  the freeze, one gated cell, (0.45, 2), was solved at N = 400 and printed
  `r − 1 = 0.0021`. Nothing in the probe was changed after that. The probe
  was committed at `39adb40` before the run.
* **The run was interrupted.** A container restart stopped it after 63 of
  its 82 solves. The 63 had been cached as they finished. The rerun took
  them from the cache and solved the remaining 19, so the archived
  `seconds`, 620, covers only the resumed part. The reproduction check
  passed on the cached solves, and every measurement in the archive comes
  from the unchanged probe.

## Scope

As frozen. This is a classical scalar on a fixed 2D surface, not a GR
solution. The tunnelling is classical wave tunnelling; nothing here is
claimed to be quantum. The residual of up to about `0.11 a⁴` is measured,
not derived. #355's second gate, the dynamical reflection axis, stays
closed.
