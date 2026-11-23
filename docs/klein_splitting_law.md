# The tunnelling law for Klein splitting: held-out result

Test of `docs/klein_splitting_law_prereg.md`, frozen at `56055ab`. Archive:
`experiments/closure_ledger/runs/20261123_klein_splitting_law/` (`probe.json`,
`probe.md`). Code: `experiments/closure_ledger/klein_splitting_law_probe.py`,
committed at `de35f55` before its results were read, and
`experiments/closure_ledger/klein_splitting_law_predictions.py`, which is part
of the freeze. Solver: PR #355's `geometrodynamics/bulk/alice_handle_ports.py`,
unchanged from `07d58de`.

## Verdict

**`CONFIRMED`.** The reproduction check and all three gates pass at both
resolutions.

* **Reproduction.** The probe reproduces #355's six in-sample splittings to
  5.1e-8 relative, against a gate of 1e-6.
* **G1, the law.** At every gated cell `0 < m(r − 1) ≤ 0.25`, where
  `r = Δ_measured / H`.
* **G2, the error ordering.**
  * c(0.2) = 0.021 < c(0.4) = 0.095.
  * Every a = 0.4 value of m(r − 1), from 0.091 to 0.100, lies in
    [0.048, 0.182].
* **G3, resolution.** The worst `|Δ₂₀₀₀/Δ₄₀₀₀ − 1|` is 7.2e-4, against 1e-3.

| a | m | predicted Δ = H | measured Δ, N=4000 | r | m(r − 1), N=2000 | m(r − 1), N=4000 | Δ₂₀₀₀/Δ₄₀₀₀ − 1 |
|---|---|---|---|---|---|---|---|
| 0.2 | 1 | 8.4421e-4 | 8.6200e-4 | 1.0211 | 0.021 | 0.021 | 6.6e-5 |
| 0.2 | 2 | 4.9980e-7 | 5.0487e-7 | 1.0101 | 0.019 | 0.020 | 7.2e-4 |
| 0.4 | 1 | 1.3262e-2 | 1.4592e-2 | 1.1003 | 0.100 | 0.100 | 4.8e-6 |
| 0.4 | 2 | 1.3728e-4 | 1.4375e-4 | 1.0471 | 0.094 | 0.094 | 4.9e-5 |
| 0.4 | 3 | 1.1766e-6 | 1.2127e-6 | 1.0306 | 0.091 | 0.092 | 1.7e-4 |
| 0.5 | 4 | 3.8280e-7 | 3.9757e-7 | 1.0386 | 0.154 | 0.154 | 1.8e-4 |

Across these cells the splitting spans more than four decades, from
1.5e-2 to 4e-7. The closed form, with nothing fitted, predicts it to
between 1% and 10%.

## What it means

The splitting of the Klein doublets on #355's handle is classical wave
tunnelling through one barrier on a ring. The barrier is the neck together
with the band's forbidden strips at the rims, and the torus twist is a Bloch
phase.

The law reproduces three things:
* the size of the splitting;
* its dependence on `m`, including the extra decay `ln(1/sin²a)` per unit `m`
  that #355's expectation omitted;
* its dependence on the aperture, through `c(a)`.

The handle's reversal time `π/(2Δ)` follows from the geometry alone. Per
cell it is in `probe.json` as `reversal_time`; for example 1.3e6 at
`a = 0.4`, `m = 3`.

This explains #355's P4 counterexample but does not change it. That verdict
stays `COUNTEREXAMPLE` under its own freeze.

## Diagnostics, not gated

**a = 0.6.** As frozen, this angle was outside the in-sample range of `c(a)`.

| m | r | m(r − 1) |
|---|---|---|
| 1 | 1.333 | 0.333 |
| 2 | 1.136 | 0.272 |
| 3 | 1.086 | 0.257 |
| 4 | 1.062 | 0.249 |
| 5 | 1.049 | 0.244 |

The law still holds to `O(1/m)`. But `m(r − 1)` exceeds the gated bound of
0.25 for `m ≤ 3`, so these cells would have failed G1. They show where the
frozen tolerance stops.

Across all angles tested, `c(a)` is 0.020, 0.049, 0.095, about 0.16 and about
0.25 at `a` = 0.2, 0.3, 0.4, 0.5 and 0.6. That is roughly `0.5–0.7 a²`. This
is an observation, not a frozen claim.

**Excited pairs.** For `n = 1, 2`, `r` uses the same construction with the
exact band state `P_{m+n}^m(cos θ)`.
* Deep below the rim barrier it is 1.02–1.20.
* It degrades as `ω` nears `m/sin a`, reaching 2.4 at `a = 0.4, m = 1, n = 1`
  and 2.8 at `a = 0.6, m = 2, n = 1`. There the rim is close to a turning
  point and the WKB step from the rim fails.
* Cells where `ω ≥ m/sin a` are not evaluated.

This is why the freeze gated `n = 0` only.

## Checks

* **Sector purity.** The archive's `purity` is the worst cosine/sine mixing
  over the *whole* spectrum. In some cells it reaches 0.50. Checked directly
  at (0.5, 4), (0.2, 1) and (0.6, 1) with N = 2000:
  * the four lowest modes have cosine fractions of 1 and below 1e-12;
  * the most mixed mode is near the top of the spectrum (index 3983 of
    4000, ω ≈ 1341), where neighbouring grid-scale modes are nearly
    degenerate.

  The lowest pair is unaffected.
* **Resolution.** The measured `Δ` agrees between N = 2000 and 4000 to
  between 1e-6 and 7e-4. Over the same range `r − 1` is 0.01–0.10.

## Relation to PR #360

PR #360 registered (`73cca3f`, after this freeze) a leading small-aperture
law, `C(m) tan^{4m}(a/2)`, and tested it at `a = 0.08–0.24` with a continuum
Wronskian solver. That law is the `a → 0` limit of `H`: there `κ_rim → m`,
`sin^{2m} a → 4^m e^{−2mL}` and `e^{−W_neck} → e^{−2mL}`, with
`L = ln cot(a/2)`. The two tests are complementary:
* #360 measures the leading law at small apertures, in the continuum;
* this test measures the finite-aperture form at `a = 0.2–0.5`, with #355's
  finite-element solver.

Two comparisons, made after this test was scored:

* **The solvers agree.** On the three cells both PRs computed, (0.4, 2),
  (0.6, 2) and (0.6, 4), this test's finite-element relative split at
  N = 4000 matches #360's continuum Wronskian to 1.6e-5, 3.5e-6 and 2.9e-5.
  This independently checks both solvers.
* **`H` is not uniformly better than the leading law.** Measured/predicted
  at N = 4000:

  | a | m | measured/H | measured/leading law |
  |---|---|---|---|
  | 0.2 | 1 | 1.021 | 1.002 |
  | 0.2 | 2 | 1.010 | 1.003 |
  | 0.4 | 1 | 1.100 | 1.019 |
  | 0.4 | 2 | 1.047 | 1.029 |
  | 0.4 | 3 | 1.031 | 1.039 |
  | 0.5 | 4 | 1.039 | 1.080 |
  | 0.6 | 1 | 1.333 | 1.029 |
  | 0.6 | 2 | 1.136 | 1.061 |
  | 0.6 | 3 | 1.086 | 1.081 |
  | 0.6 | 4 | 1.062 | 1.101 |
  | 0.6 | 5 | 1.049 | 1.121 |

  The errors have different structure. The leading law's error grows with
  `m` at fixed `a`. The error of `H` falls as `1/m`, the size of the WKB
  step through the neck at finite `ω`. So `H` is more accurate only for
  `m ≳ 3`, and the "finite-aperture" corrections it adds are not
  improvements at small `m`. A formula with the neck solved exactly at the
  mode's frequency, as #360's solver does numerically, should remove both
  errors. That is not tested here.

## Scope

As frozen. A classical scalar on a fixed 2D surface, not a GR solution. The
tunnelling is classical wave tunnelling; nothing here is claimed to be
quantum. #355's second gate, the dynamical reflection axis, stays closed.
