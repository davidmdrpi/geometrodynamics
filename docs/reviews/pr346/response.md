# PR #346 review response: preserve the failed candidate, correct its interpretation

Review: [5581410007](https://github.com/davidmdrpi/geometrodynamics/pull/346#pullrequestreview-5581410007),
at `2073e49`. The negative result remains: no sustained approach, no 25% radius
threshold crossing, and no collapse in the archived nine cases. This revision
changes interpretation and adds reproducible diagnostics. It does not replace
the frozen experiment with a different model or silently tune its parameters.

## Changes made

1. **Angular balance:** remove the claim that the logged correction demonstrates
   a physical clock support exchanging angular momentum. It is a coordinate
   canonical-balance correction, nonzero even when the ports are off and the
   mechanical momenta remain zero. The old solver and field names remain for
   archive reproducibility, not as physical terminology.
2. **Resolution:** label the throat fractions and mechanical response as
   unconverged. The committed 8π maxima are 0.1744% and 2.4496% for sphere
   cutoffs 10 and 14. This is a sensitivity comparison, not a convergence claim.
   The review's additional 4π transmission/resolution scans are attributed to
   the review; they are not rerun or mixed with these 8π figures.
3. **Radius:** replace "2.21% shrinkage" with "instantaneous dip" and report
   full radial-period averages. Baseline and dephased minimum means are
   0.995174 and 0.980937. Every case's return-window mean and half-range are
   saved separately from the frozen results.
4. **Publication:** use the published snapshot instead of a broken GIF embed.
   The original animation remains missing; its manifest entry and recovery
   identifiers are preserved in the PR.

## Reproducible numerical checks

`experiments/closure_ledger/coupled_mouth_throat_review.py` postprocesses the
existing CSVs and records their hashes in [numerical_review.json](numerical_review.json).
It does not overwrite the production manifest, solver, tests or trajectories.

With `--ports-off`, it also evolves six diagnostic free-field controls to 2π:
offsets −0.3 and −1.2 at l_max = 10, 14, 18. It compares the canonical angular
quantity, its integrated correction, and the free sphere's physical rotation
charge. The latter uses the review's current-flux expression with an overall
minus sign to match the solver's convention J = −pᵀDu:

    g_t = ∇u − ∇Θ u_τ
    L = (u_τ² − |g_t|² − u²/4)/2
    Q = −∫[(u_τ + ∇Θ·g_t)(ξ·g_t) + (∇Θ·ξ)L] dA.

This diagnostic is for the **free sphere**. It does not supply a total covariant
charge for the finite-port mechanical action. In particular, showing that a
coordinate correction closes an identity cannot remove the position-dependent
neck forces already present in that action. Those remain a model defect for
the proposed MTY interpretation.

The six controls reproduce the review's physical-charge ranges (max − min):

| Imposed offset | l_max = 10 | l_max = 14 | l_max = 18 |
| --- | ---: | ---: | ---: |
| −0.3 | 5.817 × 10⁻⁴ | 8.652 × 10⁻⁶ | 3.589 × 10⁻⁸ |
| −1.2 | 0.03755 | 0.01557 | 0.005037 |

For offset −1.2 at l_max=10, the canonical charge range is 0.03296675 and
the coordinate correction range is 0.03296675, despite zero mouth momenta.
The corrected canonical sum has range 3.40 × 10⁻⁹. Its accurate cancellation
does not cure the physical charge's truncation error. All **14 tests pass**:
the original ten and four review-diagnostic checks.

```bash
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.coupled_mouth_throat_review --ports-off
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_coupled_mouth_throat.py tests/test_coupled_mouth_throat_review.py
```

The radial average is explicitly post-hoc: a full window of length
T_r = 2π/√(k/μ), using piecewise-linear interpolation of the archived samples.
Window starts are the eligible sample times plus the last possible full window.
The averaging tests distinguish a pure oscillation from a moving mean, test
nonuniform sampling/endpoints, and reject incomplete windows. Rotation tests
check the untilted sign and a nonzero correction with no mechanical support.

## Structural correction still required: clocks from a birth history

The imposed −1.2 offset has no birth-history derivation. A future experiment
must begin at contact with synchronized clocks and evolve independent mouth
histories. The newer #345 implementation at
[`c0a1ae5`](https://github.com/davidmdrpi/geometrodynamics/commit/c0a1ae55e509efb89d7543cabefe78345d235bde),
`docs/pair_birth_timeline.md` §7, is a useful starting point. It reports zero
offset for symmetric births and transient asymmetric offsets up to 0.215 in
the examples. Its own limitations state that these clocks are bookkeeping:
there is no neck-channel feedback. Its release/re-creation steps are stated
energy/momentum-matching rules, not a derived topology change.

Three distinctions matter before implementing the review's proposal:

**Proper-time difference versus external-time offset.** Let the two clocks be
s_A(t), s_B(t), with ds_i/dt = √(1−v_i²) in the proposed ultrastatic background.
The simultaneous clock difference C(t) = s_B(t)−s_A(t) obeys

    dC/dt = √(1−v_B²) − √(1−v_A²).

But the throat's equal-clock external-time map is

    F(t) = s_B⁻¹(s_A(t)),    δ(t) = F(t)−t,
    dδ/dt = s_A'(t)/s_B'(F(t)) − 1.

Thus dC/dt is not generally dδ/dt. To first order δ ≈ −γ_B C, as #345's
history document states with the opposite clock-difference convention. An
implementation must retain both clocks and invert a monotone history, or
state and bound an approximation. It must not extrapolate silently when the
matched event lies outside the computed history. Symmetric births yield δ=0;
asymmetric birth momenta replace arbitrary imposed-offset controls.

**Action and work signs.** For a neck parameter A(t), with the other coordinates
on shell, dE/dt = −(∂L/∂A) A_dot. A row accumulating +(∂L/∂A) A_dot is the
compensating clock reservoir, not energy delivered to the field. If A is
constrained by the mouth histories, that constraint must also be varied;
postprocessing clocks on an unchanged trajectory is not a closed action or a
mechanical backreaction calculation. The free sphere may use an untilted chart,
but the current finite-port action is not established to be gauge covariant.
Simply changing its tilt does not establish physical equivalence.

**Chronology versus chart failure.** |A| = 0.9 L is a chart-safety threshold,
not a chronology horizon. For stationary mouths and a fixed traversal map,
d + L_neck − |δ| is the familiar null-loop margin. For moving mouths and a
time-dependent δ, it is a diagnostic estimate: an actual closed-null-loop
claim needs an event-to-event transit and return map evaluated on the moving
worldlines. Neither a chart stop nor an instantaneous inequality alone proves
that such a loop occurred.

## Requirements for the next mechanism test

- Derive the history-coupled neck action and a complete field/mouth/clock ledger.
  First validate equal-speed and constant-unequal-speed clock histories against
  exact equal-clock matching, and verify the clock-work sign.
- Calibrate the port against `HandleWaveSim` using matched initial fields and
  the **same observable** on both sides. The 7.5% sideways rim-momentum fraction
  is not an energy-transmission fraction and is not a direct calibration target
  for the 0.17–2.5% energy numbers.
- Establish sphere and neck resolution separately. If a tilted chart is used,
  pull initial data back from the same physical packet rather than repeating
  its coordinate profile on different slices.
- Include a predeclared radial detuning control (for example k=11.25 gives
  Ω_r=7.5) and require declining averaged radius, not merely an oscillation
  trough. Detuning from one integer sum-frequency does not prove the absence
  of all resonant responses.
- Treat contact release, antipodal refocusing and re-creation as explicit
  transfer rules until a topology-transfer law is derived. Resetting clocks at
  re-creation is part of that rule and must be disclosed.

These structural requirements remain open in #346. The corrected archive is
evidence against the particular imposed-clock scalar candidate, not a test
of the proposed complete-history MTY feedback mechanism.
