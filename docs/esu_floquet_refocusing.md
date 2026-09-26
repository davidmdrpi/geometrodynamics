# Floquet spectrum and antipodal refocusing of the four-scalar ESU

**Certified result (tensor sector).** Gravitational-wave perturbations of the
breathing four-scalar ESU are elliptic for every degree n=2..80. They
refocus at the antipode after conformal time pi, and the refocusing becomes
exact as n grows: `1-F_n ~ n^-1.952`. The support's own mass term
reduces the phase error of the bare ESU by a factor of about 5.5. The frozen
WKB prediction holds to 3.5e-5.

**Scalar sector, certified by the prospective G3 extension.** The frozen G3
left V and S UNRESOLVED; that record stands. A dated addendum
([`9a8ef99`](esu_floquet_g3_extension_prereg.md)) was published before it
was run, and under it S passes. The scalar sector has a **parametrically
unstable n=2 mode**:
- multiplier 1.2475 per transit, with period doubling against the breathing
  background;
- it lies in the RP3-compatible sector;
- degrees 3..80 are elliptic;
- the scalar refocusing fidelity saturates near .994 instead of approaching
  1.

**Vector sector: still UNRESOLVED.** It fails the extension's trace-ratio
test as well. A diagnostic, not a verdict, shows that its matrices converge
at exactly fourth order. The trace is first-order insensitive to phase error
for maps near ±I, so its errors cancel (section 3). The uncertified vector
maps are elliptic, with asymptotic refocusing matching WKB to 0.23%.

The freeze [`4e65c3e`](esu_floquet_refocusing_prereg.md) was published in
[#310](https://github.com/davidmdrpi/geometrodynamics/pull/310) before
implementation or measurement. No frozen equation, degree range, tolerance
or threshold was changed. Section 7 lists every operational choice made
after the freeze; each was fixed before any result was viewed.

## 1. Verdicts (as frozen)

| Verdict | Result |
|---|---|
| T_STABILITY | **ELLIPTIC_2_TO_80** |
| T_REFOCUSING | **ASYMPTOTIC** (p = 1.952 ± 0.001; F_n > .999 from n=6) |
| T_WKB_PREDICTION | **PASS** (measured (n+1)θ_n = 0.283349, predicted 0.283339) |
| V_STABILITY, V_REFOCUSING, V_WKB_PREDICTION | UNRESOLVED (G3) |
| S_STABILITY, S_REFOCUSING | UNRESOLVED (G3) |
| *Extension (9a8ef99):* S_STABILITY_EXT | **HYPERBOLIC_AT[2]** (also in the odd subset) |
| *Extension:* S_REFOCUSING_EXT | ASYMPTOTIC under the frozen classifier; the data show saturation (section 6) |
| *Extension:* T_STABILITY_EXT, T_REFOCUSING_EXT | ELLIPTIC_2_TO_80, ASYMPTOTIC (confirms the frozen verdicts) |
| *Extension:* V | UNRESOLVED (trace-ratio test fails; section 3) |

For the odd (RP3-compatible) subsets, T (n even) is also ELLIPTIC and
ASYMPTOTIC, with p = 1.952 ± 0.002. V and S are UNRESOLVED.

## 2. Gates and controls that passed

- **G1, explicit harmonics.** The reduced equations were imposed on the full
  coordinate linearization of the Einstein–sigma equations, including every
  Einstein component and all four field equations. Maximum relative
  residuals:
  - tensor n=2 (left-invariant σ1²−σ2²): 3.9e-16;
  - vector n=2,3 (toroidal C⁽²⁾ₙ₋₁): 1.1e-15 and 1.5e-14;
  - scalar n=2,3,4 (zonal sin((n+1)χ)/sin χ): 1.0e-14, 1.8e-14 and 6.4e-15.

  The archive is `g1_symbolic.json`.
- **G2, unused equations along solutions.** Vector ij: 7.6e-14. Scalar trace
  equation: 4.9e-16.
- **G4, determinants.** |det M − 1| ≤ 4.7e-12 in all sectors.
- **Phase independence.** Traces agree to ≤ 3.3e-13 for η₀ = .3,
  n = 2..10.
- **Controls.**

  | control | result |
  |---|---|
  | C1 free conformal field, max D_n | 2.0e-11 |
  | C2 vector with the metric coupling removed | identical to C1 |
  | C3 tensor n=2, tr M(π/2) | −0.0963065401946, against #294's −0.0963065402 (error 5.4e-12) |
  | C4 bare static tensor, phase error | ≤ 9.0e-11 |

C3 is an independent reproduction of #294's published map by a different
formulation, not a new prediction.

## 3. Why V and S are UNRESOLVED

G3 has three parts:
- DOP853 at 1e-12 and at 1e-10 must agree to 1e-7;
- RK4 at 2¹⁷ steps must agree with the primary DOP853 to 1e-7;
- the RK4 traces must converge with ratio in [8,32].

The first two pass everywhere. The largest disagreements are 9.1e-10
(DOP853 at the two tolerances) and 3.4e-12 (RK4 against primary). The ratio
condition fails:

| sector | degrees | ratio |
|---|---|---|
| V | n = 72–80 | 4.4–5.0 |
| S | n = 69–80 | 3.0–6.6 |

There, both RK4 differences from the primary are about 1e-11 to 1e-12. That
is the primary DOP853's own error, not RK4 truncation error, so the ratio
measures the reference's noise. The freeze did not name a floor below which
the ratio is not evaluated. The implementation evaluated it only when the
coarser difference exceeded 1e-11. The same failure pattern occurred in
#307.

The frozen rule makes these sectors UNRESOLVED, and they are reported as
such.

**Extension outcome.** The prospective extension
([`9a8ef99`](esu_floquet_g3_extension_prereg.md)) ran RK4 at 2^10, 2^11 and
2^12 steps against the unchanged archived maps
(`g3_extension.json`):

| sector | evaluated trace ratios | range | result |
|---|---|---|---|
| T | 76 | 18.3–31.8 | pass |
| S | 79 | 16.0–30.7 | pass |
| V | 79 | 0.84–15808 | fail (61 degrees outside [8,32]) |

Per the addendum, V stays UNRESOLVED and no further search follows.

As a diagnostic only, the V **matrix-norm** errors at the same step counts
converge at 16.0 and 16.0 for n = 18, 22, 27, 40, and at 15.7–16.0 for
n = 60, 80. The signed trace errors change sign between resolutions. The
refocusing maps are close to ±I, and the trace of a near-rotation is
stationary in its phase, so trace errors cancel at first order. The vector
maps therefore appear accurate. The registered criterion cannot show it,
and a matrix-norm criterion would need its own prospective registration.

## 4. Tensor sector (certified)

| n | F_n | θ_n |
|---|---|---|
| 2 | 0.995363 | 0.0963 |
| 3 | 0.997447 | 0.0715 |
| 5 | 0.998877 | 0.0474 |
| 10 | 0.999668 | 0.0258 |
| 20 | 0.999909 | 0.0135 |
| 40 | 0.999976 | 0.00691 |
| 80 | 0.999994 | 0.00350 |

Every multiplier lies on the unit circle; the maximum |μ| is 1 − 6e-14.

- The measured (n+1)θ_n over n = 40–80 is 0.283349. WKB, with
  `m_T² = <2R²/f> − <H²> − 1 = 0.828540 − 0.008919 − 1 = −0.180379`,
  predicts π|m_T²|/2 = 0.283339.
- The bare static ESU (m² = −1) would give θ₂ = 0.539. The supporting fields'
  gradient stress adds `2R²/f`, which nearly cancels the curvature
  detuning.

So a gravitational-wave packet launched anywhere reassembles, inverted, at
its antipode after conformal time π (Einstein proper time ≈ .968π). The
inversion is the point-caustic sign. The residual dispersion falls as
n⁻². In the RP3 quotient this is a return to the launch point.

## 5. Vector sector (maps computed, not certified)

- All multipliers lie on the unit circle for n = 2..80.
- F_n: 0.5715 at n=2, 0.9921 at n=10, 0.99987 at n=80.
- 1 − F ~ n^−1.970.
- (n+1)θ_n = 1.30440, against the WKB value π<2R²/f>/2 = 1.30147 (0.23%).

Low-degree vector perturbations of the quartet refocus poorly.

## 6. Scalar sector (certified by the G3 extension)

**n=2 is hyperbolic.** Over one transit π, the multipliers are 1.2475444,
0.8015747 and a unit-circle pair 0.4150 ± 0.9098i, with det = 1. The
half-period map has real multipliers −1.1169 and −0.8953. This is period
doubling against the π/2-periodic pump: a first-tongue parametric resonance.

The result is robust within the frozen system:
- rtol 1e-13 reproduces it;
- a start phase of 1.1 reproduces it;
- DOP853 at 1e-10 differs by 3.8e-11 in the trace;
- RK4 differs by 3e-13;
- G1 verifies the n=2 reduced equations against the full linearization to
  1e-14.

Its growth rate is ln(1.2475)/π = 0.070 per unit conformal time. The
homogeneous Eddington mode (#296) grows at √2. n=2 scalar harmonics are
RP3-compatible (odd fields, invariant metric), so the inherited antipodal
restriction does not remove this mode. Degrees 3..80 are elliptic.

This mode was not predicted. The freeze made no stability prediction. It is
consistent with the α–β coupling splitting the frequencies to about
`n+1 ± 1/√f`. At n=2 the lower branch sits near 2, which is exactly half
the pump frequency 4.

**Refocusing saturates.**
- F_n: 0.7198 at n=2, 0.9870 at n=10, 0.9926 at n=20, 0.9943 at n=80.
- At n=80 the eigenphase errors are 0.088 and 0.123 rad. Their mean,
  0.105, matches the freeze's non-binding heuristic π(<f^−1/2> − 1) = 0.1057.

The frozen ASYMPTOTIC rule needs a positive fitted exponent and F ≥ .99 for
n = 40–80. A slow residual drift satisfies it: p = 0.157 ± 0.007. Under the
extension, the frozen label is therefore ASYMPTOTIC. The data, however,
describe saturation below 1, i.e. a plateau. This weakness of the frozen
classifier is recorded here; the label is not changed.

## 7. Operational choices made after the freeze (before any result was viewed)

1. **"Converges to a limit" for PLATEAU.** The spread of F_n over n = 60–80
   must be below .01. This was added when a unit test showed that an
   alternating sequence would otherwise be called PLATEAU.
2. **RK4 ratio floor.** The ratio is evaluated only when the 2¹⁶-step
   difference exceeds 1e-11 (see section 3).
3. **G1 recording.** G1 runs as a separate module
   (`esu_floquet_symbolic.py`, about 40 min) and is archived alongside the
   probe. It is not coded into the probe's verdict function. All G1 cases
   pass.
4. **First run discarded.** A first probe run was stopped before it printed
   any sector result, because the probe source had been edited after launch
   and its recorded hashes would not have matched. The archived run is the
   complete rerun under the committed sources.

## 8. What this does and does not establish

- **Linear antipodal refocusing of gravitational waves is a derived property
  of this supported ESU.** It needs no throat. It answers where and when a
  later two-object experiment should look: at the antipode, at η = π.
- **The background is not a stable stage.**
  - The homogeneous Eddington mode multiplies perturbations by 85 per
    transit.
  - The scalar n=2 quadrupole, certified by the extension, grows by 1.25 per
    transit.
  - Any experiment lasting a transit time has to account for both.
- **Not addressed:**
  - nonlinear evolution;
  - interaction with handles or mouths;
  - momentum exchange between objects;
  - quantization.
- **No discreteness mechanism.** The integer free spectrum and the
  refocusing are kinematic consequences of compactness.

## Reproduction

```sh
python -m experiments.closure_ledger.esu_floquet_probe \
  --output /tmp/esu_floquet.json            # several minutes
python -m experiments.closure_ledger.esu_floquet_probe \
  --output /dev/null --replay experiments/closure_ledger/runs/20260926_esu_floquet/esu_floquet.json
python -m experiments.closure_ledger.esu_floquet_symbolic /tmp/g1.json   # G1, ~40 min
python -m experiments.closure_ledger.esu_floquet_g3_extension --output /tmp/g3_extension.json
pytest -q tests/test_esu_floquet.py
```

Replay checks the freeze, the source hashes and every stored matrix's
derived observables and labels, and recomputes the maps (all degrees by
default). A tampered matrix or relabelled verdict clears the evidence.
