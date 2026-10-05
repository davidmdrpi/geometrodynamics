# Specification: resonance breaking at the LRS 2/5 crossing

Date: 2026-10-05. Parent: `main` at `3304cc5` (#319 merged).
Published before any map evaluation at or near the 2/5 crossing, and before
any phase scan of the #319 loop. The code is committed with this document.
Corrections will be dated notes only.

## 1. Question

#319 found that the two-return solutions of the diagonal-circular mode form a
**closed loop** of period-two points, with multipliers identical at all 12
samples to about 1e-10. In a generic area-preserving map, a resonant invariant
circle does not survive as a continuum of periodic points. It breaks into a
Birkhoff chain: q elliptic and q hyperbolic period-q orbits, whose residues
have opposite signs. An unbroken resonant circle needs a continuous symmetry
or an extra integral, at least locally. #320 proposes a nonlinear study in a
"symmetry-reduced neighbourhood". That is well-posed for the family direction
only if this direction is an exact symmetry.

**Question.** Is a resonant circle of the same homogeneous system, with no
known protecting symmetry, a continuum of periodic points or a Birkhoff
chain, at a resolution that the #319 loop passes?

**Test case.** The LRS (linearly polarised) invariant-circle family of
`r3_extension` Part A, at its crossing of rotation **2/5** (frac(ρ) = 0.4).
- It lies between the archived circles a = .2348 (ω = 2.57394) and
  a = .2792 (ω = 2.37961).
- q = 5 is the lowest closure denominator in the family's measured range.

**Control.** The #319 diagonal-circular two-return loop (q = 2), scanned by
the same method.

## 2. Hypotheses and interpretation

| main (LRS 2/5) | control (#319 loop) | reading |
|---|---|---|
| BROKEN_CHAIN | UNBROKEN_LOOP | The LRS sector has no extra integral that protects resonant circles at this level. The #319 loop's continuum is therefore not explained by integrability of the homogeneous model. Any explanation must be specific to the diagonal-circular family (see §7c). The family tangent should not be assumed to be an exact symmetry, and a nonlinear quotient by it is not yet defined. |
| UNBROKEN_LOOP | UNBROKEN_LOOP | Resonant circles survive as continua in a second, unrelated family, at the stated resolution. This is evidence for a hidden integral or local integrability. The next step would be to identify it. |
| BROKEN_CHAIN | BROKEN_CHAIN | The scan resolves breaking of the #319 loop that #319's tests did not. That is reported as a finding about the loop. |
| any INDETERMINATE | | Reported as such; no reading. |

My prior is BROKEN_CHAIN for the main case. Both outcomes are informative.

## 3. Method (`geometrodynamics/waves/r3_breaking.py`, `experiments/closure_ledger/r3_breaking_probe.py`)

**Maps (unchanged).**
- Main: the 4D LRS section map `r3_return_map.esu_map`, z = (A, p_A, x, p_x),
  with DOP853 at 1e-12/1e-14.
- Control: the 6D diagonal map `r3_family.P` on z[:6] (12D embedding), with
  the same integrator and tolerances.

**Curves.**
- Main: the interpolated circle K* = K₁ + s(K₂ − K₁), with
  s = (ω₁ − 4π/5)/(ω₁ − ω₂). K₁ and K₂ are the archived 63-point circles
  at a = .2348 and .2792. Both are normalised by the same real first x
  harmonic, so their phases are consistent. The curve is evaluated by
  trigonometric interpolation. Node seeds are K*(φ + 4πi/5), i = 0..4.
- Control: periodic cubic splines through the 110 archived #319 loop
  points, for node 0 (v[:6]) and node 1 (v[6:]), on a shared parameter.

**Phase scan.** At each of N = 60 node-0 phases φ_j = 2πj/60, solve the
square system in X = (z_0, …, z_{q−1}, λ):

    P(z_i) − z_{i+1} = 0               i = 0..q−2
    P(z_{q−1}) − z_0 − λ g = 0
    c'(φ)·(z_0 − c(φ)) = 0

- g = Ωᵀc'(φ)/|·|, where Ω is the pairwise symplectic matrix of
  `r3_return_map.OMEGA`. For a symplectic map, g is the cokernel direction of
  the closure at a continuum of periodic points.
- Newton uses centred finite-difference Jacobians (h = 1e-7). It stops at a
  max residual of 1e-12, with at most 12 iterations. A phase converges if its
  final residual is ≤ 1e-11.
- λ(φ) = 0 exactly where a period-q orbit crosses the phase-φ slice.
  - A continuum gives λ ≡ 0.
  - A Birkhoff chain gives λ with 2q sign changes per turn and dominant
    Fourier harmonic q.
- λ is in section-coordinate units along the unit vector g.

**Noise.** At phases j = 0, 10, 20, 30, 40, 50 of each scan, re-solve with
Radau (1e-12/1e-14) by chord Newton, keeping the DOP853 Jacobian, for at
most 4 iterations. Then

    noise ν = max|λ_Radau − λ_DOP853| + max(final Radau residual)

**Isolated orbits (main scan only, if it has sign changes).**
- From every sign-change bracket, interpolate the bracketing node sets
  linearly to the zero of λ. Then solve the unconstrained closure
  P(z_i) = z_{i+1 mod 5} (λ = 0, no phase row) by least-squares Newton.
- Converged orbits (residual ≤ 1e-11) are de-duplicated up to cyclic shift
  (tolerance 1e-8).
- For each distinct orbit, report:
  - the smallest singular value of the closure Jacobian;
  - the residue R = (2 − tr C)/4 of the centre block C. C comes from
    periodic orthogonal iteration over the five one-return Jacobians, so the
    hyperbolic multipliers (≈ 85 per return) are never multiplied into the
    centre.
- Residues use finite-difference Jacobians and exact jets
  (`r3_family.DP`, dims = 4, steps (2048, 4096), with (1024, 2048) as
  comparison). They are descriptive.

## 4. Registered labels

The rule is `r3_breaking.classify` plus the probe's gates. Let Λ = max|λ| over
the converged phases, r = max(10ν, 1e-11), and S the number of cyclic sign
changes of λ.

- **BROKEN_CHAIN:**
  - Λ ≥ 10r, and S ≥ 2q with S a multiple of 2q;
  - and, for the main case, at least two distinct isolated orbits converge,
    one from an upward and one from a downward λ crossing.
- **UNBROKEN_LOOP:** Λ ≤ r, and r ≤ 1e-9.
- **INDETERMINATE:** everything else. This includes more than 6 of the 60
  phases failing to converge.

The control gets its label from the same rule with q = 2. The orbit gate does
not apply to it. Each label is reported separately; neither overrides the
other.

Also reported:
- Λ and Λ/a, where a ≈ .249 is the interpolated circle size;
- the Fourier spectrum of λ;
- residual, constraint and conditioning maxima;
- the maximum distance from the seeds;
- the residues of the distinct orbits.

**No retries.** If a stage fails, the label is INDETERMINATE and the failure is
reported. Any rerun with changed settings will be a dated, post-hoc correction.

## 5. Validation before freeze (toy only)

`tests/test_r3_breaking.py`, 8 tests. The toy is a 4D map: a hyperbolic pair
with multiplier 85, times a twist map with a q = 5 resonant kick
ε·sin(5φ).

| check | result |
|---|---|
| ε = 0 | λ ≡ 0 (< 1e-13) from seeds scaled by 1.001: UNBROKEN_LOOP |
| ε = 1e-9, 1e-6 | BROKEN_CHAIN with S = 10 and dominant harmonic 5; Λ equals the analytic 5ε/R to 1e-3 |
| chain residues | opposite signs, magnitude equal to the analytic value to 2% |
| centre block | recovers the product multipliers when each factor carries a multiplier-85 hyperbolic pair: elliptic and Jordan cases to 1e-10 |
| curve interpolants | exact on band-limited data |
| full probe pipeline on the toy | runs through to scoring |

## 6. What will not follow

- **Neither label bears on action selection, or on stability against
  inhomogeneous perturbations.**
- **UNBROKEN_LOOP is a bound at resolution r, not a proof of exactness.**
- **BROKEN_CHAIN at LRS 2/5 does not show that the #319 loop is broken.** The
  control answers that separately.
- **The scan is local to the interpolated circle.** It does not establish
  whether the family continues beyond a = .3044 (Part A, unresolved).

## 7. Disclosures

a. **Observation during the #320 review (2026-10-05).** This was computed
   before this freeze and is descriptive only. Re-integrating the 12
   archived #319 loop points shows:
   - the two-return time is identical to 10 digits (6.2716152504);
   - the first-return time varies by 3e-5;
   - A at the section varies by 7e-5.

   So the loop points are not related by a rotation of the diagonal shear
   plane. This observation motivated the experiment. It is not a test here.

b. **Correction to my #320 review note.** It said "rotation 2/5 near
   a ≈ .24, in the diagonal subsystem" as if that rotation lay on the circle
   family that holds the #319 loop. The 2/5 crossing belongs to the **LRS**
   family, whose frac(ρ) decreases from .4847 to .3587. The diagonal-circular
   family's frac(ρ) increases from .4847 through 1/2. Its next low-order
   rationals lie beyond the computed range, a ≤ .0905.

c. **Why not a second resonance of the circular family?** This is a heuristic
   and is not tested. The circular mode cycles the anisotropy through the
   three axes, so the cyclic axis permutation acts as a phase shift of 2π/3.
   Breaking at rotation p/q would then need tensor-phase harmonics that are
   multiples of lcm(3, q): 6 at q = 2, 15 at q = 5, 21 at q = 7. Higher
   resonances of that family would be at least as protected as 1/2, so they
   cannot discriminate. The LRS family has no such restriction.

d. **Pre-freeze computation on the model.** Pre-freeze computation was
   limited to:
   - timing one evaluation of each map, and of Radau, at the
     non-resonant point (A, p_A, x, p_x) = (1, 0, .01, 0);
   - reading the archived circles (Part A) and loop (#319).

   No map evaluation was made at or near the 2/5 crossing.

e. **Known archive fact.** In Part A, the first attempt at a = .256, which
   lies near this crossing, failed with "step size underflow" from a
   rescaled predictor. The retry at .2348 passed. This is not used here.

**Runtime estimate.** The scans take a few minutes each on four processes.
The Radau noise estimate takes about five minutes. The jet residues take
about four minutes per distinct orbit.

## 8. Reproduction

    python -m pytest -q tests/test_r3_breaking.py
    OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.r3_breaking_probe

Archives: `experiments/closure_ledger/runs/20261005_r3_breaking/`
(`scan_control.json`, `scan_main.json`, `noise.json`, `orbits.json`,
`result.json`). Each archive binds the SHA-256 of the seven sources. An
authenticated replay will be added with the results.
