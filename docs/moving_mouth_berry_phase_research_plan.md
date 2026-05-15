# Moving-mouth Berry phase falsification test — research plan

Opens a new research thread orthogonal to the closure-ledger and
throat-dynamics work. The closure-ledger sequence (PRs #11–22)
established the *static* eigenvalue structure of the BAM lepton
sector and the *non-local* WKB-BS reading of the inner cutoff.
This thread tests the *adiabatic* (Berry phase) content of the
throat-transport operator — a sharp falsifiable prediction that
exists independently of the closure-ledger machinery.

## The prediction

The BAM throat carries the Hopf connection `A = ½ cos(χ) dφ` on
the S² base of the Hopf bundle (`hopf/connection.py`), and the
non-orientable transport operator `T = iσ_y` with `T² = −I`
(`embedding/transport.py`). Together these encode the spin-½
double cover and the closure-quantum partnership at χ = 0.

If the throat **mouth** is adiabatically transported through a
closed loop γ in the (χ, φ) configuration space, the wavefunction
at the mouth picks up a Berry phase

```
γ_Berry  =  ∮_γ A  =  ∬_S F  =  −Ω(S) / 2
```

where Ω(S) is the solid angle enclosed by γ on the Hopf S² and
the second equality follows from Stokes + the curvature
`F = dA = −½ sin(χ) dχ ∧ dφ`. The factor of ½ is the spin-½
signature.

This is a **uniquely BAM** prediction in the following sense: the
explicit connection `A = ½ cos(χ) dφ` is the orientation-reversing
isometry of S³ that preserves the Hopf bundle (derivation in
`embedding/transport.py`). The Berry phase computed from this
specific connection must equal `−Ω/2` for any closed loop, with
errors only at the level of numerical discretization. Any
disagreement at finite resolution that does not converge as
`O(1/N²)` would falsify the BAM-derived Hopf connection.

## What "moving mouth" means concretely

The throat is a finite-radius region of S³ that, for the static
eigenvalue analyses, sat at a fixed (χ = 0, φ = 0). For this
thread, we let the mouth be adiabatically translated/rotated in
the bundle-base coordinates (χ, φ) along a closed trajectory γ(t):

  - **Trajectory γ(t)**: a closed path γ: [0, T] → S² with
    γ(0) = γ(T).
  - **Mouth state**: parametric family `|ψ(γ(t))⟩` of
    instantaneous-ground-state spinors at each point of γ.
  - **Berry connection**: pull-back of the canonical Hopf
    connection to γ. For our problem, this is just the standard
    `A = ½ cos(χ) dφ` evaluated along γ.
  - **Berry phase**: `γ_Berry = ∮_γ A = i ∮_γ ⟨ψ|dψ⟩` along the
    discretized trajectory.

The two formulas — line integral of A and overlap integral of
⟨ψ|dψ⟩ — must agree to high precision. Disagreement at the
numerical level would indicate either (a) an inconsistency in the
BAM-derived A, or (b) a bug in the spinor parametrisation. Both
are falsifying.

## Three tests with closed-form predictions

### Test 1: Equatorial loop at fixed χ = π/2

Path: φ goes from 0 to 2π at χ = π/2 (constant). The loop is a
single Hopf fibre at the equator.

Prediction: `A_φ(π/2) = ½ cos(π/2) = 0`, so `γ_Berry = 0`. The
equatorial Hopf fibre carries zero Berry phase — the equatorial
orbit is the "stable" zero-self-energy configuration in
`connection.py`.

### Test 2: Polar cap up to χ_max with full φ winding

Path: a full loop at fixed χ = χ_max (boundary of the cap from
χ = 0 to χ = χ_max).

Prediction: enclosed solid angle on the Hopf S² is
`Ω = 2π·(1 − cos χ_max)`, giving

```
γ_Berry = −π · (1 − cos χ_max)
```

Specific values:

  - χ_max = π/2 (hemisphere): γ_Berry = −π → the spin-½ sign flip.
  - χ_max = π (full sphere): γ_Berry = −2π → the closure quantum.
  - χ_max → 0 (small loop near pole): γ_Berry → −π·χ_max²/2 → 0
    quadratically.

### Test 3: 4π double-cover loop

Path: traverse a 2π loop TWICE (φ goes from 0 to 4π at χ = π/2,
or equivalently traverse the polar-cap loop twice).

Prediction: SO(3) Berry phase is doubled (mod 2π); SU(2) Berry
phase is `2 × (−π·(1 − cos χ_max))`. For χ_max = π/2 (the spin-½
configuration), the SU(2) phase is `−2π`, equivalent to `+1` in
exp(iγ) — recovering the spinor 4π-periodicity (T² = +I after
two traversals). This is the dynamical version of the static
spinor-monodromy check in `hopf/spinor.py`.

## Falsification criterion

The probe computes both numerical and analytic Berry phases for
each test. The framework PASSES the test if

```
|γ_Berry_numerical − γ_Berry_analytic|  <  10⁻⁶
```

(machine-precision agreement at high path-discretisation N ≥ 1024).

The framework FAILS if any test shows a disagreement larger than
the discretisation error `O(1/N²)`. Failure modes that would
falsify BAM's geometric content:

  - Non-zero Berry phase for the equatorial loop (test 1).
  - Wrong cosine dependence on χ_max in test 2.
  - Berry phase not doubling for the 4π loop in test 3 (would
    indicate the Hopf bundle is trivial, contradicting `T² = −I`).
  - Discretisation error scaling worse than `O(1/N²)` (would
    indicate a singular structure in A that the BAM derivation
    missed).

## Sub-targets

### (1) Numerical Berry phase via line integral of A

Discretise the trajectory γ at N steps. Compute the line integral
`∮ A·dγ` directly from the BAM connection. Compare to the
analytic prediction. This is the most direct test.

### (2) Numerical Berry phase via overlap integral of ⟨ψ|dψ⟩

Construct an explicit spinor parametrisation `|ψ(χ, φ)⟩` of the
instantaneous Hopf-fibre eigenstate. Compute the Berry phase as
`γ = i ∮ ⟨ψ|∂_φ ψ⟩ dφ + i ∮ ⟨ψ|∂_χ ψ⟩ dχ`. Verify equal to the
line integral.

### (3) Convergence with N

Verify the discretisation error scales as `O(1/N²)`. Anomalous
scaling (e.g. `O(1/N)` or constant residual) would indicate
either a singular point in A on the trajectory or a flawed
discretisation.

### (4) Holonomy ↔ closure-ledger consistency

The static Hopf-fibre holonomy `π·cos(χ)` from `connection.py` is
the closed-form integral of A over a single full-φ loop at fixed
χ. Verify the moving-mouth Berry phase reduces to this static
holonomy in the appropriate limit (constant-χ trajectory). This
is a sanity-link to the closure-ledger framework.

## Stopping condition

The thread closes when:

  (a) All four tests pass at machine precision. The BAM-derived
      Hopf connection produces the predicted Berry phases for any
      closed mouth trajectory. The spinor double cover and the
      closure-quantum 2π are confirmed dynamically.

  (b) Any test fails. The BAM derivation of `A = ½ cos(χ) dφ` is
      contradicted by the numerical Berry phase, indicating a
      flaw in the geometric foundations.

The outcome is a clean falsification test, not a derivation —
either BAM's geometric content survives or it doesn't.

## Cross-references

- `geometrodynamics/hopf/connection.py` — Hopf connection
  `A = ½ cos(χ) dφ`, curvature, and full-fibre holonomy
  `π · cos(χ)`. The probe extends this to general trajectories.
- `geometrodynamics/hopf/spinor.py` — static spinor monodromy
  (2π → −1, 4π → +1). The probe is the dynamical / adiabatic
  version.
- `geometrodynamics/embedding/transport.py` — derivation of
  `T = iσ_y` as the orientation-reversing Hopf-preserving map;
  the source of the spinor double-cover content.
- `experiments/closure_ledger/moving_mouth_berry_phase_probe.py`
  — first probe in this thread.
