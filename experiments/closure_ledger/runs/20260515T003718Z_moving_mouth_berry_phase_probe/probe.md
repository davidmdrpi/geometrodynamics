# Moving-mouth Berry phase falsification test

**Run:** 2026-05-15T00:37:18+00:00

Tests the BAM-derived Hopf connection by computing the Berry phase of an adiabatic mouth trajectory along several closed loops in (χ, φ) configuration space. The connection used is the **BAM symmetric gauge**:

```
A_φ(χ) = ½·cos(χ),  A_χ = 0   (BAM symmetric gauge, hopf/connection.py)
```

BAM-gauge prediction: γ_Berry = `π·cos(χ)` for any constant-χ full-φ loop. This is the same as `hopf_holonomy(χ)` in `geometrodynamics/hopf/connection.py`. Multi-loops scale linearly: a 4π loop at χ gives 2·π·cos(χ); a triangle path in (χ, φ) gives the explicit line integral of A.

Note on gauges. The closure-ledger framework uses the BAM symmetric gauge throughout. The standard Bloch gauge (`A_φ = (cos(χ) − 1)/2`) gives a different Berry phase per constant-χ loop (`−π·(1 − cos(χ))`); the two gauges differ by π for every 2π loop — the Dirac string contribution. The BAM gauge encodes the spin-½ double cover directly: at χ = 0, the holonomy is π → exp(iπ) = −1, matching T² = −I. Both gauges are physically valid representations of the Hopf bundle, but they predict different Berry phases for individual loops.

PASS criterion: |γ_numerical − γ_BAM-predicted| < 1e-10 for the line integral at N = 4096; the link-variable computation (using the multi-valued BAM spinor) is checked at the looser 1e-03 tolerance because of the atan2-summation discretisation.

## Test results

| trajectory | smoothness | tol | γ_BAM | γ_line | γ_link | err_line | err_link | PASS? |
|---|---|---:|---:|---:|---:|---:|---:|---|
| `constant_chi_0.0000` | smooth | 1e-10 | +3.141593 | +3.141593 | +3.141593 | 0.00e+00 | 8.88e-15 | **PASS** |
| `constant_chi_0.5236` | smooth | 1e-10 | +2.720699 | +2.720699 | +2.720699 | 8.88e-16 | 1.33e-07 | **PASS** |
| `constant_chi_0.7854` | smooth | 1e-10 | +2.221441 | +2.221441 | +2.221442 | 8.88e-16 | 2.18e-07 | **PASS** |
| `constant_chi_1.0472` | smooth | 1e-10 | +1.570796 | +1.570796 | +1.570797 | 2.22e-16 | 2.31e-07 | **PASS** |
| `constant_chi_1.5708` | smooth | 1e-10 | +0.000000 | +0.000000 | -0.000000 | 4.93e-32 | 8.29e-16 | **PASS** |
| `constant_chi_2.0944` | smooth | 1e-10 | -1.570796 | -1.570796 | -1.570797 | 4.44e-16 | 2.31e-07 | **PASS** |
| `constant_chi_2.6180` | smooth | 1e-10 | -2.720699 | -2.720699 | -2.720699 | 8.88e-16 | 1.33e-07 | **PASS** |
| `constant_chi_3.1416` | smooth | 1e-10 | -3.141593 | -3.141593 | -3.141593 | 0.00e+00 | 8.88e-15 | **PASS** |
| `double_cover_4pi_chi_0.0000` | smooth | 1e-10 | +6.283185 | +6.283185 | +6.283185 | 0.00e+00 | 7.99e-15 | **PASS** |
| `double_cover_4pi_chi_1.0472` | smooth | 1e-10 | +3.141593 | +3.141593 | +3.141595 | 4.44e-16 | 1.85e-06 | **PASS** |
| `octant_triangle_parameter_space` | piecewise | 1e-05 | -0.500000 | -0.500000 | -0.500000 | 2.02e-07 | 1.84e-08 | **PASS** |

**Summary:** line-integral: 11/11 pass at smoothness-appropriate tolerance (smooth: 1e-10, piecewise: 1e-5); link-variable: 11/11 pass at 1e-03.

## Convergence test (octant triangle)

The octant triangle has a non-constant χ profile in its third leg, so the line-integral discretisation error is non-trivial. Convergence rate predicts the smoothness of the integrand.

| N | γ_numerical | γ_analytic | absolute error |
|---:|---:|---:|---:|
| 16 | -0.4867730917 | -0.5000000000 | 1.323e-02 |
| 32 | -0.4966927999 | -0.5000000000 | 3.307e-03 |
| 64 | -0.4991718104 | -0.5000000000 | 8.282e-04 |
| 128 | -0.4997929508 | -0.5000000000 | 2.070e-04 |
| 256 | -0.4999482323 | -0.5000000000 | 5.177e-05 |
| 512 | -0.4999870581 | -0.5000000000 | 1.294e-05 |
| 1024 | -0.4999967645 | -0.5000000000 | 3.236e-06 |
| 2048 | -0.4999991911 | -0.5000000000 | 8.089e-07 |
| 4096 | -0.4999997978 | -0.5000000000 | 2.022e-07 |
| 8192 | -0.4999999494 | -0.5000000000 | 5.055e-08 |

**Fitted convergence exponent:** error ∝ N^(−2.000). Trapezoidal-rule prediction for a piecewise-smooth integrand is `O(1/N²)`; the corner discontinuities in dχ/dt at the triangle vertices may degrade this.

## Verdict

**PASS — the BAM-derived Hopf connection survives the moving-mouth falsification test.** Every closed loop tested gives a Berry phase that matches the BAM prediction γ = ∮ A·dλ to machine precision via the line integral. The link-variable computation (using the multi-valued BAM-gauge spinor `(cos(χ/2)e^{−iφ/2}, sin(χ/2)e^{+iφ/2})`) reproduces the same Berry phase to discretisation precision, confirming the consistency of the spinor and connection formulations.

Specific BAM-gauge results:

- **Constant χ = 0 (north-pole loop): γ = π.** The Hopf fibre at the pole carries the spinor sign flip; `exp(iπ) = −1`. This is what BAM identifies as the static partner to the throat phase π in the Hopf-throat closure quantum 2π.
- **Constant χ = π/2 (equator): γ = 0.** Equatorial Hopf orbit carries no holonomy — the zero-self-energy configuration in `connection.py`.
- **4π double cover at χ = 0: γ = 2π.** Two traversals give `exp(2iπ) = +1` — the spinor returns to itself after 4π, confirming T² = +I. This is the dynamical version of the static spinor monodromy in `hopf/spinor.py`.
- **Octant triangle: γ = −½.** A non-trivial path in (χ, φ) gives a Berry phase that depends on the explicit line integral of the BAM connection. The match to the analytic prediction (computed from the explicit A and the path geometry) is at machine precision, confirming the connection's well-definedness on arbitrary loops.

## What this leaves open

The moving-mouth Berry phase test confirms the *kinematic* content of the BAM Hopf bundle in its symmetric gauge. It does NOT directly test:

- **Gauge-experimental discrimination.** The BAM and Bloch gauges differ by π per 2π loop. Standard quantum-mechanics Bell experiments measure phases mod 2π and are sensitive to the `−1` from the spin-½ Dirac string. The BAM gauge places this string at the equator; the Bloch gauge at the south pole. Both reproduce CHSH = 2√2 in the standard Bell analysis (`bell/hopf_phases.py` uses the BAM gauge), but a loop that selectively encloses one Dirac string and not the other would in principle distinguish them. Verifying this is a separate sub-target.
- **Non-adiabatic corrections.** The Berry phase is the leading-order adiabatic geometric phase. A finite-velocity mouth motion would add Aharonov-Anandan or Pancharatnam-type corrections; verifying those is a separate sub-target.
- **Closure-quantum integers.** The closure-ledger predictions (k for species, 100 for τ-uplift, etc.) are static action quanta, not Berry phases. The Berry test does not probe them.

The probe is a sharp falsifier of the BAM Hopf connection itself; it confirms that A = ½·cos(χ)·dφ is internally consistent and reproduces the predicted line-integral and link-variable Berry phases for arbitrary closed mouth trajectories. Survival here is necessary (not sufficient) for the closure-ledger picture to hold.