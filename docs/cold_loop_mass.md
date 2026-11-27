# Center mass through the cold charged loop

The scalar-loaded reflection center keeps contracting in all four evolutions. In the closest-to-cold case the center crosses the reference inner radius at U = −0.01371. By the first return of a radial light ray around the circle, its mass has grown from 0.558392 to 317,351 while its areal radius is still 0.61177. This comparison uses an event defined by the periodic geometry that all three grids reach.

[Registered test](cold_loop_mass_protocol.md) · [Parent equations and exact vacuum control](charged_nariai_loop.md) · [Numerical archive](../experiments/closure_ledger/runs/20261126_cold_loop_mass/probe.json)

## Which mass is measured

In the units of #371, let U and w be the unit normal and unit spatial derivatives of areal radius. The invariant combination w²−U² is the squared gradient of r. The mass parameter used here is

$$
\mu=r^2(1+U^2-w^2-r^2/2)+Q^2/r^2.
$$

It equals the constant metric parameter in electrovac. In geometric Misner–Sharp units the unadjusted mass is proportional to r²(1+U²−w²). The cosmological subtraction and electromagnetic correction relating these two masses stay finite on the radius intervals in the experiment.

At the reflection center w vanishes. Differentiating μ with respect to its proper time and substituting the radial Einstein equation cancels the Λ and Maxwell terms, giving

$$
\dot\mu=-\tfrac13r^3 U(p^2+f^2).
$$

Here p and f are the orthonormal time and space derivatives of the neutral scalar. Scalar stress therefore increases μ during contraction. The code also integrates −r³(p²+f²) dr/3 by the trapezoid rule, starting at the exact initial mass. This checks the geometric mass using a different combination of evolved quantities from the same solver. It supplies a signed integral through a possible turning point without dividing by U.

## Replaying the Q = 0.14 acceleration spike

The table uses common areal radii on the primary packet. Both finite grids resolve the sign and the onset of large mass growth. The two finest estimates differ by 18% of the finest value at the deepest cut.

| r | U, 32,768 | μ, 16,384 | μ, 32,768 | Scalar integral, 32,768 |
|---:|---:|---:|---:|---:|
| 0.19800 | −0.299 | 0.541872 | 0.541882 | 0.541887 |
| 0.19700 | −0.602 | 0.557784 | 0.557174 | 0.557213 |
| 0.19650 | −23.86 | 24.56 | 22.52 | 22.60 |
| 0.19625 | −152.6 | 1,009 | 897.3 | 901.4 |
| 0.19600 | −974.0 | 42,419 | 36,456 | 36,666 |
| 0.19590 | −2,043 | 189,732 | 160,207 | 161,217 |

At r = 0.197 the proper-time acceleration of r² is −990.73, even though μ has risen only 2.8% since r = 0.198. Farther inward the geometric mass and scalar integral grow together. At r = 0.1959 their discrepancy on the fine grid is 0.63%, and the normalized null constraint falls from 0.00814 to 0.00283 on refinement. These observations pin down finite physical amplification. The proper clock approaches 2.06052 over this short window. Establishing exponential mass inflation at a Cauchy horizon would require identifying the horizon, normalizing an advanced clock and controlling the asymptotic flux there. The r = 0.2 vacuum inner horizon in #371 belongs to the separate Q = 0.1378695 control.

## Following the cold branch

Write a, b, c for the three squared horizon radii. Fix b = 0.4, take a = b(1−ε) and c = 2−a−b. Factoring F(r) gives Q²=abc/2 and μ=(ab+ac+bc)/2. The reference is cold as a approaches b. The charged Nariai degeneration occurs when b approaches c. All three surface gravities use |F′|/2 in the static metric normalization of #371.

Each packet starts at R = √b with amplitude 0.02ε. The Hamiltonian constraint sets P² = 6(R⁻²−1−Q²/R⁶). This path varies the scalar preparation along with Q. The initial center mass is the reference μ plus R²(0.02ε)². Values in the last two columns are measured as the packet crosses r = √a on the 2,048-point grid. These radii label comparison events in the dynamical solutions.

| ε | Q | κ₋ | κ₊ | κ cosmological | Initial μ | U at √a | μ at √a |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.50 | 0.236643 | 1.34164 | 0.39528 | 0.36221 | 0.460040 | −1.08244 | 0.694339 |
| 0.20 | 0.286217 | 0.21213 | 0.13914 | 0.29168 | 0.524806 | −0.31479 | 0.556510 |
| 0.05 | 0.304500 | 0.03586 | 0.03241 | 0.25558 | 0.551800 | −0.07037 | 0.553682 |
| 0.01 | 0.308799 | 0.00648 | 0.00636 | 0.24587 | 0.558392 | −0.01371 | 0.558466 |

The archive includes all registered radius cuts reached on each grid and the last valid stencil. None of the twelve cold-family evolutions records a negative-to-positive U crossing. The terminal stencil is a numerical limit.

Two supplementary diagnostics clarify the closest-to-cold run. First, the extremal mass at **fixed Q** follows from the smaller positive solution of Q²=s²(1−s), with μ_cold=2s−3s²/2. The initial mass exceeds μ_cold by only 0.0000040663. By r=√a the scalar has added 18.3 times that initial excess. The small reference κ₋ characterizes the preparation. The subsequent center mass demonstrably departs from it.

Second, on the uniform conformal grid radial light travels at unit coordinate speed and the circle has period 2π. A ray sent from the center on the initial slice returns at grid time 2π. This gives a common finite event even very deep in collapse.

| Points | τ at return | r | U | Geometric μ | Scalar integral | Normalized null constraint |
|---:|---:|---:|---:|---:|---:|---:|
| 512 | 0.778780 | 0.611767 | −922.46 | 318,466 | 318,395 | 6.47×10⁻⁶ |
| 1,024 | 0.778829 | 0.611772 | −921.15 | 317,574 | 317,557 | 1.63×10⁻⁶ |
| 2,048 | 0.778842 | 0.611774 | −920.83 | 317,351 | 317,346 | 4.07×10⁻⁷ |

The five archived spatial samples on that last slice span μ=275,200 to 317,351. U is negative at all five. Their masses are scalar invariants evaluated on the specified slicing. Away from the reflection center U is the normal-frame derivative.

The compact packet carries nonzero integrated scalar momentum. In conformal coordinates the scalar equation conserves the periodic integral of r³∂ₜφ, initially 2πR³P on the uniform grid. Local flux can redistribute this quantity. A decay rate at a cosmological exterior was never prescribed or measured here. Importing the stationary exterior spectral gap into this loop would leave the decisive flux assumption untested. The results disfavor a charge-only rescue along this particular constrained path. Other scalar data, including weak localized perturbations of the actual electrovac neck with zero integrated scalar momentum, remain open. [Higher-dimensional RNdS scalar perturbations](https://doi.org/10.1007/JHEP03(2019)187)

## Reproduction

From the repository root run the commands below. The worker count controls independent simulations. The archived spatial profiles supplement the registered center diagnostics.

~~~sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.gravity.cold_loop_mass --workers 4
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_*loop_mass.py'
~~~
