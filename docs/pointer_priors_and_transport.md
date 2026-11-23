# What selects the record in the two BAM measurement models?

With the elliptic mixer followed by a winding Stern–Gerlach, a localized position transported from spatial equilibrium gives the expected cosine-squared curve. A thermal mechanical latch and a phase-prepared latch give different curves. Both models retain a field in both channels. They assign the single-trial record to different physical variables.

The [protocol](pointer_and_waist_protocol.md) was published at [a50f290](https://github.com/davidmdrpi/geometrodynamics/commit/a50f2901d26a634a0f6b5b8d728277bada14952d). The separate classical-support result and its assumptions are in the [asymmetric waist analysis](dilatonic_waist_characteristics.md).

## 1. Natural preparations for the fixed latch

Use the action, damping, gain, time and classification of [#357](elliptic_pointer.md). Let $z=\cos(2\theta)$ and $w_*=(1+\rho)g/\rho=0.16981365$. The exact half-period propagator from arbitrary $(q_0,v_0)$ gives $q_T=(1+\rho)gz-\rho q_0$, $v_T=-\rho v_0$. The probe evolves both coordinates through the nonlinear latch. In particular, velocity matters for a genuine thermal preparation. The Gaussian is the Boltzmann law of the uncoupled unit harmonic oscillator before connecting the field; the load has no incident waves during the write or latch. A continuously thermalized detector would require fluctuating forces.

| Preparation | Positive fraction at z=0.5 (target 0.75) | At z=1 (target 1) | Largest absolute error, nine inputs |
|---|---:|---:|---:|
| Boltzmann, half reference rms | 0.9102 | 0.9958 | 0.1602 |
| Boltzmann, reference rms $w_*/\sqrt3$ | 0.7490 | 0.9102 | 0.0898 |
| Boltzmann, twice reference rms | 0.6301 | 0.7490 | 0.2510 |
| Uniform phase, cosine displacement, released at rest | 0.6665 | 1.0000 | 0.1050 |
| Uniform phase of freely moving harmonic oscillator | 0.6265 | 0.7813 | 0.2188 |
| Haar sphere axial projection, released at rest, width $w_*$ | 0.7500 | 1.0000 | 0.0000 |

These are 4,096 equal-weight launches per input. The 1,024-launch check gives the same decisions. All launches in these runs latched. Sobol points approximate the Gaussian joint measure; midpoint angles integrate uniform phase. The table measures a deterministic model with specified preparations, rather than statistical error in an observed experiment.

Two analytic facts delimit the negative result.

* At every positive temperature the Boltzmann distribution has positive density on an open set of initial states that enter the **wrong** basin even at z=1. Take $q_0$ slightly above $w_*$ and $v_0$ near zero. The switch position is then strictly negative and the energy below the barrier. Continuity gives an open wrong-record set, so exact certainty at the endpoint is impossible at any such temperature.
* At rest, $q_0=w_*\cos\phi$ gives $P_+=1/2+\arcsin(z)/\pi$. More generally a smooth nonconstant scalar function of one uniform circle angle cannot have an exactly uniform bounded displacement law. At a maximum its derivative vanishes. For arbitrarily small slopes near the maximum, the fraction in the top displacement interval divided by the interval width is arbitrarily large. A regular uniform density has a finite ratio. This argument concerns a smooth one-angle displacement and the #357 at-rest basin; phase-dependent velocities or a different coupling need separate dynamics.

A symmetry-based rescue exists in principle: one Cartesian projection of a Haar-isotropic unit vector on $S^2$ is exactly uniform, since the area element is $d\phi\,d(\cos\vartheta)$. It gives the Born curve for a linear, at-rest displacement of matched width. #241's ellipse supplies a circle orientation; it has yet to supply an isotropic sphere preparation, the displacement map and its gain. As a further control, averaging the actual eight-site mixer over 128 uniformly spaced orientations gave positive-force weight 0.5 for every equatorial input. Randomizing the analyzer angle removes the external phase reference.

## 2. The ellipse as a pulse before a spatial measurement

The eight-site elliptic metric hopping from [#241](sigma_z_readout_capstone.md) evolves a coherent winding superposition with equal winding intensities and relative phase $2\theta$. A pulse about the negative y axis maps its x quadrature to winding z. With ellipticity ε the duration is fixed by the leading matrix element: $T=\pi/[4(2-\sqrt2)\epsilon]$. We diagonalize the full ring Hamiltonian at ε=0.02 and 0.01. There is no angle-by-angle adjustment and no threshold drawn from an apparatus prior. The analyzer still requires a pulse area and a phase reference.

After switching off the mixer, the wave enters the co-moving force window from [#209](measurement_sector.md). Its odd-winding force is proportional to $\sin(k\pi/4)$. The quadrupole keeps the odd windings, so +1,+3 have the same spatial envelope and −1,−3 share its mirror. For each initial position we integrate the positive norm-current velocity $J/\rho$ with $\rho=\sum_k|\psi_k|^2$, $J=\sum_k\operatorname{Im}(\psi_k^*\partial_x\psi_k)$. The Gaussian envelopes are the exact accelerating and spreading solutions of this linear-force window. We keep the finite wrong-side mass at read time. This combines a spatially uniform mixer with a subsequent prescribed analyzer. It does not model a two-mouth apparatus, backreaction or an irreversible spatial amplifier.

For 2,048 deterministic position quantiles drawn from spatial equilibrium, the largest difference between the fraction arriving on the right and $\cos^2\theta$ was **0.0004883** at ε=0.02. The maximum error of the continuous right-side mass was **0.0002768**; at ε=0.01 it fell to **0.0002254**. The single-branch wrong-side Gaussian tail at the chosen read time was 0.0002092. Quantile counting can round away this finite error. At ε=0.02 the largest leakage to |k|=3 was 0.000173. The integrated trajectories preserved the cumulative norm probability to within $1.7\times10^{-10}$.

| Preparation and controls | Right fraction at θ=π/8 | At θ=π/4 |
|---|---:|---:|
| Expected $\cos^2\theta$ | 0.85355 | 0.50000 |
| Ellipse + spatial equilibrium | 0.85352 | 0.50000 |
| Same fields, initial positions at 0.6 times the equilibrium width | 0.95996 | 0.50000 |
| Mixer off, equilibrium | 0.50000 | 0.50000 |
| Pulse about x, equilibrium | 0.85352 | 0.99951 |

The narrower preparation changes the right fraction even though the fields and trajectories' velocity law are identical. Equivariance transports an equilibrium measure forward. It does not establish that every initial ensemble has that measure. The unoccupied branch here means **no position beable in that branch in a particular run**. Its field persists, and while the two envelopes overlap it still contributes to the local velocity.

### Which current is supplied by #199?

For a single winding k in the weak-field reduction of [#199](guidance_law_from_5d.md), the fiber-momentum density and current are proportional to $k\rho_k$ and $kJ_k$. The factor k cancels in the velocity. For several windings the signed sums $\sum k\rho_k,\sum kJ_k$ differ from the positive norm density and its current. As a direct counterexample, take equal k=±1 branches at x=0 during the Stern–Gerlach transit. By reflection their densities are equal and their velocities are opposite. At t=2 the signed density vanishes while the signed current divided by the positive density is **0.970588** in code units. The signed-current velocity is undefined; the positive-current velocity is zero. Moreover the ellipse breaks continuous fiber rotations during the mixing pulse.

Hermitian wave evolution gives the positive norm continuity equation. Spatial equilibrium plus #209's guidance by that current gives the observed result. Deriving the localized geometric core's mixed-channel law directly from five-dimensional dynamics remains an additional step beyond the single-winding cancellation. [#198 equivariance proof](born_rule_equivariance.md)

## Reproduction

From the repository root, with NumPy, SciPy and SymPy:

~~~bash
export OPENBLAS_NUM_THREADS=1
python -m experiments.closure_ledger.antipodal_network.natural_priors
python -m experiments.closure_ledger.antipodal_network.elliptic_transport
python -m experiments.closure_ledger.antipodal_network.dilatonic_waist
python -m unittest discover -s tests -p test_pointer_and_waist.py
~~~

[Archived results](../experiments/closure_ledger/runs/20261123_pointer_and_waist/)
