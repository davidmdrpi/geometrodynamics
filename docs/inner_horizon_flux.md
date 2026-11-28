# How much scalar flux reaches the cold inner horizon?

For the 0.5%-loaded shell of #373, the first nonzero mass correction reaches 0.6610% of the cold mass gap by 99.99% of the vacuum turning time. The correction changes very little over the last quarter of that interval. The velocity is much more sensitive near a zero. At the last sample, the first velocity correction exceeds the entire remaining inward vacuum speed. Only the fixed-background characteristic calculation reaches this sample with matter. The fully coupled run loses quantitative control earlier.

[Flux protocol](inner_horizon_flux_protocol.md) · [Nonlinear protocol](inner_horizon_nonlinear_protocol.md) · [Initial data and causal collar](subgap_neck.md) · [Weak-loading archive](../experiments/closure_ledger/runs/20261127_inner_horizon_flux/probe.json) · [Coupled archive](../experiments/closure_ledger/runs/20261127_inner_horizon_flux/coupled.json)

## The two comparison problems

On the time-oriented cover the causal collar established in #373 prevents a light signal emitted from the neck from winding around the spatial circle. An open interior gives the same evolution on a common domain of dependence if the incoming characteristic data match there. No open-exterior decay exponent was imposed here. Both scalar directions evolve from the actual #373 shell.

The clock τᵥ=1.561023291 belongs to the analytic electrovac geodesic, which arrives at the inner bifurcation sphere there. Positive outgoing vacuum velocity occurs in the analytic extension. For matter, a finite mass at a given proper time does not specify an extension. Even when an extension exists the initial globally hyperbolic evolution can leave its continuation undetermined. Examples of this distinction for characteristic Einstein–Maxwell–scalar data are available in the literature; their four-dimensional results are not used to decide this five-dimensional experiment. [Costa et al., Part III](https://arxiv.org/abs/1406.7261)

## Characteristic evolution at leading order in the shell stress

Use u=t−x, v=t+x, D±=∂ₜ±∂ₓ and K(z)=δ sinh(z). In static and trapped regions invert the signed tortoise function H(r)=sgn(r−r₊) exp(2κ₊r*) at H=−K(u)K(v). The lapse squared is (F/H)K′(u)K′(v)/κ₊². Separate logarithms of the distances to each horizon keep F/H regular at the outer bifurcation and resolve initial light feet that round to the cosmological radius in ordinary floating point. Each numerical diamond takes both parents from the prescribed initial data or from earlier valid diamonds.

Write P=r^(3/2)D₊φ and M=r^(3/2)D₋φ. On exact electrovac their linear propagation is

$$
D_-P=-\frac{3D_+r}{2r}M,\qquad D_+M=-\frac{3D_-r}{2r}P .
$$

A trapezoidal step solves the two endpoint equations simultaneously. The leading charge-renormalized mass gain follows independently along each null ray and the initially normal central geodesic. For example,

$$
 D_+\delta\mu=\frac{\kappa_+ K(v)}{3K'(v)}P^2,\qquad
 D_-\delta\mu=\frac{\kappa_+ K(u)}{3K'(u)}M^2.
$$

The null integrals start with their own initial shell mass; the central integral starts at zero. Trajectory variations are evaluated at fixed geodesic proper time. We also integrated the first variation of the exact timelike radial equation as a supplementary velocity diagnostic,

$$
\delta\ddot r =
\left(\tfrac12+3\mu_0/r^4-10Q^2/r^6\right)\delta r
-\delta\mu/r^3-r(p^2+f^2)/6 .
$$

Here p and f are the vacuum-frame derivatives of the transported scalar. The geometric mass reconstructed from δr and δṙ checks the characteristic integral to the same order. Zero initial radius and velocity variations implement the fixed outer neck in #373. Loading is proportional to scalar amplitude squared.

| Fraction of τᵥ | Standing mass gain | Initially one-way mass gain | Vacuum dr/dτ | Standing velocity correction |
|---:|---:|---:|---:|---:|
| 0.50 | 2.03555×10⁻⁸ | 2.02533×10⁻⁸ | −0.00318997 | −8.035×10⁻⁶ |
| 0.75 | 2.67729×10⁻⁸ | 2.67171×10⁻⁸ | −0.00227174 | −6.762×10⁻⁶ |
| 0.99 | 2.67732×10⁻⁸ | 2.67173×10⁻⁸ | −1.01213×10⁻⁴ | +8.027×10⁻⁷ |
| 0.9999 | 2.67732×10⁻⁸ | 2.67173×10⁻⁸ | −1.01230×10⁻⁶ | +1.138×10⁻⁶ |

The initial load is 2.02515×10⁻⁸. At half time #373's fully coupled standing stress integral is 2.03841×10⁻⁸. At the final sample the last two standing mass estimates differ by 0.0020%; changing the null-label scale by a factor of two changes the 0.004-step result by less than 0.00002%. Across sampled nonzero gains on the finest grid the relative discrepancy between the proper-time integral and the two null integrals or the geometric variation is below 0.002%. These checks address discretization and gauge dependence, not higher-order scalar backreaction.

The final radius correction is −5.2130×10⁻⁶. As a scale comparison, the shift of the *static* inner root at fixed Q for this mass gain is −δμ/(2κ₋r₋²), also about −5.21×10⁻⁶. The last vacuum radius is only 7.9×10⁻¹¹ above its own inner root. Thus evaluating a loaded turn by crossing the original r₋, or by adding the first velocity correction to a nearly zero vacuum velocity, would promote a perturbative estimate to a nonlinear event. A small fractional mass correction alone cannot justify that inference.

## Fully coupled standing shell

The second solver evolves r, the conformal factor σ, μ, and both scalar characteristics, keeping differences from exact vacuum coefficients to avoid cancellation. It integrates the affine expansions N±=exp(−2σ)D±r along their respective rays using D±N±=−exp(−2σ)P±²/(3r²). The scalar transport uses the evolved radial gradients. The lapse source is exp(2σ)[1/2+3μ/r⁴−10Q²/r⁶]−PM/(2r³). The initial longitudinal metric and its derivative come from the constraint-solved loaded mass profile, including its shifted cosmological radius. Vacuum is preserved by construction. The nonzero-load checks compare with #373 and independently reconstruct both mass and velocity.

The coupled archive uses independently sampled proper clocks. It reports two null expansions, the transported mass, the mass reconstructed from the metric, the central stress integral and the timelike integral of radial acceleration. The first invalid characteristic stencil terminates the observer's domain. At τ=0.780512, ε=0.005, step 0.003, the evolved mass gain is 2.03841×10⁻⁸ and the geometric gain is 2.03842×10⁻⁸. At τ=0.936614 the respective values are 2.31225×10⁻⁸ and 2.25824×10⁻⁸. At deeper steps in the same run the geometric mass drifts by more than the scalar increment. Those terminal values are retained to expose the loss of accuracy. At step 0.0015 the same 60%-clock mass discrepancy falls to 0.67%, with dr/dτ=−0.00305140 and its independent acceleration integral −0.00305146. This finest run reaches τ=1.0693 before the stencil guard; its terminal geometric mass error exceeds the scalar increment by two orders of magnitude.

All accepted inward null rays in the initial trapped diamond obey the same focusing sign as in #373. A geodesic turn in the loaded spacetime, an asymptotic divergence of its mass and the regularity of a Cauchy-horizon extension all remain unresolved. In particular the near-horizon table is a finite-interval leading-stress result. The experiment supplies no boundary data on the far side of the horizon.

## Replay

~~~sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.gravity.inner_horizon_flux
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.gravity.inner_horizon_coupled
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_inner_horizon_flux.py'
~~~
