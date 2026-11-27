# Sub-gap loading of an electrovac neck

At fixed charge and fixed initial black-hole neck, a localized scalar load of 0.5% or 5% of the cold mass gap produces a measurable increase in the mass seen by the falling observer. At half the exact vacuum turning time the standing and initially one-way 5% loads both leave the observer moving inward. The calculation resolves this response well before it can determine the fate of the inner horizon.

Two geometric checks change the interpretation of this experiment. The null focusing equation excludes a turn while the two past light rays originate on the inward-facing static flanks. Separately, an expanding electrovac region stays causally protected for the entire future of its cosmological observer. The first light circuit used in #372 therefore has no counterpart for these data in the globally hyperbolic development.

[Registered protocol](subgap_neck_protocol.md) · [Full archive](../experiments/closure_ledger/runs/20261126_subgap_neck/probe.json) · [Parent conventions](cold_loop_mass.md)

## Fixed neck, independently variable scalar energy

Use the untilted bifurcation slice of #371 with squared reference horizon radii (0.396, 0.4, 1.204). Units are those of #371. The invariant charge is 0.308798964 and the reference inner surface gravity is 0.00648483. At fixed charge the mass excess above the cold double root is 0.00000405030831. On every slice the geodesic starts at the same outer bifurcation sphere, with zero areal velocity. Its future begins with inward acceleration. Both Maxwell stress and geometry are evolved on the time-oriented cover.

Write s=r², a=0.396, b=0.4. Put a smooth shell on each exterior flank from s=b+0.3(b−a) to b+3(b−a). The reference mass inside is μ₀=0.558392. Across each shell the mass rises by ε times the cold gap, according to the regularized beta integral with parameters (9,9). This definition of loading compares two adjacent exact electrovac regions on the initial slice. On the periodic cover the signed increments around the circle sum to zero.

With U=0 initially, use w=∂ₗr and the unit normal and spatial scalar derivatives p and f. The Hamiltonian and momentum constraints reduce to

$$
w^2=1-\mu(r)/r^2+Q^2/r^4-r^2/2,\qquad
p^2+f^2=3r^{-3}\,d\mu/dr,\qquad h w= rpf/3.
$$

Set h=0 outside the shell. The integration constant at the neck fixes the interior electrovac mass for every load. The standing data have p=0, equal initial null energies and zero integrated scalar momentum. The globally chiral control has p=f, with the same pointwise energy and a single initially occupied null direction. Its signed longitudinal expansion fails time-reversing half-shift parity, so it is a cover control. A third preparation reverses chirality on the antipodal half, satisfying the half-shift with an odd scalar. For all three, the integrated shift charge vanishes to numerical precision. The control energies are matched locally and in the proper-volume integral.

## Proper-time comparison

We follow the initially normal timelike geodesic, including its displacement for chiral data. The exact zero-Killing-energy orbit reaches r₋=0.629285309 at τᵥ=1.561023291 and has positive velocity in the analytic extension. The table uses the common event τ=τᵥ/2 and the finest focusing grid (16,384 points). The last column integrates scalar stress along the geodesic. For a moving observer the energy flux term matters. With w measured in its spatial frame, the rate is r³[−U(p²+f²)+2wpf]/3.

| Initial data | ε | dr/dτ | Geometric mass gain | Integrated stress gain |
|---|---:|---:|---:|---:|
| Electrovac | 0 | −0.00318989 | −2.06×10⁻¹⁰ | 0 |
| Standing | 0.005 | −0.00319792 | 2.018×10⁻⁸ | 2.038×10⁻⁸ |
| Standing | 0.05 | −0.00327030 | 2.062×10⁻⁷ | 2.064×10⁻⁷ |
| Globally chiral, cover | 0.05 | −0.00326984 | 2.024×10⁻⁷ | 2.025×10⁻⁷ |
| Antipodally paired chiral | 0.05 | −0.00334970 | 4.153×10⁻⁷ | 4.156×10⁻⁷ |
| Standing | 1 | −0.00482292 | 5.249×10⁻⁶ | 5.250×10⁻⁶ |

At the 5% standing event the geometric-to-stress discrepancy falls from 3.31×10⁻⁹ to 8.34×10⁻¹⁰ to 2.10×10⁻¹⁰ over the focusing grids. The exact vacuum velocity is −0.00318997. Inside the cropped chart, edge cells supply no incoming data. The finest 5% standing grid reaches τ=1.011 (about 65% of τᵥ). The last cells already have a vacuum mass drift larger than the 5% load. Matched earlier events carry the quantitative result.

The two chiralities mix dynamically. For D±=∂ₜ±∂ₓ in conformal coordinates the exact wave equation gives

$$
D_+D_-\phi=-\frac{3}{2r}\big[(D_+r)(D_-\phi)+(D_-r)(D_+\phi)\big].
$$

On the chiral shell initially only D₊φ is occupied and D₋r is nonzero. The source generating the other component reaches magnitude 0.101 in the periodic coordinate normalization. At half time both geodesic-frame null energy components are nonzero in the cover control. The matched comparison tests initial directionality; a permanently single-stream propagating control would require a different matter model. The runs do not establish cross-flow as the unique mechanism of destruction.

## What the two evolutions can decide

We ran all six preparations on three grids in each of two formulations. The original periodic mixed-derivative scheme gives an apparent vacuum turn around τ=0.8534 on 8,192 points. The exact orbit is still contracting there. This false positive fails both the vacuum benchmark and a pointwise geometric check. The second scheme uses a cropped optical chart, integrates the affine-weighted null focusing equations and tests the mixed radial equation through the residual. For a canonical neutral scalar, radial Maxwell flux saturates the radial null energy condition and each ray obeys

$$
D_\pm(e^{-2\sigma}D_\pm r)=-\tfrac13 r e^{-2\sigma}(D_\pm\phi)^2\le0.
$$

On the initial left flank the right-moving ray points toward decreasing r. On the right flank the left-moving ray also points toward decreasing r. At every regular intersection with these initial feet, both null derivatives are strictly negative. Every future timelike radial velocity there is a positive combination of them, hence dr/dτ<0. This applies to the standing, chiral and vacuum data. It also explains why checking the proper clock and the past rays of an apparent turn is essential. The sign proof covers this causal diamond; it does not determine an extension beyond its boundary.

The stiffness can be quantified in the exact vacuum. Take the real tortoise coordinate r* with dr*/dr=1/F and r*(∞)=0. The two light rays reaching the reflection orbit at r have initial static feet with the same real r*. At 90% and 99% of τᵥ the fractions (r_c−r_foot)/(r_c−r₊) are respectively about 10⁻⁵⁹ and 10⁻¹³⁵. The optical chart resolves small coordinate distances by logarithmic parametrization. Evolved metric gradients and the finite valid stencil still limit the coupled run. The analytic orbit supplies an outgoing branch in exact electrovac only. Neither a positive critical amplitude nor its collapse to zero has been measured for the scalar family. A finite-duration bracket in one family would also leave stability under general perturbations open.

## Why a light ray cannot lap this bridge

This supplementary check uses only the untouched *outer* electrovac mass μ₀+Δμ and uniqueness in the domain of dependence. On each side the initial data are exactly electrovac for r>0.6418723. At the cosmological bifurcation sphere the future normal geodesic has zero Killing energy and expands toward infinity. Its proper time is unbounded since dr/dτ=√(−F) grows like r/√2. In regular cosmological Kruskal coordinates the two past radial light rays from the observer at radius R hit the static initial flanks at equal tortoise coordinate r*(r_foot)=r*(R). This follows by matching the magnitude of each Kruskal null label between the zero-Killing-time initial slice and the reflection orbit. With the outer electrovac roots rᵢ,

$$
r_*(r)=\sum_i\frac{1}{F'(r_i)}
\log\left|\frac{r-r_i}{r+r_i}\right|.
$$

The limiting initial foot solves r*(r_foot)=0. For ε=0.05 it is 0.85383362; across all four loads it is at least 0.85383324. Every finite-time foot lies even closer to the cosmological sphere. Thus the observer's entire past on the initial slice lies in the exterior electrovac interval, with a radius margin above the shell exceeding 0.21196. The spherical Cauchy data there determine the observer's full future neighborhood. A signal from the initial neck never reaches this worldline. On the periodic orbit space a ray returning after one winding would have to cross it. This removes the first-circuit clock for these data on the globally hyperbolic cover. It makes no identification across a future Cauchy horizon. The scalar outcome there remains unresolved.

## Reproduction

~~~sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.gravity.subgap_neck --workers 6
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p 'test_subgap_neck.py'
~~~

The archive retains common-clock samples, grid-limited endpoints, apparent turns from both schemes, constraint diagnostics, exact-orbit optical feet, and the protected cosmological intervals.
