# Classical traversable-waist research record

**Decision: close the following four tested routes to a smooth static traversable waist with causal, well-posed evolution and the stated matter content.** A background satisfying algebraic stress equations is insufficient if its physical metric perturbations have no hyperbolic evolution. All throat statements below use positive lapse, a regular radius and a strict local minimum. A prime denotes proper-radial differentiation unless the affine null parameter is specified.

| Action and scope | Necessary test | Result |
|---|---|---|
| Einstein gravity, round flaring waist, minimally coupled matter satisfying the null energy condition (NEC) | Radial null focusing | $R_{kk}=-m r''/r<0$ at a waist with an $m$-sphere orbit, whereas Einstein's equation gives $\kappa T_{kk}\ge0$. [Source audit](source_audit_prereg.md) |
| Five-dimensional constant-coupling EGB, round waist, NEC matter with no metric second derivatives in its stress and no principal metric–matter mixing for the transverse shear | Radial NEC and angular shear propagation | NEC requires $1+4\alpha/r^2\le0$. This same factor multiplies the time derivative of an angularly propagating physical shear. The strictly negative branch is elliptic in the static clock; equality degenerates. |
| Metric $f(R)$, $F=f_R>0$ smooth, regular positive spherical radius on a periodically repeating static radial null traverse, vacuum or NEC matter | Integrated affine null equation | Periodicity forces constant radius along the traverse, excluding flare-out. |
| Five-dimensional Einstein–scalar–Gauss–Bonnet, arbitrary smooth coupling, canonical massless radial scalar, static round waist without extra sources | Temporal and radial constraints plus physical radial shear | The time–radius symbol is negative definite. The lapse asymmetry required by the constraint cannot restore a hyperbolic radial mode. [Coupled calculation](dilatonic_waist_characteristics.md) |

## Why the negative constant coupling also fails the propagation test

Let $Y=-r''/r<0$, $Z=1/r^2>0$. The radial null equation is $\kappa T_{kk}=3Y(1+4\alpha Z)$ in an orthonormal frame. This recovers the threshold in the [constant-coupling audit](gauss_bonnet_prereg.md), including its erratum about exact Lovelock theory. To test the remaining branch, propagate along angular direction 2 with polarization $h_{13}$, where 1 is radial and 3 angular. At the waist this is transverse, is gauge-invariant at principal order, and sources no other principal metric component. Direct curvature variation gives

$$
2\,\delta E_{13}/h_{13}=-(1+4\alpha Z)\omega^2+k_2^2 .
$$

In the time coefficient the two complementary directions span an angular plane of curvature Z. In the spatial coefficient they span a time–angle plane of zero curvature at $r'=0$. This gives an imaginary angular frequency wherever NEC permits a strictly negative coefficient. At equality the spatial second derivative persists with no second time derivative. These are local obstructions to the requested causal evolution on a static slice. The code checks the aligned shear against all metric components at positive, zero and negative time coefficient. The negative-coupling background branch and its matter conditions are discussed by [Maeda and Nozawa](https://arxiv.org/abs/0803.1704).

## Why periodic metric f(R) fails

Let a dot denote affine differentiation on the radial ray and assume $r,F,\dot r,\dot F$ repeat after one spatial traverse. Time can stay unwrapped. The null field equation and spherical warped-product curvature give $-m\ddot r/r=\ddot F/F+\kappa T_{kk}/F$. Integration gives

$$
-m\oint(\dot r/r)^2\,d\lambda=\oint(\dot F/F)^2\,d\lambda+
\kappa\oint T_{kk}/F\,d\lambda .
$$

The signs force both squares to vanish. A single local minimum without the periodic ray is insufficient for this proof. The regular positive-F hypothesis also matters for the broader scalar–tensor literature. [Bronnikov, Skvortsova and Starobinsky](https://arxiv.org/abs/1005.3262)

The record leaves changing topology, time-dependent or non-spherical waists, nonminimal additional sources and other actions as separate cases. In particular, the prescribed network metric used in packet experiments supplies no solution to these support tests.
