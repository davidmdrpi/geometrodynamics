# A radial tensor obstruction at an asymmetric five-dimensional waist

For the action and static spherical waist below, a regular asymmetric lapse can satisfy the radial constraint. On any such background satisfying the temporal constraint as well, the physical radial shear has **negative squared characteristic speed**. This excludes a causal, strongly hyperbolic evolution of this waist in the static clock, independently of the smooth coupling function. It is a local obstruction, before a global extension or a search for long-time instabilities. [Protocol](pointer_and_waist_protocol.md) · [Covariant contraction check](../experiments/closure_ledger/antipodal_network/dilatonic_waist.py)

## Action and necessary background equations

Use signature −++++ and proper radial distance s. In units absorbing the gravitational normalization take

$$
 S=\int\sqrt{-g}\,[R-\tfrac12(\nabla\phi)^2+f(\phi)(R_{abcd}R^{abcd}-4R_{ab}R^{ab}+R^2)]\,d^5x,\qquad
 ds^2=-N(s)^2dt^2+ds^2+r(s)^2d\Omega_3^2 .
$$

Assume a smooth radial scalar, unwrapped time, 0<N<∞, r>0, r′=0 and r″>0 at the waist. The field content excludes a scalar potential, cosmological constant, material shell and other matter. Work in a local orthonormal frame, with time 0, radius 1 and angles 2,3,4. The nonzero sectional curvatures there can be written as X=−N″/N in the time–radius plane, Y=−r″/r in each radius–angle plane and Z=1/r² in each angle–angle plane. Time–angle curvature vanishes at the waist. Write p=φ′, and denote the mixed Hessian entries of the coupling by

$$
 D_0=(N'/N) f_\phi p,\qquad D_1=f_\phi\phi''+f_{\phi\phi}p^2,\qquad D_i=0\quad(i=2,3,4).
$$

With $P_{abcd}=R_{abcd}-(g_{ac}R_{bd}-g_{ad}R_{bc}-g_{bc}R_{ad}+g_{bd}R_{ac})+\tfrac12 R(g_{ac}g_{bd}-g_{ad}g_{bc})$, the metric equation is $G_{ab}+fH_{ab}+4P_{acbd}\nabla^c\nabla^d f=\tfrac12\partial_a\phi\partial_b\phi-\tfrac14 g_{ab}(\nabla\phi)^2$, with H the Lanczos tensor. The scalar equation is $\Box\phi+f_\phi\mathcal L_{\rm GB}=0$. Direct contraction gives

$$
 E^1{}_1=-3Z+12ZD_0-p^2/4=0,\qquad
 E^0{}_0=-3(Y+Z)-12fYZ+12ZD_1+p^2/4=0 .
$$

At a symmetric lapse the first equation has no real solution. Asymmetry allows it, but requires $D_0=1/4+p^2/(48Z)>0$. Smoothness thus forces $p\ne0$ and $f_\phi\ne0$. The constant-coupling result in the [original Gauss–Bonnet study](gauss_bonnet_prereg.md) and this variable-coupling test use different source assumptions.

## A physical factor of the coupled principal symbol

Take a covector with orthonormal components $(\omega,k,0,0,0)$ and perturb $h_{34}=h_{43}$, with $\delta\phi=0$. The angular off-diagonal shear is transverse and trace-free. A pure gauge perturbation $\xi_a v_b+\xi_b v_a$ has zero 34 entry. Radial propagation preserves the angular rotation group, so this shear cannot mix with the scalar. At principal order the only independent nonzero metric component is

$$
 2\,\delta E_{34}/h_{34}=-C_t\omega^2+C_s k^2,\qquad
 C_t=1+4fY-4D_1,\qquad C_s=1-4D_0 .
$$

One way to obtain these coefficients is to vary the orthonormal curvature with $\delta R_{i3j4}=\xi_i\xi_j h_{34}/2$ for $i,j=0,1$, together with curvature symmetries. The Einstein variation gives $(k^2-\omega^2)/2$. The Lanczos term contributes $-4fY\omega^2/2$; the coupling Hessian contributes $(4D_1\omega^2-4D_0 k^2)/2$. The scalar curvature-invariant variation $\delta\mathcal L_{\rm GB}$ is zero for this polarization. Terms from varying the connection in the Hessians contain at most one metric derivative. This evaluation retains the metric-to-scalar equation and all metric components sourced by this polarization.

Eliminating the Hessian entries by the two background constraints gives

$$
 C_t=\frac{p^2-12Y}{12Z}>0,\qquad C_s=-\frac{p^2}{12Z}<0,\qquad
 \boxed{\frac{\omega^2}{k^2}=-\frac{p^2}{p^2+12r''/r}<0 .}
$$

For real nonzero k the characteristic frequencies are imaginary, with growth rate proportional to |k| in the frozen-coefficient equations. Lower-order terms and gauge choices cannot remove this physical mode. A single such polarization suffices to reject the requested static-time hyperbolicity. The calculation does not classify all characteristic cones of less symmetric backgrounds.

## Checks and scope

The probe evaluates the full five-index-range Ricci, Lanczos and derivative-coupling contractions, differentiating them in a shear direction by complex step. At 32 off-shell jets varying curvature, coupling, scalar gradient, Hessian and covector, it checks both constraints, every shear-sourced metric component, and the scalar GB variation against the analytic expressions. The maximum absolute discrepancy was $2.5\times10^{-14}$. The flat Einstein control gives luminal propagation.

An example satisfying **all** metric and scalar equations at one point has $r=N=1$, $p=f_\phi=1$, $f=0$, $N'/N=13/48$, $X=2/3$, $Y=-1/2$, $\phi''=-205/48$, and $f_{\phi\phi}=35/8$. A smooth local Taylor polynomial realizes those derivatives of f. The largest background residual was $2.3\times10^{-16}$ and the directly calculated shear speed squared was −1/7. The local data alone do not assert existence of a globally regular solution.

This obstruction uses five spacetime dimensions, spherical symmetry, stationarity, the canonical kinetic sign and the specified action. Time-dependent waists, a potential or extra matter require new equations. Four-dimensional dilatonic wormhole backgrounds and their radial stability analyses concern another problem. Strong hyperbolicity results at weak coupling likewise require their own hypotheses. [Kanti, Kleihaus and Kunz](https://doi.org/10.1103/PhysRevLett.107.271101) · [Cuyubamba, Konoplya and Zhidenko](https://doi.org/10.1103/PhysRevD.98.044040) · [Kovács and Reall](https://doi.org/10.1103/PhysRevD.101.124003)
