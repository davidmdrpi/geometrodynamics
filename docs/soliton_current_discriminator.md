# Which current follows a localized soliton in the winding gradient?

The strong drive pulls the two winding components apart. Their total amplitude maximum moves ahead of both current tracers while a tail of the minority charge extends away from the peak. At the end of the drive, the peak position is **0.894819**. The laboratory energy and scalar number characteristics from the same starting point end at **0.787070** and **0.787742**. The charge median follows the number characteristic. The amplitude peak follows the deformation of the soliton. These are separately measurable features of the same nonlinear solution. The primary comparison and its resolution criterion were [registered before the run](soliton_current_protocol.md).

## A field model with an independent core

The reference vortex in [#199](guidance_law_from_5d.md) has quantized phase winding. Its current comparison used an ambient ring average in a nonrelativistic transverse condensate. Here the modeled surrogate is a charge-stabilized relativistic Q-ball. Its nonzero peak permits evaluating both currents directly at the core. This model extends the linear field in [#359](https://github.com/davidmdrpi/geometrodynamics/pull/359) by a specified nonlinear binding action. The two complex components represent the opposite-force winding channels. Their common interaction is a modeling choice. An unrestricted nonlinear scalar on the full fiber could generate higher harmonics.

In units with mass and action scale one, take the Lagrangian density
$$
 \mathcal L=\tfrac12\left(\sum_{k=\pm1}
 [|D_t u_k|^2/c^2-|\partial_xu_k|^2]-U(s)\right),\qquad
 D_t=\partial_t+ikA,\quad s=\sum_k |u_k|^2,\quad
 U(s)=c^2s-2s^2+\tfrac12s^3 .
$$
Thus $D_t^2u_k=c^2[u_{k,xx}-U'(s)u_k]$. With zero drive the exact rest solution is $u_k=\sqrt{w_k} f(x)e^{-iWt}$, where
$$
 W=c\sqrt{c^2-1},\qquad f(x)^2=\frac{2}{2+\sqrt{2}\cosh(2x)} .
$$
Integration of the static equation gives $(f')^2=f^2-2f^4+f^6/2$. This fixes the initial spatial profile and time derivative independently of any guidance law. Cubic–quintic complex scalar theories admit this family of analytic Q-balls. [Continuum reference](https://doi.org/10.1016/j.nuclphysb.2022.115964)

Initialize $w_+=0.8535130337631147$ from #359's renormalized $k=\pm1$ state at $\theta=\pi/8$, and $w_-=1-w_+$. Over $0\le t\le 2$, set $A=-x\sin^2(\pi t/2)$ near the core, with a taper to zero for $|x|\ge28$ in a periodic box $[-32,32)$. At maximum drive the two force components are +1 and −1. Evolution starts after the mixer.

## What was measured

1. Locate the maximum of $s(x,t)$ using a spline around the highest grid point. Track its curvature, height and secondary peaks. A second locator solves $s_x=0$ using a Fourier interpolant. At a nondegenerate peak its instantaneous velocity is independently available as $-s_{xt}/s_{xx}$.
2. Launch two tracers from that maximum at $t=0$. Advance one with $S/E$ and the other with $J/n$, each evaluated at its own position on the evolving fields. There is no tracer force on the soliton.
3. Independently integrate final and intermediate charge densities to find the spatial median. Initially the peak and the two medians of energy and number coincide by reflection at zero offset. Equivariance predicts that the number characteristic will continue to follow the number median as long as the scalar number density is positive in its neighborhood and boundary flux is negligible.

With $b_k=D_tu_k$, the unsigned phase density and flux are $n=\sum n_k$, $n_k=-\operatorname{Im}(u_k^*b_k)/c^2$, and $J=\sum\operatorname{Im}(u_k^*u_{k,x})$. Here “unsigned” refers to adding the two channels without a winding-sign factor. Local KG phase density can have either sign in general. The laboratory energy and flux for the specified time generator are
$$
 E=\tfrac12\!\left[\sum_k (|b_k|^2/c^2+|u_{k,x}|^2)+U(s)\right]+A\sum_k k n_k,
 \qquad S=-\sum_k\operatorname{Re}[u_{k,x}^*(b_k-ikAu_k)] .
$$
The balances are $\partial_t n+\partial_xJ=0$ and $\partial_t E+\partial_xS=(\partial_t A)\sum_k k n_k$. The potential energy of binding enters $E$. The code checks positive denominators at every sampled tracer position.

Spatial derivatives use Fourier differentiation. RK4 advances the field in a rotating frame with the exact rest frequency, retaining both relativistic frequency branches. The characteristics use a second-order predictor–corrector. The primary grid has 1,024 sites, field step 0.001 and readout step 0.01. Refinement doubles the spatial resolution and halves both steps. The matched nonrelativistic limit evolves $i\partial_t u_k=[-\partial_x^2/2+(1-4s+1.5s^2)/2+kA]u_k$ after removal of the carrier.

## Outcomes and interpretation

| Mixed preparation | Final core | Final energy tracer | Final number tracer | Maximum energy–number separation |
|---|---:|---:|---:|---:|
| $c=8$, refined | 0.894819 | 0.787070 | 0.787742 | 0.00067211 |
| $c=12$ | 0.897135 | 0.788676 | 0.788968 | 0.00029166 |
| $c=16$ | 0.897956 | 0.789237 | 0.789399 | 0.00016265 |
| Nonrelativistic limit | 0.899016 | 0.789956 | 0.789956 | 0 |

The largest primary-to-fine trajectory change was $1.28\times10^{-5}$, giving the frozen locking tolerance $6.40\times10^{-5}$. The two predictions separate by over ten times that tolerance. The maximum refined core residuals are 0.107748 for energy and 0.107076 for number. Both exceed the tolerance by more than three orders of magnitude. The refined core retains 84.6% of its initial peak height and stays single-peaked throughout the drive. At the end, the independent peak-velocity measurement gives **1.043271**, while the instantaneous energy and number velocities evaluated **at the peak itself** are **0.901542** and **0.901732**.

The energy–number separation scales approximately as $1/c^2$ across the three finite carriers. The much larger core–current separation survives the nonrelativistic limit. Its origin can be seen directly from the continuity equation in that limit. Write $J=sv$. At a maximum, $s_x=0$ and $\dot X=J_{xx}/s_{xx}=v+s v_{xx}/s_{xx}$. A velocity field with spatial curvature can move a density peak relative to the matter passing through it. The two oppositely driven components change the shape of the bound profile. Charge conservation constrains a fixed cumulative charge, while a peak is determined by a spatial derivative.

### Controls

- On the refined mixed run the charge median follows the number tracer within $1.55\times10^{-6}$. Spline and Fourier peak locations agree within $1.62\times10^{-6}$. Maximum relative charge drift is below $2\times10^{-14}$ and the energy–work residual is $1.56\times10^{-7}$ of initial energy.
- With no drive the resting peak stays fixed. A free Lorentz boost at speed 0.5 transports the core and both tracers to $x=1$ at $t=2$. Reversing the gradient reverses the mixed trajectories.
- A constant offset $A\mapsto A+0.1$ preserves the physical gradient and the core and number paths to roundoff at the endpoint. It changes the laboratory energy endpoint from 0.787070 to 0.787611 on the primary grid. This tests the choice of laboratory time generator.
- The pure-winding control isolates translation under a single force. At $c=8$ the core ends at 0.995875 and the two primary tracer endpoints are 0.996876 and 0.996908. The pure core retains over 99% of its initial peak height. In the exploratory pure-winding series, the maximum core–number residual falls from 0.001086 at $c=8$ (0.001075 refined) to 0.000494 at $c=12$ and 0.000284 at $c=16$. The nonrelativistic single-force solution translates essentially rigidly and its Fourier-located peak reaches $x=1$ within $2\times10^{-13}$. This control separates the small finite-carrier discrepancy of one winding from the much larger deformation of the mixed core.

The result is a resolved **neither** for the registered amplitude-peak definition and the specified soliton. A Q-ball is a Noether-charge-stabilized field configuration. A vortex zero, a spatially localized quantized topological charge, and a geometric wormhole waist have different definitions of position. This experiment supplies a counterexample to inferring a peak's trajectory from conservation alone. For BAM the geometric waist or topological defect must be located independently in its own evolution. Exact continuity of a KG or Dirac number current establishes equivariance for its own density under its own velocity. Predicting an exact $\cos^2\theta$ curve also requires the prepared ensemble and analyzer's branch weights at finite carrier speed.

## Where the component charge goes

The review asked for a fourth position measurement and a force sweep. These are exploratory follow-ups to the registered peak test. At $t=2$ the refined $c=8$ total number centroid, median and amplitude peak are **0.705531, 0.787742 and 0.894818**, respectively (Fourier peak locator). The individual number centroids are +0.976903 and −0.875633. Within each winding, we normalize its number density separately and integrate the fraction outside $|x-X_{\rm peak}|\le3$ using a cubic antiderivative. Initially each winding has 0.561% outside. Under full drive this reaches 11.586% for the minority and 0.570% for the majority. A single maximum thus coexists with partial spatial separation. A two-unit observation cannot determine how much charge escapes permanently.

| Relative drive, $c=8$ | Maximum peak–number-median gap through the drive | Minority outside at $t=2$ |
|---:|---:|---:|
| 0 | 0 | 0.561% |
| 0.125 | 0.000144 | 0.623% |
| 0.25 | 0.001736 | 0.827% |
| 0.5 | 0.015500 | 1.919% |
| 1 | 0.107075 | 11.586% |

Weakening the force reduces both the tail and the peak–median gap. At one-eighth drive the gap is still resolved: refining the grid and time step changes it by less than $2\times10^{-9}$. The sweep gives no sharp binding threshold. The test covers this smooth peak, this preparation and this interaction.

For the NLS limit, Ehrenfest's identity gives $\ddot{\langle x\rangle}=-\sum_k\int |u_k|^2\partial_x(kA)\,dx/\int s\,dx$. The internal nonlinear force integrates to a boundary term proportional to $[sU'(s)-U(s)]$. With an exactly linear potential and a resting ensemble, $\langle x\rangle(2)=F_0(w_+-w_-)$ because the integrated displacement of the unit pulse equals one. The numerical potential is tapered in the far tails. At full force, the refined NLS centroid is **0.707026067526**, within $2\times10^{-14}$ of the uniform-gradient value. The two component centroids are +0.979015 and −0.877733. The Fourier peak is 0.899023 and the median is 0.789956. The center of charge and the charge median obey different laws.

## Does a winding integer fix a velocity?

In a transverse plane write a complex winding field as $u=U+iV$. At a simple zero $X(t)$, differentiating $u(X(t),t)=0$ gives
$$
 D(U,V)\,\dot X=-\partial_t(U,V).
$$
The sign of the spatial Jacobian determinant fixes the local index. The field equation fixes the time derivative. This distinction already matters in the nonrelativistic mixed-channel model.

For example, take initial fields $u_+=(x+iy)e^{\alpha x-r^2/2}$ and $u_-=e^{-r^2/2}$ on the plane, with nonzero real $\alpha$. Under $i\partial_tu_k=-\Delta u_k/2+V_k u_k$, any smooth multiplicative winding potential or density interaction drops out of the equation at the zero of $u_+$. At the origin $D(U,V)=I$ and $\Delta u_+=2\alpha$, giving $\dot X=(0,-\alpha)$. The total density is positive there and $\sum_k\operatorname{Im}(u_k^*\nabla u_k)=0$. Thus even a simple **component** vortex of index +1 can move relative to the total number current. This local example uses a filled core, a density gradient, and a specified winding basis. It addresses the claimed implication from topology alone. #199's ambient-flow comparison for its prepared single-condensate vortex, and the dynamics of a common geometric waist, are separate questions.

[All centroid measurements, force levels and refinements](../experiments/closure_ledger/runs/20261123_soliton_discriminator/force_sweep.json)

~~~sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.soliton_current
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.soliton_force_sweep
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p test_soliton_current.py
~~~

[Trajectory samples, full-drive residuals and control runs](../experiments/closure_ledger/runs/20261123_soliton_discriminator/probe.json)
