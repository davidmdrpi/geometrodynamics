# A conditional classical derivation of the positive closure measure

The stated frame-restoring Hamiltonian and canonical preparation imply the
positive phase-coarea measure in the low-temperature limit. A different
restoring energy with the same closure circle gives a different measure.
Thus the closure locus is geometric, but its probability measure also
depends on a physical response and preparation. This advances the unresolved
measure derivation without closing the program-wide probability fork.

Predictions and falsifiers were committed locally in
[closure_equilibrium_prereg.md](closure_equilibrium_prereg.md), commit
76ed50e, before implementation. The original base was
92a915bfaaabd02564accd7467b81cedb1ee8c16, after PR #281. PR #282 landed
during development; its [conditioning-variable correction](conditioning_variable.md)
is incorporated unchanged. This derivation supplies the subsequent physical
response model. The preregistration's publication note distinguishes the
local freeze from the later GitHub upload.
The published preregistration is commit
[d83d46a](https://github.com/davidmdrpi/geometrodynamics/commit/d83d46aab8cfe87b0a0adc7f9674401be455d74f).
Implementation: [closure_equilibrium.py](../geometrodynamics/bulk/closure_equilibrium.py).
Independent checks: [test_closure_equilibrium.py](../tests/test_closure_equilibrium.py).
Reproducible evidence:
[probe.md](../experiments/closure_ledger/runs/20260905_closure_equilibrium_probe/probe.md)
and [probe.json](../experiments/closure_ledger/runs/20260905_closure_equilibrium_probe/probe.json).

## 1. The missing specification in zero-set conditioning

For the singlet-sign geodesic triangle, put

$$
 u=s_Aa,\qquad w=-s_Bb,\qquad
 N=x\cdot(u\times w),\qquad D=1+u\cdot w+x\cdot(u+w).
$$

The reduced holonomy is $G=\cos\theta+x\sin\theta$, with
$\theta=\operatorname{atan2}(N,D)$. The full itinerary contributes an
additional minus sign from $J^2=-1$; it does not change the classical
returned frame $R=\operatorname{Ad}_G$. The closure set is the circle
$\Gamma=\{N=0\}$, except $x=-u,-w$, where a geodesic leg is undefined.
All results below concern fixed non-collinear settings.

The circle has zero area under the uniform sphere prior. Merely saying
"condition on $N=0$" does not specify a unique limiting prescription.
On its regular arcs, two legitimate choices give

$$
 \begin{aligned}
 |\theta\bmod\pi|<\epsilon
 &\quad\Longrightarrow\quad d\mu_\Gamma\propto
       \frac{|D|}{|u\times w|}\,d\sigma,\\
 |N|<\epsilon
 &\quad\Longrightarrow\quad d\mu_\Gamma\propto
       \frac{1}{|u\times w|}\,d\sigma.
 \end{aligned}
$$

The first window means distance to the nearest multiple of $\pi$.
Using $\Omega=2\theta$, as round 5 does, divides the coarea density by
two; that common factor cancels from normalized probabilities. Round 5's
implementation used the **phase** window. Its numerical results survive;
the stronger assertion that the zero set alone forces its measure does
not. This document supplies a conditional classical response model for
that prescription.

## 2. A specified classical preparation

Compare an oriented orthonormal frame with its transported return. A
rotation through $2\theta$ has

$$
 \|R-I\|_F^2=6-2\operatorname{Tr}R=8\sin^2\theta.
$$

Choose an isotropic quadratic mismatch energy and a round rotor on $S^2$:

$$
 V_{\rm frame}(x)=\frac{K}{16}\|R(x)-I\|_F^2
                 =\frac K2\sin^2\theta(x),\qquad
 H_s(x,p)=\frac{|p|_{\rm round}^2}{2I}+V_{{\rm frame},s}(x).
$$

Here $K>0$ and $I>0$. This is an additional classical Hamiltonian on
the triangle-family configuration space. It is not derived from BAM's
field action or supplied as a local detector implementation. The energy
depends jointly on both analyzer boundary settings.

In local sphere coordinates $q^i$ with metric $h_{ij}$, canonical
Liouville measure is $d^2q\,d^2p$. Integrating
$\exp[-h^{ij}p_ip_j/(2Ik_BT)]$ over the two momenta gives
$2\pi I k_BT\sqrt{\det h}\,d^2q$. Consequently the position prior is
sphere area, provided the assumed rotor metric is round. The momentum
factor cancels between sectors only when their inertias are identical.
Canonical equilibrium is a specified preparation, not a consequence of
Hamiltonian evolution alone.

Writing $\beta=K/(k_BT)$, define the dimensionless positional partition
function and history-sector probabilities by

$$
 Z_s(\beta)=\frac1{4\pi}\int_{S^2}
       e^{-\beta\sin^2\theta_s/2}\,d\Omega,\qquad
 P_s(\beta)=\frac{\pi_s Z_s(\beta)}{\sum_t\pi_t Z_t(\beta)}.
$$

The positive coefficients $\pi_s$ are additional preparation inputs.
Equal coefficients form the baseline. Calling the resulting sector
probabilities actual detector event frequencies remains an unproved
readout assumption.

## 3. Low-temperature theorem and the punctures

For a smooth residual $F$ with a simple zero on a regular closure arc
and smooth positive stiffness $a$, let $V=KaF^2/2$.
In signed unit-normal distance $n$,
$F=|\nabla_{S^2}F|n+O(n^2)$. Gaussian integration gives

$$
 \sqrt{\frac{\beta}{2\pi}}
 \int e^{-\beta aF^2/2}f\,d\Omega
 \longrightarrow
 \int_\Gamma\frac{f}{\sqrt a\,|\nabla_{S^2}F|}\,d\sigma
$$

on compact regular arcs, for a smooth test function $f$. On closure,
$\nabla N=u\times w$ and
$|\nabla\sin\theta|=|u\times w|/|D|$. Thus the frame potential implies

$$
 \boxed{\displaystyle
 \sqrt{\frac{\beta}{2\pi}} Z_s(\beta)\longrightarrow
 \frac1{4\pi}\int_\Gamma\frac{|D_s|}{|u\times w|}\,d\sigma .}
$$

The punctures require an argument beyond regular coarea. Set
$q=|u\times w|>0$. Near either puncture, the local coordinates
$(d,n)=(D/q,N/q)$ have a nonsingular bounded area Jacobian. Indeed the
sphere-tangent gradients of $D$ and $N$ there are perpendicular and
both have norm $q$. In these coordinates
$\sin^2\theta=n^2/(d^2+n^2)$. Polar coordinates give an angular integral
$\int_0^{2\pi}e^{-\beta\sin^2\psi/2}d\psi=O(\beta^{-1/2})$;
a coordinate disc of radius $\epsilon$ therefore contributes
$O(\epsilon^2\beta^{-1/2})$ to the partition integral, uniformly for
$\beta\ge1$. After the displayed scaling it is $O(\epsilon^2)$.
Away from the circle and these discs, the potential has a positive lower
bound on compact sets and the contribution decays exponentially.

The limiting density also vanishes linearly at each puncture. Excising
arclength radius $\epsilon$ around **each** puncture removes total
unnormalized circle mass $2\epsilon^2+O(\epsilon^4)$, since
$|dD/d\sigma|=q$ there. Taking the temperature limit and then shrinking
the discs proves the whole-sphere result. No value of the holonomy at the
punctures is assigned. Constants in this argument can depend on the
fixed angle; no uniform collinear limit is claimed.

## 4. Closed form and explicit countermodels

Let $c=u\cdot w$, $\delta=\arccos c$. Parametrize the circle from the
bisector $u+w$, so
$D=1+c+\sqrt{2(1+c)}\cos\phi$. Splitting the integral at its two zeros
and integrating $|D|$ gives

$$
 M(c)=\int_\Gamma\frac{|D|}{\sqrt{1-c^2}}\,d\sigma
     =4+\frac{2(1+c)(\pi-\delta)}{\sqrt{1-c^2}}.
$$

At analyzer separation $\gamma$, $c_s=-s_As_B\cos\gamma$.
With equal priors,

$$
 E_\infty(\gamma)=
 \frac{M(-\cos\gamma)-M(\cos\gamma)}
      {M(-\cos\gamma)+M(\cos\gamma)}.
$$

This is round 5's positive correlation with the singlet partner sign.
The standard coplanar CHSH value is **2.1422831632**. This calculation
does not assert a global CHSH maximum.

| Physical response $V/K$ | Limiting density on $\Gamma$ | Consequence with equal priors |
|---|---|---|
| $\sin^2\theta/2$ | $\lvert D\rvert/q$ | Positive phase-coarea law |
| $N^2/2$ | $1/q$ | $E=0$ at every temperature |
| $(gN)^2/2$ | $1/(gq)$ | Generally changes the sector law |
| $(gN)^2/(2g^2)$ | $1/q$ | Same energy and law as $N^2/2$ |
| $(\sin^2\theta+\lambda\sin^4\theta)/2,\ \lambda\ge0$ | $\lvert D\rvert/q$ | Same leading limit; different finite-temperature corrections |

For the normal-residual model, $N_s^2$ is pointwise identical in all
four sectors. Its partition integral has the independent Gaussian oracle

$$
 Z_N=\int_0^1e^{-az^2}\,dz
 =\frac{\sqrt\pi\,\operatorname{erf}\sqrt a}{2\sqrt a},
 \qquad a=\beta(1-c^2)/2,
$$

with value 1 at $a=0$. It has the same regular closure set as the frame
model, but a different transverse stiffness.

For any positive smooth $g(x)$, changing the residual coordinate
$F\mapsto gF$ and stiffness $a\mapsto a/g^2$ preserves the energy
pointwise, hence the measure at every temperature. On the zero set
$|\nabla(gF)|=g|\nabla F|$, so the limiting density is also invariant.
Holding $a$ fixed changes the physical energy. It is not a coordinate
ambiguity of a fully specified apparatus.

The explicit control $g=1+\epsilon x\cdot(u+w)$, $|\epsilon|<1/2$,
is positive on the whole sphere and has

$$
 M_{gN}(c)=\frac{2\pi}{q\sqrt{1-2\epsilon^2(1+c)}}.
$$

At $\gamma=1,\epsilon=1/4$, its limiting correlation is
$-0.0386506750$, versus zero for the covariant-stiffness control.
The quartic frame deformation leaves the leading transverse Hessian
unchanged; its nonnegativity also bounds puncture contributions by those
of the original frame model.

## 5. Numerical evidence and reproduction

The partition integrator covers the whole sphere. In bisector coordinates
$N=qz$ and $D=1+c+\sqrt{2(1+c)}\sqrt{1-z^2}\cos\phi$.
All tested energies are even in $z$. The substitution
$z=\sinh t/\sqrt{\max(\beta,1)}$, with its full range and Jacobian,
resolves the thermal width without imposing a closure cutoff.
Gauss-Legendre integration uses 128 normal nodes and 1024 midpoint
azimuth nodes. Each coordinate is doubled separately for convergence.
The integrator does not call the limiting formula.

| $\gamma$ | $E_{\beta=4096}$ | $E_\infty$ | Largest joint-cell error | Largest scaled-mass relative error |
|---|---|---|---|---|
| $\pi/4$ | -0.5348083838 | -0.5355707908 | $1.91\times10^{-4}$ | $2.10\times10^{-3}$ |
| 1 | -0.3980052277 | -0.3984966504 | $1.23\times10^{-4}$ | $1.19\times10^{-3}$ |
| $3\pi/4$ | 0.5348083838 | 0.5355707908 | $1.91\times10^{-4}$ | $2.10\times10^{-3}$ |

These meet the frozen $2\times10^{-3}$ joint and 1% mass criteria.
Separate coordinate refinements change a joint cell by at most
$1.67\times10^{-7}$, below the frozen $10^{-4}$ threshold.
At $\gamma=1,\beta=16$, independent uniform-sphere Monte Carlo with
120,000 samples gives $E=-0.309181587\pm0.000900255$ (one standard
error), versus quadrature $-0.308680932$. It computes all four sectors
from unreduced three-vectors and estimates ratio uncertainty; half
marginals are not enforced in that control. Monte Carlo here integrates
the equilibrium law; it does not demonstrate equilibration.

The full report retains every preregistered temperature, frame-transport
identity checks, Gaussian and variable-stiffness controls, puncture
excision, the quartic deformation, and sensitivity to sector priors.
The CLI exits nonzero if any numerical criterion fails.

    python -m experiments.closure_ledger.closure_equilibrium_probe
    python -m experiments.closure_ledger.closure_equilibrium_probe --output-dir /tmp/closure-equilibrium
    python -m pytest -q tests/test_closure_equilibrium.py

## 6. What has and has not been derived

| Ingredient | Status after this calculation |
|---|---|
| Holonomy and regular closure circle | Inherited geodesic-frame geometry |
| Uniform sphere position prior | Derived from assumed round-rotor kinetics and canonical Liouville preparation |
| Positive phase-coarea law | Derived in the low-temperature limit of the stated frame energy |
| Residual-coordinate covariance | Exact when physical stiffness transforms with the residual |
| Choice of frame energy and canonical preparation | Additional classical assumptions; BAM origin and equilibration open |
| Relative sector coefficients | Still inputs; paired half marginals do not fix the like/unlike ratio |
| Local apparatus, composition, and event frequencies | No dynamical derivation supplied |
| Source-readout causality | Remains open; the conditioned model uses both future settings |
| Born rule or exhaustion of classical alternatives | Neither established |

All canonical history weights are nonnegative. This model therefore
supplies no history-by-history signed cancellation by the closure
holonomy. Positivity alone does not rule out another positive model
reproducing the integrated quantum law, nor coherent classical field
responses. The result identifies a concrete conditional measure and a
physical way to vary it; it is not a theorem selecting an aggregation rule
for the full BAM theory.
