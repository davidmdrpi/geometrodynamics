# Action completion and the first constrained kinetic test

**Result: `STATIC_ACTION_MATCH_FROZEN_MODE_REJECTED`.** A specified classical
scalar action reproduces #330's static bulk and shell stresses and its affine
surface EOS. Its nominated material sector has positive bare quadratic kinetic
coefficients. However, #330's nonzero radial oscillation with the bulk metric
and support fields frozen changes the **total conserved scalar shift charge**:
Q'_r=-r0^3 s'_stat/q is nonzero throughout the WEC/C>0 window. The trial
oscillation therefore leaves the background's fixed-charge sector. Its
linearized density-junction and scalar-current failures express this same
obstruction. It is not a physical perturbation to which a canonical energy
sign can be assigned in this action.

This conclusion assumes both **a continuous scalar whose surface value is the
trace of the frozen bulk field** and **a single-invariant barotropic surface
action P_Sigma(Y)**. It does not exclude all possible material completions.

**The full constrained radial kinetic sign is still uncomputed.** This test
identifies a precise obstruction to computing it from the old oscillator alone;
it neither demonstrates a gravitational ghost nor rejects coupled solutions of
the new action. The positive restricted radial window of #330 remains valid
within its stated prescribed-exterior model. This completion adds a material
postulate and is not a derivation of BAM matter from geometry.

The [protocol](egb_shell_kinetic_protocol.md) states the gate before the numerical
checks, after the analytic prediction. The [result archive](../experiments/closure_ledger/runs/20261013_egb_shell_kinetic/result.json)
records failed physical admission separately from successful numerical checks.
It contains no claimed physical radial norm: that field is null.

## A complete local variational candidate

Use the same two-mouth direct identification, one surface Sigma and two incident
bulk faces. The retained bulk is the S4 interval chi_angle in [a,pi-a], with
common unshifted time. There is no finite neck or Alice lift. The new neutral
scalar is denoted chi; it is continuous across the identification. In signature
-++++, set chi=sqrt(kappa5)*phi and define

\[
X=-\tfrac12 g^{AB}\partial_A\chi\partial_B\chi,\qquad
Y=-\tfrac12 h^{ab}D_a\chi D_b\chi.
\]

The candidate action is

\[
S=\frac1{2\kappa_5}\int_{M}\sqrt{-g}
 [\mathcal R+\alpha\mathcal L_{GB}+2P_b(X)]\,d^5x
+S_{\rm Myers}
+\frac1{\kappa_5}\int_\Sigma\sqrt{-h}\,P_\Sigma(Y)\,d^4x.
\]

With normals directed into the retained bulk, the two incident gravitational
boundary terms are

\[
S_{\rm Myers}=-\frac1{\kappa_5}\sum_{\pm}\int_\Sigma\sqrt{-h}
 [K_\pm+2\alpha(J_\pm-2\widehat G_{ab}K_\pm^{ab})]\,d^4x.
\]

This is the Davis convention already used in #330. Usual temporal endpoint
data/terms must be treated consistently in a Hamiltonian reduction; none are
discarded to assert a radial energy sign here. Varying the bulk metric gives
the EGB equations with the material stress below, and varying the induced metric
gives the same Davis junction. The scalar is varied continuously across Sigma,
so the surface current equation is supplied by the same action. No externally
prescribed reservoir or corrective kicks is added.

The nominated bulk law is **P_b=X-V**, with

\[
\chi_0=qt,\qquad q^2=\frac{3C}{R^2},\qquad X_0=\frac{q^2}{2},\qquad
V=\frac12\left(\frac9{R^2}+\frac{12\alpha}{R^4}\right),
\quad C=1+4\alpha/R^2>0.
\]

The scalar equation is box(chi)=0, satisfied because q and the ultrastatic
volume element are time independent. Its barred stress is
rho=X0+V=6/R^2+12 alpha/R^4 and p=X0-V=-3/R^2, exactly the required exterior.
This q and constant potential are choices tied to the prescribed background,
not a dynamically selected bulk radius.
Thus the previously unspecified support is replaced by a definite rolling
scalar, filling the entire retained bulk and carrying conserved shift charge.
It is not a field confined to the focus or shell.

Let s0,p0 be the existing barred static shell stresses and W0=s0+p0>0. For
**0<eta<=1**, take

\[
n=\frac{1+\eta}{2\eta},\qquad
P_\Sigma(Y)=p_0+\frac{W_0}{2n}
 \left[\left(\frac{Y}{X_0}\right)^n-1\right].
\]

Then s=2Y P_Y-P and p=P obey **p-p0=eta(s-s0)** exactly. At equilibrium
Y=X0 and both stresses match the junction. This is a local timelike-gradient
model (Y>0); eta=0 is not covered by this finite-power action. No claim about
large excursions through Y=0, microscopic realization, or global Hamiltonian
boundedness is made. Fine tuning of the intercept remains explicit.
At the chosen point the constant term in P_Sigma is +10.14708988/R; by itself
it contributes **negative surface vacuum energy**. The kinetic contribution
brings the total static density and pressure to the previously tested DEC
values. Positive bare kinetic coefficients and total background DEC do not
mean that every separately written material contribution has positive energy.

## What the actual material action fixes

On a fixed background, expand chi=qt+pi. The local quadratic material action is

\[
\frac1{2\kappa_5}\int\sqrt{-g}
 [Z_b\dot\pi^2-P_{b,X}|\nabla\pi|^2]
+\frac1{2\kappa_5}\int_\Sigma\sqrt{-h}
 [Z_\Sigma\dot\pi^2-P_{\Sigma,Y}|D\pi|^2],
\]

where Z_b=P_b,X+2X0 P_b,XX and Z_Sigma=P_Y+2X0 P_YY. For the nominated
completion,

\[
Z_b=P_{b,X}=1,\qquad
P_{\Sigma,Y}=\frac{W_0}{q^2}>0,\qquad
Z_\Sigma=\frac{W_0}{\eta q^2}>0,\qquad c_\Sigma^2=\eta.
\]

At (R,a,alpha,eta)=(1,.3,-.06,.1), q=1.50996688705, X0=1.14, V=4.14,
n=5.5, P_Sigma,Y=13.43072354 and Z_Sigma=134.30723541. The bulk and surface
bare kinetic signs are positive, while the old D0=-8.91048348 remains negative.
This does not establish a positive **coupled radial** kinetic norm.

Since chi=qt is stationary under a time translation combined with an internal
shift, we also differentiate the fixed-geometry local energy H-qJ, with
H=P_X dot(chi)^2-P and J=P_X dot(chi). Its quadratic time coefficient is the
same Z, after removing the linear background-charge variation. This local
material check is not the constraint-reduced gravitational Hamiltonian.
Independent finite differences of the Lagrangians and H-qJ recover all these
coefficients at two step sizes, with maximum normalized error 6.63e-10.

### Why static support does not determine kinetic sign

Consider the declared countercontrol family

\[
P_{b,\lambda}(X)=X-V+\frac{\lambda}{2X_0}(X-X_0)^2.
\]

All members have the same P_b(X0), P_b,X(X0)=1, background scalar solution,
bulk stress and frozen-bulk shell flux. The surface action is unchanged. Yet
their bulk time-kinetic coefficient is **Z_b=1+2 lambda**.

| lambda | Static stress and frozen flux | Z_b | Bulk scalar principal behavior |
|---|---|---:|---|
| -1 | Identical | -1 | Rejected: negative bare time-kinetic coefficient; elliptic scalar principal form |
| 0 | Identical | 1 | Canonical candidate; local scalar speed squared 1 |
| 1 | Identical | 3 | Positive local material coefficients; scalar speed squared 1/3 |

The negative control is not offered as a second admissible theory. It proves
that the background/flux data do not even fix the bulk scalar quadratic sector.
Nor are these asserted to be theories sharing an actual physical radial mode:
the old clamped mode fails the next constraint check. Additional derivative or
curvature couplings would be different completions and require new analysis.

## The scalar boundary equation supplies the exchange law

For the nominated canonical bulk action, variation with respect to chi gives

\[
D_a(P_{\Sigma,Y}D^a\chi)+\sum_\pm n_\pm^A\partial_A\chi=0.
\]

For a homogeneous shell, this is

\[
\frac1{r^3}\frac{d}{d\tau}(r^3P_{\Sigma,Y}\dot\chi)
=\sum_\pm n_\pm(\chi).
\]

If the bulk trial fields are held at chi=qt and the original metric, then
dot(chi)=q gamma, gamma=Delta/sqrt(f), and sum n(chi)=2qv/sqrt(f). The surface
Ward identity consequently demands

\[
\dot s+3\frac vr(s+p)
=\dot\chi\sum n(\chi)=2q^2\gamma\frac v{\sqrt f}
=\frac{6Cv\Delta}{R^2f}.
\]

Thus the old source has a definite covariant, reversible realization **if the
scalar boundary equation is satisfied**. Equal source formulas do not make
every old trajectory a solution of the new matter equations.

## The frozen radial oscillator fails a first-order constraint

For those frozen bulk fields, the induced invariant is exactly

\[
Y=X_0\left(1+\frac{v^2}{f}\right).
\]

Around r=r0, v=0 this gives **delta Y=0**, hence delta s_material=0 at first
order for every smooth P_Sigma(Y) in this class. In contrast, the gravitational
density junction requires

\[
\delta s_{\rm junction}=s'_{\rm stat}(r_0)\,\delta r,
\qquad s'_{\rm stat}(r_0)=\frac{3D_0}{r_0^2}\ne0.
\]

Since D0<0 throughout the WEC/C>0 window, the linearized junction with both
metric and scalar frozen forces **delta r=0**. The nonzero reduced oscillator
is not tangent to the solution space of this material completion. This is a
constraint failure before a kinetic sign can be assigned, not a proof of a
negative-energy physical mode and not a family-wide wormhole no-go.

Numerically s'_stat=-306.089308562. At turning-point displacements epsilon r0,
with epsilon=1e-3,1e-4,1e-5,1e-6 and v=0, the density residual divided by
(1+|s0|)epsilon approaches **4.83949**, not zero. Its ratio to the analytic
linear residual approaches one (0.99999776 at the smallest amplitude).
The failure persists in the infinitesimal limit; shrinking amplitude does not
pass the constraint. Units in this numerical normalization are R=1.

The scalar boundary residual has the independent linear form

\[
\mathcal J= q\left(\frac{3P_{\Sigma,Y}}{r_0}
-\frac2{\cos a}\right)v+O(\epsilon^2)
=-\frac{s'_{\rm stat}}q\,v+O(\epsilon^2).
\]

Its coefficient is 202.712596671. For the retained closed bulk the total shift
charge per unit S3 area, evaluated on the trial fields, is

\[
Q=qR^4\int_{\arcsin(r/R)}^{\pi-\arcsin(r/R)}\sin^3 u\,du
+r^3P_{\Sigma,Y}(Y)q\gamma.
\]

There is one shell and two excised caps. Direct differentiation gives
dot(Q)=r^3 J and at equilibrium Q'_r=-r0^3 s'_stat/q=5.23169337259. The local
boundary obstruction and closed-system charge failure agree independently.
For an oscillatory perturbation, conservation requires delta Q to be constant;
the frozen trial instead gives delta Q=Q'_r delta r(t). Changing the initial
charge sector does not cure its time dependence. A constant shift of chi also
cannot cure it, since it changes neither Y nor Q.

The two assumptions above are essential to the stated proof: trace continuity
gives delta Y=0 when the bulk is frozen, and P_Sigma(Y) then gives delta s=0.
An independent surface field or additional material invariant changes this
argument and requires a new action and all its current junctions. The review's
suggestion that a tuned two-field surface action might reopen a frozen tangent
is an **untested possibility**, not an established counterexample. Counting
adjustable Hessian entries does not establish a constraint-satisfying mode or
a positive physical kinetic norm.

Evaluating, without changing, #330's archived 20-period trajectory gives maximum
density mismatch **0.09087813** and scalar-current residual **0.08365584**.
Both current signs occur. These arrays remain solutions of the old reduced ODE;
they are not relabeled as an evolution of this new action.

## What a physical radial-mode calculation must now include

Already on a fixed metric, matching the density at first order would require
the boundary scalar trace rate to satisfy

\[
\frac{\delta\dot\chi_\Sigma}{q}
=\frac{\eta s'_{\rm stat}}{s_0+p_0}\,\delta r
=-0.999571589175\,\frac{\delta r}{R}
\quad\hbox{at the chosen point}.
\]

The near-unit numerical coefficient is not a selection rule. The exact static
Ward identity gives

\[
R\frac{\eta s'_{\rm stat}}{s_0+p_0}
=-\frac{3\eta R}{r_0}
+\frac{6\eta C}{R\cos a\,(s_0+p_0)}.
\]

At the chosen point these terms are -1.015159008547 and +0.015587419372.
At fixed geometry and static stresses the sum scales linearly with eta:
eta=.05 gives -0.499785794588 and eta=.2 gives -1.999143178350.
Nothing in this calculation selects -1.

This is a necessary boundary condition under the stated fixed-metric trial,
not a solved bulk field. For the canonical completion, bulk perturbations have
delta rho=q delta(dot chi) and radial momentum proportional to q partial_r(delta
chi). Keeping the bulk metric unperturbed would set the corresponding EGB
left-hand sides to zero, so an oscillatory scalar response cannot simply be
added while retaining the metric. The lapse, shift and spatial geometry must
respond **at the same perturbative order as the shell**. Thus the
prescribed-exterior approximation fails at first order in this completion;
neither #330's frequency nor its stability window has been reproduced by this
action. A full treatment also perturbs the normal scalar derivative;
the old prescribed flux and frequency cannot be carried over automatically.

The action now specifies the material response that #330 lacked. A coupled
calculation must retain metric/scalar perturbations on the entire compact
retained bulk, solve the bulk lapse/shift and shell density constraints, impose
the perturbed Davis and scalar-current junctions, account for total shift charge
and its constant-shift zero mode, and evaluate the resulting physical quadratic
form including gravitational boundary terms. Only then can its radial kinetic
eigenvalue and mode frequency be reported. Positive material Hessians alone
do not pass this test. Boundary-normal canonical momentum in Lovelock junction
formulas is not automatically the physical time-evolution kinetic norm.

This PR supplies that missing material action and proves why the old frozen
mode cannot answer the requested kinetic question. It does **not** claim to
have completed the coupled constraint reduction. Smooth finite thickness,
full well-posedness, nonspherical stability, focal activation and a quantum
event remain separate uncomputed requirements.

## Sources and reproduction

- [Davis, hep-th/0208205](https://arxiv.org/abs/hep-th/0208205): EGB boundary action
  and metric junction conventions. The new scalar boundary equation is derived
  here by variation, with the same normals.
- [Garriga and Mukhanov, hep-th/9904176](https://arxiv.org/abs/hep-th/9904176),
  equations (4) and (10): scalar-fluid stress and acoustic coefficient. We use
  the opposite metric signature consistently and only classical material
  identities; no quantization rule or four-dimensional cosmological perturbation
  equation is imported into this five-dimensional shell.
- [Guilleminot, Merino and Olea, 2204.00767](https://arxiv.org/abs/2204.00767):
  action and boundary-normal canonical momentum in Lovelock shell matching.

```bash
python -m experiments.closure_ledger.egb_shell_kinetic_probe
python -m pytest tests/test_egb_shell_action.py -q
```

The archive binds the new source files and original trajectory digest. Numerical
identity errors and physical admission failures are recorded separately. No
original producer, protocol or simulation archive is modified.

Validation: 88 focused tests pass with NumPy 2.3.5 and 2.2.6 (SciPy 1.17,
Python 3.12.14), including nine new action tests and the previous 79 shell,
focus/audit and negative-EGB checks. Repository-wide CI is separate.

### Independent review and follow-up

The [review of a21a9f2](https://github.com/davidmdrpi/geometrodynamics/pull/331#issuecomment-6139078212)
reports 88 passing tests with NumPy 2.4.6/SciPy 1.17.1 and a producer rerun
matching every archived field except environment. It also reports symbolic
checks of the density-slope and EOS identities and green Python 3.10/3.12 CI.
These are attributed reviewer checks; they are not additional archived runs
of this PR's producer. This follow-up leads with charge conservation, states
the two material assumptions, and expands the trace-response coefficient.
The original protocol, producer, source-bound manifest and trajectory are
unchanged.
