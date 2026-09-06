# Source readability: a parity correction and a conditional classical pointer

The causality gate remains open, but two proposed protections are insufficient.
Antipodal symmetry does **not** make the full law of an odd readout blind to
future settings. Nor does the existence of two source boundary conditions
alone prohibit a classical record. An explicit Hamiltonian pointer can leave
both source boundaries unchanged under an additional ideal preparation.

Neither result constructs an operational BAM signal. The missing connection
is still a map from the informative triangle coordinate to a local BAM field
observable, together with an admissible apparatus, its ensemble, and a theory
of independently selectable future settings. This follow-up to the
[review of #283](https://github.com/davidmdrpi/geometrodynamics/pull/283#pullrequestreview-5122111334)
also identifies a conditional local parent for its frame energy and checks
the numerical overlaps with rounds 5–8.

## 1. Odd means are blind; odd readout laws need not be

If `mu(-dx|a,b) = mu(dx|a,b)` and `F(-x) = -F(x)`, then the distribution of
`F` is symmetric and, when integrable, `E[F|a,b] = 0`. Symmetry does not make
that distribution the same at different settings. In particular its variance
can differ. Measuring an odd quantity repeatedly also measures its even
moments; the apparatus need not directly couple to `F²`.

The sign statement also needs a qualification:

\[
 P(F>0)=P(F<0)=\frac{1-P(F=0)}2.
\]

The nondegenerate continuous projections tested in round 7 have no atom at
zero and both sign probabilities are `1/2`. A projection normal to a sharp
closure circle vanishes identically; the original claim of `1/2` for *every*
axis was too broad. Historical frozen outputs are retained. Current code,
tests, and the status audit distinguish means/signs from full laws.

Here is the example frozen in the
[public preregistration](source_readout_prereg.md). Fix

\[
 a=e_z,\quad b_0=e_x,\quad b_1=e_y,\quad
 m=(1,2,0)/\sqrt5,\quad F(x)=m\cdot x.
\]

For either future choice write `x = cos(phi) a + sin(phi) b_j`. Summing the
four equal-prior sectors of the positive phase-coarea law gives

\[
 p(\phi)=\frac{1+|\sin\phi|+|\cos\phi|}{2\pi+8}.
\]

Quarter-turn invariance gives `E[sin²(phi)] = 1/2`. The two readout amplitudes
are `1/sqrt(5)` and `2/sqrt(5)`, respectively, hence

\[
 E[F|b_0]=E[F|b_1]=0,\qquad
 \operatorname{Var}(F|b_0)=\frac1{10},\quad
 \operatorname{Var}(F|b_1)=\frac4{10}.
\]

This is an exact counterexample to distributional blindness of odd
observables. It uses the same axis and the same response at both settings.

It does not imply that **every** nondegenerate odd projection separates
these settings. For a general fixed unit axis `m`, the two sharp variances
are `(m_x²+m_z²)/2` and `(m_y²+m_z²)/2`. Their difference is
`(m_y²-m_x²)/2`, which vanishes on a nonempty set of axes. In particular
`F=x_z` is continuous, odd and nondegenerate, with variance `1/2` in both
ensembles. Its **entire law** agrees: the quarter-turn about `a=e_z`
mapping `b_0` to `b_1` preserves `F` and carries one source measure to the
other. This symmetry also holds at finite beta and after the same readout
noise. Thus parity supplies no general guarantee of unreadability, but
sign tests are not the only records that can fail to distinguish a given
pair. This qualification and its two regression cases were added after
the review at `fe5016d`; they are not new preregistered results.

For finite resolution let `Y = F + Z`, with independent Gaussian noise
`Z ~ Normal(0,0.15²)`, and use the fixed event `B = {|Y|>0.6}`. Its kernel is

\[
 K(B|x)=\Phi\!\left(\frac{-0.6-F(x)}{0.15}\right)
       +\Phi\!\left(\frac{F(x)-0.6}{0.15}\right).
\]

| Ensemble | Future choice | Variance of F | P(F>0) | P(abs(Y)>0.6) |
|---|---|---:|---:|---:|
| sharp positive phase coarea | `b_0` | 0.1000000000 | 0.5 | 0.0528578345 |
| sharp positive phase coarea | `b_1` | 0.4000000000 | 0.5 | 0.5136058265 |
| #283 Gibbs, beta = 64 | `b_0` | 0.1335456636 | 0.5 | 0.1231586349 |
| #283 Gibbs, beta = 64 | `b_1` | 0.3904155247 | 0.5 | 0.4897143360 |

The finite-temperature computation integrates the whole sphere, summing
unnormalized sector densities before normalization. Independently doubling
the normal and azimuthal grids changes either variance by less than `1e-4`.
It does not assume singular supports. Thus the conditional distinction
survives both finite temperature and finite readout resolution.

There is also a noiseless antiderivative. For amplitude `A` and threshold
`t`, set `v=t/A`. When `0<v<1`,

\[
 P(|F|>t)=\frac{2[\pi/2-\arcsin v+\sqrt{1-v^2}+1-v]}{\pi+4}.
\]

At `t=0.6` the two sharp probabilities are `0` and `0.5338532795`. Midpoint
quadrature of this discontinuous indicator differs by `2.54e-5`; the smooth
Gaussian-tail calculation above avoids that indicator error.

These are **conditional record laws**, not an experiment establishing a
channel. A setting-independent constant kernel gives the same record law
for both measures. Existence of that kernel is no universal non-readability
theorem; existence of the informative Gaussian kernel is no BAM apparatus.

## 2. A canonical pointer with the backreaction included

On a regular source chart with canonical coordinates `(x^i,p_i)`, append
pointer coordinates `(Q,P)` and the classical interaction

\[
 H_{\rm tot}=H_S(x,p,t)+\frac{P^2}{2M}+h(t)P F(x).
\]

For the interaction-only impulse with `integral h dt = g`, Hamilton's
equations integrate exactly to

\[
 x'=x,\quad p'_i=p_i-gP\partial_iF,\quad
 Q'=Q+gF(x),\quad P'=P.
\]

This map is symplectic. The tests check its Jacobian and compare it with an
independent Hamilton ODE built from complex-step derivatives of `P F`,
including nonzero pointer momenta. The reciprocal force has not been dropped.

For a finite pulse, `P=0` is an invariant preparation because `H_tot` is
independent of `Q`. On it the source follows exactly its old Hamiltonian
history for **any** `H_S`, while

\[
 Q(t_1)-Q(t_0)=\int_{t_0}^{t_1}h(t)F(x(t))\,dt.
\]

Both old source endpoint values survive. If `H_S` is autonomous, the explicit
work term `partial_t H_tot = h'(t) P F` vanishes on this preparation too.
The final pointer position is a record, not prescribed final Dirichlet data.

The construction assumes an exactly zero pointer momentum and a known
initial position. Classical mechanics allows that ideal preparation; it is
additional physics here. Finite momentum produces recoil. The impulse has
`|delta p| <= |gP| |dF|`; a finite pulse also changes the source trajectory,
so its original endpoint constraints cannot simply be assumed preserved.
No claim of a robust no-disturbance apparatus with nonzero momentum width is
made. For an independent zero-mean initial momentum with standard deviation
`sigma_P`, the impulse's conditional recoil covariance is
`g² sigma_P² dF dF^T`; its root-mean-square size at fixed source position
is `|g| sigma_P |dF|`. Zero mean recoil is not zero disturbance.
A model requiring the pointer itself to return to its initial position
would impose `integral hF dt = 0`: for `h>=0` and `F=q_osc²`, that prohibits a
nonzero record. That is an *extra apparatus boundary requirement*, not a
consequence of the two source endpoints alone. We have not derived BAM's
complete apparatus boundary problem.

### A finite pulse on an actual repository source Hamiltonian

The conservative source probe
[`hamiltonian_source_eigenhistory_probe.py`](../experiments/closure_ledger/hamiltonian_source_eigenhistory_probe.py)
already implements, in `_h_red` and `_rhs_red`,

\[
 H_S=\tfrac12(P_f^2+\omega_r^2Q_f^2)
     +\tfrac12(p_s^2+\omega_0^2q_s^2)
     +\tfrac\mu4q_s^4+g_{\rm eff}q_sQ_f.
\]

We use its actual equations, `omega_0=3.2`, `mu=0.5`, with fixed
`omega_r=2.738858`, `g_eff=0.3203`, and field turning amplitude `Q_f(0)=0.4`.
Shooting the two momenta to zero at half-period gives a nontrivial periodic
coupled source/field orbit. Add `h(t)P q_s²`, with
`h(t)=(2g/T) sin²(pi t/T)`, `g=0.4`, and pointer mass one.

| Quantity | Zero-P preparation |
|---|---:|
| period T | 2.2997607422 |
| initial source displacement q_s | -0.0461461914 |
| original source closure error | 2.76e-14 |
| source endpoint change after adding pointer | 2.28e-13 |
| maximum source trajectory change | 6.04e-12 |
| pointer displacement | 0.0004258902143 |
| total energy range | 1.00e-11 |

The small nonzero discrepancies are numerical integration errors; source
preservation at `P=0` follows exactly from the equations. With `P=0.001`,
the source endpoint changes by `1.52e-5`, and the time-dependent interaction
can do work. This control does not keep energy constant by fiat.

**The oscillator `q_s` is not the quaternion frame `q`, nor is `q_s²` known
to equal an informative `F(x)`.** This is a source-local coupling in an
extended version of an existing reduced Hamiltonian. **At ideal zero pointer
momentum**, it invalidates an obstruction based solely on “measurement adds
a boundary condition.” Nonzero momentum spread generically restores recoil,
linear in the spread for the impulse's root-mean-square momentum change. It
does not complete the source-readout construction requested by the live
falsifier, because the informative field map and allowed preparation remain
unproved.

## 3. The operational test requires interventions and the instrument ensemble

Let `z` denote the physical state of a complete instrument-plus-source
model `I`. For each independently selectable future setting, the test is

\[
 P_I(Y\in B\mid a,\operatorname{do}(b))
 =\int K_I(B\mid z)\,\mu^I_{a,\operatorname{do}(b)}(dz).
\]

Operational no-signalling to the past requires equality for all early-record
events `B` and all permitted future choices, with the early instrument and
preparation fixed and no selection on a later record. One admissible event
with unequal probabilities would refute it. The expression presumes a
conditional response kernel exists for the full physical model; it does not
grant that kernel or the intervention measure to BAM.

In particular `mu^I` is the **instrument-modified** ensemble. Substituting the
unmeasured conditional posterior `rho(x|a,b)` requires a separate proof that
the physical measurement preserves that ensemble. Likewise conditioning on
`b` is not yet a theory of independently changing `b`. Neither the Gibbs
prescription in #283 nor the pointer toy supplies this missing intervention
structure. A non-readability theorem must constrain the allowed physical
observables, apparatus/preparation, modified measure, or future-setting
control. An absent implementation proves none of these constraints.

### What a finite-momentum follow-up must measure

The review proposes replacing `P=0` by a distribution. That is a useful
next experiment, but merely sampling momenta in `pointer_kick` cannot decide
the operational gate: **for the same incoming source ensemble**, its record
`Q_out=Q_in+gF(x)` is independent of `P` at every spread. The source momenta
change, and subsequent dynamics and boundary conditioning may change which
incoming histories are allowed. Reusing the old posterior would omit the
very effect being investigated.

A finite pulse has both pointer drift `P(t_1-t_0)/M` and a record of the
perturbed source trajectory. A meaningful comparison must specify the same
initial apparatus distribution for both future choices, solve the joint
boundary problem with the interaction present, determine its history weights,
and compare the full early-record laws without selecting on later outputs.
The existing Duffing example contains no analyzer-setting or frame-variable
map, so a spread experiment on that oscillator alone cannot supply the
required comparison. Such a map and a physical preparation/intervention law
are separate gaps from the ideal-pointer limitation.

There is a useful conditional continuity check on any proposed result. Let
`L_{b,sigma}` be the early-record law of that **complete instrument model**,
and suppose its zero-spread event contrast is `Delta_0>0`. If
`TV(L_{b,sigma},L_{b,0}) <= epsilon_b(sigma)`, then

\[
 |L_{b_1,\sigma}(B)-L_{b_0,\sigma}(B)|
 \geq \Delta_0-\epsilon_{b_1}(\sigma)-\epsilon_{b_0}(\sigma).
\]

If both errors tend to zero, some nonzero spreads retain the distinction.
Complete erasure at every nonzero spread would therefore require a singular
limit or failure of these hypotheses. No continuity or nonzero `Delta_0`
for an operational BAM instrument has been established here. Nor would
erasure for one particular apparatus prove that every admissible apparatus
is unreadable. These are analytic checks added in response to the review,
not a completed finite-spread dynamics experiment.

## 4. What the existing field machinery does and does not contain

The inventory below concerns the implementations at merged #283, `b0e372f`.
Scalar charges, radial amplitudes and spatial points on an `S³` must not be
identified with a Spin(3) frame merely because symbols or manifolds coincide.

| Machinery | Actual degrees of freedom and coupling | Consequence for frame energy/readout |
|---|---|---|
| [`tangherlini/dynamics.py`](../geometrodynamics/tangherlini/dynamics.py), [`radial.py`](../geometrodynamics/tangherlini/radial.py), [`operator_audit.py`](../geometrodynamics/tangherlini/operator_audit.py) | Spherical Einstein–scalar fields `phi(v,r), A, delta`; separated scalar radial master operator; radial action `integral p_r dr*` | No triangle-frame orientation in this nonlinear spherical ansatz. Keeping its actual fields fixed leaves no variable on which the proposed frame energy could depend. The radial one-form action does not supply a frame comparator. |
| [`waves/finite_throat.py`](../geometrodynamics/waves/finite_throat.py), [`finite_mouth.py`](../geometrodynamics/waves/finite_mouth.py) | Scalar tube and its DtN response; one monopole channel per mouth, with higher multipoles explicitly omitted | Supplies a positive static gradient energy usable in the conditional extension below. It does not supply transported frame channels or endpoint matching to the analyzer itinerary. |
| [`waves/throat_operator.py`](../geometrodynamics/waves/throat_operator.py), [`two_wave.py`](../geometrodynamics/waves/two_wave.py) | Complex scalar mouth charges with flux `Im(q_j* (Aq)_j)`; scalar ESU waves and their interference | Flux is a real available observable, but no frame-to-field family has been exhibited. Spatial `S³` coordinates are not automatically mouth-frame quaternions. |
| [`transaction/particles.py`](../geometrodynamics/transaction/particles.py), [`cavity.py`](../geometrodynamics/transaction/cavity.py) | Radial `ThroatMode` amplitudes, scalar `q_geom`, damped/driven cavity modes | No orientation-sensitive comparison with `Ad_G`; these `q` symbols do not define the needed frame coordinate. |
| [`transaction/derived_network.py`](../geometrodynamics/transaction/derived_network.py), `cavity.py:closure_check` | Scalar/spinor channel labels, a discrete wrap sign, and a phase-closure test with supplied `phi_spin` | Existing holonomy references are real, but neither a sign nor a phase-tolerance test supplies a continuous frame-restoring field energy or its inertia. |
| [`tangherlini/lepton_spectrum.py`](../geometrodynamics/tangherlini/lepton_spectrum.py) | Instanton-transition **surrogate** with `closure_term = CAVITY_GAMMA (1-cos(phase_per_pass*sep))/2` | A sinusoidal closure penalty already exists. It is an inserted action term in a spectral surrogate, not a field-derived energy on the triangle family. No identification of its phase with `2 theta` is supplied. |
| Conservative source Hamiltonian above | Scalar Duffing source coupled to a scalar field mode | Permits the explicit local oscillator pointer extension in section 2. Does not relate its oscillator state to `x`. |
| [`shells/junction.py`](../geometrodynamics/shells/junction.py) | Spherical Israel shell radius, density, pressure, radial potential and stiffness; vacuum regions with Birkhoff decoupling | Orientation dependence is absent in this restricted sector. Its stiffnesses are not frequencies without a kinetic metric; no rotor inertia can be read off from them. |
| [`bulk/source_audit.py`](../geometrodynamics/bulk/source_audit.py) | Scalar, Maxwell and other local stress functions | Gives matter observables to test once a frame-parametrized field solution exists. Degree in field amplitude is not parity under `x -> -x`. |

There is no implemented pullback of these field Hamiltonians onto the
triangle family. This is a scoped missing derivation, not a theorem forbidding
all nonspherical classical completions. In the spherical/monopole models an
independent internal orientation cannot acquire a nonconstant potential
while all represented fields are held fixed. Adding orientation-carrying
fields changes that premise.

The lepton surrogate is the closest existing algebraic lookalike: formally
setting `phase_per_pass*sep=2 theta` turns its closure penalty into
`CAVITY_GAMMA sin²(theta)`. That substitution would be another assumption,
and converting this surrogate action coefficient into a rotor energy would
require a dynamical normalization as well. It therefore does not ground
#283's potential. The spectral parameters and mass fits are untouched.

## 5. A conditional local parent for the Frobenius energy

The massless scalar tube's zero-frequency DtN matrix is

\[
 N_0=\frac{\mathcal A}{L}
 \begin{pmatrix}1&-1\\-1&1\end{pmatrix}.
\]

Extend this positive gradient energy to an unconstrained real matrix field
`U(s) in R^(3x3)`, with endpoint data `U(0)=I`, `U(L)=R=Ad_G`:

\[
 E[U]=\frac{\mathcal A}{2}\int_0^L\|U'(s)\|_F^2 ds.
\]

Write `U=I+s(R-I)/L+v` with `v(0)=v(L)=0`. Integration of the cross term
gives zero, so

\[
 E[U]=\frac{\mathcal A}{2L}\|R-I\|_F^2
       +\frac{\mathcal A}{2}\int_0^L\|v'(s)\|_F^2ds.
\]

The minimum is therefore exactly #283's potential with `K=8 A/L`.
Independent finite-difference elimination of the nine interior components
agrees to `4.45e-16`; the existing DtN routine at `kL=1e-6` agrees to
`7.43e-13`. Its literal `k=0` evaluation has a removable `0/0`; the control
approaches the static limit without modifying that older routine. At finite
frequency the DtN quadratic form belongs to a Helmholtz variational problem,
not this positive static energy.

This demonstrates **classical realizability after specified extensions**.
The additional channels can be viewed as three transported vectors, but the
interior matrix is not an orthonormal frame. Constraining `U(s)` to `SO(3)`
would be a different nonlinear model, generally with a different minimum.
Actual grounding still needs those channels in the action, a connection and
path identifying `R` with the analyzer holonomy, rigid endpoint matching,
and a justified length/scale. The proper throat interval has not been shown
to be the full analyzer itinerary. Static elimination alone establishes
neither a causal implementation nor equilibration.

Bi-invariance does not uniquely select the exact potential:
`||LR T-LS T||_F=||R-S||_F` for orthogonal `L,T`, and nonlinear functions of
that squared norm are invariant too. Quadratic gradient response makes this
particular chordal energy natural *within the added linear field model*.
Symmetry alone does not force it or exclude a different classical coupling.

## 6. The assumption budget is relocated and stronger, not unchanged

For a quadratic rotor with kinetic metric `G_ij`, canonical momentum
integration gives

\[
 \int d^2p\,e^{-p_iG^{ij}p_j/(2kT)}
 =2\pi kT\sqrt{\det G}.
\]

Round constant inertia `G=I_0 h` on `S²` yields uniform sphere area before
the potential is applied. Thus #283 derives the inherited position prior
**conditionally on** a canonical phase space, that kinetic metric and
canonical equilibrium. Equal sector kinetic prefactors are also needed if
the intended sector priors are to remain equal. Liouville volume by itself
does not select a canonical equilibrium state or prove that it is reached.

These are stronger sufficient premises than a Haar position prior. They
imply it, but it does not imply them. Three explicit controls show why:

| Alternative kinetic choice | Position marginal |
|---|---|
| Rotation-invariant quartic tangent momentum energy `H=|p|^4/4`, in units `kT=1` | Haar, despite non-Gaussian momenta: radial integral `2pi integral r exp(-r^4/4) dr = pi^(3/2)`, vs `2pi` for the Gaussian |
| Position-dependent scalar inertia `I(x)=1+0.5 x_z` | Density proportional to `I(x) dOmega`; `E[x_z]=1/6`, so not Haar |
| Smooth anisotropic tangent metric with determinant one | Haar despite non-round kinetics |

For the last control let `P_x=I-xx^T`, `v=P_x e_z` and
`T_x=vv^T-(v^Tv)P_x/2`. Then `G_x=exp(epsilon T_x)` has normal eigenvalue
one and tangent determinant one, is globally smooth (including the poles),
and is anisotropic where `v != 0`. The computed volume factors agree to
`1e-12`; a separate spherical quadrature verifies the variable-inertia mean.

The other load-bearing choice is **holonomy rather than `N` as the restoring
coupling**. #282's conditioning-variable choice has moved into the
Hamiltonian. #283's `N²/2` countermodel proves this relocation by giving
`E=0` at every temperature on the same regular zero set. Single-valued
`sin²(theta)` and the linear tube parent narrow the family of choices under
their added premises; they do not eliminate the choice. Neither an unchanged
assumption count nor exhaustion of classical alternatives follows.

## 7. Cross-round consistency on the shared triangle family

The sweep uses `gamma = 0.3, 0.7, 1, 1.4, 2, 2.6`, all four sectors, and
32,768 circle points. It compares actual independently implemented quantities
and normalizations, not status flags. For the singlet itinerary
`u=s_A a`, `w=-s_B b`, so `c=-s_A s_B cos(gamma)`. Round 5's raw convention
is converted with `pi-gamma`. Its full phase is `Omega=2 theta`, hence its
coarea masses carry an extra factor `1/2`; this cancels only after sector
normalization. The loop's common minus sign also cancels in a normalized
oriented table, without turning its local signed current into a positive
measure.

| Independent overlap | Maximum absolute residual |
|---|---:|
| Quaternion product vs `(D,Nx)/sqrt(D²+N²)` | 9.44e-16 |
| Holonomy trace vs minus quaternion scalar part | 1.50e-15 |
| Morse–Bott `M_0+M_pi` vs equilibrium `M(c)` | 6.62e-9 |
| `M_0-M_pi` vs `2pi(1+c)/|u×w|` | 1.43e-14 |
| Positive normalized sector laws | 7.79e-11 |
| Oriented normalized sector laws | 5.56e-17 |
| Wrapped finite-difference half-phase gradient vs `|u×w|/D` | 6.20e-10 |
| Frame-energy transverse Hessian vs `|u×w|²/D²` | 5.64e-9 |
| Source density and antipodal parity, each | 2.04e-20 |
| Source probability normalization | 2.23e-16 |
| Unequal sector-prior correlation controls | 5.56e-16 |
| Normal-energy correlation and normal-window correlation, each | 0.0 |

Derivative tests omit the singular punctures and test regular points. The
normal-window and normal-energy controls agree in their sector blindness;
their finite-window and finite-temperature distributions are not equated.

At **`c=cos(1)`**, the numerical component masses are
`11.6708058330 + 0.1695122750 = 11.8403181080`, compared with the independent
closed form `M(c)=11.8403181146`. This confirms the review's cross-round
identity to the circle integration error. No disagreements were found in
these defined overlaps. This is not a sweep of every observable in the
repository or a proof of global consistency.

The code-to-documentation check did find two stale labels in round 7. The
variance spread of its default projection is `0.0060880574`, whereas
`0.0102841815` is the spread of `E[|x.m|]`; the write-up had attributed the
latter to the former. The current history-action probe also contains 24
criteria, all passing, while the README still said 21. Both labels are
corrected without altering the underlying numerical functions or frozen
historical artifacts.

## 8. Reproduction, provenance, and remaining derivation

The new hypotheses and criteria were publicly committed at
[`d258bb1`](https://github.com/davidmdrpi/geometrodynamics/commit/d258bb14e73674fd6ecd7ff6f5d2a46edf20eeab)
before the first new implementation. A working-session reset subsequently
lost the uncommitted implementation; it was reconstructed from that public
freeze and the recorded prior results, then rerun. The current implementation
is therefore a reconstruction with the earlier outcomes known, not a second
blind test. The archived values here are from the rerun.

#283's original freeze is a separate provenance class: local `76ed50e`,
published later as preregistration-only `d83d46a`, then implementation
`7d70728`. Its public timestamps do not independently establish that the
freeze preceded computation. That limitation is unchanged by this follow-up.
Historical preregistrations and archived probe outputs have not been edited.

```bash
python -m pytest -q tests/test_source_readout.py
python -m experiments.closure_ledger.source_readout_probe \
  --output-dir experiments/closure_ledger/runs/20260905_source_readout_probe
```

The source-readout suite now has 23 tests (21 original, two post-review
blind-projection controls); the original probe has 46 numerical criteria
and exits nonzero on failure. See the
[archived report](../experiments/closure_ledger/runs/20260905_source_readout_probe/probe.md)
and [full numerical evidence](../experiments/closure_ledger/runs/20260905_source_readout_probe/probe.json).
The parity correction also updates the existing history-action tests and
probe without changing its earlier numerical values.

The next decisive construction needs one common model containing the frame
family, local field solution, apparatus and allowed future interventions.
It must compute the early-record law with that apparatus present. An
informative admissible law would establish the causality failure; a theorem
restricting every admissible law would establish non-readability. The current
work establishes neither conclusion. Composition and the Born rule remain
separate unresolved derivations.
