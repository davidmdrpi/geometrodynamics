# Finite pointer spread: the record distinction survives the coupled history test

The preregistered primary experiment passes: **nonzero pointer momentum
spread weakens the early-record distinction but does not erase it** in the
specified classical rotor history model. With source momentum width `0.1`,
pointer momentum width `0.1`, and position resolution `0.15`, the same early
event has probabilities approximately **0.158 and 0.458** for the two future
settings, a contrast of **0.300**. The actual pointer momentum marginal is
the same Gaussian in both ensembles. The source recoils, and the history
weights include that changed trajectory and its parallel transport.

This answers the finite-width objection **within an explicit extension** of
the triangle model. It does not construct the missing BAM frame-to-field
map or a physical realization of its history preparation and future-setting
control. Nor does it derive the Born rule. The rotor's local accessibility,
Hamiltonian interaction, closure weighting and preparation law are stated
inputs. A Gaussian phase-space preparation is classically admissible here;
its apparatus-level realization in BAM is still unproved.

The [public freeze](pointer_spread_prereg.md) was committed at
[`baf856d`](https://github.com/davidmdrpi/geometrodynamics/commit/baf856df06b26de6c9750697ba5df5ad521000c5)
before this implementation or any spread calculation. It fixes the primary
width, event, threshold, preparation rules, and resolution checks. No model
parameter or primary threshold was changed after seeing the results.

## 1. A complete finite-pulse experiment

Use a unit-inertia rotor on `S²`, with unit direction `x`, tangent canonical
momentum `p`, and a canonical pointer `(Q,P)`. The interaction is

\[
 H(t)=\frac{|p|^2}{2}+\frac{P^2}{2M}+h(t)P(m\cdot x),\qquad
 h(t)=\frac{2g}{\tau}\sin^2(\pi t/\tau),\quad 0\leq t\leq\tau,
\]

and `h=0` outside the pulse. Fix `M=g=tau=1` and
`m=(1,2,0)/sqrt(5)`. The embedded Hamilton equations are

\[
 \dot x=p,\qquad
 \dot p=-|p|^2x-hP[m-(m\cdot x)x],\qquad
 \dot Q=P/M+h(m\cdot x),\quad\dot P=0.
\]

The reciprocal source force is retained. The readout `Y=Q(1)` is recorded
at the end of the pulse. The source continues freely until the future
closure boundary at `t_B=2`. Later pointer drift is not added to this early
record, and no final pointer position is prescribed or postselected.

The reference source preparation is uniform area for `x0` and an isotropic
tangent Gaussian for `p0`, with standard deviation `sigma_s=0.1` in each
orthonormal component. The pointer has independent reference variables
`P ~ Normal(0,sigma_P²)` and `Q0 ~ Normal(0,0.15²)`. The main experiment
therefore does not set either source or pointer momentum identically to
zero. These widths and the closure parameter below are independent
preparation parameters; no common thermal bath or equilibration mechanism
is claimed.

Fix `a=e_z` and future choices `b_0=e_x`, `b_1=e_y`. The event is always
`B={|Y|>0.6}`. Gaussian integration over `Q0` uses the same response kernel
for both choices. It adds no outcome-dependent detector rule.

## 2. Close the actual path, including the measurement-induced motion

Let `C` be the rotor's Levi-Civita parallel transport from `x0` to its final
direction `x_B`. It obeys the quaternion equation

\[
 \dot C=\tfrac12[0,x\times p]C,\qquad C(0)=1.
\]

The full closed itinerary is the source trajectory from `x0` to `x_B`,
followed by the geodesic analyzer legs `x_B -> u -> w -> x0`, where
`u=s_A a` and `w=-s_B b`. Its lift is

\[
 G_s=\operatorname{lift}(w,x_0)\operatorname{lift}(u,w)
     \operatorname{lift}(x_B,u)C.
\]

The final leg returns to **the original source direction**. Neither the
motion before readout nor its accumulated transport is discarded. In
particular, replacing this expression by a triangle based only at `x_B`
would be a different history model. The test suite also checks that two
paths with identical endpoints and different transport receive different
closure energies.

Use the inherited positive frame penalty and equal reference sector
coefficients:

\[
 V_s=\frac{\|\operatorname{Ad}_{G_s}-I\|_F^2}{16}
     =\frac{1-G_{s,0}^2}{2},\qquad
 W_b(x_0,p_0,P)=\frac14\sum_{s_A,s_B}e^{-64V_s}.
\]

This is a **chosen history weighting**, not a probability produced by
real-time Hamiltonian evolution alone. The common Pin-loop minus sign
does not affect this frame energy. The antipodal geodesic punctures remain
undefined and have zero reference measure; the implementation raises an
error at a puncture instead of assigning it a preferred weight.

The closure is soft and finite (`beta=64`). We have not solved BAM's full
spacetime field boundary problem, nor taken an exact hard-closure limit
with a finite-width apparatus. The construction specifies which rotor
histories are weighted and how their early records are calculated.

## 3. Keep the prepared pointer marginal distinct from a reference prior

Let `pi_s` denote the source reference measure and `nu_sigma` the pointer
Gaussian. Define

\[
 Z_b(P)=\int W_b(x_0,p_0,P)\,\pi_s(dx_0\,dp_0).
\]

The **primary preparation law** fixes the actual pointer momentum marginal:

\[
 \mu_b^{\rm fixed}
 =\nu_\sigma(dP)\,\pi_s(dx_0\,dp_0)\frac{W_b}{Z_b(P)}.
\]

One can sample `P` from the specified Gaussian and then sample a source
history conditional on that `P` and the chosen boundary setting. This
defines a normalized mathematical preparation, including its correlations.
Per-`P` normalization is an explicit input. It does not derive an apparatus
or explain how BAM physically implements that source-history law.

The **joint-conditioning control** instead uses

\[
 \mu_b^{\rm joint}
 =\frac{\nu_\sigma(dP)\,\pi_s(dx_0\,dp_0)W_b}
        {\int\nu_\sigma(dP)Z_b(P)}.
\]

Its posterior pointer density is proportional to `nu_sigma(P) Z_b(P)`.
Calling its unchanged *reference* Gaussian an unchanged *actual* preparation
would be incorrect. The experiment reports the difference:

| Pointer width | Fixed-law variance of P, either setting | Joint-law variance, b_0 | Joint-law variance, b_1 |
|---:|---:|---:|---:|
| 0.1 | 0.010000 | 0.009843 | 0.010137 |
| 0.5 | 0.250000 | 0.213760 | 0.263195 |

The mean of `P` is zero in both laws by symmetry. The primary result thus
does not rely on a change in the prepared pointer's momentum distribution.
The future-setting dependence is in the conditional source histories and
their recorded interaction with that pointer.

For either law the early-event probability is

\[
 L_b(B)=\int K(B\mid Q_{\rm shift}(x_0,p_0,P))\,d\mu_b,
\]

where `K` convolves the computed pointer shift with the fixed `Q0` Gaussian.
No later pointer record is used for selection. Nevertheless, physical
realization of these conditional families under independent future-setting
interventions in BAM is not established by their mathematical definition.

## 4. Results of the frozen sweep

Define the signed contrast `Delta=L_{b_1}(B)-L_{b_0}(B)`. The frozen-posterior
control applies the **actual coupled pointer record** to weights computed
from the freely moving source path with the backreaction force removed.
It isolates the effect of omitting the instrument's change to the source
history ensemble; it is not a second valid derivation of the coupled law.

| sigma_P | Fixed-marginal contrast | Scramble SE | Joint-conditioned contrast | Frozen-posterior contrast | Reference source momentum recoil RMS |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.3430 | 0.0020 | 0.3430 | 0.3430 | about 6e-16 |
| 0.01 | 0.3424 | 0.0020 | 0.3424 | 0.3428 | 0.00816 |
| 0.025 | 0.3398 | 0.0018 | 0.3398 | 0.3419 | 0.02041 |
| 0.05 | 0.3310 | 0.0015 | 0.3310 | 0.3385 | 0.04081 |
| **0.1** | **0.2996** | **0.0012** | **0.2997** | **0.3260** | **0.08161** |
| 0.25 | 0.1690 | 0.0004 | 0.1732 | 0.2640 | 0.20382 |
| 0.5 | 0.0237 | 0.0008 | 0.0503 | 0.1598 | 0.40624 |

These are reference-grid estimates; their last digits are not claimed as
precision results. Scramble SE describes variation across four numerical
scrambles, not a rigorous confidence interval. Independent resolution
changes are reported below.

Recoil is the difference from the zero-P source trajectory at readout,
averaged over the common source and pointer reference measure before closure
weighting. It diagnoses the actual force retained in the path calculation;
it is not a separately predicted conditional recoil distribution.

At the preregistered primary width `sigma_P=0.1`, all four fixed-law
contrasts are between **0.2961 and 0.3014**, exceeding the frozen criterion
`Delta>0.1`. The mean event probabilities are **0.1584 and 0.4580**.
The RMS momentum recoil `0.0816` is substantial relative to the source's
component width `0.1`; the result is not obtained by making the reciprocal
force identically zero.

The zero-pointer-width row still has a moving source (`sigma_s=0.1`), so
its contrast `0.3430` need not equal #284's stationary-source value. In the
separate control with **both** momentum widths zero, the event probabilities
are `0.12299` and `0.48977`, compared with #284's independent sphere
quadrature `0.1231586` and `0.4897143`. That is the proper bridge to the
earlier model, rather than silently identifying the two zero-pointer rows.

Reweighting matters increasingly at larger spreads. At `sigma_P=0.5`, the
fixed-law contrast is about `0.024`, while holding the previous source
posterior fixed would report about `0.160`. Joint conditioning also differs
from preserving the actual prepared marginal. None of these distinctions
is captured by adding noise to the old readout without revisiting the paths.
No universal erasure threshold or monotonicity theorem is inferred from
seven sampled widths. Further source convergence would be needed before
making a precise disappearance claim at larger widths.

## 5. Mechanics and numerical controls

The algorithm alternates exact spherical geodesic flow and transport with
exact canonical kicks from `h(t_mid) P F`, using Strang splitting. It reads
the pointer at `t=1` and transports the source freely to `t=2` exactly.
An independent DOP853 integration uses the embedded Hamilton equations,
the transport differential equation and the external-work integral.

| Check | Result |
|---|---:|
| Nonzero-P ODE comparison, 128 pulse steps | largest component error 2.02e-6 |
| Error ratios under step halving, nonzero-P controls | 3.993–3.999 |
| Independent work balance `Delta H = integral h' P F dt` | residual below 1.44e-14 |
| Quaternion products vs existing scalar quaternion machinery | 8.89e-16 |
| Stationary augmented-loop energy vs #283's triangle energy | 1.67e-16 |
| Sphere, tangency and axial angular momentum, entire sweep | worst below 2.89e-14 |
| Full loop fixes x0, entire sweep | residual below 1.19e-12 |
| Source effective sample size at any P node/setting | at least 1,091 of 4,096 |
| Primary: 64 to 128 pulse steps | event probability movement below 6.50e-7 |
| Primary: 16 to 32 Hermite nodes | event probability movement below 2.95e-4 |
| Primary: double source base points | event probability movement below 0.00371 |
| Width 0.5: 16 to 32 Hermite nodes | event probability movement below 0.00268 |

The event-refinement bounds above are maxima across both settings, all four
scrambles and all three weighting prescriptions. The source resolution
dominates the primary numerical uncertainty, far below the margin over its
frozen criterion. Refinements hold all physical parameters fixed.

Simultaneous `(x0,p0,P)->(-x0,-p0,-P)`, with sector relabelling, preserves
the weight and reverses the record. Means vanish and smoothed sign
probabilities stay `1/2`; this does not imply equality of the full laws.
A constant kernel erases the distinction. A separate rotation-covariant
instrument with `m=a` gives a nondegenerate blind-record control for the
two future choices. The original warning against claiming that every odd
projection is informative is therefore retained with dynamics present.

## 6. What is closed and what is still open

The exact-`P=0` limitation of the earlier demonstration is removed **in this
coupled rotor extension**: finite Gaussian source and pointer momentum
widths coexist with an informative finite-resolution record, with reciprocal
force, full path holonomy, and updated conditional source weights included.
It is a counterexample to automatic erasure by any nonzero pointer spread
within the stated class of classical history laws.

It is not a universal non-readability result, and it is not a completed BAM
falsification. The distinction between a conditional statistic and a
physically available intervention remains material. The next missing objects
are an actual source-field family realizing this rotor variable, a derivation
or exclusion of this interaction and preparation, and the corresponding
operational future-setting experiment. If that realization admits this
early-record law with independently selectable future settings, it predicts
a retrocausal channel. Neither a missing implementation nor this conditional
countermodel decides whether BAM admits it.

| Ingredient | Status in this experiment |
|---|---|
| Spin-frame holonomy and geodesic transport formulas | Inherited geometry, applied to the augmented path |
| Rotor inertia, direct local access to x and pulse coupling | Chosen classical extension; no field derivation supplied |
| Gaussian phase-space preparation and its widths | Specified smooth reference data; apparatus realization open |
| Positive finite-beta frame weighting and equal reference sectors | Inherited choices, not forced by real-time Hamiltonian flow |
| Preservation of the prepared pointer marginal | Explicit primary preparation rule; compared with joint conditioning |
| Conditional readout at finite spread | Calculated with the instrument present; primary hypothesis passes |
| Operational BAM signal, Born rule, composition | Not derived |

## Reproduction and evidence

```bash
python -m pytest -q tests/test_pointer_spread.py tests/test_source_readout.py
python -m experiments.closure_ledger.pointer_spread_probe \
  --output-dir experiments/closure_ledger/runs/20260906_pointer_spread_probe
```

The new module has 14 tests; together with the preceding source-readout
suite, 37 pass. The probe passes **25/25** frozen criteria and exits nonzero
on failure. Its [report](../experiments/closure_ledger/runs/20260906_pointer_spread_probe/probe.md)
and [full numerical archive](../experiments/closure_ledger/runs/20260906_pointer_spread_probe/probe.json)
retain all replicas, posterior moments and independent refinements.
The integrated focused suite spanning this module and the preceding five
closure/history modules plus source readout passes **158/158** tests.
This experiment used the original public freeze without reconstruction or
post-result parameter changes. Historical artifacts from #283/#284 remain
unchanged.
