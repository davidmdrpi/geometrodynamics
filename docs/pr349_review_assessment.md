# Assessment of the October 25 review of PR #349

Review: [comment 6402240493](https://github.com/davidmdrpi/geometrodynamics/pull/349#issuecomment-6402240493),
against `64810a8`. This assessment corrects interpretation of the archived
results and adds targeted replications. It leaves the original manifests and
negative scalar results intact. No BAM force, source or capture law is changed.

## Findings accepted and checked

**P2-1: return cancellation is topology-specific.** The untwisted pair's inner
region is the heavy side for the prepared channel at both defects. Transmitted
waves transfer momentum toward that side on either crossing. The twisted
single-kink return exchanges channel labels, reversing which side is heavy;
its cancellation must not be generalized to an untwisted pair.

New runs use the original amplitude 0.2, speed 0.02, separation 20, dx=0.1,
dt=0.02, with records every 0.5 time units:

| Circumference | First separation below 2 | First loss of two zero crossings |
| ---: | ---: | ---: |
| 400, original archive | 161.5 | 162.5 |
| 80, new run | 119.5 | 120.5 |
| 60, new run | 100.0 | 101.0 |

The no-wave circumference-80 control continues separating. Returns strengthen
the inward branch in these cases. They do not demonstrate the behaviour of
every packet or topology. The review's coarse-sampled time 184 is not the
first approach: the original pair is already below separation 2 at 161.5 and
has no zero crossings at 164. A detector requiring exactly two roots can miss
the encounter and report a later reappearance. We record both events separately.

**P2-2: the energy/speed mismatch is substantive.** The requested small-seed
test uses radiation energy `0.006*(8/3)=0.016`, with continuum pair rest energy
8/3, speed 0.3, circumference 400, and the same packet shape and frequency.
The amplitude is solved from the actual discrete energy, not scaled by guess.
A no-wave control removes the lattice's own velocity drift. Late velocities
are fitted over t=85 to 100, after the first encounter and before exterior
circulation returns.

| dx / dt | Right velocity, with wave | No-wave velocity | Wave-induced change |
| --- | ---: | ---: | ---: |
| 0.1 / 0.02 | 0.29820408 | 0.29967077 | -0.00146669 |
| 0.05 / 0.01 | 0.29845858 | 0.29991787 | -0.00145929 |

The wave-induced change agrees to about 0.51% across refinement. The pair
keeps separating. Using the continuum translation estimator `P=(4/3) gamma v`
gives a refined momentum decrement 0.0022396 per defect, about 0.53% of its
initial 0.4193. This is an inferred translational momentum, not an independent
integral of the complete field momentum. The refined wave run's maximum
recorded relative energy error is 7.81e-8. These are analogue controls at the
requested ratio, not an exact matching of BAM's masses, sources or dispersion.

The review's order-of-hundreds estimate is a useful first-pass budget. It is
not a universal impossibility theorem for repeated conservative exchange:
crossings can recycle radiation and draw on mechanical/internal energy. The
available encounters, work and full-cycle field history must be calculated.

**P2-3: the prepared channel is unsourced.** Both psi equations are homogeneous
in psi. With psi=psi_t=0, a moving phi defect cannot generate either psi field.
The no-wave runs retain exactly zero psi, and a regression tests that invariant
subspace. Phi can radiate its own disturbances; that does not source the two
prepared channels. A linear source coupling would change the action, its
support and its deck symmetry. It is a new model requiring derivation, not a
parameter correction that this review should silently introduce.

**P3-4: exact quotient symmetry is a construction check.** The zero residual
follows from compatible initial data and the equivariant reversible update.
It checks implementation consistency, not independent physical closure or
unique selection. The useful result is the resulting history; ordinary
energy/refinement and full-field reversal checks remain separate evidence.

## Carry-over clarifications

- The original directional table already separates sources. At A's first
  antipode, A alone gives centre field **0.12828852**, rim mean **0.06132605**.
  B contributes **0.04955529** and **0.04487740** respectively. The headline
  0.17784381 is explicitly the coherent pair, not A's isolated focused field.
- Paying for the flash changes kick speed 1/3 to **0.31497141**, as independently
  solved from the existing PairBirth energy coefficients. Contact then occurs
  at **8.06928050**, rather than the fixed-speed diagnostic's 7.62478. Both
  first-flight controls capture zero fronts before contact.
- Zero captures on that first flight does not imply zero captures across
  rebirths. PR #345 at `d5921ff` explicitly records later-cycle captures. Its
  current capture and flash-payment code was inspected; its full multicycle
  resonance runs were not independently replayed here.
- The directional report already declines to choose a finite window as a
  physical impulse. The review's whole-passage numbers do not specify window
  endpoints in the comment; those exact four values are not claimed reverified.

## The proposed quotient route is useful, but does not yet select speeds

On a static round sphere times time, the free involution

\[
 J(n,t,q)=(-n,t_c-t,-q)
\]

squares to the identity for any t_c. This is different from the tested
Schwarzschild boost-offset family. For the prescribed trajectories
`theta_A=a+ut`, `theta_B=-a-ut`, taking `t_c=(pi-2a)/u` maps A to B exactly.
The twelve checked histories include resonant and nonresonant speeds; all
obey this worldline relation to rounding precision. Thus the *geometric*
quotient condition itself does not select u. It can reinterpret the two cover
mouths as one quotient history; it does not derive a force between two
independent objects or establish a birth/contact surgery law.

The elliptic de Sitter comparison is structurally relevant: [Parikh, Savonije
and Verlinde](https://arxiv.org/abs/hep-th/0209120) use antipodal CPT
identification and discuss time nonorientability. The static cylinder here is
not their de Sitter metric. In an evolving FRW metric a time-reflection
isometry additionally needs `a_FRW(t)=a_FRW(t_c-t)` on its domain; arbitrary
clock origins do not generally satisfy that condition.

Exact partner-centre focusing *at contact* gives the proposed values

\[
 u_k=\frac{\pi-2a}{(2k+1)\pi}.
\]

At a=0.3 these are 0.80901407, 0.26967136, 0.16180281 and 0.11557344.
However, the existing aperture rule admits partner capture *before* contact
over the interval

\[
 \frac{\pi-3a}{(2k+1)\pi}\le u<
 \frac{\pi-2a}{(2k+1)\pi},
\]

provided an earlier capture has not already removed the front. The actual
router confirms two interior speeds for every k=0,1,2,3. At k=1, for example,
both **0.26171361** and **0.24579811** capture the partner's flash at 3pi.
The exact centre resonance is the contact endpoint, which the present
first-flight audit excludes. The full model therefore needs a specified
simultaneous focus/contact boundary rule. Moving slightly below resonance
tests near-contact capture, not exact equality of the birth and absorption
events under J.

Consequently neither finite-aperture capture nor the coasting quotient proves
a discrete speed spectrum. A complete field/source/absorption matching could
impose additional constraints; that is the next problem to solve.

Nor does time nonorientability uniquely force a half-retarded, half-advanced
field. If J exchanges the retarded and advanced Green operators, their half
sum is an invariant particular choice for a compatible source. An invariant
homogeneous solution can still be added. For a commuting linear wave operator,
`H_invariant=(H+J*H)/2` is such a solution whenever nonzero, with the appropriate
field parity included in J*. Boundary data or an absorber condition must fix
this freedom. On the scalar round-sphere example an odd-l harmonic times
`cos(omega*(t-t_c/2))` already supplies a nonzero charge-odd homogeneous family
for every t_c. This illustrates nonuniqueness, not a new Maxwell solution.

## Neck stress and the next calculation

The transmitted part of a wave pushes toward the heavy side. That is not by
itself the sign of the **total** neck-radius force: reflection, evanescence,
interference and shape work must be included. A lower gap at a larger radius
also changes the radius derivative of the field energy. The m=0 neutral
result in #348 cannot simply be relabelled as the m=1 shape force.

The highest-priority continuation is a quotient-compatible directional field
history with explicit sources, homogeneous data and focus/contact ordering,
followed by the stress/shape-work projection on that history. The review's
charged-background dispersion idea remains a hypothesis: linear Maxwell
propagation on a fixed metric does not acquire a sign-dependent dispersion
merely by adding a static Maxwell background.

These results narrow and correct the interpretation of PR #349. The PR stays
open; no new GR solution, quantization or sustained BAM collapse is claimed.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.pr349_review_probe
python -m pytest -q tests/test_pr349_review.py
python -c 'from experiments.closure_ledger.pr349_review_probe import manifest; manifest()'
```

The seven new histories are diagnostics sampled every 0.5 time units, not full
state replay archives. The original fifteen replay archives are unchanged.
The new manifest binds all diagnostics, the producer, tests and source models.
