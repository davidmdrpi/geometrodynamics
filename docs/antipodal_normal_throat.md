# Normal reversal permits circulation through mechanical response

The antipodal normal twist admits nonzero angular momentum of a **real
wave**, and that momentum survives a mechanically responding waist. A
reflective twist forbids the same nonzero total angular momentum on its
deck-fixed fields. The distinction is reversal of the angular generator,
not the sign of the three-dimensional orientation alone.

All **11 preregistered checks pass**, across 32 trajectories at two phi
resolutions and two angular harmonics. There are 25 tests. Predictions were
published before implementation in [the freeze](antipodal_normal_throat_prereg.md),
commit `533f9257fcb28f015740b9f08557959bc3c1a83b`, on top of merged PR #352
(`d77d915`). The models, support and orienting wave are specified inputs.

The strongest limitation emerged in a **post-freeze exploration**: the
circulation is around theta, not net angular momentum transfer between the
normal faces. Exact deck symmetry forbids the latter for these fields.
That result is reported alongside the positive circulation result.

## Geometry and the meaning of the faces

On the cover T2 x [-a,a], compare

    A(phi,theta,u) = (phi+pi, theta+pi, -u)
    R(phi,theta,u) = (phi+pi, -theta, u).

| property | antipodal A | reflective R |
|---|---|---|
| surface quotient at u=0 | torus | Klein bottle |
| surface determinant | +1 | -1 |
| throat determinant | -1 | -1 |
| pushforward of d_theta | +d_theta | -d_theta |
| uniform nonzero theta velocity | descends | does not descend |
| time extension in this experiment | preserved | preserved |

Both identifications are free involutions on the cover. For A, coordinates
Phi=2phi and Theta=theta-phi show that the central surface quotient is a
torus. Its normal line reverses around a quotient loop: it is one-sided in
the non-orientable three-dimensional neighborhood. R can still carry a
nowhere-zero full flow (d_phi); a forbidden uniform angular component is
not a theorem that the entire fluid or field must stop.

The half-turn in the gluing is **monodromy**, not a prescribed connection
holonomy. Nor does it imply a conjugating internal lift, J squared equal
to -1, or time reversal. The construction is a local torus model, not an
identification with the repository's S3 geometry.

**The two normal faces are cover representatives.** A identifies u=+a
with u=-a, so the quotient of this slab has one boundary component. It is
not by itself a wormhole with two independent exterior mouth boundaries.
The two opening coordinates below are paired representatives of one
quotient opening degree of freedom on compatible histories. Independent
physical mouths require an external port geometry and attachment map.
This distinction is essential when interpreting the equal-and-opposite
normal-load bookkeeping.

## One real harmonic, with a dynamical opening

Retain a continuous theta harmonic of one real field:

    F = a(phi,u,t) cos(n theta) + b(phi,u,t) sin(n theta).

The normalization of the angular basis is absorbed into a,b. Discretize phi
periodically and retain two normal samples, coupled by kappa. With x=(a,b),

    H = 1/2 ||p||^2 + 1/2 sum_c c^T K(q) c
        + sum_s [p_qs^2/(2M) + K_s q_s^2/2],
    K(q) = L_phi + n^2 + mu^2 + kappa L_normal + V(q),
    V_s(phi) = U exp[-2 h(phi) q_s],
    M qddot_s = -K_s q_s + sum_phi h V_s (a_s^2+b_s^2).

Here h=0.6+0.4 cos(2phi), U=2, mu^2=0.25, kappa=0.75, M=1 and K_s=4.
The waist load is a derivative of this Hamiltonian, not an imposed kick.
The coupling of opening to potential is still a chosen effective model;
neither the metric nor a throat normal discretization is solved from GR.

q is an outward-opening scalar and swaps without a sign under A. Oriented
normal readouts have signs (+,-). Their cancellation is not a derivation
of translational recoil or of a bulk momentum constraint. In particular,
kappa is fixed rather than derived from a changing proper normal distance.

Spatial theta rotations give the canonical angular momentum

    L_theta = n sum(a p_b - b p_a).

It is a spatial Noether charge of the real wave, not electric charge or
spin-1/2. A commutes with its generator and preserves L_theta; R reverses
the generator and L_theta. A therefore permits nonzero total circulation;
R-fixed fields have zero total, though they can have local currents.
Both deck maps are canonical because time is preserved. Their fixed
spaces are symplectic, not the time-reversing Lagrangian halves of PR #352.

## The topology selects normal channels

For sign sector eta and normal parity epsilon_u,

    (-1)^(m+n) epsilon_u = eta.

At eta=+1 an odd m+n chooses the antisymmetric normal channel; an even
m+n chooses the symmetric channel. In particular n=1,m=0 has opposite
field amplitudes on the two faces but equal quadratic loading and equal
angular momentum. Opposite field sign is not opposite electric charge.

For a frozen waist the eigenfrequency squared is

    omega^2 = 4 sin^2(pi m/N)/(2pi/N)^2 + n^2 + mu^2 + U
              + kappa(1-epsilon_u).

The independently assembled full-grid operator matches this for all 160
resolved real modes across eta=+/-1, n=1,2 and N=16,24; maximum residual
3.91e-14. The fixed-m continuum convergence order is 1.96 to 2.00.

The frozen response of the two normal faces is the inverse of

    Z_m(Omega) = [[lambda_m+kappa-Omega^2, -kappa],
                  [-kappa, lambda_m+kappa-Omega^2]],

where lambda_m excludes the normal coupling. The deck selects the even
or odd eigenchannel by the parity rule above. The full-grid response and
this two-face response agree to 3.22e-15. Changing kappa moves only the
odd-normal pole. This is a lattice susceptibility, not a continuum
Dirichlet-to-Neumann or scattering map.

## Responding results

Start from a frozen travelling mode: a=A g, b=p_a=0, p_b=omega A g,
q=p_q=0. It is an exact travelling mode of the frozen Hamiltonian only.
Run amplitudes 0.05,0.1,0.2,0.4 for m=0,1 and n=1,2 at N=16,24 to t=40.
Velocity Verlet uses step 0.01. The following are extrema over individual
runs, not bounds inferred from medians.

| measurement | result |
|---|---:|
| largest opening | 0.0614693 |
| largest relative angular momentum drift | 6.05e-15 |
| largest relative energy excursion | 1.53e-6 |
| mouth-response amplitude slopes on [0,4] | 1.979 to 1.987 |
| field-correction amplitude slopes on [0,4] | 2.970 to 2.984 |
| coarse final-state error against DOP853 | 0.001686 |
| halved-step final-state error against DOP853 | 0.0004215 |
| error reduction | 4.0002 |

The seam, equal openings, equal face angular momenta and opposite signed
normal loads pass their separate checks. Conservation permits a deformed
circulating pattern; it does not assert uniform angular speed pointwise.
The normal grid has only two samples: the phi convergence test is not a
normal-direction continuum convergence test.

## Post-freeze exploration: circulation is not net transfer

The exact normal-link angular momentum flux into the plus face is

    I_+ = kappa n sum_phi(a_+ b_- - b_+ a_-),
    d L_+/dt = I_+,       d L_-/dt = -I_+.

For A-compatible fields, the minus face is (-1)^n times the half-shift of
the plus face. That half-shift is self-adjoint, so the two terms cancel
after summing phi. Equivalently, A keeps time and commutes with the angular
generator: PR #352's parity rule gives L_+=L_total/2 at every time.

Add a generic compatible perturbation to excite multiple phi modes, then
integrate on [0,4]: local normal flux reaches 7.26e-4, while net normal
flux stays below 9.22e-19 and face angular momentum changes by less than
1.39e-17. A non-compatible control starts with normal flux 0.12 and
changes a face's angular momentum by 0.0474. The flux routine is not
identically zero: cancellation is imposed by the exact deck symmetry.

This is a useful boundary for the proposed mechanism. Normal reversal,
circulation and mechanical response coexist; this time-preserving,
generator-preserving identification does not generate net exchange of
that angular momentum between its cover faces. To study independent
mouth transfer, specify external mouths or a different time/lift action
and derive the corresponding balance, rather than calling the conserved
circulation a transferred charge.

## Reproduce and pursue

```bash
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_antipodal_normal_throat.py
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.antipodal_normal_throat_probe \
  --output-dir /tmp/antipodal-normal-throat-run
```

Archive: `experiments/closure_ledger/runs/20261112_antipodal_normal_throat/`.
The probe exits nonzero on a failed required check and records the
post-freeze transfer measurement separately.

The next geometric task is to attach two distinct exterior ports and
compute the response there, including a normal-resolution refinement.
The next dynamical task is to decide which physical conserved generator
and time extension govern their exchange. This experiment selects allowed
mode labels but leaves amplitudes, and hence angular momentum, continuous.
It derives neither action quanta, a spin lift, electric charge, a probability
law nor an Einstein solution. It is not Stage 1 of PR #351.
