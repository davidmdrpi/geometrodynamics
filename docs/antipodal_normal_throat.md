# At fixed angular shift, normal reversal blocks net face transfer

The four-cell comparison separates two effects. **The theta action selects
whether total circulation is allowed.** Both angular-shift cells allow it;
both angular-reflection cells force zero total on compatible fields.
**At fixed shift, reversing the normal forbids the net face-to-face angular
momentum transfer that the normal-kept control allows.** Normal reversal
is compatible with circulation; it does not enable it.

This corrects the original headline at `5bb26f2`, which compared maps
that changed both actions. The reviewer supplied the missing controls and
strong-response results before this revision. The
[review amendment](antipodal_normal_throat_review_amendment.md), published
at `31cfc67e5ef92286c3893f2a8667845cf59b544d` before its implementation,
records that prior knowledge and fixes the additional runs and gates.
These are confirmations and stronger validation, not blind discoveries.

**20 amended checks and 45 tests pass.** The experiment retains the original
32 trajectories, adds 64 generic compatible trajectories, and tests a strong
opening of 1.065 against DOP853. The original [freeze](antipodal_normal_throat_prereg.md)
(`533f9257fcb28f015740b9f08557959bc3c1a83b`), runner and first archive are
unchanged; their original 11 checks and thresholds remain reproducible.

## Four geometries with the same Hamiltonian

On the cover T2 x [-a,a], every map advances phi by pi and preserves time.
Let c=+1 denote theta+pi and c=-1 denote -theta; r is the normal sign.
All four maps are free involutions and canonical symmetries of the same H.

| API name | identification (angles modulo 2pi) | surface / throat determinant | total L_theta on compatible fields | net face transfer |
|---|---|---|---|---|
| `shift_kept` | (phi+pi, theta+pi, u) | +1 / +1 | allowed | allowed |
| `antipodal` (A) | (phi+pi, theta+pi, -u) | +1 / -1 | allowed | forbidden |
| `reflective` (R) | (phi+pi, -theta, u) | -1 / -1 | zero | forbidden |
| `reflect_reversed` | (phi+pi, -theta, -u) | -1 / +1 | zero | allowed |

Thus the two non-orientable cells forbid net transfer in this family and
the two orientable cells admit it. This correspondence is **not a theorem
about arbitrary non-orientable geometries or conserved generators**.
At fixed reflection, normal reversal instead makes opposite-face transfer
possible while the total remains zero. The normal action matters jointly
with the action on the generator and the face permutation.

The shift central surface quotient is a torus (coordinates Phi=2phi,
Theta=theta-phi); the reflective central surface is a Klein bottle. A's
normal line reverses around a quotient loop, making that torus one-sided
in its three-dimensional neighborhood. Uniform nonzero theta flow descends
for either shift and fails to descend for either reflection. This does
not forbid a nowhere-zero full flow such as d_phi in the reflective case.

The half-turn gluing is monodromy, not a prescribed connection holonomy,
an internal conjugating lift, J squared equal to -1, or time reversal.
This is a local torus model, not an identification with the S3 construction.

**For either normal-reversed cell, the two normal faces are cover
representatives of one quotient boundary component.** Its compatible
opening has ONE independent scalar coordinate: q_plus=q_minus and
p_qplus=p_qminus. The normal-kept controls retain two boundary components
and may have two independent openings. None supplies an exterior
attachment map for a pair of physical mouths. In particular A does not
model two independently moving exterior mouths.

## Real wave and responding opening

Retain a continuous theta harmonic of one real field,

    F = a(phi,u,t) cos(n theta) + b(phi,u,t) sin(n theta).

With angular normalization absorbed into a,b, discretize phi periodically
and retain two normal samples coupled by kappa. For x=(a,b),

    H = 1/2 ||p||^2 + 1/2 sum_c c^T K(q) c
        + sum_s [p_qs^2/(2M) + K_s q_s^2/2],
    K(q) = L_phi + n^2 + mu^2 + kappa L_normal + V(q),
    V_s(phi) = U exp[-2 h(phi) q_s],
    M qddot_s = -K_s q_s + sum_phi h V_s (a_s^2+b_s^2).

Here h=0.6+0.4 cos(2phi), U=2, mu^2=0.25, kappa=0.75, M=1 and K_s=4.
The opening load is a derivative of this Hamiltonian. The potential,
support and two-sample normal coupling remain chosen effective inputs.
q is an outward-opening scalar, swapping without a sign when r=-1;
oriented normal readouts have signs (+,-). Their cancellation does not
solve translational recoil or a bulk momentum constraint. In particular,
kappa is fixed, not derived from a changing proper normal distance.

Spatial theta rotations of the cover give the conserved Noether charge

    L_theta = n sum(a p_b - b p_a).

This is angular momentum of the real wave. Under any of the four deck
maps, L(Dz)=c L(z), independently of r and the scalar sign sector eta.
Consequently reflection-fixed data have zero total; shift-fixed data may
carry a nonzero total. All four maps preserve time and are canonical.
Their fixed spaces are symplectic, rather than the time-reversing
Lagrangian halves studied in PR #352.

## Why the transfer table follows

After integrating over phi, the exact normal-link flux into the plus face is

    I_plus = kappa n sum_phi(a_plus b_minus - b_plus a_minus),
    d L_plus/dt = I_plus,       d L_minus/dt = -I_plus.

The two cases of the face action give the full comparison:

- **Normal kept:** each face satisfies L_s=c L_s. Reflection therefore
  makes each face value zero and constant. Shift imposes no such
  constraint, and the compatible trajectories exhibit transfer.
- **Normal reversed:** the face values satisfy L_minus=c L_plus.
  For shift, conservation of their sum makes each a constant half of the
  total. For reflection, the face values are opposite; their sum is zero
  but their magnitudes can evolve, as the control trajectories show.

For A this is also a direct cancellation: the minus-face amplitudes are
eta*(-1)^n times the phi half-shift of the plus face. The half-shift is
self-adjoint, so the two bilinear terms in the integrated flux cancel.
Local flux can remain nonzero. This is the time-kept,
generator-commuting case of PR #352's boundary pairing rule.

The original post-freeze exploration already found this obstruction for
A: local flux reached 7.26e-4, net flux stayed below 9.22e-19, and face
momentum changed by less than 1.39e-17. A non-compatible control had initial
flux 0.12 and face change 0.0474. The review supplies the missing causal
comparison: an orientable, shift-preserving control allows circulation
AND transfer under the same Hamiltonian.

## Normal channel selection and frozen response

For A only, the scalar sign sector eta selects normal parity epsilon_u by

    (-1)^(m+n) epsilon_u = eta.

At n=1,m=0,eta=+1, opposite field amplitudes on the two faces give equal
quadratic loading and equal angular momentum. Opposite field sign is not
opposite electric charge. The frozen frequency is

    omega^2 = 4 sin^2(pi m/N)/(2pi/N)^2 + n^2 + mu^2 + U
              + kappa(1-epsilon_u).

The independently assembled grid operator matches all 160 resolved real
modes in the original spectral check to 3.91e-14. The amended
`phi_symbol_convergence` check covers m=1,2,3 and n=1,2 at N=16,24,32,48,
with orders 1.963 to 1.996. It checks the phi Laplacian symbol, **not normal
refinement or convergence of nonlinear continuum dynamics**.

The frozen two-face response is the inverse of

    Z_m(Omega) = [[lambda_m+kappa-Omega^2, -kappa],
                  [-kappa, lambda_m+kappa-Omega^2]],

where lambda_m excludes normal coupling. The amended probe uses 32
one-face sources, measuring BOTH face responses in independent full-grid
solves. Its maximum discrepancy is 2.00e-15. This excites both normal
channels rather than repeating a single-eigenvector spectral check.
Changing kappa shifts only the odd-normal pole. This remains a lattice
susceptibility, not a continuum scattering or Dirichlet-to-Neumann map.

## Amended responding runs

The original 32 travelling-mode runs at amplitudes 0.05,0.1,0.2,0.4,
m=0,1, n=1,2, N=16,24 begin at q=p_q=0 and run to t=40. Their largest
opening is 0.0614693, relative momentum drift 6.05e-15, and sampled relative
energy excursion 1.53e-6. The amendment tightens their energy gate from
2e-3 to 2e-5 without rewriting the original frozen gate. The original
independent DOP853 comparison improves by 4.0002 on halving the timestep.

The new generic runs cover (N,n)=(16,1),(24,2), all four maps, both eta
signs and four independent draws per sign: **64 trajectories** to t=40,
h=.005, recording every .1. Field positions and momenta are independently
drawn and projected, then jointly normalized to norm .4. Opening positions
and momenta are random, nonzero and projected. Every draw has at least
71.9% of initial field-phase-space Fourier power outside its strongest
phi bin. The smallest initial opening momentum norm is 0.00558.

Below S=||x(0)||^2+||p(0)||^2=0.16; Delta L_s is the largest sampled
absolute change of either face value. Each row aggregates eight draws;
using the minimum Delta L_s requires transfer in EVERY allowed draw.

| N,n | cell | max abs(L_total)/S | min Delta L_s/S | max abs(I_plus)/S |
|---|---|---:|---:|---:|
| 16,1 | shift_kept | 0.19212 | 0.076958 | 0.092601 |
| 16,1 | antipodal | 0.15351 | 5.60e-16 | 3.32e-16 |
| 16,1 | reflective | 1.81e-15 | 3.62e-16 | 4.36e-16 |
| 16,1 | reflect_reversed | 2.66e-15 | 0.058403 | 0.10329 |
| 24,2 | shift_kept | 0.26451 | 0.094527 | 0.12382 |
| 24,2 | antipodal | 0.13928 | 1.62e-15 | 9.30e-16 |
| 24,2 | reflective | 4.99e-15 | 1.78e-15 | 1.14e-15 |
| 24,2 | reflect_reversed | 6.38e-15 | 0.24298 | 0.15252 |

Across all generic trajectories, maximum normalized total-momentum drift
is 6.38e-15, seam error 1.12e-14, and sampled relative energy excursion
3.19e-4 (registered bound 5e-4). Generic high-frequency data have a larger
energy error than the original low-mode runs; the two are reported
separately. Maximum forbidden face variation across all A and R draws is
4.40e-15 in these normalized units; the archive records every draw.

These symmetry-preservation checks validate the implementation and
integration, not independent dynamical evidence for a symmetry already
proved. Negative controls evaluate an actual perturbed, non-compatible
trajectory (seam error 0.684) and an artificial repeated-state history.
The seam gate rejects the former; the required transfer witness rejects
the latter. Missing cells, missing required checks and duplicate draws
cannot yield a passing verdict.

The amended geometry gate differentiates the actual wrapped coordinate
map and reconstructs its action on the real field, rather than comparing
two hard-coded Jacobians. Maximum Jacobian error is 1.03e-9, pullback
error 3.34e-16; all four state maps are canonical and preserve H.

## Strong response of the single A-compatible opening

For amplitude 2.5, N=16,n=1,m=0, q=p_q=0, integrate to t=40.

| measurement | result |
|---|---:|
| largest opening | 1.0649244 |
| relative angular momentum drift | 3.13e-15 |
| normalized seam error | 3.13e-15 |
| sampled relative energy excursion (h=.01) | 1.043e-5 |
| final-state error vs DOP853, h=.01 | 0.00147072 |
| final-state error vs DOP853, h=.005 | 0.000367828 |
| error reduction | 3.99839 |

The response is substantial within the effective Hamiltonian while the
one independent opening remains compatible with circulation and the net
transfer obstruction. Neither the two cover coordinates nor their
identical time traces constitute independent mouth motion.

## Reproduce and pursue

```bash
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_antipodal_normal_throat.py tests/test_antipodal_normal_throat_review.py
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.antipodal_normal_throat_review_probe \
  --output-dir /tmp/antipodal-normal-throat-review
# Original experiment and original thresholds:
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.antipodal_normal_throat_probe \
  --output-dir /tmp/antipodal-normal-throat-original
```

Archives: `experiments/closure_ledger/runs/20261112_antipodal_normal_throat/`
(original) and `experiments/closure_ledger/runs/20261112_antipodal_normal_throat_review/`
(amended). The amended runner includes the original calculations and exits
nonzero for a failed or missing required gate. Focused tests use the archived
generic rows to exercise gate failures; the full runner independently
recomputes all trajectories.

For BAM, the constructive lesson is that compatible circulation, selected
normal channels and mechanical response coexist. The sharper restriction
is on the proposed transfer mechanism: in this four-cell family, transfer
between faces occurs only in the orientable cells. One cannot obtain
independent mouth exchange just by strengthening the A opening response.
The next task must specify two distinct exterior ports, the conserved
generator and its time/lift action, derive their boundary balance, and
refine the normal geometry. A different attachment or time action could
change that balance; this experiment has not tested it.

Amplitudes and angular momentum remain continuous. No action quanta,
spin lift, electric charge, probability law or Einstein solution are
derived. The construction is not Stage 1 of PR #351.
