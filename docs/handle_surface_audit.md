# Trapped necks and the exterior-bulk momentum target

The archived #308 evolution contains **future-trapped neck spheres at
Einstein proper time t=0.01 in both frames**. It does not contain two
independently evolved gravitating receivers. These are separate results.
Trapped mouths remain admissible objects for the research program; the next
interaction target is propagation and response through the exterior bulk.
Traversability is not an acceptance requirement for that target.

This is a **post-hoc snapshot audit**, dated 2026-09-26 UTC, of published
commit `7a9aba6e8ac2c4520075c77f1a14662864f94da3`. The saved evolution was
not rerun. Its raw SHA256 remains
`8a858b6d7ecce0de3047f92a8baa77fd9c262ec9475d6886d539eac26cadf815`.
No old freeze, evolution equation, raw evidence or registered Boolean verdict
has been changed. New results are in
[the audit JSON](../experiments/closure_ledger/runs/20260926_handle_surface_audit/audit.json).

## Null normals, time orientation and frame

Use the repository convention K_ij=-L_n gamma_ij/2 with n future directed.
For ds_E^2=-dt^2+A^2 ds^2+B^2 dOmega^2, take e_s=A^-1 partial_s and
k_+=n+e_s, k_-=n-e_s, so k_+ dot k_-=-2. Then

    theta_E,+ = -2 l + 2 B_s/(A B),
    theta_E,- = -2 l - 2 B_s/(A B).

Both negative means future trapped; both positive means past trapped
(anti-trapped). Reversing the spatial normal only exchanges the labels.
Full time reversal reverses K and scalar normal momenta, and sends
(theta_+,theta_-) to (-theta_-,-theta_+). The antipodal spatial identification
in #308 does not reverse the timelike normal.

The physical Jordan metric is g_J=g_E/f with
f=1-(u^2+v^2)/6, f_t=-(uP+vQ)/3, f_s=-(u u_s+v v_s)/3. Using unit normals
in each frame gives

    theta_J,+ = sqrt(f) [theta_E,+ - f_t/f - f_s/(A f)],
    theta_J,- = sqrt(f) [theta_E,- - f_t/f + f_s/(A f)].

It would be wrong to assume that f>0 makes trapping conformally invariant.
The derivative terms matter. This agrees with the general frame distinction
in [Faraoni and Nielsen (2011)](https://arxiv.org/abs/1103.2089).
The implementation independently checks the area-rate identity
`theta=2(nR +/- e_s R)/R`, differentiating R_E=B and R_J=B/sqrt(f) directly
with a sixth-order stencil, against the above formula with the evolution's
fourth-order stencil. This is independent spatial evaluation of the same
saved fields, not an independent evolution.

## Measured spheres

At the neck, reflection symmetry makes the two expansions equal, up to
numerical differentiation error. Values below are for N=2048; each entry
represents both theta_+ and theta_-.

| Preparation | Time | Einstein expansions | Jordan expansions | Classification |
|---|---:|---:|---:|---|
| eta=0 | 0 | 0 within 1.2e-13 | 0 within 5.9e-14 | both marginal |
| eta=.3 | 0 | +0.001498088104 | +0.001401333105 | past trapped |
| eta=0 | .01 | -11.526001851 | -10.786309276 | future trapped |
| eta=.3 | .01 | -11.524300028 | -10.784706055 | future trapped |

All six archived runs (both preparations, N=512,1024,2048) agree on the
final neck's sign. Medium/fine differences in its expansions are at most
2.27e-9 and decrease by approximately 16 on refinement. The finest
whole-grid difference between the two area-expansion routes is below
1.13e-9. These are empirical convergence checks, not rigorous error bounds.
The classification tolerance is 1e-8, used only to avoid classifying numerical
zero as strictly trapped; it is not an event threshold or a quantum rule.

The Einstein frame has 121 future-trapped grid sections at N=2048; the
Jordan frame has 122, for either preparation. These are sampled round
spheres, **not numbers of objects or horizons**. In particular, the central
bulk sphere has positive Einstein expansions (~+0.012795) and negative
Jordan expansions (~-0.0006097). That extra Jordan sphere is associated
with the small bulk contraction; it is not evidence of a second localized
mouth. Angular dependence of these geometric diagnostics vanishes exactly
in the spherical metric ansatz.

The data therefore demonstrate trapped surfaces, not two isolated black
holes, an outermost apparent horizon, an event horizon, or a complete
Schwarzschild/Kruskal exterior. Such global identifications are not licensed
by this short-time local audit. The null energy condition of the regular
Einstein-frame scalar system does not prevent the existence of these
trapped surfaces or gravitational interaction through the exterior.

The archive retains full initial/final fields but not intermediate full
field snapshots. Its 201 diagnostic rows do not retain l, f_t and B_s.
Consequently this audit does **not** measure the trapping-onset time,
trapped-tube propagation speed, or expansions at each probe crossing.
Those stronger claims in the [#308 review](https://github.com/davidmdrpi/geometrodynamics/pull/308#issuecomment-5841520987)
come from the reviewer's separate evolution. The endpoint check corroborates
final trapping without adopting those additional measurements as our own.

## What P2=-P1 does and does not mean

There are three different constructions to keep separate.

| Construction | Surfaces and observable | What is established |
|---|---|---|
| #303 review's excised S3 complement | Two distinct boundary spheres; conformal-Killing charge integrals of symmetrized Bowen-York seed data | Reported initial constraint-charge relation C2=-C1; no Hamiltonian completion or evolution in that review |
| #304 compact mapping torus, inherited by #308 | Opposite oriented faces of one identified neck cut; radial section flux, not mouth linear momentum | Equal-and-opposite flux by smooth gluing and normal orientation |
| Desired exterior interaction experiment | Two separately tracked gravitating regions with independently evaluated surface charges and bulk flux | Not yet implemented or tested |

The provenance for the first row is the
[#303 review](https://github.com/davidmdrpi/geometrodynamics/pull/303#issuecomment-5690368230).
It reports P=(.3,-.7,.5) as a **seed parameter**, C=-P/2 for the nonzero
conformal charges, and C2=-C1. Its rotational Killing charges vanish in
the even seed sector. Calling all of this simply "P2=-P1" loses the
charge definition and conflates a seed label with a measured momentum.
The review reports separate boundary integrations but also imposes antipodal
symmetry: these are not two independently evolving object momenta. The
supporting implementation is not part of the #308 evidence audited here;
its numerical values are attributed to the review, not newly reproduced.

The same review's assertion that a closed handle has no nonzero section
charge is incorrect. An internal nonseparating surface can have a nonzero
flux, as the following exact baseline demonstrates.
For the exact compact vacuum baseline psi=1, C=1/2 and X=partial_s,

    Q_X = integral (K_ij - K gamma_ij) X^i nu^j dA = 8 pi C = 4 pi.

Cutting at that surface produces Q_+=4pi and Q_-=-4pi. The sum vanishes;
either individual flux need not. This has been added as an analytic
regression control, without changing the historical frozen gates.

For the evolved spherical Einstein-frame fields, the same unnormalised
geometric constraint flux with X=partial_s is

    Q_X(s,nu=+e_s) = -8 pi A B^2 l.

At the final eta=.3 neck the two oriented faces give
`-0.00395308850388` and `+0.00395308850388`. The audit deliberately evaluates
the same stored section with opposite normals and labels the results
`independent_surfaces=false`. This cancellation would survive even an
incorrect evolution: it cannot certify recoil, radiation flux, or exchange.
No 1/(8pi G) normalization or asymptotic particle interpretation is implied.

The magnitude is dominated by the common collapse, not by the initial eta
charge carried forward. The same archived finest-grid values give:

| Positive-normal neck flux | Value |
|---|---:|
| Initial eta=.3 (initial eta=0 is zero) | +5.14612585015e-7 |
| Final eta=0 reference | -0.00395366908817 |
| Final eta=.3 | -0.00395308850388 |
| Final difference, eta=.3 minus eta=0 | +5.80584297780e-7 |

The difference is 12.8197% above the initial eta-dependent flux, underneath
a reference contribution about 7,683 times the initial flux. Subtracting
the matched reference isolates a preparation-dependent response **in this
specified slicing and charge convention**. It does not make the response
gauge invariant or turn it into linear momentum. Both raw runs and the
subtraction must be retained, including reference terms in any balance.

There is also a distinction between radial dilation and translation. In
the round-S3 bulk chart, x_c=tanh(s) and
`gamma_round=sech(s)^2 (ds^2+dOmega^2)`, hence
`grad_round(x_c)=partial_s`. Near the two poles this acts as `+/-r partial_r`:
a dilation, with opposite signs. It remains a conformal-Killing field for
the initial conformally round metric, but its round-gradient expression is
not the gradient in the deformed physical metric. On #304's product
baseline it is instead the longitudinal Killing field of the mapping
torus; the two geometric interpretations should not be conflated.

For #306's spherical maximal initial data, A=B=psi^2,
`k=2a/B^3, l=-a/B^3`, and therefore `Q_X=8pi a(s)`. Its radial variation is
fixed by the scalar current through the momentum constraint; its constant
part is free initial data. It is a radial dilation-type constraint charge,
not a measured translational momentum. Charges transforming as vectors
under the preserved SO(3), including the transverse round-S3 embedding
coordinate gradients, vanish by symmetry on these round sections. Directed
recoil requires the nonspherical degrees of freedom in the next design.

More generally, write pi^ij=K^ij-K gamma^ij and use the convention
D_i pi^ij=j^j. For a chosen comparison vector X,

    D_i(pi^i_j X^j) = j_j X^j + pi^ij D_(i X_j).

The integrated identity includes the deformation term unless X is Killing.
If X is conformal Killing, `D_(i X_j)=lambda gamma_ij`, that term is
`-2 K lambda` in three spatial dimensions. Its usual trace-free reduction
requires **maximal slicing K=0** (or a separately justified vanishing
term); it cannot be assumed throughout #308's geodesic-slicing evolution.
At the final eta=.3 seam, K=k+2l=0.59443226538. Moreover X=partial_s need
not remain conformal Killing once the evolved A and B differ. In these
coordinates its full deformation term is

    pi^ij D_(i X_j) = -2 l A_s/A - 2(k+l) B_s/B.

Only when A_s/A=B_s/B=lambda does it reduce to -2K lambda. At the symmetric
seam the spatial derivatives vanish, but that pointwise fact does not
remove deformation terms throughout a volume. This is a
**spatial constraint identity**, not a time-evolution momentum-transfer law.
On a compact curved universe, comparing charges at different locations
also requires specifying X or a frame-transport prescription. Coordinate
components from different local frames cannot simply be added.

## Consequence for the next experiment

The next target is **two gravitating regions exchanging momentum through
the exterior bulk**. The positive result above licenses keeping trapped
mouths in the candidate model; it does not certify a two-object milestone.
The implementation sequence and required observables are in
[the exterior-bulk test design](exterior_bulk_momentum_design.md).

The #308 matter balance remains a numerical consistency check of the
coupled field evolution. Its test probes have no stress or backreaction.
At the symmetric seam A_s=f_s=0, so dp_s/dt=0 there while
`d p_hat/dt = p_hat (k + f_t/(2f))`. The leading local-frame window change
is gravitational redshift from the changing radial frame, rather than
measured recoil of another gravitating body.

`FINITE_CROSSING_IMPULSE=false` remains a measured negative for the tested
smooth crossing. The historical `DISCRETE_RECIPROCAL_MOMENTUM_EXCHANGE=false`
remains unchanged in #308's schema, but its scientific status is
**NOT_TESTED** because no independent gravitating receiver exists there.
The new audit uses that explicit status. Trapping, classical reciprocal
exchange, and action discreteness remain separate questions.

## Reproduction

From the repository root with project and test dependencies installed:

```sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.handle_surface_audit
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_handle_surfaces.py
```

The audit restores and hash-checks existing evidence and samples the archived
initial polynomial; it never calls the evolution runner. Tests prohibit
that call explicitly. Tests also cover flat/untrapped, contracting/trapped,
expanding/past-trapped and marginal controls; time and normal reversal;
conformal-frame disagreement; an independent area derivative; distinct
versus same-cut fluxes; malformed input; tampered archives; and withdrawal
of stale output after failure.

## Review follow-up, 2026-09-26 UTC

The [independent #309 review](https://github.com/davidmdrpi/geometrodynamics/pull/309#issuecomment-5843046495)
reproduced the archived audit and independently confirmed its expansions.
The reference-flux comparison above uses the already archived audit values;
no evolution or revised audit was generated. The follow-up clarifies
collapse versus preparation response, dilation versus translation, and
the maximal-slicing condition. The companion design now includes the
one-handle neck contribution, causal-route requirements and the residual
SO(2) symmetry of the specified antipodal encounter. Local trapping alone
is not promoted to a proof that all communication must use the exterior.
