# Review follow-up: separate supporting-bulk and shell failures

The [independent review of #333 at 52f3c83](https://github.com/davidmdrpi/geometrodynamics/pull/333#issuecomment-6158675664)
confirms the result and requests a better reference test. **The shell-free
supported background is already unstable.** Its conformal kinetic direction
must not be counted as a new defect introduced at the focus.

The reviewer reports 97 passing tests (NumPy 2.4.6, SciPy 1.17.1), a rerun
identical to the original archive except environment, and green Python
3.10/3.12 CI. Their independent full-tensor symbolic computation verifies the
bulk reduction. Their separate shooting and norm integrations use this PR's
junction conditions and norm formulas; they explicitly did **not** independently
rederive those boundary conditions or the symplectic form. The latter retain
this PR's internal full-action, tensor, gauge and finite-element checks.

## Corrected admission statement

Compare the new construction with the specified shell-free reference background:

| Diagnostic at C=.76 | No shell, regular S4 | Identified shell at a=.3 |
|---|---:|---:|
| Negative even kinetic directions | 1, conformal | 2: bulk block plus shell term |
| Growing physical radial modes | 1, homogeneous | 3 found in the computed spectrum |
| Real-frequency negative-norm modes | 0 | 1 found, the shell oscillation |

The **incremental shell defects** are one additional negative kinetic direction,
a negative-energy real-frequency mode, and two additional growing modes in the
computed benchmark spectrum. The structural increment is the term
3 D0 r0^4 z^2 < 0. The branch-wide statement is this extra kinetic direction;
the listed eigenvalues and growing-mode counts refer to the benchmark.

A useful focus-admission test asks whether the focus construction adds negative
kinetic directions or new instabilities relative to the same support. It does
not demand positive sign for a background conformal direction while attributing
that direction to the focus. Separately, the **complete static model remains
unstable**: adjusting the comparison does not make its Eddington mode acceptable
or remove the need to address support dynamics.

This is a post-review correction of attribution and of the focus-specific
admission criterion, not a retroactive change to the original protocol. The
original `COUPLED_RADIAL_ADMISSION_FAILED` label stands under the relative test.
Its protocol, producer and source-bound archive remain unchanged. A separate
[review result](../experiments/closure_ledger/runs/20261014_egb_coupled_review/result.json)
records `failure_reasons.background` and `failure_reasons.shell` and hashes the
original archive. It supersedes the original flat failure list for attribution.

## No-shell calculation

Remove both excisions and the surface action. Retain the same R, alpha, q and V
of the canonical bulk completion. This matches the local bulk background and
action parameters, not the total charge: restoring the caps and removing the
charged surface generally changes Q. It is a reference solution, not a
charge-conserving shell-removal trajectory. Regular scalar radial eigenfunctions on the
whole sphere are C_l^(3/2)(cos u), satisfying

\[
-\Delta_u f_l=l(l+3)f_l,\qquad
\Lambda_l=l(l+3)-4-2/C,\qquad h_l=l(l+3)-4.
\]

The homogeneous l=0 mode grows for **every C>0**, with
Lambda0=-4-2/C. The l=1 solution is the exceptional coordinate mode, not a
second physical instability. At C=.76:

| l | Lambda | Status |
|---|---:|---|
| 0 | -6.63157894737 | Homogeneous Eddington growth |
| 1 | -2.63157894737 | Gauge; zero norm |
| 2 | +3.36842105263 | Positive-norm oscillation |
| 3 | +11.3684210526 | Positive-norm oscillation |

For **unit full-sphere weighted L2 norm**, the physical kinetic coefficient is

\[
G_l=3CR^4\int_0^\pi\sin^3u[(f'_l)^2-4f_l^2]du
=3CR^4[l(l+3)-4].
\]

There is precisely one negative direction, l=0, and the l=1 null direction is
gauge. The coefficient is **3C on the full sphere**: #333's coefficient 6C
already accounts for two halves when integrating only [a,pi/2]. Using 6C on
[0,pi] would double the norm. This normalization clarification changes no signs
or frequencies. At R=1,C=.76, G0=-9.12 under the full-sphere unit norm.
The positive common Omega3/kappa5 is omitted, as in the original report.

The new probe verifies the Gegenbauer norms by quadrature and the differential
equation by polynomial derivatives. Exactly one growing physical mode is a
**benchmark statement**: l>=2 also grows whenever l(l+3)<4+2/C. In particular
l=2 is unstable for C<1/3 and marginal at C=1/3. The kinetic signature does
not change at that potential-instability threshold.

The even-lowest shell mode approaches the no-shell value along the checked
small-mouth family: Lambda=-6.37353785 at a=.3, -6.60184498 at .1 and
-6.63038020 at .02, compared with the limiting -6.63157895. This supports its
identification with the perturbed homogeneous bulk instability, independently
of simply observing a negative norm.

## Odd mode: report its density, not just its growth rate

The odd physical branch satisfies Dirichlet conditions at the mouth and equator.
At every sampled nonzero a it has h>0 and a nonzero gauge-invariant density
perturbation. It is not the exceptional h=0 coordinate mode.
As a decreases, its eigenvalue approaches the no-shell l=1 value. Domain
monotonicity bounds its spatial eigenvalue below by 4. A cutoff of cos(u) near
the shrinking pole has vanishing additional weighted gradient cost in four
spatial dimensions, providing the matching upper limit. Thus h tends to zero.
This is a singular boundary limit, not permission to discard the finite-mouth
mode as gauge.

With the original odd normalization integral_a^(pi/2) sin^3(u) f^2 du=1,

\[
\delta\bar\rho=q^2h f Q(t),\qquad
\|\delta\rho/\rho_0\|_{L^2_w}/|Q|=\frac{C}{C+1}h,
\qquad G=6CR^4h.
\]

The background density is spacetime constant, so its first variation is gauge
invariant. The density measure below is per unit shape amplitude Q, not per
unit energy; the latter would change normalization as G tends to zero.

| a | Lambda | h | Density contrast L2 / abs(Q) | Physical norm |
|---|---:|---:|---:|---:|
| .30 | -1.41649415648 | 1.215084791 | .5246957052 | 5.540786646 |
| .20 | -2.07849673388 | .5530822135 | .2388309558 | 2.522054894 |
| .10 | -2.48720360551 | .1443753419 | .06234389762 | .6583515589 |
| .05 | -2.59462018060 | .03695876677 | .01595946747 | .1685319765 |
| .02 | -2.62559942347 | .005979523895 | .002582067137 | .02726662896 |

The growth rate tends to sqrt(2/C)/R even as these density and norm measures
vanish. This does not establish uniform smallness of every derivative near the
shrinking boundary, nonlinear regularity, or absence of a finite-a instability.
The scan fixes R=1, alpha=-.06 and eta=.1. Each a has its own equilibrium-anchored
surface action; it is a comparison of static completions, **not a shrinking-shell
evolution with one fixed surface material**.

## Program consequence and reproduction

Two questions must now be separated in the ledger:

- **Bulk support:** this rolling scalar and constant potential do not give a
  stable static S4 background, even without a focus. Earlier unspecified support
  sectors are not automatically covered by this result. A future completion
  must change or stabilize the support, or demonstrate how its coupled boundary
  dynamics stabilizes the global mode; the no-shell calculation is not a theorem
  that every conceivable boundary modification leaves it unstable.
- **Focus construction:** this particular EGB-WEC shell adds a negative kinetic
  direction and benchmark growing modes relative to that background.

Both use the closed-system time-plus-scalar-shift generator H-qJ, not an ADM
energy at infinity. The conformal reference signature is not being promoted to
a new quantum ghost claim or a proof of nonlinear vacuum decay.

```bash
python -m experiments.closure_ledger.egb_coupled_background_review
python -m pytest tests/test_egb_coupled_background_review.py -q
```

The review controls are post hoc and do not alter any original numerical gate
or archive. Four new tests cover the regular-sphere spectrum/norm, the C-dependent
instability count, odd density/norm limits, and split-result provenance.

Validation: **101 focused tests pass** under NumPy 2.3.5 and 2.2.6, SciPy
1.17/Python 3.12.14: the preceding 97 plus these four review-control tests.
