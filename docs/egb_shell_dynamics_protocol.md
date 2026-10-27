# Moving EGB shell: frozen numerical protocol

This protocol is written after deriving the radial stability formula but before
running the scan and trajectories. It is a prospective numerical check, not a
blind discovery or a preregistration of the analytic result.

## Model and fine tuning

Use a five-dimensional ultrastatic round exterior with radius R, remove two
antipodal balls, and identify their S3 mouths directly. There is no finite smooth
neck. The action is Einstein--Gauss--Bonnet with constant alpha; the supporting
bulk stress is retained. This changes gravity away from the focus as well.
The distributional shell is a weak-junction regularity candidate, not bounded
pointwise curvature and not a topology-creation solution.

At a=0.3, r0=R sin(a), use p=p0+eta(s-s0), where s=kappa5 sigma and
p=kappa5 pressure. Each intercept is explicitly anchored to its equilibrium.
This is an imposed local material response, not a derived microphysical action,
an attractor, or selection of the mouth size. Let eta range from 0 to 1 as a
material-response proxy. Include the nonzero exterior flux in shell conservation.
The radial ansatz keeps the exterior geometry and supporting stress fixed;
it does not evolve the coupled bulk perturbations. No corrective kicks: budget 0.

## Checks and failure gates

1. Compare the reduced stress to a full Davis junction tensor contraction at 32
   seeded, moving, accelerating states: relative error <1e-10. Check the GR limit.
2. At 32 seeded states, check the sourced shell conservation identity by numerical
   differentiation: normalized residual <1e-7. Omitting the flux must fail at a
   moving reference state. A source-free shell must not be silently substituted.
3. Scan 101 couplings strictly inside (-R^2/4, alpha_WEC) and 101 eta values
   in [0,1]. Compare the acceleration Jacobian to the analytic omega^2 with
   normalized error <1e-6. Record the full grid. Neutral omega^2=0 is not stable.
4. Preselect alpha/R^2=-0.06, eta=0.10 as stable, and the same coupling with
   eta=0.30 as unstable. For stable motion use fractional radial displacements
   1e-3 and 1e-4, initially at rest, for 20 linear periods. Compare a tighter
   integration at 1e-3. Coordinate norm is
   sqrt(((r-r0)/r0)^2+(v/(r0*omega))^2); no coordinate invariant is implied.
   Require amplitude norm <1.1 times its initial value, period error <1e-3,
   decreasing period error at smaller amplitude, and replay error <1e-6 in this
   norm. The horizon is fixed by 20 oscillations, not a claim of indefinite
   nonlinear stability; publish the measured nonlinear frequency shift and
   cumulative phase drift. No twist-based long-time stability claim is made.
5. For instability use initial displacement 1e-5 r0 through five linear growth
   times, comparing to cosh growth. Require >50-fold growth and <2% relative
   discrepancy in displacement. Do not infer the eventual endpoint from it.
6. Integrate pressure work and exterior flux as independent ledger variables.
   Require normalized energy-plus-work-minus-flux residual <1e-8. Stable samples
   must retain shell WEC and DEC; publish minimum margins. Stop before a pole,
   equator, or acceleration-coefficient zero; any stop is unresolved, not a bounce.

If verification fails, label NUMERICALLY_UNRESOLVED. If the checks work but no
open WEC radial-stability region exists, label NO_RADIAL_WINDOW. A verified open
window is only OPEN_RADIAL_WINDOW_COUPLED_ADMISSION_OPEN. It is not an admitted
focus law. Bulk matter realization, the canonical kinetic norm, nonspherical
and coupled characteristics, finite-thickness regularity, and event activation
remain uncomputed regardless of the radial result. In particular the negative
pressure--acceleration coefficient is disclosed, not reclassified as healthy.

## Secondary scale screen

Report a one-scale Casimir/Einstein balance in native D=5 and its loop-length
dependence. State the assumptions and distinguish it from the magnetic,
near-extremal Maldacena--Milekhin--Popov mechanism. No quantum support is claimed
to derive quantum mechanics, and no D=4 Planck scale is imported into D=5.
