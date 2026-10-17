# First-transit exchange test

Frozen before this branch's numerical runs. This is an independent implementation
and replication of PR340's known reviewer results, not a blind discovery test.
The additional lambda=1/8 point and fixed-axis reference are prospective checks.

Use the unchanged PR340 Hamiltonian, total pair inertia .1, common two-packet
normalization, eta=.7 and 0, cutoff32, dt=.002. Scale prepared a,p together by
lambda=1,1/2,1/4,1/8. Run both packets and each separately, identical mouth initial
states. End at pi+1.2. Record exact .6, pi-.6, pi, pi+.6, and end; monitor full
energy and angular momentum at every integration step. Stop on the first
crossing of minimum cross-mouth angular distance .65, retaining the last valid
state and bracket. A stop makes this first-transit comparison incomplete.

Define R_h(lambda)=[L_h(both)-L_h(x)-L_h(z)]/lambda^2 and
C_h(lambda)=[4 R_h(lambda/2)-R_h(lambda)]/3. Compare the last two C estimates:
1% relative tolerance, absolute floor 1e-8, for mixed y. Leading reciprocity
requires |C_1y+C_2y| below max(.01 max|C_hy|,1e-8).

Independently solve the fixed-axis linear field/shape equations with adaptive
DOP853 (rtol=2e-10, atol=2e-12), integrating the explicit bilinear mouth torque
kappa[(q_x-U_x) grad U_z+(q_z-U_z) grad U_x]. Include reciprocal field forcing;
do not infer this reference by subtracting nonlinear moving trajectories.
Compare reference and final extrapolation with the same 1%/1e-8 budget.
Repeat the smallest two amplitudes at dt=.001 and cutoff40 separately, and the
reference at cutoff40. Require changes below the same budget. Reference torque
quadrature must satisfy the mixed field-plus-pair angular ledger within 1e-8.

For mixed retained field Jy, halving lambda should give ratio16 within20% for
the last two adjacent amplitude ratios, provided the signal exceeds numerical
uncertainty (step/cutoff changes and ledger error). Also record all vector
components and the independent-trajectory composite separation. Its smallest
two adjacent ratios should be quartic by the same criterion; no reciprocity
claim is made about the full recoil vector. Energy error must stay below1e-4
and angular ledger error below1e-8 (normalized as in PR340).

Failure: a resolved nonzero quadratic sum, disagreement with the independently
integrated bilinear torque, or a resolved non-quartic retained field defeats
the proposed leading reciprocal-exchange interpretation for this preparation.
Unconverged/underresolved results are unresolved, not success. Report all gates.

This tests classical rotational exchange in an assumed round-S2 scalar/shape
Hamiltonian. It does not test neck transport, GR, quantization, a measurement
law, or a long-time attractor. The source half-focus and target full-focus are
reported separately, without attributing the whole endpoint impulse to the
antipodal event.
