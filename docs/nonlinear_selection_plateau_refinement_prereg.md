# Prospective window-quadrature refinement after the original gate failure

Original freeze: aa7ff3d19d61c8555502eb34cb0021aa13411698.
Publish this extension before its new integrations. No selection threshold,
phase, amplitude, readout or observation window changes. The original verdict
INCONCLUSIVE_NUMERICAL_FAILURE must remain visible and the original raw archive
and diagnostics must be preserved.

The completed original scan has all 108 preparations and 15 independent
controls. One preparation, epsilon=.34171875,theta=0,phi=pi/2, exceeds the
1e-6 J_bg quadrature gate: the difference is 1.0437832688897085e-6 (DOP853)
and 1.0437832691596825e-6 (RK45). All other numerical gates pass. All five
candidate plateaus fail all their selection gates before refinement. This
motivates a numerical accuracy check, not a new selection criterion.

Repeat only this preparation with BOTH registered integrators and unchanged
tolerances, saving t=0 and .025-spaced samples in the same four late windows.
Compare .025 vs .05 window quadrature, require the original <1e-6 J_bg limit.
Also compare refined means with the original means, and check constraints,
chart, background/work checks and independent-integrator agreement with the
original tolerances. Non-refined cases retain their original passing gates.

Build a separately labeled refined analysis by replacing only the two affected
window-mean rows and their error diagnostics. Re-evaluate the same five
candidate amplitude windows and unchanged selection criteria. The original
record and verdict are not overwritten. If the refined numerical gates pass,
report the mechanism-specific candidate or negative verdict separately; if
they fail, retain INCONCLUSIVE_NUMERICAL_FAILURE. Archive the refined states,
source hashes and a replay validating which two rows were replaced.
