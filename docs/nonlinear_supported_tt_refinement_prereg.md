# Prospective amplitude-refinement extension

Date: 2026-09-13. Original scientific freeze: `a067782`.
This extension is published before any epsilon=.005 or .0025 nonlinear
trajectory is computed. It preserves the original run and its failed gates.
It does not change nonlinear_supported_tt_prereg.md.

## Observed original outcome

The original full run passes 11 of 13 gates. At the frozen finest amplitude
.01, the largest first-variation error is 0.0011758078517276156 and the
largest second-variation error is 0.0010988125538551976, against 0.001.
Thus linear_recovery and quadratic_response fail. D and C pass; N and F
are UNRESOLVED under the original dependency table, even though the exact
continuation certificate has positive margins. No threshold is relaxed.

The source and raw evidence at this decision are identified below by SHA256.
They must accompany the eventual results. Their values identify the original
run, rather than allowing the extension to be mistaken for its prehistory.

## New experiment, frozen before running

For all fifteen original U,V pairs use the original phase set: {0,pi/4,pi/2}
for the seven explicit pairs and {pi/4} for the eight random pairs. Use
matched amplitudes +/-.01, +/-.005, +/-.0025 and zero. Reuse already computed
identical trajectories where available. Use the original proper-time sample
set, integrator, tolerances, independent variational equations, coordinate
normalizations, and error norm. Do not refit a coefficient or change the
constraint-completion prescription.

The hypothesis is second-order truncation of central amplitude differences.
For both variations require finest-amplitude relative error <.001 through
t/a=2 for every pair/phase. For each pair/phase, use its maximum error over
0<t/a<=2 at each amplitude; consecutive ratios should lie in [3.5,4.5]
where both errors exceed 1e-8. Retain all ratios and floor exclusions. If
the accuracy or convergence gates fail, the refinement is UNRESOLVED.

Keep the original report and verdict unchanged. Emit a separately named
refinement report and verdict with this extension's published commit SHA.
In that verdict only, replace the two amplitude-difference gates with the
new accuracy/convergence gates; retain every other original requirement,
including independent linear-operator recovery and the future-continuation
proof. No numerical tail can substitute for that proof. The original
11/13 outcome must remain visible in the results document and PR body.

The extension supplies a sharper derivative test, not a proof that every
finite amplitude in the original grid is within its asymptotic regime.
Neither outcome permits reclassification of a failed construction as a
nonexistence theorem.

## Original source/evidence hashes

- `geometrodynamics/waves/nonlinear_supported_tt.py`: `82d52a348fd9df691a44af7edbb2596f3cccf7f78de84183e9f46ce23f58152c`
- `experiments/closure_ledger/nonlinear_supported_tt_probe.py`: `42ef10dee9f82b58a0e56ac310a0aa57a7c2e3d95921b9b403aed1284f2831f4`
- `experiments/closure_ledger/runs/20260913_nonlinear_supported_tt/probe.json`: `15da9598010dfe92aff6a16932256ccdcca8a0a2849524c8f57598922bcb0581`
