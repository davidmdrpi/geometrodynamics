# Correction note to the R3 extension specification (freeze `8584346`)

Date: 2026-09-29. Written after the run, in response to the #317 review at
`2945223`. The freeze text, the frozen probe, the archives and all four
registered labels are unchanged. This note narrows interpretation and adds
validation. It changes no threshold and no decision rule.

## 1. Specification §4 overstated closure

§4 stated that "every rational rotation p/q in the family's range has closed
histories (Poincaré–Birkhoff), at isolated actions", and that closure
actions are dense. **Both claims are withdrawn.**

Poincaré–Birkhoff needs an invariant annulus of an area-preserving map,
with twist conditions on its boundaries. None of that was established for
the four-dimensional return map over the measured range. The centre
manifold, its induced area form and the boundary twist were not checked.

Part D enumerates rationals lying between sampled rotation numbers. It does
not solve P^q(z) = z. Its output is therefore reported as **candidate
closure rationals, periodic histories unverified**.

## 2. Part A's endpoint is a continuation failure

The family "ending" (§3) is the registered stopping rule, applied to
numerical failures:
- Newton or integrator failure at a = .3620;
- then "no real clock velocity" at a Newton trial state for the retry at
  a = .3320.

The retry's initial predictor has q'^2 between 2.114 and 2.244 at every
grid point, and the predictor for the a = .3620 step has q'^2 between 1.907
and 2.140. The failure therefore does not show that section data or
circles cease to exist. NO_TURN_IN_FAMILY holds for the 27 accepted
circles through a = .3044. Beyond that it is unresolved.

## 3. Replay must not trust stored acceptance flags

The frozen `score_a` reads each circle's stored `ok` flag. Corrupting an
accepted circle's residual, or deleting its K, leaves the summary
unchanged. The new `experiments/closure_ledger/r3_extension_replay.py`
closes that gap. It:
- authenticates the three archive files by pinned SHA-256;
- validates required fields, and the shape and finiteness of each K;
- recomputes each acceptance decision from the stored numbers and from K
  (action, Fourier tail, amplitude condition) against the registered
  thresholds;
- checks the ladder and retry schedule;
- only then re-scores with the frozen functions.

With `--full`, it also re-evaluates every accepted circle's invariance
residual on the full system. The frozen probe is not edited, because its
source hash is bound into the archives.

Result: `python -m experiments.closure_ledger.r3_extension_replay --full`
returns VERIFIED, with both labels unchanged. The largest re-evaluated
invariance residual over the 27 accepted circles is 4.9e-12.
