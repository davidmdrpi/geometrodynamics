# Boundary-pairing baseline and next experiment

**The existing supported Maxwell throat passes the boundary-joining
qualifications. The responding conjugate-mouth experiment is preregistered,
not yet implemented.**

[Protocol and prior-PR review](conjugate_boundary_pairing_protocol.md) were
published at [`836bb278`](https://github.com/davidmdrpi/geometrodynamics/commit/836bb278074030262819c4bbfd91c5b791068be1)
before the qualification code or new numerical output. Baseline main is
`3293d5bad3c8e103c8514ad85afc0128050004e4` (PR #350).

The source, support, readout and composition distinctions from PRs #279–#288
and the Maxwell/support findings of #345/#349 change this experiment's design.
In particular, the preceding discussion's antisymmetric spinor pairing remains
a conditional benchmark. It has not been extracted from the Maxwell action.

## Qualification measurements

The probe covers lapse depths 0, 0.4 and 0.8 at frequencies 20 and 60. Each
bridge is split once at its waist and once at two off-centre surfaces, producing
12 joining comparisons. Whole and segmented transfer integrations do not
impose the transfer determinant or DtN symmetry. Shared boundary fields are
eliminated by outward-flux continuity. Two independent load comparisons retain
complex phase information. Errors use `norm(delta)/max(1,norm(reference))`
where applicable; they are not relative errors in tiny transmission fractions.

| Check | Maximum error | Registered limit |
| --- | ---: | ---: |
| Whole versus joined DtN | 7.9764e-12 | 2e-7 |
| DtN reciprocity | 2.0579e-11 | 2e-7 |
| Whole map versus inherited DtN | 1.4551e-11 | 2e-7 |
| DtN-derived versus outgoing-solution scattering | 1.0545e-10 | 2e-7 |
| Three-piece joining order | 1.0407e-15 | 2e-7 |
| Coherent-field load versus polarized load matrix | 3.9449e-16 | 2e-6 |
| Internal load versus scattering-phase derivative | 4.0255e-11 | 2e-6 |
| Fixed-background third-order interference | 2.5251e-16 | 1e-12 |
| Conditional spinor identity and transport controls | 8.8818e-16 | 1e-12 |

All **12 mandatory qualification gates pass**. The positive outgoing-power
Gram and its vanishing third-order interference are structural calibration
identities. They do not establish a physical history measure.

Two counterexamples remain visible:

- Four individually positive, normalized context-specific Gram matrices have
  unbiased marginals and **CHSH=4**. No common joint history functional is
  supplied. Thus separate detector-configuration Gram tests cannot establish
  the Tsirelson bound.
- Unequal paired channel efficiencies change the conditional spinor
  correlation at angle 1 by **0.2878701**, despite the invariant epsilon
  identity. A global Pin sign cancels from the squared pairing. Physical
  channel response and a relative-path sign experiment remain necessary.

## What the next implementation must discriminate

The first action expansion predicts a useful separation. In the specified
effective support comparator, the leading interaction is `q psi^2`. With an
unpumped equilibrium and zero initial mechanical displacement, mouth response
starts at field amplitude squared. The field's first backreaction correction
is cubic and a quadratic field readout receives quartic corrections, at fixed
finite time in the smooth regime. A zero-field linearization therefore cannot
test whether dynamic geometry selects the readout.

The protocol next calls for a solved responding background, the full coupled
Jacobi response, and coherent cut-and-rejoin tests retaining mechanical work.
The retarded response must be distinguished from a symmetric on-shell-action
Hessian and from a fluctuation/noise kernel. Any actual spinor boundary map
must come from matched physical degrees of freedom, not a relabeling of the
two endpoint indices. A stable detector transition and its preparation are
separate requirements for an event-frequency interpretation.

Current verdict: `BASELINE_QUALIFIED_DYNAMIC_EXPERIMENT_PENDING`.
Dynamic mouth response is `NOT_IMPLEMENTED`; physical spinor matching,
detector-frequency law and common all-setting history kernel are `NOT_DERIVED`.
This is a qualified starting point for the next calculation, not a no-go for
the proposed geometric mechanism or a claim to have derived quantum theory.

## Reproduce

```bash
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.boundary_pairing_baseline_probe \
  --output /tmp/boundary_pairing_baseline.json
OPENBLAS_NUM_THREADS=1 python -m pytest -q \
  tests/test_boundary_pairing.py tests/test_maxwell_throat_response.py \
  tests/test_mouth_spin_frame.py tests/test_derivations.py
```

The targeted suite passes **66 tests**. It includes analytic free-interval
oracles, deformed/charged bridge cuts, an intentionally wrong interface sign,
singular-coordinate rejection, an imaginary cross-term sign check, nonlinear
readout control and missing/nonfinite-gate and CLI failure-path tests. The
full repository suite was not rerun; no existing solver is modified.

[Machine-readable probe](../experiments/closure_ledger/runs/20261108_boundary_pairing_baseline/probe.json)
and [source/evidence checksums](../experiments/closure_ledger/runs/20261108_boundary_pairing_baseline/manifest.json)
retain the unrounded values and exact source identities.
