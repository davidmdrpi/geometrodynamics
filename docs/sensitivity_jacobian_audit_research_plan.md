# The sensitivity audit: Jacobian rank, the forced core, the isolation dimension (PR #173)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The inverse problem

What remains open is no longer static topology or a static equation of
state. This probe runs the **dynamical inverse problem**: let the continuous
geometry vary and **measure** how the observables respond, via the Jacobian

```
J_ij = ∂(observable_i) / ∂(input_j)
```

evaluated at the lock. Its singular-value decomposition yields three numbers:

- **rank(J) = the isolation dimension** — the number of independent
  observable directions the free inputs can dial (the *fitted* dimension);
- **the forced core = n_obs − rank(J)** — the observable combinations no
  input can move, *forced* by the rigid structure at zero input cost;
- **the compensator redundancy = n_inputs − rank(J)** — input directions
  that produce no independent observable change (degenerate knobs).

The explicit goal: the largest observable set the rigid core forces at zero
input cost — and how over-complete the parametrization is.

## What is run

**Observables** (the live, repo-reproduced set, 14 total): 4 quark mass
ratios (`s,c,b,t / d`), the 5 CKM magnitudes, the Jarlskog `J` and the
angles `β, γ` (from `LOCKED_QUARK_PARAMS_V4` via `extract_physical_spectrum`
/ `extract_ckm_matrix`), and the 2 charged-lepton mass ratios (`μ/e, τ/e`).

**Inputs**: the free (fitted) continuous knobs — explicitly **not** the
k₅-derived locks (`φ_h = π/k₅`, `χ = k₅(k₅−1)`, `uplift = 1−1/k₅²`,
`action = π`, `winding = max`), which are zero-cost by construction. So the
forced core is what survives variation of every *fitted* knob.

## The measured result (not predetermined)

| sector | observables | free knobs | rank | forced core | redundancy |
|---|---:|---:|---:|---:|---:|
| quark | 12 | 15 | 8 | **4** | 7 |
| lepton | 2 | 5 | 2 | 0 | 3 |
| **total** | **14** | **20** | **10** | **4** | **10** |

The quark singular-value spectrum has a **clean gap**: eight values from
`22.6` down to `1.2e-2`, then a drop of `~5e-6` to numerical zero — so
`rank = 8` is robust.

### The forced core = CKM unitarity

The 4 forced combinations are **entirely CKM** (forced weight on CKM
`> 0.99`, on masses `~1e-16`): they are the **CKM unitarity relations**.
BAM's `V = U₊†U₋` is *exactly* unitary (`‖V†V−I‖ ~ 1e-16`), so the 8 CKM
observables lie on the 4-parameter unitary manifold, forcing exactly
`8 − 4 = 4` relations. This is the largest set the rigid core forces at zero
input cost — a genuine structural prediction (a non-unitary mixing model
would not have it), though it is the *standard* unitarity, not a
BAM-specific numerical relation.

### The masses are fitted

The quark and lepton mass ratios carry **no** weight in the forced core —
the geometric ladder sets *where* they sit, but the knobs span them. There
is **no forced mass relation** in the live observable set.

### The compensator redundancy

`n_inputs − rank = 10` input directions move no observable independently. In
the quark sector the redundancy is **dominated by the mass-preserving
diagonal shifts** (kernel share ~0.8) — the compensators introduced in
#161/#164 to move the CKM at fixed masses. This is exactly the *loose knob /
compensator* structure the program flagged qualitatively (n_part, the
diag-shifts), now **measured**: the v4 quark parametrization is substantially
over-complete (15 knobs spanning an 8-dimensional observable subspace).

### CP at zero cost — the honest test

Adding `φ_h` as a 16th input leaves the rank unchanged (`8 → 8`): the
CP-phase observable direction is **already spanned** by the magnitude
couplings, so deriving `φ_h = π/k₅` saves **0** effective inputs in the
Jacobian sense. Honest reading: *"CP at zero parameters"* is a **counting**
statement (no dedicated CP number is fit to data), not a reduction in the
Jacobian rank — and the over-completeness makes this a weak test, since
essentially any extra knob is redundant here.

## What this is

An **audit**: a measurement of the predictive content, honest where it is
not flattering. The forced core is real (CKM unitarity) but modest; the
masses are fitted; the parametrization is over-complete; and the derived CP
phase is a counting economy, not a Jacobian reduction. The method is general
— it extends to any sector the repo reproduces — and is the natural tool for
letting continuous geometry search the spectrum/signature landscape rather
than adding more static proofs.

## Reproduce

```bash
python -m experiments.closure_ledger.sensitivity_jacobian_audit_probe
# Verdict: JACOBIAN_AUDIT_FORCED_CORE_4_CKM_UNITARITY_MASSES_FITTED_REDUNDANCY_10
```
