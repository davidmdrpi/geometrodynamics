# Attacking the fine-tuning: the electron near-zero — stabilized or dialed? (PR #194)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The question #192/#193 left sharp

#192 found the locked surrogate's electron eigenvalue is a **near-zero**
(0.1996 in action units vs μ 41.26, τ 694.98) sitting 1.4% of fiber squash
from a zero crossing. #193 showed the near-zero has no counterpart in the
actual wave operator (E₁ = 2 + 1/λ² ≥ 2): it enters only with the
instanton dynamics. The remaining question: is the near-cancellation
**stabilized** by anything — a symmetry, an index, a seesaw structure, an
attractor — or **genuinely dialed**?

Either outcome is a real result: a mechanism would be a discovery; an
exclusion sharpens the hierarchy problem the surrogate carries.

## The anatomy

The near-zero is a cancellation between two a-priori unrelated quantities:

```
h₁₁ = 6.8754  (the k=1 diagonal: 2π base action + resistance)
− repulsion = 6.6758  (transport coupling to k=3,5: t₁₃ = 14.85, t₁₅ = 10.46)
= E_e = 0.1996   (a 2.9% residue)
```

The zero locus `det H = 0` is codimension 1 in the 7-parameter space — a
zero eigenvalue needs no symmetry, only one tuned combination.

## Candidates tested and excluded

| candidate | test | result |
|---|---|---|
| chiral / sublattice | `SHS = −H` for all 8 sign matrices | none exists (needs zero diagonal; tr H = 736) |
| spectral reflection | eigenvalues symmetric about 0 | absent |
| index-like zero-mode | overlap of the near-zero eigenvector with structured vectors | best 0.917 < 0.999 — pinned to nothing |
| seesaw suppression | sign stability (multiplicative smallness is sign-stable) | E_e **flips sign** under ±2% transport (+0.43 → +0.20 → −0.03): cancellation type |

## The tuning, quantified and localized

Barbieri–Giudice sensitivities Δ_p = |d ln E/d ln p|:

| parameter | Δ(E_e) | Δ(E_μ) | Δ(E_τ) |
|---|---:|---:|---:|
| transport_strength | **−57.1** | 0.26 | 0.00 |
| action_base | **+31.5** | 0.15 | 0.01 |
| action_slope | **+30.9** | −0.14 | −0.00 |
| hard_pinhole_gamma | **+17.9** | 0.46 | 0.03 |
| resistance_scale | +7.4 | 0.13 | 0.05 |
| k_uplift_beta | +1.2 | 0.00 | 0.90 |

**Only the near-zero is tuned** — every μ- and τ-level sensitivity is
below 1. The heavy rungs are natural outputs of the geometry-scale
parameters; the whole fine-tuning is concentrated in the electron's
cancellation (the surrogate-side fingerprint of the #192 result).

The one dialed combination (codimension 1): n̂ ∝ ∇ln E_e is **76%
transport** opposing **42% base action + 41% slope + 24% pinhole**, with
global tuning Δ = |∇ln E_e| = **74.7** — a 1% move along n̂ changes E_e by
~75%.

**Cross-check against #192:** contracting ∇ln E_e with the fiber-map
direction (base action + slope + resistance + uplift, each ×λ) gives
**+71.1** — exactly the Berger λ-sensitivity #192 measured. The Berger
squash found the fine-tuning because it has an 87% overlap with the
dialed direction: the two probes see one and the same dial.

## The Monte Carlo null

20 000 samples, log-uniform ±25% around the locked point (fixed seed):

- `P(|E_e| ≤ 0.1996) = 0.077` vs the naive codimension-1 linear-measure
  estimate `2/(Δ·2ln1.25) = 0.060` — same order;
- the E_e histogram is **flat through zero**: no accumulation at 0 (which
  an attractor/mechanism would produce), no gap (which a protective
  repulsion would);
- `P(μ/e ≥ 206.7) = 0.040` — the observed hierarchy is reached only
  inside the tuned sliver.

The near-zero is exactly as rare as generic linear measure predicts — a
~7% accident under a ±25% prior, with nothing enhancing or suppressing it.

## Origin: the calibration imports the hierarchy

The surrogate is calibrated to the observed masses. With the matrix's
natural scale O(10), demanding μ/e = 206.77 **forces** |E_e| = E_μ/206.77
= 0.1996 — the near-zero is not an internal accident, it is the observed
hierarchy transferred from the data to the parameters by the fit. The
surrogate *carries* the hierarchy problem; it does not solve it.

**Numerology guardrail** (the #165 anti-rigging rule): the det-zero roots
were measured against round constants and **not** matched — e.g. the
transport root 25.538 vs 8π = 25.133 differs by 1.6% (reported, rejected
as a derivation).

## What a real solution would require

New structure that pins a k=1 zero mode: an index, a chiral grading of
the winding sectors, or a geometric identity tying the 2π base action to
the transport repulsion. Finding or refuting such structure in the throat
geometry is the follow-up.

## Honest scope

- The attack is on the locked surrogate (where #192 measured the
  fine-tuning); the #193 operator enters only as the contrast (no
  near-zero to tune).
- The Monte Carlo prior (log-uniform ±25%) is a declared choice; the
  linear-measure agreement is prior-consistent, not prior-independent.
- Exclusions cover the standard protection mechanisms expressible at the
  3×3 matrix level; a mechanism living below the surrogate's resolution
  (in the throat field theory) is exactly what the follow-up must probe.

## Reproduce

```bash
python -m experiments.closure_ledger.electron_near_zero_naturalness_probe
# Verdict: ELECTRON_NEAR_ZERO_IS_UNPROTECTED_CANCELLATION_TUNING_IMPORTED
#          _BY_CALIBRATION_NO_STABILIZING_SYMMETRY_FOUND
```
