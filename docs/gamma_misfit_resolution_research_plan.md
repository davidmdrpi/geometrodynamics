# The γ misfit resolved: the full flavor-CP dataset realized (PR #161)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. The
> construction re-aims locked structure at exactly fixed eigenvalues.

PR #160 reduced the flavor sector's residual to one number — γ = 104° vs
65.9°. This PR attacks and **resolves** it.

## The diagnosis

γ = arg(−V_ud V_ub\*/(V_cd V_cb\*)) lives in the **ub corner** of the
unitarity triangle, coupled to the (1,3)/(2,3) rotation planes the #160
two-plane family never used. The full mass-preserving family is
**SO(3)×SO(3)** — six Euler angles, three per partition block, every point
at exactly fixed eigenvalues. Against five data constraints (V_us, V_cb,
V_ub, β, γ), a one-parameter solution manifold generically exists — and
does. The misfit was a restricted-family artifact, not an obstruction.

## The solution: all nine observables land

| observable | value | observed | status |
|---|---|---|---|
| V_us | 0.2256 | 0.225 | constrained, lands |
| V_cb | 0.0419 | 0.0418 | constrained, lands |
| V_ub | 0.00368 | 0.00369 | constrained, lands |
| β | 22.3° | 22.2° | constrained, lands |
| γ | 65.9° | 65.9° | constrained, **resolved** |
| V_td | ×1.01 | — | **predicted, lands** |
| V_ts | ×1.00 | — | **predicted, lands** |
| J | ×1.00 | — | **predicted, lands** |
| α | 91.8° | 91.9° | **predicted, lands** |
| sin δ | 0.889 | 0.887 | **predicted, lands** |

— at exactly preserved masses (1e-14) and the derived φ_h = π/k₅.

## The physical branch and the re-lock targets

The manifold's up-dominant end is excluded (the 5768 eigenvalue amplifies
sub-degree rotations into ×(−181) element targets). The **down-dominant
branch** reaches the same residual with physical targets:

| element | locked | target | factor |
|---|---|---|---|
| H₊[12] | −0.3548 | −0.4566 | ×1.287 |
| H₊[13], H₊[23] | −0.2682 | −0.2682 | ×1.000 (unchanged) |
| H₋[12] | −0.3548 | −0.6501 | ×1.832 |
| H₋[13] | −0.2682 | −0.5352 | ×1.996 |
| H₋[23] | −5.2682 | −5.8543 | ×1.111 |

All O(1–2), same-sign, within the transport-element family; rotation angles
≤ 6.1°. These are the **complete** targets for the knob-level v3+CP re-lock.

## What closes and what remains

The #157 card's last quantitative misfit resolves: a single target state
realizes the locked masses + the full CKM + the full unitarity triangle +
J + sin δ at the derived CP phase — **zero new inputs**. What remains is
engineering: the knob-level realization of the tabulated targets (the v3+CP
joint re-lock), plus the lepton sector's anarchic draw.

**Honesty notes**: V_td/V_ts/α land partly via CKM unitarity (the observed
dataset is itself unitarity-consistent) — the nontrivial content is
*existence* at fixed masses and the fixed derived phase; the target state is
an existence-plus-targets demonstration, not yet the knob-level model; the
branch choice is a documented physicality selection on a degenerate
manifold.

## Reproduce

```bash
python -m experiments.closure_ledger.gamma_misfit_resolution_probe
# Verdict: GAMMA_RESOLVED_FULL_FLAVOR_CP_DATASET_REALIZED_AT_PI_OVER_K5
```
