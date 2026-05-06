# Quark β: focused boundary-correction probe

**Run:** 2026-05-06T00:03:13+00:00
**Targets:** N_l = 100, N_q = 466, ΔN = 366; k_5 = 5.
**Natural δ set:** [-10, -6, -5, -4, -2, -1, 0, 1, 2, 4, 5, 6, 10].

## All structural-unit fits, ranked by joint cleanness

| unit | u | N_l = m·u+δ | N_q = m·u+δ | ΔN = m·u+δ | Σ|δ| | δ ∈ natural? | δ_q == δ_ΔN? | δ_l = 0? |
|---|---:|---|---|---|---:|---|---|---|
| `k_5` | 5 | 20·5 +0 | 93·5 +1 | 73·5 +1 | 2 | ✓ | ✓ | ✓ |
| `k_5 + 1` | 6 | 17·6 -2 | 78·6 -2 | 61·6 +0 | 4 | ✓ | — | — |
| `2 k_5` | 10 | 10·10 +0 | 47·10 -4 | 37·10 -4 | 8 | ✓ | ✓ | ✓ |
| `k_5² + 1` | 26 | 4·26 -4 | 18·26 -2 | 14·26 +2 | 8 | ✓ | — | — |
| `(k_5 − 1)²` | 16 | 6·16 +4 | 29·16 +2 | 23·16 -2 | 8 | ✓ | — | — |
| `k_5(k_5+1)/2` | 15 | 7·15 -5 | 31·15 +1 | 24·15 +6 | 12 | ✓ | — | — |
| `(k_5 + 1)²` | 36 | 3·36 -8 | 13·36 -2 | 10·36 +6 | 16 | — | — | — |
| `k_5²` | 25 | 4·25 +0 | 19·25 -9 | 15·25 -9 | 18 | — | ✓ | ✓ |
| `k_5² − 1` | 24 | 4·24 +4 | 19·24 +10 | 15·24 +6 | 20 | ✓ | — | — |
| `k_5(k_5 + 2)` | 35 | 3·35 -5 | 13·35 +11 | 10·35 +16 | 32 | — | — | — |
| `(2 k_5)²` | 100 | 1·100 +0 | 5·100 -34 | 4·100 -34 | 68 | — | ✓ | ✓ |
| `4 (k_5² + 1)` | 104 | 1·104 -4 | 4·104 +50 | 4·104 -50 | 104 | — | — | — |

## Best fit

**Unit:** `k_5` = 5.

**Decomposition:**

- `N_lepton  = 20 · 5 +0` ⇒ residue δ = **0**
- `N_quark   = 93 · 5 +1` ⇒ residue δ = **1**
- `ΔN_q-l   = 73 · 5 +1` ⇒ residue δ = **1**

**Structural cleanness criteria met:** δ_lepton = 0 (lepton sector closes exactly on a multiple of k_5) AND δ_quark = δ_gap = 1 (the quark sector carries a single fixed boundary correction +1 across both targets).

This is consistent with reading the quark closure-quantum count as `m · k_5 + 1`, where m is the integer winding number that absorbs all per-species mass-perturbation drift, and the +1 is a single sector-level boundary correction tied to the orientation-reversing throat closure (parallel to the lepton-sector exact reading `N_l = 20 · k_5`).

## Robustness re-read under documented N-drifts

From `docs/quark_axioms.md` §8 N-stability ablation, restated under the (m, δ) decomposition with u = `k_5`:

| perturbation | N_q | m | δ | δ-drift |
|---|---:|---:|---:|---:|
| baseline (anchor=d, PDG, min_eig) | 466 | 93 | +1 | +0 |
| PDG × 1.10 (uniform scale) | 466 | 93 | +1 | +0 |
| PDG × 0.90 (uniform scale) | 466 | 93 | +1 | +0 |
| anchor = s | 476 | 95 | +1 | +0 |
| anchor = c | 474 | 95 | -1 | -2 |
| anchor = b | 474 | 95 | -1 | -2 |
| anchor = t | 482 | 96 | +2 | +1 |
| c × 1.10 | 432 | 86 | +2 | +1 |
| b × 1.10 | 494 | 99 | -1 | -2 |
| t × 1.10 | 494 | 99 | -1 | -2 |
| t × 0.90 | 440 | 88 | +0 | -1 |
| all ±5% (deterministic) | 510 | 102 | +0 | -1 |

**δ-invariance rate under documented drifts:** 33% of perturbations leave δ at its baseline value of +1.

**Only 33% of documented drifts leave δ at baseline**; under per-species and anchor perturbations, both m and δ wander. The structural reading is descriptively useful for the BASELINE locked value but is NOT a robust topological invariant — consistent with §8's headline conclusion that N (and therefore both m and δ) is a fit compensator. See `quark_beta_robustness_audit.py` for the full audit.

## Verdict

**The cleanest structural reading is `m · k_5 + 1`** with the lepton sector exact (δ = 0) and the quark sector carrying a single fixed boundary correction δ = +1. The pattern is suggestive but not yet derived: `+1` could be a Z₂ partition residue (orientation-flip closure), an l = 0 s-wave closure quantum, or a single-mouth boundary contribution. Identifying which is the next probe target.