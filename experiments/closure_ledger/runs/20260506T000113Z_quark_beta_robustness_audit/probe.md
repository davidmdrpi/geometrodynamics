# Quark β: §8 robustness audit sub-probe

**Run:** 2026-05-06T00:01:13+00:00
**k_5 = 5**, **baseline N = 466** (m = 93, δ_signed = +1).

Audits the structural reading `N_q = m·k_5 + δ` against the EXACT §8 N-ablation values logged in `docs/quark_axioms.md`. Earlier probes used extrapolated N values for some perturbation classes; this audit drops them and reports only what the §8 log actually documents.

## Per-ablation (m, δ) decompositions

| ablation | class | N | m | δ_signed | drift | fit residual |
|---|---|---:|---:|---:|---:|---:|
| baseline (anchor=d, PDG, min_eig) | baseline | 466 | 93 | +1 | +0 | 1.6% |
| PDG × 1.10 (uniform scale) | uniform_scale | 466 | 93 | +1 | +0 | 1.6% |
| PDG × 0.90 (uniform scale) | uniform_scale | 466 | 93 | +1 | +0 | 1.6% |
| anchor = s | anchor | 476 | 95 | +1 | +0 | 1.1% |
| anchor = c | anchor | 474 | 94 | -1 | -2 | 1.0% |
| anchor = b | anchor | 474 | 94 | -1 | -2 | 1.9% |
| anchor = t | anchor | 482 | 96 | +2 | +1 | 1.6% |
| c × 1.10 | single_species | 432 | 86 | +2 | +1 | 1.8% |
| b × 1.10 | single_species | 494 | 98 | -1 | -2 | 4.2% |
| t × 1.10 | single_species | 494 | 98 | -1 | -2 | 5.5% |
| t × 0.90 | single_species | 440 | 88 | +0 | -1 | 3.7% |
| all ±5% (deterministic) | all_species | 510 | 102 | +0 | -1 | 4.6% |

## δ distribution across all ablations

| signed δ | count | fraction |
|---:|---:|---:|
| -1 | 4 | 33% |
| +0 | 2 | 17% |
| +1 | 4 | 33% |
| +2 | 2 | 17% |

## m envelope

baseline m = 93, m range = [86, 102], drift envelope ±9 units.

## Claim verdicts

### Claim 1: δ = +1 in baseline / uniform-scale runs
- Under perturbations that do NOT change the physics (just the MeV anchor or a uniform mass rescaling), δ should equal the baseline +1.
- **Expected:** all 3 runs at δ = +1
- **Observed:** 3/3 at δ = +1 (rate 100%)
- **Verdict:** ✓ **TRUE**

### Claim 2: δ = +1 across single-species perturbations
- Under single-species ±10% perturbations (c, b, t each ×1.10 or ×0.90), δ should remain at +1 if the structural reading is robust.
- **Expected:** all 4 runs at δ = +1
- **Observed:** 0/4 at δ = +1; observed δ values: [-1, -1, 0, 2]
- **Verdict:** ✗ **FALSE**

### Claim 3: δ = +1 across anchor-species changes
- Anchor-species choice (s, c, b, t instead of d) is a convention; structural reading expects δ = +1 throughout.
- **Expected:** all 4 runs at δ = +1
- **Observed:** 1/4 at δ = +1; observed δ values: [-1, -1, 1, 2]
- **Verdict:** ✗ **FALSE**

### Claim 4: overall δ-invariance rate ≥ 50%
- Aggregated over every logged ablation, at least half should sit at δ = +1 if the boundary correction is a structural invariant.
- **Expected:** ≥ 50% of runs at δ = +1
- **Observed:** 33% of runs at δ = +1
- **Verdict:** ✗ **FALSE**

## Verdict summary

- **TRUE:** 1/4
- **PARTIAL:** 0/4
- **FALSE:** 3/4


## Headline correction to prior probes

**The earlier boundary-correction probe overstated the +1 invariance.** Using the actual §8 ablation N values, only **33%** of runs sit at δ = +1; the remainder are distributed across δ ∈ {-1, 0, 1, 2}. The δ-distribution is approximately uniform on the integers in (−k_5/2, k_5/2] = (−2, 2], consistent with the §8 conclusion that 'N is a compensator, not a topological invariant'.

**Implications for the structural decomposition:**

- The structural reading `N_q = ((k_5−1)·k_5 + 2·k_5(k_5+2) + N_c)·k_5 + (N_c − 2)` correctly fits the *baseline* value 466 and the *uniform-scale* invariance, but does NOT predict the perturbation behaviour that §8 documents. Under per-species or anchor changes, both m and δ drift.
- The +1 = N_c − 2 boundary-correction interpretation is consistent at the baseline but is NOT a robust topological invariant. It survives uniform-mass-scale perturbations (those don't change the underlying physics, just the MeV anchor) but breaks under anchor-species and per-species perturbations.
- This is consistent with the §8 conclusion: **N (and therefore m and δ both) is a fit compensator. The structural decomposition is a useful descriptor of the BASELINE locked value, not a derivation of N.**

## Conclusion

The structural decomposition `m_q = (k_5−1)·k_5 + 2·k_5(k_5+2) + N_c` and boundary correction `δ = N_c − 2` remain **descriptively useful for the baseline** and the uniform-scale runs — these read off the locked N = 466 in a geometrically meaningful way. They are NOT a derivation of N in the strong sense: under per-species and anchor perturbations, both m and δ wander, exactly as §8 reports. The repo's own conclusion stands: **β_quark is phenomenological**, and the cleanest structural reading is a post-hoc identification of the locked value, not a first-principles prediction.

The next concrete sub-target is to ask whether **a different perturbation-robust quantity** (e.g. the NN_l + N_c piece alone, or the lepton-block sub-piece) is what's actually topologically locked, with the rest being the compensator.