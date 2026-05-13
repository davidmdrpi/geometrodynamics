# Resistance-reading disambiguation probe

**Run:** 2026-05-13T00:44:03+00:00

Disambiguates the two within-1 % candidates for `resistance_scale` identified by the opening transport/resistance origin probe:

- **(A) closure_quantum**: `resistance = 7π/100 ≈ 0.2199` (+0.94 %), R-independent.
- **(B) eigenfrequency**: `resistance = 4·(ω(1,0;R) − 1) ≈ 0.2189` (+0.48 %), R-dependent (varies with the bisection).

Each reading is paired with `transport = 8π = 4·(2π)` (the transport closure-quantum reading) and re-runs the R_OUTER self-consistency loop from probe 8 — bisecting R*_μ and R*_τ independently for each reading.

## Results

| reading | description | R*_μ | R*_τ | |ΔR*|/R*_μ | γ at R* | transport | resistance | ω(1,0) | err μ | err τ |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `locked_baseline` | transport = 25.1, resistance = 0.217869 (locked values) | 1.262239 | 1.262338 | 0.0078% | 22.4994 | 25.1000 | 0.2179 | 1.0537 | -0.0000% | +0.1614% |
| `closure_quantum` | transport = 8π, resistance = 7π/100 (both R-independent) | 1.262636 | 1.262662 | 0.0021% | 22.5076 | 25.1327 | 0.2199 | 1.0535 | +0.0000% | +0.0436% |
| `eigenfrequency` | transport = 8π, resistance = 4·(ω(1,0;R) − 1) (R-dependent) | 1.258316 | 1.258284 | 0.0025% | 22.4167 | 25.1327 | 0.2220 | 1.0555 | -0.0000% | +0.1453% |
| `closure_quantum_locked_transport` | transport = 25.1 (locked), resistance = 7π/100 | 1.258067 | 1.258217 | 0.0120% | 22.4113 | 25.1000 | 0.2199 | 1.0556 | +0.0000% | +0.2602% |
| `eigenfrequency_locked_transport` | transport = 25.1 (locked), resistance = 4·(ω(1,0;R) − 1) | 1.259967 | 1.259920 | 0.0038% | 22.4520 | 25.1000 | 0.2190 | 1.0547 | +0.0000% | +0.2146% |

## Verdict

Cross-species agreement — locked: 0.0078%, closure: 0.0021%, eig: 0.0025%.
R*-match to locked baseline — closure: 0.0314%, eig: 0.3109%.
γ at R* vs canonical 22.5 — closure: 0.0337%, eig: 0.3704%.

**closure-quantum reading PREFERRED.** Both readings close the loop at high cross-species precision, but closure-quantum's R* (1.262636) matches the locked baseline (1.262239) to 0.0314%, while eigenfrequency drifts by 0.3109%. The closure-quantum reading also lands γ at R* on 22.5 to 0.0337%, preserving the pinhole-origin anchor; eigenfrequency drifts γ by 0.3704%.

### Reading-by-reading interpretation

**Locked baseline** (control). Reproduces probe 8: R* ≈ 1.2622, cross-species agreement 0.0078 %, γ = 22.499 ≈ 22.5. ✓

**Closure-quantum reading**. `transport = 8π`, `resistance = 7π/100`, both R-independent. The bisection selects R*_μ = 1.262636, R*_τ = 1.262662, agreement 0.0021 %. γ at R* = 22.508.

**Eigenfrequency reading**. `transport = 8π`, `resistance = 4·(ω(1,0;R) − 1)`. resistance varies with R during the bisection — this is the FULLY GEOMETRIC loop (no constants other than m_e and the closure-quantum integers). R*_μ = 1.258316, R*_τ = 1.258284, agreement 0.0025 %. ω(1,0) at R* = 1.0555 (resistance at R* = 0.2220).

## Implication

**Cross-species agreement** — Locked: 0.0078% · Closure-quantum: 0.0021% · Eigenfrequency: 0.0025%.
**R*-match to locked** — Closure-quantum: 0.0314% · Eigenfrequency: 0.3109%.

**Closure-quantum reading wins.** Both readings produce tighter cross-species agreement than the locked baseline (0.0078 %), but the closure-quantum reading's R* lands on the locked-baseline R* to 0.0314 %, while the eigenfrequency reading drifts by 0.3109 %. γ at R* under closure-quantum matches the canonical 22.5 to 0.0337 %, preserving the pinhole-origin anchor identified in `pinhole_origin_probe`. The structural reading is:

  `transport = 8π = 4·(2π)` — 4th closure quantum.
  `resistance = 7π / 100` — closure-quantum fraction.

Both are pure closure-quantum invariants of the BAM framework. With this reading the R_OUTER self-consistency loop closes on principled inputs alone — no transport or resistance constants are required as external inputs.
