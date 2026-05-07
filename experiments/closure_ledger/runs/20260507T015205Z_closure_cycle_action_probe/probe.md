# Closure-cycle action quantum probe

**Run:** 2026-05-07T01:52:05+00:00

First concrete probe of the ℏ-origin research plan (`docs/hbar_origin_research_plan.md`). Tests prediction P1: the per-species Layer-1 closure-cycle phase is an integer multiple of 2π, and that integer N is the geometric content of the closure structure.

## Layer-1 ledger contributions per species

All values in units of π. `Φ_total / 2π` is the integer-quantum count if Layer 1 is exactly quantized.

| species | k | antipodal | Hopf | throat | uplift | Φ_total | Φ_total / 2π | integer? |
|---|---:|---:|---:|---:|---:|---:|---:|:---:|
| electron | 1 | 2.0π | 1.0π | 1.0π | 0.0π | 4.0π | 2.000000 | **✓** |
| muon | 3 | 6.0π | 1.0π | 1.0π | 0.0π | 8.0π | 4.000000 | **✓** |
| tau | 5 | 10.0π | 1.0π | 1.0π | 200.0π | 212.0π | 106.000000 | **✓** |

**Prediction P1:** PASS. Per-species integer quantum counts: `electron: N = 2`, `muon: N = 4`, `tau: N = 106`.

Each species' Layer-1 ledger sum is an integer multiple of `2π` (the action quantum candidate). The N counts grow as `{k+1, k+1, k+1+100·(k=5)} = {2, 4, 106}` — the first two are simply (k + 1) closure passes (electron and muon have no β·max(0, k-3)² contribution), and the τ row picks up the closure-quantum integer 100 from the 4β = 100·(2π) lock.

## R_MID = ℏ/(m_e c) identification (prediction P2)

Reduced Compton wavelength of the electron: λ_e_reduced = ℏ/(m_e c) = 3.861593e-11 cm.

Under P2, BAM's geometric `R_MID = 1` corresponds to λ_e_reduced. The closure cycle in physical units:

| species | N quanta | cycle (geom units) | cycle (cm) | cycle (Compton wavelengths) |
|---|---:|---:|---:|---:|
| electron | 2 | 12.5664 | 4.8526e-10 | 12.5664 |
| muon | 4 | 25.1327 | 9.7052e-10 | 25.1327 |
| tau | 106 | 666.0176 | 2.5719e-08 | 666.0176 |

By construction, the cycle length in Compton wavelengths equals N · 2π (since R_MID = λ_e_reduced is the conversion). This is **self-consistent** — but it does not yet PREDICT ℏ; it just identifies that under the canonical R_MID convention, the closure cycle traverses N · 2π reduced Compton wavelengths. Predicting ℏ requires R_MID to be geometrically determined (THESIS.md "Self-consistent throat radius" target), which is open.

## Layer-2 effects on integer quantization

Each row tests whether adding a candidate Layer-2 contribution preserves integer-quantum closure. (None of the prior probes' candidates did — but the table records by how much they break it.)

| Layer 2 candidate | Δφ per species (π) | new φ/2π per species | preserves integer? | new spread (rad) |
|---|---|---|:---:|---:|
| Layer 2 absent (ledger universal at 0 mod 2π) | [0.000, 0.000, 0.000] | [2.000, 4.000, 106.000] | ✓ | 0.0000 |
| C1 (eigenvector-weighted B1 modes) | [0.864, 0.788, 0.761] | [2.432, 4.394, 106.380] | — | 0.3257 |
| D1 (operator-valued V_j-V_i Hermitian matrix element) | [0.098, 1.914, 1.988] | [2.049, 4.957, 106.994] | — | 0.5767 |

Layer 1 alone preserves integer quantization (P1 PASS). C1 and D1 — the best Layer-2 scalar/operator candidates from the closure-ledger sweep — both BREAK integer quantization (they shift the per-species N by non-integer amounts). This is consistent with their FAIL on Layer-2 universality: a candidate that broke the per-species integer count but still made the residue universal would imply a NEW action quantum different from `2π`, which the lemma's structure does not allow.

## Verdict

**P1 PASS.** Each species' Layer-1 ledger sum is exactly `N · 2π` for integer N ∈ {2, 4, 106}. The geometric action quantum candidate is `2π` per closure unit; summing the four wired channels in the locked baseline lands every species on an integer multiple, with no fractional residue.

**P2 self-consistent, not yet predictive.** Identifying BAM's geometric R_MID = 1 with the reduced Compton wavelength of the electron (λ_e_reduced ≈ 3.86 × 10⁻¹¹ cm) makes the conversion `ℏ_SI = m_e R_MID c` a tautology. Predicting ℏ requires R_MID to be geometrically determined from a self-consistency condition (e.g. equilibrium throat radius for the locked mass spectrum), which is open.

**Layer 2 status.** No Layer-2 candidate from the closure-ledger sweep preserves integer quantization at the per-species level. This refines the earlier Layer-2 verdict: any Layer-2 closure form must contribute exactly an integer multiple of `2π` per species (NOT just be universal mod 2π) to preserve the action-quantum reading. The current candidates do neither.

## What's next

The ℏ-origin research plan identifies four sub-targets in decreasing tractability. P1 PASS narrows the next sub-probe scope:

- **(1) Layer 2 closure as integer-quantum contribution.** The next probe should test whether a worldline-integral Layer-2 form (rather than the per-mode WKB sum used by C1/D1) yields integer-multiple-of-2π per species. The natural candidate is the full `∮p·dq` over the closed orbit on the Tangherlini grid, with the closure cycle treated as a single integral rather than a sum of channels.
- **(2) Aharonov-Bohm form.** Compute the AB phase explicitly for a Hopf-fibre loop and verify it equals `2π` per spinor double-cover closure. This is the conceptually-clean geometric realization of the action quantum and should be verifiable in 50 lines of code given `geometrodynamics/hopf/connection.py`.
- **(3) Tangherlini eigenvalue ω_0 = m_e c² / ℏ?** Test whether `ℏ ω_0 = m_e c²` is consistent across species. If the lowest Tangherlini eigenvalue equals the electron rest energy in geometric units, that's an ℏ relationship.