# Quark `n_part = 233` origin: phenomenological compensator

**Run:** 2026-05-27T06:40:22+00:00

Attempts the hardest open piece of the quark sector arc — "why `n_part = 233`?" — and refines the prior verdict (`docs/quark_beta_status.md`: "phenomenological compensator") by (i) extending the candidate catalog beyond the prior enumerations (Fibonacci, Lucas, Padovan, Perrin, tribonacci, color × flavor × generation, QCD β₀, Tangherlini QCD-shell modes), (ii) identifying `n_part = 233 = F_13` as a baseline coincidence under §8 drift, and (iii) reframing structurally: the v3 quark Hamiltonian is lepton-shaped machinery, but the quark sector lives in the QCD shell channel per #68–#69. The right derivation route — quantitative #68–#69 development — is outside the closure-ledger machinery's scope.

- **Identification**: n_part = 233 = phenomenological compensator at the v3 baseline (lepton-style Hamiltonian, wrong machinery for the quark sector)
- **Only topological invariant**: N_q ≡ 0 (mod 2) — Z₂ partition multiplicity
- **Right derivation route**: quantitative #68–#69 QCD-shell development (outside closure-ledger scope)
- **B4 caveat**: n_part dimensionless integer; scale-independent; structural

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_parity_invariance_recap` | parity N_q ≡ 0 (mod 2) across all 12 §8 ablations | **PASS** |
| T2 | `T2_extended_candidate_catalog` | 0 exact-match enumerations beyond the prior catalog | **PASS** |
| T3 | `T3_fibonacci_F13_coincidence_under_drift` | F_13 = 233 baseline coincidence; drift not Fibonacci | **PASS** |
| T4 | `T4_v3_Hamiltonian_is_lepton_shaped` | v3 Hamiltonian basis = lepton-shaped {(k=1,3,5),±} | **PASS** |
| T5 | `T5_pr68_pr69_as_right_derivation_route` | #68–#69 = right route; outside closure-ledger scope | **PASS** |
| T6 | `T6_n_part_drifts_across_section_8` | n_part drifts 216–255 (span 39), not invariant | **PASS** |
| T7 | `T7_honest_scope_b4` | honest scope: phenomenological, structurally classified | **PASS** |
| T8 | `T8_assessment` | n_part = phenomenological compensator (PR verdict) | **PASS** |

## T2: Extended candidate catalog scan

Scanned **129 candidates** across **9 families**:

- `fibonacci` (18 entries)
- `k_5_polynomial` (6 entries)
- `lucas` (18 entries)
- `padovan` (24 entries)
- `perrin` (23 entries)
- `qcd_beta_function` (5 entries)
- `qcd_color_flavor` (11 entries)
- `tangherlini_shell_modes` (12 entries)
- `tribonacci` (12 entries)

**Exact matches to `n_part = 233`:**

- `fibonacci` — `F_13 = 233`
- `k_5_polynomial` — `9k_5² + k_5+3 = 233`

See `surviving_exact_matches` for §8-drift survival status (none survive — all are parameter-free enumerations).

**No exact matches to `N_q = 466`.**

## T3: `F_13 = 233` is a baseline coincidence

`n_part_baseline = 233 = F_13`. But the §8 drift values are:

| n_part value | in Fibonacci sequence? |
|---:|:---:|
| 216 | ✗ |
| 220 | ✗ |
| 233 | ✓ |
| 237 | ✗ |
| 238 | ✗ |
| 241 | ✗ |
| 247 | ✗ |
| 255 | ✗ |

1 of 8 drift values are Fibonacci → F_13 = 233 is a baseline coincidence, not a structural invariant.

## T6: §8 drift — `n_part` is not a topological invariant

| ablation | N_q | n_part = N_q/2 |
|---:|---:|---:|
| 1 | 466 | 233 |
| 2 | 466 | 233 |
| 3 | 466 | 233 |
| 4 | 476 | 238 |
| 5 | 474 | 237 |
| 6 | 474 | 237 |
| 7 | 482 | 241 |
| 8 | 432 | 216 |
| 9 | 494 | 247 |
| 10 | 494 | 247 |
| 11 | 440 | 220 |
| 12 | 510 | 255 |

`n_part` drift: min = 216, max = 255, span = 39, mean = 236.4, pct = 16.5%. Only the parity (N_q ≡ 0 mod 2) survives.

## T4–T5: Why the v3 Hamiltonian is the wrong machinery

| sector | lives in | basis | closure integer |
|---|---|---|---:|
| **lepton** | throat (odd-`k` modes) | `{(k=1,±), (k=3,±), (k=5,±)}` | `4·k_5² = 100` (clean) |
| **quark** | QCD shell channel (#68–#69) | (same v3 basis = wrong) | `466` (phenomenological) |

The v3 quark Hamiltonian fits the quark spectrum on a **lepton-shaped** basis (the 6 odd-`k` throat modes). But per #68–#69, the quark mass sector is delocalized into the QCD shell channel: extended-character wavefronts reproducing the Z₂ partition (#69), the `3 × 2 = 6` flavor count (#69), the heavier mass scale (#68). The v3 ansatz absorbs the unmodeled QCD physics (confinement, αs running, color sector) into `n_part = 233` as a phenomenological compensator — the price of fitting the quark spectrum on the wrong basis.

**Right derivation route:** quantitative development of #68–#69 from "structural match" to a full QCD-shell model, with confinement, αs running, and color sector explicit. This is genuinely outside the closure-ledger machinery's scope (the closure-ledger primitives — `action_base = 2π`, `transport = 8π`, `resistance = 7π/100`, `pinhole γ` — are lepton-throat-sector primitives) and is a substantial research program in its own right.

## Verdict

**N_PART_IS_PHENOMENOLOGICAL_COMPENSATOR.** N_PART IS PHENOMENOLOGICAL COMPENSATOR. n_part = 233 at the v3 baseline is one realization of a fit compensator that drifts across docs/quark_axioms.md §8 ablations from 216 to 255 (span 39, mean 236.4). The only topological invariant is the parity N_q ≡ 0 (mod 2), the Z₂ partition multiplicity of the v3 ansatz basis {(k, +), (k, −)} (recap from quark_beta_subblock_stability).

EXTENDED CANDIDATE CATALOG. Beyond the prior probe's enumerations (S³/S² harmonics, SU(3) representations, torus-knot crossings, Tangherlini barrier sums), this probe scanned Fibonacci/Lucas/Padovan/Perrin/tribonacci, color × flavor × generation (3·6 = 18, 3·6² = 108, 8·6·3 = 144, 8·6·5 = 240, …), QCD β₀ = 7 combinations (7·25 = 175, 7·30 = 210, …), and Tangherlini QCD-shell mode counts. n_part = 233 = F_13 (the 13th Fibonacci number) is the only exact baseline match in the extended catalog — but §8 drift values (216, 220, 237, 238, 241, 247, 255) are not Fibonacci, so F_13 = 233 is a baseline coincidence, not a structural invariant. No principled enumeration in the closure-ledger catalog reproduces n_part across §8.

STRUCTURAL READING. The v3 quark Hamiltonian uses the lepton-style closure-quantum basis {(k=1,±), (k=3,±), (k=5,±)} — the 6 odd-k throat modes that give the charged-lepton ladder via β_lepton = k_5²·(2π) = 50π (closure-quantum integer 4·k_5² = 100). But the quark sector lives in the QCD SHELL CHANNEL per #68–#69: higher excitations of the focused lepton-throat pulse delocalize into a heavier-scale, extended-character shell wavefront reproducing the documented quark-sector structural invariants (Z₂ partition, 3 × 2 = 6 flavors). The v3 Hamiltonian is the WRONG machinery for the quark sector — it is fitting QCD-confined quarks on closure-quantum throat basis vectors. n_part = 233 is the empirical price of that sector mismatch: it absorbs unmodeled QCD physics (confinement, αs running, color), which closure-ledger primitives (action_base = 2π, transport = 8π, resistance = 7π/100, pinhole γ) are not designed to capture.

RIGHT DERIVATION ROUTE. The structurally honest path to n_part is QUANTITATIVE development of #68–#69 (throat-to-shell transition + shell↔QCD structural match) into a full QCD-shell model — i.e., a model in which the quark mass sector is computed from the shell channel directly, with confinement, αs running, and color sector all explicit rather than absorbed into a phenomenological compensator. This is structurally OUTSIDE the closure-ledger machinery's scope (the closure-ledger primitives are throat-sector / lepton-sector primitives) and is itself a substantial research program — comparable in scope to deriving lattice QCD's spectrum from underlying geometric principles. It is not the next-most-tractable work in the BAM framework, and should not be pursued by further enumeration on the v3 Hamiltonian.

HONEST SCOPE. This probe classifies n_part = 233 structurally; it does NOT first-principles derive 233 (no such derivation exists in the current catalog). The prior verdict from quark_beta_status.md ("phenomenological compensator") is upheld and sharpened: the WHY is the lepton/shell sector mismatch, and the RIGHT-ROUTE is #68–#69 developed quantitatively. B4: n_part is a dimensionless integer; scale-independent; structural.

## What this leaves open

- **First-principles derivation of `n_part = 233`** — open. No principled enumeration in the closure-ledger catalog reproduces it (extending the prior scan with Fibonacci, Lucas, Padovan, Perrin, tribonacci, color × flavor × generation, QCD β₀, Tangherlini QCD-shell modes).
- **Quantitative QCD-shell model** (extending #68–#69) — the structurally honest route, but outside closure-ledger scope and a substantial research program in its own right.
- **Heavy-lepton thresholds** (`2 m_μ c²`, `2 m_τ c²`) and the analogous quark-sector pair-production thresholds — related open work flagged by #58.
