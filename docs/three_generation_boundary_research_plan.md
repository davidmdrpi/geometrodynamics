# Three-generation boundary: β-uplift cutoff + throat-shell availability

Combines the two mechanisms identified in PRs #67–#69 to pin the sharp
three-generation boundary of the charged-lepton sector at `k ≤ 5`:

  - **β-uplift (closure-quantum growth)** — within the locked lepton
    surrogate, the closure-quantum sum is
    `Φ_avail(k) = 2π(k+1) + 50π·max(0, k−3)²` (the odd-k lemma). The
    second term grows **quadratically**: `25·(k−3)²` closure quanta added
    starting at `k = 5` (tau). It scales the muon/tau masses correctly
    within the three generations and quantifies the steep closure-quantum
    cost of any would-be higher generation.
  - **Throat-shell mode availability (#68)** — the throat-localized
    fermionic ladder has exactly **three** modes (`n = 0, 1, 2`); from
    `n = 3` the modes are shell-saturated (participation `→ 2/3`, the
    QCD/quark sector, #69). The hard cutoff in mode availability.

The combination pins `k ≤ 5`: closure arithmetic alone closes for *any*
integer `k` (#67), so the cutoff is **not arithmetic** — it is the
mode-availability boundary, with the β-uplift's quadratic growth
quantifying the within-lepton scaling.

## The k ↔ n mapping (B2 reading)

The closure-ledger B2 candidate (`sk_bridge.py`) maps each odd depth `k`
to the radial overtone `n = (k − 1) / 2` at `l = 1`:

| k | n | species (#68 + #69) | β-uplift `50π·(k−3)²` | N_total = (k+1) + 25(k−3)² | throat-localized? |
|---:|---:|---|---:|---:|---|
| 1 | 0 | electron | 0 | 2 | **yes** (lepton) |
| 3 | 1 | muon | 0 | 4 | **yes** (lepton) |
| 5 | 2 | tau | 100·(2π) | 106 | **yes** (lepton, last) |
| 7 | 3 | (would-be 4th gen) | 400·(2π) | 408 | **NO** — shell-saturated (#68); QCD sector (#69) |
| 9 | 4 | (would-be 5th gen) | 900·(2π) | 910 | NO — shell (QCD) |

The would-be 4th charged lepton at `k = 7` maps to `n = 3` — the
shell-saturated mode (#68), which carries the quark/QCD structural
invariants (#69), not the charged-lepton throat-localized character. So
there is no `k = 7` charged lepton: **the lepton sector ends with the
third throat-localized mode (τ at `k = 5, n = 2`)**.

## The two mechanisms together

  1. **β-uplift quadratic growth** (50π·(k−3)²) — the closure-quantum
     cost: 0 (e, μ), 100 (τ), then jumps to 400, 900, … beyond. The
     within-lepton mass-scaling structure of the locked surrogate
     (`compute_knotted_lepton_spectrum`, `TAU_BETA_50PI`), giving the
     correct μ, τ ratios at `k = 3, 5`.
  2. **Throat-shell availability** (#68) — only `n = 0, 1, 2` are
     throat-localized; `n ≥ 3` is shell-saturated (QCD, #69). The
     hard boundary.

The closure arithmetic by itself does NOT cut (#67: `Φ_avail(k) ≡ 0 mod
2π` for every integer `k`). The lepton/quark boundary is the
mode-availability one; the β-uplift quantifies the within-lepton scaling
that the surrogate uses to fit the muon and tau masses correctly given
the three generations.

## Honest scope

  - **Is:** an explanation of the sharp `k ≤ 5` boundary as the join of
    the β-uplift quadratic growth (within-lepton scaling) and the
    throat-shell availability (#68's mode-availability cutoff at `n = 3`),
    with the `k ↔ n` map (B2) showing the would-be 4th generation lives
    in the shell/QCD sector (#69).
  - **Is not:** a first-principles derivation of `β_lepton = 50π` itself,
    nor of `n_part = 233` for the quark sector. The closure-quantum
    `β_lepton` is the locked surrogate's structure (well-documented in
    `docs/hbar_origin_status.md`); the structural origins (`8π = 4·(2π)`
    transport, `7π/100` resistance) are established in the ℏ-origin
    thread.

## B4 accounting

`k` is a **dimensionless integer**; the closure-quantum integers
`N_layer1(k), N_uplift(k), N_total(k)` are integer counts; the
mode-availability cutoff at `n = 3` is geometric. The boundary is
topological/structural — independent of the single anchor `m_e`. The
mass *values* at `(e, μ, τ)` carry the scale; the count of generations
does not.

## Tests

  T1. **β-uplift quadratic growth.** `Φ_uplift(k) = 50π·(k−3)²` activates
      at `k = 5`; the closure-quantum integer is `25·(k−3)²` → 0, 0,
      100, 400, 900, … (quadratic explosion).
  T2. **Closure arithmetic closes for any k (recap #67).** Cutoff is not
      arithmetic.
  T3. **k ↔ n mapping (B2).** `n = (k − 1) / 2`; the three leptons
      `(k=1,3,5) ↔ (n=0,1,2)`; would-be 4th gen `k=7 ↔ n=3`.
  T4. **Throat-shell availability cutoff (recap #68 + #69).** `n = 0,1,2`
      throat-localized (leptons); `n ≥ 3` shell-saturated (QCD).
  T5. **The would-be 4th generation is in the shell sector.** `k=7 ↔
      n=3` shell-saturated → no charged lepton at k=7; the
      three-generation cutoff IS the lepton/quark sector boundary.
  T6. **Combined mechanism.** β-uplift quadratic cost + throat-shell
      availability pin `k ≤ 5`.
  T7. **Falsification / B4.** a 4th charged lepton would falsify; BAM
      passes (none observed; `n = 3` is shell). `k` dimensionless;
      structural.
  T8. **Assessment.**

## Verdict structure

  - **THREE_GENERATIONS_PINNED** (expected): the sharp `k ≤ 5`
    charged-lepton boundary is pinned by the join of (i) the β-uplift's
    quadratic closure-quantum growth `50π·(k−3)²` (within-lepton scaling)
    and (ii) the throat-shell mode-availability cutoff (#68) at `n = 3`,
    where the B2 mapping `n = (k−1)/2` sends the would-be 4th generation
    `k = 7` to the shell-saturated sector (the QCD/quark structural
    sector, #69), so no charged lepton exists there. Closure arithmetic
    by itself does not cut (#67). The three observed generations are
    exactly the three throat-localized odd-k fermionic modes.

  - **NO_SHARP_BOUNDARY**: the mechanisms do not combine to pin `k ≤ 5`.

## What this leaves open

  - **First-principles `β_lepton = 50π`.** The structural origins are in
    the ℏ-origin thread (`8π = 4·(2π)`, `7π/100`); a clean derivation of
    the `50π` from these primitives sharpens further.
  - **Why exactly the THREE modes** (the n=3 cutoff in the cavity) — the
    cavity geometry sets it (`R_OUTER − R_MID` and the centrifugal
    barrier); a deeper structural argument tied to the bulk geometry is
    future work.

## Cross-references

  - `docs/odd_k_closure_lemma.md` — the closure arithmetic.
  - `docs/even_k_absence_research_plan.md` — odd-k = fermionic (#67).
  - `docs/throat_to_shell_transition_research_plan.md` — the throat-shell
    crossover (#68).
  - `docs/shell_to_qcd_match_research_plan.md` — the shell ↔ QCD match
    (#69).
  - `docs/hbar_origin_status.md` — the closure-quantum origins
    (`8π = 4·(2π)`, `7π/100`).
  - `geometrodynamics/tangherlini/lepton_spectrum.py` — `TAU_BETA_50PI`,
    `LEPTON_BASELINE_DEPTHS = (1, 3, 5)`.
  - `experiments/closure_ledger/three_generation_boundary_probe.py` —
    this probe.
