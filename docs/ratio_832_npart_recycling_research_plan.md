# Is 832 = N_q + ΔN an independent scale-ratio selection, or recycled n_part? (PR #107)

PR #106 left the lepton/QCD scale ratio `√σ/m_e ≈ 830` as the one
UNDERIVED open dimensionless residual whose derivation would reduce BAM to
a single anchor (the bulk-gravity scale `G`). A tempting candidate
"derivation":

```
N_lepton = 100,   N_q = 466,   ΔN = N_q − N_lepton = 366,
N_q + ΔN = 832  ≈  √σ/m_e ≈ 830   (0.2% match).
```

This probe tests, skeptically, whether the bulk shell-stress integral
INDEPENDENTLY selects 832 — or whether 832 just recycles the
phenomenological compensator `n_part`. **Answer: it recycles `n_part`. The
match is a baseline coincidence; the channel-normalisation derivation via
this route fails.**

## 832 is built from n_part

`N_q + ΔN = 2·N_q − N_lepton = 2·(2·n_part) − 4·k_5² = 4·n_part − 4·k_5²
= 4·233 − 100 = 832`. A linear function of `n_part` — the quark closure
integer PR #76/#97 established is a PHENOMENOLOGICAL COMPENSATOR (it
absorbs the quark flavor puzzle; drifts 216–255 across the `quark_axioms`
§8 ablations; only its parity is invariant). 832 inherits that status.

## The decisive test: §8 drift

If 832 independently selected the FIXED observed ratio 830, it would be
§8-stable. It is not. Propagating `n_part ∈ {216, …, 255}` through
`4·n_part − 100`:

```
832-analogue ∈ [764, 920]   (span 156, ≈ ±9%),
```

so the quantity drifts ~±10% while `√σ/m_e = 830.3` is fixed. The baseline
`n_part = 233` merely lands it at 832 (0.2% from 830) — a **baseline
coincidence**, the same kind as `50π·k_5 = 785` (PR #106) and
`F_13 = 233` (PR #76), not a stable selection.

## No independent bulk shell-stress integral selects 832

The genuine bulk shell-stress integrals over the Tangherlini geometry are
`O(10–70)`, nowhere near 466/832:

| bulk shell-stress integral | value |
|---|---:|
| `Σ ω²(l=1, n=3..5)` | 69.8 |
| `Σ (n+1)π` (Bohr–Sommerfeld closure) | 47.1 |

The number 466 enters ONLY through the v3 Hamiltonian closure count
`4β_quark/(2π) = 2·n_part` — i.e. through the fit. No independent integral
yields ~466 or ~832.

## The circularity

`n_part` was FIT to reproduce the quark spectrum (which already encodes
the physical scales). Recovering a scale ratio from `n_part` is therefore
circular — you get back (an unstable version of) what was put in. So
`832 ≈ 830` is the compensator echoing the spectrum it was fit to, not a
derivation of the lepton/QCD hierarchy.

## What this means for the ledger

The channel-normalisation derivation via `N_q + ΔN` FAILS. The PR #106
status stands unchanged: `√σ/m_e ≈ 830` remains an UNDERIVED open
dimensionless residual, and BAM remains at one foundational scale (`G`) +
one open ratio. A genuine derivation would need an INDEPENDENT bulk
shell-stress integral (not the v3-fit closure count) that selects ~830
AND is §8-stable — which is not available.

## Tests

| # | test | finding |
|---|---|---|
| T1 | observation | `N_q+ΔN = 2N_q−N_lepton = 832 ≈ 830` (0.2%) — tempting |
| T2 | built from n_part | `832 = 4·n_part − 4·k_5²` |
| T3 | n_part = compensator | §8-drifts 216–255, parity-only (PR #76/#97) |
| T4 | §8-drift (decisive) | `4·n_part−100` drifts 764–920 (±9%) vs fixed 830 |
| T5 | no independent integral | shell integrals `O(10–70)`, never ~466/832 |
| T6 | circularity | `n_part` fit to the spectrum it would "predict" |
| T7 | ledger unchanged | `√σ/m_e` underived; one scale `G` + one ratio |
| T8 | assessment | `RATIO_832_RECYCLES_N_PART_NOT_INDEPENDENT` |

## Established and open

  - **Established (BAM-native):** `832 = 2N_q − N_lepton = 4·n_part −
    4·k_5²` is built from the compensator `n_part`; it §8-drifts 764–920
    (so the `832 ≈ 830` match is a baseline coincidence); and no
    independent bulk shell-stress integral selects ~466/832. We are
    recycling `n_part` — the channel-normalisation derivation fails.

  - **Open (unchanged from PR #106):** a genuine derivation of `√σ/m_e ≈
    830` from an INDEPENDENT, §8-stable bulk integral. `√σ/m_e` stays an
    underived open residual; BAM stays at one scale `G` + one open ratio.

## Cross-references

  - `docs/scale_count_anchors_research_plan.md` — PR #106, the underived
    `√σ/m_e ≈ 830` ratio this candidate tried (and failed) to derive.
  - `docs/quark_npart_origin_research_plan.md` / `docs/quark_beta_status.md`
    — PR #76/#97, the `n_part` compensator and its §8 drift.

## Run

```
python -m experiments.closure_ledger.ratio_832_npart_recycling_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_ratio_832_npart_recycling_probe/`.
Expected verdict: `RATIO_832_RECYCLES_N_PART_NOT_INDEPENDENT`, 8/8 PASS.
