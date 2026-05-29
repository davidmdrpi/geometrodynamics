# Neutrinoless double-beta (0νββ) effective-mass probe (PR #95)

The neutrino arc (#85–#94) fixed the *structure* of the BAM neutrino
sector. This probe turns it into a concrete, falsifiable prediction for
the one observable that directly tests the Majorana nature — the
effective Majorana mass

```
m_ββ = | Σ_i U_ei² m_i |
     = | c12²c13² m1 + s12²c13² m2 e^{iα21} + s13² m3 e^{iα31} |.
```

This combines everything the arc established:

  - **0νββ occurs at all** because the neutrino is **Majorana** (`c₁=0`,
    C-invariant, PR #86) — a Dirac neutrino would forbid it;
  - the **normal ordering** (PR #91: generations = cavity overtones,
    `m_ν ∝ m_D`) selects the NO band of `m_ββ`;
  - the **anarchic Majorana phases** (PR #94) populate the *whole* band,
    including the cancellation trough where `m_ββ → 0`;
  - the **light absolute scale** (PR #90: `m_ν ~ few meV`) places us in
    the deep-cancellation region.

## The prediction

Using the observed `Δm²` and mixing angles (BAM fixes the *ordering*,
*Majorana nature*, *phase distribution*, and *scale*, not `Δm²`
themselves), scanning the Majorana phases uniformly:

| m_lightest (meV) | m_ββ min (meV) | m_ββ max (meV) |
|---:|---:|---:|
| 0 | 1.45 | 3.68 |
| 2 | 0.15 | 5.11 |
| 3 | 0.02 | 5.88 |
| 5 | 0.02 | 7.49 |

| | m_ββ (meV) |
|---|---|
| **BAM (normal, light scale)** | ≲ 8 (cancellation to ~0) |
| inverted ordering (contrast) | 19–48 (floor ~19, no cancellation) |
| current bound (KamLAND-Zen) | 28–122 |
| next-gen reach (LEGEND-1000 / nEXO) | ~9–20 |

So BAM predicts `m_ββ ≲ 8 meV` — below the current bound (consistent with
the null result) and largely below even next-generation reach. **This is
a sharp falsifier:** a 0νββ discovery with `m_ββ ≳ 19 meV` (the inverted-
ordering floor) would imply inverted ordering or a quasi-degenerate
scale, contradicting the BAM normal-ordering + light-scale prediction.

## Tests

| # | test | finding |
|---|---|---|
| T1 | setup | `m_ββ = |Σ U_ei² m_i|`; needs Majorana + ordering + phases + scale |
| T2 | occurs | 0νββ occurs ⟸ neutrino Majorana ⟸ `c₁=0` (PR #86) |
| T3 | ordering | normal (PR #91) ⟹ NO band, below the IO floor ~19 meV |
| T4 | phases | anarchic (PR #94) ⟹ full band incl. cancellation → ~0 |
| T5 | scale | light (PR #90) ⟹ `m_ββ ≲ 8 meV` |
| T6 | experiment | below current (28–122 meV) & next-gen (~9–20); falsifiable |
| T7 | honest scope | qualitative firm; exact `m_ββ` a band |
| T8 | assessment | `ZERONUBB_OCCURS_NORMAL_ORDERING_M_BB_FEW_MEV` |

## Established and open

  - **Established (BAM-native):** 0νββ occurs (the neutrino is Majorana,
    PR #86); BAM selects the normal-ordering band (PR #91), below the IO
    floor; the anarchic Majorana phases (PR #94) populate the full band
    (cancellation to ~0 allowed); at the BAM light scale (PR #90)
    `m_ββ ≲ 8 meV` — below current bounds and a target/falsifier for
    next-gen experiments.

  - **Open:** a single `m_ββ` value — it is a band, because the lightest
    neutrino mass is unmeasured and the Majorana phases are anarchic
    (uniform). The exact spectrum (the PR #91 `χ_n`-corrected ratios, the
    absolute scale) and the specific phases are the residuals.

## Run

```
python -m experiments.closure_ledger.zeronubb_effective_mass_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_zeronubb_effective_mass_probe/`.
Expected verdict: `ZERONUBB_OCCURS_NORMAL_ORDERING_M_BB_FEW_MEV`, 8/8 PASS.
