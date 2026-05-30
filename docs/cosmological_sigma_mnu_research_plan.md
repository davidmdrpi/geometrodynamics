# Cosmological Σm_ν prediction probe (PR #96)

The same BAM neutrino spectrum behind the 0νββ prediction (PR #95) also
fixes the **sum of neutrino masses** `Σm_ν = m1 + m2 + m3` — the quantity
probed by cosmology (CMB lensing + large-scale structure / BAO). Where
0νββ tests the Majorana nature, Σm_ν tests the ordering and the absolute
scale, and the two come from one spectrum.

## The prediction

BAM fixes the **normal ordering** (PR #91: generations = cavity
overtones, `m_ν ∝ m_D`) and a **light absolute scale** (PR #90: lightest
~ few meV). With the observed `Δm²`:

```
Σm_ν(m_lightest) = m1 + √(m1² + Δm²_21) + √(m1² + Δm²_31)
```

| m_lightest (meV) | Σm_ν (meV) |
|---:|---:|
| 0 | 58.7 (the NO floor) |
| 2 | 60.9 |
| 5 | 65.2 |
| 10 | 74.2 |

So **BAM predicts `Σm_ν ≈ 59–65 meV`** — pinned near the normal-ordering
floor (`√Δm²_21 + √Δm²_31 ≈ 58.7 meV`), the light scale keeping it out of
the quasi-degenerate regime. The inverted-ordering floor (the contrast)
is `≈ 99 meV`.

## Cosmological comparison (the sharp, topical test)

| | Σm_ν (meV) |
|---|---|
| **BAM (normal, light scale)** | ≈ 59–65 |
| inverted-ordering floor (contrast) | ≈ 99 |
| Planck 2018 + BAO (95%) | < 120 |
| DESI DR1 + CMB (95%) | < 72 |
| DESI DR2 + CMB (95%, tightest) | ≲ 60–64 |

The BAM prediction sits exactly where cosmology is now probing — the
normal-ordering floor. It is **comfortably below Planck**, **just inside
DESI DR1 + CMB**, and **right at the DESI DR2 + CMB frontier**.

**Sharp falsifiers:**
  - a robust cosmological `Σm_ν < 58.7 meV` (below the NO floor) would
    exclude normal ordering ⟹ BAM fails (such a result would also be in
    tension with the oscillation `Δm²` themselves — a deep consistency
    test);
  - a quasi-degenerate detection (`Σm_ν ≳ 100 meV`) would contradict the
    BAM light scale.

## One spectrum, two observables

`Σm_ν` (~59–65 meV) and the 0νββ effective mass `m_ββ` (≲ 8 meV, PR #95)
are both consequences of the **same** light, normal-ordered, Majorana
spectrum (with anarchic phases) — a joint, cross-checkable prediction.

## Tests

| # | test | finding |
|---|---|---|
| T1 | setup | `Σm_ν = m1+m2+m3`; needs ordering + scale + `Δm²` |
| T2 | NO floor | normal ordering ⟹ floor ≈ 58.7 meV (IO floor ≈ 99 meV) |
| T3 | light scale | `Σm_ν ≈ 59–65 meV` (pinned near the floor) |
| T4 | cosmology | below Planck (<120), at the DESI DR2+CMB frontier (~60–64) |
| T5 | falsifiable | `< 58.7` ⟹ NO excluded; `≳ 100` ⟹ not light |
| T6 | consistency | one spectrum: `Σm_ν` + `m_ββ` (PR #95) cross-checkable |
| T7 | honest scope | pinned near the NO floor; exact value the absolute-scale residual |
| T8 | assessment | `SIGMA_MNU_AT_NORMAL_ORDERING_FLOOR_60MEV` |

## Established and open

  - **Established (BAM-native):** the light, normal-ordered spectrum gives
    `Σm_ν ≈ 59–65 meV` — at the NO floor, below the IO floor (~99 meV),
    consistent with Planck and at the DESI-DR2 + CMB frontier; a
    near-term-falsifiable target, from the same spectrum as the PR #95
    0νββ prediction.

  - **Open:** the exact `Σm_ν` — a narrow band (59–65 meV) because the
    lightest neutrino mass is unmeasured; the absolute scale is the PR #90
    residual (the BAM light scale keeps the band narrow).

## Run

```
python -m experiments.closure_ledger.cosmological_sigma_mnu_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_cosmological_sigma_mnu_probe/`.
Expected verdict: `SIGMA_MNU_AT_NORMAL_ORDERING_FLOOR_60MEV`, 8/8 PASS.
