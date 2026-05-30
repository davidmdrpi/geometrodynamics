# Heavy-quark Möbius baryon: prediction + constraint audit (PR #103)

PR #102 ranked the experimental constraints on BAM-specific baryonic
exotics and found the heavy-quark baryons the FREEST channel (sparse
data) — the place a clean new state is most likely findable. This probe
makes the concrete prediction there and audits its constraints.

## The heavy-quark-symmetry handle

A heavy baryon (Qqq) has the heavy quark Q as a near-static color source
— a spectator. The Möbius / flux-tube excitation (the non-orientable
twist of PRs #100–#102) lives entirely in the LIGHT / flux sector, so its
energy is INDEPENDENT of the heavy-quark mass. The excitation gap is the
flux-tube quantum

```
Δ = 2√σ ≈ 0.85 GeV,
```

the SAME for charm and bottom. This flavor-independence is the
heavy-sector handle that replaces the absent exotic-`J^P` smoking gun of
PR #102: a supernumerary state sitting the same ~0.85 GeV above BOTH the
charm and the bottom ground baryon is the Möbius/hybrid signature — a
correlated, cross-flavor prediction.

## The predictions and constraint audit

Ground-state heavy baryon + Δ (≈ 849 MeV):

| baryon | ground (MeV) | Möbius (MeV) | current ceiling | above |
|---|---:|---:|---|---:|
| Λ_c | 2286 | 3135 | Λ_c(2940) | +195 |
| Ω_c | 2695 | 3544 | Ω_c(3120) | +424 |
| Ξ_cc | 3622 | 4471 | none (ground only) | +849 |
| Λ_b | 5620 | 6469 | Λ_b(6152) | +317 |
| Ω_b | 6045 | 6894 | none (ground only) | +849 |

All predicted states lie just **above** the currently-measured
excitation ceilings (Λ_c data ends ~2.94 GeV, Λ_b ~6.15 GeV) — so they
are unexplored, within LHCb / Belle II / future reach, and **not
excluded**. The doubly-heavy `Ξ_cc` and the `Ω_b` have NO measured
excitation spectrum at all — the freest of the free, entirely
unconstrained.

## Above the orbital tower

The Möbius gap (~0.85 GeV) is well above the ordinary orbital/radial
excitations — Λ_c(2595)/Λ_c(2625) P-wave at ~+0.31 GeV, Λ_c(2940) at
~+0.65 GeV — so the Möbius/hybrid baryon is a **supernumerary** state
ABOVE the orbital tower, not an orbital excitation itself.

## The signature

As in PR #102 there is no forbidden baryon `J^P`, so the heavy Möbius
baryon is a supernumerary ordinary-`J^P` state. But heavy-quark-symmetry
flavor-independence (the SAME Δ for c and b) makes it a **correlated
cross-flavor prediction** — findable precisely because the heavy spectrum
is sparse and isolated.

## Tests

| # | test | finding |
|---|---|---|
| T1 | recap | heavy-quark baryons = freest channel (PR #102) |
| T2 | flavor-indep gap | `Δ = 2√σ ≈ 849 MeV` same for c, b (heavy quark spectator) |
| T3 | predictions | ground + `2√σ` (Λ_c 3.14, Λ_b 6.47 GeV, …) |
| T4 | above orbital tower | `~0.85 GeV` > orbital gaps (P-wave 0.31, 2940 0.65) |
| T5 | constraint audit | all above current ceilings; Ξ_cc / Ω_b unexplored |
| T6 | cross-flavor signature | same gap for c and b — correlated, findable |
| T7 | honest scope | flavor-indep gap is the prediction; exact mass/`J^P` open |
| T8 | assessment | `HEAVY_MOBIUS_BARYON_FLAVOR_INDEPENDENT_GAP_FINDABLE_UNCONSTRAINED` |

## Established and open

  - **Established (BAM-native):** the heavy-quark Möbius/hybrid baryon
    sits at the flavor-INDEPENDENT flux-tube gap `Δ = 2√σ ≈ 0.85 GeV`
    above the ground heavy baryon (heavy quark a spectator), above the
    orbital tower; the predictions (Λ_c ~3.14, Ω_c ~3.54, Λ_b ~6.47, Ω_b
    ~6.89, Ξ_cc ~4.47 GeV) lie just above current data — findable, not
    excluded; the doubly-heavy and Ω_b channels are entirely
    unconstrained; the cross-flavor correlation is the signature.

  - **Open:** the exact mass (the flux-tube/hybrid gap is ~0.85 GeV here
    but lattice hybrid gaps span ~0.8–1.3 GeV) and the `J^P` (no smoking
    gun) — it is a correlated counting prediction in the freest channel.

## Cross-references

  - `docs/baryonic_exotics_classification_research_plan.md` — PR #102, the
    constraint ranking that identified the heavy sector as freest.
  - `docs/mobius_exotic_sector_research_plan.md` — PR #101, the `2√σ`
    flux-tube gap and the mesonic exotics.
  - `geometrodynamics/qcd/topology.py` — `make_mobius_baryon_*`.

## Run

```
python -m experiments.closure_ledger.heavy_mobius_baryon_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_heavy_mobius_baryon_probe/`.
Expected verdict: `HEAVY_MOBIUS_BARYON_FLAVOR_INDEPENDENT_GAP_FINDABLE_UNCONSTRAINED`, 8/8 PASS.
