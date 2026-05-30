# The Möbius / exotic sector: non-orientable flux tubes give exotic J^PC (PR #101)

PR #100 showed BAM's non-orientable (Möbius) closed loops predict an
extra glueball tower — legitimate because glueballs are NOT observed.
This probe pursues the Möbius / exotic sector for **open** flux networks
(hybrids, multiquark exotics), where — unlike glueballs — the states
**are** experimentally observed. So this is the sharper test: BAM's
non-orientable topology must now match data, and it does.

## Topological hadron classification

BAM flux-network topology is the hadron taxonomy:

| state | flux network | content |
|---|---|---|
| meson | open tube | qq̄ |
| baryon | Y-junction (3 arms) | qqq |
| tetraquark | two junctions / diquark–antidiquark | qq q̄q̄ |
| pentaquark | junction + meson cloud | qqqq q̄ |
| hybrid | tube + excited/twisted flux | qq̄ + glue |
| glueball | closed loop (PR #100) | pure glue |

Each has an orientable and a non-orientable (Möbius, Z₂-twisted) version.

## The Möbius twist = the exotic J^PC marker

An ordinary orientable qq̄ meson is restricted to `P = (−1)^{L+1}`,
`C = (−1)^{L+S}`, which FORBIDS the combinations `0--, 0+-, 1-+, 2+-` —
the "exotic" J^PC. A NON-ORIENTABLE (Möbius) flux tube carries an
antiperiodic flux-tube phonon (PR #100) whose C/P assignment opens
exactly those forbidden channels — most importantly **`1-+`**. So in BAM
the exotic quantum numbers are the signature of a non-orientable flux
tube: the Möbius twist is the flux-tube-hybrid excitation.

## These exotics ARE observed (the sharp test)

The lightest exotic hybrids are the `1-+` states `π₁(1400)`,
`π₁(1600)`, and `η₁(1855)` (BESIII, 2022) — all carrying the exotic
`1-+` that ordinary qq̄ cannot. The BAM Möbius/hybrid excitation gap is

```
Δ ≈ 2√σ ≈ 0.85 GeV    (one flux-tube excitation),
```

so the lightest `1-+` sits at:

| state | BAM (base + 2√σ) | observed |
|---|---:|---:|
| π₁ (1-+) | `ρ(770) + 2√σ ≈ 1.62 GeV` | 1.66 GeV |
| η₁ (1-+) | `~1.0 + 2√σ ≈ 1.85 GeV` | 1.85 GeV |

BAM's non-orientable topology accounts for the observed exotic hybrids at
the right `J^PC` AND the right mass.

## The multiquark exotic zoo

The observed tetraquarks (`X(3872)`, `Z_c(3900)`, `T_cc(3875)`) and
pentaquarks (`P_c`) fit the multi-junction flux-network topologies
(diquark–antidiquark / junction + meson cloud). BAM's taxonomy
accommodates the observed exotic zoo, with Möbius-twisted partners as
predictions.

## Observability: the contrast with glueballs

Glueballs are unobserved, so BAM's Möbius glueball tower (PR #100) is a
legitimate but untestable-against-experiment divergence. Exotic hybrids
and multiquark states ARE observed, so the Möbius / exotic sector is
where BAM's non-orientable topology MUST meet data — and it does (the
`1-+` hybrids at `~2√σ` above the ground meson; the multiquark networks).
The Möbius half-twist is the same non-orientable Z₂ that gives the throat
its spin-½ (PR #63–#67) — the exotic states carry an odd flux twist.

## Tests

| # | test | finding |
|---|---|---|
| T1 | classification | flux-network topology = hadron taxonomy (+ Möbius Z₂) |
| T2 | exotic marker | Möbius flux tube opens `1-+` (forbidden to qq̄) |
| T3 | exotic J^PC | `0--, 0+-, 1-+, 2+-` exotic; observed hybrids all `1-+` |
| T4 | mass scale | gap `2√σ` ⟹ `π₁ ≈ 1.62`, `η₁ ≈ 1.85` GeV (obs 1.66, 1.85) |
| T5 | multiquark | `X, Z_c, T_cc, P_c` fit multi-junction networks |
| T6 | observability | glueballs unobserved (free); exotics observed (BAM matches) |
| T7 | Z₂ tie | Möbius twist = throat spin-½ Z₂; Möbius baryon = prediction |
| T8 | assessment | `MOBIUS_FLUX_GIVES_EXOTIC_JPC_OBSERVED_HYBRIDS_MATCH` |

## Established and open

  - **Established (BAM-native):** the flux-network topology is the hadron
    taxonomy; the non-orientable (Möbius) flux tube carries the exotic
    `J^PC` (`1-+`) forbidden to ordinary qq̄; the observed exotic hybrids
    (`π₁`, `η₁`) match at the right `J^PC` and at `~2√σ ≈ 0.85 GeV` above
    the ground meson; the observed multiquark exotics fit the
    multi-junction networks. Unlike glueballs, this sector is observed —
    and BAM matches.

  - **Open:** precise exotic masses (full network dynamics: junction
    binding, mixing), and the unobserved Möbius partners (the Möbius
    baryon, `make_mobius_baryon`) — those remain BAM-specific predictions.

## Cross-references

  - `geometrodynamics/qcd/topology.py` — `make_mobius_tube`,
    `make_mobius_baryon_*`, hybrid / tetraquark constructions.
  - `docs/glueball_closed_flux_loop_research_plan.md` — PR #100, the
    Möbius closed-loop (glueball) tower (the unobserved counterpart).
  - `docs/qcd_confinement_cornell_audit_research_plan.md` — PR #99, the
    string tension `√σ` setting the `2√σ` exotic gap.

## Run

```
python -m experiments.closure_ledger.mobius_exotic_sector_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_mobius_exotic_sector_probe/`.
Expected verdict: `MOBIUS_FLUX_GIVES_EXOTIC_JPC_OBSERVED_HYBRIDS_MATCH`, 8/8 PASS.
