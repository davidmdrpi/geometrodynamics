# BAM non-orientable hadron sector: experimental note (PR #110)

PRs #100–#109 built BAM's non-orientable (Möbius / closed-flux-loop) hadron
predictions across four sub-sectors — mesonic `1⁻⁺` hybrids, glueballs, and
heavy Möbius baryons with their decays. This PR COMPILES them into a single
compact **experimental note** of the kind an LHCb / Belle II / BESIII
analyst can read off: predicted masses, Q-values, preferred/suppressed
modes, and analysis handles. Every number is recomputed from the one input
(`√σ`) so the note cannot drift from the probes it summarises.

The compiled note is committed standalone as
`docs/bam_nonorientable_experimental_note.md`.

## The single input

The whole sector is set by one QCD scale, the string tension `√σ ≈ 0.424
GeV` (the B4 confinement anchor). The flux-tube excitation quantum is `Δ =
2√σ ≈ 0.849 GeV`; the lowest closed-flux-loop (glueball) level is `√(4πσ) ≈
1.50 GeV`. No new parameters — the note is a pushforward of `√σ`.

## What the note collects

  1. **Mesonic `1⁻⁺` hybrids (PR #101):** ground meson + `2√σ` — π₁ ≈ 1.62,
     η₁ ≈ 1.85 GeV — matched to the observed `π₁(1600)` / `η₁(1855)`; the
     exotic `1⁻⁺` (forbidden to ordinary qq̄) is the Möbius/twist marker and
     the one place the sector has a smoking-gun `J^PC`.
  2. **Glueballs as closed flux loops (PR #100):** `0⁺⁺` ground `√(4πσ) ≈
     1.50 GeV` (benchmarks the lattice `0⁺⁺` `√σ` scale to ~13%); unobserved
     — the freest channel.
  3. **Heavy Möbius baryons — masses (PR #103):** ground heavy baryon +
     `2√σ` — Λ_c ~3135, Ω_c ~3544, Ξ_cc ~4471, Λ_b ~6469, Ω_b ~6894 MeV —
     supernumerary, above the orbital tower, just above current data.
  4. **Heavy Möbius baryons — decays (PR #109):** twist-unwinding ⟹ the
     hybrid selection rule (single-S-wave-π-to-ground SUPPRESSED; `Σ_Q π` /
     isoscalar dipion / P-wave+π PREFERRED) and the cross-flavor Q-match
     (`Λ_Q ππ` 569, `Λ_Q η` 301 MeV identical for c and b).

## Analysis handles (the note's payload)

  - **Cross-flavor Q-match** — the all-light release energies are
    flavor-independent: the SAME dipion spectrum (Q = 569 MeV) and η recoil
    (Q = 301 MeV) above the charm and the bottom ground baryon.
  - **Hybrid selection rule** — the naive single-S-wave-π-to-ground
    transition is suppressed; a radial excitation would show the opposite.
    The branching PATTERN is the discriminator.
  - **Isoscalar dipion** — the preferred `Λ_Q(ππ)_S` channel peaks at high
    `m(ππ)` (coherent twist unwinding).
  - **Broad ⟹ amplitude analyses** — widths ~tens–150 MeV; resolve in
    Dalitz / amplitude fits, not as sharp peaks.
  - **Exotic `J^PC` where available** — the mesonic sector has the `1⁻⁺`
    smoking gun; the baryon sector (supernumerary ordinary-`J^P`) leans on
    the cross-flavor correlation instead.

## Tests

| # | test | finding |
|---|---|---|
| T1 | scope | consolidate the #100–#109 non-orientable sector into one note |
| T2 | single input | `√σ` 424, `2√σ` 849, `√(4πσ)` 1504 MeV |
| T3 | mesonic `1⁻⁺` | π₁ ~1.62, η₁ ~1.85 GeV — matched |
| T4 | glueball | `0⁺⁺` ground `√(4πσ)` ~1.50 GeV |
| T5 | heavy masses | Λ_c 3135, Ω_c 3544, Ξ_cc 4471, Λ_b 6469, Ω_b 6894 MeV |
| T6 | decays | channels + Q; cross-flavor match (569/301) |
| T7 | modes/handles | hybrid selection rule + analysis handles |
| T8 | assessment | `NONORIENTABLE_SECTOR_EXPERIMENTAL_NOTE_COMPILED` |

## Established and open

  - **Established (BAM-native):** the four non-orientable sub-sectors are
    compiled into one compact, internally-consistent experimental note
    (masses, Q-values, preferred/suppressed modes, analysis handles), all
    pushed forward from the single input `√σ`.

  - **Open (unchanged):** exact masses within the lattice hybrid-gap band
    (~0.8–1.3 GeV); absolute branching fractions and total widths; the
    baryon `J^P`. The note is a reference card, not new physics — it carries
    the established content and the open items unchanged.

## Cross-references

  - `docs/bam_nonorientable_experimental_note.md` — the compiled note (the
    deliverable).
  - `docs/mobius_exotic_sector_research_plan.md` — PR #101, the mesonic
    `1⁻⁺` hybrids and the `2√σ` gap.
  - `docs/glueball_closed_flux_loop_research_plan.md` — PR #100, the
    `√(4πσ)` glueball ground.
  - `docs/heavy_mobius_baryon_research_plan.md` — PR #103, the heavy Möbius
    baryon masses.
  - `docs/heavy_mobius_baryon_decay_research_plan.md` — PR #109, the decays
    and the hybrid selection rule.

## Run

```
python -m experiments.closure_ledger.nonorientable_experimental_note_probe
```

Writes `probe.json`, `probe.md`, and the standalone `experimental_note.md`
under
`experiments/closure_ledger/runs/<UTC timestamp>_nonorientable_experimental_note_probe/`.
Expected verdict: `NONORIENTABLE_SECTOR_EXPERIMENTAL_NOTE_COMPILED`, 8/8 PASS.
