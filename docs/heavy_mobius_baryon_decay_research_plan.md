# Heavy-quark Möbius baryon: decay channels + search strategy (PR #109)

PR #103 predicted the heavy Möbius/hybrid baryon at the flavor-INDEPENDENT
flux-tube gap `Δ = 2√σ ≈ 0.85 GeV` above the ground heavy baryon (Λ_c
~3135, Λ_b ~6469 MeV, …), findable and unconstrained in the sparse heavy
sector. But it left the `J^P` open and named no way to TELL the state apart
from an ordinary radial/orbital excitation once a bump is seen. This probe
completes that: **how the Möbius baryon decays, and where and how to search
for it.**

## The decay mechanism: twist-unwinding

The Möbius excitation is the non-orientable (orientation `−1`) flux-tube
sector of PRs #100–#102; the ground heavy baryon is orientable (`+1`). To
reach the ground state the half-twist must **unwind**, shedding the stored
`2√σ ≈ 0.85 GeV` as light, isoscalar hadrons. This is the heavy-baryon
analog of a gluonic hybrid de-exciting — the flux/topological degree of
freedom is radiated, not the heavy quark (a spectator).

## The hybrid selection rule (the falsifiable handle)

This is the prediction that distinguishes a Möbius/hybrid baryon from an
ordinary radial excitation. The flux-tube model's hybrid selection rule
forbids a hybrid from decaying into two ground-state (both-S-wave) hadrons —
the symmetric flux configuration has zero overlap with the antisymmetric
(excited) tube. Inherited by the Möbius BARYON:

  - **SUPPRESSED:** Möbius → (ground heavy baryon) + (single S-wave π). The
    naive, most phase-space-favored channel is the one the topology
    suppresses.
  - **PREFERRED:** `Σ_Q π` (the light diquark in its spin-1 config), the
    isoscalar S-wave dipion `Λ_Q(ππ)_S` (the twist unwinding coherently,
    like `ψ(2S) → J/ψ ππ`), and decays to a P-wave heavy baryon + π.

An ordinary radial excitation does the OPPOSITE — single π to the ground
state. So the branching **pattern**, not the mass, is the test.

## The open channels and Q-values

`Möbius_c = 3135 MeV`, `Möbius_b = 6469 MeV`. Release energies:

| channel | charm Q (MeV) | bottom Q (MeV) | role |
|---|---:|---:|---|
| `Λ_Q π⁺π⁻` (isoscalar S-wave) | 569 | 569 | twist-unwinding — PREFERRED |
| `Σ_Q π` | 542 | 515 | spin-1 diquark — PREFERRED |
| `Σ_Q* π` | 477 | 496 | spin-1 diquark — PREFERRED |
| `Λ_Q η` | 301 | 301 | isoscalar — allowed |
| `D N` / `B N` | 332 | 251 | open-flavor — threshold-sensitive |

## The cross-flavor Q-match (the clincher)

Because BOTH the gap (`2√σ`) and the light-meson thresholds (`ππ`, `η`) are
flavor-independent, the all-light release energies are **identical** for
charm and bottom: `Λ_Q ππ` at `Q = 569 MeV` and `Λ_Q η` at `Q = 301 MeV`,
the same dipion mass spectrum scaled only by the heavy-baryon recoil.
Observing the same Q-value structure above both the charm and the bottom
ground baryon — with the hybrid branching pattern — is the Möbius
signature. The `Σ_Q π` channels differ only by the `Σ_Q − ground`
hyperfine splitting (167 MeV for c, 194 for b), itself a checkable offset.

## Width and difficulty (honest)

With several open channels and `Q ≈ 0.5 GeV`, the Möbius baryon is NOT
narrow: lattice and flux-tube hybrid widths run `~tens–150 MeV`. The
selection-rule suppression of the single-pion-to-ground channel partially
protects it (narrowing it relative to naive phase space), but the state is
still broad — best resolved in amplitude (Dalitz) analyses, not as a sharp
peak. This is the honest cost of sitting `0.85 GeV` up.

## Search strategy

  - **LHCb** (the workhorse): amplitude analyses of `Λ_c⁺π⁺π⁻` and
    `Λ_b⁰π⁺π⁻` (prompt and from b-hadron decays), `Σ_Q π`, and the
    open-flavor `D⁰p` / `B⁻p` final states; the doubly-heavy `Ξ_cc⁺⁺π` and
    `Ω_b` channels are wide open.
  - **Belle II:** `e⁺e⁻ → cc̄`, `Λ_c⁺π⁺π⁻` and `Σ_c π` near threshold.
  - **Discriminators:** (i) the suppressed single-π-to-ground branch; (ii)
    the isoscalar S-wave dipion peaked at high `m(ππ)`; (iii) the
    cross-flavor Q-match (`569` / `301 MeV` for `ππ` / `η` above both c and
    b).

## Tests

| # | test | finding |
|---|---|---|
| T1 | recap | #103 prediction; this probe = decay channels + search |
| T2 | mechanism | twist-unwinding: Möbius (−1) → ground (+1) sheds `2√σ` as light hadrons |
| T3 | channels + Q | `Λ_Q ππ` 569, `Σ_Q π` 542/515, `Σ_Q* π` 477/496, `Λ_Q η` 301 MeV |
| T4 | selection rule | single-π-to-ground SUPPRESSED; `Σ_Q π` / dipion / P-wave+π PREFERRED (≠ radial) |
| T5 | cross-flavor Q | `Λ_Q ππ` 569, `Λ_Q η` 301 MeV identical for c and b |
| T6 | width | broad (~tens–150 MeV); best in amplitude analyses |
| T7 | search/scope | LHCb / Belle II channels; pattern + Q-structure are the predictions |
| T8 | assessment | `HEAVY_MOBIUS_BARYON_TWIST_UNWINDING_HYBRID_SELECTION_RULE_CROSS_FLAVOR_FINDABLE` |

## Established and open

  - **Established (BAM-native):** the decay proceeds by twist-unwinding
    (non-orientable → orientable), so the Möbius baryon inherits the hybrid
    selection rule — single-S-wave-π-to-ground SUPPRESSED, `Σ_Q π` /
    isoscalar dipion / P-wave+π PREFERRED — which TELLS it apart from a
    radial excitation; and the all-light Q-values (`Λ_Q ππ` 569, `Λ_Q η`
    301 MeV) are cross-flavor identical, a correlated signature, with
    concrete LHCb / Belle II search channels.

  - **Open:** absolute branching fractions and the total width (need the
    flux-tube decay amplitudes, not computed here; the state is expected
    broad), and the `J^P` (open since PR #102). The predictions are the
    branching PATTERN and the Q-structure, not partial rates.

## Cross-references

  - `docs/heavy_mobius_baryon_research_plan.md` — PR #103, the mass
    prediction (ground + `2√σ`) this probe gives the decays/search for.
  - `docs/baryonic_exotics_classification_research_plan.md` — PR #102, the
    constraint ranking (heavy sector = freest; no exotic-`J^P` smoking gun).
  - `docs/mobius_exotic_sector_research_plan.md` — PR #101, the `2√σ`
    flux-tube gap and the hybrid/Möbius mesonic exotics.
  - `geometrodynamics/qcd/topology.py` — `make_mobius_baryon_*`.

## Run

```
python -m experiments.closure_ledger.heavy_mobius_baryon_decay_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_heavy_mobius_baryon_decay_probe/`.
Expected verdict:
`HEAVY_MOBIUS_BARYON_TWIST_UNWINDING_HYBRID_SELECTION_RULE_CROSS_FLAVOR_FINDABLE`,
8/8 PASS.
