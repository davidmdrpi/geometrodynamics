# Heavy-quark Möbius baryon: prediction + constraint audit (PR #103)

**Run:** 2026-05-30T16:16:34+00:00

Makes the concrete BAM-specific baryonic-exotic prediction in the FREEST channel (PR #102's heavy-quark baryons). **Key handle:** by heavy-quark symmetry the heavy quark is a spectator, so the Möbius/flux excitation gap `Δ = 2√σ ≈ 0.85 GeV` is **flavor-independent** (same for charm and bottom) — the cross-flavor signature that replaces the absent exotic-`J^P` smoking gun. The predictions land just **above** current data — findable, not excluded.

- **Identification**: heavy-quark Möbius/hybrid baryon at the flavor-independent flux-tube gap 2√σ ≈ 0.85 GeV above the ground heavy baryon (Λ_c ~3.14, Λ_b ~6.47 GeV, …), above the orbital tower and just above current data — findable, not excluded; doubly-heavy / Ω_b entirely unconstrained
- **Gap**: 2√σ ≈ 849 MeV, flavor-independent (heavy quark spectator)
- **Predictions**: Λ_c ~3.14, Ω_c ~3.54, Λ_b ~6.47, Ω_b ~6.89, Ξ_cc ~4.47 GeV
- **Signature**: supernumerary state at the SAME gap for c and b (cross-flavor)
- **Audit**: just above current ceilings ⟹ findable; Ξ_cc / Ω_b entirely unexplored
- **Open**: exact mass (lattice hybrid gap 0.8–1.3 GeV); J^P (no smoking gun)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_recap_freest_channel` | heavy-quark baryons = freest channel (PR #102) | **PASS** |
| T2 | `T2_flavor_independent_flux_gap` | heavy-quark symmetry ⟹ gap 2√σ flavor-independent (c=b) | **PASS** |
| T3 | `T3_heavy_mobius_predictions` | predictions: ground + 2√σ (Λ_c 3.14, Λ_b 6.47 GeV, …) | **PASS** |
| T4 | `T4_above_orbital_tower` | gap ~0.85 GeV above the orbital tower (supernumerary) | **PASS** |
| T5 | `T5_constraint_audit` | all above current ceilings ⟹ findable; Ξ_cc/Ω_b unexplored | **PASS** |
| T6 | `T6_cross_flavor_signature` | cross-flavor: same gap for c and b (correlated signature) | **PASS** |
| T7 | `T7_honest_scope` | flavor-indep gap is the prediction; exact mass/J^P open | **PASS** |
| T8 | `T8_assessment` | HEAVY_MOBIUS_BARYON_FLAVOR_INDEPENDENT_GAP_FINDABLE_UNCONSTRAINED | **PASS** |

## Predictions and constraint audit

| baryon | ground (MeV) | Möbius (MeV) | current ceiling | above ceiling |
|---|---:|---:|---|---:|
| Lambda_c | 2286 | 3135 | Λ_c(2940) | +195 MeV |
| Omega_c | 2695 | 3544 | Ω_c(3120) | +424 MeV |
| Xi_cc | 3622 | 4471 | none (ground only) | +849 MeV |
| Lambda_b | 5620 | 6469 | Λ_b(6152) | +317 MeV |
| Omega_b | 6045 | 6894 | none (ground only) | +849 MeV |

All predictions sit above the current excitation ceilings — findable, not excluded. The doubly-heavy `Ξ_cc` and `Ω_b` have no measured excitations at all: entirely unconstrained.

## Verdict

**HEAVY_MOBIUS_BARYON_FLAVOR_INDEPENDENT_GAP_FINDABLE_UNCONSTRAINED.** THE HEAVY-QUARK MÖBIUS BARYON SITS AT A FLAVOR-INDEPENDENT FLUX-TUBE GAP, JUST ABOVE CURRENT DATA — FINDABLE AND UNCONSTRAINED. PR #102 found the heavy-quark baryons the freest channel for BAM-specific baryonic exotics. This probe makes the concrete prediction there and audits its constraints.

THE HEAVY-QUARK-SYMMETRY HANDLE. A heavy baryon (Qqq) has the heavy quark as a near-static color source — a spectator. The Möbius / flux-tube excitation (the non-orientable twist of PRs #100–#102) lives entirely in the LIGHT / flux sector, so its energy is INDEPENDENT of the heavy-quark mass: the gap is the flux-tube quantum Δ = 2√σ ≈ 0.85 GeV, the SAME for charm and bottom. This flavor-independence is the heavy-sector handle that replaces the absent exotic-J^P smoking gun of PR #102 — a supernumerary state sitting the same ~0.85 GeV above BOTH the charm and bottom ground baryon is the Möbius/hybrid signature.

THE PREDICTIONS. Ground heavy baryon + Δ: Λ_c ~3.14, Ω_c ~3.54, Λ_b ~6.47, Ω_b ~6.89, Ξ_cc ~4.47 GeV. The Möbius gap (~0.85 GeV) is well above the ordinary orbital/radial excitations (Λ_c P-wave ~+0.31 GeV, Λ_c(2940) ~+0.65 GeV), so the Möbius/hybrid baryon is a SUPERNUMERARY state ABOVE the orbital tower, not an orbital excitation.

CONSTRAINT AUDIT. All predicted states lie just above the currently-measured excitation ceilings (Λ_c data ends ~2.94 GeV, Λ_b ~6.15 GeV) — unexplored, within LHCb / Belle II / future reach, and NOT excluded. The doubly-heavy Ξ_cc and the Ω_b have no measured excitation spectrum at all — BAM's prediction there is entirely unconstrained, the freest of the free. The signature is correlated, not exotic-J^P: as in PR #102 there is no forbidden baryon J^P, so the heavy Möbius baryon is a supernumerary ordinary-J^P state — but the heavy-quark-symmetry flavor-independence (same Δ for c and b) makes it a correlated cross-flavor prediction, findable precisely because the heavy spectrum is sparse and isolated.

HONEST SCOPE. ESTABLISHED (BAM-native): the heavy-quark Möbius/hybrid baryon sits at the flavor-INDEPENDENT flux-tube gap Δ = 2√σ ≈ 0.85 GeV above the ground heavy baryon, above the orbital tower; the predictions lie just above current data — findable, not excluded; the doubly-heavy and Ω_b channels are entirely unconstrained; the cross-flavor correlation is the heavy-sector signature. NOT established: the exact mass (the flux-tube/hybrid gap is ~0.85 GeV here but lattice hybrid gaps span ~0.8–1.3 GeV) or the J^P (no smoking gun) — it is a correlated counting prediction in the freest channel.

## What this leaves open

- **The exact mass** — the flux-tube/hybrid gap is `2√σ ≈ 0.85 GeV` here, but lattice hybrid gaps span ~0.8–1.3 GeV; the prediction is the flavor-independent *gap*, not a sub-percent mass.
- **The `J^P`** — no smoking gun (PR #102); a supernumerary ordinary-`J^P` state. The cross-flavor correlation (same gap for c and b) is the testable signature.
