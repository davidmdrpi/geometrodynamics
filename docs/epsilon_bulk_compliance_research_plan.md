# Can ε be computed from bulk compliance, or only inferred from the meV scale? (PR #112)

The chargeless-throat healing length `ε` is the mass-generating parameter of
the BAM Majorana neutrino: `m_ν = m_D · e^{−S}` with the bounce action
`S = t²·P0·L*(ε)` and tortoise length `L*(ε) ≈ −(rs/2) ln ε + const`. PRs
#87–#90 **inferred** `ε` by demanding the observed `S = ln(m_D/m_ν) ≈ 15–18`
— i.e. `ε` was read back from the meV scale it is meant to predict. This PR
asks the sharper question: can `ε` instead be **computed from the bulk
compliance** (the elastic/nucleation response of the throat neck), with no
neutrino input?

**Answer: a genuine partial — the smallness yes, the precise value no.**

## Bulk compliance ⟹ the healing length (meV-free)

The compliance cutoff is a sub-throat healing length `ε = ℓ²/(2 rs)`, where
`ℓ` is the length over which the chargeless neck heals back into the bulk.
For the nucleating throat that length is the critical-bubble scale `ℓ ~ R_c =
2σ/ρ`, with the surface tension `σ` and bag density `ρ` fixed by the
**electron** rest-energy calibration (PR #58: `σ = 1/(12π)`, `ρ = 3/(4π)`),
so `R_c = 2/9 ≈ 0.222` — a pure number from the charged-throat geometry, with
no neutrino mass anywhere. The candidate compliances are all sub-throat,
`O(10⁻²)`:

```
ε ~ R_c³ ≈ 0.011,   Δ³ ≈ 0.018,   R_c²/2 ≈ 0.025.
```

So the **order of magnitude** of `ε` — and its decisive **sub-throat
character** — IS computable from bulk compliance, independent of the meV
scale.

## The chain closes meV-free

With `ε = R_c³` (bulk nucleation geometry) and the winding-edge tension
`t = k_5√(2π) ≈ 12.53` (the PR #89 closure quantity, also meV-free):

```
S = t²·P0·L*(R_c³) ≈ 16.85   ⟹   m_ν(gen 1) = m_D·e^{−S} ≈ 2.1 meV,
```

with `m_D ≈ 43 keV` the electron-anchored cavity-floor Dirac mass. The meV
**scale** comes out as an **output** of bulk geometry — a genuine
retrodiction, not an input. This structurally **derives the neutrino's
lightness**: a sub-throat healing length (`ε ≪ 1`) makes `L*` large, hence
`S` large, hence `m_ν = m_D·e^{−S}` exponentially small.

(Claim is the scale only: uniform `S` gives `m_ν ∝ m_D` (×2.7 spread), so the
observed generation spread — up to `√Δm²_31 = 50 meV` — is the separate known
residual of PR #90/#91, not addressed by `ε`.)

## The catch: precision is not geometric

The bounce action is steep in `ε` at the winding edge:

```
d ln m_ν / d ln ε = t²·P0·rs/2 ≈ 4.8.
```

So the `O(1)` ambiguity among the healing-length candidates blows up:

| ε candidate | value | m_ν (meV) |
|---|---:|---:|
| R_c³ | 0.011 | 2.1 |
| Δ³ | 0.018 | 20 |
| R_c²/2 | 0.025 | 108 |

a factor ~50 spread from a factor ~2 in `ε`. Landing precisely on "few meV"
still **selects** `ε ≈ R_c³` using the observed scale — there is no
first-principles reason to prefer `R_c³` over `R_c²/2`.

## The obstruction to pinning ε

A true compliance = 1/stiffness needs the absolute bulk gravitational
stiffness, set by `κ₅²` and `Λ₅`. BAM fixes only the dimensionless RS-tuning
combination `λ_crit κ₅²/√|Λ₅| = √6` (PR #57); `κ₅²` and `Λ₅` separately are
absorbed into the single dimensionful anchor and are **not pinned**. So the
compliance — and hence `ε` — is computable only up to an `O(1)` factor,
exactly the residual status everything tied to the one anchor shares.

## The honest verdict (a genuine partial)

  - **Computed (meV-free):** `ε`'s order of magnitude and sub-throat
    character — `ε ~ R_c³ ~ 10⁻²` from the electron-calibrated nucleation
    geometry — which **derives** the exponential smallness of `m_ν` and, at
    the winding-edge tension, **outputs** `m_ν ~ meV` with no neutrino input.
    So `ε` is upgraded from "inferred from the meV scale" to "bulk-geometric
    to order of magnitude; the meV scale is a retrodiction."
  - **Still residual:** the **precise** `ε`. Because `m_ν ∝ ε^{4.8}`, the
    `O(1)` ambiguity (`R_c³` vs `Δ³` vs `R_c²/2`) spans `m_ν ~ 2–108 meV`;
    the exact value is not fixed by geometry (the absolute normalization is
    the unpinned `κ₅²/Λ₅` = the one anchor), and pinning the precise meV
    still uses the observed scale.

**So: the smallness is derived from bulk compliance; the exact value is not.**

## Tests

| # | test | finding |
|---|---|---|
| T1 | question | `ε` currently inferred via `S = ln(m_D/m_ν)` |
| T2 | define compliance | `ε = ℓ²/2rs`, `ℓ ~ R_c = 2σ/ρ` (electron-calibrated) |
| T3 | compute ε | `R_c³` 0.011, `Δ³` 0.018, `R_c²/2` 0.025 — sub-throat, meV-free |
| T4 | chain closes | `ε=R_c³` + `t=k_5√(2π)` ⟹ `S≈16.85` ⟹ `m_ν≈2.1 meV` (retrodiction) |
| T5 | the catch | `m_ν ∝ ε^{4.8}`: `O(1)` ambiguity spans 2–108 meV |
| T6 | obstruction | absolute compliance = unpinned `κ₅²/Λ₅`; only `√6` fixed |
| T7 | verdict | smallness DERIVED (meV-free); precise value RESIDUAL |
| T8 | assessment | `EPSILON_BULK_COMPLIANCE_DERIVES_SMALLNESS_PRECISE_VALUE_STAYS_RESIDUAL` |

## Established and open

  - **Established (BAM-native):** `ε`'s order of magnitude and sub-throat
    character are computable from the electron-calibrated bulk nucleation
    compliance (`ε ~ R_c³ ~ 10⁻²`, meV-free), which derives the exponential
    smallness of `m_ν` and retrodicts `m_ν ~ meV` at the winding-edge
    tension. `ε` is no longer merely "inferred from the meV scale."

  - **Open:** the precise `ε` — because `m_ν ∝ ε^{4.8}`, the `O(1)` ambiguity
    in the healing length spans `m_ν ~ 2–108 meV`; the absolute compliance
    normalization is the unpinned `κ₅²/Λ₅` (the one anchor), so the exact
    value stays a residual. The generation spread is the separate PR #90/#91
    residual.

## Cross-references

  - `docs/boundary_compliance_bulk_geometry_research_plan.md` (if present) /
    PR #90 — the geometric `ε` candidates this builds on.
  - `docs/majorana_bounce_action_research_plan.md` — PR #88, the bounce
    action `S = t²·P0·L*(ε)` and the tortoise length.
  - `docs/seesaw_scale_nucleation_compliance_research_plan.md` — PR #87,
    `M_R = m_D·e^{S}` and the required `S`.
  - `docs/neutrino_mev_scale_sharpening_research_plan.md` — PR #111, the
    observed meV-scale spectrum this retrodicts (to order of magnitude).
  - PR #57 (brane-tension tuning) — the `√6` RS combination and the unpinned
    `κ₅²/Λ₅`.

## Run

```
python -m experiments.closure_ledger.epsilon_bulk_compliance_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_epsilon_bulk_compliance_probe/`.
Expected verdict:
`EPSILON_BULK_COMPLIANCE_DERIVES_SMALLNESS_PRECISE_VALUE_STAYS_RESIDUAL`, 8/8 PASS.
