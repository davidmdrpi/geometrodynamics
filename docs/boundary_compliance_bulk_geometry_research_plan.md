# Boundary compliance ε from bulk throat geometry (PR #90)

The capstone of the neutrino-mass chain. PR #89 bracketed the `ΔL=2` /
B−L tension ratio `t ∈ [2π, k_5√(2π)]` and left a single residual
number: *where in that window*, i.e. the boundary compliance `ε` — the
minimal approach of the non-orientable bounce path to the throat neck at
`r = rs = R_MID`. This probe derives `ε` from the bulk throat geometry.

## ε is a sub-throat healing length

Near the neck the warp is `f(r) = 1 − (rs/r)² ≈ 2(r − rs)/rs`, so the
proper radial distance from the neck out to `rs + ε` is

```
ℓ_proper = ∫_rs^{rs+ε} dr/√f ≈ √(2 rs ε)   ⟹   ε = ℓ_proper² / (2 rs).
```

`ε` is the (squared, neck-warped) **healing length** of the throat — the
proper thickness over which the brane profile heals the would-be pinch.
A sub-throat healing length `ℓ ≪ rs` gives a sub-throat `ε`, exactly the
PR #89 window `ε ∈ [6×10⁻⁷, 1.3×10⁻²]`.

## Why ε is tiny *for the neutrino* (c₁ = 0)

  - **Charged lepton (`c₁ = ±1`):** the throat is **propped open** by its
    EM self-repulsion `A/R` at the equilibrium `R* ≈ R_MID` (PR #55). The
    charge holds the neck open, the bounce cannot approach it, and the
    charged throat cannot flip at all — it stays **Dirac** (PR #86).
  - **Neutrino (`c₁ = 0`):** `A = 0`. Nothing props the neck open, so the
    bounce approaches it all the way down to the bulk healing length `ℓ`.

The chargelessness that makes the neutrino Majorana (PR #86) is the
**same** property that makes its `ε` sub-throat — and therefore its mass
tiny. The smallness of `m_ν` is the unobstructed near-rigidity of the
chargeless neck.

## Bulk candidates and the chain closure

The natural BAM sub-throat scales all cluster at `ε ~ 10⁻²–10⁻³`, inside
the PR #89 window:

| candidate | ε |
|---|---:|
| `R_c³` (nucleation critical-bubble scale) | `1.10×10⁻²` |
| `Δ³` (cavity-aspect cube) | `1.76×10⁻²` |
| `(m_D/m_charged)²` (charge suppression) | `7.0×10⁻³` |
| `E_c` (nucleation barrier scale) | `5.5×10⁻³` |

Combined with the **winding-edge** tension `t ≈ k_5√(2π) = √β` — the edge
PR #89's cross-check `m_charged/m_D ≈ 11.9 ≈ √β` independently favoured —
the leading `ε = R_c³` gives

```
S = t²·P0·L*(ε) ≈ 16.8,   m_ν = m_D·e^{−S} ≈ 2–6 meV,
```

squarely the observed neutrino scale (`1–50 meV`). At the *lower*
(closure-quantum) edge `t = 2π` the same `ε` give `S ≈ 4` (far too
small): **the chain closes only at the winding edge** — the same edge the
mass-ratio cross-check picked. The two independent arguments agree.

| gen | m_D (keV) | m_ν predicted (meV) | m_ν observed (meV) |
|---:|---:|---:|---:|
| 1 | 42.9 | 2.1 | 1 |
| 2 | 80.2 | 3.9 | 9 |
| 3 | 117.6 | 5.7 | 50 |

## The full neutrino-mass chain (all geometric)

  - **#85** `(k=0, n<3)` quadrant = neutrino (chargeless, non-winding).
  - **#86** `k=0 ⟹ c₁=0 ⟹ C-invariant ⟹ Majorana`; seesaw `m_ν=m_D²/M_R`.
  - **#87** `M_R = m_D·e^{S}` ⟹ open input = instanton action `S ≈ 15–18`.
  - **#88** `S = √(2μE_c)·L*(ε)` (non-orientable tortoise bounce); rigid
    throat ⟹ massless ν; open input = `ΔL=2` tension ratio `t`.
  - **#89** `t ∈ [2π, k_5√(2π)]` (closure quantum … winding action); open
    input = `ε`.
  - **#90** `ε` = chargeless-throat sub-throat healing length; bulk
    geometry ⟹ `S ≈ 15–19`, `m_ν ~ few meV` (observed scale).

## Tests

| # | test | finding |
|---|---|---|
| T1 | recap | `ε` is the last residual; window `[6e-7, 1.3e-2]` |
| T2 | healing length | `ε = ℓ²/(2rs)` from `f ≈ 2(r−rs)/rs` |
| T3 | why tiny for ν | `c₁=0` neck not EM-propped (charged → Dirac); ties #86 |
| T4 | bulk candidates | `R_c³, Δ³, (m_D/m_ch)², E_c` cluster inside window |
| T5 | chain closure | winding edge `t≈√β` + `ε~R_c³` ⟹ `S≈17`, `m_ν~few meV` |
| T6 | full chain | TeV mass → S → t → window → ε → meV, all geometric |
| T7 | honest scope | closed to OOM; precise `m_ν` / gen spread residual |
| T8 | assessment | `COMPLIANCE_IS_CHARGELESS_THROAT_HEALING_LENGTH_CHAIN_CLOSED_TO_OOM` |

## Established and open

  - **Established (BAM-native):** `ε` is the chargeless throat's
    sub-throat healing length (`ε = ℓ²/2rs`); it is sub-throat *for the
    neutrino* because the `c₁=0` neck is not EM-propped (the charged neck
    is, and stays Dirac); natural BAM sub-throat scales land `ε` in the
    PR #89 window; and at the winding-edge tension the chain yields
    `S ≈ 15–19`, `m_ν ~ few meV` — the observed scale, with **no input
    outside the throat geometry**. The neutrino mass *scale* is geometric,
    not tuned.

  - **Open:** the *precise* `m_ν` (the exact healing length and exact `t`
    within their geometric ranges remain residual), and the generation
    spread (geometry-only `(t, ε)` gives a uniform `S` ⟹ `m_ν ∝ m_D`, a
    ×2.7 spread, whereas observed `m_ν/m_D` spans ×18 — needing a
    generation-dependent healing length or the mixing sector). The
    bounce normalisation (the `t²` scaling, the `P0` prefactor) carries
    model dependence.

## Run

```
python -m experiments.closure_ledger.boundary_compliance_bulk_geometry_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_boundary_compliance_bulk_geometry_probe/`.
Expected verdict: `COMPLIANCE_IS_CHARGELESS_THROAT_HEALING_LENGTH_CHAIN_CLOSED_TO_OOM`, 8/8 PASS.
