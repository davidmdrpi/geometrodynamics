# Majorana instanton-bounce action probe (PR #88)

Follows PR #87, which reframed the neutrino's Majorana seesaw scale as
an instanton/bounce action `M_R = m_D·e^{S}`, `S = ln(m_D/m_ν) ≈ 15–18`,
but left `S` itself OPEN (verdict
`SEESAW_SUPPRESSION_IS_THROAT_ANTITHROAT_TUNNELING_S_OPEN`). This probe
constructs a **reduced Euclidean bounce action** for the `ΔL=2`
throat ↔ antithroat transition from first BAM ingredients — the PR #58
nucleation potential, the PR #86 chargeless (`k=0, c₁=0`) Majorana
sector, and the non-orientable tortoise geometry — and asks whether the
intrinsic throat tension and boundary compliance naturally yield
`S ≈ 15–18`.

## The non-orientable tortoise path

The throat sits at `r = rs = R_MID`. The 5D tortoise coordinate

```
r*(r) = r + (rs/2)·ln|(r − rs)/(r + rs)|
```

diverges **logarithmically** as `r → rs` (slope `rs/2`): the throat is
infinitely deep in tortoise distance. The `ΔL=2` Majorana transition is
the throat ↔ antithroat flip — the inner/outer swap `C` (`c₁ → −c₁`,
PR #63), realised in `solve_radial_modes` as the **odd extension across
the throat** (`u_full = [−u[::-1], u]`, the orientation-reversing /
non-orientable identification). The bounce is the Euclidean path of the
throat configuration along this non-orientable tortoise path.

## The reduced bounce action

Flipping the throat orientation costs the PR #58 nucleation/pinch energy
`E_c = (16π/3)σ³/ρ²` and carries the throat's breathing inertia
`μ = 4πσ rs²` (the brane mass at the throat scale). The reduced bounce
(thin-wall momentum × path length) is

```
S = ∫ √(2 μ E_c) dr* = √(2 μ E_c) · L*(ε),
   L*(ε) = r*(R_OUTER) − r*(rs + ε) = −(rs/2) ln ε + const,
```

with `ε` a **boundary-compliance cutoff** (how close to the throat the
path approaches). `S ∝ ln(1/ε)` — the action is a tortoise logarithm.

## Two structural results (BAM-native)

  1. **Rigid throat ⟹ massless neutrino.** A perfectly rigid
     (incompressible) throat has `ε → 0`, so `L* → ∞`, `S → ∞`, and
     `m_ν = m_D·e^{−S} → 0`. The boundary compliance `ε > 0` is exactly
     what gives the neutrino mass; the smallness of `m_ν` is the
     near-rigidity of the throat.

  2. **`S` is a tortoise logarithm.** `S ∝ ln(1/ε)` is naturally `O(10)`
     and varies only slowly — matching the `O(15)`, generation-stable
     `S` PR #87 required. The keV→TeV gap is a logarithm, not a
     hierarchy.

## The magnitude test (honest)

With the throat tension fixed by the **electron** throat (PR #58/#87:
`σ = 1/(12π)`, `ρ = 3/(4π)`, geometric units `mc²=1`, `R*=1`), the
under-barrier momentum is `√(2 μ E_c) ≈ 0.060`. Then:

| compliance ε | tortoise L* | bounce S |
|---:|---:|---:|
| 1e-2 | 1.82 | 0.110 |
| 1e-4 | 4.13 | 0.250 |
| 1e-8 | 8.74 | 0.528 |
| 1e-16 | 17.9 | 1.085 |

Even an absurdly rigid `ε = 1e-16` gives `S ≈ 1.1` — **~40× short** of
the required `~16`. Pure compliance cannot close the gap (it would need
`ε ~ 10⁻²³⁰`). Matching `S ∈ [15, 18]` at a *sane* near-rigid compliance
(`ε ~ 10⁻³–10⁻⁶`) requires the `ΔL=2` (lepton-number-violating / B−L)
throat tension to be a factor `t ≈ 6–12` stiffer than the EM-throat
tension (barrier `~10²–10³×` larger, since `E_c ∝ σ³`):

| compliance ε | required tension ratio t | barrier boost t³ |
|---:|---:|---:|
| 1e-2 | 12.0 | 1750× |
| 1e-3 | 9.4 | 840× |
| 1e-6 | 6.4 | 260× |

## Tests

| # | test | finding |
|---|---|---|
| T1 | required `S` | `S = ln(m_D/m_ν) ≈ 14.7–17.6`, gen-stable (PR #87 target) |
| T2 | tortoise path | odd extension across throat; `r*` log-diverges, slope `rs/2` |
| T3 | rigid limit | `ε→0 ⟹ S→∞ ⟹ m_ν=0`; compliance = mass-generating parameter |
| T4 | EM-calibrated bounce | `S ≲ 1` even near-rigid → ~40× short |
| T5 | tension ratio | matching `S∈[15,18]` needs `t ≈ 6–12` (B−L/EM) |
| T6 | generation | uniform `S ⟹ m_ν∝m_D`; `×18` spread needs gen-compliance/mixing |
| T7 | honest scope | structure established; open input → `O(10)` tension ratio |
| T8 | assessment | `MAJORANA_BOUNCE_IS_TORTOISE_LOG_S_OPEN_AS_TENSION_RATIO` |

## Established and open

  - **Established (BAM-native):** the bounce coordinate is the
    non-orientable tortoise path (odd extension, `c₁ → −c₁`); the rigid
    throat gives an exactly massless Majorana neutrino, so the compliance
    is the mass-generating parameter; the bounce action is a tortoise
    logarithm `S ∝ ln(1/ε)`, hence naturally `O(10)` and coarsely
    generation-stable — the form PR #87's `S` required.

  - **Open:** `S ≈ 15–18` does **not** fall out of the EM-throat tension
    (under by `~40×`); matching needs a `ΔL=2` tension `~6–12×` stiffer
    (a dimensionless ratio, not yet derived), the compliance `ε` (a
    sub-throat length, not yet derived from the bulk), and the detailed
    `m_ν` spectrum (the `×18` `m_ν/m_D` spread needs a generation-
    dependent compliance or the mixing sector).

**Progressive localisation of the open input:** a mysterious `~TeV` mass
(PR #86) → an `O(15)` instanton action `S` (PR #87) → an `O(10)`
dimensionless B−L/EM tension ratio (this probe). Each step ties the
unknown more tightly to BAM geometry.

## Run

```
python -m experiments.closure_ledger.majorana_bounce_action_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_majorana_bounce_action_probe/`.
Expected verdict: `MAJORANA_BOUNCE_IS_TORTOISE_LOG_S_OPEN_AS_TENSION_RATIO`, 8/8 PASS.
