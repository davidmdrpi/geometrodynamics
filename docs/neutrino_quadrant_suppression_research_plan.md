# Neutrino-quadrant suppression: Majorana seesaw from `c₁ = 0` (PR #86)

PR #85 mapped the unified `(k, n)` operator into four quadrants and
found the `(k = 0, n < 3)` quadrant gives the lightest states —
candidate neutrinos — but with `m_ν / m_charged ≈ 0.07`, far above the
observed `< 10⁻⁶`. So the bare quadrant is ~10⁵–10⁶ too heavy for
neutrinos. This probe identifies the BAM-native suppression mechanism.

## The structural key: `k = 0 ⟹ c₁ = 0 ⟹ Majorana`

The throat winding number `k` IS the Hopf charge: a winding mode
carries `c₁ = ±1` (the charged-lepton sector, PR #71), a non-winding
mode carries `c₁ = 0`. Under charge conjugation `C` (the inner/outer
swap `c₁ → −c₁`, PR #63):

| sector | k | c₁ | under C | nature | seesaw? |
|---|---:|---:|---|---|:---:|
| neutrino | 0 | 0 | 0 → 0 (invariant) | **Majorana** | ✓ |
| charged lepton | 2g−1 | ±1 | +1 → −1 (e⁻→e⁺) | Dirac | ✗ |

The neutrino is its own C-conjugate, hence **necessarily Majorana** —
not by assumption, but because the chargeless (`k = 0`) sector is the
C-invariant sector.

## The suppression: Majorana seesaw

A Majorana mass term violates lepton number by `ΔL = 2` and admits the
seesaw

```
m_ν  =  m_D² / M_R,
```

where `m_D` is the bare neutrino-quadrant Dirac mass (the cavity floor
`√(ω²(0, n))`) and `M_R` is the heavy Majorana mass — the
lepton-number-violating scale. In BAM, `M_R` is the throat ↔
antithroat coupling scale (a Majorana mass flips throat to antithroat,
`ΔL = 2`).

Because `M_R ≫ m_D`, the seesaw produces `m_ν ≪ m_D` automatically —
the anomalous lightness of neutrinos is **generic** to the mechanism.
The specific value needs `M_R`, but the smallness does not.

## Why only neutrinos are suppressed

The seesaw is available ONLY to the Majorana (C-invariant, `c₁ = 0`)
sector. Charged leptons (`k ≠ 0`, `c₁ = ±1`) are Dirac — under `C`,
`e⁻ → e⁺` (distinct) — so they cannot have a `ΔL = 2` Majorana mass,
get no seesaw, and keep their full winding mass `β·k²`. **This is the
BAM-native explanation of why the neutrino is anomalously light and the
charged lepton is not**: the two differ precisely by their Hopf charge
(winding) `c₁`.

## Numbers

With the electron anchored at `(k=1, n=0)`, the bare neutrino Dirac
masses are the cavity floors:

| gen | n | m_D (keV) | m_ν obs (eV) | required M_R (GeV) |
|---:|---:|---:|---:|---:|
| 1 | 0 | 42.9 | 0.001 | 1836 |
| 2 | 1 | 80.2 | 0.009 | 715 |
| 3 | 2 | 117.6 | 0.05 | 277 |

The required seesaw scale `M_R ≈ 0.3–1.8 TeV` — a single new heavy
scale, roughly generation-uniform.

## The scale `M_R` is open

No BAM scale in the current catalog lands at ~TeV:

| candidate | value | matches ~TeV? |
|---|---|:---:|
| pair-production `2 m_e c²` | ~1 MeV | ✗ |
| leptoquark gen-1 (raw operator) | sub-MeV | ✗ |
| bulk/cosmological `ℏc/R_cosmo` | ~10⁻³³ eV | ✗ |

So `M_R` — the lepton-number-violating (throat↔antithroat / B−L
breaking) scale — is a **new heavy input**, not yet BAM-derivable. The
TeV value is read off from the observed neutrino mass, not predicted.

## What this establishes (and does not)

  - **Establishes (BAM-native):** the neutrino is Majorana because the
    chargeless `k = 0` quadrant is the C-invariant sector
    (`c₁ → −c₁ = 0`); the Majorana seesaw is therefore the natural
    suppression mechanism; and the seesaw is available ONLY to the
    `c₁ = 0` sector, explaining why neutrinos (not charged leptons) are
    anomalously light.

  - **Does not establish:** the value of `M_R` (a new heavy input,
    ~TeV from observed `m_ν`, not predicted; no current BAM scale
    matches); absolute neutrino masses (need `M_R` + the L_eff
    unification still open from PR #83 + the B4 anchor); the neutrino
    mass ordering / PMNS mixing.

## B4 accounting

`m_D` is the cavity floor scaled by the electron anchor; `M_R` is open;
absolute `m_ν` needs `M_R` + the L_eff unification (PR #83) + the
single B4 anchor (PR #53). The Majorana/Dirac structural distinction
(`c₁ = 0` vs `c₁ = ±1`) is scale-free.

## Tests

  T1. Suppression gap (recap PR #85): bare ν/charged ≈ 0.07;
      observed < 10⁻⁶; ~10⁵–10⁶ suppression needed.
  T2. `k = 0 ⟹ c₁ = 0 ⟹ C-invariant ⟹ Majorana` (BAM-native).
  T3. Charged lepton `k ≠ 0 ⟹ c₁ = ±1 ⟹` Dirac, no seesaw (why only
      ν suppressed).
  T4. Bare Dirac masses = cavity floors (~43, 80, 118 keV).
  T5. Required seesaw scale `M_R = m_D²/m_ν ≈ 0.3–1.8 TeV`.
  T6. BAM heavy-scale candidates; none matches ~TeV → `M_R` open.
  T7. Honest scope: mechanism BAM-native, scale `M_R` open.
  T8. Assessment.

## Verdict structure

  - **`NEUTRINO_SUPPRESSION_IS_MAJORANA_SEESAW_SCALE_OPEN`**
    (expected): the suppression mechanism is the Majorana seesaw,
    forced by `c₁ = 0` C-invariance and available only to the
    chargeless sector; the seesaw scale `M_R` is a new heavy input,
    open.
  - **`NEUTRINO_SUPPRESSION_INCONCLUSIVE`**: a structural test fails.

## What this leaves open

  - **The seesaw scale `M_R`** — the lepton-number-violating
    (throat↔antithroat / B−L) scale; ~TeV from observed `m_ν`, not
    predicted; no current BAM scale matches.
  - **Absolute neutrino masses** — need `M_R` + L_eff unification
    (PR #83) + B4 anchor.
  - **Neutrino mass ordering / PMNS mixing** — beyond this scope.

## Cross-references

  - `docs/winding_shell_quadrant_research_plan.md` — PR #85, the
    four-quadrant map that surfaced the candidate neutrino quadrant
    and its mass-scale caveat.
  - `docs/charge_conjugation_swap_research_plan.md` — PR #63, the
    inner/outer swap `C: c₁ → −c₁` that makes the `c₁ = 0` sector
    C-invariant (Majorana).
  - `docs/beta_lepton_derivation_research_plan.md` — PR #71, the
    charged-lepton winding sector (`c₁ = ±1`, Dirac).
  - `docs/throat_shell_mass_operator_unification_research_plan.md` —
    PR #83, the unified operator providing the bare cavity-floor
    Dirac masses.
  - `docs/pati_salam_throat_shell_bridge_research_plan.md` — PR #82,
    which flagged BAM-native neutrinos (the Majorana option) as one of
    three open extensions.
  - `experiments/closure_ledger/neutrino_quadrant_suppression_probe.py`
    — this probe.
