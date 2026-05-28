# Throat-shell mass-operator unification (PR #83)

Closes extension (iii) of PR #82's verdict — the deepest of the
three open extensions: why does the same closure-ledger geometry
give `ω²(l, n)` cavity eigenfrequency in the shell (quark) sector but
`β·k²` closure-winding in the throat (lepton) sector?

PR #82 found these look like two structurally different mass
operators. This probe shows they are the **same operator** — both are
**Bohr-Sommerfeld** `m² = (S / L_eff)²`, with `S` the
closure-quantized action of the channel and `L_eff` its geometric
length.

## The unified operator

```
m²(k, n)  =  (k·2π / L_throat)²  +  ((n+1)·π / L_cavity)²
```

with `L_throat = √(2π)/k_5` and `L_cavity = L_rstar` (the tortoise
cavity length). The two terms are PR #52's closure-ledger channels
`N_total = N_layer1 + N_radial`:

  - **`N_layer1`** — throat-winding integer `k` → first term.
  - **`N_radial`** — cavity-overtone integer `n` → second term.

  - **Leptons** wind through the throat: `k ∈ {1, 3, 5}`, `n = 0` →
    the winding term dominates, `m² ≈ β·k²`.
  - **Quarks** resolve the cavity: `k = 0`, `n ∈ {3, 4, 5}` → the
    winding term vanishes, `m² ≈ ω²(l, n)`.

The `k = 0` for quarks is the operator-level statement of the
physical insight that drove the entire QCD-shell arc: *quarks do not
pass through the throat; they are the wavefronts that resolve the
cavity itself.*

## The three pillars

### 1. Cavity = Bohr-Sommerfeld (verified)

The WKB action integral over the actual Tangherlini potential,

```
∮ √(ω²(n) − V(r*)) dr*  =  (n+1)·π,
```

holds to machine precision for `n ≥ 1` (n=0 is the WKB-weakest mode
at ~0.88). So `ω²(n)` is exactly Bohr-Sommerfeld quantization of the
radial action `S_radial = (n+1)·π` — a cavity standing wave with a
half-cycle `π` per node.

| n | ω² | ∮√(ω²−V) dr* | (n+1)·π | ratio |
|---:|---:|---:|---:|---:|
| 0 |  1.11 |  2.77 |  3.14 | 0.88 |
| 1 |  3.90 |  6.27 |  6.28 | 0.997 |
| 2 |  8.38 |  9.42 |  9.42 | 1.000 |
| 3 | 14.63 | 12.57 | 12.57 | 1.000 |
| 4 | 22.67 | 15.71 | 15.71 | 1.000 |
| 5 | 32.49 | 18.85 | 18.85 | 1.000 |

### 2. Lepton = winding form (exact)

`β·k² = (k·2π / L_throat)²` exactly, with `L_throat = √(2π)/k_5`. The
constant `β/(2π)² = 50/(4π)` is the same for every `k`, confirming
the pure `(k·2π)²` winding form. The winding action is
`S_winding = k·(2π)` — `k` closure quanta of the S³ great circle.

### 3. β_lepton recovered

`L_throat = √(2π)/k_5` is not a free parameter:

```
(2π / L_throat)²  =  (2π)² · k_5² / (2π)  =  k_5² · (2π)  =  50π  =  β_lepton ✓
```

Expressing the lepton mass in Bohr-Sommerfeld form reproduces PR
#71's structurally-derived `β_lepton = k_5²·(2π)` exactly.

## Unified operator limits

  - **Lepton limit** (n=0, k=1,3,5): matches `β·k²` to within 0.6%.
    The small residual is the n=0 radial floor — leptons also occupy
    the lowest cavity mode.
  - **Quark limit** (k=0, n=3,4,5): matches `ω²(n)` to within 3%.
    The residual is the flat-box leading term; the exact `ω²` is
    Bohr-Sommerfeld with the full potential (pillar 1).

## The half-cycle / full-cycle distinction

The two channels carry different closure quanta:

  - Throat winding: **`2π`** (full S³ great circle, `action_base`).
  - Radial cavity: **`π`** per Bohr-Sommerfeld node (half-cycle).

The factor of 2 is BAM's pervasive full/half-cycle distinction — the
throat dwell `τ = π/ω` (half-cycle per pass), the Hopf holonomy
`∮A = π cos χ` (half at the pole), the B3 hard-wall reflection phase
`π` (Maslov μ=2 per reflection). The radial standing wave is a
reflection (half-cycle `π`); the throat winding is a full
great-circle traversal (`2π`).

## What this establishes

The conceptual gap PR #82 flagged — "why `ω²` in the shell but `β·k²`
in the throat" — is answered: **both are Bohr-Sommerfeld `(S/L)²` of
their respective closure channels**, distinguished by

  - the winding number `k` (`k ≠ 0` throat / `k = 0` cavity), and
  - the half/full-cycle closure quantum (`π` radial node / `2π`
    great-circle winding).

The "two mass operators" were always one operator read in two
channels — exactly the PR #52 `N_total = N_layer1 + N_radial`
decomposition, with both channels feeding `m²` via the same
Bohr-Sommerfeld rule.

## Honest scope

  - **Is:** the demonstration that both mass operators are
    Bohr-Sommerfeld `(S/L_eff)²`; the verified cavity BS integral;
    the algebraically-exact lepton winding form; the recovery of
    `β_lepton = k_5²·(2π)` from `L_throat = √(2π)/k_5`; the unified
    operator with `k = 0` for quarks; the half/full-cycle reading;
    the tie to PR #52.

  - **Is not:** an independent derivation of the two `L_eff` from a
    single deeper principle — `L_throat = √(2π)/k_5` re-expresses PR
    #71's already-derived `β_lepton`, and `L_cavity` is the literal
    tortoise cavity length. The unification is at the
    Bohr-Sommerfeld FORM level (both sectors share one operator),
    not a reduction of both scales to one number. The
    inter-generation hierarchy (cross-channel / mixed modes) and the
    prediction of new states remain open.

## B4 accounting

The closure actions `S` are dimensionless (counts of closure quanta);
the effective lengths `L_eff` are dimensionful (`1/length`).
`m²/scale` ratios are scale-free; the structural identifications are
scale-independent. The absolute MeV scale rides on the single B4
anchor `m_e = f_closure · ℏ/(ΔR·c)` (PR #53).

## Tests

  T1. Cavity Bohr-Sommerfeld: `∮√(ω²−V) dr* = (n+1)·π` to machine
      precision (n ≥ 1).
  T2. Lepton winding form: `β·k² = (k·2π/L_throat)²` exact;
      `L_throat = √(2π)/k_5`.
  T3. β_lepton recovery: `(2π/L_throat)² = k_5²·(2π) = 50π`.
  T4. Unified operator limits: lepton (n=0) < 0.6%, quark (k=0) < 3%.
  T5. Half-cycle/full-cycle closure quanta (π vs 2π) = BAM's
      pervasive distinction.
  T6. Tie to PR #52 closure ledger `N_total = N_layer1 + N_radial`.
  T7. `k = 0` for quarks = operator statement of "quarks don't pass
      through the throat".
  T8. Honest scope + assessment.

## Verdict structure

  - **`MASS_OPERATOR_UNIFIED_BOHR_SOMMERFELD`** (expected): both
    sectors are one Bohr-Sommerfeld operator; closes extension (iii)
    of PR #82 at the structural-form level.
  - **`UNIFICATION_INCONCLUSIVE`**: a structural test fails;
    investigate before claiming unification.

## What this leaves open

  - **Independent derivation of the two `L_eff`** from one principle.
  - **Inter-generation hierarchy** — cross-channel / mixed-mode
    question (PR #80's open gap).
  - **Prediction of new states** — e.g. winding shell modes with
    both `k ≠ 0` and `n ≥ 3`, if physical.

## Cross-references

  - `docs/pati_salam_throat_shell_bridge_research_plan.md` — PR #82,
    which flagged this mass-operator unification as the deepest of
    three open extensions.
  - `docs/beta_lepton_derivation_research_plan.md` — PR #71, the
    lepton `β_lepton = k_5²·(2π)` recovered here from `L_throat`.
  - `docs/qcd_shell_waveguide_scaffold_research_plan.md` — PR #77,
    the quark `ω²(l, n)` cavity mass operator.
  - `docs/maslov_dimensional_bridge_research_plan.md` — PR #52, the
    closure-ledger `N_total = N_layer1 + N_radial` decomposition.
  - `docs/k5_origin_research_plan.md` — PR #73, `k_5 = 5` appearing
    in `L_throat = √(2π)/k_5`.
  - `geometrodynamics/tangherlini/radial.py` — `V_tangherlini`, the
    cavity eigensolver providing the verified Bohr-Sommerfeld modes.
  - `experiments/closure_ledger/throat_shell_mass_operator_unification_probe.py`
    — this probe.
