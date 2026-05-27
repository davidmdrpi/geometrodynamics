# Boundary stress derivation of `χ_n` and singlet constraint (PR #79)

Derive or constrain the generation-dependent Z₂ partition splitter
`χ_n` from boundary stress across the inner/outer shell mouths, and
add the singlet (color-neutral) projection constraint on physical
states. Third PR of the four-PR QCD-shell arc (#77 scaffold → #78
mass-ordering audit → #79 boundary stress → #80 color algebra).

## The boundary stress derivation

The Tangherlini radial equation `H ψ_n = ω² ψ_n` on `[R_MID, R_OUTER]`
with Dirichlet at both ends defines a scalar field whose boundary
stress at each mouth is `T_rr|_boundary ∝ (∂ψ_n/∂r*)²`. For
normalized `ψ_n` (`∫|ψ_n|² dr* = 1`):

```
T_inner(n)  =  (∂ψ_n/∂r*)²  at  r = R_MID + ε       (throat-side mouth)
T_outer(n)  =  (∂ψ_n/∂r*)²  at  r = R_OUTER − ε     (cavity-edge mouth)
```

The inner/outer swap `r ↦ 2·R_MID − r` is PR #63's charge-conjugation
involution `C`. Decompose under it:

```
T_even(n)  =  (T_inner + T_outer) / 2   (Z₂-symmetric, mass² shift)
T_odd(n)   =  (T_inner − T_outer) / 2   (Z₂-antisymmetric, χ_n)
```

The Z₂-antisymmetric piece is the structurally derived `χ_n` — no
free parameter once cavity geometry is fixed.

## Numerical findings (l = 1, n = 0..5)

| n | ω | ω² | T_inner | T_outer | T_odd/T_even |
|---:|---:|---:|---:|---:|---:|
| 0 | 1.055 |  1.11 | 6.7×10⁻¹ | 3.6×10⁻¹ | 0.304 |
| 1 | 1.974 |  3.90 | 2.26     | 1.86     | 0.097 |
| 2 | 2.894 |  8.38 | 4.94     | 4.57     | 0.039 |
| 3 | 3.825 | 14.63 | 8.70     | 8.34     | 0.021 |
| 4 | 4.761 | 22.67 | 13.53    | 13.18    | 0.013 |
| 5 | 5.700 | 32.49 | 19.44    | 19.10    | 0.009 |

Three structural findings:

  1. **`T_inner > T_outer` for every mode** — uniform-positive sign of
     `T_odd`. The 5D Tangherlini centrifugal + throat-curvature
     potential shifts the radial profile toward the throat side.
  2. **Asymmetry decreases monotonically with `n`** — from 0.30
     (electron-focused) to 0.009 (shell-saturated). Shell modes feel
     the mouth asymmetry only weakly because the standing wave fills
     the cavity nearly uniformly.
  3. **No sign flip between `n=3` and `n=4`** — PR #78's existence
     proof (`χ_3 < 0, χ_4, χ_5 > 0`) is **structurally incompatible**
     with the boundary-stress derivation.

## Magnitude audit

For the within-generation mass-ordering ratio
`m²_light / m²_heavy = (ω² − χ)/(ω² + χ) = r`, the required
`χ/ω² = (1 − r)/(1 + r)`. With observed quark ratios:

| splitting | required χ/ω² | derived χ_n/ω² | deficit |
|---|---:|---:|---:|
| u/d at n=3 | 0.647 | 0.021 | 31× |
| s/c at n=4 | 0.986 | 0.013 | 76× |
| b/t at n=5 | 0.954 | 0.009 | 106× |

Boundary stress alone gives `χ_n/ω²` 30–100× too small to drive the
observed within-generation splittings.

## The honest reading

**Boundary stress IS the right structural slot for `χ_n`.** It
provides:

  - The correct Z₂-antisymmetric origin (PR #63's `C` involution).
  - The correct qualitative scaling (largest for throat-focused
    leptons, smallest for shell-saturated quarks).
  - No free parameter — once cavity geometry is fixed, `χ_n` is
    determined by the Tangherlini eigensolver.

**But boundary stress is INSUFFICIENT for the quark splittings.** The
within-generation inversion (and the inter-generation hierarchy)
must come from **PR #80's color algebra** populating `H_couple`.

## The v3 species ↔ partition map: flagged

The v3 convention `(k=1, +) = u` (with u lighter than d) is
inconsistent with the natural boundary-stress reading: `χ_n > 0`
uniformly ⟹ + heavier than − at every generation. Two reconciliation
paths:

  - **(a) Revise the v3 map.** Under "+ = heavier" uniformly:
    `(n=3, +) = d, (n=3, −) = u`; `(n=4, +) = c, (n=4, −) = s`;
    `(n=5, +) = t, (n=5, −) = b`. All generations have the same
    partition convention; uniform-sign `χ_n` is consistent with the
    structural reading.
  - **(b) Keep the v3 map; attribute the n=3 inversion to a
    different mechanism.** The within-generation inversion at n=3
    (u<d) would then arise from PR #80's color sector or another
    mechanism, and `χ_n` from boundary stress contributes only a
    small additional splitting.

PR #80 settles this.

## Singlet constraint

The 6-state shell waveguide basis is at the **flavor** level
(`u, d, s, c, b, t`), not at the color level. Color is implicit;
physical states are color singlets. The singlet projector `P_S` is
the identity on the flavor basis — a placeholder until PR #80
identifies the color algebra. `P_S` commutes with the diagonal
flavor Hamiltonian; the spectrum is unchanged under projection.

The non-trivial singlet content arrives in PR #80, where the color
algebra acting on `(l, n, p)` is identified and the projection onto
color-singlet operators becomes meaningful.

## What PR #79 fixes

  1. **Structural origin of `χ_n`** = Z₂-antisymmetric piece of the
     cavity-mouth boundary stress.
  2. **No free parameter** once cavity geometry `[R_MID, R_OUTER]`
     is fixed.
  3. **Uniform-sign / shell-suppressed** structural pattern.
  4. **PR #78 sign-flipping ansatz overruled** — not derivable from
     boundary stress.
  5. **Singlet projector placeholder** — identity on flavor basis;
     awaits PR #80's color algebra.

## What PR #80 must contribute

  1. **Color algebra** acting on `(l, n, p)` — SU(3), SU(2) × Z₂,
     Pati-Salam SU(4), or another candidate.
  2. **`H_couple` inter-mode mixing** populated by the color-algebra
     transformation rules.
  3. **Within-generation splittings** beyond the small boundary-stress
     `χ_n`.
  4. **Inter-generation hierarchy** spanning ~9 orders of magnitude
     in mass².
  5. **v3 species ↔ partition map settlement**.

After PR #80 populates the operator, the `n_part` audit can be re-run
honestly. Until then, `n_part = 233` remains a residual
phenomenological compensator with sharpened scope.

## Honest scope

  - **Is:** structural derivation of `χ_n` with no free parameter;
    identification of the uniform-sign / shell-suppressed pattern;
    overrule of PR #78's sign-flipping ansatz; quantitative magnitude
    audit (50–100× too small); singlet projector placeholder.
  - **Is not:** a derivation of quark masses; a resolution of the
    within-generation inversion or inter-generation hierarchy; an
    identification of the color algebra; a re-audit of `n_part` (all
    pending PR #80).

## B4 accounting

`T_inner, T_outer` are dimensionful (`1/length²` for the normalization
chosen). The ratio `χ_n/ω²` is **scale-free**; the structural pattern
(uniform sign, shell suppression) is scale-independent. The absolute
MeV scale rides on the single B4 anchor `m_e = f_closure · ℏ/(ΔR·c)`
(PR #53).

## Tests

  T1. Compute `T_inner(n), T_outer(n)` for l=1, n=0..5; verify
      positivity (Dirichlet boundary conditions ⟹ nonzero slope).
  T2. Z₂ decomposition under the inner/outer swap (PR #63's `C`):
      `T_even, T_odd`. Verify `T_odd > 0` uniformly; asymmetry
      decreases with `n`.
  T3. `χ_n = T_odd(n)` derivation; report values and `χ_n/ω²`
      ratios.
  T4. Sign-pattern audit vs PR #78's existence proof: no sign flip
      → PR #78 ansatz overruled.
  T5. Magnitude audit: required vs derived `χ/ω²`; 30–100× deficit
      across all three blocks.
  T6. Singlet projector placeholder = identity on the 6-state flavor
      basis; commutes with flavor Hamiltonian.
  T7. Structural interpretation: PR #79 fixes the `χ_n` slot; PR #80
      owes the color sector.
  T8. Honest scope + verdict.

## Verdict structure

  - **`CHI_N_DERIVED_BOUNDARY_STRESS_INSUFFICIENT`** (expected): T1–T7
    pass; structural slot derived without free parameter; magnitude
    too small and sign uniform; PR #78 sign-flipping ansatz overruled;
    PR #80's color sector flagged as the necessary contributor.
  - **`CHI_N_NOT_DERIVED`**: a structural test fails; investigate
    before proceeding to PR #80.

## What this leaves open

  - **PR #80** — color algebra acting on `(l, n, p)`; `H_couple`
    inter-mode mixing.
  - **v3 species ↔ partition map** — flagged for revision; PR #80
    settles.
  - **`n_part` audit re-run** — after PR #80 populates the operator.

## Cross-references

  - `docs/qcd_shell_waveguide_scaffold_research_plan.md` — PR #77,
    the basis + operator scaffold (`H_kin + H_Z2 + H_couple`).
  - `docs/shell_mass_ordering_audit_research_plan.md` — PR #78, the
    mass-ordering audit that posed the sign-flipping `χ_n` question.
  - `docs/charge_conjugation_swap_research_plan.md` — PR #63, the
    inner/outer swap `r ↦ 2·R_MID − r` (`C` involution).
  - `docs/throat_dirac_spinor_research_plan.md` — PR #66, the SUSY
    factorization of the radial Dirac operator with the two-mouth
    structure.
  - `geometrodynamics/tangherlini/radial.py` — `V_tangherlini`,
    `r_to_rstar`, `rstar_to_r` (the cavity eigensolver).
  - `experiments/closure_ledger/qcd_shell_waveguide_scaffold_probe.py`
    — PR #77 basis + operator that this probe extends.
  - `experiments/closure_ledger/boundary_stress_chi_n_probe.py` —
    this probe.
