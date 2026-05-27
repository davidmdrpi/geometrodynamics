# QCD shell waveguide: basis + operator scaffold

> *Quarks do not pass through the throat; they are the wavefronts that
> resolve the cavity itself.*

Foundation for the four-PR quantitative QCD-shell arc the user laid
out:

  - **PR #77 (this PR)** — QCD shell waveguide basis/operator scaffold
  - **PR #78** — shell Hamiltonian mass-ordering / `n_part` audit
  - **PR #79** — boundary stress tensor and singlet constraint
  - **PR #80** — color algebra candidate: SU(3), SU(2) × Z₂, or other

## The physical insight

The user's reframe of the quark sector: quarks are not throat
traversals (as the lepton ladder is — e, μ, τ are odd-`k` modes that
pass through the non-orientable throat). They are the **shell-saturated
standing waves on the cavity** `[R_MID, R_OUTER]` — wavefronts whose
⟨r⟩ has plateaued at the shell center and whose participation ratio
has locked to the uniform-standing-wave value `2/3`. The mode
resolves the cavity rather than focuses on the throat.

This is the quantitative development of PR #68 (focused pulse →
extended wavefront) + PR #69 (shell ↔ QCD structural match) that
PR #76's verdict identified as the **right derivation route** for
the quark sector. The v3 Hamiltonian was the wrong machinery; the
shell waveguide is the right machinery.

## The shell waveguide basis

A **shell waveguide basis state** is a tuple `(l, n, p)`:

  - `l` — S³ angular momentum (Casimir `l(l+2)`, the Hopf-bundle's
    angular base, BAM's foundational primitive per PR #73).
  - `n` — radial overtone index, shell-saturated branch. For `l = 1`,
    `n ≥ 3` is the saturated regime (per PR #68's metrics).
  - `p ∈ {+, −}` — Z₂ partition from B2 (the non-orientable throat
    `T = iσ_y`, `T² = −I`).

The lowest 6 shell-saturated states constitute the **3 × 2 = 6 flavor**
structural count documented by PR #69.

Two natural enumerations of the 6 states:

  - **`n_varied`**: fix `l = 1`, vary `n ∈ {3, 4, 5}`, both `p`.
  - **`l_varied`**: fix `n = 3` (shell-saturated), vary `l ∈ {1, 2, 3}`,
    both `p`.

PR #78 will choose between them based on which reproduces the quark
mass ordering. PR #77 scaffolds both.

## The operator scaffold

The shell Hamiltonian acts on the 6-state basis as a 6×6 Hermitian
matrix:

```
H  =  H_kin  +  H_Z2  +  H_couple
```

  - **`H_kin`** — diagonal cavity-mass operator. Per-mode entry:
    `ω²(l, n)` from the Tangherlini radial eigensolver. This is the
    "energy resolves the cavity" content; distinct from the lepton
    Hamiltonian's diagonal `β·k²·(2π)` (winding cost on the throat
    traversal, PR #71).
  - **`H_Z2`** — block-`σ_z` partition splitter on the `p = ±` index
    per `(l, n)` block. Strength `χ` is a slot for PR #79–#80
    structural input (cavity-wall stress tensor / color algebra). PR
    #77 leaves `χ = 0`.
  - **`H_couple`** — inter-mode coupling matrix between different
    `(l, n)` blocks. Zero placeholder; PR #78–#80 populate.

With `χ = 0` and `H_couple = 0`, the eigenvalues are `{ω²(l, n)}` each
doubled by the Z₂ partition.

## Distinctness from the lepton / v3 machinery

| sector | lives in | basis | diagonal |
|---|---|---|---|
| **lepton** (PR #71) | throat (odd-`k` traversal) | `{(k=1,±), (k=3,±), (k=5,±)}` | `β·k²·(2π) = 50π·k²` |
| **shell waveguide** (this PR) | cavity wavefronts (shell-saturated) | `{(l, n, p)}_{n≥3}` | `ω²(l, n)` cavity eigenfrequency squared |

The lepton ladder rides the throat with the closure-quantum winding
cost; the shell waveguide rides the cavity wavefronts with the
cavity-mode eigenfrequency squared. Per PR #76's diagnosis, fitting
quarks on the lepton basis absorbed the unmodeled QCD physics into
the phenomenological `n_part = 233`. The shell waveguide basis is
structurally distinct and is the right machinery to host that
physics natively.

## Honest scope

  - **Is:** the 6-state shell waveguide basis (both enumerations), the
    6×6 operator scaffold `H_kin + H_Z2 + H_couple` with `χ` and
    `couple_matrix` as named slots, distinctness verification (lepton
    ladder = throat-to-shell transition, shell modes on the plateau),
    and named hooks for PR #78–#80.

  - **Is not:** a derivation of quark masses (PR #78), an audit of
    `n_part` on the shell basis (PR #78), a definition of the boundary
    stress tensor or singlet constraint (PR #79), or an identification
    of the color algebra (PR #80). The 6×6 operator is constructed but
    deliberately unpopulated.

## B4 accounting

Cavity eigenfrequencies `ω(l, n)` are dimensionful (`1/length`); mass
**ratios** between shell modes are scale-free. The scaffold is
structurally scale-independent — `χ`, `couple_matrix` are
dimensionless coupling constants. The absolute scale rides on the
single B4 anchor `m_e = f_closure · ℏ/(ΔR·c)` (PR #53).

## Tests

  T1. **Shell waveguide basis constructed.** 6 states (l=1, n=3,4,5,
      ±); each saturated (⟨r⟩ on plateau ≥ 0.040, PR within 0.02 of
      2/3).
  T2. **Lepton-to-shell transition.** Lepton ladder (n=0,1,2) is the
      throat-to-shell transition (⟨r⟩ rises 0.021 → 0.039 → 0.044,
      below plateau); shell modes (n≥3) saturated on plateau; electron
      clearly below plateau.
  T3. **Operator scaffold.** `H_kin` diagonal from ω²(l,n); `H_Z2`
      zero placeholder (`χ = 0`); `H_couple` zero placeholder;
      `H_total` Hermitian; eigenvalues equal `{ω²}` each doubled.
  T4. **Basis flexibility.** Both `n_varied` and `l_varied`
      enumerations supported; eigenvalues differ — PR #78 will choose.
  T5. **Six-flavor structural count.** 3 generations × 2 partitions =
      6 flavors (PR #69).
  T6. **Hooks for PR #78–#80.** Named slots: PR #78 populates
      `H_couple` and chooses the enumeration; PR #79 sets `χ` from the
      boundary stress tensor + adds singlet projection; PR #80
      identifies the color algebra acting on `(l, n, p)`.
  T7. **Honest scope / B4.** Scaffold only; no masses, no `n_part`
      audit. Cavity eigenfrequencies dimensionful; scaffold scale-free
      in ratios.
  T8. **Assessment.**

## Verdict structure

  - **`SHELL_WAVEGUIDE_SCAFFOLD_CONSTRUCTED`** (expected): the 6-state
    shell basis is built, the 6×6 operator scaffold has the named
    slots, the basis is distinct from lepton/v3 throat modes, both
    enumerations are supported, and the hooks for PR #78–#80 are in
    place.
  - **`SCAFFOLD_INCOMPLETE`**: a structural test fails; investigate
    before proceeding to PR #78.

## What this leaves open

  - **PR #78** — populate `H_couple` with inter-mode mixing, choose
    between `n_varied` and `l_varied`, and test whether the shell
    basis reproduces the quark mass ordering (m_u < m_d, m_c > m_s,
    etc.) AND derives `n_part` structurally rather than as a
    phenomenological compensator.
  - **PR #79** — define the boundary stress tensor on the cavity wall
    and set `χ` from it; add a singlet (colorless) projection as a
    physical-state constraint.
  - **PR #80** — identify the color algebra acting on `(l, n, p)`:
    SU(3) (color triplet, the canonical QCD reading), SU(2) × Z₂
    (partition-flavored, smaller), Pati-Salam SU(4) (the `n` overtone
    as a 4th leptocolor), or another candidate.

## Cross-references

  - `docs/throat_to_shell_transition_research_plan.md` — PR #68, the
    focused-pulse → extended-wavefront metric framework this scaffold
    builds on.
  - `docs/shell_to_qcd_match_research_plan.md` — PR #69, the
    3 × 2 = 6 flavor structural match.
  - `docs/quark_npart_origin_research_plan.md` — PR #76, the diagnosis
    that the v3 Hamiltonian is the wrong machinery; this PR is the
    right machinery.
  - `docs/beta_lepton_derivation_research_plan.md` — PR #71, the
    lepton-sector analogue with the clean closure-quantum integer
    `4·k_5² = 100`.
  - `docs/k5_origin_research_plan.md` — PR #73, the `k_5 = D_bulk =
    dim(S³) + 2 = 5` primitive and the `l(l+2)` S³ Casimir on the
    Hopf-bundle's angular base.
  - `geometrodynamics/tangherlini/radial.py` — `solve_radial_modes`
    and `V_tangherlini`, the cavity eigensolver.
  - `experiments/closure_ledger/throat_to_shell_transition_probe.py`
    — `radial_ladder` and the saturation metrics this scaffold reuses.
  - `experiments/closure_ledger/qcd_shell_waveguide_scaffold_probe.py`
    — this probe.
