# Pair-annihilation crossing — research plan

Follow-on to PR #36 (Breit-Wheeler cross-process validation: the
closed-form Compton vertex factor F²(x, c) crossed via Mandelstam
substitution reproduces γγ → e⁺e⁻ exactly). PR #36 closed one edge of
the Compton ↔ BW ↔ annihilation crossing triangle. This probe closes
the remaining two edges:

  - Compton → annihilation directly via the (s_C → u_ann, t_C → s_ann,
    u_C → t_ann) crossing.
  - BW ↔ annihilation via T-invariance (the two processes share the
    same |M̄|²(s, t, u) function but differ in phase space / flux).

and verifies that the **loop closes** — that the BAM kernel produces a
consistent amplitude regardless of which crossing path is taken.

## The triangle

```
          Compton  (γe⁻ → γe⁻)
            / \
           /   \
          /     \
         /       \
        v         v
       BW ◀────▶ annihilation
    (γγ → e⁺e⁻)   (e⁺e⁻ → γγ)
       (T-reversal)
```

Each edge is a crossing of the same underlying QED tree amplitude.
The amplitude is a single algebraic function `|M̄|²(s, t, u)`; each
process is its value in a different Mandelstam region. The BAM
construction provides a closed-form expression for this function in
Compton variables `(x, c)`. The question is whether the same closed
form, analytically continued, lives consistently in all three
kinematic regions.

## Crossing maps

  - **Compton → BW** (PR #36):
    `k → k₁`, `p → −p₊`, `k' → −k₂`, `p' → p₋`
    → `s_C → u_BW`, `t_C → s_BW`, `u_C → t_BW`.

  - **Compton → annihilation** (this PR):
    `k → −k₂`, `p → p₋`, `k' → k₁`, `p' → −p₊`
    → `s_C → u_ann`, `t_C → s_ann`, `u_C → t_ann`.

Both crossings use the same Mandelstam permutation `(s, t, u) →
(u, s, t)`. The physical distinction between BW and annihilation is
which two legs are incoming — but Mandelstam invariants are
time-reversal invariant, so

    |M̄|²_BW(s, t, u) = |M̄|²_ann(s, t, u)

as functions of the invariants. The cross sections differ only in
the kinematic conversion (flux factor and final-state phase space).

## Annihilation kinematics in CM

CM frame with `e⁻` along `+z`:

    p₋ = (E, 0, 0, p),   p₊ = (E, 0, 0, −p),   p = √(E² − m²) = E·β
    k₁ = (E, E sinθ, 0, E cosθ),
    k₂ = (E, −E sinθ, 0, −E cosθ),

with √s = 2E, β = √(1 − 4m²/s). Mandelstam:

    s_ann = 4E² = s
    t_ann = m² − 2E²(1 − β·cosθ)
    u_ann = m² − 2E²(1 + β·cosθ)

— identical functional form to BW kinematics (with the role of θ
reinterpreted as the photon-to-electron angle rather than the
electron-to-photon angle).

## Predictions

### P1. Mandelstam crossing identity for annihilation

    |M̄|²_ann(β, cosθ) = −|M̄|²_KN evaluated at
                          (s_C = u_ann, t_C = s_ann, u_C = t_ann)

PASS at machine precision is the standard QED prediction.

### P2. BAM kernel crossed to annihilation reproduces |M̄|²_ann

Applying the conversion factor `−2/x_⊗²` (Compton x²/2 relation +
fermion-crossing sign) to the BAM product evaluated at the
annihilation-crossed variables yields

    |M̄|²_ann_BAM = −2·(f_baseline_C · F²_C)/x_⊗²

with `x_⊗ = (m² − t_ann)/(u_ann − m²)` and
`c_⊗ = 1 + 2·s_ann·m²/[(u_ann − m²)(m² − t_ann)]`. PASS at machine
precision means the BAM closed-form F is process-general under the
annihilation crossing.

### P3. T-invariance: BAM-BW = BAM-ann as functions of (s, t, u)

Both BW and annihilation use the same Mandelstam permutation; their
crossed BAM kernels must therefore be the same function. Verify
numerically that the BAM-predicted `|M̄|²_BW` and `|M̄|²_ann` at the
same `(β, cosθ)` agree at machine precision.

### P4. Triangle loop closure

Going around the triangle should be the identity. Define the
Mandelstam permutations

    π_C→BW :  (s, t, u) → (u, s, t)     (Compton → BW)
    π_BW→ann :  identity                  (T-reversal, same Mandelstam)
    π_ann→C :  inverse of π_C→ann = (s, t, u) → (t, u, s)  (ann → Compton)

The composition `π_ann→C ∘ π_BW→ann ∘ π_C→BW` should be the identity
permutation. Verify by composing the corresponding variable
substitutions and checking that round-tripping a generic Compton
amplitude `|M̄|²_KN(s, t, u)` lands back on itself.

### P5. Differential and total annihilation cross section

The annihilation differential cross section in CM frame:

    dσ_ann/dΩ_CM = (1/(128π²·s·β)) · |M̄|²(s, t, u)
                  (with 1/2 for identical final photons)

vs the BW differential:

    dσ_BW/dΩ_CM = β·|M̄|²/(128π²·s)

so `dσ_ann = dσ_BW/β²` at the same kinematic point. Integrate to
get the Dirac annihilation cross section

    σ_ann(s) = (π·r_e²/(2β²))·(1−β²)·[(3−β⁴)·log((1+β)/(1−β))
                                       − 2β·(2−β²)]

The BAM-predicted differential, integrated over `cosθ`, must
reproduce the Dirac formula to machine precision. The threshold
divergence `σ_ann ∼ π/β` as `β → 0` (s-wave Coulomb-like
enhancement) and the high-energy logarithmic falloff must be
reproduced.

## Verdict structure

  - **TRIANGLE_CLOSES**: P1–P5 all pass at machine precision. The
    same closed-form Compton vertex factor F², expressed in
    Lorentz invariants and analytically continued via the
    appropriate Mandelstam permutations, reproduces Compton (PR #35),
    BW (PR #36), and annihilation (this PR). The three pairwise
    crossings are consistent; the loop is identity. The BAM tree
    kernel is process-general across the full QED two-photon-two-
    fermion triangle.

  - **PARTIAL_CLOSURE**: P1, P2 pass but P4 (loop closure) fails —
    indicating an inconsistency in the algebraic structure
    of the crossing permutations.

  - **ANNIHILATION_BREAKS**: P2 fails. The BAM kernel reproduces
    Compton and BW but not annihilation — a Compton-/BW-specific
    algebraic feature that doesn't survive the second crossing.

## What this leaves open

  - **Bhabha** (`e⁺e⁻ → e⁺e⁻`) and **Møller** (`e⁻e⁻ → e⁻e⁻`)
    scattering: two-channel processes (s+t or t+u) requiring
    interference of two crossed copies of the elementary vertex.
  - **Loop corrections**: still tree-level only.
  - **BAM first-principles derivation of F²**: still open.

## Cross-references

  - PR #35: `compton_vertex_resummation_probe.py` — closed-form F².
  - PR #36: `breit_wheeler_cross_process_probe.py` —
    Compton → BW edge.
  - `experiments/closure_ledger/pair_annihilation_crossing_probe.py`
    — this probe.
