# Color algebra: `SU(2) × Z₂` is the BAM-native choice (PR #80)

Final PR of the four-PR QCD-shell arc (#77 scaffold → #78
mass-ordering audit → #79 boundary stress → #80 color algebra).
Identifies the color algebra acting on the shell waveguide basis,
populates `H_couple`, settles the v3 species ↔ partition map
question, and re-audits the `n_part = 233` compensator.

## The candidates

| algebra | generators | BAM-derivable? | verdict |
|---|---:|:---:|---|
| SU(3) (canonical QCD color) | 8 | ✗ | no triplet structure with right algebra in current scaffold |
| **SU(2) × Z₂** | 4 | **✓** | **natural BAM-derivable choice** |
| Pati-Salam SU(4) | 15 | partial | requires BAM throat↔shell algebra map (open, beyond #68) |
| U(1) Cartan-only | 1 | ✓ | too weak; abelian — no inter-mode mixing strength |

**Decision: `SU(2) × Z₂`** is the BAM-native color algebra.

  - **SU(2)** from B2 / Hopf holonomy (the spin-½ structure derived
    across PRs #59–#66; `T = iσ_y`, `T² = −I`).
  - **Z₂** from PR #63's inner/outer swap (the charge-conjugation
    involution `C`).
  - Both factors derive from established BAM primitives.

## Why not SU(3)?

The natural BAM "triplet" candidates all yield SU(2)/SO(3)
algebras, **not** SU(3):

  - **3 generations from `(k_5+1)/2 = 3`** (PR #72) — gives `S_3`
    permutation / SO(3), not SU(3).
  - **Three Hopf fibrations of S³** (i, j, k quaternion axes) — SO(3)
    permutation, not SU(3).
  - **S³ isometries** — SO(4) = SU(2) × SU(2), not SU(3).
  - **Hopf bundle structure group** — U(1), not SU(3).
  - **Bulk 5D = time × radial × S³** — no natural SU(3) substructure.

The 4 generators of SU(2)×Z₂ vs SU(3)'s 8 generators is a
**substantive structural difference, not a contradiction with QCD**.
BAM's color algebra is a model of the shell-waveguide internal
symmetry, not a derivation of canonical QCD color from underlying
geometry. SU(3) would require an additional structural input outside
the current scaffold — most plausibly a quantitative Pati-Salam SU(4)
extension unifying throat-leptons and shell-quarks via a BAM-native
throat↔shell algebra map (beyond PR #68's structural transition).

## Construction on the 6-state basis

Basis ordering: `[(n=3, +), (n=3, −), (n=4, +), (n=4, −), (n=5, +), (n=5, −)]`.

SU(2) generators (Pauli matrices on partition index, per generation
block):

```
T_1 = I_3 ⊗ σ_x      T_2 = I_3 ⊗ σ_y      T_3 = I_3 ⊗ σ_z
```

Z₂ generator (n=3 ↔ n=5 swap, n=4 fixed):

```
            ⎡0 0 1⎤
P_gen   =   ⎢0 1 0⎥                Z₂ = P_gen ⊗ I_2
            ⎣1 0 0⎦
```

Algebra verification (T2 of the probe):

  - `[T_1, T_2] = 2i T_3`, `[T_2, T_3] = 2i T_1`, `[T_3, T_1] = 2i T_2`
    (SU(2)).
  - `Z₂² = I` (involution).
  - `[T_i, Z₂] = 0` for all i (product algebra structure).

## H_couple population

```
H_couple  =  α · T_3  +  β · T_1  +  γ · (Z₂ − I)
```

with `α, β, γ` real coupling constants. Hermitian by construction.
Combined with the kinetic and partition-splitter pieces:

```
H  =  H_kin(ω²(l, n))  +  H_Z2(χ_n)  +  H_couple(α, β, γ)
```

where `χ_n` is the boundary-stress derivation from PR #79.

## Singlet projection

The fully-singlet state under SU(2)×Z₂ is the **symmetric sum over
all 6 flavor states** `(1/√6)(u+d+s+c+b+t)`. The projector

```
P_S  =  |singlet⟩⟨singlet|
```

satisfies `P_S² = P_S`, Hermitian, `Tr(P_S) = 1` (1-D fully-singlet
subspace).

**Important structural point.** `P_S` does **not** commute with the
diagonal `H_kin` — the singlet is not a mass eigenstate. Mass
eigenstates are individual flavors. This is structurally expected:
physical **observables** (e.g., scattering amplitudes, total
cross-sections) are color-singlet, but the mass spectrum is
diagonalized on the flavor states themselves.

## v3 species ↔ partition map: settled

Under PR #79's uniform-sign `χ_n` reading (`+ = heavier` from
cavity-mouth boundary stress) and the SU(2)×Z₂ structure, the
revised species map is:

```
(n=3, +) = d, (n=3, −) = u
(n=4, +) = c, (n=4, −) = s
(n=5, +) = t, (n=5, −) = b
```

Each generation block has `+` heavier than `−`, consistent with the
boundary-stress derivation. This **revises** v3's `(k=1, +) = u`
(where u was lighter) but is consistent with the rest of the
shell-waveguide arc.

## `n_part` re-audit

With the full Hamiltonian populated, the eigenvalue range factor
saturates at single-digit / modest-two-digit values even for large
illustrative couplings. The observed inter-generation mass² range
factor is `~6.4·10⁹`. **No bounded algebra acting on the 6-state
shell basis can span this range** — algebra matrix elements are
bounded by their norm, so the spectrum range is bounded by `O(|α| +
|β| + |γ|)`. Even pushing to `α, β, γ ~ 10⁹` would be a
phenomenological tuning, not a structural derivation.

**The inter-generation mass hierarchy is therefore outside the scope
of BAM's color algebra on the shell waveguide.** It must come from
elsewhere — coupling to physics outside the shell sector (deeper
Tangherlini bulk modes, EW symmetry breaking, Yukawa-like couplings
to a separate sector). `n_part = 233` (PR #76) remains a
phenomenological compensator, but its **scope is now sharply
identified**: it absorbs the inter-generation hierarchy that no BAM
color algebra acting on the shell basis can naturally produce.

## Four-PR arc closure summary

| PR | what it established |
|---|---|
| #77 | shell waveguide basis + operator scaffold; quarks reframed as cavity wavefronts |
| #78 | shell basis structurally better than v3; uniform χ insufficient; `n_part` not closed at #78 alone |
| #79 | `χ_n` derived structurally from cavity-mouth boundary stress; uniform-sign, shell-suppressed; sign-flipping ansatz overruled; magnitude 30–100× too small |
| #80 | BAM-native color algebra = SU(2)×Z₂; `H_couple` populated; v3 species map settled (`+ = heavier`); `n_part` remains compensator with sharply identified scope |

**Closed:**

  - Shell waveguide basis is the right machinery.
  - Shell basis structurally distinct from v3 lepton-shaped basis.
  - `χ_n` has a no-free-parameter structural origin.
  - BAM-native color algebra identified.
  - v3 species ↔ partition map revised consistently.

**Remaining open:**

  - Inter-generation mass hierarchy (~9 orders of magnitude in mass²).
  - `n_part = 233` as a residual phenomenological compensator with
    sharpened scope.
  - Pati-Salam SU(4) throat↔shell unification (potential extension
    beyond this arc).
  - Absolute MeV scale (B4 anchor, single dimensionful input).

## Honest scope

  - **Is:** identification of BAM-native color algebra (SU(2)×Z₂);
    construction and verification of generators on the 6-state basis;
    population of `H_couple`; settlement of v3 species ↔ partition
    map; populated singlet projector; final `n_part` re-audit showing
    inter-generation hierarchy is outside BAM color scope.

  - **Is not:** a derivation of quark masses; an identification of
    SU(3) from current primitives; a closure of the inter-generation
    mass hierarchy. The four-PR arc closes structurally; the residual
    compensator has sharply identified scope rather than being
    resolved.

## B4 accounting

Algebra generators are **dimensionless**; coupling constants `α, β, γ`
are dimensionful (mass² units when paired with the `H_kin` ω² scale).
The structural identifications are scale-free. Absolute MeV scale
rides on the single B4 anchor `m_e = f_closure · ℏ / (ΔR·c)` (PR #53).

## Tests

  T1. Candidate color algebras evaluated against BAM primitives;
      SU(2)×Z₂ is the only candidate with all factors BAM-derivable.
  T2. SU(2)×Z₂ generators constructed on the 6-state basis;
      algebra verified (SU(2) commutators, Z₂ involution, product
      structure).
  T3. `H_couple = α·T_3 + β·T_1 + γ·(Z₂ − I)` populated with
      illustrative couplings; Hermitian; bounded eigenvalue spread.
  T4. Singlet projector `P_S` constructed; `P_S² = P_S`, Hermitian,
      `Tr(P_S) = 1`; does NOT commute with `H_kin` (structurally
      expected).
  T5. v3 species map revised under uniform-sign `χ_n` reading
      (`+ = heavier`); each block has + heavier than − (consistent
      with observed).
  T6. `n_part` re-audit: coupling scan from α=0 to α=100; range
      factor saturates at single-digit / modest-two-digit values;
      far from observed ×6.4·10⁹.
  T7. Why SU(3) is not derivable from current scaffold (all natural
      triplet candidates yield SU(2)/SO(3), not SU(3)).
  T8. Four-PR arc closure summary.

## Verdict structure

  - **`COLOR_ALGEBRA_SU2_Z2_BAM_NATIVE_MASS_HIERARCHY_OPEN`**
    (expected): BAM-native color algebra is SU(2)×Z₂ (from
    established primitives); the inter-generation mass hierarchy is
    outside the scope of any BAM color algebra acting on the shell
    basis; `n_part = 233` remains a phenomenological compensator with
    sharply identified scope.
  - **`COLOR_ALGEBRA_INCONCLUSIVE`**: a structural test fails;
    investigate before declaring arc closure.

## What this leaves open

  - **Inter-generation mass hierarchy.** Outside BAM color algebra
    scope. Requires Pati-Salam SU(4) extension (with throat↔shell
    algebra map) or coupling to a separate sector.
  - **`n_part = 233`** as a residual compensator with sharply
    identified scope (absorbs the inter-generation hierarchy).
  - **Absolute MeV scale** — single B4 anchor (PR #53), independent.

## Cross-references

  - `docs/qcd_shell_waveguide_scaffold_research_plan.md` — PR #77,
    the basis + operator scaffold.
  - `docs/shell_mass_ordering_audit_research_plan.md` — PR #78, the
    mass-ordering audit that posed the sign-flipping `χ_n` question.
  - `docs/boundary_stress_chi_n_research_plan.md` — PR #79, the
    cavity-mouth boundary-stress derivation of `χ_n`.
  - `docs/charge_conjugation_swap_research_plan.md` — PR #63, the
    inner/outer swap that provides the Z₂ factor.
  - `docs/throat_dirac_spinor_research_plan.md` — PR #66, the
    spin-½ structure that provides the SU(2) factor.
  - `docs/quark_npart_origin_research_plan.md` — PR #76, the
    diagnosis that motivated the four-PR arc.
  - `docs/three_throat_modes_research_plan.md` — PR #72, the
    `(k_5+1)/2 = 3` generation count.
  - `experiments/closure_ledger/color_algebra_shell_probe.py` — this
    probe.
