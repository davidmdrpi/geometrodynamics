# BAM effective-action scaffold — barrier closure status

Tracks the closure programme for the covariant BAM effective-action
scaffold. The scaffold (`bam_effective_action_scaffold_probe`) proposed
a single 5D variational principle unifying three targets — the Compton
vertex `F²(x, c)`, the Hopf-bundle U(1) connection `A_φ = ½cos χ`, and
the 5D Tangherlini bulk boundaries `ΔR` — and identified **five
mismatch terms** (B1–B5) blocking full closure. This document records
how each barrier has since been addressed.

## The candidate action

```
S_BAM = ∫_{M₅} √(−g₅)[ (R₅−2Λ₅)/2κ₅ − ¼ F_{MN}F^{MN}
                       + ψ̄(iΓ^M D_M − m)ψ + L_throat ]
      + S_∂[hard walls] + S_closure[∮A = 2πn]
```

Three sectors close the three targets:

| sector | target | status |
|---|---|---|
| (A) U(1) gauge | `A_φ = ½cos χ` | closes outright (`c₁ = 1`, no input) |
| (T) throat | `F² = K²·Q` | closes given the topological sector |
| (G) gravity | `ΔR = R_OUTER − R_INNER` | metric bulk-derived; boundaries from BCs |

## Barrier ledger

| barrier | type | original status | now | by |
|---|---|---|---|---|
| **B1** closure quantum `∮A = 2πn` | topological | imposed constraint | **CLOSED** | winding θ-term `S_top = 2π·n` (topological-discrete sector probe) |
| **B2** antipodal Z₂ `T = iσ_y` | discrete | imposed identification | **CLOSED** | `RP³ = S³/Z₂` + non-trivial spin structure (topological-discrete sector probe) |
| **B3** hard-wall BC (Dirichlet at throat) | boundary | imposed by hand | **CLOSED** | single-valuedness under `T² = −I` forces `ψ(throat) = 0` (hard-wall boundary derivation probe) |
| **B5** 5D→4D reduction producing F² | reduction | unconstructed | **CLOSED** | one separable master functional `ℳ = G_C ⊗ 𝒢_{S³}` on the warped product `C × S³` yields masses (poles), K (throat boundary), Q (S³ Hopf); vertex residue = `F²=K²·Q` (master integral probe) |
| **B4** dimensional bridge `ℏ = m_e·R_MID·c` | scale | one external anchor | **IRREDUCIBLE** | the closure-ledger/Maslov machinery is scale-free, so exactly one dimensionful anchor is mathematically required; `m_e` is the minimal choice (Maslov dimensional-bridge probe) |

## How each barrier was addressed

### B1 + B2 → the topological/discrete sector (CLOSED)

Both barriers are data of a single topological/discrete sector:

```
RP³ = S³/Z₂  +  non-trivial spin structure (T² = −I)  +  2π winding θ-term
```

  - **B2** is the deck transformation of the double cover `S³ → RP³`
    (`σ: p → −p`, a free involution; `π₁(RP³) = Z₂`), and `T² = −I`
    selects the non-trivial of RP³'s two spin structures
    (`H¹(RP³, Z₂) = Z₂`; antiperiodic 4π-spinors). Topological data,
    not an imposed symmetry.
  - **B1** is the winding number of the phase map `S¹ → U(1)`:
    `∮dφ = 2π·n`, a topological total-derivative term `S_top = 2π·n`
    (integer-quantized, metric-independent, doesn't modify the local
    EOM → variationally consistent).
  - The two are unified by the double cover: a great circle on `S³`
    (length 2π = closure quantum, B1) is a double traversal of the
    non-contractible `RP³` loop (the `π₁` generator, B2).

With both as action data, `K(x) = 2x/(1+x)` and `Q(x, c)` follow from
the topological sector + stationary action, no longer imposed.

### B3 → consequence of the spin structure (CLOSED)

The hard-wall (Dirichlet) BC at the throat is **not** an independent
imposition once the non-trivial spin structure (B2) is in place:

```
T = iσ_y, T² = −I  (eigenvalues ±i, no +1 eigenvector)
   →  single-valuedness at the throat fixed point: ψ = T·ψ
   →  T²ψ = Tψ = ψ but T²ψ = −ψ  ⟹  ψ = −ψ  ⟹  ψ(throat) = 0
```

Realized concretely in the Tangherlini radial solver, whose modes are
extended antisymmetrically across the throat
(`u_full = [−u_reflected, u]`, odd to machine precision), producing the
Dirichlet node. Neumann/Robin are ruled out (no nonzero `T`-invariant
spinor). B3 is absorbed into the topological sector.

### B5 → master integral on the warped product (CLOSED)

The 5D → 4D reduction factorizes into three channels of one action on
the shared internal geometry:

| channel | integrate over | produces |
|---|---|---|
| radial | `r ∈ [R_MID, R_OUTER]` | KK masses `ω(l,n)` |
| S³ angular | `Ω ∈ S³` | gauge `c₁ = 1` + propagator `1/q²`; helicity `Q(x,c)` |
| throat | `r → R_MID` pinch | impedance `K(x)` |

All three share `R_MID`, the closure quantum `2π`, and `T² = −I`,
connecting the mass and amplitude sub-threads.

**Central finding (radial reduction bridge):** `F²(x, c)` is **NOT** a
radial overlap integral — radial overlaps are kinematics-independent
constants (`δ_mn`), while `F²(x, c)` varies strongly with the
scattering kinematics. So `F²` is the throat-channel form factor, not a
KK overlap; the naive "F² from radial integration" is falsified.

**Resolution (master integral).** The three channels are unified in a
single separable functional on the warped-product internal geometry
`M_int = C × S³` (`C` = radial cavity `[R_MID, R_OUTER]`):

```
ℳ(ω; x, c) = G_C(r, r′; ω) ⊗ 𝒢_{S³}(Ω, Ω′)
```

read three ways from one object:

  - **poles in ω** → the mass spectrum `ω(l,n)` (radial ladder `n` ×
    S³ Casimir `l`, the latter the centrifugal term of the warp);
  - **throat boundary of `G_C`** → `K(x) = 2x/(1+x)` (dwell-time
    impedance `Z(ω)=π/ω` in series);
  - **S³ Hopf reduction of `𝒢_{S³}`** → `Q(x,c) = x²+x(1−x)²/(1+c²)`
    (Hopf-fibre helicity spinor; `(1+c²)/2 = cos⁴(θ/2)+sin⁴(θ/2)`).

The **vertex residue** reproduces `F²(x,c) = K²·Q` to machine precision
while the **poles** give the masses — masses AND the F² vertex from one
functional. The `F²=K²·Q` factorization is **not** a failure to unify:
it is the direct consequence of the product internal geometry
(separation of variables `Ψ = Σ u_{l,n}(r)𝒴_l(Ω)` makes the Green
function a sum of factor products). B5′ is **closed**.

### B4 → dimensional bridge (IRREDUCIBLE)

The action fixes dimensionless structure (`F²`, `c₁`, the closure
quantum, `K`, `Q` are all scale-free); the absolute MeV scale enters
only through the single anchor `ℏ = m_e·R_MID·c`. The Maslov
dimensional-bridge probe audits whether this last anchor can be derived
away, and shows it **cannot** — provably:

  - **The closure-ledger/Maslov machinery is scale-free.** Every
    quantity it produces is dimensionless: the winding integers `k, n`,
    the Maslov index `μ`, action ratios `S/2π`, `ω·R_MID`, mass *ratios*,
    and the geometric residuals (`R_OUTER`, `ε = 7π/(100·5⁴)`,
    transport `= 8π`, resistance `= 7π/100`, the `1.054` factor). This
    is shown concretely by **scale invariance**: rescaling
    `R_MID → λ·R_MID` leaves every dimensionless output unchanged to
    machine precision and shifts only the absolute scale.

  - **A scale-free theory cannot produce a dimensionful scale.** By
    dimensional analysis, exactly **one** external dimensionful anchor
    is mathematically required; `m_e` (equivalently `R_MID` via
    `ℏ = m_e·R_MID·c`, dimensionally consistent) is the minimal choice.

So B4 is **irreducible by dimensional necessity** — a structural
feature, not an unsolved gap (cf. SI fixing `c, ℏ, e` by convention).
The Maslov reading also ties B4 to B3: the radial ledger integer `n+1`
is the Maslov index of the doubly-Dirichlet throat cavity (`μ = 4`),
whose throat half (`μ = 2`, reflection phase `π`) is the B3 hard wall —
and that `π` is the closure half-quantum and the dwell phase `ω·τ`
(`τ = π/ω`). The companion ℏ-origin thread
(`docs/hbar_origin_status.md`) had already reduced every dimensionless
residual to closure-quantum form; this audit closes the interpretation.
What remains genuinely open is `m_e` itself, which by the audit cannot
come from scale-free geometry.

**The anchor as a geometric invariant (ΔR).** The single required
anchor need not be the particle mass `m_e` — it can be the **invariant
bulk separation** `ΔR = R_OUTER − R_INNER = 0.52·R_MID`, *provided* ΔR
is a proper (cosmologically fixed) length. The ΔR scale-modulus probe
shows it is: the throat is a static bound vacuole (discrete spectrum +
vacuum + dimensionless BC), so it is decoupled from Hubble flow
(`ΔR/R_cosmo ~ 10⁻³⁹`), and a comoving throat (`rs ∝ a`) is
observationally excluded — it would redshift particle masses as `(1+z)`
against quasar bounds `≲10⁻⁵`. The bridge then reads
`m_e = f_closure·ℏ/(ΔR·c)` with `f_closure = ΔR/R_MID = 0.52`, making
`m_e` a consequence of a fixed bulk length and predicting that local
throat ratios (lepton mass ratios included) are constant in cosmic time
while only `ΔR/R_cosmo(t) ∝ 1/a` drifts. This **relocates** the anchor
to a geometric invariant; it does not evade the scale-modulus theorem
(ΔR is still the one external dimensionful number, its value underived).

**The anchor as a self-energy equilibrium.** The
`self_consistent_throat_radius_probe` recasts the (previously imposed)
throat radius as a **finite-self-energy stable equilibrium**: the throat
caps the EM field at `R_MID` (no `r < R_MID`), making the self-energy
finite (`U_EM/(m c²) = α/2`, no UV divergence — unlike a point charge),
and `E(R) = A/R + B·R²` (EM repulsion vs cohesion) has a unique stable
minimum `R* = (A/2B)^{1/3}`. Consistent with the theorem, `R*` rides on
one dimensionful coupling (`B → B/λ³` ⟹ `R* → λ R*`); the BAM-native
balance `m c² = U_EM` is `R`-independent and instead fixes `g = 2/α`,
relocating the scale question to `α`. The chain is therefore: imposed
`R_MID` → invariant geometric length `ΔR` → finite-self-energy
equilibrium — each step more physical, none deriving the absolute value.

The cohesive `B·R²` (posited in the equilibrium) is then **derived** by
`cohesive_tension_derivation_probe`: it is the throat **brane tension**
`E = σ·Area = 4πσR²` (so `B = 4πσ`), the `R²` power being the unique
signature of a constant surface tension by power-counting — distinct
from the induced Tangherlini junction tension (`R¹`, computed from
`f(r)=1−(rs/r)²`), Einstein–Hilbert (`R¹`), and the cosmological bag
(`R³`). The tension is set by the bulk gravity sector (`σ ∝ √|Λ₅|/κ₅`,
Randall–Sundrum-like), so its *value* is the single dimensionful anchor:
the derivation fixes the cohesive term's form and identity, not the
absolute scale.

`brane_tension_tuning_probe` sharpens that bulk-gravity relation to the
**exact** RS fine-tuning. The `Z₂` Israel junction for a pure-tension
brane gives `K_μν = −(κ₅²λ/6) h_μν`; the bulk `AdS₅` equation
`Λ₅ = −6k²` and staticity (`K_μν = k h_μν`) then fix
`λ_crit = 6k/κ₅² = √(6|Λ₅|)/κ₅²` — the **dimensionless tuning factor is
√6**. The tuning is the flat / static-throat condition (`Λ₄ ∝ λ²−λ_crit²`
vanishes at `λ_crit`; over-/under-tension give dS/AdS throats), tying the
critically-tuned brane to the static equilibrium of `#55`. The fine-tuning
is one condition among `(λ, Λ₅, κ₅)`, so a net one dimensionful
combination remains — the single anchor (the bulk gravitational scale
`k = √|Λ₅/6|`); `√6` and the flatness condition are derived, the scale is
not. The `AdS₅` warp over the bulk depth `ΔR` gives an RS exponential
hierarchy `e^{−kΔR}`.

Finally, `pair_production_threshold_probe` shows the **pair-production
threshold** falls out as the lowest stable configuration: a throat
carries one Hopf charge (`|c₁| = 1`), so conservation (`Σ c₁ = 0`) forces
creation as a C-conjugate throat–antithroat pair (the antipodal `Z₂`,
B2), with threshold `E_thr = 2 E(R*) = 2 m_e c² = 1.022 MeV`. A
bubble-nucleation barrier (`R_c = 2σ/ρ`, the brane tension `σ` as the
surface cost) gives the *disperse-below / persist-above* dichotomy, and
the Schwinger critical field `E_S = m_e²c³/(eℏ)` (where
`e E_S R_MID = m_e c²`) ties the throat scale to the threshold. The
factor 2 and the structure are derived; the absolute `2 m_e c²` rides on
the single anchor `m_e c² = ℏc/R_MID`.

## Summary

The scaffold began with five barriers and is now **complete**. Four
(B1, B2, B3, B5) are **closed**: B1, B2, B3 promoted to a topological/
discrete action sector (`RP³ + spin structure + winding θ-term`) whose
spin structure also forces the hard-wall BC; B5 closed by the master
integral (one `C × S³` functional yielding masses, `K`, and `Q`, with
vertex residue `F²=K²·Q` — the factorization being the product-geometry
consequence). The fifth (B4) is **irreducible by dimensional
necessity**: the closure-ledger/Maslov machinery is scale-free, so
exactly one external dimensionful anchor is mathematically required,
and `m_e` is the minimal choice.

The programme is therefore complete in a precise sense: four barriers
derived from geometry, the fifth shown to be the single mandatory
dimensionful unit. The only genuinely open input is `m_e` itself, which
the B4 audit shows cannot come from scale-free geometry.

## Probe ledger

| probe | barriers | verdict |
|---|---|---|
| `bam_effective_action_scaffold_probe` | B1–B5 (map) | SCAFFOLD_WITH_BARRIERS |
| `topological_discrete_sector_probe` | B1, B2 | PROMOTION_SUCCEEDS |
| `hard_wall_boundary_derivation_probe` | B3 | HARD_WALL_DERIVED |
| `radial_reduction_bridge_probe` | B5 | BRIDGE_FACTORIZED |
| `bulk_boundary_interaction_probe` | B5′ (radial+throat) | BULK_BOUNDARY_FORMULATED |
| `master_integral_probe` | B5′ (+ S³ Q) | MASTER_INTEGRAL_COMPLETE |
| `maslov_dimensional_bridge_probe` | B4 + Maslov ledger | B4_IRREDUCIBLE |
| `delta_r_scale_modulus_probe` | B4 anchor (ΔR) | DELTA_R_INVARIANT |
| `self_consistent_throat_radius_probe` | B4 anchor (self-energy) | SELF_CONSISTENT_THROAT_EQUILIBRIUM |
| `cohesive_tension_derivation_probe` | B4 anchor (cohesive term) | COHESIVE_TENSION_DERIVED |
| `brane_tension_tuning_probe` | B4 anchor (bulk-gravity tuning) | BRANE_TUNING_DERIVED |
| `pair_production_threshold_probe` | B4 anchor (pair threshold) | PAIR_THRESHOLD_DERIVED |
| `stable_moving_throat_probe` | throat = particle (Lorentz) | MOVING_THROAT_COVARIANT |
| `spin_wigner_rotation_probe` | throat = spin-½ (Wigner) | SPIN_WIGNER_COVARIANT |
| `gyromagnetic_ratio_probe` | throat g = 2 (magnetic moment) | G_FACTOR_DERIVED |
| `throat_vertex_loop_probe` | throat g−2 = α/2π (one loop) | SCHWINGER_RECONSTRUCTED |
| `charge_conjugation_swap_probe` | C = inner/outer swap (c₁→−c₁) | C_IS_INNER_OUTER_SWAP |
| `cpt_assembly_probe` | CPT = C·P·T (throat histories) | CPT_ASSEMBLED |
| `cpt_dirac_operator_probe` | CPT operator Θ ∝ γ⁵ (Dirac spinor) | CPT_OPERATOR_CONSTRUCTED |
| `throat_dirac_spinor_probe` | throat 4-spinor from S_BAM (Dirac factorization) | THROAT_DIRAC_DERIVED |
| `even_k_absence_probe` | even-k absence (spin-statistics) | EVEN_K_EXCLUDED_BY_SPIN_STATISTICS |
| `throat_to_shell_transition_probe` | lepton throat → QCD shell channel | THROAT_TO_SHELL_TRANSITION_CONFIRMED |
| `shell_to_qcd_match_probe` | shell ↔ QCD structural invariants | SHELL_REPRODUCES_QCD_STRUCTURE |
| `three_generation_boundary_probe` | sharp `k ≤ 5` three-generation boundary | THREE_GENERATIONS_PINNED |
| `beta_lepton_derivation_probe` | `β_lepton = k_5²·(2π) = 50π` | BETA_LEPTON_DERIVED |
| `three_throat_modes_probe` | `#gen = (k_5+1)/2 = 3` from `k_5` | THREE_GENERATIONS_FROM_K5 |
| `k5_origin_probe` | `k_5 = D_bulk = dim(S³)+2 = 5` | K_5_FROM_BULK_DIMENSION |
| `s_bam_loop_measure_probe` | `1/(2π)` in `a = α/(2π)` = BAM closure quantum | LOOP_MEASURE_IDENTIFIED |
| `quark_npart_origin_probe` | `n_part = 233` (quark) = phenomenological compensator | N_PART_IS_PHENOMENOLOGICAL_COMPENSATOR |
| `qcd_shell_waveguide_scaffold_probe` | shell waveguide basis + operator scaffold (PRs #77–#80 arc) | SHELL_WAVEGUIDE_SCAFFOLD_CONSTRUCTED |
| `shell_mass_ordering_audit_probe` | shell mass-ordering / `n_part` audit on PR #77 basis | SHELL_BASIS_STRUCTURALLY_BETTER_N_PART_NOT_YET_RESOLVED |
| `boundary_stress_chi_n_probe` | `χ_n` derived from cavity-mouth boundary stress; singlet placeholder | CHI_N_DERIVED_BOUNDARY_STRESS_INSUFFICIENT |

## Cross-references

  - `docs/bam_effective_action_scaffold_research_plan.md` — the original
    5-barrier map.
  - `docs/topological_discrete_sector_research_plan.md` — B1 + B2.
  - `docs/hard_wall_boundary_derivation_research_plan.md` — B3.
  - `docs/radial_reduction_bridge_research_plan.md` — B5.
  - `docs/bulk_boundary_interaction_research_plan.md` — B5′ (radial+throat).
  - `docs/master_integral_research_plan.md` — B5′ closed (+ S³ Q).
  - `docs/maslov_dimensional_bridge_research_plan.md` — B4 audit +
    Maslov closure-ledger (B4 irreducible).
  - `docs/delta_r_scale_modulus_research_plan.md` — ΔR invariant under
    S³ expansion; the B4 anchor as a geometric invariant.
  - `docs/self_consistent_throat_radius_research_plan.md` — the B4 anchor
    as a finite-self-energy equilibrium.
  - `docs/cohesive_tension_derivation_research_plan.md` — the cohesive
    `B·R²` term derived as the throat brane tension.
  - `docs/brane_tension_tuning_research_plan.md` — the RS-like
    bulk-gravity fine-tuning (factor `√6`).
  - `docs/pair_production_threshold_research_plan.md` — the
    pair-production threshold `2 m_e c²` (lowest stable throat pair).
  - `docs/stable_moving_throat_research_plan.md` — the boosted throat /
    Lorentz-covariance falsifier (throat = particle).
  - `docs/spin_wigner_rotation_research_plan.md` — the Hopf-spin / Wigner-
    rotation falsifier (throat = spin-½ particle).
  - `docs/gyromagnetic_ratio_research_plan.md` — `g = 2` from the
    Pauli/SU(2) + Hopf monopole (the magnetic moment).
  - `docs/throat_vertex_loop_research_plan.md` — the one-loop Schwinger
    anomaly `a = α/2π` from the throat-vertex loop (reconstruction).
  - `docs/charge_conjugation_swap_research_plan.md` — C = the inner/outer
    swap (`c₁ → −c₁`); charge conjugation as geometry.
  - `docs/cpt_assembly_research_plan.md` — C·P·T assembled into the
    geometric CPT symmetry on throat histories.
  - `docs/cpt_dirac_operator_research_plan.md` — the explicit CPT operator
    `Θ = γ⁰γ¹γ²γ³ = −iγ⁵` on the throat Dirac spinor.
  - `docs/throat_dirac_spinor_research_plan.md` — the throat Dirac 4-spinor
    derived from the radial Dirac/SUSY factorization of `S_BAM`.
  - `docs/even_k_absence_research_plan.md` — even-k absence as a
    spin-statistics selection rule (charged leptons = odd-k fermions).
  - `docs/throat_to_shell_transition_research_plan.md` — higher
    excitations delocalize from the lepton throat into the QCD shell
    channel (focused pulse → wavefront).
  - `docs/shell_to_qcd_match_research_plan.md` — shell modes reproduce
    the documented structural invariants of the quark sector (Z₂
    partition, 3×2=6 flavors, heavier scale, extended character).
  - `docs/three_generation_boundary_research_plan.md` — the sharp
    `k ≤ 5` boundary from β-uplift quadratic growth + throat-shell
    availability (combining #67–#69).
  - `docs/beta_lepton_derivation_research_plan.md` — `β_lepton =
    k_5²·(2π) = 50π` from closure-quantum primitives + topological
    charge; closes the PR #70 follow-on.
  - `docs/three_throat_modes_research_plan.md` — `#generations =
    (k_5+1)/2 = 3` from the same `k_5` primitive (closes the "why 3
    throat modes" follow-on).
  - `docs/k5_origin_research_plan.md` — `k_5 = D_bulk = dim(S³)+2 = 5`
    (reduces "why k_5 = 5" to "why the Hopf bundle / S³").
  - `docs/s_bam_loop_measure_research_plan.md` — the `1/(2π)` in the
    Schwinger anomaly `a = α/(2π)` identified as the BAM closure quantum
    (same `2π` as `action_base`, closure ledger, `β_lepton`, Hopf, throat
    dwell, `ε` integer); closes the structural piece of PR #62's open
    follow-on. Full covariant `(2π)^d` path-integral derivation remains
    future work.
  - `docs/quark_npart_origin_research_plan.md` — `n_part = 233` (quark)
    classified as a phenomenological compensator at the v3 baseline;
    extended catalog (Fibonacci, Lucas, Padovan, Perrin, tribonacci,
    color × flavor × generation, QCD β₀, Tangherlini QCD-shell modes)
    yields no exact match surviving §8 drift; structural reading is
    "v3 Hamiltonian is lepton-shaped, quark sector lives in QCD shell
    channel per #68–#69". Right derivation route (quantitative #68–#69)
    is outside closure-ledger scope.
  - `docs/qcd_shell_waveguide_scaffold_research_plan.md` — PR #77, the
    foundation of the four-PR QCD-shell arc (#77 scaffold → #78 mass-
    ordering / `n_part` audit → #79 boundary stress tensor + singlet
    constraint → #80 color algebra). Quarks reframed as cavity
    wavefronts that resolve the shell, NOT throat traversals. 6-state
    `(l, n, p)` basis with 6×6 operator scaffold `H = H_kin + H_Z2 +
    H_couple`.
  - `docs/shell_mass_ordering_audit_research_plan.md` — PR #78, the
    mass-ordering / `n_part` audit on the PR #77 shell basis. Finds:
    (i) shell basis is structurally better than v3 (cavity wavefronts;
    `ω²(l, n)` kinetic; Z₂ partition slot for the within-generation
    inversion); (ii) uniform `χ·σ_z` cannot reproduce the inversion
    (best 2/3 blocks); (iii) sign-flipping `χ_n` can (existence
    proof); (iv) shell kinetic spans ×2.2 in mass² vs observed
    ×6.4·10⁹ — ~9 orders unaccounted for; (v) `n_part` NOT resolved
    at PR #78 alone — depends on PR #79's `χ_n` derivation and PR
    #80's `H_couple` population.
  - `docs/boundary_stress_chi_n_research_plan.md` — PR #79, derives
    `χ_n = T_odd(n) = (T_inner − T_outer)/2` from the Z₂-antisymmetric
    piece of cavity-mouth boundary stress (PR #63's inner/outer swap).
    Findings: (i) `χ_n` structurally pinned with no free parameter
    once cavity geometry fixed; (ii) sign is uniform-positive across
    all `n` (no sign flip), overruling PR #78's sign-flipping ansatz;
    (iii) magnitude is shell-suppressed (`χ_n/ω² ~ 0.01–0.02` for
    `n ≥ 3`), 30–100× too small for observed within-generation
    splittings; (iv) within-generation inversion and inter-generation
    hierarchy ⟹ PR #80 color sector; (v) v3 species ↔ partition map
    flagged for revision; (vi) singlet projector placeholder (identity
    on flavor basis), awaits PR #80 color algebra.
  - `docs/odd_k_closure_lemma.md` — the closure arithmetic this upgrades.
  - `docs/hbar_origin_status.md` — B4 (the m_e anchor).
  - `docs/tree_qed_status.md` — the tree-QED result the F² target
    summarises.
