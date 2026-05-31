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
| `color_algebra_shell_probe` | BAM-native color algebra = SU(2)×Z₂; H_couple populated; v3 species map settled | COLOR_ALGEBRA_SU2_Z2_BAM_NATIVE_MASS_HIERARCHY_OPEN |
| `pati_salam_throat_shell_bridge_probe` | throat ↔ shell n+3 Z₂ bridge built; 3 open extensions for full SU(4) | PATI_SALAM_THROAT_SHELL_Z2_BUILT_FULL_SU4_REQUIRES_EXTENSIONS |
| `throat_shell_mass_operator_unification_probe` | lepton β·k² and quark ω²(l,n) unified as one Bohr-Sommerfeld operator m²=(S/L)² | MASS_OPERATOR_UNIFIED_BOHR_SOMMERFELD |
| `winding_shell_quadrant_probe` | (k≠0, n≥3) quadrant = leptoquark sector; complete four-quadrant interpretation | WINDING_SHELL_QUADRANT_IS_LEPTOQUARK_SECTOR |
| `neutrino_quadrant_suppression_probe` | neutrino = Majorana (k=0 ⟹ c₁=0 ⟹ C-invariant); seesaw mechanism, M_R scale open | NEUTRINO_SUPPRESSION_IS_MAJORANA_SEESAW_SCALE_OPEN |
| `seesaw_scale_nucleation_compliance_probe` | M_R grounded in PR #58 throat↔antithroat nucleation; Σc₁=0 = only-neutrino rule; barrier-height M_R falsified; suppression = tunnelling, M_R = m_D·e^{S}, S≈15–18 open | SEESAW_SUPPRESSION_IS_THROAT_ANTITHROAT_TUNNELING_S_OPEN |
| `majorana_bounce_action_probe` | reduced Euclidean bounce on the non-orientable tortoise path; rigid throat ⟹ massless ν; S ∝ ln(1/ε) (O(10), gen-stable); EM-throat tension under-produces S ~40×; S≈15–18 needs ΔL=2 tension ratio t≈6–12 | MAJORANA_BOUNCE_IS_TORTOISE_LOG_S_OPEN_AS_TENSION_RATIO |
| `b_minus_l_tension_ratio_probe` | ΔL=2/B−L tension ratio t = global-closure enhancement of local EM tension; bracketed parameter-free by closure quantum 2π (lower) and winding action k_5√(2π)=√β (upper): t∈[6.28,12.53], matching PR #88's 6–12; residual = compliance ε | B_MINUS_L_TENSION_BRACKETED_BY_CLOSURE_AND_WINDING_ACTIONS |
| `boundary_compliance_bulk_geometry_probe` | ε = chargeless-throat sub-throat healing length (ε=ℓ²/2rs; c₁=0 neck not EM-propped, charged → Dirac); bulk scales (R_c³,Δ³) land ε in window; winding-edge t≈√β ⟹ S≈15–19, m_ν~few meV (observed scale, untuned); precise spectrum residual | COMPLIANCE_IS_CHARGELESS_THROAT_HEALING_LENGTH_CHAIN_CLOSED_TO_OOM |
| `epsilon_bulk_compliance_probe` | is ε computed from bulk compliance or inferred from meV? COMPUTED (meV-free): healing length ℓ~R_c=2σ/ρ from the ELECTRON calibration (PR #58, R_c=2/9) ⟹ ε~R_c³≈0.011 sub-throat O(10⁻²), no neutrino input; with t=k_5√(2π)=√β_lepton ⟹ S≈16.85 ⟹ m_ν≈2.1 meV (scale output/retrodiction), deriving the exponential smallness (ε≪1⟹large S⟹tiny m_ν). RESIDUAL: precise ε — m_ν∝ε^4.8 ⟹ O(1) ambiguity (R_c³→2, Δ³→20, R_c²/2→108 meV) spans ×50; absolute compliance normalization = unpinned κ₅²/Λ₅ (only √6 fixed, PR #57). Smallness derived from bulk compliance; exact value not | EPSILON_BULK_COMPLIANCE_DERIVES_SMALLNESS_PRECISE_VALUE_STAYS_RESIDUAL |
| `generation_spread_pmns_mixing_probe` | generations = cavity overtones ⟹ bare m_ν ∝ m_D (normal ordering 1:1.87:2.74); spread widened by overtone-dependent neck coupling (PR #79 χ_n ↓ with n ⟹ higher-n less suppressed ⟹ heavier); large PMNS = cross-channel (charged k≠0 × neutrino k=0), small CKM = intra-channel (shell × shell) ⟹ PMNS ≫ CKM; angles/spectrum open | PMNS_CROSS_CHANNEL_CKM_INTRA_CHANNEL_SPREAD_FROM_OVERTONE_COUPLING |
| `cross_channel_pmns_overlap_probe` | naive radial overlap → near-permutation (small); lepton gens in different coordinates (closure-winding k vs radial-overtone n) ⟹ anarchic map; observed PMNS typical of Haar U(3) (30th/57th/4th pct), CKM extremely atypical (aligned, joint p≈0); specific angles open (θ13 mild tension) | PMNS_ANARCHIC_CROSS_COORDINATE_CKM_ALIGNED_INTRA_COORDINATE |
| `theta13_residual_alignment_probe` | θ13=U_e3 is the corner / two-hop element (gap |g−i|=2); residual nearest-neighbour alignment (throat↔shell coupling local in (k,n)) suppresses it ⟹ θ13 robustly smallest, observed θ13 moves 4th→~21st percentile (PR #92 tension resolved), θ12/θ23 stay typical; exact θ13 (μ one param, median saturates ~14–16°) open | THETA13_SUPPRESSED_BY_RESIDUAL_NEAREST_NEIGHBOUR_ALIGNMENT |
| `cp_majorana_phase_probe` | CP violation generic (winding amplitudes Hopf-complex e^{ikχ}, PR #60; CP-conservation measure-zero); Jarlskog dichotomy: \|J_PMNS\|≈0.026 typical of anarchy (51st/81st pct), \|J_CKM\|≈3e-5 extremely atypical (aligned, suppressed); two Majorana phases EXIST ⟸ neutrino Majorana ⟸ c₁=0 (PR #86), 0νββ; specific phase values anarchic/not pinned | CP_GENERIC_MAJORANA_PHASES_EXIST_VALUES_ANARCHIC |
| `zeronubb_effective_mass_probe` | 0νββ occurs ⟸ Majorana (c₁=0, PR #86); BAM normal-ordering band (PR #91) below IO floor ~19 meV; anarchic Majorana phases (PR #94) populate full band incl. cancellation→~0; light scale (PR #90) ⟹ m_ββ ≲ 8 meV — below current (28–122 meV) & next-gen (~9–20 meV) reach; falsifiable (discovery ≳19 meV ⟹ IO/degenerate); exact m_ββ a band | ZERONUBB_OCCURS_NORMAL_ORDERING_M_BB_FEW_MEV |
| `cosmological_sigma_mnu_probe` | same light, normal-ordered spectrum ⟹ Σm_ν ≈ 59–65 meV (NO floor 58.7 meV, below IO floor 99 meV); consistent with Planck (<120), at the DESI DR2+CMB frontier (~60–64 meV); falsifiable (Σ<58.7 ⟹ NO excluded; Σ≳100 ⟹ not light); cross-checks PR #95 (m_ββ ≲ 8 meV) | SIGMA_MNU_AT_NORMAL_ORDERING_FLOOR_60MEV |
| `neutrino_mev_scale_sharpening_probe` | sharpens the #96 band into a PINNED spectrum: NuFIT 6.0 fixes m₂=8.65, m₃=50.34 meV (NO floor 59.0); DESI DR2+CMB (≲60–64) corners m₁≲3 meV ⟹ Σm_ν∈[59.0,62.6] (tightened from 59–65, toward the floor); pinned spectrum m=(≲3,8.65,50.34) meV; m_β≈8.8–9.3 meV; m_ββ NONZERO floor [1.5,3.7] meV (no full cancellation in NO: s12²c13²m₂=2.60 > s13²m₃=1.10); honest reachability — only Σm_ν near-term testable (DESI at floor now), m_β ~4–5× below Project 8, m_ββ ~3–10× below LEGEND-1000/nEXO; flag: some 2025 DESI+CMB fits prefer Σ at/below floor ⟹ tension for all NO models; open: m₁ band (0–3 meV) + anarchic Majorana phases | NEUTRINO_MEV_PREDICTIONS_SHARPENED_SPECTRUM_PINNED_SIGMA_TESTABLE |
| `npart_dynamical_hierarchy_probe` | n_part=233 revisited: a huge hierarchy CAN be geometric (neutrino e^{−S}), so size isn't the obstruction; the quark hierarchy is IRREGULAR (c/u≈588 vs t/c≈136, up/down asym) — the QCD-RG signature; geometric shell span ×2.2 vs observed ×6.4×10⁹; quark is the program's one dynamical sector; gap N_q−N_lepton=366 = dynamical excess; PR #76 upheld+sharpened | QUARK_HIERARCHY_DYNAMICAL_N_PART_COMPENSATES |
| `quark_hierarchy_flavor_puzzle_probe` | refines #97: quark mass RATIOS are RG-invariant (γ_m flavor-universal) ⟹ hierarchy is NOT α_s running but the FLAVOR PUZZLE (Yukawas); quark Yukawas overflow the compressed shell capacity (×1.49) by ~×5×10⁴; BAM captures quark STRUCTURE (6=3×2, Z₂, k=0, 3 gens) but not the magnitudes; #97 core (dynamical/non-geometric) stands | QUARK_HIERARCHY_IS_FLAVOR_PUZZLE_NOT_RG_RUNNING |
| `qcd_confinement_cornell_audit_probe` | confinement geometry audit: Cornell V(L)=σL−A·ℏc/L (linear=flux-tube wormhole bridge, Coulomb=throat/gluon exchange); string breaking = Schwinger exp(−πm_q²/(σL)) = the PR #58 throat-pair mechanism with eE→σ; BAM σ reproduces Regge α'=1/(2πσ)=0.884 GeV⁻² and L_break; √σ≈0.42 GeV = the one QCD scale anchor (B4 analogue), calibrated not derived | CONFINEMENT_GEOMETRIC_STRING_BREAK_IS_SCHWINGER_SCALE_IS_QCD_ANCHOR |
| `glueball_closed_flux_loop_probe` | glueballs = pure-confinement benchmark (closed flux loops, no quark/flavor input); BAM orientable ground √(4πσ)≈1.50 GeV (3.5√σ) benchmarks lattice 0++ (4.1√σ) to ~13%; closed-string glueball Regge slope = half the meson; BAM non-orientable Möbius sector ⟹ extra glueball tower (half-int modes, +πσ in M²) interleaving the orientable one (≈2× states); legitimate vs lattice not experiment (glueballs unobserved) | GLUEBALL_SCALE_BENCHMARKS_LATTICE_MOBIUS_TOWER_IS_BAM_TOPOLOGY |
| `mobius_exotic_sector_probe` | flux-network topology = hadron taxonomy (meson/baryon/tetraquark/pentaquark/hybrid/glueball + Möbius Z₂); non-orientable (Möbius) flux tube carries the exotic J^PC (1-+) forbidden to ordinary qq̄; observed 1-+ hybrids π₁(1600), η₁(1855) match at right J^PC and at ρ/ω+2√σ (≈1.62, 1.85 GeV); multiquark exotics (X,Z_c,T_cc,P_c)=multi-junction networks; unlike glueballs, exotics observed ⟹ BAM topology meets data and matches | MOBIUS_FLUX_GIVES_EXOTIC_JPC_OBSERVED_HYBRIDS_MATCH |
| `baryonic_exotics_classification_probe` | BAM baryonic exotics (Möbius/hybrid baryon) have NO exotic-J^P smoking gun (any J^P ordinary for qqq, no C) ⟹ supernumerary ordinary-J^P states (signature = counting); land in light N*/Δ* (~1.79, 2.08 GeV = base+2√σ); constraint ranking light N*/Δ* > strange hyperons > heavy baryons; Möbius doubling must coincide w/ observed states or decouple else excluded; MOST-constrained corner (opposite of glueballs) | BARYONIC_EXOTICS_LACK_EXOTIC_JP_LIGHT_SECTOR_MOST_CONSTRAINED |
| `heavy_mobius_baryon_probe` | heavy-quark baryons = freest channel; heavy quark spectator ⟹ Möbius/flux gap 2√σ≈0.85 GeV FLAVOR-INDEPENDENT (same c,b); predictions Λ_c~3.14, Ω_c~3.54, Λ_b~6.47, Ω_b~6.89, Ξ_cc~4.47 GeV — all above current data (findable, not excluded), above orbital tower; Ξ_cc/Ω_b entirely unexplored; cross-flavor correlation = signature (no exotic J^P); exact mass/J^P open | HEAVY_MOBIUS_BARYON_FLAVOR_INDEPENDENT_GAP_FINDABLE_UNCONSTRAINED |
| `heavy_mobius_baryon_decay_probe` | completes #103 with decays + search: decay = twist-unwinding (non-orientable −1 → orientable +1 ground sheds 2√σ as light isoscalar hadrons) ⟹ inherits flux-tube HYBRID SELECTION RULE — single-S-wave-π-to-ground SUPPRESSED, Σ_Q π / isoscalar dipion Λ_Q(ππ) / P-wave+π PREFERRED (the branching PATTERN distinguishes Möbius from a radial excitation); open channels Λ_Q ππ Q=569, Σ_Q π 542/515, Σ_Q* π 477/496, Λ_Q η 301, DN/BN 332/251 MeV; CROSS-FLAVOR Q-MATCH all-light Q identical c=b (Λ_Q ππ 569, Λ_Q η 301; Σ_Q π offset only by hyperfine 167/194); broad (~tens–150 MeV) ⟹ amplitude analyses at LHCb (Λ_Q ππ, Σ_Q π, DN/BN; Ξ_cc/Ω_b wide open) / Belle II; branching fractions/width/J^P open | HEAVY_MOBIUS_BARYON_TWIST_UNWINDING_HYBRID_SELECTION_RULE_CROSS_FLAVOR_FINDABLE |
| `nonorientable_experimental_note_probe` | compiles the PRs #100–#109 non-orientable (Möbius/closed-flux-loop) hadron sector into one compact LHCb/Belle II/BESIII-style experimental note (deliverable: docs/bam_nonorientable_experimental_note.md), every number a pushforward of √σ: single input √σ 424 / 2√σ 849 / √(4πσ) 1504 MeV; mesonic 1⁻⁺ π₁~1.62, η₁~1.85 GeV (matched); glueball 0⁺⁺ √(4πσ)~1.50 GeV (unobserved); heavy Möbius baryons Λ_c 3135, Ω_c 3544, Ξ_cc 4471, Λ_b 6469, Ω_b 6894 MeV; decays via twist-unwinding (single-π-to-ground SUPPRESSED; Σ_Q π/isoscalar dipion/P-wave+π PREFERRED) with cross-flavor Q-match (569/301); analysis handles (branching pattern, isoscalar high-m(ππ) dipion, broad⟹amplitude fits, 1⁻⁺ smoking gun); open: exact masses 0.8–1.3 GeV band, branching fractions/widths, baryon J^P | NONORIENTABLE_SECTOR_EXPERIMENTAL_NOTE_COMPILED |
| `program_synthesis_probe` | capstone: classifies all results into 5 epistemic tiers and counts the input budget — 2 dimensionful anchors (B4: m_e, √σ; the mandatory minimum) + 2 localized open dimensionless residuals (neutrino ε, quark n_part) + 1 universal flavor puzzle; the rest ~22 derived-geometry + 6 non-orientable topological predictions (matched→free) | BAM_TWO_ANCHORS_LOCALIZED_RESIDUALS_DERIVED_GEOMETRY_TOPOLOGICAL_PREDICTIONS |
| `alpha_G_ledger_classification_probe` | places α and G in the #104 ledger: G = dimensionful ANCHOR (GR-foundational scale, root of m_e/√σ via the RS tuning λ_crit=√(6\|Λ₅\|)/κ₅², PR #57); α = UNIVERSAL dimensionless RESIDUAL (used as input A_EM=αℏc/2, a=α/2π; structure derived not value; only running derived — the 137 problem; sits with flavor puzzle); ℏ = geometric (closure quantum, ℏ=m_e·R_MID·c); c = units | G_IS_ANCHOR_ALPHA_IS_UNIVERSAL_RESIDUAL |
| `scale_count_anchors_probe` | m_e and √σ NOT independent — both descend from the one bulk-gravity scale G (PR #57), so dimensionful-anchor count reduces 2→1; but the ratio √σ/m_e≈830 is UNDERIVED (no clean closure match; nearest 50π·k_5=785, 5.4% off — a near-coincidence like F_13=233), so it becomes a new open dimensionless residual; a repackaging (dimensionful→dimensionless), total irreducible inputs unchanged; cleaner "one fundamental scale G" picture | M_E_SQRT_SIGMA_NOT_INDEPENDENT_ONE_G_PLUS_UNDERIVED_RATIO |
| `ratio_832_npart_recycling_probe` | tests N_q+ΔN=832≈√σ/m_e≈830 (0.2%) as a derivation of the #106 ratio: 832=2N_q−N_lepton=4·n_part−4·k_5² is BUILT from the n_part compensator; §8-drift test decisive — 4·n_part−100 drifts 764–920 (±9%) while 830 is fixed ⟹ baseline coincidence; independent bulk shell-stress integrals O(10–70) (Σω²≈70, Σ(n+1)π≈47), never ~466/832; circular (n_part fit to the spectrum); #106 ratio stays UNDERIVED | RATIO_832_RECYCLES_N_PART_NOT_INDEPENDENT |
| `lepton_qcd_ratio_legitimate_search_probe` | the fit-independent search #107 called for: scans quantities built ONLY from fixed geometry (k_5=5, β_lepton=50π, 2π) against √σ/m_e=830.3 under 4 criteria — C1 fit-independent, C2 §8-stable, C3 <1%, C4 principled (no ad-hoc factor); C2 is AUTOMATIC for geometric candidates (no quark-ablation dependence) but C3∧C4 fail — best principled 2π·k_5³=β_lepton·k_5=785.4 (−5.4%); every sub-% match needs a reverse-engineered factor (π·265, (4/3)·k_5⁴, k_5⁵/3.77); exponential route ln(830)=6.72 vs clean 2π=6.28 (7% off); cavity integrals O(10–350) select nothing near 830; √σ/m_e stays UNDERIVED, now plausibly IRREDUCIBLE like α; BAM does NOT collapse to a single anchor (one scale G + ratio + α + flavor puzzle) | LEPTON_QCD_RATIO_NO_PRINCIPLED_GEOMETRIC_SELECTION_STAYS_UNDERIVED |

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
  - `docs/color_algebra_shell_research_plan.md` — PR #80, identifies
    the BAM-native color algebra as **SU(2) × Z₂**: SU(2) from
    B2/Hopf holonomy (PRs #59–#66, `T = iσ_y`, `T² = −I`); Z₂ from PR
    #63's inner/outer swap (C involution). SU(2) acts on the partition
    index per generation block; Z₂ swaps n=3 ↔ n=5. SU(3) NOT
    derivable from the current scaffold (all natural triplet
    candidates give SU(2)/SO(3) algebras). Findings: (i) `H_couple`
    populated with SU(2)×Z₂ generators; (ii) singlet projector built
    (1-D fully-singlet subspace = symmetric sum over 6 flavors);
    (iii) v3 species ↔ partition map revised under uniform `+ =
    heavier` reading: `(n=3, +) = d, (n=3, −) = u`, etc.; (iv) `n_part`
    re-audit: eigenvalue range factor of full Hamiltonian saturates
    at single-digit / modest-two-digit values, while observed
    inter-generation mass² hierarchy is ~6.4·10⁹ — **outside the scope
    of any BAM color algebra on the shell basis**. n_part = 233
    remains a phenomenological compensator with sharply identified
    scope. Four-PR QCD-shell arc (#77→#80) closes structurally; the
    inter-generation hierarchy remains genuinely open and most
    plausibly requires Pati-Salam SU(4) extension with a quantitative
    throat↔shell algebra map.
  - `docs/pati_salam_throat_shell_bridge_research_plan.md` — PR #82,
    builds the BAM-native throat ↔ shell `n + 3` Z₂ bridge (each
    generation has a lepton at `n = g - 1` and a quark-pair at
    `n = g + 2`; the shift = PR #68's shell-saturation threshold; no
    free parameter). Constructs the unified 12-state radial-overtone
    basis `(l=1, n=0..5, p=±)`. Mass-ratio audit under cavity-ω²
    convention: Gen 3 within 17%, Gen 1 off by factor 2.5, Gen 2 has
    **wrong sign** (BAM predicts quark heavier than lepton;
    observation has them ~equal). Identifies three open extensions
    required for full SU(4) PS: (i) BAM-native neutrinos (candidate
    channels: opposite-chirality Weyl, sterile Majorana, separate
    radial mode); (ii) 3-fold quark color (PR #80's open gap); (iii)
    **lepton-quark mass-operator unification** — v3 leptons use
    `β·k²` closure-winding (PR #71), PR #77 quarks use `ω²(l, n)`
    cavity eigenfrequency. Cavity-ω² alone cannot give the observed
    `(τ/e)² ~ 10⁷` lepton hierarchy (throat-region spread is only
    ~7.5). PR #82 sharpens the PS extension scope; does not close it.
  - `docs/throat_shell_mass_operator_unification_research_plan.md` —
    PR #83, closes extension (iii) of PR #82 at the structural-form
    level: the lepton `β·k²` (PR #71) and quark `ω²(l, n)` (PR #77)
    mass operators are the **same Bohr-Sommerfeld operator**
    `m² = (S/L_eff)²`. Unified form `m²(k, n) = (k·2π/L_throat)² +
    ((n+1)·π/L_cavity)²` with `L_throat = √(2π)/k_5`,
    `L_cavity = L_rstar`. Pillars: (1) cavity `∮√(ω²−V)dr* = (n+1)·π`
    Bohr-Sommerfeld verified to machine precision (n≥1); (2) lepton
    `β·k² = (k·2π/L_throat)²` exact; (3) `(2π/L_throat)² = k_5²·(2π)
    = 50π = β_lepton` recovered (PR #71). The two channels are PR
    #52's `N_total = N_layer1 + N_radial`; the closure quanta `2π`
    (throat full great circle) vs `π` (cavity half-cycle node) are
    BAM's pervasive full/half-cycle distinction; and `k = 0` for
    quarks is the operator-level statement of "quarks don't pass
    through the throat". Open: independent derivation of the two
    `L_eff` from one principle; the inter-generation hierarchy
    (cross-channel/mixed modes); prediction of new states.
  - `docs/winding_shell_quadrant_research_plan.md` — PR #85, maps the
    full `(k, n)` lattice of the unified operator into four quadrants
    (one sector each per generation): neutrino candidate `(0, g−1)`,
    quark `(0, g+2)`, charged lepton `(2g−1, 0)`, and the
    previously-empty **leptoquark `(2g−1, g+2)`**. The `(k≠0, n≥3)`
    quadrant flagged by PR #83 carries BOTH throat-winding (lepton) and
    cavity-resolution (quark) character — both mass terms add, so it is
    the heaviest state per generation, and it is the operator-level
    realization of the Pati-Salam `SU(4)/SU(3)` coset (quark↔lepton
    converters, PR #82). The complementary `(k=0, n<3)` quadrant is a
    candidate neutrino sector (light, non-winding) — partially closing
    PR #82's missing-neutrino extension, with the honest caveat that
    the BAM ν/charged-lepton mass ratio ~0.07 is far above observed
    `< 10⁻⁶` (needs extra suppression). Structural map only; absolute
    masses need the L_eff unification (PR #83 open) + B4 anchor.
  - `docs/neutrino_quadrant_suppression_research_plan.md` — PR #86,
    identifies the neutrino-quadrant suppression mechanism. The
    `(k=0, n<3)` quadrant has `c₁ = 0` (no winding ⟹ no Hopf charge);
    under `C` (`c₁ → −c₁`, PR #63) it is invariant, so the neutrino is
    **necessarily Majorana**. A Majorana mass admits the seesaw
    `m_ν = m_D²/M_R` with `m_D` the bare cavity-floor Dirac mass
    (~43–118 keV) and `M_R` the lepton-number-violating
    (throat↔antithroat) scale. Because `M_R ≫ m_D`, the smallness of
    `m_ν` is generic. Only the chargeless `c₁=0` sector gets the
    seesaw — charged leptons (`c₁=±1`) are Dirac and keep their full
    winding mass — explaining why only neutrinos are anomalously
    light. Required `M_R ≈ 0.3–1.8 TeV`, a new heavy input not yet
    BAM-derivable (no current BAM scale matches ~TeV). Mechanism
    BAM-native; scale open.
  - `docs/seesaw_scale_nucleation_compliance_research_plan.md` — PR #87,
    grounds PR #86's open `M_R` in the PR #58 throat-nucleation
    framework. A `ΔL=2` Majorana mass IS a throat↔antithroat (antipodal
    `Z₂`, inner/outer swap `C`) transition; PR #58's `Σc₁=0` on a
    *single* state reproduces PR #86's only-neutrino selection rule
    (`k=0` flips `0→0`, allowed; `k≠0` gives `Σc₁=∓2`, forbidden). The
    literal `M_R = `nucleation-barrier-height hypothesis is **falsified**
    (with the electron-throat `σ, ρ`, `E_c ≈ 2.8 keV`, ~10⁸ too small).
    Instead the suppression is the **tunnelling amplitude through** the
    barrier, `m_ν = m_D·e^{−S}`, so `M_R = m_D²/m_ν = m_D·e^{S}`: the
    ~TeV scale is the keV Dirac floor exponentially lifted, and the open
    input becomes a modest, generation-stable bounce action `S ≈ 15–18`
    — the instanton-rate follow-on PR #58 flagged. Mechanism + selection
    rule BAM-native; `S` (hence absolute `m_ν`) open.
  - `docs/majorana_bounce_action_research_plan.md` — PR #88, builds the
    reduced Euclidean bounce for the `ΔL=2` flip and sharpens PR #87's
    open `S`. The bounce runs along the **non-orientable tortoise path**
    (the odd extension across the throat, `c₁ → −c₁`); the tortoise
    coordinate diverges logarithmically at the throat, so a **rigid
    throat ⟹ massless neutrino** and the boundary compliance `ε` is the
    mass-generating parameter. The action is a tortoise logarithm
    `S = √(2 μ E_c)·L*(ε) ∝ ln(1/ε)` — naturally `O(10)` and coarsely
    generation-stable, the form PR #87 required. But with the EM-throat
    tension (PR #58/#87 `σ, ρ`) it **under-produces** by `~40×`
    (`S ≲ 1`); matching `S ≈ 15–18` at a sane compliance needs the
    `ΔL=2` (B−L) throat tension `~6–12×` stiffer than the EM-throat
    tension. Progressive localisation of the open input: `~TeV` mass
    (PR #86) → `O(15)` action `S` (PR #87) → `O(10)` tension ratio
    (PR #88).
  - `docs/b_minus_l_tension_ratio_research_plan.md` — PR #89, constrains
    PR #88's open tension ratio `t`. Since the `ΔL=2` flip reverses the
    throat's orientation (`c₁ → −c₁`), it is a **global** operation, so
    `t` is a global-closure enhancement of the **local** EM surface
    tension (PR #56). It is bracketed, parameter-free, by the two basic
    BAM action scales: the **closure quantum `2π`** (a single
    great-circle orientation reversal, lower bound) and the **winding
    action `k_5√(2π) = √β_lepton`** (a full throat winding, upper bound),
    so `t ∈ [2π, k_5√(2π)] ≈ [6.28, 12.53]` — exactly PR #88's required
    `6–12` (the computed `[6.41, 12.05]` sits inside). The residual is
    "where in the window" = the compliance `ε` (`t=2π ⟹ ε≈6e-7`,
    `t=√β ⟹ ε≈1.3e-2`); cross-check `m_charged/m_D ≈ 11.9 ≈ √β` lands at
    the winding edge. A constraint + identification, not a unique
    derivation (the `(t,ε)` degeneracy + bounce-normalisation caveats
    remain). Localisation: `~TeV` (PR #86) → `O(15)` `S` (#87) → `O(10)`
    `t` (#88) → the `[2π, k_5√(2π)]` window + compliance (#89).
  - `docs/boundary_compliance_bulk_geometry_research_plan.md` — PR #90,
    the capstone: derives PR #89's residual compliance `ε` from the bulk
    throat geometry. Near the neck `f ≈ 2(r−rs)/rs`, so `ε = ℓ²/(2rs)` is
    the throat's (neck-warped) **healing length**. It is sub-throat *for
    the neutrino* because the chargeless (`c₁=0`) neck has no EM term to
    prop it open (the charged `c₁=±1` neck is propped open and stays
    Dirac, PR #86) — the same chargelessness that makes the neutrino
    Majorana makes its `ε` tiny, hence its mass tiny. Natural BAM
    sub-throat scales (`R_c³`, `Δ³`, `(m_D/m_charged)²`, `E_c`) land `ε`
    inside the PR #89 window; at the **winding-edge** tension `t ≈ √β`
    (the edge PR #89's `m_charged/m_D ≈ 11.9 ≈ √β` cross-check favoured)
    the chain yields `S ≈ 15–19`, `m_ν ~ few meV` — the observed scale,
    with no input outside the throat geometry. At the `2π` edge the same
    `ε` give `S ≈ 4` (too small): the chain closes only at the winding
    edge, the same one the cross-check picked. So the whole chain (`~TeV`
    mass → `S` → `t` → window → `ε` → bulk healing length) is closed at
    order-of-magnitude — the neutrino mass *scale* is geometric, untuned;
    the precise `m_ν` and the generation spread (`×18` vs the geometric
    `×2.7`) are the residual.
  - `docs/epsilon_bulk_compliance_research_plan.md` — PR #112, sharpens the
    question "is `ε` computed from bulk compliance, or inferred from the meV
    scale?" The healing length `ℓ ~ R_c = 2σ/ρ` uses `σ,ρ` from the
    **electron** calibration (PR #58: `R_c = 2/9`), so `ε ~ R_c³ ≈ 0.011` —
    sub-throat, `O(10⁻²)` — is computed with NO neutrino input; with
    `t = k_5√(2π) = √β_lepton` the chain gives `S ≈ 16.85`, `m_ν ≈ 2.1 meV`,
    so the meV **scale** is an output (retrodiction) and the exponential
    smallness (`ε ≪ 1 ⟹ large S ⟹ tiny m_ν`) is DERIVED. But because
    `m_ν ∝ ε^{4.8}`, the `O(1)` ambiguity (`R_c³`→2, `Δ³`→20, `R_c²/2`→108
    meV) spans ×50, and the absolute compliance normalization is the
    unpinned `κ₅²/Λ₅` (only `√6` fixed, PR #57); the PRECISE `ε` stays a
    residual. So `ε` is upgraded from "inferred from the meV scale" to
    "bulk-geometric to order of magnitude" — the smallness derived, the
    exact value not.
  - `docs/generation_spread_pmns_mixing_research_plan.md` — PR #91,
    addresses PR #90's two residuals (the generation spread and the large
    PMNS mixing). Generations are the cavity radial overtones `n`, so the
    bare prediction is **normal ordering** with `m_ν ∝ m_D` (cavity-floor
    ratios `1:1.87:2.74`). The spread is widened in the right direction
    by the **overtone-dependent neck coupling**: PR #79's boundary stress
    `χ_n` decreases with `n` (0.304, 0.097, 0.039), so higher-`n`
    neutrinos are less throat-coupled ⟹ more compliant ⟹ less suppressed
    ⟹ relatively heavier (lifting `m₃` toward the `Δm²`-implied value).
    The headline: large PMNS vs small CKM is the **cross-channel**
    (leptons: charged throat-winding `k≠0` × neutrino cavity-resolving
    `k=0`) vs **intra-channel** (quarks: up & down both cavity-shell
    `k=0`) distinction — the BAM-native reason `PMNS ≫ CKM`. The spread
    direction and the mixing dichotomy are structural; the precise
    spectrum (`ε_n(χ_n)` is `O(1)`, absolute scale unmeasured) and the
    explicit angles (cross-channel overlap integrals) are open.
  - `docs/cross_channel_pmns_overlap_research_plan.md` — PR #92, computes
    the cross-channel overlap. A literal same-coordinate radial overlap
    (winding-imprint `sin(kπs)` × cavity overtones) is a near-permutation
    matrix ⟹ **small** mixing — so large PMNS is *not* a literal radial
    overlap. The real structure: the lepton generation labels live in
    DIFFERENT coordinates — charged leptons in the closure-winding
    `k=1,3,5` (Hopf fibre), neutrinos in the radial-overtone `n=0,1,2`
    (cavity) — so their map has no preferred alignment ⟹ **anarchic**
    (Haar-random) PMNS. Quantitatively the observed PMNS angles
    (33.4°, 49°, 8.6°) are typical of a Haar `U(3)` (30th/57th/4th
    percentile), while CKM (13°, 2.4°, 0.2°) is extremely atypical (joint
    `p ≈ 0`) = aligned, consistent with up & down sharing the
    radial-overtone (shell) coordinate. So PMNS ∈ anarchy class
    (cross-coordinate), CKM ∈ aligned class (intra-coordinate). The
    class-level separation is BAM-native; the specific angles are not
    pinned (anarchy is statistical; θ13 sits at the 4th percentile, the
    one mild tension).
  - `docs/theta13_residual_alignment_research_plan.md` — PR #93, resolves
    PR #92's θ13 tension. θ13 = `|U_e3|` is the corner element, linking
    the lowest winding (`k=1`) to the highest overtone (`n=2`) — the most
    coordinate-distant (two-hop) pair (gap `|g−i|=2`); θ12, θ23 are
    adjacent (gap 1). Because the throat↔shell coupling (PR #82 `+3`
    shift, PR #83 operator) is **local** in the `(k,n)` lattice, the
    `g=1↔g=3` corner is a suppressed two-hop amplitude — a residual
    nearest-neighbour alignment. A structured-anarchy model (corner
    variance `exp(−μ)`, μ=0 = PR #92 pure anarchy) with `μ≈3` shifts the
    θ13 distribution down (median 33°→~16°), makes θ13 robustly the
    smallest angle (frac 0.50→0.72), and moves observed θ13=8.6° from the
    4th to the ~21st percentile — resolving the tension — while θ12
    (~44th) and θ23 (~70th) stay typical. The mechanism robustly explains
    θ13-smallest; the exact value (μ; θ13 median saturates ~14–16°) and
    the BAM origin of the locality are open.
  - `docs/cp_majorana_phase_research_plan.md` — PR #94, the CP-phase
    sector. CP violation is **generic**: the winding amplitudes carry the
    Hopf holonomy `e^{ikχ}` (PR #60), so the cross-channel overlaps are
    intrinsically complex and `δ_CP ≠ 0, π` with probability 1 (CP
    conservation is measure-zero — no BAM symmetry forces real
    amplitudes). The Jarlskog invariant mirrors the angle dichotomy:
    `|J_PMNS| ≈ 0.026` is typical of anarchy (51st/81st percentile, large
    CP violation), `|J_CKM| ≈ 3×10⁻⁵` is extremely atypical (~0.1th, =
    aligned ⟹ CP suppressed). And the **two Majorana phases exist**
    because the neutrino is Majorana (`c₁=0`, PR #86) — CP phases of the
    ΔL=2 throat↔antithroat sector (PRs #87–#90), observable in 0νββ; a
    Dirac neutrino would have none. The specific phase values are anarchic
    (uniform), set by the Hopf phases of the overlaps and the bounce —
    not pinned (`δ_CP` is itself poorly measured, consistent with
    uniform).
  - `docs/zeronubb_effective_mass_research_plan.md` — PR #95, turns the
    arc into a falsifiable 0νββ prediction. The effective Majorana mass
    `m_ββ = |Σ U_ei² m_i|` combines: 0νββ **occurs** (neutrino Majorana ⟸
    `c₁=0`, PR #86; a Dirac neutrino would forbid it); **normal ordering**
    (PR #91) selects the NO band; **anarchic Majorana phases** (PR #94)
    populate the whole band incl. the cancellation trough (`m_ββ → ~0`);
    and the **light scale** (PR #90, ~few meV) gives `m_ββ ≲ 8 meV`. That
    is below the current bound (KamLAND-Zen 28–122 meV, so the null result
    is expected) and largely below next-gen reach (LEGEND-1000/nEXO,
    ~9–20 meV), and below the inverted-ordering floor (~19 meV). Sharp
    falsifier: a 0νββ discovery with `m_ββ ≳ 19 meV` would imply inverted
    ordering or a quasi-degenerate scale, contradicting BAM. The exact
    `m_ββ` is a band (lightest mass unmeasured + anarchic phases).
  - `docs/cosmological_sigma_mnu_research_plan.md` — PR #96, the
    cosmological companion to PR #95. The same light, normal-ordered
    spectrum fixes `Σm_ν = m1+m2+m3`: the NO floor is
    `√Δm²_21 + √Δm²_31 ≈ 58.7 meV` (the IO floor ≈ 99 meV), and the light
    scale (PR #90) gives `Σm_ν ≈ 59–65 meV` — pinned near the floor, not
    quasi-degenerate. This is consistent with Planck (<120 meV), just
    inside DESI DR1+CMB (<72 meV), and right at the DESI DR2+CMB frontier
    (~60–64 meV). Sharp falsifiers: a robust `Σm_ν < 58.7 meV` excludes
    normal ordering (and is in tension with the oscillation `Δm²`); a
    quasi-degenerate `Σm_ν ≳ 100 meV` contradicts the light scale. `Σm_ν`
    and `m_ββ` (PR #95) are one spectrum's two observables — a joint,
    cross-checkable prediction; the exact `Σm_ν` is a narrow band (the
    lightest mass is unmeasured).
  - `docs/neutrino_mev_scale_sharpening_research_plan.md` — PR #111,
    sharpens the PR #96 band into a PINNED meV-scale spectrum. Updating to
    NuFIT 6.0 fixes `m₂ = 8.65`, `m₃ = 50.34 meV` (NO floor `Σm_ν = 59.0`),
    and the 2025 DESI DR2 + CMB bound (`≲ 60–64 meV`) corners `m₁ ≲ 3 meV`
    ⟹ `Σm_ν ∈ [59.0, 62.6] meV` (tightened from 59–65, toward the floor).
    The pinned spectrum gives the laboratory effective masses: `m_β ≈
    8.8–9.3 meV` and a NONZERO 0νββ floor `m_ββ ∈ [1.5, 3.7] meV` (in NO
    the solar term `s12²c13²m₂ = 2.60` exceeds the reactor term `s13²m₃ =
    1.10 meV`, so the contributions cannot fully cancel). Honest
    reachability: only `Σm_ν` is near-term testable (DESI, at the floor
    now); `m_β` sits ~4–5× below Project 8 and `m_ββ` ~3–10× below
    LEGEND-1000 / nEXO. Flag: some 2025 DESI + CMB fits already prefer
    `Σm_ν` at/below the floor ⟹ tension for all normal-ordered models. Open:
    `m₁` within its cornered band (0–3 meV) and the anarchic Majorana
    phases (which set `m_ββ` within the floor band).
  - `docs/npart_dynamical_hierarchy_research_plan.md` — PR #97, revisits
    the `n_part = 233` quark compensator (PR #76) with the now-complete
    lepton/neutrino sectors. The neutrino arc proved a huge hierarchy can
    be geometric (the `e^{S}` tortoise bounce, ~10⁶), so *size* is not the
    obstruction. The quark inter-generation hierarchy is non-geometric for
    a diagnosable reason: it is **irregular** (up-type `c/u≈588` vs
    `t/c≈136` ⟹ not exponential; up/down asymmetric ⟹ not a single power
    law) — the signature of QCD-RG running (`α_s` logarithmic). The
    geometric shell `ω²(1,n=3,4,5)` carries only ×2.2 of the ×6.4×10⁹
    observed mass² span. So the quark sector is the program's ONE
    DYNAMICAL hierarchy; the quark closure integer is the only one that
    §8-drifts, and the lepton↔quark gap `N_q − N_lepton = 366` quanta is
    the dynamical (QCD) excess. PR #76's compensator verdict is upheld and
    sharpened — `n_part` compensates a dynamical hierarchy; the right
    route is a QCD-shell model *with* `α_s` running. Not a derivation
    (none should exist in the geometric machinery).
  - `docs/quark_hierarchy_flavor_puzzle_research_plan.md` — PR #98,
    refines #97 by taking the first step on its "right route" and testing
    the mechanism. QCD's mass anomalous dimension `γ_m` is
    flavor-universal, so quark mass *ratios* are RG-invariant — `α_s`
    running sets the overall scale, not the hierarchy. So the hierarchy is
    NOT QCD running (as #97 said) but the **flavor puzzle** (the irregular
    Yukawa couplings, free SM inputs). The quark Yukawas overflow the
    compressed shell-overtone capacity (`ω(1,n=3,4,5)` range ×1.49) by
    ~×5×10⁴, which is why `n_part` compensates — whereas the charged
    leptons (also a flavor puzzle) are fit by the winding ladder
    (`k∈{1,3,5}`, PR #71) that has the range. BAM captures the quark
    sector's STRUCTURE (counting / quantum numbers) geometrically; the
    Yukawa MAGNITUDES are the flavor puzzle, open across all physics.
    #97's core (dynamical / non-geometric, `n_part` compensates) stands;
    the mechanism is sharpened from "QCD-RG" to "flavor puzzle".
  - `docs/qcd_confinement_cornell_audit_research_plan.md` — PR #99, pivots
    from the quark *mass* sector to the QCD *confinement* sector (the
    geometric part of QCD in BAM) and audits the Cornell potential
    `V(L)=σL−A·ℏc/L`. The linear `σL` is the flux tube = a 1D
    wormhole-bridge of constant tension; the Coulomb `−A·ℏc/L` is
    short-distance throat/gluon exchange. String breaking is Schwinger
    pair nucleation `exp(−πm_q²/(σL))` — the QED Schwinger form with
    `eE→σ`, i.e. the **PR #58 throat-pair mechanism** (`e E_S R_MID =
    m_e c²`) in the QCD sector: the string snaps when `σL ≈ 2 m_q`. The
    BAM `σ` reproduces the Regge slope `α'=1/(2πσ)=0.884 GeV⁻²` (observed
    ~0.88–0.93) and the string-breaking length (`L≈1.4 fm` vs lattice
    1.35). `√σ ≈ 0.42 GeV` is the single QCD scale anchor (the B4
    analogue: lepton `m_e` ↔ QCD `√σ`); the Cornell form + Schwinger
    break + Regge slope are geometric, the `σ` value calibrated to lattice
    (not derived).
  - `docs/glueball_closed_flux_loop_research_plan.md` — PR #100, uses
    closed flux loops (glueballs) as a pure-confinement benchmark vs
    lattice QCD — the cleanest confinement probe (no valence quarks, no
    flavor puzzle). The BAM orientable closed-loop ground state
    `√(4πσ) ≈ 1.50 GeV` (3.5√σ) benchmarks the lattice 0++ `√σ` scale
    (4.1√σ) to ~13%, and the closed-string glueball Regge slope is half
    the meson. **Where BAM's topology diverges:** the machinery has both
    orientable (`make_glueball_ring`, periodic) and **non-orientable**
    (`make_mobius_tube`, antiperiodic) closed loops; the Möbius
    antiperiodic boundary condition shifts the modes integer → half-
    integer, giving an extra **Möbius glueball tower** shifted by `πσ` in
    `M²` that interleaves the orientable one (≈2× the states). Since
    glueballs are *not experimentally observed* (they mix with qq̄ mesons),
    this is a legitimate BAM-vs-lattice difference for a non-observable —
    testable against lattice (pure-glue states), not contradicted by
    experiment. The √σ scale + the topological doubling are robust; the
    exact `M/√σ` coefficients need the full closed-loop dynamics.
  - `docs/mobius_exotic_sector_research_plan.md` — PR #101, pursues the
    Möbius topology into the **open** flux-network exotics (hybrids,
    multiquark), where — unlike glueballs — the states ARE observed. The
    flux-network topology is the hadron taxonomy (meson tube, baryon
    Y-junction, tetraquark / pentaquark multi-junction, hybrid
    tube+twist, glueball loop). A non-orientable (Möbius) flux tube
    carries the antiperiodic phonon that opens the **exotic J^PC** (`1-+`,
    forbidden to ordinary qq̄ with `P=(−1)^{L+1}`, `C=(−1)^{L+S}`), and
    the observed exotic hybrids `π₁(1600)`, `η₁(1855)` (all `1-+`) match
    at the right J^PC and at `ρ/ω + 2√σ ≈ 1.62, 1.85 GeV` (the `2√σ`
    flux-tube gap). The observed tetraquarks (`X, Z_c, T_cc`) and
    pentaquarks (`P_c`) fit the multi-junction networks. The contrast
    with PR #100: glueballs unobserved (BAM free to differ), exotics
    observed (BAM's non-orientable topology must — and does — meet data).
    The Möbius twist is the same Z₂ giving the throat spin-½ (PR #63–#67);
    the Möbius baryon is a BAM-specific prediction.
  - `docs/baryonic_exotics_classification_research_plan.md` — PR #102,
    classifies the BAM-specific baryonic exotics (Möbius / hybrid baryon)
    and ranks the channels by experimental constraint. The key subtlety:
    unlike mesons (where `1-+` is a smoking-gun exotic via `C`), baryons
    have NO forbidden `J^P` (`P=(−1)^L`, `S∈{½,3/2}`, no `C` ⟹ every
    half-integer `J^P` ordinary), so BAM's Möbius/hybrid baryons are
    *supernumerary ordinary-`J^P`* states — identifiable only by counting.
    They land in the light N*/Δ* region (`nucleon/Δ + 2√σ ≈ 1.79, 2.08
    GeV`), the densest, best-measured part of the spectrum — the MOST
    experimentally constrained corner of BAM's non-orientable predictions
    (opposite extreme from the unobserved glueballs). The Möbius doubling
    must either coincide with observed-but-unexplained resonances (fill
    missing resonances) or decouple from `πN` (the standard
    missing-resonance mechanism), else be excluded. Constraint ranking:
    light N*/Δ* > strange hyperons > charm/bottom baryons (the freest).
  - `docs/heavy_mobius_baryon_research_plan.md` — PR #103, the concrete
    prediction in that freest channel. Heavy-quark symmetry (the heavy
    quark is a spectator) makes the Möbius/flux excitation gap
    `Δ = 2√σ ≈ 0.85 GeV` FLAVOR-INDEPENDENT — the same above the charm and
    bottom ground baryons — which is the cross-flavor signature replacing
    the absent exotic-`J^P` smoking gun. Predictions: Λ_c ~3.14, Ω_c
    ~3.54, Λ_b ~6.47, Ω_b ~6.89, Ξ_cc ~4.47 GeV — all just above the
    currently-measured excitation ceilings (findable at LHCb/Belle II, not
    excluded) and above the orbital tower (a supernumerary state). The
    doubly-heavy Ξ_cc and Ω_b have no measured excitations at all —
    entirely unconstrained. A correlated counting prediction; exact mass
    (lattice hybrid gap 0.8–1.3 GeV) and `J^P` open.
  - `docs/heavy_mobius_baryon_decay_research_plan.md` — PR #109, the decay
    channels + search strategy completing PR #103. The Möbius excitation is
    the non-orientable (orientation −1) flux sector and the ground heavy
    baryon is orientable (+1), so the decay proceeds by UNWINDING the
    half-twist, shedding `2√σ ≈ 0.85 GeV` as light isoscalar hadrons (a
    hybrid de-excitation; heavy quark a spectator). This inherits the
    flux-tube HYBRID SELECTION RULE — single-S-wave-π-to-ground SUPPRESSED;
    `Σ_Q π` / isoscalar S-wave dipion `Λ_Q(ππ)` / P-wave-baryon+π PREFERRED
    — the branching PATTERN that distinguishes the Möbius baryon from an
    ordinary radial excitation (which does the opposite). Because the gap
    `2√σ` and the light-meson thresholds are both flavor-independent, the
    all-light Q-values are CROSS-FLAVOR IDENTICAL (`Λ_Q ππ` 569, `Λ_Q η`
    301 MeV for both c and b; `Σ_Q π` offset only by the `Σ_Q − Λ_Q`
    hyperfine splitting 167/194 MeV). Honest: with several open channels at
    `Q ≈ 0.5 GeV` the state is broad (~tens–150 MeV), best resolved in LHCb
    / Belle II amplitude analyses of `Λ_Q ππ`, `Σ_Q π`, `DN`/`BN` (`Ξ_cc`,
    `Ω_b` wide open). Absolute branching fractions, total width, and `J^P`
    remain open — the predictions are the branching pattern and the
    Q-structure, not partial rates.
  - `docs/nonorientable_experimental_note_research_plan.md` — PR #110, the
    compact experimental note compiling the PRs #100–#109 non-orientable
    sector for an LHCb / Belle II / BESIII reader (deliverable:
    `docs/bam_nonorientable_experimental_note.md`). Every number is a
    pushforward of the single input `√σ`: the mesonic `1⁻⁺` hybrids (π₁
    ~1.62, η₁ ~1.85 GeV, matched), the `0⁺⁺` glueball (`√(4πσ)` ~1.50 GeV,
    unobserved), the heavy Möbius baryon masses (Λ_c 3135 … Ω_b 6894 MeV)
    and their decays (twist-unwinding ⟹ single-π-to-ground suppressed, `Σ_Q
    π` / isoscalar dipion / P-wave+π preferred; cross-flavor Q-match 569 /
    301 MeV). Analysis handles: the branching pattern (vs a radial
    excitation), the isoscalar high-`m(ππ)` dipion, broad widths ⟹
    amplitude fits, and the `1⁻⁺` smoking gun in the mesonic sector. A
    reference card — established masses + decay pattern; exact masses
    (0.8–1.3 GeV band), branching fractions / widths, and baryon `J^P` open.
  - `docs/program_synthesis_research_plan.md` — PR #104, the capstone
    synthesis. Classifies every major result into five epistemic tiers and
    counts the input budget: BAM's entire DIMENSIONFUL content reduces to
    **two B4 anchors** — `m_e = ℏc/R_MID` (QED/lepton) and `√σ ≈ Λ_QCD`
    (confinement) — the irreducible B4 minimum (one scale per sector,
    PR #52). The genuinely-open DIMENSIONLESS inputs are localized to two —
    the neutrino compliance `ε` (bracketed `[2π, k_5√(2π)]`) and the quark
    `n_part` (a flavor-puzzle compensator). Beyond these there is one
    UNIVERSAL open problem, the flavor puzzle (the quark Yukawa hierarchy,
    derivable by no theory — not BAM-specific). Everything else is derived
    geometry (~22 results) or a non-orientable topological prediction (~6,
    spanning matched → falsifiable → constrained → findable → free). In one
    line: two mandatory anchors + a couple of localized residuals + the
    universal flavor puzzle, the rest geometry and falsifiable predictions.
  - `docs/alpha_G_ledger_classification_research_plan.md` — PR #105,
    places the fundamental constants in the #104 ledger. **G** is the
    dimensionful ANCHOR — the GR-foundational scale (the throat is a
    gravitational wormhole; its size, the one B4 length, is set by the
    bulk gravity via `λ_crit = √(6|Λ₅|)/κ₅²`, PR #57) and the root the
    #104 sector anchors (`m_e`, `√σ`) descend from. **α** is a UNIVERSAL
    dimensionless RESIDUAL — used as a numerical input throughout
    (`A_EM = α·ℏc/2`, `a = α/2π`); BAM derives the charge unit (`|c₁|=1`),
    the `1/2π` measure, and α's running, but the value 1/137 is a free
    input as in the SM (the "137 problem"), so it sits with the flavor
    puzzle, not the BAM-specific residuals. **ℏ** is geometric (the
    closure quantum, `ℏ = m_e·R_MID·c`); **c** is units. Refines #104: α
    was a silent residual input to its "derived geometry" tier, and G is
    the root of its two sector anchors.
  - `docs/scale_count_anchors_research_plan.md` — PR #106, settles the
    scale-count question #105 raised. `m_e` and `√σ` are NOT independent —
    both are brane scales of the one bulk geometry, descending from the
    bulk gravity `G` (PR #57: `m_e=ℏc/R_MID` with `R_MID` from `λ_crit=
    √(6|Λ₅|)/κ₅²`; `σ∝√|Λ₅|/κ₅²`) — so the dimensionful-anchor count
    reduces **2→1** (`G`). But their dimensionless ratio `√σ/m_e≈830` (the
    lepton-throat/QCD-confinement hierarchy) is **underived** — no clean
    closure number (nearest `50π·k_5=785`, 5.4% off, a near-coincidence
    like `F_13=233`). So it is a **repackaging, not a free reduction**: a
    dimensionful anchor becomes a dimensionless residual (joining `ε`,
    `n_part`, `α`), total irreducible inputs unchanged. The gain is the
    GR-foundational cleanliness — the sole fundamental *scale* is `G`,
    everything else dimensionless. Deriving the ~830 ratio (the channel
    normalisation) would reduce BAM to a single irreducible input.
  - `docs/ratio_832_npart_recycling_research_plan.md` — PR #107, tests the
    tempting `N_q + ΔN = 832 ≈ √σ/m_e ≈ 830` (0.2%) as a derivation of the
    #106 ratio, and rejects it: `832 = 2N_q − N_lepton = 4·n_part − 4·k_5²`
    is built from the `n_part` compensator. The decisive §8-drift test —
    propagating `n_part ∈ {216..255}` through `4·n_part − 100` gives
    `[764, 920]` (±9%) while the observed 830 is fixed — shows it is a
    baseline coincidence (like `50π·k_5=785`, `F_13=233`), not a stable
    selection. No independent bulk shell-stress integral yields ~466/832
    (the natural ones, `Σω²≈70`, `Σ(n+1)π≈47`, are `O(10–70)`); 466 enters
    only via the v3-fit closure count `4β_quark/(2π)=2·n_part`. It is
    circular (n_part was fit to the spectrum). So `√σ/m_e` stays
    UNDERIVED; the PR #106 ledger is unchanged.
  - `docs/lepton_qcd_ratio_legitimate_search_research_plan.md` — PR #108,
    the fit-independent search #107 called for. Scans quantities built ONLY
    from fixed geometry (`k_5=5`, `β_lepton=50π`, `2π`) against
    `√σ/m_e=830.3` under four criteria (C1 fit-independent, C2 §8-stable,
    C3 <1%, C4 principled). C2 is automatic for geometric candidates (they
    never touch the quark ablations that made `n_part` drift), so the
    binding bars are C3 ∧ C4 — and no candidate clears both: the best
    *principled* candidate `2π·k_5³ = β_lepton·k_5 = 785.4` is `−5.4%` off,
    and every sub-% match needs a reverse-engineered factor (`π·265`,
    `(4/3)·k_5⁴`, `k_5⁵/3.77`). The exponential route fails too
    (`ln(830)=6.72` vs the clean action `2π=6.28`, 7% off), and the
    Tangherlini cavity integrals (`O(10–350)`) select nothing near 830. So
    `√σ/m_e` stays UNDERIVED and is now plausibly IRREDUCIBLE, like `α`;
    BAM does NOT collapse to a single anchor (one scale `G` + this ratio +
    `α` + the flavor puzzle).
  - `docs/odd_k_closure_lemma.md` — the closure arithmetic this upgrades.
  - `docs/hbar_origin_status.md` — B4 (the m_e anchor).
  - `docs/tree_qed_status.md` — the tree-QED result the F² target
    summarises.
