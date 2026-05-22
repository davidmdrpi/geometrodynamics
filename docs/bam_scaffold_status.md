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
| **B4** dimensional bridge `ℏ = m_e·R_MID·c` | scale | one external anchor | **OPEN** | (the single `m_e` anchor; documented in `docs/hbar_origin_status.md`) |

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

### B4 → dimensional bridge (OPEN)

The action fixes dimensionless structure (`F²`, `c₁`, the closure
quantum, `K`, `Q` are all scale-free); the absolute MeV scale enters
only through the single anchor `ℏ = m_e·R_MID·c`. The 2π winding (B1)
is metric-independent, but the physical great-circle length `2πR` and
the mass scale still need `m_e`. This is the last external input; it is
the subject of the separate ℏ-origin thread
(`docs/hbar_origin_status.md`), which reduced the lepton surrogate to
this one anchor.

## Summary

The scaffold began with five barriers. Four (B1, B2, B3, B5) are now
**closed**. B1, B2, B3 were promoted to a topological/discrete action
sector (`RP³ + spin structure + winding θ-term`) whose spin structure
also forces the hard-wall BC. B5 is closed by the master integral: a
single separable functional `ℳ = G_C ⊗ 𝒢_{S³}` on the warped product
`C × S³` yields the masses (poles), `K` (throat boundary), and `Q` (S³
Hopf), with vertex residue `F²=K²·Q` — masses and the vertex from one
object; the factorization is the product-geometry consequence. One
(B4) remains **open** — the single `m_e` dimensional anchor.

The sole survivor is the most fundamental: the external mass scale
(B4), the subject of the separate ℏ-origin thread.

## Probe ledger

| probe | barriers | verdict |
|---|---|---|
| `bam_effective_action_scaffold_probe` | B1–B5 (map) | SCAFFOLD_WITH_BARRIERS |
| `topological_discrete_sector_probe` | B1, B2 | PROMOTION_SUCCEEDS |
| `hard_wall_boundary_derivation_probe` | B3 | HARD_WALL_DERIVED |
| `radial_reduction_bridge_probe` | B5 | BRIDGE_FACTORIZED |
| `bulk_boundary_interaction_probe` | B5′ (radial+throat) | BULK_BOUNDARY_FORMULATED |
| `master_integral_probe` | B5′ (+ S³ Q) | MASTER_INTEGRAL_COMPLETE |

## Cross-references

  - `docs/bam_effective_action_scaffold_research_plan.md` — the original
    5-barrier map.
  - `docs/topological_discrete_sector_research_plan.md` — B1 + B2.
  - `docs/hard_wall_boundary_derivation_research_plan.md` — B3.
  - `docs/radial_reduction_bridge_research_plan.md` — B5.
  - `docs/bulk_boundary_interaction_research_plan.md` — B5′ (radial+throat).
  - `docs/master_integral_research_plan.md` — B5′ closed (+ S³ Q).
  - `docs/hbar_origin_status.md` — B4 (the m_e anchor).
  - `docs/tree_qed_status.md` — the tree-QED result the F² target
    summarises.
