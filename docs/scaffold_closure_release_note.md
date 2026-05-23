# BAM effective-action scaffold — closure release note (through PR #53)

A milestone summary of the BAM effective-action scaffold programme,
which began with five mismatch terms (B1–B5) blocking a single covariant
5D action and is now **complete**: four barriers are derived from
geometry, and the fifth (B4, the single dimensionful anchor) is shown to
be irreducible by dimensional necessity — and relocatable to the
invariant bulk separation ΔR.

This note collects the closure arc (PRs #49–#53). The living tracker is
`docs/bam_scaffold_status.md`; this is the point-in-time release summary.

## Barrier ledger — final

| barrier | type | original status | final status | by (PR) |
|---|---|---|---|---|
| **B1** closure quantum `∮A = 2πn` | topological | imposed constraint | **CLOSED** — winding θ-term `S_top = 2π·n` | #49 |
| **B2** antipodal `Z₂` (`T = iσ_y`) | discrete | imposed identification | **CLOSED** — `RP³ = S³/Z₂` + non-trivial spin structure | #49 |
| **B3** hard-wall throat BC | boundary | imposed by hand | **CLOSED** — single-valuedness under `T² = −I` ⟹ `ψ(throat) = 0` | (hard-wall derivation) |
| **B5** 5D→4D reduction → F² | reduction | unconstructed | **CLOSED** — one `C × S³` master functional yields masses **and** `F²=K²·Q` | #50, #51 |
| **B4** dimensional bridge `ℏ = m_e·R_MID·c` | scale | one external anchor | **IRREDUCIBLE** — scale-free machinery needs exactly one anchor; relocatable to ΔR | #52, #53 |

## The closure arc (PRs #49–#53)

  - **#49 — topological/discrete sector** (`topological_discrete_sector_probe`).
    B1 + B2 promoted to action data: `RP³ + spin structure + winding
    θ-term`. The spin structure (`T² = −I`) also forces the B3 hard wall.

  - **#50 — radial reduction bridge** (`radial_reduction_bridge_probe`).
    The 5D→4D reduction factorizes into three channels (radial → masses,
    S³ → gauge+propagator, throat → F²). Central honest finding: `F²` is
    the throat-channel form factor, **not** a radial overlap. B5 reduced
    to the residual "one master integral producing masses and F²
    together."

  - **#51 — bulk-boundary interaction + master integral**
    (`bulk_boundary_interaction_probe`, `master_integral_probe`).
    The same throat cavity yields the mass spectrum (Green-function
    poles) and `K` (throat dwell-time impedance); then the S³ Hopf
    Q-channel is integrated, completing the master functional

    ```
    ℳ(ω; x, c) = G_C(r, r′; ω) ⊗ 𝒢_{S³}(Ω, Ω′)
    ```

    on the warped product `C × S³`: ω-poles → masses, throat boundary →
    `K(x)=2x/(1+x)`, S³ Hopf → `Q(x,c)`; vertex residue → `F²=K²·Q` to
    `2e-14`. **Masses and the full vertex from one object.** The
    factorization is the product-geometry consequence (separation of
    variables), not a failure to unify. **B5 closed.**

  - **#52 — B4 audit + Maslov closure-ledger**
    (`maslov_dimensional_bridge_probe`). The closure-ledger Bohr–
    Sommerfeld integer `n+1` per radial mode is the **Maslov index** of
    the doubly-Dirichlet throat cavity (`μ=4`); its throat half (`μ=2`,
    reflection phase `π`) is the B3 wall = closure half-quantum = dwell
    phase. The machinery is **scale-free** (rescaling `R_MID → λ·R_MID`
    leaves every dimensionless output invariant to machine precision), so
    a dimensionful scale cannot emerge — exactly **one** external anchor
    is mathematically required. **B4 irreducible by dimensional
    necessity** (a structural feature, cf. SI fixing `c, ℏ, e`).

  - **#53 — ΔR scale-modulus** (`delta_r_scale_modulus_probe`). The one
    required anchor need not be the particle mass `m_e` — it can be the
    invariant bulk separation `ΔR = R_OUTER − R_INNER = 0.52·R_MID`,
    *provided* ΔR is a proper (cosmologically fixed) length. It is: the
    throat is a static bound vacuole (discrete spectrum + vacuum +
    dimensionless BC) decoupled from Hubble flow (`ΔR/R_cosmo ~ 10⁻³⁹`),
    and a comoving throat is observationally excluded (it would redshift
    masses as `(1+z)` against quasar bounds `≲10⁻⁵`). The bridge becomes
    `m_e = f_closure·ℏ/(ΔR·c)` with `f_closure = 0.52` — `m_e` a
    consequence of a fixed bulk length, with local mass ratios predicted
    constant in cosmic time.

## Headline

The scaffold is **complete** in a precise sense:

  - **Four barriers (B1, B2, B3, B5) are derived** from the geometry —
    the topological/discrete sector, the spin-structure hard wall, and
    the `C × S³` master integral.
  - **The fifth (B4) is the single mandatory dimensionful unit.** It is
    not an unsolved gap: a scale-free theory provably requires exactly
    one external dimensionful anchor. That anchor is identified as a
    cosmologically-invariant geometric length (ΔR), with `m_e` a
    consequence.

## What remains genuinely open

  - **The value of the one anchor.** ΔR (equivalently `R_MID`, `m_e`) is
    still one external dimensionful number — relocated to a geometric
    invariant, not derived. Pinning it to a second fixed scale (e.g. a
    closure-quantum relation to the Planck length) would be a genuine
    derivation.
  - **A first-principles covariant action.** The channel reductions
    (radial Sturm–Liouville, throat dwell-time, Hopf helicity, Maslov
    indices) are read from the WKB/boundary structure; writing them as
    one explicit 5D Lagrangian density with the throat boundary term is
    the standing follow-on.
  - **Loop corrections.** The amplitude thread is tree-level throughout.

## Cross-references

  - `docs/bam_scaffold_status.md` — the living barrier tracker.
  - `docs/topological_discrete_sector_research_plan.md` — B1 + B2 (#49).
  - `docs/radial_reduction_bridge_research_plan.md` — B5 factorized (#50).
  - `docs/bulk_boundary_interaction_research_plan.md` — B5′ radial+throat (#51).
  - `docs/master_integral_research_plan.md` — B5 closed (#51).
  - `docs/maslov_dimensional_bridge_research_plan.md` — B4 audit (#52).
  - `docs/delta_r_scale_modulus_research_plan.md` — B4 anchor as ΔR (#53).
  - `docs/hbar_origin_status.md` — the dimensionless residuals closed.
  - `docs/tree_qed_status.md` — the tree-QED result the F² target summarises.
