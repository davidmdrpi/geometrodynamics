# Spin-½ Wigner rotation falsifier probe

The companion to the moving-throat probe (PR #59), which showed a boosted
throat carries the correct relativistic **energy–momentum**
(`E²−(pc)²=(mc²)²`, invariant mass = static eigenvalue). This probe tests
the other half of "throat = relativistic particle": does the throat's
**spin** transform correctly under boosts — i.e. does the Hopf-holonomy
spin (the Berry phase `∮A = π cos χ`) reproduce the relativistic
**Wigner rotation**? A genuine falsifier: if the geometric-phase spin
does not match the Wigner SU(2) holonomy (wrong spin value, wrong double
cover, wrong solid-angle law), BAM fails.

## The structures

  - **Hopf holonomy (BAM spin).** The Hopf connection `A_φ = ½ cos χ`
    (`geometrodynamics/hopf/connection.py`) is the spin-½ monopole
    connection; its holonomy around a fibre is `∮A = π cos χ`. For a cap
    at polar angle `χ` the solid angle is `Ω = 2π(1 − cos χ)`, so the
    holonomy is `−½Ω + π` — **the spin-½ Berry phase = ½ × solid angle**.

  - **Wigner rotation (relativistic spin).** Two non-collinear Lorentz
    boosts compose into a boost *and* a rotation: in `SL(2,C)`,
    `B₂B₁ = U·P` (polar decomposition), with the unitary part
    `U ∈ SU(2)` the Wigner rotation. The boost is the spinor
    `B(ζ,n̂) = cosh(ζ/2) + sinh(ζ/2) n̂·σ` — the spin-½ representation
    (the `½` in `ζ/2`).

## The claim (the falsifier)

Both are the **same spin-½ SU(2) holonomy**:

  - both carry the **spin-½ factor ½** (the Hopf monopole charge `½`; the
    spinor `ζ/2` in the boost);
  - both live in the **spinor double cover** (`2π → −1`, `4π → +1`; the
    Hopf/RP³ double cover, B2, `T² = −I`);
  - both obey the **geometric-phase law** "rotation = ½ × enclosed solid
    angle" for spin-½ — the Wigner rotation for a closed boost loop and
    the Hopf Berry phase for a closed Bloch-sphere loop.

If the throat's Hopf spin reproduces the Wigner rotation, BAM's throat is
a genuine relativistic spin-½ particle. If not — if `A_φ = c cos χ` with
`c ≠ ½` (wrong spin), or no double cover (wrong statistics) — BAM fails.

## B4 accounting

The spin / Wigner structure is purely **geometric and dimensionless**
(angles, `SU(2)` elements, solid angles) — no scale dependence. Spin-½ is
a topological/geometric property, independent of the single dimensionful
anchor `m_e` (PRs #55–#58). Consistent with the whole arc: structure
derived, scale is the one anchor.

## Tests

  T1. **Hopf holonomy = spin-½ Berry phase.** `∮A = π cos χ = −½Ω + π`
      with `Ω = 2π(1−cos χ)` — the spin-½ monopole (½ × solid angle).
  T2. **Wigner rotation from `SL(2,C)`.** `B₂B₁ = U·P`; the Wigner angle
      from `U`'s `SU(2)` trace matches the closed-form
      `tan(ω/2)` formula; collinear boosts → no rotation.
  T3. **Spinor double cover.** The Wigner rotation is `SU(2)`:
      `R(2π) = −I`, `R(4π) = +I` — the same double cover as Hopf/RP³ (B2).
  T4. **Thomas precession (infinitesimal Wigner).** Small-rapidity
      `ω ≈ ½ ζ₁ζ₂ sin θ` (the ½ = spin-½); the Thomas factor
      `γ²/(γ+1)`.
  T5. **Common spin-½ holonomy.** Both the Hopf Berry phase and the
      Wigner rotation are the spin-½ `SU(2)` holonomy with the
      characteristic `½` and the geometric-phase law.
  T6. **Falsification criterion.** `c = ½` → spin-½ (not `c ≠ ½`); double
      cover present → fermion; BAM passes.
  T7. **B4 accounting.** Spin/Wigner structure dimensionless/geometric;
      scale-independent; spin-½ independent of the anchor.
  T8. **Assessment.**

## Verdict structure

  - **SPIN_WIGNER_COVARIANT** (expected): the throat's Hopf-holonomy spin
    reproduces the relativistic Wigner rotation — both are the spin-½
    `SU(2)` holonomy with the factor `½`, the spinor double cover
    (`2π → −1`), and the geometric-phase law "rotation = ½ × solid
    angle." Combined with PR #59 (energy–momentum), the boosted throat is
    a genuine relativistic spin-½ particle. BAM survives the falsifier.

  - **WIGNER_FALSIFIED**: the Hopf spin gives the wrong spin value, lacks
    the double cover, or does not match the Wigner `SU(2)` holonomy — the
    throat is not a relativistic spin-½ particle.

## What this leaves open

  - **The throat spinor from the full action.** The Wigner rotation here
    is the generic spin-½ `SL(2,C)` result; constructing the explicit
    boosted throat *spinor* of `S_BAM` and confirming its transformation
    is the follow-on.
  - **g − 2.** The geometric gyromagnetic ratio (`g = 2` from the Hopf
    monopole) and its loop corrections.
  - **Exact Wigner ↔ hyperbolic-area match.** The boost holonomy lives on
    the velocity hyperboloid; relating its hyperbolic area to the Bloch
    solid angle beyond the leading `½` factor.

## Cross-references

  - `docs/stable_moving_throat_research_plan.md` — energy–momentum (#59).
  - `geometrodynamics/hopf/connection.py` — `hopf_connection`,
    `hopf_holonomy` (`A_φ = ½ cos χ`, `∮A = π cos χ`).
  - `docs/topological_discrete_sector_research_plan.md` — the RP³ double
    cover / `T² = −I` (B2).
  - `docs/maslov_dimensional_bridge_research_plan.md` — the B4 theorem.
  - `experiments/closure_ledger/spin_wigner_rotation_probe.py` — this
    probe.
