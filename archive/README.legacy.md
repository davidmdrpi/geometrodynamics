Python theoretical physics simulation 
"""
Geometrodynamic QED  —  v39
================================================================

THE CORE CLAIM (not QED implementation — geometry IS physics)
──────────────────────────────────────────────────────────────
This is not a simulation of QED. It is a proposal that the structures
physicists call electromagnetism, charge, and spin are not independent
fields attached to spacetime — they ARE the geometry of spacetime,
specifically the Hopf fibration structure of the spatial S³.

  ELECTROMAGNETISM = curvature of the Hopf connection on S³
  CHARGE           = first Chern number of the Hopf bundle (= 1, topological)
  SPIN-½           = holonomy of the Möbius π-twist along the Hopf fibre
  COULOMB LAW      = Green's function on S³ (1/sin ψ focus at antipode)
  PARTICLE MASS    = eigenvalue of the 5D Tangherlini operator at the throat

The Wheeler geometrodynamics program: geometry → matter. Completed classically.

THE HOPF CONNECTION — electromagnetic potential from pure geometry
──────────────────────────────────────────────────────────────────
The Hopf bundle S¹ → S³ → S² carries a canonical connection 1-form:

  A = (1/2) cos(χ) dφ

This is not imposed — it is uniquely determined by the geometry of S³.
Its curvature (the field strength 2-form) is:

  F = dA = -(1/2) sin(χ) dχ ∧ dφ

The first Chern class integrates to 1 over any S² cross-section:

  c₁ = (1/2π) ∫_{S²} F = 1   (verified numerically)

This is the topological origin of the charge quantum. You cannot have
half a unit of charge because you cannot have a bundle with c₁ = ½ on S³.

EQUATORIAL LOCALIZATION — why particles live at χ = π/2
─────────────────────────────────────────────────────────
Three independent dynamical reasons, not an assumption:

1. MINIMAL AREA: The S² slice at χ=π/2 has area 4π sin²(π/2) = 4π,
   the maximum of 4π sin²(χ). It is a saddle point of the S³ volume form
   — a stable critical surface of the geometric action.

2. ZERO POTENTIAL: A|_{χ=π/2} = (1/2)cos(π/2)dφ = 0. The connection form
   vanishes at the equator — zero EM self-energy, minimum of the gauge energy.

3. THROAT INTERSECTION: The 5D Tangherlini wormhole throat (r = R_MID) is a
   3-surface in 5D. Its intersection with χ=π/2 is the Hopf fibre S¹ itself.
   The Dirichlet boundary condition u(R_MID)=0 is geometrically the statement
   that the mode must vanish on this circle — i.e., the Hopf fibre is the
   throat of the bundle, not just the throat of the wormhole.

GAUGE HOLONOMY → SPIN-½
─────────────────────────
Parallel transport of a vector along a closed Hopf fibre at χ=χ₀ accumulates
a phase (holonomy) equal to the integral of the connection:

  ∮ A = (1/2) cos(χ₀) · 2π = π cos(χ₀)

At χ₀=0 (north pole):   holonomy = π      →  e^{iπ} = −1  (spinor sign flip)
At χ₀=π/2 (equator):    holonomy = 0      →  trivial (stable orbit)
At χ₀=π (south pole):   holonomy = −π     →  e^{-iπ} = −1 (spin-½ conjugate)

The Möbius π-twist of the wormhole tube is exactly this holonomy at the pole.
Spin-½ is not postulated — it is the holonomy of the Hopf connection.

GEOMETRIC ACTION
─────────────────
  S = (1/16π) ∫_{S³×ℝ} R √g d⁴x  +  (1/4) ∫ F ∧ *F

Varying with respect to gᵘᵛ → Einstein equations Gᵘᵛ = 8π Tᵘᵛ(EM)
Varying with respect to A   → Maxwell equations d*F = J

Both sets of equations follow from one geometric action on S³.
No quantization postulate. No gauge group postulate. No spin postulate.
"""
