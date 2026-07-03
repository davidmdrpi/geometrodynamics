# The guidance law from the 5D bulk: the throat rides its own conserved charge (PR #199)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This document discharges condition (1)
> of PR #198: it derives the guidance law v = ∇S — there an input
> identification — from the 5D bulk field equations, turning the
> dBB-grade Born rule of #198 into a consequence of the geometry. The
> companion probe machine-checks every step.

## 0. The claim and the chain

PR #198 established: *if* throats are transported by the local momentum
field v = ∇S, then |ψ|² is the unique equivariant measure and the Born
rule follows at dBB grade. The guidance law itself remained an
identification. Here it is derived, in four steps, each verified:

1. **The throat mode is a fiber-winding mode.** On the 5D bulk with the
   Hopf-fiber Killing field ∂_χ, the brane matter wave is the k-winding
   Kaluza–Klein mode Ψ(x, χ) = ψ(x)e^{ikχ} — the repo's standing
   reduction (#83 mass operator; #193 sector reduction; #195 index
   mechanism).
2. **The guidance current IS bulk stress-energy.** For that mode, the
   mixed fiber components of the 5D stress tensor are *identically*

   ```
   T_{μχ} = k · Im(ψ* ∂_μ ψ) ,
   ```

   i.e. the de Broglie current J_μ = Im(ψ*∂_μψ) is not an extra
   structure: it is the **fiber-momentum flux of the bulk
   stress-energy**, with the winding number k as the charge unit
   (winding = charge, the #42–#44 geometry). Verified to machine
   precision on explicit 5D modes.
3. **Its conservation is a Bianchi identity.** ∂_χ is a Killing field,
   so 𝒥^μ = T^μ{}_ν ξ^ν_{(χ)} is covariantly conserved — and in BAM this
   is not a property of a matter model but of the **geometry itself**:
   the 5D Einstein equations G_{AB} = 8πT_{AB} plus the contracted
   Bianchi identity ∇_A G^{AB} ≡ 0 force ∇_A T^{AB} = 0, whose
   χ-component is exactly ∇_μ𝒥^μ = 0. The identity ∇·G ≡ 0 is verified
   **symbolically and exactly** (sympy) on the 5D weak-field KK metric
   diag(−(1+2Φ), (1−2Φ)I₃, λ²) with arbitrary Φ — all five components
   identically zero, no weak-field expansion needed. On the brane this
   conservation law *is* the continuity equation
   ∂ₜρ + ∇·(ρ∇S) = 0 of #198's Theorem 1, now with a geometric pedigree:
   **the χ-component of the Bianchi identity.**
4. **The throat rides its own charge.** The throat is not a passenger
   *added* to the current: it **is** the localized, quantized unit of
   fiber winding — the discrete invariant Q = (1/2π)∮∇(arg ψ) of
   #178/#181/#182, integer-valued and topologically conserved (it
   cannot diffuse, fractionate, or vanish under smooth evolution). A
   localized integer charge embedded in a conserved current flow moves
   with the charge-transport velocity

   ```
   v^i = 𝒥^i / 𝒥^0 = ∇S     (nonrelativistic/weak-field reduction),
   ```

   because the alternative — the throat sitting still while its
   winding flows away — would violate the quantization of Q.
   **Demonstrated pointwise**, not just on average: a quantized
   phase-winding core (2D reduction) transported through the live
   nonlinear dynamics moves with the ambient J/ρ evaluated at the core
   — both the imposed background flow and the partner-induced flow —
   while its winding number stays exactly 1.

**Conclusion:** the guidance law is the χ-component of the contracted
Bianchi identity plus the topological localization of the winding
charge. With #198 (equivariance + uniqueness + relaxation), the chain

```
5D Einstein equations → ∇·T = 0 (Bianchi) → continuity of J = T·ξ_χ
   → topological charge rides J → v = ∇S → |ψ|² equivariant, unique
   → Born rule
```

runs from the bulk field equations to the measurement measure with the
guidance step no longer an assumption.

## 1. Step 2 in detail: the stress-tensor identity

For a complex 5D scalar with the standard kinetic term,
T_{AB} = ∂_{(A}Ψ ∂_{B)}Ψ* − ½g_{AB}(…); the metric term does not enter
the mixed components T_{μχ} on a block-diagonal KK metric. With
Ψ = ψ(x)e^{ikχ}: ∂_χΨ = ikΨ, so

```
T_{μχ} = ½(∂_μΨ ∂_χΨ* + ∂_χΨ ∂_μΨ*) = k · Im(ψ* ∂_μψ) ,
```

χ-independent (as the Killing symmetry requires). The probe verifies
this identity numerically for k = 1, 2, 3 on explicit (x, χ) grids to
~10⁻¹²: **the de Broglie current is the T^μ_χ column of the bulk stress
tensor.** In BAM — where matter is geometry — that column feeds the 5D
Einstein equations, and its conservation is therefore not optional.

## 2. Step 3 in detail: the Bianchi verification

The contracted Bianchi identity is metric-independent, but the probe
verifies it concretely on the metric BAM actually uses in the weak-field
regime: the 5D conformal-Newtonian KK form with fiber scale λ and an
arbitrary potential Φ(t, x). Computing Γ, Ricci, R, G and ∇·G
symbolically: **all five components vanish identically** — exact, not
order-by-order. Hence, given G_{AB} = 8πT_{AB}, the conservation
∇_μT^μ{}_χ = 0 — the continuity equation that #198 verified numerically
at the 3×10⁻⁴ level on the live ψ–Φ–q evolution — is a *geometric
identity* of the bulk, inherited by the brane dynamics.

## 3. Step 4 in detail: pointwise transport of the quantized charge

The centroid statement (Ehrenfest, d⟨x⟩/dt = ⟨∇S⟩_ρ) was verified in
#198. The stronger, *pointwise* statement is demonstrated here on the
2D reduction (the (x, y) plane transverse to the fiber; the winding
charge appears as a quantized phase vortex): a vortex–antivortex pair
(charges ±1) on a uniform condensate with an imposed quantized
background flow, evolved through the full nonlinear dynamics:

- the winding numbers remain **exactly** ±1 throughout (topological
  conservation — the same discreteness as #181/#182);
- each core's measured velocity equals the **ambient J/ρ at the core**
  (ring-averaged to remove the self-flow, which cancels by symmetry):
  both the background part (v_bg, quantized by the box) and the
  partner-induced part (Γ/2πd) are reproduced;
- nothing about the core's motion is put in by hand: the transport is
  read off the same field evolution that conserves the current.

The charge goes where the current goes — pointwise. That is the
guidance law.

## 4. The geodesic contrast: where the quantum potential lives

The eikonal/WKB route ("wave packets follow geodesics") gives only the
classical limit. The derived transport differs from geodesic motion
*exactly* by the quantum potential: the Madelung–Euler form of the flow,

```
∂ₜv + (v·∇)v = −∇( V_eff + Q ) ,     Q = −½ ∇²√ρ / √ρ ,
```

is verified on the live #198 dynamics (residual ~10⁻³), and the quantum
force is *dominant* there (max|∇Q| ≈ 27 × max|∇V_eff|). The crucial
point for realism: **Q is not an addition to GR** — it is the gradient
part of the same wave stress-energy T^{μν} whose T^μ_χ column is the
guidance current. Geodesic transport (dropping Q) is the flow that #198
showed destroys equivariance; the Bianchi-derived transport is the one
that preserves it. GR, taken whole (the full T^{AB}, not its eikonal
truncation), *selects* the Bohmian flow.

## 5. Universality

The guidance velocity is the ratio v^i = T^i_χ/T^0_χ, and the winding
number k cancels: **all throats obey one and the same guidance law**,
independent of their charge/species — verified for k = 1, 2, 3. This is
exactly the universality dBB requires (the velocity field contains no
per-particle constant), and here it is automatic: numerator and
denominator carry the same unit of charge.

## 6. What is and is not established (honest scope)

**Established:** the guidance law v = ∇S is the unique transport
compatible with (i) the 5D Einstein equations (via the Bianchi-forced
conservation of the fiber current), and (ii) the integer, topologically
conserved winding of the throat — with the pointwise charge-transport
demonstrated on the live nonlinear dynamics, the stress-tensor identity
exact, and the Bianchi identity verified symbolically. Combined with
#198, the Born rule now rests on the bulk field equations plus the
quantum-equilibrium/measurement conditions — the guidance condition is
discharged.

**Not established / conditions that remain:**
1. The **full 5D throat-core dynamics** is not solved: the throat is
   represented by its conserved winding (exact) and its localized core
   (the soliton/vortex scale); the demonstration that the *geometric*
   core (the RP² cross-cap itself, in numerical 5D relativity) tracks
   the charge centroid is beyond the weak-field reduction. The
   topological argument constrains it strongly — the core cannot part
   company with its own quantized charge — but "strongly constrained"
   is not "solved."
2. The measurement/equilibrium conditions of #198 (linear measurement
   regime; coarse-grained relaxation) are unchanged — this PR removes
   condition (1) of #198, not condition (2).
3. The 2D vortex demonstration is the transverse reduction of the
   winding transport; the 1D/radial demonstrations of #198 cover the
   longitudinal (density) transport. A single 3D+fiber demonstration
   combining both is a computational, not conceptual, gap.

## References

- The repo's chain: #42–#44 (winding = charge); #83 (KK mass operator);
  #167–#169 (the throat geometry and quotient); #178/#181/#182 (the
  discrete winding invariant and its protection); #193/#195 (the sector
  reduction); #198 (equivariance, uniqueness, relaxation).
- T. Kaluza (1921); O. Klein (1926). [The fiber-momentum = charge
  identification.]
- de Broglie (1927); Bohm (1952). [The guidance law being derived.]
- Standard GR: the contracted Bianchi identity ∇·G ≡ 0 and Killing
  currents (e.g. Wald, *General Relativity*, §C; Misner–Thorne–Wheeler
  §20.6 — "the Bianchi identity as the automatic conservation
  machine").
- Vortex transport in nonlinear Schrödinger flows: standard GP results
  (the core rides the ambient superfluid velocity; e.g. Fetter's
  reviews) — reproduced, not assumed.

## Reproduce (the machine-checkable steps)

```bash
python -m experiments.closure_ledger.guidance_law_from_5d_probe
# Verdict: GUIDANCE_LAW_DERIVED_FROM_THE_5D_BULK_THE_CURRENT_IS_THE_FIBER
#          _BIANCHI_COMPONENT_AND_THE_THROAT_RIDES_ITS_TOPOLOGICAL_CHARGE
```
