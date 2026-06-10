# Electric and magnetic form factors: the EM gauge-arc capstone (PR #147)

> **Framing (to avoid a category error).** QFT on the *fixed classical* throat
> geometry — geometry → fields, **not** quantum gravity. The form factors are
> those of the matter QFT's dressed charge and moment on the classical
> antipodal cavity.

The EM gauge arc is complete through one loop except for one assembly: #141
built the vertex, #142 the Ward identity, #143 the α ledger, #144 the photon Π
and the running, #145 the exact charge (`Z₁ = Z₂`, `F₁(0) = c₁`), #146 the
electric form factor and its geometric radius. The remaining face is the
relativistic decomposition

    Γ^μ(q) = γ^μ F₁(q²) + (i σ^{μν} q_ν / 2m) F₂(q²),

whose two scalars are the electric (charge) and magnetic (moment) form
factors, with `g/2 = G_M(0)/c₁ = 1 + F₂(0)`. This PR assembles it on the
cavity, re-verifies the arc's magnetic keystones (#61 g = 2, #62 a = α/2π)
alongside the electric ones (#145/#146) in the #131 capstone convention, and
exhibits the single Dirac-algebra fact behind the arc's asymmetry: **why the
charge is exact while the moment is dressed.**

## The Gordon decomposition (the algebraic heart)

`ū(p′)γ^μu(p) = ū(p′)[(p+p′)^μ + iσ^{μν}q_ν]u(p)/2m` — verified with explicit
Dirac spinors over random on-shell momenta to **~1e-15** (Clifford algebra
exact). The electric/magnetic split is not an ansatz; it is the Dirac algebra
of the #141 minimal coupling.

## Why F₁ is pinned and F₂ is free — one identity

The Ward contraction `q_μΓ^μ` kills the F₂ term twice over:
`q_μσ^{μν}q_ν = 0` identically (antisymmetry — **exact**, verified) and
`ū(p′)q̸u(p) = 0` on shell (**~1e-16**). So the #142/#145 Ward identity
constrains F₁ only: the charge `F₁(0) = c₁` is exact and coupling-independent,
while F₂ — the anomalous moment — is gauge-free and dresses at every loop. One
identity explains both the exact charge and the dressed moment.

## The keystones, re-verified together

- **Tree, g = 2 (#61):** the Pauli/Hopf operator identity
  `(σ·D)² = D² − σ·B` verified by finite-difference application to random
  spinor test functions in a constant-B field (**~1e-6**). The SU(2)
  anticommutator factor of 2 IS `g_s = 2`; `F₂(0) = 0` at tree level.
- **Loop, F₂(0) = α/2π (#62):** the Schwinger Feynman-parameter simplex
  integral `∫ dz dy 2z(1−z)/(1−z)² = ∫₀¹ 2z dz = 1` (numerically 0.9999998)
  ⟹ `a = α/2π = 0.00116141` vs the measured `a_e = 0.00115965` (+0.15%, the
  α² Sommerfield term and beyond); `g = 2(1 + α/2π) = 2.0023228` vs
  `2.0023193`.

## The Sachs assembly

`G_E(q)` is the #146 dressed-density transform (geometric radius
`r_E = 0.2649` tortoise units). The magnetization density rides the same
charged-mode profile (the spin is carried by the charged line; the neutral φ
is spinless), so `r_M = r_E` exactly and `G_M(q)/G_M(0) = G_E(q)/G_E(0)` —
form-factor scaling in the minimal model. At q = 0 the arc's asymmetry is
explicit: `G_E(0) = c₁ = 1` **exact** (Ward-pinned), `G_M(0) = 1 + α/2π =
1.0011614` **dressed** (Ward-free).

## The arc, one primitive

Every face — minimal coupling + Z₂ selection (#141), Ward/masslessness (#142),
the α ledger (#143), Π and the running (#144), exact universal charge (#145),
G_E and the geometric radius (#146), the Gordon split + Ward F₂-freedom +
g = 2 + α/2π (#147) — derives from the **unitary antipodal throat carrying the
integer Hopf charge** (#129/#58). The single EM input is the value `α(μ₀)`
(#143).

## Scope and epistemic ledger

- **Derived:** the Gordon split (exact), the Ward F₁/F₂ asymmetry (exact), the
  g = 2 operator identity, the α/2π simplex value, the geometric radii.
- **Modelled:** the F₂(q) shape (flat ⟹ scaling); the #136-posture couplings.
- **Input:** the boundary value `α(μ₀) ≈ 1/137` (#143).
- **Open:** the α² Sommerfield term; the measurable `r_E − r_M` splitting;
  recoil / O(q²/m²) F₁–F₂ mixing; the absolute normalisation (#133); the
  flavor residuals (#134).

## Reproduce

```bash
python -m experiments.closure_ledger.em_form_factor_decomposition_probe
# Verdict: EM_FORM_FACTORS_F1_WARD_PINNED_F2_ALPHA_OVER_2PI_RADII_GEOMETRIC_CAPSTONE
```
