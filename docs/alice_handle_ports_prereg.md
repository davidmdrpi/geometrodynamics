# Pre-registration: what crosses an Alice handle, booked at its two mouths

Publish this file before implementing or running anything. Baseline: `main`
at `137c780`. No random numbers are used.

**Prior state.** Nothing below has been computed. The derivation was done by
hand in the session that wrote this file; freeze and implementation are by
one author. Four recent results leave the same gap:

| PR | has | lacks |
|---|---|---|
| #338 | a scalar wave crossing the necks of the handled sphere | a ledger; a torus/Klein comparison of what crosses |
| #340, #354 | mouths that respond; energy flux out of a neighbourhood | any propagation through a handle |
| #352 | exchange rules for a Klein handle | an exterior; it is a 1D ring network |
| #353 | the angular action decides whether circulation is allowed | "an exterior attachment map for a pair of physical mouths" |

## Question

A real scalar wave crosses one handle joining antipodal mouths of a sphere.
Its two rims are the ports. How much energy and how much axial angular
momentum pass each port, and how does that depend on whether the handle is
glued orientation-preservingly (torus-type) or orientation-reversingly
(Klein-type, the "Alice" handle)?

## Model

The surface of #338 with ONE handle on the z-axis, unit sphere, mouth angle
`a`. Sphere band `a ≤ θ ≤ π − a` with azimuth `φ`. Neck `|s| ≤ S`,
`S = sin a cos a`, angle `ψ`, metric `ds² + r(s)² dψ²`, `r = sqrt(s² + b²)`,
`b = sin² a`. Rim A is `θ = a ~ s = −S`; rim B is `θ = π − a ~ s = +S`. The
field obeys `u_tt = Δu` on the glued surface, continuous with continuous
normal derivative across both rims. Geometry is fixed: nothing moves.

Gluing at rim A is `ψ = φ` in every case. At rim B:

* torus-type, twist `c`: `ψ = φ + c`;
* Klein-type: `ψ = −φ`.

## Derivation

**Orientation.** With the outward normal on the band, the induced boundary
direction is `−φ̂` at rim A and `+φ̂` at rim B; on the neck it is `−ψ̂` at
`s = −S` and `+ψ̂` at `s = +S`. `ψ = φ + c` matches direction at both rims,
so the glued surface is orientable. `ψ = −φ` matches at A and reverses at B,
so it is not.

**Symmetry.** A rotation `φ → φ + α` extends to the neck as `ψ → ψ + α` and
preserves `ψ = φ + c`. It does not preserve `ψ = −φ`. So torus-type gluing
keeps the axial `U(1)` for every twist `c`, and Klein-type gluing breaks it
to the reflection `φ → −φ, ψ → −ψ`. A constant in `ψ = −φ + c` is removed by
rotating band and neck together by `c/2`: the Klein handle has no twist
modulus, only a reflection axis.

**Reduction.** Along the closed loop `ℓ` (band, then neck, total length
`P = π − 2a + 2S`) the circle radius `ρ(ℓ)` is `sin θ` on the band and `r(s)`
on the neck; it is `C¹` at both rims. Write

    u = Σ_m [ A_m(ℓ, t) cos mφ + B_m(ℓ, t) sin mφ ],

with `ψ` continued from rim A on the neck. Each of `A_m`, `B_m` obeys
`w_tt = (1/ρ)(ρ w')' − m² w/ρ²` on the loop, with one junction condition at
rim B:

* torus-type: `(A + iB)` is multiplied by `e^{imc}` across rim B;
* Klein-type: `A` is continuous and `B` changes sign.

**Angular momentum.** With density `−u_t u_φ` and flux `u_ℓ u_φ`,

    L[ℓ₁, ℓ₂] = π m ∫ (A Ḃ − Ȧ B) ρ dℓ,      Φ(ℓ) = π m ρ (A′B − B′A),

`Φ` being the flux in the `+ℓ` direction. Let `Φ_A` be the flux into the
neck at rim A and `Φ_B` the flux out of the neck at rim B, both in neck
variables. Then `dL_neck/dt = Φ_A − Φ_B`, and on the band the flux entering
at rim B is `Φ_B` for torus-type and `−Φ_B` for Klein-type, because `B`
changes sign there. Hence

    d(L_band + L_neck)/dt = 0            torus-type, any c,
    d(L_band + L_neck)/dt = −2 Φ_B       Klein-type.

The second line is bookkeeping. Its content is whether `∫Φ_B dt` is nonzero.

## Predictions

Separable solver: for each `m`, linear finite elements on the loop with `N`
uniform nodes and the junction condition at rim B; generalized eigenproblem;
time evolution by mode sum, so there is no time-stepping error. Geometries
`a = 0.3` (the repository default) and `a = 0.5`. `m = 0, 1, 2, 3`.
`N = 2000`, repeated at `N = 4000`.

**P1. `m = 0` cannot tell the gluings apart.** The `m = 0` spectra of
torus-type (any `c`) and Klein-type are identical to `1e-10` relative, and
so is the energy that crosses each rim for `m = 0` initial data.

**P2. Torus-type conserves axial angular momentum.** For `c = 0` and
`c = 1`, and for generic initial data in `m = 1, 2, 3`, `L_band + L_neck` is
constant to `1e-9` of its scale, and every `m ≥ 1` eigenvalue is doubly
degenerate to `1e-10` relative.

**P3. Klein spectrum is a union of two torus spectra.** For each `m ≥ 1`,
the Klein-type spectrum is the `c = 0` torus spectrum together with the
`c = π/m` torus spectrum, each taken once, to `1e-10` relative. The cosine
sector carries the first and the sine sector the second.

**P4. Klein doublets split, and the split is the reversal rate.** Pair the
`n`-th cosine-sector mode `(ω₊, v₊)` with the `n`-th sine-sector mode
`(ω₋, v₋)`, and write `Δ = |ω₊ − ω₋|`, `ω̄ = (ω₊ + ω₋)/2`. For each `m ≥ 1`
the lowest six pairs are split, `Δ > 1e-6 ω̄`. For every such pair with
`Δ ≤ 0.2 ω̄`, take `A = v₊` at rest and `B = 0` with velocity `ω₋ v₋`. Then
`L(0) ≠ 0`, and `L_band + L_neck` first changes sign at `t = π/(2Δ)`, to
within a relative `Δ/ω̄` (the bound that the fast term
`−(Δ/2) cos(2ω̄t)` allows).

**P5. The identity, and it is not vacuous.** For Klein-type runs,
`L(t) − L(0) + 2∫₀ᵗ Φ_B dt` vanishes to `1e-8` of the scale of `L`, with
`Φ_B` the discrete flux for which the semi-discrete continuity equation is
exact. In every run of P4, continued to `t = π/Δ`, `|∫Φ_B dt|` reaches at
least `10%` of `|L(0)|`.

**P6. What comes through the handle arrives reversed.** `a = 0.5`, `m = 1`
and `m = 2`. Initial data: a compactly supported packet on the band,
`θ ∈ [a + 0.05, a + 0.45]`, with a `cos²` envelope, moving toward rim A,
carrier wavenumber `14` along `ℓ` (above the neck's barrier `m/b ≤ 8.7`),
with `L(0) > 0`. Region `R`: the band within `0.25` of rim B.
Window: from first arrival through the neck until `t = π − 2a − 0.7`, the
earliest time anything can reach `R` round the sphere at unit speed. In that
window, `L_R` of the Klein-type run equals minus `L_R` of the `c = 0`
torus-type run, to `1e-3` of the largest `|L_R|` in the window, and that
largest value exceeds `1e-4 |L(0)|`.

**P7. The #338 mesh agrees, with no fitted parameter.** `build_surface`
with one handle, `n_sphere = 6000`, `n_ring = 48`, both gluings.

* Geometry, checked before any spectrum: across glued rim vertices, the
  far-rim azimuth minus the near-rim azimuth is constant for the torus-type
  mesh, and that constant is its `c`; their sum is constant for the
  Klein-type mesh.
* The lowest 20 nonzero eigenvalues `ω²` of the mesh operator agree with the
  separable solver's, all `m` merged and sorted, to `2%`, and the largest
  discrepancy decreases at `n_sphere = 12000`, `n_ring = 64`.

**Expectation, not a gate.** Below the neck's barrier, `ω < m/b`, the
splitting should be of order `exp(−2m asinh(cot a))`: about `2%` for `m = 1`
at `a = 0.3` and larger at `a = 0.5`. Reversal is therefore slow at the
default geometry.

## Verdicts

`CONFIRMED` if P1 to P7 hold. `COUNTEREXAMPLE` if any structural statement
fails at both resolutions: P1, the degeneracy in P2, P3, the split in P4, or
the sign in P6. `UNRESOLVED_NUMERICALLY` if a tolerance is missed but the
miss shrinks at the finer resolution; one further doubling of `N` is allowed
and no other change. `INCOMPLETE` if P7's geometry check fails, since then
the mesh is not the modelled surface. Tolerances are not adjusted after the
fact; a missed tolerance is reported as missed.

## Scope

* A scalar on a fixed 2D surface. Not a GR solution and not on `S³`.
* The reversal of `m` is a consequence of a broken symmetry. It is classical
  and is not claimed to be quantum.
* A scalar has no normal component, so this addresses the angular-action
  column of #353 and its missing exterior, not its normal-sign column.
* The resemblance to Alice strings and Cheshire charge is an analogy.

## A second gate, to be frozen separately

With a fixed Klein handle the angular momentum that P5 books is lost to the
rigid gluing. If the reflection axis `χ` is made dynamical, with the gluing
`ψ = −φ + 2χ` and a kinetic term `½ I χ̇²`, then rotating band, neck and axis
together is a symmetry, and Noether's theorem gives a conserved
`L_band + L_neck + I χ̇`. The handle would then carry what the field loses.
This is stated here as the intended next step. It is not part of this
freeze: it gets its own pre-registration, written only after the verdict
above, and only if that verdict is `CONFIRMED`.
