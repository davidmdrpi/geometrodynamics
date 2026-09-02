# The finite mouth, its gluing, and where the transport lives

*Derivation document for the finite-mouth topology round. Pre-registered in
`docs/finite_mouth_topology_prereg.md` (commit `d9d85bc`), implemented in
`geometrodynamics/bulk/mouth_topology.py`, checked by
`tests/test_mouth_topology.py` and `experiments/closure_ledger/finite_mouth_topology_probe.py`.
Nothing here uses a singlet, a Born rule, a projector, a tensor product, CHSH
or a QED amplitude.*

## The question

`docs/qft_emergence_audit.md` localised the missing chain

```
finite physical mouth  →  orientation/bundle transport  →  antipodal field boundary condition
```

The repository has the three objects: the PR #277 handle, `J = iσ_y K` in
`embedding/transport.py`, and the PR #129 condition
`Φ(U,V,Ω) = Φ(−U,−V,Ω̄)`. This document asks whether the first forces the
other two, and answers with the gluing data written out.

**Verdict:** `FINITE_MOUTH_ADMITS_BUT_DOES_NOT_SELECT_THE_BAM_LIFT`, with two
sharpenings: the antipodally glued two-mouth handle is *orientable* in the
bulk, and `J` is a *rotation*, not a reflection, whose antilinearity is a
choice of complex structure.

## 1. The handle and its gluing data (H1)

Start from PR #277: `Σ = S⁴_R`; excise geodesic 4-balls of angular radius `a`
at `P_A, P_B`; join their `S³` boundaries by the tube `T = [−S,S] × S³`,
`dΣ² = ds² + f²dΩ₃²`, `f = √(s²+b²)`, `b = R sin²a`, `S = R sin a cos a`.
Let `C = S⁴_R ∖ (B_A ∪ B_B)`. Then `M = C ∪ T` is a closed 4-manifold.

**The seam maps are free.** Darmois matching at a mouth is two scalars,
`f_m = R sin a` and `f′ = cos a`, both invariant under the isometry group
`O(4)` of the round mouth `S³`. So the identification of the tube's `S³` with
the ball boundary at each seam is an arbitrary `ψ_A, ψ_B ∈ O(4)`. Isometries
of `T` act by a single `g ∈ O(4)` on every slice, and isometries of `C` fixing
both balls act by the same `g` on both boundaries, so only the loop monodromy
`m = ψ_B⁻¹ψ_A` matters, up to conjugation. Nothing in the geometry fixes it.

**Topology.** `C ≅ S³ × [0,1]`, so `M` is the mapping torus of `m`. Since
`π₀ Diff(S³) = π₀ O(4) = Z₂` (Cerf), the diffeomorphism type is fixed by
`det m`: `S³ × S¹` for `det m = +1`, the non-orientable `S³ ×̃ S¹` for
`det m = −1`. `w₁(M)` on the handle loop is `det m` (`mapping_torus_w1`).

**The antipodal map is orientation-preserving.** `−I₄` on `R⁴ ⊃ S³` has
`det = +1`. So an antipodally glued two-mouth handle is **orientable**
(`classify_gluing(antipodal_gluing()).det_bulk = +1`). Whatever the
repository's "non-orientable throat" means, it is not this.

**What the antipode does twist.** The brane is the equatorial `S³ ⊂ S⁴`; if
the ball centres lie on it, one of the four mouth directions is the brane
normal and the other three span the brane mouth `S²`. A brane-preserving
monodromy is `m = (m₃, ε) ∈ O(3) × O(1)`, giving four classes with

```
w₁(bulk handle)          = det m₃ · ε
w₁(brane handle)         = det m₃         (S² × S¹ or S² ×̃ S¹)
w₁(brane normal bundle)  = ε
```

The antipodal class is `(−, −)`: the brane `S²`-handle is non-orientable
(`−I₃` reverses `S²`), its normal line bundle is twisted, and their product
leaves the bulk orientable. This is the `RP^{d}` versus `RP^{d−1}` parity
that `viz/hyperspherical.py` tabulates, now located in the gluing.

**Antipodal versus identity is geometry, not topology, on the handle.** `−I₄`
is isotopic to `I` inside `SO(4)`, so `M_{−I}` and `M_I` are diffeomorphic.
They are not isometric: the twist is a modulus, and its only static signature
is the action of `m` on harmonics. That action is *computed*
(`harmonic_representation`, a least-squares fit of `p ∘ m` on the
`(ℓ+1)²`-dimensional space of degree-`ℓ` harmonics): for `m = −I` it is the
scalar `(−1)^ℓ` with residual `< 1e-15`; for a reflection it is not a scalar
at all (trace `3` on the 9-dimensional `ℓ = 2` space). So the `(−1)^ℓ` the
repository uses is the antipodal parity of harmonics and nothing more.

**One traversal, two traversals.** Around the handle loop the transverse `S³`
is transformed by `m` once and by `m²` twice. For `m = ±I`, `m² = I`. For the
map in `transport.py`, `σ = L_{−j}`, `σ² = −I₄` is the antipodal map and
`σ⁴ = I`: two passes of `σ` are the antipode, not the identity.

## 2. The involution (H2)

Suppose `P_B = −P_A` on `S⁴`. This is **not** specified in PR #277 and is a
choice. Then the antipode `x ↦ −x` of `S⁴` preserves `C` and swaps the
balls. At the seams it acts as `Ω ↦ −Ω`, and the only isometry of `T` of the
form `(s, Ω) ↦ (±s, gΩ)` matching `−I₄` at both ends is

```
ι(s, Ω) = (−s, −Ω).
```

*Uniqueness.* An involution in `O(5)` is conjugate to `diag(±1, …)`; it has a
fixed point on `S⁴` iff some eigenvalue is `+1`. Only `−I₅` is free
(`free_involutions_of_o5`). So `ι` is the unique free isometric involution
of `M` that swaps the mouths.

*Properties* (`antipodal_involution`): free (`−I₄` has no eigenvalue `+1`,
so no fixed point at `s = 0`); an isometry for both lapses because `f` and
`N` are even in `s`; `det dι = (−1)(+1) = −1` on the bulk handle;
`det dι|_brane = (−1)(−1) = +1` on the brane slice `[−S,S] × S²`. This is
exactly the PR #196 computation for the `RP³` geon, now with its bulk
counterpart: **the brane quotient is orientable, the bulk quotient is not.**

*The quotient.* `C/ι = RP⁴ ∖ B⁴`; `T/ι` is the twisted `I`-bundle over the
neck `S³/± = RP³`, which is again `RP⁴ ∖ B⁴`; glued along the single mouth,
`M/ι = RP⁴ # RP⁴`. Its brane slice is `RP³ # RP³` (orientable, spin), its neck
is `RP³`, its brane neck is `RP² = S²/±`. The words "`RP²` mouth",
"non-orientable throat" and "Möbius wrap" all refer to **this quotient**, not
to the two-mouth handle.

*Pin types, corrected after review.* With `w(RPⁿ) = (1+a)^{n+1}`: `RP²` has
`w₁ = a`, `w₂ = a²`, so it is Pin⁻ only; `RP⁴` has `w₁ = a`, `w₂ = 0`,
`w₁² = a² ≠ 0`, so it is Pin⁺ only, and so is `RP⁴ # RP⁴`. **A first draft
of this document called that a mismatch. It is not; it is the standard
mechanism by which an ambient Pin⁺ structure induces an intrinsic Pin⁻
structure on a surface whose normal bundle is twisted.** The deck generator
of the neck `RP²` reverses both normal directions (`∂_s` and the brane
normal `n`), so `ν = λ ⊕ λ`, `w(ν) = 1 + a²`, and

```
w(TM|_N) = w(TN) w(ν) = (1 + a + a²)(1 + a²) = 1 + a         (a³ = 0)
```

so `w₂(TM|_N) = 0` (ambient Pin⁺ restricts) while `w₂(TN) + w₁(TN)² = 0`
(intrinsic Pin⁻). See §5b and `docs/mouth_pin_holonomy_prereg.md` (R1).
*Evidence class:* topological theorem, conditional on `P_B = −P_A`.

## 3. The scalar boundary condition (H3)

A scalar on `M/ι` is a scalar on `M` with `Φ(ιx) = ηΦ(x)`, `η = +1` for a
function and `η = −1` for a section of the orientation line bundle. With
`Φ = Σ ψ_ℓ(s) Y_ℓ(Ω)` and the computed `Y_ℓ(−Ω) = (−1)^ℓ Y_ℓ(Ω)`:

```
ψ_ℓ(−s) = η (−1)^ℓ ψ_ℓ(s)
```

so at the neck `s = 0`:

| `η (−1)^ℓ` | neck condition |
|--|--|
| `+1` | `ψ′_ℓ(0) = 0` (Neumann) |
| `−1` | `ψ_ℓ(0) = 0` (Dirichlet) |

For `η = +1` this is PR #129's even-Neumann / odd-Dirichlet structure,
obtained at the finite ultrastatic neck, `N = 1`, with no horizon anywhere.
For `η = −1` the two swap: the "twisted/Möbius field" PR #129 mentions in
passing is the other line bundle. **`η` is a choice the geometry does not
make.**

**The half-tube admittance.** On `M/ι` the two mouths are one mouth, and the
Dirichlet-to-Neumann map is a scalar per `ℓ`: the outward flux at the mouth
for unit data there with the neck condition. It is the `(1, η(−1)^ℓ)`
eigenvalue of PR #277's two-mouth `Y_ℓ(0)`:

```
Y_ℓ^N = 2π²F² [ k tanh(kX) − cos a ]        Y_ℓ^D = 2π²F² [ k coth(kX) − cos a ]
```

(`k = ℓ+1`, `F = R sin a`, `X = arcosh(1/sin a)`; from
`coth 2x ∓ csch 2x = tanh x, coth x`). `Y_0^N = 0` exactly, because
`tanh X = cos a` and the Neumann monopole is the constant; `Y_0^D = 2G`.

An independent conservative tridiagonal solve on `[0, S]` — Neumann by a
half-cell balance at the neck, Dirichlet by dropping the neck unknown, never
using the `s = b sinh x` reduction — reproduces the pre-registered targets at
`R = 1, a = 0.3` to `≤ 5e-7` relative at 2000 steps and converges at second
order (ratios `4.00, 4.00, 4.00`):

| `ℓ` | `Y^N` target | `Y^D` target |
|--|--|--|
| 0 | `0.000000000` | `0.157587622` |
| 1 | `1.797266559` | `1.804461992` |
| 2 | `3.524607516` | `3.524854051` |
| 3 | `5.248595411` | `5.248602920` |

**Controls.** *Antipodal:* replace `Ω ↦ −Ω` by the identity. The involution
`(s, Ω) ↦ (−s, Ω)` fixes the whole neck (not free; the quotient is a manifold
with boundary) and every `ℓ` gets the same condition: the grading is gone.
*Geometry:* across `a ∈ {0.05, 0.3, 0.8, 1.2}` the labels are invariant and
`Y_1^D` runs from `0.049` to `47.3`. *Topology:* identity gluing `w₁ = +1`,
reflection gluing `w₁ = −1`, antipodal gluing `w₁ = +1`.

## 4. The horizon "limit" (H4)

There is none to take. The Tangherlini branch is the same spatial handle with
`N = |s|/f`; `f` and `N` are both even in `s`, so the static radial operator
commutes with `s ↦ −s` for both lapses (residual `< 3e-11`), and the neck
sector labels are lapse-independent. The neck is the bifurcation surface of
the Tangherlini extension, and `s ↦ −s` at fixed `t` there is the Kruskal map
`(U, V) ↦ (−U, −V)`; `Ω ↦ −Ω` is `Ω ↦ Ω̄`. So PR #129's condition **is**
`ι`-invariance with `η = +1`.

One datum PR #129's map carries that `ι` does not: `(U,V) ↦ (−U,−V)` is a
rotation by `π` in the Lorentzian `(U,V)` plane, `det = +1`, which reverses
the time orientation `U + V`. The reflection `(U,V) ↦ (V,U)`, `det = −1`,
also swaps the exteriors and also restricts to `ι` on `t = 0`, but keeps the
time orientation. The PR #129 extension is PT-type (5D `det = +1`, not
time-orientable quotient); the alternative is P-type (5D `det = −1`). Which is
physical is a further choice, recorded as the datum `T`. For the static
problem the two agree.

## 5. The spinor lift, and what `J` is (H5)

**Built from the quaternions, not from Pauli matrices.** `R⁴ = H` with
`L_q, R_q` left and right multiplication (`quaternion_left`,
`quaternion_right`, from the multiplication table). `Spin(4) = SU(2)_L × SU(2)_R`
acts on vectors by `v ↦ g_L v g_R⁻¹` and on the two half-spinors `S₊ = H`,
`S₋ = H` by `L_{g_L}` and `L_{g_R}`; the half-spinor complex structure is
`R_i`, the one that commutes with the left action.

**What the map in `transport.py` is.** Reading `σ(z₁, z₂) = (z̄₂, −z̄₁)` in
the coordinates `q = z₁ + z₂ j`, `hopf_transport_matrix` builds its real
`4×4` matrix from the function itself and finds it equal to `L_{−j}`:
left multiplication by `−j`, an isoclinic rotation by `π/2`, `det = +1`,
`σ² = −I₄`. With `h(q) = q⁻¹ i q` (the left-fibre Hopf map), `h(σq) = −h(q)`
and `σ(e^{iφ}q) = e^{−iφ}σ(q)`. So `σ` is the base antipode composed with
fibre reversal, orientation-preserving on `S³`, as the audit found.

**Which complex structure makes it `K`.** `L_{−j}` anticommutes with `L_i`
and commutes with `R_i` (`complex_structure_commutation`). Therefore:

* in the **Hopf** complex structure `L_i` (fibre phase = scalar
  multiplication, the reading BAM uses when it identifies a state with a
  point of `S³ ⊂ C²`), `σ` is antilinear: `iσ_y K`;
* in the **half-spinor** complex structure `R_i`, `σ` is the complex-linear
  matrix `[[0, 1], [−1, 0]] = iσ_y ∈ SU(2)`.

The `K` is a statement about the state space's complex structure, not about
the geometry. `J = iσ_y K = (spin rotation by π about y) ∘ (reversal of the
Hopf `U(1)`)`. The second factor is charge conjugation on the `U(1)` bundle,
a `Z₂` bundle datum; the first is `−j ∈ SU(2)_L`, a rotation.

**What the gluing maps lift to.** In `Cl(4)`, built as its own regular
representation (`clifford_regular`):

* the neck reflection `dι = diag(−1, +1, +1, +1)` lifts to `±e_s`, with
  `e_s² = +1` in Pin⁺ and `−1` in Pin⁻: **four** lifts, each a real matrix
  anticommuting with the volume element, hence exchanging the two chiralities
  (`pin_lifts_of_reflection`);
* the antipode `−I₄ ∈ SO(4)` lifts to `±e₀e₁e₂e₃ = ±` volume element, with
  eigenvalues `±1`: the chirality sign (`spin_lifts_of_antipode`).

**Corrected after review (R2 in `docs/mouth_pin_holonomy_prereg.md`).** A
first draft said "no lift of any element of `O(4)` is `J`". That was too
broad, and wrong: `σ = L_{−j} ∈ SO(4)` has the Spin(4) lift `(−j, 1)` (up to
the central sign), which acts on one half-spinor as the linear quaternion
`−j`, i.e. as `iσ_y` in the `R_i` basis. So `J` **is** the spin lift of a
geometric rotation, and that rotation is an admissible seam monodromy (the
fibre-reversing Hopf-isometry class). What remains true and narrower:

* `J` is not a lift of `ι` (whose lifts are the odd, chirality-exchanging
  `±e_s`) nor of `−I₄` (whose lifts are `±` the volume element);
* the geometry does not select `m = σ` over `m = I` or `m = −I`;
* the `K` is a statement about complex structure: one real map, antilinear
  for `L_i`, linear for `R_i`. "The Pin matrices are real, therefore
  complex-linear" does not decide which complex structure is physical.

Whether the *mouth* — rather than a chosen seam monodromy — produces an
element of this class is the question of §5b.

**Lift control.** The four lifts of `ι` are pairwise distinct within a Pin
type and differ in square across types; the two lifts of `−I₄` differ by
sign. Uniqueness of "the" transport is false even before `J` is considered.

## 5b. The mouth holonomy, and where the order-four structure comes from

*Pre-registered in `docs/mouth_pin_holonomy_prereg.md` (`7f46fff`) after
review; results in the corresponding section of
`finite_mouth_topology_probe` and `tests/test_mouth_topology.py`.*

**What the deck generator does to the frame (P1).** At a neck point
`Ω ∈ S² ⊂ S³_{s=0}` take the frame `(∂_s, n, t₁, t₂)`, `n` the brane normal
inside the mouth `S³`, `t₁, t₂` tangent to the neck `S²`. The generator of
`π₁(RP²)` lifts to the great semicircle from `Ω` to `−Ω` in the `t₁`
direction, followed by the deck map `ι`. Parallel transport along the
semicircle, integrated from `v′ = −(v·γ′)γ` on the round neck
(`transport_along_great_semicircle`), followed by `dι` (`∂_s ↦ −∂_s`, and
every `S³`-tangent vector `v ↦ −v`), gives

```
H = diag(−1, −1, +1, −1)      in (∂_s, n, t₁, t₂),      det H = −1
```

Both normals reverse; the tangent block is a reflection (the Möbius
monodromy of `RP²`). So `ν = λ ⊕ λ`. *Class:* analytic identity, computed.

**Stiefel–Whitney (P2).** In `Z₂[a]/(a³)`: `w(ν) = 1 + a²`,
`w(TN) = 1 + a + a²`, `w(TM|_N) = 1 + a`. So `w₂(TM|_N) = 0` — the ambient
Pin⁺ structure restricts to the neck — while `w₂(TN) + w₁(TN)² = 0` and
`w₂(TN) ≠ 0`: intrinsically Pin⁻, not Pin⁺. This is the retraction R1: the
two Pin types are the same structure seen from the ambient and intrinsic
sides.

**The conversion (P3).** With `e_s, e_n` the normal generators of `Cl(4)`
(`e² = +1`), the twisted tangent generators `ẽ_t = e_t e_s e_n` satisfy
`ẽ_t² = −1` and `ẽ_{t₁}ẽ_{t₂} = −ẽ_{t₂}ẽ_{t₁}`: they generate `Cl⁻(2) ≅ H`
(`twisted_tangent_generators`, in the regular representation). The normal
volume element has `(e_s e_n)² = −1`: the Spin(2) lift of the normal
`π`-rotation, of order four. **Not inserted.**

**The lift (P4).** `H` is the product of the reflections in `∂_s`, `n`, `t₂`,
so its Pin lift is `±e_s e_n e_{t₂} = ±ẽ_{t₂}`: odd in `Cl(4)`, hence
chirality-exchanging, with

```
H̃² = −1   in Pin⁺(4),        H̃² = +1   in Pin⁻(4).
```

**Consistency with the cover (P5).** `H̃²` is the spin holonomy of the closed
great circle `γ ∪ ι(γ)` on the double cover. Computed independently by
transporting a tangent vector around latitude circles of the neck `S²` and
unwrapping the rotation angle continuously to the equator
(`closed_loop_spin_holonomy`): the angle follows `2π(1 − cos θ₀)` to
`2.6e-11` and reaches `2π` at the equator, whose spin lift is `−1`. So the
Pin⁺ convention is the one coherent with the cover's spin structure; Pin⁻
would give `+1` and is excluded twice over (by `w₁² ≠ 0` and by this).
**The order-four structure is the spin holonomy of a `2π` tangent rotation
on the round neck, and the deck-generator holonomy is its square root.**

**Comparison with `J` (P6), at three levels.**

| level | mouth holonomy | `J` | equivalent? |
|--|--|--|--|
| vector, `O(4)` | `H`: `det −1`, eigenvalues `−1,−1,−1,+1` | `σ = L_{−j}`: `det +1`, isoclinic | no |
| ambient spinor, `Pin(4)` | `±e_s e_n e_{t₂}`, odd, chirality-exchanging | `(−j, 1)`, even | no |
| intrinsic, `Pin⁻(2) ⊂ SU(2)` | `L_j` (with `ẽ_{t₁} ↦ i`, `ẽ_{t₂} ↦ j`, `ẽ_{t₁}ẽ_{t₂} ↦ k`) | `L_{−j}` | **yes, up to `Spin(2)` conjugation (angle `π/2`) and sign** |

`Cl⁻(2) ≅ H` canonically, acting on itself by left multiplication;
`Pin⁻(2)` is generated by the unit vectors `cos α·i + sin α·j` and
`Spin(2) = {e^{θk}}` — exactly the normaliser `Pin(2) ⊂ SU(2)_L` that §5
found as the Hopf-fibre-preserving isometry group. The mouth holonomy is
left multiplication by a unit quaternion in the `(i, j)` plane, of order
four; `σ = L_{−j}` is one element of that component
(`compare_mouth_holonomy_with_transport`). **Outcome B**: the mouth supplies
`J² = −1` and the component; it leaves a `U(1)` direction (the tangent
direction of the generator path, a `Spin(2)` conjugation) and a sign.

**The `K` (P7).** In the module `H`, `L_j` anticommutes with the `Spin(2)`
generator `L_k` and commutes with `R_i` (`intrinsic_pin2_module`). So the
mouth holonomy is antilinear for the complex structure `L_k` and linear for
`R_i`. If the Hopf fibre phase of the mouth is identified with the `Spin(2)`
rotation of the mouth's tangent plane, the transport is antilinear
**canonically**: the `K` is the statement that `Pin⁻(2) ∖ Spin(2)`
anticommutes with the `Spin(2)` generator. That identification (fibre phase
= mouth tangent rotation) is physical and unproved: *conditional*.

**Selection (P8).** `RP²` carries two Pin⁻ structures, `RP⁴ # RP⁴` four Pin⁺
structures, `S³ × S¹` two spin structures. The sign of `H̃` is not fixed by
the metric or by `ι`; it is the choice of Pin structure. This is the
quantity `embedding/topology.py` stores as `wrap_parity`: now identified,
still not selected.

**What this changes in the audit table.** "Is `J` a valid lift?" — yes: of
the rotation `σ`, and the induced Pin⁻ holonomy of the mouth lies in the
same component with the same square. "Is it unique?" — no: up to a `U(1)`
and a sign. "Where does the order-four structure come from?" — from the
ambient Pin⁺ structure and the two twisted normal lines, i.e. from the spin
holonomy of the round neck; forced, not inserted.

*Dependency ledger:*

```
mouth holonomy H̃  =  H̃( ι [derived given P_B = −P_A], ν = λ⊕λ [derived from ι],
                         Pin⁺ [derived: w₁² ≠ 0 and P5], t₂ [gauge: Spin(2)],
                         sign [chosen: Pin structure] )
"H̃ = J"           =  ( H̃, Cl⁻(2) ≅ H [canonical], H ≅ Hopf C² [chosen complex structure],
                        t₂ = −j [chosen within U(1)], sign [chosen] )
```

## 5c. The Hopf bundle is the spin-frame bundle of the mouth, and the sign pairs

*Pre-registered in `docs/mouth_spin_frame_prereg.md` (`6bc4306`);
module `geometrodynamics/bulk/mouth_spin_frame.py`; tests
`tests/test_mouth_spin_frame.py`; probe section H7.*

**The `U(1)` direction was gauge.** The unfixed direction of §5b is the
choice of a unit imaginary quaternion `u` orthogonal to the fibre generator;
fibre rotations conjugate any such `u` to any other. The geometry selects the
conjugacy class `{u ∈ Pin⁻(2) ∖ Spin(2) : u² = −1}`, which is all an invariant
statement can select. Removed from the unselected column before this round.

**T1 — the frame map.** With the brane mouth `S² = {x ∈ Im H : |x| = 1}`,
the bulk mouth `S³ = {q ∈ H : |q| = 1}`, the brane normal at the mouth centre
identified with `1 ∈ H` (geometric: the pole of the mouth `S³` in the
brane-normal direction) and a reference brane direction `i` (gauge),

```
q ↦ (x, e₁, e₂) = (q⁻¹iq, q⁻¹jq, q⁻¹kq)
```

lands in the oriented orthonormal frame bundle `F_SO(S²) ≅ SO(3)`
(orthonormality, tangency and `det(e₁,e₂,x) = +1` to `1e-15`), with `q` and
`−q` coinciding and distinct `q ≠ ±q′` giving distinct frames: it is
`Spin(3) → SO(3)`. **The bulk mouth `S³` is the spin-frame bundle of the
brane mouth `S²`.**

**T2 — the Hopf fibre is `Spin(2)`.** Under `q ↦ e^{iφ}q`, `x` is fixed and
`(e₁, e₂)` rotates by `2φ` (residual `9e-16`). The Levi-Civita form of the
round `S²`, pulled back through the frame map, is `ω = ⟨ė₁, e₂⟩ = −2A` with
`A = ⟨iq, q̇⟩` the Hopf connection (residual `1.7e-9` on random paths;
analytically `ė₁ = q⁻¹[j, q̇q⁻¹]q = −2w₁e₂ + 2w₃x`, `A = w₁`). Consequently
the repository's `c₁(Hopf) = 1` stands against `e(TS²) = χ(S²) = 2` from
Gauss–Bonnet: the factor two is `Spin(2) → SO(2)`. So the identification
that §5b listed as "physical and unproved" — Hopf fibre phase = mouth
tangent-plane `Spin(2)` rotation — is standard spin geometry. Given it, the
antilinearity of the mouth holonomy for the fibre complex structure is
canonical.

**T3 — the deck lift.** The antipode `A: x ↦ −x` reverses the orientation
of `S²`, so `dA` sends oriented frames to anti-oriented ones. Exactly the
unit quaternions `u ⊥ i` cover `A` on the spin-frame bundle (tested on
`−j, j, k, (j+k)/√2`; `i` and `(i+j)/√2` do not); the frame image of `L_{−j}`
in the frame coordinates is `diag(1, −1)`, i.e. `dA` composed with one
reflection, orientation restored; and `L_u² = −1`, the central element of
`Spin(3)`, the `2π` rotation. Fibre rotations conjugate the whole circle
`u ⊥ i` into itself: one class. The BAM map `σ = L_{−j}` is the member
`u = −j`; §5b's intrinsic holonomy is the same class. On the two-sheeted
Pin⁻(2) bundle of `S²` (the two orientations), `Ã = (L_u on +, L_{−u} on −)`
satisfies `Ã² = 1` and is free: the deck involution. Its single-sheet
iterate `L_u² = −1` is why the holonomy of the generator squares to the
`2π` rotation.

**T4–T5 — the two Pin⁻ structures and their ABK invariants.** The sign of
`u` on the `+` sheet, `ε = ±1`, is the Pin⁻ structure; neither is preferred
by the metric or by `ι`. The `Z₄`-quadratic enhancement on
`H₁(RP²; Z₂) = {0, g}` is `q(g) = ε mod 4` (quadratic, with `g·g = 1`); the
Gauss sum `1 + i^{ε} = √2 e^{iεπ/4}` gives

```
ABK(RP², ε = +1) = 1,        ABK(RP², ε = −1) = 7        (mod 8)
```

and the `H¹` action shifts `q(g)` by `2`, exchanging them.

**T6 — pairing.** With `Ω₂^{Pin⁻} ≅ Z₈` generated by `RP²` (imported:
Kirby–Taylor; Anderson–Brown–Peterson), a Pin⁻ 3-manifold `W` with
`∂W = RP²_a ⊔ RP²_b` exists iff `ABK_a + ABK_b ≡ 0 (mod 8)`:

| sectors | `ABK` sum | bounds |
|--|--|--|
| `(+, −)`, `(−, +)` | `0` | **yes** |
| `(+, +)` | `2` | no |
| `(−, −)` | `6` | no |
| single `RP²` | `±1` | no |

**The topology does not select a sign. It forces mouths to be created in
opposite Pin⁻ sectors, and forbids a single neck from being created at
all.** The sign is a conserved pair label, not a nuisance parameter. This
is conditional on modelling pair creation as a Pin⁻ bordism of the mouth
surfaces (their Pin⁻ structures induced from one on the worldvolume), which
is a modelling choice.

**Verdict of this section:**
`HOPF_IS_MOUTH_SPIN_BUNDLE_AND_PAIRING_FIXES_OPPOSITE_PIN_SECTORS`,
conditional on (i) the antipodal quotient construction and (ii) the bordism
modelling of pair creation, both stated in the pre-registration before the
computation.

*Dependency ledger:*

```
"Hopf = P_Spin(mouth S²)"  =  ( brane normal ↔ 1 ∈ H [geometric], reference direction i [gauge],
                                orientation class of H [chosen], Ad_{q⁻¹}(i, j, k) [definition] )
deck lift class            =  ( A [derived from ι], u ⊥ i [gauge], sign ε [chosen: Pin⁻ structure] )
opposite sectors           =  ( Ω₂^{Pin⁻} = Z₈ [imported], ABK(RP², ε) = ε [computed],
                                pair creation = Pin⁻ bordism of mouths [chosen] )
```

## 6. The five discrete objects, kept apart

| symbol | what it is | domain | class |
|--|--|--|--|
| `(−1)^ℓ` | antipodal parity of degree-`ℓ` harmonics on `S³` | scalar field on the neck | analytic identity (computed) |
| `η` | line-bundle sign on `M/ι` | field on the quotient | chosen |
| `C` | reversal of the Hopf `U(1)` fibre | `U(1)` bundle | chosen |
| `T` | time-orientation reversal in `(U,V) ↦ (−U,−V)` | spacetime extension of `ι` | chosen (PT versus P) |
| `η_wrap` / orientation | `det m` on the handle, `det dι` on the quotient | manifold | derived from the chosen gluing class |

They live on different objects. Multiplying them into one sign, as the
repository's `wrap_parity × orientation × (−1)^ℓ` bookkeeping does, is
legitimate only after each has been fixed, and the geometry fixes only the
first and (given the gluing class) the last.

## 7. Audit table

| Question | Result | Evidence class |
|--|--|--|
| Is the physical mouth non-orientable? | Two-mouth handle with antipodal gluing: **no** (bulk `w₁ = 0`; the brane `S²`-handle and the brane normal bundle are twisted). Quotient `M/ι`: **yes** (`RP⁴ # RP⁴`). | topological theorem |
| What object carries `w₁ ≠ 0`? | On the handle: the brane `S² ×̃ S¹` and the brane's normal line bundle. On the quotient: the bulk and its brane neck `RP²`. Never the bulk mouth `S³`. | topological theorem |
| Is `J = iσ_y K` a valid lift? | **Yes**: the Spin(4) lift `(−j, 1)` of the rotation `σ = L_{−j}`; and the induced Pin⁻ holonomy of the `RP²` mouth lies in the same fibre-reversing component of `Pin(2) ⊂ SU(2)`, with the same square `−1`. Not a lift of `ι` (`±e_s`) or of `−I₄` (`±` volume). | analytic identity |
| Is it unique? | **No**: determined up to `Spin(2) = U(1)` conjugation (the tangent direction of the generator path) and a sign (the Pin structure). Outcome B. | analytic identity |
| Where does `J² = −1` come from? | From the ambient Pin⁺ structure and the two twisted normal lines `ν = λ ⊕ λ`: `(e_s e_n)² = −1`, and `H̃²` is the spin holonomy `−1` of the `2π` tangent rotation around the neck's great circle. Forced by the quotient, not inserted. | analytic identity; numerically converged (`2.6e-11`) |
| Is PR #129's antipodal BC derived? | **Conditionally**: it is the `η = +1` sector of the unique free involution, given `P_B = −P_A` and the quotient rather than the double cover. Neither is forced. | conditional on an unproved physical identification |
| Does its horizon limit reproduce even/odd BCs? | There is no limit: the labels are lapse-independent and hold at the finite ultrastatic neck; reproduced against the PR #277 oracle to `5e-7`, second order. | numerically converged |
| Is the Hopf `U(1)` the mouth's `Spin(2)`? | **Yes, exactly**: the bulk mouth `S³` with the brane normal as identity is the spin-frame bundle of the brane mouth `S²`; fibre angle `φ` = frame angle `2φ`; `ω_LC = −2A`; `c₁ = 1` against `e = 2`. The `K` of the mouth holonomy is therefore canonical. | analytic identity |
| Do the two mouth sectors pair? | Opposite Pin⁻ sectors bound a worldvolume (`ABK 1 + 7 ≡ 0`); equal sectors and single necks do not. The sign is a conserved pair label. | imported theorem + computed invariants; conditional on the bordism modelling |
| Which inputs remain postulates? | `P_B = −P_A`; quotient versus cover; `η`; `T` (P or PT); the bordism modelling of pair creation; which sector is called the particle. The `U(1)` direction is gauge and the Hopf/`Spin(2)` identification is derived. | definition / chosen |

## 8. Dependency ledger

```
antipodal BC   =  BC( P_B = −P_A [chosen], quotient-not-cover [chosen], η [chosen],
                      ι [derived given P_B = −P_A], (−1)^ℓ [identity], a, R [geometry],
                      P-or-PT [chosen] )
J = iσ_y K     =  J( Hopf complex structure on C² [chosen], C [chosen],
                     −j ∈ SU(2)_L [gauge/convention] )
transport of ι =  ±e_s ( Pin type [chosen], sign [chosen], neck frame [gauge] )
handle topology = mapping torus of m ( m ∈ O(4) [chosen]; (det m₃, ε) [derived from m] )
```

Every headline above with a `chosen` entry in its tree is stated
conditionally.

## 9. Verdict

`FINITE_MOUTH_ADMITS_BUT_DOES_NOT_SELECT_THE_BAM_LIFT`, with the following
sharpening after review (§5b):

* **What is not selected** is the seam gluing of the two-mouth handle, the
  choice of quotient over cover, `η`, the `U(1)` direction of the mouth
  holonomy, and its sign.
* **What is supplied, and was not before.** Given `P_B = −P_A` and the
  quotient, the neck `RP²` has normal bundle `λ ⊕ λ`; the ambient Pin⁺
  structure (forced by `w₁² ≠ 0` and by coherence with the cover's spin
  holonomy) induces the intrinsic Pin⁻ structure; the deck-generator
  holonomy on the restricted spinor bundle is `±ẽ_{t₂}` with `H̃² = −1`,
  the square root of the spin holonomy `−1` of a `2π` tangent rotation on
  the round neck. In the intrinsic module `Cl⁻(2) ≅ H` this is left
  multiplication by a unit quaternion of the `(i, j)` plane — the component
  of `Pin(2) ⊂ SU(2)` that contains `σ = L_{−j}`. **The order-four structure
  `J² = −1` is geometric.** Outcome B of the follow-up: the BAM transport is
  reached up to a `U(1)` and a sign.
* **What the first draft got wrong**, recorded as R1 and R2 in
  `docs/mouth_pin_holonomy_prereg.md`: the Pin types are not mismatched,
  and `J` is the spin lift of a rotation.
* **After the spin-frame round (§5c)** the defensible statement about `J`
  is: *the finite non-orientable mouth geometrically forces the Pin⁻
  order-four transport class; the Hopf fibre is the mouth's `Spin(2)`, so
  the antilinearity is canonical; the two signs are conjugate topological
  sectors that pair creation must produce in opposite pairs.* What remains
  chosen is the antipodal quotient construction, `η`, the time extension,
  and the bordism modelling of creation. It is no longer accurate to say
  that BAM "inserts `iσ_y K`".
* The antipodally glued two-mouth handle remains orientable in the bulk;
  the non-orientable object remains the quotient; PR #129's condition
  remains the `η = +1` scalar sector of `ι`, obtained here without a horizon.

## 10. What the repository should change

* `embedding/transport.py`: corrected in this round (docstrings; the fibre
  check now tests fibre-to-fibre mapping instead of `|w| = 1`).
* README claim rows citing "derived in `embedding/transport.py` without
  ansatz" and "the unique orientation-reversing isometry": to be reworded to
  "`L_{−j}`, base antipode × fibre reversal, orientation-preserving".
* `embedding/topology.py`: the Boolean `wrap_parity`/`orientation` fields
  should be replaced by, or documented as shorthand for, the gluing class
  `(det m₃, ε)` and the quotient sign `η`.
* PR #129's docstring: "fixed by the BAM antipodal postulate" is correct and
  should stay; it may now cite `ι` as the finite-mouth realisation of that
  postulate, with the two choices it depends on.
