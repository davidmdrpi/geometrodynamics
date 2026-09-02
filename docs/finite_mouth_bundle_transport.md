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

*A Pin mismatch.* With `w(RPⁿ) = (1+a)^{n+1}`: `RP²` has `w₁ = a`, `w₂ = a²`,
so it is Pin⁻ only; `RP⁴` has `w₁ = a`, `w₂ = 0`, `w₁² = a² ≠ 0`, so it is
Pin⁺ only, and so is `RP⁴ # RP⁴` (`w₂ = 0`, `w₁² = a₁² + a₂² ≠ 0`). The Pin⁻
the repository assigns to the `RP²` mouth is the structure induced on `RP²` by
the spin structure of the orientable brane slice; it is not the Pin structure
of the bulk quotient, which is of the opposite type
(`pin_structures_rp`). *Evidence class:* topological theorem (Stiefel–Whitney
classes of projective spaces), conditional on `P_B = −P_A`.

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

Every lift is a real matrix in the regular representation, hence
complex-linear in any complex structure. `J` is antilinear in the Hopf
structure. **No lift of any element of `O(4)` is `J`.** The lift of `ι` is a
chirality-exchanging reflection; the lift of `−I₄` is a chirality sign; `J`
is a rotation with an antiunitary factor bolted on. The central falsifier of
the round is therefore answered: `iσ_y K` cannot be obtained from the mouth
gluing without inserting (i) a two-component complex structure on the state
space — the Hopf one, in which the fibre is the phase — and (ii) the
antilinear `C`. Both are inserted by hand in `embedding/transport.py` and
`bell/bulk_identity.py`. The derivation there is circular in exactly the
sense the pre-registration anticipated.

**Lift control.** The four lifts of `ι` are pairwise distinct within a Pin
type and differ in square across types; the two lifts of `−I₄` differ by
sign. Uniqueness of "the" transport is false even before `J` is considered.

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
| Is `J = iσ_y K` a valid lift? | **No**, of any gluing map: all Pin lifts are complex-linear; `J` is antilinear for the Hopf structure. `J = (spin rotation by π) ∘ (U(1) reversal)`. | analytic identity |
| Is it unique? | Not as a lift. The lifts of `ι` are four, of `−I₄` two; among Hopf-fibre-reversing isometries `J` is one point of a two-parameter component. | analytic identity |
| Is PR #129's antipodal BC derived? | **Conditionally**: it is the `η = +1` sector of the unique free involution, given `P_B = −P_A` and the quotient rather than the double cover. Neither is forced. | conditional on an unproved physical identification |
| Does its horizon limit reproduce even/odd BCs? | There is no limit: the labels are lapse-independent and hold at the finite ultrastatic neck; reproduced against the PR #277 oracle to `5e-7`, second order. | numerically converged |
| Which inputs remain postulates? | `P_B = −P_A`; quotient versus cover; `η`; `T` (P or PT); the Hopf complex structure on the state space; `C`; the Pin type and sign. | definition / chosen |

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

`FINITE_MOUTH_ADMITS_BUT_DOES_NOT_SELECT_THE_BAM_LIFT`

The trichotomy was well posed. Its middle option needs one sharpening: the
geometry does not merely fail to select the BAM lift, it contradicts two of
the words attached to it. The antipodally glued two-mouth handle is orientable
in the bulk, and the non-orientable object — the quotient — is Pin⁺, not the
Pin⁻ assigned to its `RP²` neck. `J` is a rotation whose `K` is a choice of
complex structure. What the finite mouth *does* supply, and supplies
uniquely once `P_B = −P_A` is chosen, is the involution `ι`; PR #129's
boundary condition is its `η = +1` scalar sector, and this document is the
first place that sector is obtained without a horizon.

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
