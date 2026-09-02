# Pre-registration: does the finite mouth force the transport and the antipodal boundary condition?

**Status: frozen before any code exists.** This document is committed on its
own, ahead of the module, solver, probe and tests it constrains. Its success
criteria may be changed afterwards only by an explicit correction note that
says what the original criterion got wrong.

The question, from `docs/qft_emergence_audit.md` (revision `470a8b6`):

```
finite physical mouth  →  orientation/bundle transport  →  antipodal field boundary condition
```

The repository has the three boxes. Does the geometry of the first force the
other two? Allowed answers: yes, no, underdetermined.

---

## 0. What is taken as given (and its evidence class)

| input | source | class |
|--|--|--|
| bulk spatial slice `Σ = S⁴_R` with the observed `S³` its totally geodesic equator | PR #277 | chosen (physical identification) |
| two geodesic 4-balls of angular radius `a < π/2` excised at `P_A, P_B` | PR #277 | chosen; **the centre placement is not specified anywhere in PR #277** and is recorded here as a free input |
| tube `T = [−S,S] × S³`, `dΣ² = ds² + f²dΩ₃²`, `f = √(s²+b²)`, `b = R sin²a`, `S = R sin a cos a` | PR #277, derived from `⁴R = 0` + Darmois | analytic identity |
| closed-form static admittance `Y_ℓ(0)` of the tube | PR #277 (`static_admittance`) | analytic identity, regression-checked |
| the PR #129 condition `Φ(U,V,Ω) = Φ(−U,−V,Ω̄)` with even-ℓ Neumann / odd-ℓ Dirichlet | PR #128/#129 | imposed postulate on the Tangherlini horizon |
| `J = iσ_y K` with `J² = −1`, `h(Jz) = −h(z)`, `J(e^{iφ}z) = e^{−iφ}Jz`, `det_R J = +1` | `embedding/transport.py`, corrected by the audit | analytic identity; the "orientation-reversing on S³" claim is **false** and is not used |

Nothing below may use a singlet, a Born rule, a projector, a tensor product,
CHSH, or a QED amplitude.

---

## 1. Hypotheses, with what would falsify each

### H1 — the seam gluing is unconstrained by the geometry

The Darmois conditions at a mouth are the two scalars `f_m = R sin a` and
`f′ = cos a`, both `O(4)`-invariant on the mouth `S³`. So each seam map
`ψ_A, ψ_B ∈ O(4)` is free, and only the loop monodromy
`m = ψ_B⁻¹ψ_A` (up to conjugation) is meaningful.

*Predictions.*

* **T1.** The closed 4-manifold `M = C ∪ T` is, up to diffeomorphism, the
  mapping torus of `m`, and its diffeomorphism type is fixed by
  `det m ∈ {±1}`: `S³ × S¹` (orientable) or the non-orientable
  `S³ ×̃ S¹`. Its `w₁` evaluated on the handle loop is `det m`.
* **T2.** The antipodal map `−I₄` on the mouth `S³` has `det = +1`.
  Therefore an antipodally glued two-mouth handle is **orientable** in the
  bulk. *Falsifier:* `w₁(M_{−I}) ≠ 0`.
* **T3.** With the brane (equatorial `S³`) preserved, `m = (m₃, ε) ∈ O(3) × O(1)`
  block-diagonally, and there are exactly four topological classes
  `(det m₃, ε)`. Then
  `w₁(bulk handle) = det m₃ · ε`,
  `w₁(brane handle S² × S¹ or S² ×̃ S¹) = det m₃`,
  `w₁(normal line bundle of the brane) = ε`.
  The antipodal class is `(−, −)`: brane handle non-orientable, normal
  bundle twisted, bulk orientable.
* **T4.** Within `SO(4)`, `−I₄` is isotopic to `I`, so on the two-mouth handle
  "antipodal versus identity gluing" is a *geometric* (twist-modulus)
  distinction, not a topological one; its only spectral signature is the
  sign `(−1)^ℓ` on the off-diagonal tube coupling. *Falsifier:* the two
  gluings give the same closed-manifold static spectrum at odd `ℓ`.

*Evidence class of H1 if it holds:* topological theorem (T1–T3), analytic
identity (T4).

### H2 — there is exactly one free involution, and it is the antipode

If, and only if, `P_B = −P_A` on `S⁴`, the `S⁴` antipode `x ↦ −x` preserves
`C`, swaps the two mouths, and extends through the tube as
`ι(s, Ω) = (−s, −Ω)`.

*Predictions.*

* **I1.** `ι` is an isometry of `M` for both lapses (`N = 1` and
  `N = |s|/f`), is free, and is orientation-reversing: `det dι = −1`.
* **I2.** `ι` is the **unique** free isometric involution of `M` that swaps
  the mouths: the only fixed-point-free involution in `O(5)` is `−I₅`, and
  the only extension through the tube that is an isometry and matches
  `−I₄` at both seams is `(s, Ω) ↦ (−s, −Ω)`. *Falsifier:* exhibit another.
* **I3.** `M/ι` is non-orientable; its brane slice is the orientable `RP³`
  geon of PR #169/#196 (`det = (−1)(−1) = +1` there); its neck cross-section
  is `S³/± = RP³` (orientable) and the brane neck cross-section is
  `S²/± = RP²` (non-orientable). This is the object the repository's
  phrases "`RP²` mouth", "non-orientable throat" and "Pin⁻" refer to, and
  it lives on the **quotient**, not on the two-mouth handle.
* **I4.** The PR #129 identification is precisely `ι`-invariance:
  `s ↦ −s` at fixed `t` is `(U, V) ↦ (−U, −V)` in the Tangherlini Kruskal
  chart, and `Ω ↦ −Ω` is `Ω ↦ Ω̄`. It is available for the ultrastatic
  lapse too: **no horizon limit is needed**.

*If H2 holds it is a topological theorem, conditional on the chosen centre
placement `P_B = −P_A` (chosen) and on taking the quotient rather than the
double cover (chosen).*

### H3 — the scalar boundary condition on the quotient

A scalar on `M/ι` is an `ι`-invariant scalar on `M`, or a section of the
orientation line bundle: `Φ(ιx) = η Φ(x)`, `η = ±1`. With
`Φ = Σ ψ_ℓ(s) Y_ℓ(Ω)` and `Y_ℓ(−Ω) = (−1)^ℓ Y_ℓ(Ω)`:

```
ψ_ℓ(−s) = η (−1)^ℓ ψ_ℓ(s)
```

*Predictions.*

* **B1.** `η = +1`: even `ℓ` ⇒ `ψ′_ℓ(0) = 0` (Neumann), odd `ℓ` ⇒ `ψ_ℓ(0) = 0`
  (Dirichlet). This is PR #129's structure, obtained at the finite neck with
  `N = 1`. `η = −1` swaps the two. The sign `η` is a discrete choice the
  geometry does not fix.
* **B2.** The half-tube DtN admittance (one mouth, neck condition) equals the
  `(1, η(−1)^ℓ)` eigenvalue of the PR #277 two-mouth oracle:

  ```
  Y_ℓ^N = 2π²F² [ k tanh(kX) − cos a ]        (Neumann at the neck)
  Y_ℓ^D = 2π²F² [ k coth(kX) − cos a ]        (Dirichlet at the neck)
  ```

  with `k = ℓ+1`, `F = R sin a`, `X = arcosh(1/sin a)`. Numeric targets at
  `R = 1, a = 0.3`, to be matched by an independent tridiagonal solve to
  better than `1e-5` relative and with second-order convergence
  (successive-refinement ratios `≥ 3.5`):

  | `ℓ` | `Y_ℓ^N` | `Y_ℓ^D` |
  |--|--|--|
  | 0 | `0.000000000` | `0.157587622` |
  | 1 | `1.797266559` | `1.804461992` |
  | 2 | `3.524607516` | `3.524854051` |
  | 3 | `5.248595411` | `5.248602920` |

  `Y_0^N = 0` exactly (the constant is the Neumann monopole; `tanh X = cos a`)
  and `Y_0^D = 2G` with `G = π²R²sin⁴a/cos a = 0.078793811`.
* **B3.** *Antipodal control.* Replace `Ω ↦ −Ω` by the identity
  (`ι₀(s,Ω) = (−s, Ω)`): the involution is no longer free (it fixes the whole
  neck), the quotient is a manifold with boundary, and the condition becomes
  `ψ′ = 0` for **all** `ℓ` (`η = +1`) or `ψ = 0` for all `ℓ` (`η = −1`). The
  `(−1)^ℓ` grading must disappear. *Falsifier:* it survives.
* **B4.** *Geometry control.* Across `a ∈ {0.05, 0.3, 0.8, 1.2}` the sector
  labels (which `ℓ` are Neumann) are invariant; the admittances vary.

*Evidence class if B1 holds:* analytic identity, **conditional** on H2's two
choices and on `η`.

### H4 — the horizon limit is not a limit

The Tangherlini branch is the same spatial handle with `N = |s|/f`. The neck
is the bifurcation surface; the condition `ψ_ℓ(−s) = η(−1)^ℓ ψ_ℓ(s)` is a
statement about the spatial profile and does not see the lapse.

* **L1.** The neck sector labels are identical for both lapses.
* **L2.** In the Tangherlini lapse, `ι` maps exterior I to exterior III and
  reverses the time orientation of the bifurcate Killing horizon
  (`U + V ↦ −(U + V)`), so PR #129's map contains a **time reversal** that
  the spatial `ι` alone does not. This is recorded as a distinct datum `T`.

### H5 — the spinor lift, and where `J` actually lives

* **S1.** In the neck frame `(∂_s, e_i)` with `e_i` left-invariant on
  `S³ = SU(2)`, `dι = diag(−1, +1, +1, +1)`: a single reflection. Its lifts
  to `Pin^±(4)` are `±γ_s`. There are two lifts per Pin type; they are
  complex-**linear** on the Dirac spinor and exchange the two Weyl
  chiralities of `Spin(4)`.
* **S2.** `J = iσ_y K` is antilinear on `C²` with the complex structure in
  which the Hopf fibre acts as scalar multiplication (`L_i`, left
  multiplication by `i` on `H`), and it is the complex-**linear** matrix
  `iσ_y ∈ SU(2)` with respect to the other complex structure (`R_i`). Its
  real `4×4` matrix anticommutes with `L_i` and commutes with `R_i`.
  *Falsifier:* either commutation fails.
* **S3.** Therefore `J` is not a lift of `ι`, of `m = −I₄`, or of any
  element of `O(4)`: no Pin lift is antilinear. `J` factorises as
  (spin rotation by `π` about `y`) ∘ (reversal of the fibre `U(1)`). The
  second factor is charge conjugation `C` on the `U(1)` bundle — a separate
  `Z₂` datum on the bundle, not on the manifold — and the first is a
  geometric rotation `−j ∈ SU(2)_L` whose lift to `Spin(4)` is `(−j, 1)`,
  not a reflection. The central falsifier of the task is thereby answered:
  **`J` cannot be obtained from the mouth gluing without inserting, by
  hand, (i) a two-component complex structure and (ii) the antilinear
  `C`.** If a construction produces `iσ_y K` from `ι` alone, S3 is false.
* **S4.** *Lift control.* The lifts of `ι` (`±γ_s` in `Pin⁺` and `Pin⁻`) and
  the lifts of `−I₄` (`±(1, −1) ∈ SU(2) × SU(2)`) are enumerated; they are
  pairwise distinct operators on the Dirac spinor, so "the" transport is
  not unique.

Nothing in H5 uses Pauli matrices as a starting point: the representation is
built from the Clifford algebra of the neck frame.

### Bookkeeping of the five discrete objects

`(−1)^ℓ` is the antipodal parity of degree-`ℓ` harmonics on `S³` (analytic
identity). `η` is the line-bundle sign on the quotient (chosen). `C` is the
`U(1)` fibre reversal (chosen, bundle datum). `T` is the time reversal
contained in `(U,V) ↦ (−U,−V)` (part of PR #129's map, absent from the spatial
`ι`). `η_wrap`/orientation is `det m` or `det dι` (topological, from the
chosen gluing class). They are not multiplied together here.

---

## 2. Verdict rule, fixed in advance

* `FINITE_MOUTH_FORCES_TRANSPORT_AND_ANTIPODAL_BC` — only if the geometry
  fixes the gluing class **and** forces the quotient **and** forces
  `η` **and** `J` is the unique lift of `ι`. Expected: **no**, because H1
  predicts the gluing is free.
* `FINITE_MOUTH_ADMITS_BUT_DOES_NOT_SELECT_THE_BAM_LIFT` — if `ι` exists and
  is unique given the centre placement (H2), the PR #129 condition is its
  `η = +1` sector (H3), and `J` is compatible with but not produced by the
  gluing (H5). Expected outcome.
* `FINITE_MOUTH_INCOMPATIBLE_WITH_CURRENT_BAM_TRANSPORT` — if `M` admits no
  free involution swapping the mouths, or the quotient sector fails to
  reproduce B2's oracle, or `J` is inconsistent with every allowed lift.

A narrower verdict is permitted if the trichotomy is shown to be mis-posed;
that must be stated with the reason.

---

## 3. Required controls (each must be able to fail)

| control | what changes | what must happen |
|--|--|--|
| topology | identity gluing → reflection gluing (`ε = −1`) | `w₁` on the loop flips `+1 → −1` |
| antipodal | `Ω ↦ −Ω` → `Ω ↦ Ω` in the involution | involution stops being free; `(−1)^ℓ` grading disappears |
| lift | enumerate `±γ_s` in `Pin^±`, `±(1,−1)` | distinct operators; uniqueness false |
| geometry | `a ∈ {0.05, 0.3, 0.8, 1.2}` | sector labels invariant, admittances vary |
| convergence | steps `500 → 1000 → 2000 → 4000` | relative error ratios `≥ 3.5` |

---

## 4. Dependency ledger to be printed by the probe

```
antipodal BC  =  BC( P_B = −P_A [chosen], quotient-not-cover [chosen], η [chosen],
                     ι [derived given the first], (−1)^ℓ [identity], a, R [geometry] )
J = iσ_y K    =  J( C² complex structure [chosen], C [chosen], −j ∈ SU(2)_L [gauge/convention] )
transport of ι on spinors  =  ±γ_s  ( Pin type [chosen], sign [chosen], frame [gauge] )
```

Any headline whose tree contains a `chosen` load-bearing entry is to be
stated conditionally.
