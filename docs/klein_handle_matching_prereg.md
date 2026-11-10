# Pre-registration: the matching map from the recorded BAM geometry to the seam cells

Publish this file before implementing or running anything. Baseline: branch
head `562af71` (round 14). Seed `2026111004`.

**Prior state.** Nothing below has been computed. The predictions are
analytic, derived in this file from the objects recorded in
`docs/finite_mouth_bundle_transport.md`. As in rounds 12–14, freeze and
implementation are written in one session by one author; the only ordering
evidence is that this commit precedes the code.

## Question

Rounds 13–14 found that the pairing between the two quadrature sectors
exists only on a seam that reverses both time and the fibre of a charged
field, and left `bam_cell = NOT_ESTABLISHED`. The repository records two
candidate inputs:

* PR #129's antipodal condition, the PT-type extension of
  `ι(s, Ω) = (−s, −Ω)`: **time reversal**;
* PR #279's mouth transport, `σ = L_{−j}`, whose spin lift is
  `J = iσ_y K`: **fibre reversal**.

Do they act on a fibre-charged field as the doubly reversed seam does? And
is the sign of round 14's trapped flux the Pin sector label of PR #279, as I
guessed at the end of round 14?

## Objects, as recorded

Bulk mouth `S³ = {q ∈ H, |q| = 1}`, `q = z₁ + z₂ j`. Hopf fibre
`q ↦ e^{iφ} q`, i.e. `(z₁, z₂) ↦ e^{iφ}(z₁, z₂)`; Hopf base map
`h(q) = q⁻¹ i q` (`mouth_spin_frame.spin_frame`). A **fibre-charged
function** of charge `k` satisfies `f(e^{iφ} q) = e^{ikφ} f(q)`; these are
the polynomials of bidegree `(p, r)` in `(z, z̄)` with `p − r = k`. The `k`-th
fibre harmonic of a real function on `S³` is such an `f`, with the `−k`
harmonic its conjugate.

* `−I₄`: the angular part of `ι`.
* `σ(z₁, z₂) = (z̄₂, −z̄₁)` (`embedding.transport.orientation_reversal_on_s3`).
* The spatial brane mouth: the equatorial `S² = {Re q = 0}`, the brane normal
  being `1 ∈ H` (`finite_mouth_bundle_transport.md` §1, §5c).

## Analytic predictions

### M1. The antipodal map is a fibre translation

`−q = e^{iπ} q`. So `−I₄` commutes with the fibre rotation, fixes every point
of the Hopf base, and acts on charge-`k` functions as the scalar `(−1)^k`,
complex-linearly. Predicted for all bidegrees with `p + r ≤ 4`: pullback
equals `(−1)^k` to `1e-10`; `h(−q) = h(q)` to `1e-12`.

### M2. The spatial brane mouth is not the Hopf base

On the equator `E = {Re q = 0}`:

* `−I₄` preserves `E` and acts on it as the antipode;
* `h` restricted to `E` is two-to-one, `h(p) = h(−p)`, with
  `h(p) = 2(p·i) p − i`;
* `σ` does **not** preserve `E`: predicted `max |Re(σp)| > 0.5` over random
  `p ∈ E`;
* `σ` covers the antipode of the Hopf base, `h(σq) = −h(q)`, and reverses
  the fibre, `σ(e^{iφ} q) = e^{−iφ} σ(q)`.

So the antipodal map is the antipode of the spatial mouth and the identity
on the Hopf base, while the transport is the antipode of the Hopf base and
does not map the spatial brane to itself. The two recorded maps are attached
to two different 2-spheres.

### M3. The transport on fibre-charged functions

Pullback by `σ` sends charge `k` to charge `−k`. With `c` complex
conjugation, `A = c ∘ σ*` maps each bidegree space to itself, antilinearly,
and

    A² = (−1)^k,        σ² = −I₄.

On the charge-1 doublet `span{z₁, z₂}`, `A(a z₁ + b z₂) = −b̄ z₁ + ā z₂`:
the repository's `J = iσ_y K` up to the sign that is the Pin structure.
Predicted to `1e-10` for all bidegrees with `p + r ≤ 4`.

### M4. PR #129's map on a bulk scalar is the fibre-preserving cell

By M1, the PT-type condition `Φ(−t, ι x) = η Φ(t, x)` on the fibre-charged
components of a bulk scalar reads

    Ψ_k(−t, −s, q) = η (−1)^k Ψ_k(t, s, q):

complex-linear and time-reversing, round 13's `(reversed, kept)` cell with
twist `η(−1)^k`. By rounds 13–14 the charge then vanishes identically on
admissible fields, there is no pairing, and flux through the handle is
inadmissible. Checked on the round-13 lattice for both twists: total charge
on the fixed space zero to `1e-12`.

**As recorded, PR #129's condition does not give the doubly reversed seam.**

### M5. The transport seam is consistent only with the matching loop holonomy

Take the handle traversal to act on the field by an internal map `B` with
`B² = ±1`, and let `W = ±1` be the holonomy of the field around the doubled
handle loop (`W = −1`: antiperiodic). On the lattice cover with one
sign-twisted seam edge the spatial deck map is `P_ε = P·diag(ε)`, `ε = ±1` on
the two sheets, with `P_ε² = W`. The deck map `D = T ∘ (B ⊗ P_ε)` satisfies

    D² = B² · W.

Predicted on the lattice, with `H` invariant under `D` in every case:

| field | `B` | `W = +1` | `W = −1` |
|---|---|---|---|
| singlet | conjugation, `B² = +1` | involution, fixed space half | `D² = −1`, fixed space `{0}` |
| doublet | `J = iσ_y K`, `B² = −1` | `D² = −1`, fixed space `{0}` | involution, fixed space half |

`σ² = −I₄` is the `2π` frame rotation, whose spin lift is the `−1` that
`finite_mouth_bundle_transport.md` §5b computes as the holonomy of the neck's
great circle. On the lattice that `−1` is `W`. The spinor doublet is
admissible on the transport seam exactly because the two signs cancel.

### M6. On the transport seam with time reversal, charge survives and spin vanishes

Doublet, `B = J`, `W = −1`, `s_t = −1`, mouths responding. For a generator
`X` of a symmetry of `H`, the Noether form obeys `N_X ∘ D = s_t κ_X N_X`,
`κ_X = ±1` as `X` commutes or anticommutes with `D`. `J` anticommutes with
`i` and commutes with `su(2)`. Predicted on the fixed space, both sizes:

| seam | time | charge (`X = i`) | spin (`X = iσ_a`) |
|---|---|---|---|
| `J`, `W = −1` | reversed | nondegenerate, balanced signature | **zero** |
| `J`, `W = −1` | kept | **zero** | nondegenerate |
| linear, `W = +1` | reversed | **zero** | **zero** |
| linear, `W = +1` | kept | nondegenerate | nondegenerate |

Zero means `≤ 1e-12`; balanced signature at the primary size is `(64, 64)`
for the charge. On the first row: the fixed space is a Lagrangian half,
responding histories satisfy `z(τ) = D z(−τ)` to `1e-9` at amplitudes
`1/64 … 1/8`, charge and all three spin components are conserved to `1e-9`,
and the total spin of admissible histories is zero.

So the doubly reversed structure of rounds 13–14 is reached by the recorded
transport, for the spinor doublet, **if** the field is a spin-frame field
whose complex structure is the fibre rotation and **if** the time extension
is PT.

### M7. The trapped-flux sign is not the Pin label

In the antilinear cells the Pin sign is the sign of the antilinear deck map
(`σ` against `σ⁻¹`, `J` against `−J`). In round 14's model:

* the quarter rotation `R: Ψ → iΨ` is a symmetry of `H`, leaves the gauge
  field and the handle flux untouched, and conjugates the deck map to the
  one of opposite sign;
* charge conjugation `C` is a symmetry of `H`, flips the handle flux, and
  commutes with the deck map.

Predicted: both identities to `1e-12`, and both symmetries map admissible
histories to admissible histories of the transformed labels (`1e-9`). The
two signs are flipped by different symmetries, so they are **independent**
labels, and for a single handle neither is an invariant. My round-14 guess
that the flux sign "should be tied to the Pin sector label" is predicted to
be **wrong**.

## Verdict fields

    numerical_validity         PASS | UNRESOLVED_NUMERICALLY | INCOMPLETE
    antipodal_map              FIBRE_TRANSLATION_FIXING_THE_HOPF_BASE
    two_spheres                SPATIAL_BRANE_MOUTH_IS_NOT_THE_HOPF_BASE
    transport_on_fibre_modes   ANTILINEAR_SQUARE_MINUS_ONE_ON_ODD_CHARGE_J_ON_THE_DOUBLET
    pr129_cell                 FIBRE_PRESERVING_NO_PAIRING
    spinor_seam_consistency    TRANSPORT_SEAM_NEEDS_ANTIPERIODIC_DOUBLED_LOOP
    transport_seam_charges     CHARGE_SURVIVES_SPIN_VANISHES_WITH_TIME_REVERSAL
    flux_sign_vs_pin_label     INDEPENDENT_LABELS
    bam_cell                   DEPENDS_ON_FIELD_TYPE_AND_TIME_EXTENSION
    support_dynamics           EFFECTIVE_SUPPORT_ASSUMED
    positive_response          NOT_DERIVED
    event_law                  NOT_DERIVED

Status vocabulary as in rounds 12–14. A prediction contradicted by a
resolved measurement is reported as wrong in the first paragraph of the
write-up.

## What may not be concluded

* **Which field BAM's matter is.** M4 and M6 give different cells for a bulk
  scalar and for a spin-frame field. Nothing here decides between them, and
  `P` against `PT` remains the chosen datum `T`.
* **A physical pairing.** M3's identification of `A` with `J` on the doublet
  is a representation identity. In the boundary-pairing protocol's
  vocabulary the spinor matching stays `NOT_DERIVED`: no action-derived
  boundary map produces it.
* The lattice doublet keeps only the lowest fibre harmonic with scalar
  hopping. It is not a monopole-harmonic field theory on the handle, and no
  gauge field is coupled to the doublet here.
* The identification of the lattice sign `W` with the spin holonomy of
  §5b is by matching `σ² = −I₄`; the neck geometry is not rebuilt.
* No norm, detector law or probability follows from any of this.
