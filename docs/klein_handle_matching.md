# The matching map: the two recorded reversals live on different spheres

**All seven frozen predictions hold, 30/30 required checks, with nothing
retuned — and the answer to the question as I posed it is no.** The
repository's time-reversing antipodal map (PR #129) and its fibre-reversing
mouth transport (PR #279) do not act together on one fibre-charged field as
the doubly reversed seam of rounds 13–14.

* The antipodal map is a **Hopf fibre translation**. On the fibre-charged
  components of a bulk scalar it is complex-linear, so PR #129's condition is
  the *fibre-preserving* cell. By rounds 13–14 that cell has no charge, no
  pairing and no admissible handle flux.
* The fibre reversal belongs to a **different map on a different sphere**.
  The antipodal map is the antipode of the spatial brane mouth and the
  identity on the Hopf base; the transport is the antipode of the Hopf base
  and does not map the spatial brane mouth to itself.
* The doubly reversed structure **is** reached by the recorded transport,
  for the spinor doublet, under two conditions the repository records as
  choices: the field is a spin-frame field whose complex structure is the
  fibre rotation, and the time extension is PT. It then works only because
  two minus signs cancel.
* My guess at the end of round 14, that the sign of the trapped flux is the
  Pin sector label, is **wrong**. The two signs are independent.

So `bam_cell` moves from `NOT_ESTABLISHED` to
`DEPENDS_ON_FIELD_TYPE_AND_TIME_EXTENSION`: the dependence is now explicit,
and it is not resolved.

Part I is representation theory on `S³` and Part III is symmetry algebra; the
numerics confirm derivations. Only Part II integrates anything.

*Pre-registered in [klein_handle_matching_prereg.md](klein_handle_matching_prereg.md)
(`3d3ee38`), committed before the module and unchanged. Module
`geometrodynamics/bulk/klein_handle_matching.py`; probe
`klein_handle_matching_probe.py`; tests `tests/test_klein_handle_matching.py`
(24, one added after the round-14 review). Archive `experiments/closure_ledger/runs/20261110_klein_handle_matching/`.*

```bash
python -m experiments.closure_ledger.klein_handle_matching_probe   # 30/30
python -m pytest -q tests/test_klein_handle_matching.py             # 24 tests
```

| verdict field | value |
|---|---|
| numerical_validity | `PASS` |
| antipodal_map | `FIBRE_TRANSLATION_FIXING_THE_HOPF_BASE` |
| two_spheres | `SPATIAL_BRANE_MOUTH_IS_NOT_THE_HOPF_BASE` |
| transport_on_fibre_modes | `ANTILINEAR_SQUARE_MINUS_ONE_ON_ODD_CHARGE_J_ON_THE_DOUBLET` |
| pr129_cell | `FIBRE_PRESERVING_NO_PAIRING` |
| spinor_seam_consistency | `TRANSPORT_SEAM_NEEDS_ANTIPERIODIC_DOUBLED_LOOP` |
| transport_seam_charges | `CHARGE_SURVIVES_SPIN_VANISHES_WITH_TIME_REVERSAL` |
| flux_sign_vs_pin_label | `INDEPENDENT_LABELS` |
| bam_cell | `DEPENDS_ON_FIELD_TYPE_AND_TIME_EXTENSION` |
| positive_response, event_law | `NOT_DERIVED` |

## What was matched

The objects are the ones `docs/finite_mouth_bundle_transport.md` records, used
through the repository's own functions:

* bulk mouth `S³ = {q ∈ H, |q| = 1}`, `q = z₁ + z₂ j`; Hopf fibre
  `q ↦ e^{iφ} q`; Hopf base map `h(q) = q⁻¹ i q`
  (`mouth_spin_frame.spin_frame`);
* `−I₄`, the angular part of `ι(s, Ω) = (−s, −Ω)`;
* `σ(z₁, z₂) = (z̄₂, −z̄₁)` (`embedding.transport.orientation_reversal_on_s3`),
  whose spin lift is `J = iσ_y K`;
* the spatial brane mouth, the equatorial `S² = {Re q = 0}`.

A **fibre-charged function** of charge `k` satisfies
`f(e^{iφ} q) = e^{ikφ} f(q)`: the polynomials of bidegree `(p, r)` in
`(z, z̄)` with `p − r = k`. No bundle is postulated. The `k`-th fibre harmonic
of any function on `S³` is such an `f`, and a map of `S³` acts on it by
pullback. Whether that action is linear or antilinear is then a computation.

## Part I. The two recorded maps

| | antipode `−I₄` | transport `σ` |
|---|---|---|
| on the fibre | translation by `π` | reversal |
| on the Hopf base | identity | antipode |
| on the spatial brane mouth | antipode | does not preserve it |
| on charge-`k` functions | the scalar `(−1)^k`, linear | antilinear, square `(−1)^k` |

**M1.** `−q = e^{iπ} q` (`1.6e-16`). The antipode commutes with the fibre
rotation, fixes every point of the Hopf base, and acts on all 15 bidegrees
with `p + r ≤ 4` as `(−1)^k` (`4.1e-15`).

**M2.** On the spatial brane mouth the antipode is the antipode, and the
Hopf map is two-to-one there, `h(p) = h(−p) = 2(p·i)p − i` (`5.6e-16`). So
the same map is the antipode of one sphere and the identity on the other.
The transport covers the antipode of the Hopf base (`3.3e-16`) and reverses
the fibre, but moves the spatial brane mouth off itself: `Re(σp)` reaches
`0.998` on it. It agrees with the repository's recorded matrix `L_{−j}`,
`det = +1`, `σ² = −I₄`.

**M3.** Pullback by `σ` sends charge `k` to `−k` (the residual against the
same charge is `0.99`). With `c` complex conjugation, `A = c ∘ σ*` preserves
each bidegree, is antilinear, and `A² = (−1)^k` (`4.7e-15`). On the charge-1
doublet in the basis `(z₁, z₂)` its matrix is

    [[0, −1], [1, 0]],        A(a z₁ + b z₂) = −b̄ z₁ + ā z₂,

which is `−iσ_y K`: the repository's `J` up to the sign that is the Pin
structure. So `J` is not only the spin lift of `σ`; it is literally how `σ`
acts on the lowest fibre-charged functions.

**M4.** By M1, PR #129's PT-type condition on the fibre-charged components
of a bulk scalar is

    Ψ_k(−t, −s, q) = η (−1)^k Ψ_k(t, s, q),

complex-linear and time-reversing. On the round-13 lattice, for both twists,
the charge is deck-odd and vanishes on the fixed space (`1.1e-15`).

## Part II. When the transport seam admits a field

`σ` has order four, so it cannot be an involutive deck map by itself. On the
handle the question is what happens after two traversals. Let the traversal
act on the field by `B`, `B² = ±1`, and let `W = ±1` be the field's holonomy
around the doubled handle loop. On the lattice cover with one sign-twisted
seam edge, the deck map `D = T ∘ (B ⊗ P_ε)` has `D² = B²·W`.

**M5.** With `H` invariant under `D` in every case:

| field | `B` | `W` | `D²` | fixed space (primary) |
|---|---|---:|---:|---:|
| singlet | conjugation | `+1` | `+1` | 64 of 128 |
| singlet | conjugation | `−1` | `−1` | 0 |
| doublet | `J` | `+1` | `−1` | 0 |
| doublet | `J` | `−1` | `+1` | 128 of 256 |

The spinor doublet is admissible on the transport seam **only** on the
antiperiodic loop. `σ² = −I₄` is the `2π` frame rotation, and
`finite_mouth_bundle_transport.md` §5b computes its spin lift, the holonomy
of the neck's great circle, as `−1`. On the lattice that `−1` is `W`. The
`J² = −1` that the repository has long treated as the signature of spin-½ is
here a consistency condition: it must be cancelled by the spin holonomy for
any field to exist on the seam at all.

**M6.** Which conserved quantities survive on the fixed space follows from
one rule, `N_X ∘ D = s_t κ_X N_X`, with `κ_X = ±1` as the generator commutes
or anticommutes with `D`. `J` anticommutes with `i` and commutes with
`su(2)`:

| seam | time | charge | spin |
|---|---|---|---|
| `J`, `W = −1` | reversed | `(64, 64)` | zero |
| `J`, `W = −1` | kept | zero | `(64, 64)` |
| linear, `W = +1` | reversed | zero | zero |
| linear, `W = +1` | kept | `(64, 64)` | `(64, 64)` |

Zero means `≤ 6.4e-15`. Time reversal keeps the charge and removes the total
spin; without it the reverse. On the first row, with the mouths responding:
the fixed space is a Lagrangian half (130 of 260), histories satisfy
`z(τ) = D z(−τ)` to `2.4e-15`, charge and the three spin components are
conserved to `1.5e-12`, total spin stays zero (`9.9e-15`), and the sheet
charge varies while obeying `Q_1(τ) + Q_1(−τ) = Q`.

So the doubly reversed structure of rounds 13–14, with the charge as the
pairing of the two sectors, is realised by the recorded transport for the
spinor doublet. "Total spin zero" is a statement about the cover: the spin
on one sheet is free and even in time.

## Part III. The flux sign is not the Pin label

In the antilinear cells the Pin sign is the sign of the antilinear deck map,
`σ` against `σ⁻¹`, `J` against `−J`. In round 14's model:

| symmetry of `H` | Pin sign | handle flux |
|---|---|---|
| quarter rotation `Ψ → iΨ` | flipped | unchanged |
| charge conjugation | unchanged | flipped |

Both identities hold exactly. Both symmetries carry admissible histories to
admissible histories of the transformed labels (`4.7e-16`), with the flux
history unchanged in one case and negated in the other. As a control, the
rotated history tested against the *original* Pin sign misses by `0.30`.

Two signs flipped by two different symmetries are independent labels. For a
single handle neither is an invariant: all four combinations are equivalent.
Each becomes meaningful only relative to a second handle. That agrees with
the repository's own reading of the Pin sign as a conserved **pair** label,
and it means the flux sign is a second, separate pair label. My round-14
sentence that the flux sign "should be tied to the Pin sector label" was a
guess, and it is retracted.

### Note added after the review of round 14

That review ([klein_handle_electric.md](klein_handle_electric.md), R3) showed
that the real and imaginary parts of the field are gauge dependent in the
gauged model. Checked against this round:

* The quarter rotation `Ψ → iΨ` **is** a constant gauge transformation in
  round 14's model (`5.6e-17`). So there the Pin sign of a single handle is
  a gauge choice, not only a convention. The conclusion above, that it is
  not an invariant of one handle and is independent of the flux sign, stands
  and is stronger than stated. Charge conjugation is a genuine symmetry, so
  the flux sign distinguishes two physically distinct mirror configurations.
* Whether a **relative** Pin sign between two handles is gauge invariant in
  a gauged model is not analysed. The repository's own setting differs: its
  `U(1)` is the frame rotation, with the Levi-Civita connection fixed by the
  geometry and no independent gauge field.
* Parts I and II are unaffected: Part I is about maps on `S³`, and the
  charges in Part II are Noether forms of global symmetries of an ungauged
  model.

Round 14's chain entry also changed: the orienting datum is trapped flux
**or a matter current**, not trapped flux alone. The table below is updated.

## What this does and does not show

Shown: how each recorded map acts on fibre-charged functions; that they are
attached to different spheres; that PR #129's condition on a bulk scalar is
the fibre-preserving cell; that the transport seam needs the antiperiodic
loop for the doublet; which charges survive in each cell; and that the Pin
sign and the flux sign are independent.

Not shown:

* **Which field BAM's matter is.** M4 and M6 give different cells for a bulk
  scalar and for a spin-frame field. The repository uses both: a scalar with
  PR #129's condition for its exchange kernel, and `J` for its spinor
  structure. This round shows those two usages sit in different cells. It
  does not say which one is matter, and `P` against `PT` remains the chosen
  datum `T`.
* **A physical pairing.** M3 is a representation identity. In the
  boundary-pairing protocol's vocabulary the spinor matching stays
  `NOT_DERIVED`: no action-derived boundary map produces `J`, and the
  protocol's warning against identifying two-dimensional spaces because
  their sizes agree is not lifted by this.
* **A field theory on the handle.** The lattice doublet keeps only the
  lowest fibre harmonic with scalar hopping. The identification of `W` with
  the spin holonomy is by matching `σ² = −I₄`; the neck geometry is not
  rebuilt, and no gauge field is coupled to the doublet.
* **Anything about probability.** The pairing is indefinite. No norm,
  detector law or event frequency follows.

## The chain after four rounds

| step | what it gives | what it assumes |
|---|---|---|
| 12 | two complementary quadrature sectors | a seam that reverses time |
| 13 | the sectors are `Re` and `Im` of one field; charge is their pairing | a charged field; a seam that also reverses the fibre |
| 14 | the sheet pairing is the electric flux through the seam | a conjugation-odd initial datum: trapped flux, or a matter current |
| 15 | the recorded transport supplies the fibre reversal for the spinor doublet | a spin-frame field, not a bulk scalar; PT; the `−1` spin holonomy |

Round 15 removes one free input and names two. The fibre reversal no longer
has to be assumed for a spin-frame field: it is how the recorded transport
acts. In exchange, the chain now visibly requires that matter be a
spin-frame field and that the time extension be PT, and it shows that the
scalar construction the repository actually uses for PR #129 does not
supply the pairing.

## Next

1. **Decide the field type from an action.** The cleanest test is whether a
   first-order (Dirac-type) field on the handle, with the recorded Pin
   structure, has the transport seam forced on it, where a scalar does not.
   That would replace "a spin-frame field" in the table by a derivation.
2. **The datum `T`.** Nothing in rounds 12–15 selects PT over P. The table
   in M6 says what each choice costs: PT loses total spin, P loses total
   charge. A constraint from the repository's pair-creation picture, where
   charge is conserved and spin is not fixed, would bear on it.
3. **Two handles.** The Pin sign and the flux sign are both pair labels.
   Whether the bordism rule of PR #279, which forces opposite Pin sectors,
   constrains the relative flux sign needs two handles on one network.
