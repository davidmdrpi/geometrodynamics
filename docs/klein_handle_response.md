# A localized, responding Klein handle: the structure survives, and the pairing is the charge

**All eight frozen predictions hold, 32/32 required checks at both network
sizes, with nothing retuned.** A handle whose seam reverses time keeps the
complementary-quadrature structure of round 12 when it is localized and when
its mouths respond, exactly and at every amplitude. The mouths' mechanical
response supplies **no** pairing between the two sectors. A pairing exists
only for a charged field, and only when the seam reverses the fibre as well
as time: it is the conserved charge, and the variable that couples to it is
the electric potential, which in this calculation is prescribed, not
responding. It is an indefinite pairing, not a norm.

Most of these results are consequences of one reversing symmetry. The
numerics confirm that the derivation and the implementation agree; they are
not independent evidence. No frozen prediction failed, so there is nothing
to report as wrong, but one sentence of the freeze was imprecise and three
controls were added afterwards (C1–C4 below).

*Pre-registered in [klein_handle_response_prereg.md](klein_handle_response_prereg.md)
(`5d15963`), committed before the module and unchanged. Module
`geometrodynamics/bulk/klein_handle_response.py`; probe
`klein_handle_response_probe.py`; tests `tests/test_klein_handle_response.py`
(26). Archive `experiments/closure_ledger/runs/20261110_klein_handle_response/`.*

```bash
python -m experiments.closure_ledger.klein_handle_response_probe   # 32/32
python -m pytest -q tests/test_klein_handle_response.py             # 26 tests
```

| verdict field | value |
|---|---|
| numerical_validity | `PASS` |
| support_dynamics | `EFFECTIVE_SUPPORT_ASSUMED` |
| electric_response | `EXTERNAL_POTENTIAL_PRESCRIBED` |
| localized_structure | `LAGRANGIAN_HALF_SURVIVES_A_LOCALIZED_RESPONDING_HANDLE` |
| local_bracket | `SHEET_DATA_FREE_AT_MIRROR_WRONSKIAN_NOT_CONSERVED` |
| neutral_sector_pairing | `NOT_SUPPLIED_BY_MECHANICAL_RESPONSE` |
| charged_sector_pairing | `CHARGE_PAIRS_THE_SECTORS_ONLY_WITH_BOTH_REVERSALS` |
| positive_response, spinor_matching, event_law | `NOT_DERIVED` |
| bam_cell | `NOT_ESTABLISHED` |

## Scope

This answers the question the review of round 12 left:

> whether a localized, responding mouth pair preserves that structure — and
> whether its physical response supplies a pairing between the
> complementary sectors.

It is **not** Stage 1 of
[conjugate_boundary_pairing_protocol.md](conjugate_boundary_pairing_protocol.md).
It uses a lattice network, not the supported Maxwell bridge, and closes none
of that protocol's gates. It follows the protocol's rules where they apply:
the mouth coordinate is solved from the same Hamiltonian as the field, the
support is an explicit oscillator, and no bilinear form is called a norm or
a spinor pairing.

## The model

A ring of exterior sites with two ports joined by a tube. The tube's seam
may reverse time. The time-orientation cover is two copies of the exterior
with two tubes cross-connected; in cover time it is an ordinary static
network and the sheet swap `σ` is a spatial symmetry of it. Sizes: 32 and 52
cover sites.

    H = ½|p|² + ½ XᵀK(q,q')X + Σ_i A_i ρ_i + Σ_tubes (p_q²/2M + K_s q²/2),
    K = graph Laplacian + m² + V,   V_i = U_0 exp(−2 q h_i) on tube sites.

`V` is the lattice form of the protocol's `λN ℓ(ℓ+1)/r_q²`. Each tube has a
waist coordinate with `M q̈ = −K_s q + Σ h_i V_i |Ψ_i|²`: a waist-opening
load, solved, never prescribed. Integration is by implicit midpoint, which
is symplectic, symmetric and preserves quadratic invariants.

The seam is described by two signs, `s_t = −1` if it reverses time and
`s_f = −1` if it reverses the fibre of a charged field `Ψ = x + iy` (an
antilinear action), acting through the deck map `D`. `D` scales the
symplectic form by `s_t` and the charge density by `s_t s_f`.

## Results

Numbers are for the primary size unless a bound is given; bounds cover
both sizes. The replicate values are in the archive.

### 1. The localized handle keeps the Lagrangian half (H1)

For the neutral field with frozen mouths, the time-reversing deck map is an
antisymplectic involution and its fixed space is exactly half the phase
space (32 of 64; 52 of 104), with the symplectic form zero on it. A
time-preserving seam gives a symplectic subspace instead.

The two sectors have a concrete meaning on a handle. The cover operator's
`σ`-even spectrum is the spectrum of the **ordinary** handle, and its
`σ`-odd spectrum that of the handle with one **sign-twisted** tube edge
(agreement `≤ 5.4e-15` against independently built networks; the two
spectra differ by up to `0.256`). So an admissible field is a time-even
solution on the ordinary handle plus a time-odd solution on the twisted
one. Round 12's cosine and sine quadratures become these two networks.

### 2. Locally nothing is halved; the sheet Wronskian is not conserved (H2)

At the mirror time the Cauchy data on one sheet are free (rank 32 of 32; 52
of 52). An observer in the exterior has a full phase space. What the seam
changes is the bracket's conservation: for two admissible histories the
sheet-1 Wronskian satisfies `Ω_1(τ) = Ω_1(−τ)` (`≤ 1.6e-16`) and is **not
constant** — it varies by `0.0197` and `0.0221` for unit-norm pairs, while
the cover Wronskian stays zero. Symplectic flux crosses the seam.

This sharpens round 12. The halving is a statement about the cover. On the
handle itself, exterior evolution across an interaction with the mouths is
not a canonical transformation, though it is symmetric in time.

### 3. Response preserves the structure exactly (H3, H4)

The Hamiltonian is even in the field and the momenta, so `D` is an exact
reversing symmetry of the coupled nonlinear system: `D Φ D = Φ⁻¹`
(`1.2e-16`). Admissible histories satisfy `z(τ) = D z(−τ)`, to at worst
`6.1e-16` across both sizes:

| amplitude | seam defect | mouth excursion | energy drift |
|---:|---:|---:|---:|
| 1/64 | `2.5e-16` | `6.3e-06` | `2.6e-10` |
| 1/32 | `3.3e-16` | `2.5e-05` | `1.0e-09` |
| 1/16 | `3.1e-16` | `1.0e-04` | `4.1e-09` |
| 1/8 | `3.9e-16` | `4.0e-04` | `1.6e-08` |
| 1/2 | `2.5e-16` | `6.4e-03` | `2.6e-07` |

The tangent space to the admissible histories along a history with moving
mouths is Lagrangian at every sampled time (`1.1e-14`) and has half the
dimension (34 of 68; 54 of 108); the linearized scheme agrees with finite
differences to `1.8e-11`. Admissible nonlinear histories form a Lagrangian
submanifold.

The response is real and follows the protocol's order counting: mouth
excursion slopes `2.0000, 2.0000, 1.9999`, field-correction slopes
`3.0000, 2.9999, 2.9998`. The two mouths respond differently
(`max|q − q'|` up to `8.7e-05`) and mirror each other, `q'(τ) = q(−τ)`.

### 4. Four seam types for a charged field (H5, H6)

| time | fibre | fixed space | charge on it | `i` preserves it |
|---|---|---|---|---|
| kept | kept | symplectic | nondegenerate, `(32, 32)` | yes |
| kept | reversed | symplectic | **zero** | no |
| reversed | kept | Lagrangian | **zero** | yes |
| reversed | reversed | Lagrangian | nondegenerate, `(32, 32)` | no |

Either reversal alone twists the charge density and the total charge of
every admissible field vanishes identically. Both together untwist it. In
the doubly reversed cell the real part of `Ψ` lies in the untwisted sector
and the imaginary part in the twisted sector (`≤ 5.5e-16`): the two twist
sectors of round 12 are the real and imaginary parts of one field.

The seam admits a static potential only of the matching parity,
`A_{σ i} = s_t s_f A_i`. With it the seam relation holds to `1.7e-16`;
with the opposite parity it fails by `0.63` of the field scale.

### 5. Where the pairing is, and what does not supply it (H7, H8)

In the doubly reversed cell the charge `Q = ⟨x, p_y⟩ − ⟨y, p_x⟩` **is** the
round-12 Wronskian pairing between the two sectors (identity to `1.4e-17`).
It is nondegenerate, conserved with the potential on and the mouths moving
(drift `2.2e-14`), and it is carried through the seam: the sheet charge
varies by `0.035` while obeying `Q_1(τ) + Q_1(−τ) = Q` (`2.2e-14`).

* **The mechanical response supplies no coupling.** From admissible data
  with `Im Ψ = 0`, the imaginary part stays exactly `0.0` at every
  amplitude, with the mouths moving by up to `7.1e-03`. This is the `U(1)`
  invariance of the coupling `q|Ψ|²`.
* **The potential couples the sectors at first order**: slopes
  `1.0000, 1.0000` in the potential strength.

So the located variable is the electric potential across the handle. This
is the protocol's "which physical variable carries the antisymmetric
pairing", answered for this model: not the waist.

It is a pairing and nothing more. `Q` has signature `(32, 32)`, and
multiplication by `i` maps no admissible field to an admissible field.

## Corrections and controls added after the first run

### C1. An imprecise sentence in the freeze

H5 says a twist `η` in the antilinear cells "is a rotation of `Ψ` by `π`".
The correct statement, which is what the check tests, is that the `η = −1`
deck map is the `η = +1` deck map **conjugated by the quarter rotation**
`Ψ → iΨ` (error `0.0` in the antilinear cells, `2.0` in the linear ones,
where the twist remains a genuine label). The conclusion drawn in the
freeze is unchanged.

### C2. The frozen "stress" amplitude barely moved the mouths

At amplitude `1/2` the waist excursion is `6.4e-03`. That is not a stress
test. Supplementary runs, not frozen and gating nothing:

| amplitude | seam defect | mouth excursion | barrier factor | energy drift |
|---:|---:|---:|---:|---:|
| 2 | `3.3e-16` | `0.115` | `0.816` | `2.2e-06` |
| 4 | `4.8e-16` | `0.446` | `0.454` | `9.3e-06` |

With the barrier more than halved the seam relation is still exact.

### C3. The potential's coupling is not specific to the seam

On the ordinary handle (no reversal) the same potential generates an
imaginary part with the same first-order slopes, `1.0000, 1.0000`. Minimal
coupling rotates `Re` into `Im` everywhere. What is specific to the doubly
reversed seam is that `Re` and `Im` are its two Lagrangian sectors, that the
charge exists at all, and the parity condition on the potential.

### C4. Sheet data away from the mirror time

"Free at the mirror time" holds by construction. Measured at
`τ = 1, 2, 4` the sheet-1 data of admissible fields still have full rank.
The conditioning of that map depends on the units of `x` against `p` and is
not reported as a physical number.

## What this does and does not show

Shown, for this lattice model: the complementary-quadrature structure is a
property of a time-reversing seam that survives localization and nonlinear
mouth response; the waist response does not pair the sectors; a charged
field on a seam reversing both time and fibre carries both sectors, with the
charge as their pairing and the electric potential as the variable coupled
to it.

Not shown:

* **Which seam BAM has.** The repository records a time-reversing antipodal
  map (PR #129) and a fibre-reversing mouth transport (PR #279). If both
  carried over to a fibre-charged field on a handle, that is the doubly
  reversed cell. No matching map from the Hopf and spin-frame data to this
  field is constructed, and the protocol forbids identifying structures
  because their descriptions agree. `bam_cell` is `NOT_ESTABLISHED`.
* **A responding electric variable.** The potential is static and
  prescribed. No Wheeler-flux solution is built, so "the response supplies
  the pairing" is **not** established; only that the mechanical response
  does not and that the electric sector is where to look.
* **Other couplings.** The "no mechanical coupling" result is for the
  `U(1)`-invariant coupling `q|Ψ|²`. A coupling that breaks `U(1)` is not
  examined.
* **Probability.** An indefinite conserved pairing is not a norm, a
  detector law or a Hilbert-space inner product. Friedman and Higuchi's
  obstruction to global states stands.
* **GR support.** The oscillator is an effective comparator.
* **A continuum statement.** The symmetry results are exact at any size;
  nothing here is an existence or uniqueness claim.
* A charged field is additional field content. The protocol's Maxwell
  master channel is neutral.

## For the boundary-pairing protocol

Three things Stage 1 can use, as inputs and not as results of that protocol:

1. **A symmetry oracle.** Any discretization of a time-reversing seam whose
   action is time-reversal invariant must satisfy `D Φ D = Φ⁻¹` and
   `z(τ) = D z(−τ)` to solver precision, at every amplitude. A failure is a
   discretization or coupling error, not physics.
2. **The order counting is confirmed in a comparator.** Mouth response at
   `A²` and field correction at `A³`, as the protocol predicts.
3. **A located candidate for the pairing variable.** The waist does not
   carry the antisymmetric pairing. If any variable does, this calculation
   points at the static-flux sector the protocol sets aside as "separate",
   coupled to a fibre-charged field. Testing that needs a responding
   electric degree of freedom and the matching map; neither exists yet.
