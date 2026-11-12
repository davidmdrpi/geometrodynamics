# A responding electric field on the Klein handle: the sheet pairing is an electric flux, and two inputs can drive it

**The seven frozen predictions hold, and three of the conclusions I drew
from them were too strong.** Review of `562af71` found all three; they are
reproduced and corrected here (R1–R3). The corrected probe passes 45/45
required checks at both network sizes.

What stands:

1. Gauss's law gives the pairing between the two quadrature sectors a
   **gauge-invariant realisation**: its value on one sheet is the electric
   flux through the seam.
2. The complementary-quadrature structure survives electric response, in
   the ambient phase space and, as now tested, in the Gauss-reduced one.
3. Flux trapped around the handle loop is admissible in exactly the two seam
   types where round 13 found a charge.

What was wrong:

* I wrote that the pairing "becomes dynamically active only through trapped
  handle flux". **False.** Trapped flux is sufficient, not necessary.
  Admissible data with no charge, no electric field and no flux, but with a
  matter phase structure, generate sheet charge while the handle-loop flux
  stays zero throughout (R1).
* I described the tangent test as done "with Gauss's law imposed". **It was
  not.** Its perturbations left the constraint surface (R2).
* I presented the generated imaginary part of the field as a physical
  orientation. **It is gauge dependent** (R3).

The corrected statement of the mechanism: the sheet pairing is driven by
charge-conjugation-odd initial data, and there are two kinds — trapped
electric flux, and a matter current. Neither is derived here.

Also still true from the first write-up: I expected the trapped flux to be
screened and it is not reliably (S1), and charge separates across the seam
only when the seam reverses time (S2).

*Pre-registered in [klein_handle_electric_prereg.md](klein_handle_electric_prereg.md)
(`4ab8bcd`), committed before the module and unchanged. Reviewed at
`562af71`. Module `geometrodynamics/bulk/klein_handle_electric.py`; probe
`klein_handle_electric_probe.py`; tests `tests/test_klein_handle_electric.py`
(34). Archives: the run as first written in
`experiments/closure_ledger/runs/20261110_klein_handle_electric/`, the
corrected run in `runs/20261110_klein_handle_electric_review/`.*

```bash
python -m experiments.closure_ledger.klein_handle_electric_probe   # 45/45
python -m pytest -q tests/test_klein_handle_electric.py             # 34 tests
```

| verdict field | value |
|---|---|
| numerical_validity | `PASS` |
| support_dynamics | `EFFECTIVE_SUPPORT_ASSUMED` |
| electric_dynamics | `LATTICE_GAUGE_FIELD_SOLVED_FROM_THE_ACTION` |
| structure | `AMBIENT_LAGRANGIAN_HALF_SURVIVES_ELECTRIC_RESPONSE` (narrowed) |
| constrained_structure | `LAGRANGIAN_IN_THE_GAUSS_REDUCED_SPACE` (new) |
| flux_admissibility | `HANDLE_FLUX_ADMISSIBLE_IFF_CHARGE_UNTWISTED` |
| spontaneous_pairing | `NOT_GENERATED_FROM_REAL_DATA_WITHOUT_FLUX` (narrowed) |
| matter_phase_route | `MATTER_CURRENT_GENERATES_SHEET_PAIRING_WITHOUT_FLUX` (new) |
| flux_pairing | `TRAPPED_FLUX_DRIVES_GAUGE_INVARIANT_CURRENT_ODD_IN_FLUX` (narrowed) |
| gauge_dependence | `IMAGINARY_PART_IS_GAUGE_DEPENDENT` (new) |
| flux_dynamics | `FLUX_RESPONDS_AT_SECOND_ORDER` |
| pairing_carrier | `SHEET_PAIRING_IS_THE_SEAM_FLUX_BY_GAUSS` |
| seam_specificity | `FLUX_COUPLING_ALSO_ON_THE_ORDINARY_HANDLE` |
| orienting_datum | `CONJUGATION_ODD_INITIAL_DATA_FLUX_OR_MATTER_CURRENT_NOT_DERIVED` |
| positive_response, spinor_matching, event_law | `NOT_DERIVED` |
| bam_cell | `NOT_ESTABLISHED` |

## The model

The cover network of round 13: a ring with one handle, two sheets, two
cross-connected tubes; 32 and 52 sites. One Hamiltonian, lattice scalar
electrodynamics in temporal gauge, with the two waist coordinates:

    H = ½ Σ|Π_i|² + ½ Σ_{e: i→j} |Ψ_j − e^{i g a_e} Ψ_i|²
        + ½ Σ (m² + V_i(q)) |Ψ_i|² + ½ Σ E_e² + Σ_tubes (p_q²/2M + K_s q²/2).

`a_e` is the link potential and `E_e` its conjugate electric field. The
Gauss function `G = div E − gρ` is conserved, and physical states have
`G = 0`. Nothing is prescribed. Phase-space dimensions 200 and 320.

On a graph a gauge field has no radiation; its free content is one mode per
independent cycle. Electric flux around the handle loop is source-free flux
held by the topology — Wheeler's charge without charge — here a dynamical
variable with a conjugate holonomy.

Temporal gauge leaves the time-independent gauge transformations. The
gauge-invariant observables are the electric field, the link currents
`J_e = Im(conj(U_e Ψ_i) Ψ_j)`, the charge density, loop holonomies, and
quantities built from them. The real and imaginary parts of `Ψ` separately
are not gauge invariant (R3).

Under the deck map `a → s_f a` and `E → s_t s_f E` on the swapped edge. All
four seam types are symmetries of the one `H` (to `2.6e-16`).

## What was corrected after review

### R1 (high). Trapped flux is sufficient, not necessary

E3 showed that a **real** field with no gauge field stays exactly so, and I
generalised it into "response alone does not couple the sectors". E3 is
correct for its preparation and no further. Real data are precisely the data
with no link current.

The reviewer's preparation is admissible Klein data with `x` swap-even, `y`
swap-odd and `Π = a = E = 0`. Its charge density, electric field and trapped
flux are all zero. Its link current is not (`8.2e-04`, `4.6e-04`). Evolved:

| size | amplitude | max sheet charge | max handle flux | max field |
|---|---:|---:|---:|---:|
| primary | 1/64 | `2.6e-06` | `2.1e-22` | `1.6e-05` |
| primary | 1/8 | `1.7e-04` | `1.1e-20` | `1.0e-03` |
| replicate | 1/64 | `4.5e-06` | `6.4e-23` | `9.3e-06` |
| replicate | 1/8 | `2.9e-04` | `5.4e-21` | `5.9e-04` |

Sheet charge appears at second order in the field amplitude (slopes
`2.0000, 2.0001, 2.0002`), the histories are admissible and obey Gauss's
law (`≤ 1.9e-15`), and the handle-loop flux and holonomy stay at zero. The
flux stays zero for a reason: these data are fixed by plain time reversal as
well as by the deck map, which together force the electric field to be
swap-odd at all times (measured exactly), and the handle loop is swap-even.

Three further facts:

* **Charge conjugation reverses the sheet charge, exactly.** The sign of the
  matter current is an orienting datum, as the sign of trapped flux is.
* **It does not need the electric coupling.** With `g = 0` the same data
  give nearly the same sheet charge (`1.6900e-04` against `1.6895e-04`
  with the coupling on) and no field. The sheet charge is charge transport by the
  matter current; with `g ≠ 0` Gauss's law makes the seam flux follow it.
  Round 13's own archive already showed a varying sheet charge for complex
  data. I did not connect the two.
* **It is not a numerical effect.** Halving the step changes the maximum by
  `1.4e-03` and `5.4e-05`; an independent DOP853 integration agrees to
  `4.7e-04` and `4.7e-05`, with Gauss's law held to `6.0e-14`.

So there are two routes to a nonzero sheet pairing, both starting from
conjugation-odd initial data: trapped electric flux with a real field, or a
matter current with no flux. They are distinguishable: the first carries
handle-loop flux, the second does not.

### R2 (high). The tangent test left the Gauss surface

The frozen E1 test perturbed along every deck-fixed direction. At a generic
admissible state those directions violate linearized Gauss by `0.7071`. It
tested an ambient property and said so nowhere.

Done on the constraint surface:

| | primary | replicate |
|---|---:|---:|
| phase space | 200 | 320 |
| independent Gauss constraints | 32 | 52 |
| ambient deck-fixed directions | 100 | 160 |
| of which tangent to `G = 0` | **84** | **134** |
| gauge directions among those | 16 | 26 |
| reduced space `{G = 0}/gauge` | 136 | 216 |
| deck-fixed directions in it | **68** | **108** |

The constrained dimensions are the reviewer's. On the surface the symplectic
form is degenerate, and its kernel is exactly the gauge directions (nullity
32 and 52, equal to the rank of the Hamiltonian fields of `G`). The
constrained fixed space is isotropic (`1.3e-15`), and after removing its
gauge directions it is exactly half the reduced space: Lagrangian there.
Propagated along a responding history the constrained tangents stay on the
surface (`1.4e-12`) and isotropic (`4.7e-12`).

So the structure does survive in the reduced space, and that is now a
tested statement. The frozen label is narrowed to the ambient property and a
separate field records the reduced one.

Not covered: states fixed by the deck map only up to a gauge transformation.
The count above shows the strictly fixed states already fill a Lagrangian
subspace of the reduced space, so nothing is missing locally, but other
components exist. Round 15 meets one, the sector with holonomy `π` around
the doubled loop.

### R3 (medium). The `Re → Im` rotation is gauge dependent

A time-independent gauge transformation with swap-odd phases preserves the
strict Klein deck condition (`≤ 1.2e-16`). Applied to the real initial data of
E4 it changes the imaginary norm from `0` to `0.0326` (`0.0225` at the
replicate size) and changes the whole imaginary history, while the electric
field, link currents, charge density, sheet charge, loop holonomy and energy
are unchanged (`≤ 1.8e-15`).

So E4's statement, an imaginary part first order in the flux and exactly odd
in it, is a result **in temporal gauge with zero initial link potential**.
The gauge-invariant version of the same experiment:

| observable | slopes in the flux | odd in the flux |
|---|---|---|
| link current | `0.9999, 0.9998` | exactly |
| sheet charge | `1.0000, 1.0000` | exactly |
| loop holonomy | `1.0000, 1.0000` | exactly |

And the gauge-invariant version of E3: from real data with no flux, the
link currents, charge density and electric field all stay exactly zero.

What the trapped flux orients is therefore a current and a charge
separation, not a rotation of real into imaginary. The sign of electric flux
relative to an oriented cycle remains meaningful.

## Results that stand as first written

### The structure survives electric response, in the ambient space (E1)

Generic complex admissible data on the Gauss surface, with nonzero link
potentials and moving mouths, in all four cells at four amplitudes:

| quantity | worst value |
|---|---:|
| seam relation `z(τ) = D z(−τ)` (time-reversing cells) | `8.5e-15` |
| fixed set invariant under the flow (time-preserving cells) | `3.5e-16` |
| Gauss function | `2.2e-15` |
| total charge | `6.9e-15` |
| energy drift | `2.4e-06` |

The seam relation is a statement about states on the Gauss surface and is
unaffected by R2. The fixed space is half the ambient phase space and its
ambient tangent space is Lagrangian (`1.5e-11`); the reduced statement is R2.

### Trapped flux is admissible exactly where the charge exists (E2)

| time | fibre | admissible flux dimension | flux through the tubes |
|---|---|---:|---|
| kept | kept | 2 | allowed |
| kept | reversed | 1 | forbidden (`4.5e-16`) |
| reversed | kept | 1 | forbidden (`4.5e-16`) |
| reversed | reversed | 2 | allowed |

### The flux responds (E5)

The handle-loop flux changes at second order in field amplitude (slopes
`2.0000, 2.0001, 2.0005`) and first order in the flux.

### The sheet pairing is the flux through the seam (E6)

Summing Gauss's law over sheet 1, whose only boundary is the two seam edges,

    g Q_1 = E_{θ, last} − E_{θ', last},

to `1.2e-17`. Total charge is zero — a closed network is neutral — and
`Q_1(τ) = −Q_1(−τ)`. This is Gauss's law, true in any gauge theory, and
every quantity in it is gauge invariant. The content is the identification:
what round 13 called the pairing of the two quadrature sectors on one sheet
is an electric flux through the mouths. This is the result the review
singles out as surviving, and R1 strengthens it: both routes produce a sheet
charge, and in both it equals the seam flux.

### Flux coupling is not specific to the time-reversing seam (E7)

The frozen control gives the same first-order response on the ordinary
handle.

## Supplementary measurements, not frozen

### S1. The flux is not reliably screened

Over twelve real initial fields per case, at the time of largest change, the
flux decreased in 4 and 5 of 12 on the Klein handle and 3 and 7 of 12 on the
ordinary one: 19 of 48. My expectation of screening came from the
low-frequency response of a gapped field and does not apply to these
non-stationary fields.

### S2. Charge separates only across the time-reversing seam

With trapped flux and a real field, the sheet charge reaches `1.5e-05` and
`5.4e-06` on the Klein handle and is zero (`≤ 4.6e-19`) on the ordinary one,
where swap-even data stay swap-even. Not a frozen prediction; recorded as an
observation with a test.

## What this does and does not show

Shown, for this lattice model: the pairing between the two quadrature
sectors has a gauge-invariant electric realisation, the flux through the
seam; it is driven by conjugation-odd initial data of two kinds, trapped
flux or a matter current; and the Lagrangian structure survives in the
Gauss-reduced space.

Not shown:

* **Where the orienting datum comes from.** Trapped flux and matter current
  are both initial data. Nothing derives either, or its sign.
* **Which route BAM uses, if either.** Deciding that needs the matching to
  the geometry, which round 15 addresses and does not settle.
* **A complex structure or a norm.** The pairing is indefinite and a closed
  network forces its total to zero. Nothing here is a probability.
* **Gauge-equivalent admissible states beyond the strictly fixed ones.**
* **Maxwell theory.** A graph gauge field has no radiation. This is not the
  boundary-pairing protocol's radiative channel, and not its Stage 1.
* **GR support.** The oscillator is an effective comparator.

## The assumption budget so far

| step | what it gives | what it assumes |
|---|---|---|
| 12 | two complementary quadrature sectors | a seam that reverses time |
| 13 | the sectors are `Re` and `Im` of one field; charge is their pairing | a charged field; a seam that also reverses the fibre |
| 14 | the sheet pairing is the electric flux through the seam | a conjugation-odd initial datum: trapped flux, or a matter current |

None of the three is derived from general relativity here, and none yields a
positive norm, a detector law or an event frequency.

## Next

1. **Compare the two routes on equal terms.** For matched energy, which
   gives the larger and longer-lived seam flux, and do they interfere when
   both are present? The review suggests this, and it is the direct way to
   ask what actually supplies BAM's boundary pairing.
2. **The matter-current route on the ordinary handle.** S2 used trapped
   flux. Whether a matter current separates charge across a time-preserving
   seam is untested here.
3. **The relative Pin sign of two handles.** Round 15's quarter rotation is
   a constant gauge transformation in this model, so a single handle's Pin
   sign is a gauge choice there; that is now noted in round 15. Whether a
   relative sign between two handles is gauge invariant is not analysed.
