# A responding electric field on the Klein handle: response does not make the pairing, trapped flux does

**All seven frozen predictions hold, 27/27 required checks at both network
sizes, with nothing retuned.** With the electric field solved from the same
Hamiltonian as the charged field and the mouths:

1. the complementary-quadrature structure survives;
2. electric response to the field, by itself, couples the two sectors
   **not at all** — exactly zero, by charge-conjugation symmetry;
3. what couples them is electric flux **trapped around the handle loop**,
   at first order, and the sign of that flux sets the sense of the
   `Re → Im` rotation. Such flux is admissible in exactly the two seam
   types where round 13 found a charge. By Gauss's law the pairing's value
   on one sheet is the electric flux through the seam.

So the answer to "does a responding electric degree of freedom supply the
pairing" is **no**: response carries and transports it, but does not create
it. It enters as an initial datum, the trapped flux, whose sign breaks
charge conjugation. Nothing here derives that datum.

Two things did not go as I expected, and one was not predicted:

* I expected the trapped flux to be **screened**. The sign was deliberately
  left unfrozen, and it is not robust: over 48 initial fields the flux
  decreased in 19 and increased in 29 (S1 below).
* **Charge separates across the seam only when the seam reverses time.** On
  the ordinary handle the sheet charge is exactly zero. This is the one
  effect in this round that is specific to the time-reversing seam (S2).

Most results are consequences of symmetries of one Hamiltonian. The numerics
confirm that derivation and implementation agree; they are not independent
evidence.

*Pre-registered in [klein_handle_electric_prereg.md](klein_handle_electric_prereg.md)
(`4ab8bcd`), committed before the module and unchanged. Module
`geometrodynamics/bulk/klein_handle_electric.py`; probe
`klein_handle_electric_probe.py`; tests `tests/test_klein_handle_electric.py`
(24). Archive `experiments/closure_ledger/runs/20261110_klein_handle_electric/`.*

```bash
python -m experiments.closure_ledger.klein_handle_electric_probe   # 27/27
python -m pytest -q tests/test_klein_handle_electric.py             # 24 tests
```

| verdict field | value |
|---|---|
| numerical_validity | `PASS` |
| support_dynamics | `EFFECTIVE_SUPPORT_ASSUMED` |
| electric_dynamics | `LATTICE_GAUGE_FIELD_SOLVED_FROM_THE_ACTION` |
| structure | `LAGRANGIAN_HALF_SURVIVES_ELECTRIC_RESPONSE` |
| flux_admissibility | `HANDLE_FLUX_ADMISSIBLE_IFF_CHARGE_UNTWISTED` |
| spontaneous_pairing | `NOT_GENERATED_BY_RESPONSE_WITHOUT_FLUX` |
| flux_pairing | `TRAPPED_FLUX_COUPLES_SECTORS_AND_SETS_ORIENTATION` |
| flux_dynamics | `FLUX_RESPONDS_AT_SECOND_ORDER` |
| pairing_carrier | `SHEET_PAIRING_IS_THE_SEAM_FLUX_BY_GAUSS` |
| seam_specificity | `FLUX_COUPLING_ALSO_ON_THE_ORDINARY_HANDLE` |
| flux_origin | `INITIAL_DATA_NOT_DERIVED` |
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
held by the topology — Wheeler's charge without charge — and here it is a
dynamical variable with a conjugate holonomy.

Under the deck map `a → s_f a` and `E → s_t s_f E` on the swapped edge.
With the gauge field dynamical there is no parity condition left to impose:
all four seam types are symmetries of the one `H` (to `2.6e-16`).

Implementation checks, not predictions: the vector field is Hamiltonian
against a finite-difference gradient, gauge transformations commute with the
flow, the Gauss function is conserved off the constraint surface, and at
`g = 0` the model reduces to round 13's.

## Results

Bounds cover both sizes. Tables are the primary size unless stated.

### 1. The structure survives electric response (E1)

Generic complex admissible data on the Gauss surface, with nonzero link
potentials and moving mouths, in all four cells at four amplitudes:

| quantity | worst value |
|---|---:|
| seam relation `z(τ) = D z(−τ)` (time-reversing cells) | `8.5e-15` |
| fixed set invariant under the flow (time-preserving cells) | `3.5e-16` |
| Gauss function | `2.2e-15` |
| total charge | `6.9e-15` |
| energy drift | `2.4e-06` |

The fixed space is half the phase space (100 of 200; 160 of 320) and the
tangent space along a responding history is Lagrangian (`1.5e-11` from
finite-difference tangents). The histories were not trivial: local charge
density up to `8.4e-04`.

### 2. Trapped flux is admissible exactly where the charge exists (E2)

| time | fibre | admissible flux dimension | flux through the tubes |
|---|---|---:|---|
| kept | kept | 2 | allowed |
| kept | reversed | 1 | forbidden (`4.5e-16`) |
| reversed | kept | 1 | forbidden (`4.5e-16`) |
| reversed | reversed | 2 | allowed |

The cover has three independent cycles. Flux around the handle loop is
swap-even and survives only when `s_t s_f = +1`, the same two cells in which
round 13 found a nonzero charge. In the mixed cells the only admissible flux
is the two rings opposed, with nothing through the handle.

### 3. Response alone does not couple the sectors (E3)

From real admissible data with no gauge field, at every amplitude including
`1/2`, with the gauge field free to respond and the mouths moving by up to
`7.2e-03`: the imaginary part, the link potentials and the electric field
all stay exactly `0.0`, on the Klein handle and on the ordinary one.

Charge conjugation is a symmetry of `H`, and a real field with no gauge
field is its fixed set. No amount of response leaves it.

### 4. Trapped flux couples the sectors and sets the orientation (E4, E7)

| size | seam | `Im` slopes in flux | oddness error |
|---|---|---|---:|
| primary | Klein | `1.0000, 0.9999` | `0.0` |
| primary | ordinary | `1.0000, 0.9999` | `0.0` |
| replicate | Klein | `1.0000, 0.9999` | `0.0` |
| replicate | ordinary | `1.0000, 0.9999` | `0.0` |

With flux `Φ` around the handle loop, `ȧ = E` winds the holonomy and the
hopping phases rotate the real part into the imaginary part, at first order
in `Φ`. Reversing `Φ` reverses the imaginary part and leaves the real part
unchanged, exactly. The sign of the trapped flux is the sense of rotation.

The frozen control E7 shows the same on the ordinary handle. The coupling of
`Re` to `Im` by flux is **not** specific to the time-reversing seam. What is
specific is which flux is admissible (result 2) and the charge separation
below.

### 5. The flux responds (E5)

The handle-loop flux is not conserved once charged matter is present:
its change is second order in field amplitude (slopes `2.0000, 2.0001,
2.0005`) and first order in the flux (`1.0000, 1.0000`). The sign is the
subject of S1.

### 6. The sheet pairing is the flux through the seam (E6)

Summing Gauss's law over sheet 1, whose only boundary is the two seam edges,

    g Q_1 = E_{θ, last} − E_{θ', last},

holds along the flux histories to `1.2e-17`. Total charge is zero to
`2.9e-18` — a closed network is neutral — and `Q_1(τ) = −Q_1(−τ)` to
`2.9e-18`, with `Q_1` reaching `6.7e-07`. Round 13's sheet pairing
`Ω_1(Re, Im)` is therefore the difference of the electric fields on the two
seam edges, divided by `g`.

This is Gauss's law, true in any gauge theory. The content is the
identification: the quantity round 13 called the pairing of the two
quadrature sectors on one sheet is an electric flux through the mouths.

## Supplementary measurements, not frozen

### S1. The flux is not reliably screened

The freeze expected screening, `ΔΦ·Φ < 0`, and left the sign unfrozen. The
frozen runs disagreed with each other: on the Klein handle the flux
increased at the primary size and decreased at the replicate. A survey over
twelve real initial fields per case, at the time of largest change:

| size | seam | flux decreased |
|---|---|---:|
| primary | Klein | 4 of 12 |
| primary | ordinary | 3 of 12 |
| replicate | Klein | 5 of 12 |
| replicate | ordinary | 7 of 12 |

On `[0, 4]` the sign depends on the field's state. My expectation came from
the low-frequency response of a gapped field; these are non-stationary
fields at finite time, where it does not apply. No statement about screening
is made.

### S2. Charge separates only across the time-reversing seam

Over the same fields, the largest sheet charge:

| size | seam | largest | smallest per-field maximum |
|---|---|---:|---:|
| primary | Klein | `1.5e-05` | `1.4e-06` |
| primary | ordinary | `4.6e-19` | `2.5e-20` |
| replicate | Klein | `5.4e-06` | `5.2e-07` |
| replicate | ordinary | `3.2e-19` | `2.8e-20` |

On the ordinary handle the sheet swap is a symmetry and swap-even data stay
swap-even, so the two seam edges carry equal field and `Q_1 ≡ 0`. On the
time-reversing seam, admissible data are swap-even only at the mirror time;
`Q_1` grows, odd in time, in every field tried. A sheet of the Klein handle
acquires a net charge that the same sheet carries with opposite sign at the
mirrored time, with the total always zero.

This was not a frozen prediction. It follows from the symmetries once seen,
and it is recorded as an observation with a test, not as a confirmed
prediction.

## What this does and does not show

Shown, for this lattice model: the pairing between the two quadrature
sectors is an electric quantity; response transports it and does not create
it; it becomes dynamically active only through trapped handle flux, which is
admissible only where the charge exists and whose sign fixes an orientation.

Not shown:

* **Where the flux comes from.** It is an initial datum. E3 shows response
  cannot produce it from a real field. The repository posits one Hopf
  charge per throat; whether that supplies this datum, and its sign, is not
  addressed.
* **Which seam BAM has.** As in round 13, `bam_cell` is `NOT_ESTABLISHED`.
  No matching map from the Hopf and spin-frame data to a fibre-charged field
  on a handle is constructed.
* **A complex structure or a norm.** An orientation for the `Re → Im`
  rotation is not a complex structure on the admissible set, the pairing is
  indefinite, and a closed network forces its total to zero. Nothing here is
  a probability.
* **The reduced phase space.** Statements are on the unreduced space with
  Gauss's law imposed. States fixed by the deck map only up to a gauge
  transformation are not classified.
* **Maxwell theory.** A graph gauge field has no radiation. This is not the
  protocol's radiative master channel, and not Stage 1.
* **GR support.** The oscillator is an effective comparator.

## The assumption budget so far

Each step of rounds 12–14 added an input. In order:

| step | what it gives | what it assumes |
|---|---|---|
| 12 | two complementary quadrature sectors | a seam that reverses time |
| 13 | the sectors are `Re` and `Im` of one field; charge is their pairing | a charged field; a seam that also reverses the fibre |
| 14 | the pairing is an electric flux; its coupling has an orientation | trapped flux around the handle, with a sign |

None of the three is derived from general relativity here, and none yields a
positive norm, a detector law or an event frequency. What has been located
is narrower and concrete: **which classical variable carries the
antisymmetric pairing, and what initial datum orients it.**

## Next

1. **The matching map.** Decide whether the repository's time-reversing
   antipodal map (PR #129) and fibre-reversing mouth transport (PR #279)
   act on a fibre-charged field as the doubly reversed seam does. That is
   the one input that would connect this chain to BAM's geometry.
2. **The origin and sign of the trapped flux.** If one Hopf charge per
   throat is the flux datum, its sign should be tied to the Pin sector label
   of PR #279, where opposite sectors are produced in pairs. That is a
   checkable correspondence.
3. **The reduced phase space**, where gauge-equivalent admissible states are
   identified.
