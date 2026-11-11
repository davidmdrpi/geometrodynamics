# Pre-registration: how the attachment selects the temporal symmetry of exchange

Publish this file before implementing or running anything. Baseline: branch
head `3104e37`. Seed `2026111102`.

**Prior state.** Round 16's review established, for the gauged singlet,
`Q_1(s_t t) = −s_t s_f Q_1(t)` on zero-total-charge histories. Round 15
established which conserved quantities survive on each seam for the spinor
doublet. Nothing below has been computed. Freeze and implementation are
written in one session by one author.

## Question

Given a seam, what temporal symmetry does the exchange of each conserved
quantity across it have? And what does that say about the two recorded BAM
attachments?

## Derivation

Let `X` generate a symmetry of `H`, with Noether density `n_X` and sheet
value `N_1`. Round 15 showed `N_X ∘ D = s_t κ_X N_X`, `κ_X = ±1` as `X`
commutes or anticommutes with the deck map `D`; site by site this is
`n_i(Dz) = s_t κ_X n_{σ i}(z)`. With `N_tot = N_1 + N_2` conserved:

* seam reverses time (`z(τ) = D z(−τ)`): `N_1(τ) = −κ_X (N_tot − N_1(−τ))`.
  - `κ_X = −1`: `N_1(τ) + N_1(−τ) = N_tot`. The sheet value is **odd** in
    time about `N_tot/2`.
  - `κ_X = +1`: `N_tot = 0` and `N_1(τ) = N_1(−τ)`: **even**.
* seam keeps time (`D z(τ) = z(τ)`): `N_1 = κ_X N_2`.
  - `κ_X = +1`: `N_1 = N_tot/2`, **constant**: no exchange.
  - `κ_X = −1`: `N_tot = 0`, `N_1` **unconstrained**.

No zero-charge assumption is needed. This contains round 16's relation as
the case `κ = s_f`.

## Predictions

Round 15's ungauged doublet model (`klein_handle_matching.Multiplet`) with
responding mouths, generic deck-fixed data, both network sizes, 6 draws,
histories on `[−4, 4]`. Charge is `X = i`; spin is `X = iσ_a`. The transport
`J` anticommutes with `i` and commutes with spin.

| seam | time | charge | spin |
|---|---|---|---|
| transport `J`, antiperiodic loop | reversed | odd | even |
| transport `J`, antiperiodic loop | kept | unconstrained | constant |
| linear, periodic loop | reversed | even | even |
| linear, periodic loop | kept | constant | constant |

**S1.** Each stated relation holds to `1e-9` of the field scale, in every
draw.

**S2.** The symmetry is not vacuous: where exchange is allowed (odd, even,
unconstrained) the sheet value varies by more than `1e-6` of the field
scale in every draw; and the excluded parity component is absent (`1e-9`).
For "unconstrained", both the even and the odd part in time are present
(`> 1e-7`).

**S3.** On the transport seam with time reversal, charge exchange and spin
exchange have **opposite** temporal symmetry.

## Reading for the recorded attachments

No computation; a table assembled from rounds 15–17, conditional on round
15's matching, which has not been reviewed.

| field | recorded map | time extension | charge exchange | spin exchange |
|---|---|---|---|---|
| bulk scalar | antipode (PR #129) | PT | even | — |
| bulk scalar | antipode | P | none | — |
| spin-frame doublet | transport (PR #279) | PT | odd | even |
| spin-frame doublet | transport | P | unconstrained | none |

## Verdict fields

    numerical_validity     PASS | UNRESOLVED_NUMERICALLY | INCOMPLETE
    parity_rule            EXCHANGE_PARITY_FIXED_BY_TIME_REVERSAL_AND_COMMUTATION
    non_vacuous            ALLOWED_EXCHANGE_OCCURS_EXCLUDED_PARITY_ABSENT
    charge_against_spin    OPPOSITE_PARITY_ON_THE_TRANSPORT_SEAM_WITH_TIME_REVERSAL
    selection              BY_FIELD_TYPE_AND_TIME_EXTENSION_NOT_BY_GEOMETRY_ALONE
    bam_cell               DEPENDS_ON_FIELD_TYPE_AND_TIME_EXTENSION
    positive_response      NOT_DERIVED
    event_law              NOT_DERIVED

## What may not be concluded

* The attachment does not select the parity by itself. The field type and
  the chosen datum `T` do. This round turns that dependence into an
  observable; it does not remove it.
* The doublet model is ungauged and keeps only the lowest fibre harmonic.
* The parities are properties of deck-compatible histories. Nothing is said
  about other data.
* The reading for BAM rests on round 15's matching, which is unreviewed.
* A parity in time is not a norm, a detector law or a probability.
