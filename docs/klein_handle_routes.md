# The two routes on equal terms: the matter current wins, and the flux lasts

**Nine of ten frozen predictions hold as frozen. One is contradicted, and
review of `76c2067` found three conclusions I drew from the others too
strong** (R1–R3 below). At matched total energy, over 12 draws at each of
two network sizes:

* **The matter current gives the larger sheet charge**, by a median factor
  of 23 to 430 at short times and of about three over a long window. The trapped-flux route was larger in 1 of 96 long-window
  comparisons.
* **The trapped flux is the longer-lived input.** Over a window of 200 time
  units it never changes sign in any of 96 histories. Its median excursion
  is about 1 %; individual histories move by up to 6.3 %. The matter
  current's initial pattern loses half its self-overlap in half a time unit.
* **Their contributions add at leading order in the matter amplitude**, and
  the relative sign of the two orienting data is a physical label.

Corrected after review:

* I wrote that the matter route "needs only the fibre reversal". Fibre
  reversal is **sufficient, not necessary**. Matter data on the seam that
  reverses time and keeps the fibre also develop sheet charge. Only the
  ordinary handle forbids it (R1).
* I wrote "no first-order interference". The symmetry forbids a term in
  `x·y·Φ` and allows `x²y²Φ`; the odd part of the additivity defect is
  nonzero and quartic in the field (R2).
* I wrote that the flux "never moves by more than 1.5 %". That was a bound
  on **medians** (R3).

The structural finding survives and is sharper: the pairing and the
Lagrangian halving are distinct properties, and the seam selects the
**parity in time** of the sheet charge.

The contradicted prediction is P8. I froze a formula for the early change of
the trapped flux with a 5 % tolerance; one draw misses by 6.4 %. The sign
rule in it holds in all 22 resolved draws, and the formula is the correct
leading order. The tolerance was mis-set. I then told the user the cause was
a neglected rate term; that term vanishes identically by parity. Both errors
are recorded below. P8 was also offered as the explanation of round 14's
unexplained screening result, and it is **not** one.

Most of the confirmed predictions are symmetry or power-counting statements
about one Hamiltonian. The numerics confirm the derivations; the ranking and
the long-window statistics are the measured content.

*Pre-registered in [klein_handle_routes_prereg.md](klein_handle_routes_prereg.md)
(`fd78693`), committed before the module and unchanged. Module
`geometrodynamics/bulk/klein_handle_routes.py`; probe
`klein_handle_routes_probe.py`; tests `tests/test_klein_handle_routes.py`
(32). Archives: the run as first written in
`experiments/closure_ledger/runs/20261111_klein_handle_routes/`, the
corrected run in `runs/20261111_klein_handle_routes_review/`.*

```bash
python -m experiments.closure_ledger.klein_handle_routes_probe   # about 10 minutes
python -m pytest -q tests/test_klein_handle_routes.py             # 32 tests
```

| verdict field | value |
|---|---|
| numerical_validity | `PASS` |
| flux_route_requirement | `FLUX_ROUTE_NEEDS_FIELD_MOMENTUM` |
| energy_scaling | `MATTER_LINEAR_FLUX_THREE_HALVES` |
| short_time_ranking | `MATTER_CURRENT_DOMINATES_BY_INVERSE_SQRT_ENERGY` |
| flux_fraction | `OPTIMUM_AT_ONE_THIRD_IN_THE_FLUX` |
| coupling_dependence | `MATTER_ROUTE_IS_TRANSPORT_FLUX_ROUTE_IS_ELECTRIC` |
| electric_sectors | `ROUTES_OCCUPY_DIFFERENT_ELECTRIC_SECTORS` |
| interference | `ADDITIVE_AT_LEADING_ORDER_IN_MATTER_AMPLITUDE_RELATIVE_SIGN_PHYSICAL` (narrowed) |
| screening_rule | **`COUNTEREXAMPLE`** (frozen tolerance) |
| screening_leading_order | `BOND_ORDER_FIXES_THE_SIGN_AND_LEADING_ORDER_OF_EARLY_FLUX_CHANGE` (post-freeze) |
| reversal_requirements | `AT_REST_MATTER_ROUTE_DOES_NOT_NEED_TIME_REVERSAL` (narrowed) |
| seam_cells | `MATTER_PAIRING_IN_EVERY_CELL_EXCEPT_THE_ORDINARY_HANDLE` (new) |
| temporal_symmetry | `TIME_PARITY_OF_SHEET_CHARGE_IS_SET_BY_THE_SEAM` (new) |
| long_time_ranking | `FLUX_ROUTE_GAINS_AT_LONG_TIMES` |
| flux_retention | `SIGN_RETAINED_IN_EVERY_DRAW_EXCURSIONS_UP_TO_SEVEN_PERCENT` (new) |
| long_window_numerics | `STATISTICS_STABLE_UNDER_STEP_AND_INTEGRATOR` |
| orienting_datum, positive_response, event_law | `NOT_DERIVED` |

## Setup

Round 14's model unchanged: lattice scalar electrodynamics with two
responding mouths on the doubly reversed seam, `g = 1` unless stated. Only
gauge-invariant observables are used: sheet charge `Q_1`, seam flux `g Q_1`,
handle-loop flux `Φ_h`, link currents, loop holonomy, energy. Every
preparation is deck-fixed, has zero total charge, and lies on the Gauss
surface; each is scaled to the same total energy `ε`.

| | field | electric field | gauge-invariant content |
|---|---|---|---|
| **M**, matter current | complex, moving | loop-free Gauss solution | link currents; no trapped flux |
| **M0**, at rest | complex, `Π = 0` | zero | link currents only |
| **F**, trapped flux | real, moving | flux around the handle loop | flux; no link current |
| **F0** | real, at rest | flux | flux only |
| **B**, both | as M | Gauss solution plus flux | both |

F puts a fraction `f` of the energy in the flux, default `1/2`. Energies
`2.5e-4, 1e-3, 4e-3, 1.6e-2`. Short window `[0, 4]`, long window `[0, 200]`.

## What holds

### P1. A flux with the field at rest gives no sheet charge at all

F0 data are fixed by the deck map, which makes `Q_1` odd in time, and by
charge conjugation composed with plain time reversal, which makes it even.
Measured: `Q_1 ≤ 1.0e-17` at every energy and draw, while the flux itself
does evolve. With a moving real field the sheet charge is at least
`3.0e-09`. **The flux route needs field momentum.** Round 14's flux runs had
it without my noticing that it mattered.

### P2, P3. Scaling and ranking on the short window

| size | energy | matter current | at rest | flux route | ratio F/M: median (range) |
|---|---:|---:|---:|---:|---|
| primary | `2.5e-4` | `3.7e-06` | `5.7e-06` | `3.6e-08` | `0.0054` (`0.0023` to `0.033`) |
| primary | `1e-3` | `1.5e-05` | `2.3e-05` | `2.9e-07` | `0.0108` (`0.0045` to `0.066`) |
| primary | `4e-3` | `5.9e-05` | `9.1e-05` | `2.3e-06` | `0.0215` (`0.0091` to `0.131`) |
| primary | `1.6e-2` | `2.4e-04` | `3.6e-04` | `1.8e-05` | `0.0426` (`0.0181` to `0.261`) |
| replicate | `2.5e-4` | `4.4e-06` | `5.2e-06` | `7.8e-09` | `0.0023` (`0.0006` to `0.010`) |
| replicate | `1.6e-2` | `2.8e-04` | `3.4e-04` | `4.0e-06` | `0.0188` (`0.0048` to `0.079`) |

Log-log slopes of the medians: `1.0000` and `1.0001` for the matter
current, `1.4998` and `1.5000` for the flux route, `0.4972` and `0.4998`
for their ratio. Every per-draw slope, 24 for each route, lies within
`0.005` of the predicted value. No draw at any energy has the flux route ahead.

The sheet charge is bilinear in the field. The matter route needs nothing
else, so it is linear in energy. The flux route needs one power of the
flux, which costs `sqrt(ε)`.

The fitted ratio reaches `1` at energies of `9` and `46` in lattice units,
far outside the tested range and outside the small-amplitude regime the
scaling assumes. I draw no conclusion from that extrapolation.

### P4. One third of the energy in the flux is optimal

| `f` | `0.1` | `1/3` | `0.5` | `0.9` |
|---|---:|---:|---:|---:|
| predicted `sqrt(3f)(1−f)·3/2` | `0.7394` | `1` | `0.9186` | `0.2465` |
| measured, primary | `0.7395` | `1` | `0.9185` | `0.2464` |
| measured, replicate | `0.7395` | `1` | `0.9185` | `0.2464` |

`f = 1/3` is best in all 24 draws.

### P5. The matter route is transport; the flux route is electric

Over `g = 0.5, 1, 2` with the same data, the matter current's sheet charge
changes by parts in `10⁵` (slopes `1e-05`, `4e-05`); the flux route's is
first order in `g` (`0.9999`, `0.9997`). The seam flux is `g Q_1` in both
(`8.1e-17`). So the matter route moves charge across the seam whether or not
there is an electric field, and Gauss's law reports it; the flux route is an
electric effect.

### P6. Different electric sectors

For the matter current at rest the electric field is swap-odd at all times:
handle-loop flux, loop holonomy, loop current and the swap-even field are
all zero (`≤ 6.1e-20`), while the sheet charge is `5.7e-05`. The flux route
lives in the swap-even loop mode and produces a swap-odd seam field only as
a response.

This is a property of the at-rest subclass. A moving matter current does
develop a small loop flux (`6.0e-06`), so the two sectors are not strictly
separated in general.

### P7. Additive at leading order; the relative sign is physical

Under charge conjugation `Q_1`, the imaginary part and the flux are all odd,
so a term in `x·y·Φ` is forbidden. Measured: the defect of

    Q^B[x, y, Φ] = Q^M[x, y] + Q^F[x, Φ] + Q^F[y, Φ]

is second order in the flux (slopes `2.012, 2.006` and `1.993, 1.997`), at
`1.4e-10` against a matter-current signal of `4.7e-05`. Reversing the flux
alone changes the history at first order; reversing both orienting data is
an exact symmetry.

So when both are present their contributions add at leading order in the
matter amplitude, with a relative sign that is gauge invariant and
conjugation invariant. Neither sign means anything alone; their product
does. The additivity is not exact: see R2.

### P9. The matter route at rest does not need the time reversal

*The table first printed here said the matter route needs the fibre
reversal. That was not tested and is false; see R1.*

The matter current at rest is fixed by the deck map of the time-preserving
fibre-reversing seam as well, exactly, and its history stays in that fixed
set while separating charge (`5.7e-05`). On the ordinary handle swap-even
data keep `Q_1` at zero (`1.6e-17`) despite a local charge density of
`7.8e-04`, with or without trapped flux. Trapped flux is admissible only
when both reversals or neither are present (round 14), and separates charge
only with both.

This narrows what round 12's time reversal contributes to the pairing: the
fibre reversal alone suffices. It is not the only sufficient condition (R1).

### P10 and the long window

| size | energy | rms M | rms F | rms B | ratio, long | ratio, short | F ahead | flux retained (min) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| primary | `2.5e-4` | `7.2e-06` | `2.4e-06` | `3.9e-06` | `0.318` | `0.0054` | 0 of 12 | `0.9997` |
| primary | `1e-3` | `3.0e-05` | `1.0e-05` | `1.5e-05` | `0.282` | `0.0108` | 0 of 12 | `0.9974` |
| primary | `4e-3` | `1.2e-04` | `5.4e-05` | `6.0e-05` | `0.422` | `0.0215` | 0 of 12 | `0.9986` |
| primary | `1.6e-2` | `4.7e-04` | `2.3e-04` | `2.8e-04` | `0.485` | `0.0426` | 0 of 12 | `0.9901` |
| replicate | `2.5e-4` | `6.2e-06` | `1.6e-06` | `2.8e-06` | `0.320` | `0.0023` | 0 of 12 | `0.9998` |
| replicate | `1e-3` | `2.5e-05` | `8.0e-06` | `1.2e-05` | `0.341` | `0.0047` | 0 of 12 | `0.9996` |
| replicate | `4e-3` | `1.0e-04` | `3.5e-05` | `4.5e-05` | `0.347` | `0.0094` | 0 of 12 | `0.9992` |
| replicate | `1.6e-2` | `3.6e-04` | `1.4e-04` | `2.0e-04` | `0.352` | `0.0188` | 1 of 12 | `0.9928` |

**P10 holds**: the long-window ratio exceeds the short-window ratio at every
energy, by a factor of 11 to 136. The reason given in the freeze is the
right one in outline: the loop holonomy winds steadily, so the flux route's
response is not confined to first order in the flux.

The rest of this table is measurement, declared as such in advance:

* **The long-window ratio is about one third and nearly independent of
  energy**, `0.28` to `0.48`. The short-window square-root law does not
  persist. The flux route catches up to the same order as the matter route
  and stops there. The matter route still leads, in 95 of 96 comparisons.
* **The flux is long-lived.** It never changes sign in any draw. The
  "retained" column is the **median** over draws of each history's minimum;
  per-draw extrema are in R3. In advance I estimated a slow oscillation that
  could reverse it within the window at the highest energy. That estimate
  was wrong by a wide margin.
* **The matter current's pattern is short-lived.** The initial link-current
  pattern drops below half its self-overlap at `t = 0.46` (`0.48` at the
  replicate size), at every energy. The sheet charge it drives keeps
  oscillating long after.
* **Neither is a steady pairing.** The one-sided mean divided by the rms
  sits between `0.15` and `0.55` for every route: slow beats, neither a
  constant nor a zero-mean oscillation on this window.
* **Both together are not more than the matter current alone** at matched
  energy: B has part of its energy in the flux, and its rms lies between
  those of F and M.

Numerics: halving the step changes the long-window rms by at most `5.1e-04`;
an independent DOP853 integration agrees to `1.7e-04` in rms and `0.44 %`
pointwise over the full window, with Gauss's law held to `3.2e-14`.

## What was corrected after review

### R1 (high). Fibre reversal is sufficient, not necessary

P9 tested the at-rest matter current on two seams and I wrote a table
claiming the matter route needs the fibre reversal. The cell with time
reversed and the fibre kept was never tested. Admissible complex matter data
in every cell, with zero total charge, zero sheet charge at the mirror time
and no trapped flux, at energy `1e-3`:

| time | fibre | median peak `Q_1` (range), primary | replicate | parity in time |
|---|---|---|---:|---|
| kept | kept | `≤ 1.4e-15` | `≤ 2.4e-15` | zero |
| kept | reversed | `2.4e-05` (`5.5e-06` to `6.7e-05`) | `1.1e-05` | unconstrained |
| reversed | kept | `1.5e-05` (`3.4e-06` to `5.1e-05`) | `9.7e-06` | **even** |
| reversed | reversed | `1.8e-05` (`8.5e-06` to `3.4e-05`) | `1.2e-05` | **odd** |

Every draw in the three reversed cells develops sheet charge; an independent
DOP853 integration agrees to `7.9e-04`. Only the ordinary handle forbids it.
So the matter route needs **at least one** reversal, and the comparable
sizes mean no cell is singled out as minimal by strength.

The reviewer's symmetry relation organises this. For zero-total-charge,
deck-compatible histories,

    Q_1(s_t t) = − s_t s_f Q_1(t),

verified in all four cells to `7.3e-15`. With no reversal it forces
`Q_1 = 0`. With only the fibre reversed it is no constraint. With time
reversed it fixes the parity: **even in time if the fibre is kept, odd if it
is reversed.** The seam does not decide whether a pairing exists so much as
what temporal symmetry the charge exchange has.

### R2 (medium). Additivity is leading order in the matter amplitude

Charge conjugation forbids `x·y·Φ` and allows `x²y²Φ`, and the separately
evolved components do not see each other's back-reaction. The part of the
additivity defect that is odd in the flux,
`(D(Φ) − D(−Φ)) / 2Φ`:

| size | energy | endpoint, first draw | median peak | change on halving `Φ` |
|---|---:|---:|---:|---:|
| primary | `1e-3` | `−4.034e-09` | `2.2e-09` | `3.7e-04` |
| primary | `1.6e-2` | `−1.033e-06` | `5.5e-07` | `3.7e-04` |
| replicate | `1e-3` | `−1.444e-09` | `5.7e-10` | `7.2e-04` |
| replicate | `1.6e-2` | `−3.716e-07` | `1.5e-07` | `7.2e-04` |

It is nonzero in every draw, independent of the flux, and quartic in the
field (slope `2.000` in energy at both sizes). The primary values agree with
the reviewer's `−4.03776e−9` and `−1.03383e−6` to the accuracy of the step.
In the frozen range this odd part is under 1 % of the defect, which is why
P7's slopes came out near `2`; the prediction passed without the claim built
on it being exact.

### R3 (medium). Flux retention: medians are not bounds

The long-window table reports the median over draws. Per-draw extrema of
`Φ_h(t)/Φ_h(0)` on `[0, 200]`:

| size | energy | median min | worst min | worst max | sign changes |
|---|---:|---:|---:|---:|---:|
| primary | `2.5e-4` | `0.9997` | `0.9931` | `1.0057` | 0 |
| primary | `1e-3` | `0.9974` | `0.9917` | `1.0127` | 0 |
| primary | `4e-3` | `0.9986` | `0.9565` | `1.0248` | 0 |
| primary | `1.6e-2` | `0.9901` | **`0.9373`** | `1.0520` | 0 |
| replicate | `2.5e-4` | `0.9998` | `0.9963` | `1.0057` | 0 |
| replicate | `1e-3` | `0.9996` | `0.9905` | `1.0164` | 0 |
| replicate | `4e-3` | `0.9992` | `0.9796` | `1.0303` | 0 |
| replicate | `1.6e-2` | `0.9928` | `0.9788` | **`1.0433`** | 0 |

The worst history loses 6.27 % and the largest increase is 4.33 %, the
reviewer's figures. What is supported is that the sign is retained in all 96
histories and that excursions grow with energy.

## What does not hold: P8

### The prediction

For a real field at `t = 0` the link currents and their rate vanish, and

    Φ̈_h(0) = − g² Φ B_0 / |C|,        B_0 = Σ_{e ∈ C} x_tail x_head,

so I predicted `ΔΦ_h = −½ g² Φ (B_0/|C|) t²` within 5 % at `t = 0.1` for
every draw with a resolved bond order.

### The result

| | primary | replicate |
|---|---:|---:|
| resolved draws | 11 | 11 |
| sign as predicted | 11 | 11 |
| worst error at `t = 0.1` | `3.9 %` | **`6.4 %`** |

The frozen tolerance fails at one replicate draw. The verdict field is
`COUNTEREXAMPLE` for the claim as frozen. This is a statement about the
prediction, not the numerics, and `numerical_validity` stays `PASS`.

### What is and is not true

I first attributed the miss to the rate of the bond order, a term of
relative size `Ḃ_0 t / 3B_0` that I had left out for a moving field. That
was wrong. `Ḃ_0` vanishes identically, because the field is swap-even, its
momentum swap-odd and the handle loop swap-invariant (measured `≤ 5.9e-16`).
Including the term changes nothing.

The correction that matters is second order in `t`, with a coefficient that
grows when the bond order is small compared with the corresponding sum over
momenta. Tested directly, after the freeze:

| `t` | `0.04` | `0.1` | `0.2` |
|---|---:|---:|---:|
| median error, primary | `0.26 %` | `1.7 %` | `6.6 %` |
| median error, replicate | `0.24 %` | `2.1 %` | `8.5 %` |
| worst error, primary | `0.57 %` | `3.9 %` | `15.6 %` |
| worst error, replicate | `0.74 %` | `6.4 %` | `26.1 %` |

The error falls as `t²` (slopes `2.00`, `2.22`). So the frozen formula is
the leading order and the sign rule is exact at early times; the 5 %
tolerance at `t = 0.1` was an underestimate on my part, by a factor that
depends on the draw.

### It does not explain round 14's screening result

The freeze offered P8 as the explanation of S1, where 19 of 48 random fields
screened. Half of that survives: a random field has no preferred sign of
bond order (5 of 12 positive at each size), so there is no built-in
screening even at the earliest times. But S1 measured the sign at the time
of largest change over `[0, 4]`, and the early sign predicts that in 10 of
12 draws at the primary size and **5 of 12** at the replicate, which is
chance. The later sign is not set by the initial bond order. S1 remains
unexplained beyond "there is no reason to expect screening".

## What this does and does not show

Shown, for this lattice model at matched total energy: the matter current is
the stronger route at every energy and time scale tested; the trapped flux
is the more durable input, retaining its sign in every history; a matter
current separates charge on every seam except the ordinary handle, with a
parity in time fixed by the seam; and the two routes add at leading order in
the matter amplitude, with a physical relative sign.

Not shown:

* **Which route BAM uses, or which seam.** Round 15 left the seam cell
  dependent on the field type and on P against PT. R1 shows the matter
  route exists in three of the four cells, so its existence does not
  discriminate between them; the parity in time does.
* **A minimal setting.** I called the fibre-reversing, time-preserving seam
  the minimal setting for the stronger route. With three cells admitting it
  at comparable strength, that is withdrawn.
* **A different definition of equal terms.** Matched total energy is one
  choice. Matching the conjugation-odd content, or the action, is not tried
  and could rank the routes differently.
* **A lifetime.** A finite lattice without dissipation has no decay. The
  durability statements are statistics on `[0, 200]`.
* **A physical energy scale.** Everything is in lattice units of a toy
  model.
* **A steady pairing.** No route produces a sheet charge with a stable sign
  over the long window.
* **Anything about probability.** No norm, detector law or event frequency
  follows.

## Consequence for the chain

| step | what it gives | what it assumes |
|---|---|---|
| 12 | two complementary quadrature sectors | a seam that reverses time |
| 13 | the sectors are `Re` and `Im` of one field; charge is their pairing | a charged field; a seam that reverses the fibre |
| 14 | the sheet pairing is the electric flux through the seam | a conjugation-odd initial datum |
| 15 | the recorded transport supplies the fibre reversal for the spinor doublet | a spin-frame field |
| 16 | the stronger datum is a matter current; any reversal admits it, and the seam fixes its parity in time | a matter current, not derived |

Round 16 changes the weight of the earlier steps. The pairing and the
halving, which rounds 12 and 13 presented as one structure, are distinct
properties. The halving needs the time reversal. A matter-driven pairing
needs at least one reversal of either kind. Where both properties are
present the seam fixes how they combine: time reversal alone gives a sheet
charge even in time, both reversals give one odd in time.

## Next

1. **How the attachment selects the temporal symmetry of charge exchange.**
   This is the reviewer's suggestion and the sharper question. An even
   sheet charge and an odd one are physically different exchanges: one is
   symmetric about the mirror time, the other reverses through it. Round 15
   tied the cell to the field type and to P against PT; this ties the cell
   to an observable parity.
2. **A steady pairing.** No preparation tried gives a sheet charge of stable
   sign. Whether one exists, for instance with a degenerate pair of even and
   odd modes, is untested.
3. **Another definition of equal terms.**
