# Pre-registration: the two routes to the sheet pairing, on equal terms

Publish this file before implementing or running anything. Baseline: branch
head `dc91d89` (round 14 corrected, round 15). Seed `2026111101`.

**Prior state.** The round-14 review established that two kinds of
charge-conjugation-odd initial data drive the sheet pairing on the Klein
handle: flux trapped around the handle loop, and a matter current. Numbers
from that review are known: a matter-current preparation at rest reaches a
sheet charge of order `1e-4` at field amplitude `1/8`, and trapped flux of
`4e-3` with a moving real field reaches `1e-6` to `1e-5`. Nothing else below
has been computed. Freeze and implementation are written in one session by
one author; the only ordering evidence is that this commit precedes the
code.

## Question

At matched total energy, which route gives the larger and the longer-lived
seam flux, and do they interfere when both are present?

## Lessons applied from the two reviews

* Only gauge-invariant observables are used: sheet charge `Q_1`, seam flux
  `g Q_1`, handle-loop flux `Φ_h`, link currents, loop holonomy, energy.
* Every preparation lies on the Gauss surface.
* A result for one preparation is not generalised: each route is sampled
  over 12 draws at each of two network sizes, and special subclasses are
  named as such.
* Controls isolate what the time reversal contributes.
* Quantities I cannot predict are listed as measurements, in advance.

## Model and preparations

Round 14's model unchanged (`bulk/klein_handle_electric.py`): lattice scalar
electrodynamics with two responding mouths, `g = 1` unless stated, doubly
reversed seam, implicit midpoint with step `0.02`. Mouths start at rest.

All preparations are deck-fixed and satisfy Gauss's law. With `ε` the total
energy and `|C|` the number of edges of the handle loop:

* **M, matter current.** Complex admissible field (`x` swap-even, `y`
  swap-odd, `p_x` swap-odd, `p_y` swap-even), total charge zero, electric
  field the loop-free solution of Gauss's law, no trapped flux. Scaled to
  `H = ε`.
* **M0, matter current at rest.** The review's subclass: `x` swap-even, `y`
  swap-odd, `Π = a = E = 0`.
* **F, trapped flux.** Real admissible field (`x` swap-even, `p_x`
  swap-odd) with energy `(1 − f) ε`, plus flux `Φ = sqrt(2 f ε / |C|)`
  around the handle loop. Default `f = 1/2`.
* **F0, trapped flux with the field at rest.** As F with `p_x = 0`.
* **B, both.** The field of M with trapped flux added.

Energies `ε ∈ {2.5e-4, 1e-3, 4e-3, 1.6e-2}`. Short window `[0, 4]`; long
window `[0, 200]`. Twelve draws per size from the frozen seed.

## Analytic predictions

### P1. A flux with the field at rest gives no sheet charge at all

F0 data are fixed by the deck map and also by `C ∘ T_0`, charge conjugation
composed with plain time reversal (`a → −a`, `E → E`, `y → −y`,
`p_x → −p_x`). The first makes `Q_1` odd in time, the second even. So

    Q_1 ≡ 0     for F0,

to `1e-12`, at every energy and draw. **The flux route needs field
momentum.** M0, fixed by `T_0` alone, is under no such constraint.

### P2. Scaling with energy on the short window

`Q_1` is bilinear in the field. For M it is zeroth order in everything else,
so `max |Q_1| ∝ ε`. For F it is first order in `Φ ∝ sqrt(ε)`, so
`max |Q_1| ∝ ε^{3/2}`. Predicted log-log slopes over the four energies, from
the median over draws: `1 ± 0.1` for M and M0, `1.5 ± 0.15` for F.

### P3. Matter current dominates at matched energy, by a factor that shrinks as `sqrt(ε)`

`R(ε) = max|Q_1|_F / max|Q_1|_M` has slope `0.5 ± 0.2`, and the median `R`
is below `1` at all four energies. The energy at which the fitted line
reaches `R = 1` is a measurement.

### P4. The best split of energy between field and flux is one third in the flux

At fixed `ε`, `max|Q_1|_F ∝ sqrt(f) (1 − f)`, maximal at `f = 1/3`.
Predicted ratios to the `f = 1/3` value at `ε = 1e-3`, each within `10 %`
(median over draws):

| `f` | `0.1` | `0.5` | `0.9` |
|---|---:|---:|---:|
| ratio | `0.739` | `0.919` | `0.246` |

### P5. The matter route is transport; the flux route is electric

At `ε = 1e-3`, over `g ∈ {0.5, 1, 2}` with the same field data: `max|Q_1|`
is independent of `g` for M0 (slope `0 ± 0.1`) and first order in `g` for F
(`1 ± 0.1`), where F's flux is held at the `g = 1` value. The seam flux is
`g Q_1` in both.

### P6. The two routes occupy different electric sectors

For M0 the electric field is swap-odd at all times, so the handle-loop flux,
the loop holonomy and the loop current `Σ_C J_e` are zero to `1e-12`. F
carries handle-loop flux by construction. The sheet charge of either is the
swap-odd seam field.

### P7. No first-order interference; the relative sign is physical

Under charge conjugation `Q_1` is odd, `y` is odd and `Φ` is odd. A term
proportional to `x·y·Φ` would be even, so it is absent. Therefore, with
`Q^M[x, y]` the matter-current history, and `Q^F[x, Φ]`, `Q^F[y, Φ]` the
flux histories of its two real components,

    Q^B[x, y, Φ] = Q^M[x, y] + Q^F[x, Φ] + Q^F[y, Φ] + O(Φ²).

Predicted: the additivity defect has slope `2 ± 0.3` in `Φ` over
`Φ ∈ {1e-3, 2e-3, 4e-3}`. Reversing `Φ` at fixed matter current changes the
history; the simultaneous reversal of both is a symmetry (`1e-12`). So the
**relative** sign of the two orienting data is a gauge-invariant,
conjugation-invariant label.

### P8. Early flux change is set by the loop bond order

For a real field, moving or not, at `t = 0`: `J_e = 0`, `J̇_e = 0` and

    Φ̈_h(0) = − g² Φ B_0 / |C|,       B_0 = Σ_{e ∈ C} x_tail x_head.

Predicted for F at `t = 0.1`: `ΔΦ_h = −½ g² Φ (B_0/|C|) t²` within `5 %`,
for every draw with `|B_0| > 0.02 ‖x‖²`, and with the predicted sign in all
of them. The flux is screened where neighbouring field values on the loop
are aligned and anti-screened where they alternate. This is offered as the
explanation of round 14's unexplained S1, where 19 of 48 random fields
screened: for a random field `B_0` has no preferred sign.

### P9. Which reversal each route needs

* The M0 preparation is also fixed by the deck map of the **time-preserving**
  fibre-reversing seam, and its history is the same history. So the matter
  route's sheet charge needs the fibre reversal, not the time reversal.
  Predicted: fixed to `1e-12`, fixed set invariant along the history to
  `1e-9`.
* On the ordinary handle, swap-even complex data with zero charge keep
  `Q_1 ≡ 0` (`1e-12`), with or without trapped flux.
* Trapped flux is inadmissible on either singly reversed seam (round 14).

So the flux route needs both reversals to exist and to separate charge; the
matter route needs only the fibre reversal.

### P10. At long times the flux route gains

On `[0, 200]` the loop holonomy winds, `∮a ≈ Φ |C| t`, so the flux route's
response is not confined to first order. Predicted, qualitatively: the ratio
of root-mean-square sheet charges, `R_long = rms_F / rms_M`, exceeds the
short-window `R` at each of the four energies (medians).

## Measurements, declared in advance and not predicted

* The values of `R` and `R_long`, and the energy where the fitted `R`
  reaches `1`.
* Long-window rms sheet charge for M, M0, F and B, and its one-sided mean
  divided by its rms (a persistent pairing has this near `1`, an oscillating
  one near `0`).
* Flux retention `Φ_h(t)/Φ_h(0)` on the long window, and whether the flux
  changes sign.
* How long the matter current's initial link-current pattern survives: the
  first time its normalised overlap with itself falls below `1/2`.
* Agreement of the long-window statistics under a halved step and under an
  independent DOP853 integration, for one draw.

## Verdict fields

    numerical_validity         PASS | UNRESOLVED_NUMERICALLY | INCOMPLETE
    flux_route_requirement     FLUX_ROUTE_NEEDS_FIELD_MOMENTUM
    energy_scaling             MATTER_LINEAR_FLUX_THREE_HALVES
    short_time_ranking         MATTER_CURRENT_DOMINATES_BY_INVERSE_SQRT_ENERGY
    flux_fraction              OPTIMUM_AT_ONE_THIRD_IN_THE_FLUX
    coupling_dependence        MATTER_ROUTE_IS_TRANSPORT_FLUX_ROUTE_IS_ELECTRIC
    electric_sectors           ROUTES_OCCUPY_DIFFERENT_ELECTRIC_SECTORS
    interference               NO_FIRST_ORDER_INTERFERENCE_RELATIVE_SIGN_PHYSICAL
    screening_rule             EARLY_FLUX_CHANGE_SET_BY_LOOP_BOND_ORDER
    reversal_requirements      MATTER_ROUTE_NEEDS_ONLY_FIBRE_REVERSAL
    long_time_ranking          FLUX_ROUTE_GAINS_AT_LONG_TIMES
    support_dynamics           EFFECTIVE_SUPPORT_ASSUMED
    orienting_datum            NOT_DERIVED
    positive_response          NOT_DERIVED
    event_law                  NOT_DERIVED

Status vocabulary as in rounds 12–15. A prediction contradicted by a
resolved measurement is reported as wrong in the first paragraph of the
write-up.

## What may not be concluded

* **Which route BAM uses.** Round 15 left the seam cell dependent on the
  field type and on `P` against `PT`. P9 bears on that and does not settle
  it.
* **A physical energy scale.** The crossover energy is in lattice units of a
  toy model.
* Matched total energy is one definition of equal terms. Matching the
  conjugation-odd content, or the action, could rank the routes
  differently; it is not tried.
* A finite lattice without dissipation has no decay. "Longer-lived" here
  means statistics on a fixed window, not a lifetime.
* No norm, detector law or probability follows.
