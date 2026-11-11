# How the attachment selects the temporal symmetry of exchange

**The frozen predictions hold, 33/33 required checks at both network sizes.**
The answer to the question has two halves.

* **There is a rule.** For any conserved quantity, the parity in time of its
  value on one sheet is fixed by two things: whether the seam reverses time,
  and whether the symmetry generating that quantity commutes with the
  seam's deck map. On the transport seam with time reversal this makes
  charge exchange **odd** in time and spin exchange **even**: opposite
  symmetries for the two quantities carried by one field.
* **The attachment alone does not select it.** The rule takes the field type
  and the time extension as inputs. For the two attachments the repository
  records, that gives four different answers. This round turns that
  dependence into an observable; it does not remove it.

The rule is a symmetry consequence, and the numerics confirm that derivation
and implementation agree. The reading for BAM rests on round 15's matching,
which has not been reviewed.

*Pre-registered in [klein_handle_parity_prereg.md](klein_handle_parity_prereg.md)
(`95d7a80`), committed before the module and unchanged. Module
`geometrodynamics/bulk/klein_handle_parity.py`; probe
`klein_handle_parity_probe.py`; tests `tests/test_klein_handle_parity.py` (6).
Archive `experiments/closure_ledger/runs/20261111_klein_handle_parity/`.*

```bash
python -m experiments.closure_ledger.klein_handle_parity_probe   # 33/33
python -m pytest -q tests/test_klein_handle_parity.py             # 6 tests
```

## The rule

Let `X` generate a symmetry of `H`, with sheet value `N_1` and conserved
total `N_tot`. Round 15 showed the Noether density picks up `s_t κ_X` under
the deck map, `κ_X = ±1` as `X` commutes or anticommutes with it. Then:

| seam | `X` and the deck map | sheet value |
|---|---|---|
| reverses time | anticommute | `N_1(τ) + N_1(−τ) = N_tot`: **odd** about `N_tot/2` |
| reverses time | commute | `N_tot = 0`, `N_1(τ) = N_1(−τ)`: **even** |
| keeps time | commute | `N_1 = N_tot/2`: **constant**, no exchange |
| keeps time | anticommute | `N_tot = 0`, `N_1` **unconstrained** |

No zero-charge assumption is needed. Round 16's relation
`Q_1(s_t t) = −s_t s_f Q_1(t)` is the case `κ = s_f` with `N_tot = 0`.

## The test

Round 15's ungauged spinor doublet with responding mouths, generic
deck-fixed data, 6 draws at each size, histories on `[−4, 4]`. Charge is
`X = i`; spin is `X = iσ_a`. The transport `J` anticommutes with `i` and
commutes with spin; a linear seam commutes with both. Parts are in units of
the squared field norm, primary size:

| seam | time | quantity | predicted | even part | odd part |
|---|---|---|---|---|---|
| transport | reversed | charge | odd | `≤ 1.1e-13` | `4.4e-03` to `1.4e-02` |
| transport | reversed | spin | even | `1.1e-02` to `4.5e-02` | `≤ 6.7e-16` |
| transport | kept | charge | unconstrained | `9.1e-03` to `7.2e-02` | `6.5e-03` to `3.8e-02` |
| transport | kept | spin | constant | variation `≤ 1.1e-15` | |
| linear | reversed | charge | even | `8.0e-03` to `5.2e-02` | `≤ 5.5e-16` |
| linear | reversed | spin | even | `6.3e-03` to `5.0e-02` | `≤ 1.5e-15` |
| linear | kept | charge | constant | variation `≤ 2.2e-16` | |
| linear | kept | spin | constant | variation `≤ 4.7e-16` | |

Every relation holds in every draw, for all three spin components, at both
sizes. Where exchange is allowed it occurs in every draw, and the excluded
parity is absent to rounding. The histories satisfy the seam relation to
`5.2e-15` and conserve each total to `1.3e-12`.

On the transport seam with time reversal the sheet charge is odd about half
the total charge and the sheet spin is even with zero total. Charge that
has crossed the seam by time `τ` is returned by time `−τ`; spin is exchanged
symmetrically about the mirror time.

## Reading for the recorded attachments

Assembled from rounds 15–17. No further computation.

| field | recorded map | time extension | charge exchange | spin exchange |
|---|---|---|---|---|
| bulk scalar | antipode (PR #129) | PT | even | — |
| bulk scalar | antipode | P | none | — |
| spin-frame doublet | transport (PR #279) | PT | odd | even |
| spin-frame doublet | transport | P | unconstrained | none |

So what selects the temporal symmetry is the pair (field type, time
extension). The geometry of the attachment supplies the maps; which map acts
on the field, and whether time is reversed with it, are the two inputs round
15 identified as chosen. The gain is that the choice now has an observable
signature: PT against P is the difference between a parity in time and none;
scalar against spin-frame is the difference between even and odd charge
exchange.

## What this does and does not show

Shown: a general parity rule for exchange across a seam, verified for charge
and spin of the doublet on four seam types.

Not shown:

* **A selection by geometry.** Recorded in the verdict as
  `BY_FIELD_TYPE_AND_TIME_EXTENSION_NOT_BY_GEOMETRY_ALONE`.
* **That either attachment acts as tabulated.** That is round 15's matching,
  unreviewed. Rounds 12, 14 and 16 each had conclusions that outran their
  tests; this table should be read with that in mind.
* **Anything for data that are not deck-compatible**, a gauged doublet, or
  higher fibre harmonics.
* **A norm, detector law or probability.**

## Next

1. Have round 15 reviewed. Two rounds now rest on it.
2. Find something in the repository's pair-creation picture that fixes the
   time extension: for instance, whether creation and annihilation of a
   mouth pair require charge exchange to be odd about the mirror time.
3. A gauged doublet, so the parity statements can be made for the electric
   flux through the seam, as round 14 did for the singlet.
