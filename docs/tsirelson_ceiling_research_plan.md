# Probe I — what caps BAM's nonlocality at Tsirelson? (PR #236)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## 0. The question the arc had not asked

#206 states the program's sharpest ontological constraint in its own
words:

> A single classical field on 3-space with local dynamics and local
> readout is a local-hidden-variable model, and cannot violate CHSH.
> That is Bell's theorem itself.

So BAM spends its one legitimate nonlocal resource — the bridge — and
reaches CHSH = 2√2. But **nonlocality is not a dial that stops where
quantum mechanics stops.** Popescu–Rohrlich boxes are perfectly
no-signaling and reach the algebraic maximum CHSH = 4. Tsirelson's number
appears in eight documents of this repository, always as the value being
matched, never as something the geometry *explains*.

> **Can the bridge be driven past Tsirelson — and since it cannot, which
> structure forbids it: the geometry, or the imported formalism?**

The two answers are not equivalent. If the geometry forbids it, the
classical-origin conjecture has derived a quintessentially quantum bound
from classical ingredients. If the formalism forbids it, then #232's
blanket caveat *"quantization is imported, not derived"* has a precise —
and much smaller — residue, and that residue is exactly what the
conjecture still owes.

**It is the formalism. And the residue is one identity.**

## 1. The floor is bridge-carried

| bridge coupling `gb` | raw `\|c₊₋\|` | raw `\|c₋₊\|` |
|---:|---:|---:|
| 0.8 (committed) | 3.747e-01 | 3.747e-01 |
| 0.4 | 1.999e-01 | 1.999e-01 |
| 0.1 | 5.101e-02 | 5.101e-02 |
| **0.0 (cut)** | **4.140e-16** | **4.671e-16** |

and without the handle the model is a local-hidden-variable model: the 16
deterministic local strategies cap at CHSH = **2** exactly. This is
#206's result, re-derived as the baseline.

> **A reporting trap worth recording.** `extract_pair_state` normalizes,
> so the *cut* bridge still yields a unit-norm "pair state" — made
> entirely of numerical noise. The raw amplitudes above are what actually
> exhibits the cut.

## 2. The geometry's settings are exactly sufficient — and no more

BAM's local setting is a **fiber-frame rotation**, which on the `k = ±1`
qubit is `diag(e^{iθ}, e^{−iθ})` — a z-rotation. That is a **one-parameter
U(1)**, not SU(2). Conjugating the fixed transverse readout sweeps the
x–y plane and nothing else.

It is exactly enough, for a reason worth stating: **CHSH's optimum lies in
a plane.** The lattice-derived state is the singlet at fidelity
`1.000000000000`; the x–y block of its T-matrix has singular values
`(1, 1)`; the plane-restricted Horodecki value is therefore `2√2`
analytically — and is attained at every grid resolution tested:

| grid `n` | 25 | 49 | 97 | 193 |
|---|---|---|---|---|
| max CHSH | 2.828427125 | 2.828427125 | 2.828427125 | 2.828427125 |

Grid-independent, so a saturation and not a discretization artifact.
**Neither a sub-quantum deficit nor any excess** — precisely the
one-parameter family needed to reach the ceiling.

## 3. It never exceeds

Maximizing over every bridge the model has — all 8 gluing holonomies `s`
(the π-holonomy handle is only one of them), all preparations `(β, α)`:

| max CHSH | Tsirelson | excess |
|---:|---:|---:|
| 2.828427125 | 2.828427125 | **+8.9e-16** |

So BAM does not predict super-quantum correlations. On this axis it agrees
with experiment **for a reason**, and it is falsifiable in the other
direction.

## 4. But locality is not what caps it

The obvious guess is that the two mouths' operations act on commuting
algebras — and they do, overwhelmingly:

```
‖[P_A, P_B]‖ = 4.7e-113        mouth-mode overlap = 2.3e-113
```

because the committed bridge is a **link in `H` between distinct lattice
sites**, not an identification of degrees of freedom. (Worth stating
plainly: #206's prose — "one object", "two frames of one bulk defect" —
reads like an identification, and the implementation is not one.)

**Dropping commutativity does not help.** For *any* dichotomic Hermitian
observables on **one** Hilbert space, with no tensor split at all,

```
(B₀ + B₁)² + (B₀ − B₁)² = 4·I          [residual 5.3e-15, 800 draws]
```

and Cauchy–Schwarz then gives `CHSH ≤ 2√2` **without invoking locality
anywhere**. Measured, taking the optimal `A` for each `B`: the best
attainable value is **2.828427095** — Tsirelson again.

> So the ceiling is not supplied by the bridge, nor by the mouths being
> far apart, nor by anything else in the geometry. It is supplied by the
> **readout being a ±1-valued observable on a Hilbert space**.

## 5. The residue of "quantization is imported"

| ingredient of a Bell test | supplied by |
|---|---|
| the pair state (the singlet) | **geometry** — field + bulk gluing, fidelity 1.000000000000 |
| the nonlocality (violating Bell at all) | **geometry** — bridge-carried; 3.7e-01 with the handle, 4.1e-16 without |
| the measurement settings | **geometry** — a U(1), exactly sufficient |
| **the ceiling (why not 4?)** | **not the geometry** — `B² = I` alone |

Three of four ingredients are genuinely geometric. The fourth is not, and
the *entire* missing ingredient is one algebraic identity: `B² = I`. Not
Hilbert space in general, not the Born rule, not linearity — just the
readout being two-valued.

**That identity is what a classical-origin derivation still owes, and
naming it is the progress here. Deriving it is not attempted and is not
claimed.** The narrowing is real but it is a narrowing, not a closure. An
argument that a classical throat readout *must* be two-valued would finish
it; the repo's Pin⁻ / `T² = −1` machinery (#195/#196, #227–#231) is where
such an argument would have to come from, and none of it currently
addresses the readout algebra.

## 6. No-signaling is an equal-time statement

#206's T6 checks that Bob's reduced state is invariant under Alice-side
unitaries. That is a **partial-trace identity true of every bipartite
state** — it tests nothing about the geometry.

Tested on the lattice instead — Alice's setting applied as an operator
supported at mouth A, then evolution, then Bob's channel populations:

| `gb` | `dt = 0` | `0.5` | `2` | `8` |
|---:|---:|---:|---:|---:|
| 0.8 (committed) | 0 | 2.273e-02 | 1.019e-01 | **2.315e-01** |
| 0.2 | 0 | 5.262e-03 | 1.344e-02 | 3.600e-02 |
| 0.0 (cut) | 0 | 1.5e-15 | 1.2e-15 | 1.5e-15 |

At `dt = 0` the shift is exactly zero — the supports are disjoint. But the
committed handle is **always on**, so once any evolution separates the
setting from the readout, Alice's choice propagates to Bob, scaling with
the coupling and vanishing when the handle is cut.

This substantiates with numbers the caveat #206 states only in prose —
*"the physical bridge is non-traversable; the lattice handle is a modeling
stand-in"* — and converts it into a **requirement**: the bridge must be
dynamically inert at measurement time, which the committed model is not.
A constraint on the model, not a refutation of the program — but the kind
of thing that has to be discharged rather than asserted.

## 7. Tests

| # | test | outcome |
|---|---|---|
| T1 | goal | posed |
| T2 | the floor is bridge-carried, and 2 without it | 3.7e-01 → 4.1e-16; local max 2 |
| T3 | the U(1) of settings is exactly sufficient | 2√2 at every grid `n` |
| T4 | no bridge exceeds Tsirelson | excess +8.9e-16 |
| T5 | the cap is `B² = I`, not locality | 2.828427095 without commutativity |
| T6 | the residue of imported quantization | one identity |
| T7 | no-signaling is equal-time only | 2.3e-01 at `dt = 8` |
| T8 | assessment | 8/8 |

## 8. Open

  - **Whether a classical throat readout can be argued to be two-valued.**
    That would close the gap this probe locates. The Pin⁻ / `T² = −1`
    machinery is the natural source and does not currently address the
    readout algebra.
  - **The handle must be made dynamically inert at measurement time**, or
    the model signals (§6).
  - This is **one** 2-outcome, 2-setting scenario on one lattice. Nothing
    here bounds multipartite or higher-input scenarios, where the quantum
    set is not characterized by a single inequality.

## 9. Reproduce

```bash
python -m experiments.closure_ledger.tsirelson_ceiling_probe
```

Expected verdict:
`BAM_REACHES_TSIRELSON_EXACTLY_AND_NEVER_EXCEEDS_IT_BUT_THE_CEILING_IS_NOT_GEOMETRIC_IT_FOLLOWS_FROM_DICHOTOMIC_READOUT_ALONE_SO_THE_IMPORTED_QUANTIZATION_REDUCES_TO_B_SQUARED_EQUALS_I_AND_NO_SIGNALING_IS_EQUAL_TIME`, 8/8 PASS (~2 min).

## 10. Cross-references

  - `docs/configuration_space_emergence.md` — #206, whose lattice, state
    and no-go this builds on and whose T6 this replaces with a real test
  - `docs/tensor_product_emergence.md` — #232, whose "quantization is
    imported, not derived" this prices
  - `docs/measurement_sector.md` — #209, the Born-rule/pointer sector,
    untouched here
