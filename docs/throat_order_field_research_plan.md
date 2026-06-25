# The throat-order field q(t,r,θ): the throat as a topological defect (PR #178)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## One field for three discrete facts

The antipodal-mechanics arc has established three discrete facts about the
throat, each in a different language:

- **#174** — the throat sector is an **odd-k winding ladder**, forced and
  rigid to the non-orientable 5D geometry (the discrete charge k).
- **#175** — a smooth continuous geometry reaches the discrete sector **only
  through an amplitude-zero node** forced at the antipodal focus (the gate).
- **#176/#177** — real (weak-field) GR self-gravity reproduces the focusing
  **disperse/collapse threshold** (`M_c ∝ 1/G`) — the nucleation of a throat.

This probe introduces a single field that makes all three the **defect data
of one object**: a complex Ginzburg–Landau order parameter

```
q(t, r, θ) = |q| e^{iφ},     V(q) = (λ/4)(|q|² − q₀²)² .
```

The **ordered** vacuum `|q| = q₀` is the orientable bulk; the **topological
defects** of `q` are the throats.

## What each test shows (measured)

### T2 — The two phases (the Landau potential)

The Mexican-hat `V(q) = (λ/4)(|q|² − q₀²)²` has two phases:

| point | V | V″ | phase |
|---|---:|---:|---|
| `q = 0` | 0.25 | **−1** | disordered, **unstable** maximum (symmetric) |
| `|q| = q₀` | 0.00 | **+2** | ordered, **stable** minimum (broken symmetry) |

The bulk sits in the ordered minimum (`|q| = q₀`); the degeneracy over the
phase φ is the U(1) that a defect winds.

### T3 — The throat is a vortex defect

The radial Ginzburg–Landau profile `f(r) = |q|(r)` solving

```
f″ + f′/r − k² f / r² = λ f (f² − q₀²),   f(0) = 0,  f(∞) = q₀
```

exists for each winding `k`:

| k | core `|q|` | bulk `|q|` | core size `r(|q|=q₀/2)` |
|---:|---:|---:|---:|
| 1 | 0.00 | 1.00 | 0.80 |
| 3 | 0.00 | 1.00 | 1.40 |
| 5 | 0.00 | 1.00 | 2.02 |

`|q| = 0` at the **core** (the field must vanish where the phase is
undefined), healing to the bulk vacuum `q₀` away from it — a **localized
defect**, the throat itself. The core widens with `k`: a higher-charge
throat costs a wider core.

### T4 — Winding = the discrete k (odd, by the #174 grading)

The topological charge `∮∇φ/2π` is the integer winding (measured `1, 3, 5`
for `k = 1, 3, 5`). It is the homotopy class in `π₁(S¹) = ℤ`; it cannot
change continuously while `|q| > 0`, so the discrete throat ladder **is** the
defect's winding number. The realized sector is **odd-k** — the #174
orientability grading (the non-orientable 5D geometry forces the odd rungs).

### T5 — The defect core IS the antipodal node (#175)

At the core `|q| = 0`. The phase is undefined where the amplitude vanishes,
so a defect's core is exactly the **forced amplitude-zero node** of #175.
The gate is topological: a loop enclosing no node has winding `0` (the
continuous sector), and acquiring the discrete charge requires the amplitude
to pass through zero — there is **no continuous path between winding sectors
that keeps `|q| > 0`** everywhere. The antipodal focus of #175 is where the
order field is driven to zero, opening the core.

### T6 — Nucleation IS the focusing threshold (#176/#177)

The disordered `q = 0` state is unstable (`V″ = −1 < 0`). Under the
Ginzburg–Landau gradient flow `∂_t|q| = −λ|q|(|q|² − q₀²)`, `q = 0` is a
repeller and `q₀` the attractor (a `10⁻³` seed → `q₀`, a finite seed →
`q₀`), so once a region is driven off the symmetric point an ordered domain
with a fixed-winding core **nucleates**. The **trigger** — what drives the
region off `q = 0` — is the self-gravitating focusing of #176/#177: the
`M_c ∝ 1/G` collapse concentrates the field and crosses the disorder→order
transition. Nucleating a throat-defect is the crossing of that threshold.

## Honest scope

This is the **effective** Ginzburg–Landau / Landau order-parameter level: q
is introduced as the coarse-grained order field whose defects are the
throats, and the arc's three discrete facts are shown to be its defect data.
It does **not** derive `V(q)` — the coefficients `λ, q₀` from the 5D bulk
action — nor does it dynamically couple `q` to the self-gravitating metric of
#176/#177 (here the focusing is the trigger, not yet a solved-together
field–metric system). Those two — the microscopic `V(q)` from the geometry,
and the q–metric coupling — are the follow-ups. This probe establishes the
**unifying field and its defect structure**, not its first-principles
derivation from the bulk.

## What this establishes

The throat stops being three separate facts and becomes **one object**: a
vortex defect of the throat-order field `q(t,r,θ)`. Its core is the antipodal
node (#175), its winding is the discrete odd-k charge (#174), and its
nucleation is the self-gravitating focusing threshold (#176/#177).

## Reproduce

```bash
python -m experiments.closure_ledger.throat_order_field_probe
# Verdict: THROAT_ORDER_FIELD_INTRODUCED_DEFECTS_ARE_THROATS_UNIFYING_THE_ARC
```
