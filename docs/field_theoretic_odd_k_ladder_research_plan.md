# The field-theoretic odd-k ladder on the Berger sphere (PR #193)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The follow-up #192 promised

PR #192 deformed the locked lepton Hamiltonian — an instanton-transition
**surrogate** — through a declared fiber/base ingredient map, and found a
finite window (λ_break = 0.986) with a metric-fine-tuned electron
near-zero. Its stated follow-up was the field-theoretic version: the
**actual deformed wave operator on S³_λ**, no surrogate and no ingredient
map. This probe is that version.

## The operator and the sector reduction

The scalar Laplacian on the unit Berger sphere is the genuine SU(2) result
(#165, imported and re-validated): `Δ(j,m) = 4j(j+1) + 4m²(λ⁻²−1)`. A
mode's Hopf-fiber winding is `k = 2m` (fiber phase `e^{ikθ}`, period 2π),
so the k-winding sector is `m = k/2`, `j ≥ k/2`, and

```
E(j, k; λ) = 4j(j+1) − k² + k²/λ² ,
sector ground (j = k/2):   E_k(λ) = 2k + (k/λ)² .
```

This is a **Kaluza–Klein split derived from the spectrum, not assumed**:

- `(k/λ)²` is the fiber winding momentum on a fiber of length ∝ λ —
  exactly the `(k·2π/L_throat)²` throat term of the unified mass operator
  (THESIS PR #83), now derived;
- `2k` is the base part: the k-sector reduces to a **monopole problem on
  the base S² with charge q = k/2** (the winding *is* the charge on the
  base — the #42–#44 Hopf⟷charge geometry), whose ground eigenvalue
  `l(l+1) − q²` at `l = q` gives `4q = 2k`. Half-integer q for odd k is
  the double-valued (Pin-twisted) monopole bundle — the odd-k sector is
  twisted already on the base.

The monopole reduction is verified by an **independent numerical field
solve** — the Wu–Yang charged Laplacian on S², cell-centered finite volume
with the natural no-flux pole boundaries built in (faces at 0 and π where
sin θ = 0), N = 6000:

| k | q | ground (numeric) | exact | first excited | exact |
|---:|---:|---:|---:|---:|---:|
| 1 | ½ | 0.49999999 | 0.5 | 3.500000 | 3.5 |
| 3 | 3/2 | 1.49999994 | 1.5 | 6.499999 | 6.5 |
| 5 | 5/2 | 2.49999988 | 2.5 | 9.499999 | 9.5 |

## Result 1: absolute structural protection

In closed form, for **every** λ ∈ (0, ∞):

```
E_k(λ) = 2k + (k/λ)²  ≥  2k  ≥ 2 ,
E_{k+2}(λ) − E_k(λ) = 4 + (4k+4)/λ²  >  4 .
```

The squash limit stiffens the fiber momenta (E → ∞); the stretch limit
floors at the base monopole part 2k; the sectors can never cross or
vanish. Numeric sweep λ ∈ [0.05, 20]: min E₁ = 2.0025, min gaps 4.02 and
4.04. The {1,3,5} ladder of the actual wave operator has **no λ_break at
all** — where the #192 surrogate's electron level crossed zero at a 1.4%
squash, the operator's k=1 sector is bounded below by 2 at *any* squash.

The deck grading is λ-independent for a structural reason: the antipodal
point lies **on the Hopf fiber** (the fiber translation θ → θ + π), so a
winding-k mode picks up `e^{ikπ} = (−1)^k` under the deck map, and the
deck map is an isometry of **every** Berger metric (the squash only
rescales the fiber direction). Odd k = the antipodally-odd (Pin-twisted)
sector of RP³; even k = honest RP³ functions. The #183 algebra, realized
spectrally at every λ.

## Result 2: the hierarchy is not kinematic

The same operator's mass ratios are pinned to O(1) at **every** λ
(conformal ω = √(E+1)):

- μ/e-analog `ω₃/ω₁` ∈ [1.528, 3.0] over all λ (limits √7/√3 stretch, 3
  squash); at the round point 2.0 — vs observed **206.8** (factor ≥ 69
  away at closest approach).
- τ/μ-analog `ω₅/ω₃` ∈ [1.254, 1.667] — vs observed **16.8**.

No Berger deformation of the bare operator produces the lepton hierarchy.
Combined with #192 (the hierarchy *is* reached by the instanton dynamics,
at the price of a metric-fine-tuned near-cancellation), the division of
labor is now **measured**:

> **Structure from kinematics/topology** (absolutely protected — three
> gapped generations on every Berger sphere); **hierarchy from dynamics**
> (the instanton near-cancellation, metric-fine-tuned, with no counterpart
> in the wave operator).

## Controls

- **Coupling:** minimal (ω = √E) instead of conformal — ordering holds at
  every sampled λ and the μ/e-analog stays bounded by 3.
- **Even-k:** the untwisted sectors obey the same closed form (k=0 the
  constant mode, E₂ = 4 + 4/λ²) — the grading separates towers, it does
  not gap the untwisted one.
- **The #192 side-by-side:** the surrogate's λ_break = 0.98598 is
  reproduced by import; the field-level E₁ infimum is 2.

## Honest scope

- Scalar wave operator on S³_λ; the sectoring by fiber winding is exact
  (m is a good quantum number at every λ), but the throat is represented
  by its winding sector, not by a solved throat geometry — the leptons
  here are sector ground states, not self-gravitating solitons.
- The Berger family is one deformation axis (the one that separates fiber
  from base); this is protection over that family, not over all metrics.
- Mass ratios via conformal (and minimal, as control) frequencies; code
  units, unit S³ radius.

## Reproduce

```bash
python -m experiments.closure_ledger.field_theoretic_odd_k_ladder_probe
# Verdict: FIELD_THEORETIC_ODD_K_LADDER_ABSOLUTELY_PROTECTED_ON_EVERY_BERGER_SPHERE
#          _STRUCTURE_KINEMATIC_HIERARCHY_DYNAMICAL
```
