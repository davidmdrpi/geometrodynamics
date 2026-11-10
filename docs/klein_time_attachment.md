# The Klein time attachment: where a classical field starts to look quantum

**All four frozen predictions are confirmed, 20/20 required checks, with no
parameter changed after the freeze.** A classical real field on a
time-nonorientable quotient of the Einstein static universe has half a phase
space, a deck map that acts as complex conjugation, a sharp 0-or-π phase
closure between two caps, and — in a generic closed history — no free
radiation at all. Three of these are things the repository's ledger currently
records as chosen. None of them is a probability rule.

*Pre-registered in [klein_time_attachment_prereg.md](klein_time_attachment_prereg.md)
(`bb6ff2c`), committed before the module. Module
`geometrodynamics/bulk/klein_time.py`; probe `klein_time_probe.py`; tests
`tests/test_klein_time.py` (52). Archive
`experiments/closure_ledger/runs/20261110_klein_time/`.*

```bash
python -m experiments.closure_ledger.klein_time_probe   # 20/20
python -m pytest -q tests/test_klein_time.py             # 52 tests
```

## The construction

`docs/finite_mouth_bundle_transport.md` §4 records that PR #129's antipodal
map is PT-type — it reverses time orientation — while the alternative P-type
map does not, and adds: *"Which is physical is a further choice, recorded as
the datum `T`. For the static problem the two agree."* This round does the
dynamic problem, where they do not agree.

Background: `R × S³`, radius 1, conformally coupled real scalar, zonal
sector `φ = u(t,χ)/sin χ`. Then `u_tt = u_χχ` with `u = 0` at both poles,
exactly; modes `sin(kχ)`, `ω_k = k`; antipode `χ ↦ π − χ`. On the lattice
`dt = dχ = π/24` the leapfrog rule is exact propagation with exactly these
integer frequencies (period error `0.0`).

    P-type (control):   (P φ)(t,χ)   = φ(t, π−χ)
    PT-type cap at c:   (T_c φ)(t,χ) = η φ(2c − t, π−χ),   η = ±1

The PT quotient contains a one-sided spacelike `RP³` at `t = c`: a cross-cap
in time. It is locally the same Einstein solution as the cover, since a
fluid at rest has a stress tensor even under `u^a ↦ −u^a`.

## Results

### 1. Half a phase space

`T` reverses the conserved Wronskian, so its fixed fields satisfy
`Ω(a,b) = −Ω(a,b) = 0`.

| quotient | fixed dimension (of 46) | Wronskian on the fixed space |
|---|---:|---|
| P-type | 24 | rank 24: a symplectic subspace |
| PT-type, every cap and twist tested | 23 | zero to `7e-15`: a Lagrangian subspace |

The P-type quotient keeps half the modes with both quadratures: an ordinary
classical field theory on `RP³`. The PT-type quotient keeps **every** mode
with **one** quadrature. Its admissible fields have no canonical pairs among
themselves; the conjugate partners are exactly the `η = −1` fields (pairing
smallest singular value `≥ 0.183`). Marching the leapfrog by explicit loop
confirms the identification level by level (`≤ 2.4e-14`; generic data miss
by `≥ 0.40`).

A classical state is a point of phase space: two real functions. A field on
the PT quotient is specified by one. That is the data of a polarisation —
of a wavefunction in geometric quantisation — and it is forced by the
topology, not imposed.

### 2. The deck map is a conjugation

With `j` the quarter turn of each mode plane in the sense of time evolution:
`j² = −1`, `[j, M] = 0`, and **`T j = −j T`** (`2.5e-14`), while the P-type
map commutes with `j`. Further `j L_+ = L_−`, `Ω(·, j·)` is symmetric and
definite, and the Hermitian form `g + iΩ` is real on `L_+`.

So the complex mode space is `L_+ ⊗ C`: the untwisted admissible fields are
its real vectors, the twisted fields its imaginary ones, and the imaginary
unit is the map between the two line bundles. The repository's throat
transport `J = iσ_y K` carries an antilinear `K` whose origin
`docs/finite_mouth_bundle_transport.md` traces to the identification of the
Hopf phase with the mouth's `Spin(2)`. Here antilinearity has a second,
independent source: any time-reversing deck map is antilinear for the
positive-frequency structure. Whether the two coincide is not shown.

### 3. Two caps force phase closure, and the twist picks the branch

Caps at `0` and `τ` compose to a time translation by `2τ`: compact time, a
closed history. Predicted and found, mode by mode:

    2kτ ≡ 0 (mod 2π)  if η_0 η_1 = +1,      2kτ ≡ π (mod 2π)  if η_0 η_1 = −1.

| `τ` | untwisted survivors | twisted survivors |
|---|---|---|
| `π` | all 23 | none |
| `π/2` | even `k` (11) | odd `k` (12) |
| `π/3` | `3 | k` (7) | none |
| `π/4` | `4 | k` (5) | `k ≡ 2 mod 4` (6) |
| `7π/48` | none | none |

The surviving sets match the frozen table exactly, only the product of the
twists matters, and an independent construction (the kernel of the sourced
operator on a periodic time lattice) gives the same counts.

Two ledger entries are touched:

* `history/closure.py:11` postulates *"total phase around every loop ≡ 0 or
  π (mod 2π)"*. For a linear field between two cross-caps this is a theorem,
  it is sharp, and the branch is the holonomy of the field's line bundle
  around the time loop.
* `odd_k_closure_lemma.md` Part 1 says even `k` is orientation-preserving
  closure, odd `k` orientation-reversing, and records the selection as a
  choice. At `τ = π/2` — half the antipodal light-crossing time — the
  untwisted history carries exactly the even harmonics and the twisted one
  exactly the odd. The choice is **relocated** into the twist `η_0 η_1`; it
  is not removed.

**Caveat, frozen in advance.** The harmonic index here is not shown to be
the lepton depth `k`, and the mode phase `2kτ` is not the frame holonomy
`Ω/2` of rounds 5–10. These are matches of form.

### 4. A generic closed history has no free field

At `τ = 7π/48` no homogeneous field survives, so the sourced problem has
exactly one solution (smallest singular value `0.0171`). Its Green function
is symmetric (`4.3e-15`) and deck-invariant, and the response to a source on
one time level has before/after norm ratio `1.000000`: it is exactly
time-symmetric about the source. At a resonant `τ = π/2` the problem is
solvable if and only if the source is orthogonal to the 22 kernel modes
(residual `4.7e-13` against `4.1230`).

In such a history the field is not an independent system. It is fixed by its
sources through a reciprocal, time-symmetric kernel — the direct-action
structure of Wheeler and Feynman, obtained from the boundary-value problem
rather than assumed. This supports audit item F4's global-constraint reading
and gives no support to traversable causal exchange.

## How the frozen text was met

* **P4, the retardation ratio.** The freeze predicted "of order one, pass
  threshold `> 0.1`". The measured value is exactly `1`, and for a reason
  not stated in the freeze: the periodic operator is invariant under time
  reflection about any level. The freeze also described an invariant source
  "confined to one time level"; such a source must sit on the cap, and that
  is the primary case. The off-cap single-site source is not deck-invariant
  and is reported as supplementary.
* **Explicit marching** was added as a construction check. It is not one of
  the four numbered predictions; it is grouped under P1.
* No tolerance, lattice size, cap position or table entry was changed.

## What this does and does not show

Shown, for one linear scalar on a hand-built global quotient: the PT-type
identification produces a polarisation, a conjugation, a sharp closure
spectrum and slaving of field to source. Quantum theory has all four. A
classical field on an orientable spacetime has none.

Not shown:

* **No probability.** The Born rule, the counting function `Φ`, the sector
  coefficients and the readout are untouched. Verdict fields `hbar`,
  `born_rule`, `counting_function`, `lepton_index` are `NOT_DERIVED` /
  `NOT_SELECTED` / `NOT_IDENTIFIED`.
* **No new mathematics in results 1–2.** They are the lattice form of Kay's
  F-locality (1992), Friedman (1995), Gibbons (1986), Chamblin–Gibbons
  (1995); the `T`-invariant universe of Boyle, Finn and Turok (2018) is the
  one-cap case. The classical mode space with `j` being the one-particle
  space of the quantised field is standard and is not a derivation of
  quantisation.
* **Not a localized handle.** Two mouths glued with time reversal inside one
  exterior are not a global quotient. That is the object the research record
  asks for.
* **The caps and the twist are chosen.** Nothing here selects `τ` or `η`.
* **No spinor.** With `T² = −1` the structure is quaternionic, not real, and
  no nonzero invariant section exists; that case needs its own round.

## Next

1. Repeat result 1 for a localized time-reversing handle, where the deck map
   exists only on the two-sheeted orientation cover.
2. Do the spinor case and compare its closure condition with
   `Φ_throat = π` per passage.
3. Ask whether the Lagrangian half supplies the missing measure. Two
   Lagrangian subspaces have a canonical half-density pairing; whether that
   pairing bears on the choice of `Φ` is open and is the first place this
   construction could meet the probability problem.
