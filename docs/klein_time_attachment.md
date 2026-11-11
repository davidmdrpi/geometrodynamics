# The Klein time attachment: complementary quadratures and a twist-dependent resonance

**Two results stand and one interpretation is withdrawn.** A classical real
field on a time-nonorientable quotient of the Einstein static universe keeps
one quadrature of every mode — a Lagrangian half of the cover's phase space —
and two such attachments compose to a resonance condition whose branch is set
by the field's twist. The frozen claim that a generic closed history has **no
free field** is false in the continuum: the review of `5b79975` supplied an
explicit counterexample, reproduced here. It was an artefact of the frozen
cutoff. The time-symmetric response is also not a consequence of the
attachment. None of this is a probability rule.

*Pre-registered in [klein_time_attachment_prereg.md](klein_time_attachment_prereg.md)
(`bb6ff2c`), committed before the module and left unchanged. Reviewed at
`5b79975`; corrections R1–R5 below. Module
`geometrodynamics/bulk/klein_time.py`; probe `klein_time_probe.py`; tests
`tests/test_klein_time.py` (80). Archives: the run as first written in
`experiments/closure_ledger/runs/20261110_klein_time/`, the corrected run in
`runs/20261110_klein_time_review/`.*

```bash
python -m experiments.closure_ledger.klein_time_probe   # 28/28
python -m pytest -q tests/test_klein_time.py             # 80 tests
```

| verdict field | value |
|---|---|
| phase_space | `ADMISSIBLE_FIELDS_ARE_A_LAGRANGIAN_HALF` |
| complex_structure | `TIME_DECK_IS_A_CONJUGATION_OF_THE_COVER_MODE_SPACE` (narrowed) |
| phase_closure | `TWO_CAPS_FORCE_PHASE_0_OR_PI_FIXED_BY_TWIST` |
| field_freedom | `UNIQUE_SOURCED_SOLUTION_IN_THE_FROZEN_TRUNCATION` (narrowed) |
| continuum_free_field | `COUNTEREXAMPLE` |
| irrational_spacing | `OPEN_SMALL_DIVISORS` |
| attachment_contribution | `QUADRATURE_SELECTION_NOT_THE_TIME_SYMMETRIC_RESPONSE` |
| sector_pairing | `NOT_SUPPLIED` |
| hbar, born_rule | `NOT_DERIVED` |
| counting_function, lepton_index | `NOT_SELECTED`, `NOT_IDENTIFIED` |

## The construction

`docs/finite_mouth_bundle_transport.md` §4 records that PR #129's antipodal
map is PT-type — it reverses time orientation — while the alternative P-type
map does not, and adds: *"Which is physical is a further choice, recorded as
the datum `T`. For the static problem the two agree."* This round does the
dynamic problem, where they do not agree.

Background: `R × S³`, radius 1, conformally coupled real scalar, zonal
sector `φ = u(t,χ)/sin χ`. Then `u_tt = u_χχ` with `u = 0` at both poles,
exactly; modes `sin(kχ)`, `ω_k = k`; antipode `χ ↦ π − χ`. On the lattice
`dt = dχ = π/N` the leapfrog rule is exact propagation with exactly these
integer frequencies, for the harmonics `k < N` the lattice carries.

    P-type (control):   (P φ)(t,χ)   = φ(t, π−χ)
    PT-type cap at c:   (T_c φ)(t,χ) = η φ(2c − t, π−χ),   η = ±1

The PT quotient contains a one-sided spacelike `RP³` at `t = c`: a cross-cap
in time. It is locally the same Einstein solution as the cover, since a
fluid at rest has a stress tensor even under `u^a ↦ −u^a`.

## What stands

### 1. Complementary quadratures: a Lagrangian half

`T` reverses the conserved Wronskian, so its fixed fields satisfy
`Ω(a,b) = −Ω(a,b) = 0`.

| quotient | fixed dimension (of 46) | Wronskian on the fixed space |
|---|---:|---|
| P-type | 24 | rank 24: a symplectic subspace |
| PT-type, every cap and twist tested | 23 | zero to `7e-15`: a Lagrangian subspace |

The P-type quotient keeps half the modes with both quadratures. The PT-type
quotient keeps every mode with one quadrature, and the twist decides which.
The untwisted and twisted sectors are complementary: the Wronskian vanishes
inside each and pairs them with each other (smallest singular value
`≥ 0.183`). Marching the leapfrog by explicit loop confirms the
identification level by level (`≤ 2.4e-14`; generic data miss by `≥ 0.40`).

Both quotients carry the same amount of data on their own spatial section:
on the cap `RP³` the PT field is an even function with an odd normal
derivative, the P field an even function with an even one. What the PT
quotient lacks is not data but a bracket. The Wronskian current is odd under
the deck map, so its flux through a one-sided surface is not defined, and no
symplectic form descends.

### 2. Two attachments compose to a twist-dependent resonance

Caps at `0` and `τ` compose to a time translation by `2τ`. Mode `k`
survives only if

    2kτ ≡ 0 (mod 2π)  if η_0 η_1 = +1,      2kτ ≡ π (mod 2π)  if η_0 η_1 = −1.

| `τ` | harmonics, untwisted / twisted (`k < 24`) | orientable closed time | two caps |
|---|---|---:|---:|
| `π` | all 23 / none | 46 / 0 | 23 / 0 |
| `π/2` | even `k` / odd `k` | 22 / 24 | 11 / 12 |
| `π/3` | `3 | k` / none | 14 / 0 | 7 / 0 |
| `π/4` | `4 | k` / `k ≡ 2 mod 4` | 10 / 12 | 5 / 6 |
| `7π/48` | none / none below `k = 24` | 0 / 0 | 0 / 0 |

The surviving sets match the frozen table, and only the product of the
twists matters. The last two columns are dimensions and are the new control
(R2): a **time-orientable** closed time of the same period, periodic or
antiperiodic, selects exactly the same harmonics with both quadratures. So
the resonance condition belongs to temporal closure. What the two caps add
is that closure is produced by composing them rather than imposed, and that
each surviving mode keeps one quadrature instead of two.

Two ledger entries are touched, with that qualification:

* `history/closure.py:11` postulates *"total phase around every loop ≡ 0 or
  π (mod 2π)"*. For a linear field in a closed history this is a theorem,
  and the branch is the holonomy of the field's line bundle around the time
  loop.
* `odd_k_closure_lemma.md` Part 1 records the even/odd-`k` selection as a
  choice. At `τ = π/2` the untwisted history carries exactly the even
  harmonics and the twisted one exactly the odd. The choice is **relocated**
  into the twist; it is not removed.

**Caveat, frozen in advance.** The harmonic index here is not shown to be
the lepton depth `k`, and the mode phase `2kτ` is not the frame holonomy
`Ω/2` of rounds 5–10. These are matches of form.

## What was corrected

### R1 (high). "No free field" was a property of the cutoff

The freeze took `τ = 7π/48` as its generic case and predicted no surviving
mode. That holds for `k < 24` only. In the continuum, with `τ/π = p/q` in
lowest terms, the untwisted condition is `q | k` and the twisted condition
is `q` even with `k` an odd multiple of `q/2`. Holding the physical cap
separation fixed and refining:

| `N` | untwisted survivors | twisted survivors |
|---:|---|---|
| 24 | none | none |
| 48 | none | 24 |
| 72 | 48 | 24 |
| 96 | 48 | 24, 72 |

This is the reviewer's table, reproduced. The explicit continuum mode

    u(t,χ) = sin(48χ) sin(48t)

satisfies both untwisted caps (defects `4.9e-14`, `4.7e-14`, by marching on
`N = 96`), and `sin(24χ) sin(24t)` satisfies the twisted pair. More
generally **every rational spacing admits the untwisted mode `k = q`**
(checked for `1/3, 2/5, 3/7, 5/11, 4/9`). So at rational spacing a
compatible source keeps homogeneous freedom and an incompatible one has no
solution.

The frozen label `GENERIC_CLOSED_HISTORY_HAS_NO_FREE_FIELD` is withdrawn.
What the frozen checks establish is `UNIQUE_SOURCED_SOLUTION_IN_THE_FROZEN_TRUNCATION`.

**Irrational spacing is open.** No mode satisfies the closure condition
exactly, so a solution is unique if it exists. Existence is a small-divisor
question. The continuum operator is diagonal with symbol `(q/α)² − k²`,
`α = τ/π`, and its smallest value near the diagonal is about
`2k·|kα − q|/α`. A diagnostic over `k ≤ 10⁵` (not a required check):

| `τ/π` | smallest gap | at `k` |
|---|---:|---:|
| `7/48` | 0 | 48 |
| golden ratio conjugate | 1.0000 | 1 |
| `√2 − 1` | 1.0000 | 1 |
| `π − 3` | 0.0481 | 113 |

*Corrected after review of `c433f91`:* the diagnostic first rounded `kα` to
the nearest integer, which is not the minimiser of a difference of squares.
It now tries both neighbouring integers, zero included. The golden-ratio
entry was `1.3820` at `k = 2`; the true minimum is `1` at `k = 1`, `q = 0`.
The other entries are unchanged. Corrected run:
`runs/20261111_klein_time_review2/`.

The gap stays bounded away from zero exactly when `α` is badly
approximable, which is a set of measure zero; for almost every spacing it is
expected to come arbitrarily close to zero. That is a statement about the
diagonal continuum operator from standard Diophantine approximation. It is
not an existence theory, it says nothing about which source regularity
restores solvability, and it is not the lattice limit.

### R2 (medium). The time-symmetric response is not the attachment's doing

The symmetric operator and its inverse exist before any deck constraint is
imposed. On the same period with the same operator, a single-site source
that is **not** deck-invariant already has a unique response, with a
symmetric Green function and a response time-symmetric about the source to
`1.8e-15`. The first write-up's sentence that this is "the direct-action
structure of Wheeler and Feynman, obtained from the boundary-value problem"
attributed to the attachment something temporal closure supplies alone. It
is withdrawn, and so is the claim of support for audit item F4.

What the attachment does add, measured against the orientable control:

| | orientable closed time | with the cap constraint |
|---|---:|---:|
| homogeneous dimension at `τ = π/2` | 22 | 11 |
| admissible source dimension, `m = 7` | 161 | 81 |
| solvability conditions at resonance | 22 | 11 |

It halves each: one quadrature per mode, deck-invariant sources only, and
half the Fredholm conditions, the other half holding automatically.

### R3 (medium). The complex structure belongs to the cover

`T j = −j T` and `j L_+ = L_−` hold (`≤ 4.6e-14`), with `Ω(·, j·)`
symmetric and definite. But `j L_+ = L_−` means `j` maps **no** nonzero
admissible field of a fixed twist to an admissible field of that twist
(`dim(L_+ ∩ jL_+) = 0`, every cap and twist tested). Time translation moves
the cap, so it does not act within a sector either.

So the statement is: putting the two twist sectors together reconstructs the
cover's complex mode space, on which `T` is a conjugation. Neither sector
carries a complex structure from this construction, and nothing here
justifies combining the two bundles physically. The first write-up's
sentence that "the complex mode space is `L_+ ⊗ C`" was true as linear
algebra and misleading as physics; the label is narrowed to
`TIME_DECK_IS_A_CONJUGATION_OF_THE_COVER_MODE_SPACE`.

This matches a known obstruction rather than evading it. Friedman and
Higuchi (Phys. Rev. D 52, 5687, 1995) show that on time-nonorientable
spacetimes of this kind local quantum states do not extend consistently to a
global one. A route from this construction to quantum structure has to solve
that problem, not cite the classical correspondence around it.

### R4 (medium). A missed tolerance is not a refutation

The freeze's failure rule labelled a missing check `REFUTED`. The verdict
now separates `INCOMPLETE` (a required check was not run),
`UNRESOLVED_NUMERICALLY` (a check ran and missed its tolerance) and
`COUNTEREXAMPLE` (an explicit counterexample is recorded). Only the last is
a mathematical statement. The frozen file is unchanged; this is a recorded
departure from its rule.

### R5 (found while correcting). An overstated data count and a citation

* The first write-up said a PT field "is specified by one" real function
  where a classical state needs two, "the data of a wavefunction". That
  compared the quotient with its cover. On its own spatial section the PT
  quotient has as much data as the P quotient (result 1). The sentence is
  withdrawn.
* "Friedman (1995)" in the freeze and first write-up is Friedman and
  Higuchi, as cited in R3; that reference was checked against its arXiv
  record. The freeze's other citations were written from memory. The
  attribution to Kay (1992) in particular is unverified and is dropped here.

### As first reported

The retardation ratio was frozen as "of order one"; it is exactly `1`, by
the time-reflection symmetry R2 identifies. Explicit marching was added as a
construction check under P1. No tolerance, lattice size, cap position or
frozen table entry has been changed, before or after review.

## What this does and does not show

Shown, for one linear scalar on a hand-built global quotient: a
time-reversing attachment constrains complementary quadratures, and two of
them compose to a twist-dependent resonance.

Not shown:

* **No probability.** The Born rule, the counting function `Φ`, the sector
  coefficients and the readout are untouched.
* **No pairing between the sectors is supplied.** The Wronskian pairs the
  untwisted and twisted sectors on the cover; nothing physical here couples
  them. `docs/conjugate_boundary_pairing_protocol.md` (PR #351) is explicit
  that a symplectic form is not a positive norm and that the variational
  pairing, the conjugate-mouth transport and a detector response must not be
  conflated. That applies to this round.
* **No new mathematics in result 1.** It restates on a lattice what is
  known for time-nonorientable spacetimes (Gibbons 1986; Chamblin–Gibbons
  1995; Friedman–Higuchi 1995); the `T`-invariant universe of Boyle, Finn
  and Turok (2018) is the one-cap case.
* **Not a localized handle.** Two mouths glued with time reversal inside one
  exterior are not a global quotient.
* **The caps and the twist are chosen.** Nothing here selects `τ` or `η`.
* **No spinor.** With `T² = −1` the structure is quaternionic, not real.

## Next

The question this round leaves is the one the review names and PR #351's
protocol is built to answer: does a localized, dynamically responding mouth
pair preserve the complementary-quadrature structure, and does its physical
response supply a pairing between the two sectors? Concretely, for that
protocol:

1. Test whether a time-reversing seam on the two-sheeted orientation cover
   of a handle leaves a Lagrangian half, as the global quotient does.
2. Ask whether the mouths' boundary response couples the untwisted and
   twisted quadratures, and if so whether that coupling is the variational
   pairing or only a coordinate one (`COORDINATE_PAIRING_ONLY` in the
   protocol's vocabulary).
3. Settle the irrational-spacing problem before using any closed-history
   uniqueness statement.
