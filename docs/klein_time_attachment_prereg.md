# Pre-registration: what a time-nonorientable (Klein) attachment does to a classical field

Publish this file before implementing or running anything. Baseline: main
`3293d5b` (PR #350 merged). No random seed is needed for the frozen checks;
where random vectors are used as generic controls the seed is `2026111001`.

**Prior state.** Nothing below has been computed. The predictions are
analytic, derived in this file before any code exists. This freeze and its
implementation are written in one session by one author; the only ordering
evidence is that this commit precedes the code.

## Why this round

Three items are open in the repository's own ledger:

* the v1.3.5 research record lists "establish the Klein/time attachment" as a
  requirement for the next model;
* `history/closure.py:11` states phase closure "≡ 0 or π (mod 2π)" as an
  axiom, and rounds 5–8 record the conditioning on phase as **chosen**;
* `docs/odd_k_closure_lemma.md` records the even/odd-`k` split as a **choice
  of boundary condition**;
* audit item F4 leaves traversable causal exchange unsupported and develops
  the global-constraint reading instead.

`docs/gr_handle_orientation.md` already separates the time-preserving
antipodal map `J_plus` from the time-reversing `J_minus` on Kruskal. This
round asks what the time-reversing identification does to a **classical real
field**, on the repository's standard closed background, with no quantum
input.

## The model

Einstein static universe `R × S³`, radius `a = 1`, conformally coupled real
scalar. Restrict to zonal fields about a pole, `φ = u(t,χ)/sin χ`. Then

    u_tt = u_χχ,   u(t,0) = u(t,π) = 0,

exactly, with modes `sin(kχ)`, `ω_k = k`, `k = n + 1 ≥ 1`. The antipode is
`χ ↦ π − χ` and `sin(k(π−χ)) = (−1)^{k+1} sin(kχ)`.

Deck maps on the cover, with twist `η = ±1` (the field is a section of the
trivial or the nontrivial real line bundle):

    orientable control:   (P φ)(t,χ)   = φ(t, π−χ)
    Klein cap at time c:  (T_c φ)(t,χ) = η φ(2c − t, π−χ)

`T_c` is a free involutive isometry; the quotient is time-nonorientable and
contains the one-sided spacelike hypersurface `{c} × RP³` (a cross-cap in
time). A perfect fluid at rest has `T_ab` even under `u^a ↦ −u^a`, so the
quotient is locally the same Einstein solution as the cover. Nothing in this
round derives the caps; they are placed by hand.

**Discretisation.** Lattice `χ_i = iπ/N`, `dt = dχ`, leapfrog
`u^{j+1} = L u^j − u^{j−1}` with `L` the nearest-neighbour sum. At unit
Courant number this is exact d'Alembert propagation and its discrete
frequencies are exactly `ω_k = k`. State `s_j = (u^j, u^{j+1})`, one-step
matrix `M`, conserved Wronskian `Ω(s,s') = u^j·u'^{j+1} − u^{j+1}·u'^j`.
A cap at `c = m·dt/2` (integer `m`, so caps may sit on or between time
levels) acts on states as `T_m = η F S M^{m−1}` with `F` the spatial flip
and `S` the swap of the two time levels. Frozen size `N = 24`
(`M_dim = 23` interior points, phase-space dimension `46`).

## Analytic predictions

### P1. Admissible fields form a Lagrangian half; the orientable quotient does not

`T` reverses time, so `T*Ω = −Ω`. For `a, b` fixed by `T`,
`Ω(a,b) = Ω(Ta,Tb) = −Ω(a,b) = 0`. Both eigenspaces `L_± = ker(T ∓ 1)` are
isotropic and they span, so each has exactly half the dimension.

Predicted, for caps `m ∈ {0, 1, 7, 24}` and both twists:

* `T_m² = 1`, `T_mᵀ Ω T_m = −Ω`;
* `dim L_+ = dim L_− = 23`;
* `Ω` restricted to `L_+` and to `L_−` vanishes;
* the pairing `Ω : L_+ × L_− → R` is nondegenerate.

Control. `P` preserves time: `PᵀΩP = +Ω`, and its fixed space is a
**symplectic** subspace (restricted `Ω` of full rank). For `N = 24` the
flip-even interior vectors have dimension `12`, so the fixed space has
dimension `24` and restricted rank `24`.

The distinction to be reported is structural, not a dimension count: the
orientable quotient keeps half the modes with both quadratures; the Klein
quotient keeps every mode with one quadrature.

### P2. The time deck is a conjugation; admissible plus twisted is the complex mode space

Let `j` be the real operator that rotates each mode plane by a quarter turn
in the sense of time evolution (`M = exp(dt·A)`, `j = A|A|^{-1}`). Predicted:

* `j² = −1`, `[j, M] = 0`;
* `T j = −j T` (the deck map is antilinear for the positive-frequency
  complex structure);
* `j L_+ = L_−`;
* `g(a,b) = Ω(a, j b)` is symmetric and definite, and the Hermitian form
  `h = g + iΩ` is **real** on `L_+`.

So `(V, j) = L_+ ⊗ C` canonically, with `T` as complex conjugation: the
untwisted admissible fields are the real vectors, the `η = −1` fields the
imaginary ones.

### P3. Two caps force phase closure, with the branch fixed by the twist

Caps at `0` and `τ`, twists `(η_0, η_1)`. `T_τ T_0` is translation by `2τ`
times `η_0 η_1`, so the quotient has compact time and is a closed history.
Mode `k` with spatial parity `p_k`: cap 0 selects `cos(kt − α)` with
`2α = 0` if `η_0 p_k = +1`, else `π`; cap `τ` then requires

    2kτ ≡ 0  (mod 2π)   if η_0 η_1 = +1,
    2kτ ≡ π  (mod 2π)   if η_0 η_1 = −1.

The phase accumulated around the time loop is `0` or `π` mod `2π`, sharply,
and which one is the holonomy of the field's line bundle around that loop.

Predicted dimension of the joint admissible space, `N = 24`, `τ = mπ/48`:

| `m` | `τ` | `(+,+)` | `(+,−)` |
|---:|---|---:|---:|
| 48 | `π` | 23 (all `k`) | 0 |
| 24 | `π/2` | 11 (even `k`) | 12 (odd `k`) |
| 16 | `π/3` | 7 (`3 | k`) | 0 |
| 12 | `π/4` | 5 (`4 | k`) | 6 (`k ≡ 2 mod 4`) |
| 7 | `7π/48` | 0 | 0 |

`(−,−)` must equal `(+,+)` and `(−,+)` must equal `(+,−)`. The surviving
modes must be exactly the listed `k`, not only the right number of them.

Corollary to be reported with its caveat: at `τ = π/2`, half the antipodal
light-crossing time, the untwisted closed history carries exactly the even
harmonics and the twisted one exactly the odd harmonics. This is the same
`Z_2` rule as `odd_k_closure_lemma.md` Part 1 (accumulated phase `kπ` must
match the loop's orientation class). **It is a structural match only**: the
identification of the harmonic index with the lepton depth `k` is not made
here and must not be claimed.

### P4. A generic closed history has no free field: sources determine it

Write the sourced leapfrog on the periodic time lattice (untwisted,
`m` steps): `D u = s`, `D = (shift_+ + shift_−) ⊗ 1 − 1 ⊗ L`, symmetric.

* Independent route to P3: `dim ker D = 2 × (P3 count)` and the
  `T_0`-invariant part of `ker D` has dimension equal to the P3 count, for
  every row of the table.
* For `m = 7`: `D` is invertible. Every invariant source has exactly one
  admissible field; `G = D^{-1}` is symmetric (emitter–absorber reciprocity)
  and commutes with `T_0`.
* That unique field is **not retarded**: for a source confined to one time
  level, the field is nonzero on the levels before it. Report the ratio of
  field norm on the three levels before the source to the three after;
  predicted to be of order one, pass threshold `> 0.1`.
* For `m = 24` (resonant): the sourced problem is solvable if and only if
  the source is orthogonal to `ker D` (Fredholm alternative).

## Verdict fields and failure rule

    phase_space        ADMISSIBLE_FIELDS_ARE_A_LAGRANGIAN_HALF
    complex_structure  TIME_DECK_IS_A_CONJUGATION_OF_THE_MODE_SPACE
    phase_closure      TWO_CAPS_FORCE_PHASE_0_OR_PI_FIXED_BY_TWIST
    field_freedom      GENERIC_CLOSED_HISTORY_HAS_NO_FREE_FIELD
    hbar               NOT_DERIVED
    born_rule          NOT_DERIVED
    counting_function  NOT_SELECTED
    lepton_index       NOT_IDENTIFIED

If any required check fails, the corresponding field is `REFUTED` with the
measured numbers and the write-up says so in its first paragraph. A missing
check counts as a failure. Tolerance for "vanishes": `1e-9` in the natural
norm; for "nondegenerate"/"definite": smallest singular value or eigenvalue
`> 1e-6`. Ranks use relative threshold `1e-9`.

## What may not be concluded

* None of this is quantum mechanics. A Lagrangian subspace, a conjugation
  and a discrete closure spectrum are kinematic structures quantum theory
  also has; no `ħ`, no probability and no measurement rule follows.
* The classical solution space of a linear field with its
  positive-frequency complex structure *is* the one-particle space of the
  quantised field; that correspondence is standard and is not a derivation
  of quantisation.
* The mathematics of P1–P2 is known: Kay's F-locality (1992), Friedman's
  Lorentzian universes-from-nothing (1995), Gibbons (1986) and
  Chamblin–Gibbons (1995) on time-nonorientability and the complex
  structure, Boyle–Finn–Turok (2018) on the `T`-invariant universe. The
  contribution here is to compute it on this repository's background and
  connect it to three items its ledger records as chosen.
* A global quotient is not a localized Klein handle. Two mouths glued with
  time reversal inside one exterior are not treated.
* One linear zonal scalar is not the BAM field content, and the spinor case
  (`T² = −1`, quaternionic rather than real structure) is not treated.
* The frame-holonomy phase `Ω/2` of rounds 5–10 is not the mode phase `2kτ`
  of P3. A shared "0 or π" form is not an identification.
