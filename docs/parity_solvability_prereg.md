# Pre-registration: does constraint solvability force antipodal parity?

Publish this file before implementing or running anything. Baseline: main
`080c1cc` (PR #289 merged), which supplies `waves/reciprocal_scalar_tt.py`
and `waves/scalar_tt_constraints.py`. Seed `2026090713`.

This round follows item 5 of the COMMENT review on #289. During that review
a single measurement was made and is recorded here as **prior state, not a
finding of this round**: a mixed-parity scalar built from degrees 2 and 3
has dipole overlap `0.787` against Monte-Carlo noise `~0.002` for
parity-pure data. Nothing else about the question below has been computed.

## The question

#289 solves the linearized Hamiltonian constraint
`(Delta + 3/a^2) u = -kappa delta_rho/4` on the round `S^3`. The operator
has a genuine four-dimensional `l=1` kernel, so a solution exists only when

    P^A[phi] = int rho[phi] x^A dV = 0,   A = 0,1,2,3,

with `x^A` the four ambient coordinates restricted to `S^3`. #289 satisfies
this because its scalar is a single odd multiplet, making `rho` even.

**Is antipodal parity purity forced by that solvability, or merely
sufficient?** If forced, general relativity alone would select the parity
structure that BAM currently imposes as a boundary condition (audit item
F6). If merely sufficient, that route to F6 is closed.

"Forced" must be made precise or the question is empty, because `P^A` is
four real conditions and mixed-parity data satisfying them accidentally
will always exist. The question this round answers is therefore stated at
the level of **linear subspaces**:

> Are the two antipodal parity eigenspaces the maximal linear subspaces of
> admissible scalar data on which the dipole obstruction vanishes
> identically, for all amplitudes?

A subspace answer is the right notion because a physical sector must be
closed under superposition and rescaling. Any *isolated* mixed-parity
configuration with zero dipole is not a counterexample and must not be
reported as one.

## Analytic predictions, derived before implementation

These are pre-calculation predictions to be verified or refuted, not
numerical discoveries. If any is wrong, record the correction explicitly.

### P1. Only the parity-odd part of the source can obstruct

`x^A` is antipodally odd, so only the odd part of `rho` contributes to
`P^A`. Writing `phi = phi_e + phi_o` by antipodal parity, every term of the
inherited improved stress is quadratic, so `rho` splits into an even part
(from `phi_e` with `phi_e` and `phi_o` with `phi_o`) and an odd cross part.
Hence

    P^A[phi] = P^A_cross(phi_e, phi_o),

a **bilinear** pairing between the even and odd sectors. In particular
`P^A` vanishes identically on each parity eigenspace, which is why definite
parity is sufficient. This much restates #289.

### P2. The obstruction reduces to a degree-1 triple overlap

For `phi_n` a pure degree-`n` multiplet and `phi_n'` a pure degree-`n'`
multiplet, using `Delta phi_n = -lambda_n phi_n`,
`lambda_n = n(n+2)/a^2`, and
`grad(phi_n).grad(phi_n') = [Delta(phi_n phi_n')
+ (lambda_n + lambda_n') phi_n phi_n']/2`, the cross energy density is

    rho_cross = phidot_n phidot_n' + Delta(phi_n phi_n')/6
                + [ (lambda_n + lambda_n')/2 + 1/a^2 ] phi_n phi_n'.

Pairing with `x^A` and using `int Delta(f) x^A dV = -(3/a^2) int f x^A dV`,

    P^A = int phidot_n phidot_n' x^A dV
        + [ (n(n+2) + n'(n'+2) + 1) / (2 a^2) ] int phi_n phi_n' x^A dV.

So the entire obstruction is carried by the triple overlaps
`int Y_n Y_n' Y_1 dV`, in the field and momentum sectors separately.

### P3. The selection rule: only ADJACENT degrees obstruct

A product of harmonic polynomials of degrees `n` and `n'`, restricted to
`S^3`, decomposes into degrees `|n-n'|, |n-n'|+2, ..., n+n'`. It contains
degree 1 only if `n + n'` is odd and `|n - n'| <= 1 <= n + n'`. Since
opposite parity forces `|n - n'|` odd, this requires

    |n - n'| = 1.

**Prediction: the dipole obstruction couples only adjacent multiplets.**
Equal degrees never obstruct (`n + n'` even), which is P1 again.

### P4. The predicted answer is NEGATIVE

If P3 holds, definite antipodal parity is **sufficient but not necessary**.
The obstruction vanishes identically on any span of pure multiplets whose
degree set contains no adjacent pair. Definite parity is one such condition
— degrees of one parity differ by at least 2 — but strictly stronger than
needed. Explicit predicted counterexamples, both mixed parity:

    degrees {1, 4},   degrees {2, 5},   degrees {3, 6}.

Each should give an identically zero dipole at every amplitude, while
`{1, 2}`, `{2, 3}` and `{3, 4}` should not.

**The expected outcome of this round is therefore that constraint
solvability does NOT force antipodal parity, closing this route to F6.**
Recording that expectation in advance is deliberate: a negative result is
the anticipated result and is fully acceptable. It must not be softened if
it holds, and P3 must not be quietly widened if it fails.

### P5. What the parity sectors are, stated as interpretation only

Even harmonics descend to the antipodally identified `S^3/Z_2`; odd
harmonics are sections of the twisted real line bundle over it. So the two
parity eigenspaces are exactly the two `Z_2` sectors of the identified
space. This is an observation about which data are single-valued after
identification. It is **not** a derivation of the antipodal boundary
condition, and this round may not present it as one.

## What this round does not address

The nonlinear constraint, orders beyond `phi^2`, the supporting matter that
#289 sets to zero, evolution and constraint propagation, and any
non-subspace (variety) characterization of the zero set. A mixed-parity
configuration that satisfies `P^A = 0` by accident is expected to exist and
is not evidence either way. No counting function, Born law, triangle map or
readout follows from anything here.

## Frozen checks

Dimensionless `a = kappa = 1`. Use the existing `harmonic_multiplet`
machinery and its analytic moments; do not introduce a new quadrature rule
where an exact moment is available. Seed `2026090713` for random data.

1. **Bilinearity (P1).** For 200 random `(phi_e, phi_o)` pairs, the even-even
   and odd-odd contributions to `P^A` vanish to `1e-12`, and `P^A` equals
   its cross part to `1e-12`.
2. **Reduction (P2).** For all degree pairs with `n, n' <= 6`, the closed
   coefficient `(n(n+2) + n'(n'+2) + 1)/2` reproduces `P^A` computed from
   the full inherited improved stress to relative `1e-10`, at nonzero
   momentum as well as `p = 0`.
3. **Selection rule (P3).** Compute `int Y_n Y_n' Y_1 dV` for all
   `1 <= n, n' <= 6` and report the full table. Adjacent pairs must be
   nonzero by at least `1e-3` in operator norm; every non-adjacent pair must
   vanish to `1e-12`. Both must be checked; a table of zeros alone is not
   evidence of a selection rule.
4. **Counterexamples (P4).** For degree sets `{1,4}`, `{2,5}`, `{3,6}`,
   200 random mixed-parity states at three amplitudes each give `|P^A|`
   below `1e-12`. For `{1,2}`, `{2,3}`, `{3,4}`, the same sampling gives a
   maximum `|P^A|` above `1e-3`, and halving the amplitude of one parity
   component halves `|P^A|` to `1e-8` — confirming the bilinearity rather
   than a generic nonzero number.
5. **Maximality.** For the truncation `V_1 + V_2`, verify that no nonzero
   `u` in `V_2` is annihilated by all of `V_1`, so no subspace
   `V_1 + U` with `U` nonzero evades the obstruction. Report whether any
   graph subspace `{v + Tv}` with `T: V_1 -> V_2` nonzero has an identically
   vanishing dipole; search by solving the symmetrized quadratic conditions,
   and report `UNRESOLVED` rather than asserting absence if the search is
   inconclusive.
6. **Momentum sector, secondary.** Report whether the six Killing charges
   and the four gradient conformal Killing charges obey the same
   adjacent-degree rule. No prediction is frozen for this; report what is
   measured, with its own verdict field.
7. A failed required check yields `UNRESOLVED` on every verdict field, names
   the failures, writes valid JSON, overwrites a stale passing report, and
   exits nonzero.

## Verdicts

Separate fields, never one label:

- `dipole_obstruction_structure`: `BILINEAR_CROSS_PARITY_ADJACENT_DEGREE_ONLY`,
  or `UNRESOLVED`.
- `antipodal_parity_status`: `SUFFICIENT_NOT_NECESSARY`, `FORCED_AT_SUBSPACE_LEVEL`,
  or `UNRESOLVED`.
- `f6_consequence`: `CONSTRAINT_SOLVABILITY_DOES_NOT_DERIVE_THE_ANTIPODAL_CONDITION`,
  `SUPPORTS_A_DERIVATION`, or `UNRESOLVED`.
- `momentum_sector`: reported, with `NOT_PREDICTED_IN_ADVANCE` recorded.
- `triangle_map`, `readout`: `NOT_DERIVED`.

## Deliverables

One reusable module, independent tests, an archived report, a write-up
separating prior results from new ones, and an update to
`docs/qft_emergence_audit.md` recording the F6 consequence. Preserve every
earlier freeze and archive. Pin this commit in the implementation.
