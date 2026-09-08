# Does constraint solvability force antipodal parity?

**No.** It is sufficient and strictly stronger than necessary.

Pre-registered in [`parity_solvability_prereg.md`](parity_solvability_prereg.md)
at [`495f1f1`](https://github.com/davidmdrpi/geometrodynamics/commit/495f1f185411a80cfd688b0b33291d1986f9702a),
committed and published before any of the code or numbers below. Baseline is
main `080c1cc` (PR #289 merged). This follows item 5 of the COMMENT review on
#289.

| verdict field | value |
|---|---|
| dipole obstruction structure | `BILINEAR_CROSS_PARITY_ADJACENT_DEGREE_ONLY` |
| antipodal parity status | `SUFFICIENT_NOT_NECESSARY` |
| F6 consequence | `CONSTRAINT_SOLVABILITY_DOES_NOT_DERIVE_THE_ANTIPODAL_CONDITION` |
| momentum sector | Killing charge is independent of parity (not predicted in advance) |
| triangle map, readout | `NOT_DERIVED` |

## The question and why it needed sharpening

#289 solves `(Delta + 3/a^2) u = -kappa delta_rho/4` on the round `S^3`. That
operator has a genuine four-dimensional `l=1` kernel, so a solution exists
only if `P^A = int rho x^A dV = 0`. #289 satisfies this because its scalar is
a single odd multiplet, making `rho` even. The review asked whether general
relativity was therefore *deriving* the antipodal parity that BAM imposes as
a boundary condition (audit item F6).

Asked naively the question is empty: `P^A` is four real conditions, so
mixed-parity data satisfying them accidentally always exists. The freeze
therefore posed it at the level of **linear subspaces** — a physical sector
must be closed under superposition and rescaling — and asked whether the two
parity eigenspaces are the maximal such subspaces.

## What the obstruction actually is

Every term of the inherited improved stress is quadratic, and `x^A` is
antipodally odd, so only the parity-odd cross part of `rho` survives the
pairing. Writing `phi = phi_e + phi_o`, the obstruction is a **bilinear**
pairing between the two parity sectors. That is why definite parity works.

For pure multiplets of degrees `n` and `n'`, using
`grad(phi_n).grad(phi_n') = [Delta(phi_n phi_n') + (lambda_n + lambda_n')
phi_n phi_n']/2` and `int Delta(f) x^A dV = -(3/a^2) int f x^A dV`,

    P^A = int phidot_n phidot_n' x^A dV
        + [ (n(n+2) + n'(n'+2) + 1) / (2 a^2) ] int phi_n phi_n' x^A dV.

The whole obstruction is carried by the degree-1 triple overlap
`int Y_n Y_n' Y_1 dV`. A product of harmonic polynomials of degrees `n` and
`n'` decomposes into degrees `|n-n'|, |n-n'|+2, ..., n+n'`, so it contains
degree 1 only when `n + n'` is odd and `|n - n'| <= 1`. Opposite parity makes
`|n - n'|` odd, so

> **the dipole obstruction couples only adjacent multiplets, `|n - n'| = 1`.**

Measured over every pair up to degree 6: the minimum adjacent overlap norm is
`3.464102`, the maximum non-adjacent norm is `1.39e-14`.

## The answer

Definite antipodal parity is sufficient because degrees of one parity differ
by at least 2, so no adjacent pair exists. But that is strictly stronger than
required: **any** degree set without an adjacent pair is equally unobstructed,
including mixed-parity ones.

| degrees | mixed parity | adjacent pair | max \|P^A\| over 600 samples |
|---|---|---|---:|
| `{1, 4}` | yes | no | `5.62e-15` |
| `{2, 5}` | yes | no | `1.52e-14` |
| `{3, 6}` | yes | no | `2.86e-14` |
| `{1, 2}` | yes | yes | `5.076` |
| `{2, 3}` | yes | yes | `6.933` |
| `{3, 4}` | yes | yes | `9.032` |
| `{1, 3}`, `{2, 4}`, `{1,3,5}`, `{2,4,6}` | no | no | `0` exactly |

So constraint solvability does not derive the antipodal condition. **This
route to F6 is closed.** The freeze recorded that expectation in advance
precisely so it could not be softened once it held.

**The parity eigenspaces are not maximal even within an adjacent truncation.**
Two narrower statements do hold, and are retained: for `V_n + V_{n+1}` no
nonzero `u` in `V_{n+1}` is annihilated by all of `V_n` (smallest singular
value `6.928203` at `n=1`), and no nonzero `T: V_n -> V_{n+1}` gives a graph
subspace with identically vanishing dipole (nullity `0` on a `40 x 36`
system). But both families enlarge or map *all* of `V_n`, so neither contains
a subspace with a smaller projection — and such subspaces evade:

    S = span{ x_0, x_1 x_2 }  in  V_1 + V_2

is mixed parity, sits inside an adjacent pair, and has
`int x^A x_0 x_1 x_2 dV = 0` for every ambient coordinate. Its dipole vanishes
at every amplitude on both routes (`4.9e-15` and `3.3e-14` over 200 random
amplitudes), against `4.67` for generic data in the same adjacent pair. So the
parity answer fails locally as well as globally, which strengthens the
negative conclusion rather than weakening it.

## The momentum sector is not a parity condition at all

No prediction was frozen here. Measured: parity-pure data always kills the
Hamiltonian dipole exactly, yet carries large `SO(4)` Killing charges.

| degree | parity | Hamiltonian dipole | Killing charge |
|---:|---|---:|---:|
| 1 | odd | `0` | `1.000000` |
| 2 | even | `0` | `1.071142` |
| 3 | odd | `0` | `3.745197` |
| 4 | even | `0` | `6.992044` |

Completing the audit to all six Killing and all four gradient conformal
charges exposes their structure. **Killing charges are diagonal in degree**: a
Killing field is divergence free, the improved-stress correction integrates
away against it, and the canonical charge `-p^T G q` needs field and momentum
in the same multiplet. **Gradient conformal charges are not**: `grad(x^A)` has
divergence `-lambda_1 x^A`, the improvement survives, and the charge is a
field-momentum pairing obeying the same adjacency rule as the Hamiltonian
dipole. With field in `V_2` and momentum in `V_3` all six Killing charges and
the Hamiltonian dipole vanish exactly while a gradient charge reaches `0.267`;
at equal or non-adjacent degrees the gradient charges vanish (`2.9e-17`,
`4.4e-15`).

So there are three independent obstructions, and none of them is a parity
condition.

## Interpretation, offered as observation only

Even harmonics descend to the antipodally identified `S^3/Z_2`; odd harmonics
are sections of the twisted real line bundle over it. The parity eigenspaces
are therefore exactly the two `Z_2` sectors of the identified space, and data
that is single-valued after identification is automatically unobstructed.
That is a statement about which data descend, **not** a derivation of the
antipodal boundary condition, and the result above shows general relativity
does not force the choice.

## Implementation correction

**C1.** The overlap route summed ordered pairs `(n,m)` and `(m,n)`, each with
the full frozen `P2` coefficient, double counting the total cross term. It was
caught by the frozen improved-stress cross-check, which disagreed by exactly
`2.000` in every case while agreeing exactly in direction. One overall factor
of `1/2` fixes it, and the two routes then agree to `~1e-15`. No frozen
prediction changed: a global factor cannot move a zero, so the selection rule
and both counterexample classes are unaffected and only reported magnitudes
halve. The bilinearity check was also moved onto the independent stress route,
because measured through the overlap route it is bilinear by construction and
could not fail.

## Scope

Linearized Hamiltonian constraint at order `phi^2` on the round `S^3`, zero
supporting-matter perturbation, linear subspaces of pure-multiplet data. Not
the nonlinear constraint, not evolution, not a variety-level characterization
of the zero set. Isolated mixed-parity configurations with accidental zero
dipole exist and are not evidence either way. No counting function, Born law,
triangle map or readout follows.

## Reproduce

```bash
python -m experiments.closure_ledger.parity_solvability_probe   # 11 required checks
python -m pytest -q tests/test_parity_solvability.py             # 61 tests
```


## Review corrections (C2-C6)

Recorded, not applied silently. None changes the adjacency selection rule or
the counterexample classes; C4 strengthens the negative conclusion.

**C2 - a stale factor of two.** `subspace_maximality` retained the `2 *` that
C1 removed from the dipole normalization, so its reported smallest singular
value was `13.856406` where the consistent value is `6.928203`. It rescales
every singular value and cannot change the kernel, so the conclusion stands.

**C3 - a missing check could not fail the verdict.** `verdict` gated only on
`all(checks.values())`, so `verdict({"unrelated": True}, ...)` returned the
full affirmative result and a probe that silently dropped every real check
would still have passed. A test asserting `verdict({"a": True}, ...)` had
locked that in. The frozen required names are now declared in
`REQUIRED_CHECKS`, presence is required as well as truth, the probe asserts
its check set matches exactly, and a test drops each required check in turn
and drives the failing CLI end to end.

**C4 - the "fails only globally" claim was false.** See above. The two
maximality results are narrower than the interpretation placed on them; the
explicit counterexample is now a required check and a test.

**C5 - the momentum audit was one third of what was frozen.** It computed
three of six Killing charges and none of the four gradient conformal charges.
Constructed data makes the omission concrete: the three originally reported
charges come to `7.1e-16` while the complete six reach `0.180` on the same
data. Both sectors are now computed, and their structure is reported above.

**C6 - the archived verification under-delivered against the freeze.** The
reduction check sampled seven hand-picked degree sets with field data only,
where frozen check 2 demands every pair through degree 6 at nonzero momentum
as well as `p = 0`; and the bilinearity scan used 8 samples where check 4
demands 200. Both now match the freeze: 21 unordered pairs times three data
kinds is 63 reduction cases, and the bilinearity scan runs 200 samples per
pair. The high-sample scan uses the coarse sphere rule, licensed by a new
check that the rules are exact for these polynomial integrands (fine and
coarse grids agree to better than `1e-12` relative).