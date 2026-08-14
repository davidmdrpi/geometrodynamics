# Does a solved field reproduce the branch ledger, with its phases?

`geometrodynamics/waves/field_solve.py` · renderer
`scripts/geometrodynamics_v54_field_solve.py` · probe
`experiments/closure_ledger/field_solve_probe.py`

## 0. The question

PR #253 built a **ray** ledger — short way, long way, winding — and finished by
conceding that rank counting had reached its limit: it could not supply a
quantity that *vanishes* when a source is removed rather than merely becoming
underdetermined. Before building that quantity, one thing has to be checked.

> **Does a solved wave field, with the throat imposed as the same
> identification, reproduce the exact branch ledger — its arrival times, its
> amplitudes, and phases the ray picture never had?**

**Yes**, and more sharply than a stationary-phase argument would give.

## 1. Scope, stated before the result

* a **linear scalar field** on a **fixed** background (the Einstein static
  universe `S³ × R`), `c = 1`;
* the throat is still an **identification map**, not a solution.
  `shells/junction.py` (PR #249) priced it and the bill is inherited, unpaid;
* **no backreaction, no stress tensor, no topology change, no rate**;
* **no two-source invariant yet** — the cross-quantity that vanishes when a
  source is removed is the next step, not this one;
* the self-consistent sum over **repeated** throat traversals is PR #251's fixed
  point and is **not** redone here. Only single-traversal arrivals are computed.

## 2. Why the answer is exact rather than asymptotic

On `S³` the scalar Laplacian has eigenvalues `−n(n+2)`, and `R = 6`, so the
conformal term `ξR = 1` gives

```
ω² = n(n+2) + 1 = (n+1)²        ⟹        ω_n = n + 1
```

**Integer frequencies.** The retarded Green function is therefore exactly
periodic, and Poisson summation turns the mode sum into a sum over images:

```
G(χ,t) = 1/(4π sin χ) [ Σ_k δ(t − χ − 2πk) − Σ_k δ(t + χ − 2πk) ]
```

So the geometric-optics branches are the **exact support** of the field. Nothing
has to be argued asymptotically, and there is no small parameter.

**The solve against the closed form.** A truncated eigenmode sum — which never
sees an image — against the image sum — which never sees a mode: **`8.3e-13`**.

## 3. The field's support is the ray ledger

Peaks read off the **solved** field land on PR #253's branch times:

| | |
| --- | ---: |
| worst time error | **`3.0e-04`** |
| grid spacing | `6.3e-04` |
| branches matched by a peak | **all** |
| spurious peaks | **none** |

The error is half a grid cell, so it is grid-limited rather than a
disagreement. `t = χ + 2πk` is the short way with `k` windings; `t = 2πk − χ` is
the long way — precisely the branch set that round enumerated by hand.

## 4. The amplitude is PR #251's shell law, derived rather than imposed

That round set `A ∝ 1/sin χ` by conserving energy across a shell of area
`4π sin²χ`. Here nothing is imposed: the peak of the solved field, multiplied by
`sin χ`, is the same constant at every `χ`.

| | |
| --- | ---: |
| relative spread of `peak × sin χ` | **`7.0e-16`** |
| constant | `0.634936` |
| predicted `1/(4π w √(2π))` | `0.634936` |

## 5. And the field supplies phases the ray ledger could not

**This is the part that is new information rather than a re-derivation.**

Every arrival carries a sign, and it is `(−1)^m` with `m` the number of **focal
crossings** — the antipode at `t = π`, the source point again at `t = 2π`, and so
on. That is the **Maslov index**.

| `t` | branch | crossings | field sign | Maslov |
| ---: | --- | ---: | :---: | :---: |
| `0.700` | short, `k=0` | 0 | `+` | `+` |
| `5.583` | long, `k=0` | 1 | `−` | `−` |
| `6.983` | short, `k=1` | 2 | `+` | `+` |
| `11.866` | long, `k=1` | 3 | `−` | `−` |

**12 of 12** signs agree across the `χ` values checked.

A path-length ledger gives arrival times and has **no way** to produce a sign.
This is the first quantity in the arc that the ray picture could not in
principle have carried — which is exactly the kind of thing the move to a field
was supposed to buy.

## 6. And the ledger belongs to the *conformal* field specifically

Conformal coupling is not a convenience. The **minimally** coupled massless
field has `ω = √(n(n+2))`, irrational, so its Green function has no image
structure at all:

| coupling | `ω` | worst `\|φ\|` between branches, relative to peak |
| --- | --- | ---: |
| conformal | `n+1` | **`4.0e-08`** |
| minimal | `√(n(n+2))` | **`0.63`** |

Sixty-three percent of the peak amplitude sits *between* the arrivals. Huygens'
principle holds for one and fails for the other.

PR #253 never said which field its ledger described, **because rays cannot tell
the two apart** — both have the same geodesics. That is a scope fact the ray
picture was structurally unable to notice, and it only appears once the field is
solved.

## 7. The throat reproduces the closure condition

Imposed as the same identification, now on a field:

```
φ(M⁺, t) = η · φ(M⁻, t + Δ)
```

A through-throat contribution therefore lands at `ℓ₁ + Δ + ℓ₂`, summed over
branch pairs. Setting `Δ` to minus a branch-pair sum — **exactly PR #253's
closure condition** — puts an arrival back on the emission event, at `t = 0`.
**9 of 9** branch-pair choices close.

And the field adds what that round could not: the closing arrival's **sign**, `η`
times the two Maslov factors. Flipping the orientation `η` flips it.

## 8. What this closes, and what it does not

**Closes:** that the ray ledger of PR #253 and the shell law of PR #251 are both
consequences of a solved linear field rather than modelling choices; that the
branches are exact support rather than a stationary-phase approximation; that
the ledger carries a Maslov phase the ray picture omitted; and that the sharp
ledger belongs to the conformally coupled field specifically.

**Does not close:** everything the roadmap has after this. No two-source
invariant — the quantity that *vanishes* when a source is removed is still not
built. No stationary-action test, so "past and future jointly satisfy least
action" remains interpretation. No backreaction, no topology change, and the
throat is still an identification map with an unpaid bill.

## 9. What this changes for the next step

The two-source invariant `𝒞 = A_A² A_B² (k_A · k_B)²` now has a concrete
obstacle, and it is the branch structure established here rather than a
technicality. At a generic point the field arrives on **several branches at
once**, each with its own propagation direction and its own Maslov sign. So
`k_A = ∇S_A` is genuinely multivalued, and a `𝒞` written as though it were
single-valued would silently average over sheets — the field-level version of
the scope error PR #253 was corrected for.

The branch labels and signs built here are what a branch-resolved `𝒞` will be
written in terms of.
