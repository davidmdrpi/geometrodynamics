# The local geometry of the quark fit manifold

*Module `geometrodynamics/qcd/response_jacobian.py`, probe
`experiments/closure_ledger/quark_response_jacobian_probe.py` (7/7), tests
`tests/test_quark_response_jacobian.py`.*

---

## Why a matrix and not more scalars

PR #272 measured the elasticity of one knob at a time. It found the three quark
residuals individually right and jointly wrong, and named the gap a *missing
correlation* — guessing it lived between `transport` and `resistance`, since
`transport` appears both alone and inside `resistance`.

**That guess was wrong**, and more importantly the object was wrong. A
collection of scalar residuals cannot express what is happening. The right
object is the response map

```
J_ia = ∂ ln(m_i / m_d) / ∂ ln p_a ,        i ∈ {s, c, b, t}
```

and its singular value decomposition, which settles in one calculation what no
number of per-knob percentages could.

---

## The null structure: three different objects

Three directions show zero first derivative at the lock. They are
**mathematically different** and putting them in one category was the previous
round's second error.

| direction | status | mechanism | verified |
|--|--|--|--|
| `action_base` | **exact invariance, all orders** | `H(a) = H(0) + a·I` | `1.8e-13` |
| `phase` | **Z₂-even, quadratic** | enters only as `cos(φ·dk)` | `0.0` |
| `partition_mixing` | **Z₂-even, quadratic** | `H(−p) = D H(p) D†` | `0.0` |

### `action_base` — an exact gauge

`_diagonal_entry` adds it identically to all six diagonal entries and nowhere
off-diagonal, so `H(a) = H(0) + a·I`, every eigenvalue shifts by `a`, and the
`min_eigenvalue` spectrum-zero subtraction removes it identically. A **3×**
change moves the anchored masses by `1.6e-12` — flat at finite displacement,
not merely stationary.

One correction to the natural expectation: this is **not the `d`-anchor doing
the work.** The cancellation happens in the zero-point subtraction, which is
*upstream* of the anchor, and the alternative `spectrum_zero_mode="action_base"`
kills it too. So `action_base` does **not** reappear in an absolute-scale
observable of this model — it is a gauge of the model as defined, and it is
dropped from the first-order parameter space for that reason.

### `phase` — even, but only here

It enters the same-partition off-diagonal through `cos(phase·dk)`, so
`H(−φ) = H(φ)` exactly. But the *different*-partition element carries
`e^{i·phase·k}`, which is not even; it is switched off only because
`partition_mixing = 0` at the v3 lock. Turn mixing on and the evenness breaks
(`0.096` at `φ = 0.37`).

> **The Z₂ is a property of the lock, not of the model.**

### `partition_mixing` — a discrete gauge symmetry

It appears only on `+`/`−` off-diagonal elements, so flipping its sign is
conjugation by `D = diag(σ(p)) = diag(+1,−1,+1,−1,+1,−1)`:

```
H(−p) = D H(p) D†        exact, to 0.0
```

hence isospectral, hence the eigenvalues are even in `p`. That is stronger than
"even function" — it is an exact discrete symmetry of the spectrum.

### The chart

Both quadratic directions are handled by `q = x²`, under which `∂ln m/∂q` is
finite and constant. The licence is empirical and exact: a 10× step in `x`
gives **100×** the response, to 2%. Their first derivatives in `x` are not
informative and are never reported as such.

---

## Result 1 — the fit manifold is four-dimensional

```
rank J = 4        against 8 first-order knobs
```

The rank cannot exceed the number of independently scored masses — `u` is zero
by construction under `min_eigenvalue` and `d` is the anchor, so the ladder
scores exactly four numbers. **Four of the eight first-order directions
therefore move no observable at all**, before any accidental degeneracy, and
that is before counting the exact gauge and the two quadratic directions.

| `a` | `σ_a` | `σ_a/σ_1` |
|--|--|--|
| 1 | 19.2190 | 1.0000 |
| 2 | 7.5739 | 0.3941 |
| 3 | 2.5349 | 0.1319 |
| 4 | 0.8519 | 0.0443 |

**Condition number 22.6.** This is *not* a sloppy model in the Sethna sense —
sloppy spectra span 10³–10⁶. The identifiable subspace is well conditioned. The
degeneracy is **dimensional**: eleven knobs against four observables.

> `pinhole`, `transport`, `resistance` and `N` are **not independently
> constrained physical observables.** The masses fix at most four combinations
> of everything. Every compensation seen since PR #76 is this one fact.

The stiffest direction `v1` is `uplift_asymmetry` (−0.988) and `u1` is almost
purely `b` (0.9992): the single largest entry in the matrix is
`∂ln m_b/∂ln uplift_asymmetry = −18.96`. **The shell-index axiom
`ε = 1 − 1/k₅² = 24/25` is by far the most load-bearing number in the quark
ladder** — and it is the one PR #272 classified as *exactly invariant*, which
it is, since it reads no radial operator. Being safe from the operator
correction is not the same as being unimportant. Its elasticity had never been
measured.

---

## Result 2 — the correction moves *away* from what the masses want

The decisive test. Project the operator-induced parameter displacement
`δx_geom` onto the mass-optimal direction `δx_min = J⁺r`.

| operator | `\|δx_geom\|` | `cos Θ` (parameter) | `cos Θ` (observable) | invisible fraction |
|--|--|--|--|--|
| legacy | 0.014819 | +0.2870 | **+0.4638** (62.4°) | 0.6630 |
| corrected | 0.007906 | −0.2785 | **−0.6163** (128.0°) | 0.6922 |

This is the third case of the trichotomy. The legacy geometry had partial
overlap with what the data want; **the corrected geometry is actively opposed
to it.**

> PR #272's three per-knob improvements were not merely uninformative about the
> ladder. They were **misleading** about it. Every residual moved toward its
> locked value while the displacement they jointly produce turned from
> partly-helpful to actively-harmful.

A scalar residual cannot see this. Only the projection can — which is the
argument for the matrix in one line.

And about **two thirds** of each geometric displacement lies in the null space
of `J` and moves no observable at all, so even the part that points somewhere
mostly points nowhere.

---

## Result 3 — the `0.018%` repair is local flexibility

The minimum-norm correction `δx_min = J⁺r` has `|δ ln p| = 0.0229` — every knob
moving by ~1% or less — and drives the ladder from the lock's `1.61%` to
**`0.018%`** (exact nonlinear re-run, not the linearisation).

That is important, but not as a success of the physics. It proves something
more basic:

> **The frozen quark Hamiltonian has enough local freedom to absorb essentially
> the entire observed residual.** The old `1.61%` floor was not imposed by the
> functional form of the spectrum. It was a consequence of where the chosen
> lock sits in parameter space.

Two measurements say it is a compensator rather than structure.

### Which directions do the repairing

| `a` | `σ_a` | `c_a = u_aᵀr` | `c_a/σ_a` | share of `\|δx\|²` |
|--|--|--|--|--|
| 1 | 19.2190 | −0.007223 | −0.000376 | 0.0003 |
| 2 | 7.5739 | +0.015230 | +0.002011 | 0.0077 |
| 3 | 2.5349 | −0.005377 | −0.002121 | 0.0086 |
| 4 | 0.8519 | +0.019377 | **+0.022745** | **0.9835** |

Reporting `σ_a` alone would have said nothing. **98.35% of the repair rides the
weakest direction** and the stiffest contributes 0.03% — the model's
least-constrained freedom doing essentially all the work.

### Leave-one-species-out

Fit the correction on three masses; evaluate the fourth with **nothing
readjusted**.

| held out | `\|δx\|` | fitted three (exact) | **held-out (exact)** |
|--|--|--|--|
| `s` | 0.022498 | 0.015% | **+1.905%** |
| `c` | 0.004813 | 0.041% | **+2.562%** |
| `b` | 0.022198 | 0.059% | **−10.410%** |
| `t` | 0.008997 | 0.002% | **−2.519%** |

The fitted species land at `0.002–0.06%` every time; the withheld one does not
come along. **Local flexibility, not structure.**

---

## What this retroactively qualifies

PR #272 reported that the locked ladder "wants" `pinhole = 22.228`. That was a
single-axis scan. With all knobs free, the minimum-norm solution wants
`22.0225` instead — close to the *legacy* derived `22.008`, which is part of
why the legacy set looked better.

> **"What the ladder wants" is not a number.** It depends on what else is
> allowed to move, and it ranges over at least `22.02 … 22.23` depending on
> that choice. Quoting a single value implies an identifiability the model does
> not have.

---

## The ledger

| claim | verdict | evidence |
|--|--|--|
| `action_base` is a free parameter of the mass spectrum | **WITHDRAWN — EXACT GAUGE** | `H(a) = H(0) + a·I`, removed upstream of the anchor; flat to all orders |
| `phase` and `partition_mixing` are null directions | **MISCLASSIFIED** | Z₂-even at the lock, quadratic away from it |
| `pinhole`, `transport`, `resistance` independently constrained by the masses | **WITHDRAWN** | `rank J = 4` against 8 knobs; four directions move no observable |
| the v3 ladder's `1.61%` floor is set by the functional form | **WITHDRAWN** | a `0.0229` displacement reaches `0.018%` |
| that `0.018%` fit is evidence for the model | **NOT ESTABLISHED** | 98.3% rides the weakest direction; leave-one-out misses by 10.4% |
| PR #272's per-knob improvements moved the ladder toward the data | **REFUTED** | `cos Θ = −0.616` corrected against `+0.464` legacy |
| the quark model is "sloppy" in the Sethna sense | **NOT SUPPORTED** | condition number 22.6 — the degeneracy is dimensional |
| `N = 466` drifting is a defect of `N` | **REFRAMED** | symptom of a four-dimensional observable space; any knob would drift |
| the missing correlation is a scalar relation `R = f(p, T)` | **REFUTED** | the nearest pair is `gamma_q`/`transport` at 178.9°, not `transport`/`resistance` |

---

## What would settle it

**More observables, not more knobs.** The rank is capped by the observable
count, so no amount of further geometric derivation can make `pinhole`,
`transport` and `resistance` separately identifiable from four masses.

The v4 flavor-CP layer already produces CKM angles and the Jarlskog invariant
from the *same* Hamiltonian, and by construction inherits the v3 eigenvalues —
so those observables are free of the mass sector's parameters in the sense that
matters here: they add rows to `J` without adding columns. Extending the
Jacobian to the joint mass + CKM response is the natural successor, and it is
the only route that can raise the achievable rank above four and make these
directions identifiable for the first time.

Two smaller open items: the near-degenerate pair `gamma_q`/`transport` at
`178.9°` is essentially one handle with two names — an overall scale of the
`d`-anchored ladder — and deserves an algebraic explanation of the kind the
three nulls now have; and `chi_q_k3`/`eta_k3k5_minus` at `13.6°` are nearly
parallel and probably should not both exist.
