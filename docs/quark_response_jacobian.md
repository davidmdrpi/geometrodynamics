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
| `phase` | **antiunitary `Z₂` of the spectrum, quadratic** | `H(−φ,p) = H(φ,p)*` for *arbitrary* `p` | `0.0` |
| `partition_mixing` | **unitary-conjugation `Z₂`, quadratic** | `H(−p) = D H(p) D†` | `0.0` |

The two `Z₂`s are of **different kinds** — one antiunitary, one unitary — and
that distinction is what the first draft of this round got wrong.

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

### `phase` — an antiunitary `Z₂`, model-wide

It enters the same-partition off-diagonal through `cos(phase·dk)`, so at the v3
lock `H(−φ) = H(φ)` exactly. The first draft stopped there and concluded that
because the *different*-partition element carries `e^{iφk}` — which is not even
— the symmetry was **a property of the lock**. Turn mixing on and the matrix
equality does break: `0.096` at `φ = 0.37, p = 0.05`.

**That was wrong.** The Hamiltonian satisfies

```
H(−φ, p) = H(φ, p)*        for arbitrary p        exact, to 0.0
```

— same-partition entries are real `cos(φ·dk)`, different-partition entries
`e^{iφk}` conjugate — and since `H` is Hermitian, `H*` = `Hᵀ` is isospectral.
The adiabatic labeller uses overlap magnitudes, which conjugation preserves.
So the **spectrum** is even in `φ` for every `p`: verified `0.0` on the
eigenvalues *and* on the anchored masses at `p = 0.05, 0.3`.

> **Matrix inequality is not spectral asymmetry.** The `Z₂` is a property of
> the model, not of the lock, and it is *antiunitary* — a complex-conjugation
> symmetry, distinct in kind from `partition_mixing`'s unitary one.

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

## Result 1 — the identifiable tangent space is four-dimensional

```
rank J = 4        against 8 first-order knobs
```

The rank cannot exceed the number of independently scored masses — `u` is zero
by construction under `min_eigenvalue` and `d` is the anchor, so the ladder
scores exactly four numbers. **Four of the eight first-order directions
therefore move no observable at all.** The *identifiable tangent space* is
4-dimensional; the model's local non-identifiability is larger still, since the
exact gauge and the two quadratic coordinates are excluded from this count.

| `a` | `σ_a` | `σ_a/σ_1` |
|--|--|--|
| 1 | 19.2190 | 1.0000 |
| 2 | 7.5739 | 0.3941 |
| 3 | 2.5349 | 0.1319 |
| 4 | 0.8519 | 0.0443 |

The four **nonzero** singular values span only `22.6×`, so there is no long
hierarchy among the identifiable directions. That is *not* the same as "not a
sloppy model": `JᵀJ` carries **four exact zero eigenvalues** in these
coordinates, and the full eleven-coordinate model carries more. The dominant
pathology is **structural non-identifiability**, not ill-conditioning.

> `pinhole`, `transport`, `resistance` and `N` are **not independently
> constrained physical observables.** The masses fix at most four combinations
> of everything. Every compensation seen since PR #76 is this one fact.

**Scope.** Every right-space quantity in this document — `δx_min`, the share
decomposition, the singular vectors, the parameter-space cosine, the "invisible
fraction" — is Euclidean in the eight positive knobs in **unit log
coordinates**, with the quadratic coordinates held fixed. Rescaling any column
changes them. The chart-independent statements are the rank, the zero count of
`JᵀJ`, and anything computed in observable space.

**And which knobs drift is itself structure.** The share of each unit
coordinate lying in the identifiable row space runs from `1.0000`
(`uplift_asymmetry`, fully identifiable) through `0.7616` (`gamma_q`), `0.7027`
(`beta`), `0.6970` (`chi_q_k3`), `0.6766` (`pinhole`), `0.1373` (`transport`),
`0.0244` (`resistance`) to `0.0003` (`eta_k3k5_minus`, essentially entirely
null). **"Any knob would drift" is false.**

The stiffest direction `v1` is `uplift_asymmetry` (−0.988) and `u1` is almost
purely `b` (0.9992): the single largest entry in the matrix is
`∂ln m_b/∂ln uplift_asymmetry = −18.96`. **The shell-index axiom
`ε = 1 − 1/k₅² = 24/25` is by far the most load-bearing number in the quark
ladder** — and it is the one PR #272 classified as *exactly invariant*, which
it is, since it reads no radial operator. Being safe from the operator
correction is not the same as being unimportant. Its elasticity had never been
measured.

---

## Result 2 — which way does the correction move?

Two different questions live here, and the first draft of this round ran them
together.

**(a) Which candidate displacement from the lock lands nearer the data?**

| operator | `\|δx_geom\|` | `cos(J·g, r)` | L2 log-residual | max rel err | null-space fraction |
|--|--|--|--|--|--|
| legacy | 0.014819 | **+0.4638** | 0.0548 | 3.19% | 0.6630 |
| corrected | 0.007906 | **−0.6163** | 0.0433 | 3.80% | 0.6922 |

**(b) Which way does the *correction itself* move?** That is
`Δg = g_corrected − g_legacy`, compared against the residual the legacy triple
leaves, `r − J·g_legacy`:

| `\|Δg\|` | `\|r − J·g_legacy\|` | `cos(J·Δg, r − J·g_legacy)` |
|--|--|--|
| 0.015537 | 0.054813 | **+0.8734** (29.1°) |

> **Withdrawn.** The first draft's headline — *"the corrected geometry moves
> away from what the masses want"*, from `cos(J·g_corrected, r) = −0.616` — is
> a true statement about `g_corrected` **as a displacement from the lock**, and
> an invalid basis for a claim about what correcting the operator does. Taken
> as the move it actually makes, the correction points *toward* the residual
> the legacy triple leaves.

**And the two objectives disagree**, which is the substance rather than a
technicality:

| objective | legacy | corrected | |
|--|--|--|--|
| L2 log-residual `\|r − Jg\|` | 0.0548 | 0.0433 | **improves** |
| max relative error | 3.19% | 3.80% | **worsens** |

Both are true. The direction of "improvement" is objective-dependent and **no
metric-free claim is available**. The repository's historical score is the
max-relative-error one, which is why the residual round saw a worsening; the
L2 log-residual, which weights all four species, sees the opposite.

**What survives, and needs no cosine.** Per-knob proximity to a fitted value
carries no information about the sign or size of the joint effect on the
spectrum. That was the point of the round, and it does not depend on which
objective or which displacement one picks.

About two thirds of each geometric displacement lies in the null space of `J`
and moves no observable at all.

## Result 3 — the `0.018%` repair is local flexibility

The minimum-norm correction `δx_min = J⁺r` has `|δ ln p| = 0.0229`, its largest
single knob (`beta`) moving `1.78%`, and drives the ladder from the lock's
`1.61%` to **`0.0179%`** — computed by applying `exp(δ)` to the lock and
re-solving the spectrum, not quoted from prose.

That is important, but not as a success of the physics. It proves something
more basic:

> **The frozen quark Hamiltonian has enough local freedom to absorb essentially
> the entire observed residual.** The old `1.61%` floor was not imposed by the
> functional form of the spectrum. It was a consequence of where the chosen
> lock sits in parameter space.

One measurement below is metric-dependent and one is not.

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

*This decomposition is metric-dependent (see the scope note above). What is
not: the rank deficiency guarantees that some small displacement removing the
residual exists, whichever chart it is measured in.*

### The holdout measures the regulariser, not the Hamiltonian

The first draft fit the correction on three masses, evaluated the fourth, and
read the miss as the Hamiltonian failing a prediction. **That was an
overclaim.**

For each holdout `J_keep` is `3×8` with a five-dimensional kernel, and because
the full `J` has rank 4 the held-out row is *not* in the span of the other
three. So a `z ∈ ker(J_keep)` that moves the withheld species always exists,
and `δ + λz` fits it while holding the other three exact:

| held out | `dim ker` | `\|j_h·Z\|` | min-norm miss | after a kernel shift | norm cost |
|--|--|--|--|--|--|
| `s` | 5 | 4.2520 | +1.912% | **+6.9e-16%** | 1.02× |
| `c` | 5 | 1.1095 | +2.519% | **−3.5e-16%** | 4.77× |
| `b` | 5 | 17.6989 | −9.705% | **−3.6e-15%** | 1.03× |
| `t` | 5 | 1.2102 | −2.521% | **+0.0e+00%** | 2.55× |

The pseudoinverse's miss is therefore a property of the **minimum-log-norm
regulariser**, not of the model. Nothing in the physics selects that
regulariser.

The rank deficiency already established the local flexibility. The holdout
added an overclaim on top of it and is now reported as a regulariser
diagnostic.

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
| four scored masses can identify the current quark parameterisation | **WITHDRAWN** | `rank J = 4` against 8 first-order knobs; `JᵀJ` has 4 exact zeros, the full 11-knob model more |
| `action_base` is a free parameter of the mass spectrum | **WITHDRAWN — EXACT GAUGE** | `H(a) = H(0) + a·I`, removed upstream of the anchor; flat to all orders |
| `phase` and `partition_mixing` are null directions | **MISCLASSIFIED** | two different `Z₂`s — antiunitary and unitary — both quadratic |
| the `phase` `Z₂` is a property of the lock, not the model | **WITHDRAWN — THIS ROUND'S OWN ERROR** | `H(−φ,p) = H(φ,p)*` for arbitrary `p`; the spectrum is even at every mixing |
| per-knob proximity determines the effect on the spectrum | **REFUTED** | the three radial residuals do not compose linearly into the ladder |
| the corrected geometry moves AWAY from what the masses want | **WITHDRAWN — THIS ROUND'S OWN ERROR** | that used `g_corrected` as a displacement from the lock; the move is `Δg`, and `cos = +0.873` toward the residual |
| one objective settles whether the correction helped | **NOT AVAILABLE** | L2 improves `0.0548 → 0.0433` while max-rel-err worsens `3.19% → 3.80%` |
| the v3 ladder's `1.61%` floor is set by the functional form | **WITHDRAWN** | a displacement whose largest knob moves `1.78%` reaches `0.0179%` |
| that repair is evidence for the model | **NOT ESTABLISHED** | 98.3% rides the weakest direction in the chosen chart, and rank deficiency guarantees such a displacement exists |
| leave-one-species-out shows the Hamiltonian fails to predict | **WITHDRAWN — THIS ROUND'S OWN OVERCLAIM** | `ker(J_keep)` is 5-D and the held-out row is reachable; a kernel shift fits it to `~1e-15` at `1.02–4.77×` the norm |
| the quark model is "sloppy" in the Sethna sense | **NOT THE RIGHT DIAGNOSIS** | nonzero singular values span only `22.6×`, but the rank deficiency is exact — the pathology is rank, not conditioning |
| `N = 466` drifting is a defect of `N`, and any knob would drift | **REFRAMED, SECOND HALF FALSE** | identifiable share runs `1.0000` (`uplift_asymmetry`) to `0.0003` (`eta_k3k5_minus`) |
| the missing correlation is a scalar relation `R = f(p, T)` | **REFUTED** | the nearest pair is `gamma_q`/`transport` at 178.9°, not `transport`/`resistance` |

---

## What would settle it

**More observables, not more knobs.** The rank is capped by the observable
count, so no amount of further geometric derivation can make `pinhole`,
`transport` and `resistance` separately identifiable from four masses.

The v4 flavor-CP layer produces CKM angles and the Jarlskog invariant from the
*same* Hamiltonian and inherits the v3 eigenvalues, so it **adds observable
rows**. But it is **not automatically column-free**: v4 `QuarkParams` carries
`phi_h`, targeted `eta` couplings on three elements, and per-shell diagonal
shifts. Whether those count as new response columns depends on which are
externally derived and fixed — `phi_h = π/k₅` is derived (PR #159), the `eta`s
and diagonal shifts are fitted. **A joint identifiability audit has to make
that call explicitly before claiming the rank rises.**

Two smaller open items, both from the SVD rather than added by hand: the
near-degenerate pair `gamma_q`/`transport` at `178.9°` is essentially one
handle with two names — an overall scale of the `d`-anchored ladder — and
deserves an algebraic explanation of the kind the three nulls now have; and
`chi_q_k3`/`eta_k3k5_minus` at `13.6°` are nearly parallel, with the latter
almost entirely in the null space, so it is not clear both should exist.
