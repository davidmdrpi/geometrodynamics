# Is the v4 CKM realization a prediction, or a locally flexible fit?

*Module `geometrodynamics/qcd/flavor_identifiability.py`, probe
`experiments/closure_ledger/flavor_identifiability_probe.py` (9/9), tests
`tests/test_flavor_identifiability.py`.*

---

## The answer

**`rank K = 4`. The CKM realization is a fit, not a prediction.**

---

## The construction

Not another pseudoinverse. Build two response maps over **one common parameter
chart** `x`:

```
J_M = ∂y_M/∂x ,   y_M = ln(m_{s,c,b,t} / m_d)
J_F = ∂y_F/∂x ,   y_F = (θ₁₂, θ₂₃, θ₁₃, δ)
```

take the mass-preserving tangent space `N_M = ker J_M`, and form

```
K = J_F N_M
```

`rank K` counts the physically independent CKM directions reachable **without
disturbing the masses**. Rank is invariant under any invertible
reparameterisation of `x`, which is exactly why it is the right object here —
no metric, no chart, no regulariser is being smuggled in. That was the flaw in
the previous round's pseudoinverse geometry, and it is structurally absent
from this one.

### Four observables, not ten

`y_F = (θ₁₂, θ₂₃, θ₁₃, δ)` in the PDG convention, built **only** from `|V_ij|`
and the Jarlskog invariant — both invariant under the arbitrary per-column
phases `eigh` returns. Verified to `1e-16` under random rephasing of both
eigenbases. Reading the raw complex entries would have made every Jacobian
column a function of LAPACK's convention.

`J`, `α`, `β`, `γ`, `|V_td|`, `|V_ts|` and the remaining moduli are **validation
outputs, not observable rows**. They are functions of these four.

| observable | value | deg | measured |
|--|--|--|--|
| `θ₁₂` | 0.227526 | 13.036° | ~0.2265 |
| `θ₂₃` | 0.041932 | 2.403° | ~0.0424 |
| `θ₁₃` | 0.003681 | 0.211° | ~0.00365 |
| `δ` | 1.150760 | 65.934° | ~1.144 |
| `J` | 3.0936e-5 | — | ~3.08e-5 |

The chart `x` carries 18 coordinates: the eight v3 knobs, the three v4 targeted
`η`s, `φ_h`, and the six diagonal-shift components.

---

## The result

```
rank J_M = 4        dim ker J_M = 14        rank K = 4
singular values of K = [2.611, 0.468, 0.00827, 0.00689]
```

Spread `379×`, no near-degeneracy, kernel verified to `5.3e-15`. Stable across
a 3×3 grid of finite-difference step (`1e-5…1e-7`) and null-space cutoff
(`1e-6…1e-10`) — rank 4 in all nine cells.

**A tangent-space construction** — solving for an arbitrary target `δy_F` using
only mass-preserving directions. This is algebra on the matrices above, *not*
independent evidence; the nonlinear check is below:

| trial | flavor miss | mass disturbance | `\|δx\|` |
|--|--|--|--|
| 0 | 3.6e-14 | 4.1e-14 | 67.0 |
| 1 | 1.2e-14 | 5.7e-15 | 79.2 |
| 2 | 6.3e-15 | 2.0e-14 | 52.0 |

> **The mass-preserving parameter freedom spans the entire physical flavor
> space.** Fitting the CKM is a successful realisation, not a prediction.

**And therefore there is no prediction to extract.** `rank K = 4` leaves zero
left-null vectors, so there is no `w` with `wᵀK = 0` and no first-order
invariant `wᵀδy_F = 0` to compare against experiment. The round was built to
extract one if it existed.

---

## The restricted question — which is the decisive one

*Added after review; the first draft did not compute it.*

The headline `rank K` lets **every** coordinate move, including the eight v3
knobs. But the v4 construction was explicitly **additive over a frozen v3
lock**, so the question that bears on predictive surplus is the restricted one:

```
K_v4 = J_F[:,G] · ker(J_M[:,G])        with the v3 knobs held fixed
```

The first draft reported `rank J_F[:,G]` — which lets the masses move — beside
an unrestricted `K` over all 18 coordinates. **Neither is the restricted
mass-preserving map**, and neither establishes that the *actual* v4 calibration
freedoms span all four CKM directions while preserving the frozen masses.

| group | coordinates | `dim ker J_M` | **rank K_v4** | smallest σ |
|--|--|--|--|--|
| v4 additions only, `φ_h` fixed | 9 | 5 | **4** | 8.03e-06 |
| v4 additions only, `φ_h` released | 10 | 6 | **4** | 1.21e-05 |
| **with the `η₃₅` retune, `φ_h` fixed** | **10** | 6 | **4** | **3.10e-03** |
| with the `η₃₅` retune, `φ_h` released | 11 | 7 | **4** | 6.72e-03 |

> **They do span it.** `rank K_v4 = 4` with the frozen v3 lock *and* `φ_h` both
> held fixed — so the headline is **stronger** than the first draft claimed,
> not weaker. The surplus claim fails on the restricted question, not merely on
> the permissive one.

**The group includes the `eta_k3k5_minus` retune** (`5.0 → 5.586`). The round's
own rule was to count every numerical quantity selected using flavor data
regardless of whether the symbol already existed, and the first draft's
nine-symbol group broke that rule by omitting it. Including it also **improves
the conditioning by two and a half orders** — the smallest singular value rises
`8.0e-06 → 3.1e-03` — so counting it honestly makes the rank-4 claim more
robust, not less.

---

## And it is nonlinearly true

*Added after review; the first draft's evidence here was circular.*

The tangent-space construction solves `K c = δy_F` and reports `J_M N c ≈ 0`.
That is **algebra on the same first-order matrices whose rank was just
computed** — not independent evidence, and the first draft presented it as
though it were.

The real test: scale one of those directions by `ε`, re-run the **actual
Hamiltonian**, and check that both errors fall as `O(ε²)`.

| `ε` | mass error | flavor miss |
|--|--|--|
| 1.00e-04 | 1.200e-06 | 2.052e-06 |
| 5.00e-05 | 3.000e-07 | 5.119e-07 |
| 2.50e-05 | 7.501e-08 | 1.278e-07 |
| 1.25e-05 | 1.876e-08 | 3.194e-08 |

Halving `ε` quarters both — mass ratios `[3.999, 4.000, 3.999]`, miss ratios
`[4.009, 4.005, 4.002]`. Clean second order, so the reachability is a property
of the model and not only of its Jacobian.

**What it licenses.** The claim is about **infinitesimal** CKM directions
reachable locally, not about arbitrary finite CKM matrices. The excursion
needed for a unit `δy_F` is `|δx| ≈ 50–80`, far outside the regime this scaling
tests. Wherever the first draft said *"the same Hamiltonian could have
reproduced any other CKM"*, read *"any infinitesimal CKM direction, locally"*.

---

## The `φ_h` A/B test

The library treats `φ_h = π/k₅` as *derived* and as the source of CP structure.
If holding it fixed dropped the flavor rank from 4 to 3, and the missing
direction were the observed CP relation, that would be meaningful evidence that
topology removes one calibration freedom.

| case | coordinates | dim ker `J_M` | **rank K** | singular values |
|--|--|--|--|--|
| `φ_h` released | 18 | 14 | **4** | [2.611, 0.468, 0.00827, 0.00689] |
| `φ_h` fixed | 17 | 13 | **4** | [0.542, 0.366, 0.00793, 0.00357] |

**Fixing it does not lower the rank.** The other fitted matrix elements absorb
arbitrary CKM data on their own, so the derived phase produces no flavor
prediction by itself.

**Is it actually a CP handle?** Yes, and now shown rather than asserted: its
`J_F` column is **0.99978** along `δ` (`θ₁₂ −0.0208`, `θ₂₃ +0.0001`,
`θ₁₃ −0.0003`, `δ −0.9998`). That share **is** chart-independent — rescaling
the `φ_h` coordinate scales its column without rotating it.

It is also the most *efficient* CP handle: releasing it multiplies the leading
singular value by **4.8**. But **that ratio is chart-dependent**, unlike the
rank and unlike the direction: singular-value magnitude is Euclidean in the
chosen linear absolute coordinates, and rescaling any column changes it. Kept
as a scoped numerical diagnostic, not a physical efficiency statement.

Efficiency and identifiability are different, and only the second is invariant.

This confirms and sharpens **PR #173**, which found that *adding* `φ_h` as an
input left the observable rank unchanged and concluded that *"CP at zero
parameters" is a counting economy, not a Jacobian reduction*. That was a
statement about spanning. This is the stronger one: with `φ_h` held at its
derived value, the mass-preserving freedom still covers all four physical
flavor coordinates.

---

## The calibration-DOF census, measured rather than counted

"New symbol count" and "calibration degree-of-freedom count" are different
numbers. Measuring the second as `rank J_F` restricted to each knob group:

| group | symbols | measured flavor rank |
|--|--|--|
| v3 knobs | 8 | 4 |
| v4 targeted `η`s | 3 | **3** |
| `eta_k3k5_minus` retune | 1 | **1** |
| `diag_shift_plus` | 3 | **2** |
| `diag_shift_minus` | 3 | **2** |
| all diagonal shifts | 6 | **3** |
| `φ_h` | 1 | **1** |
| v4 additions, `φ_h` fixed | 9 | **4** |
| v4 additions, `φ_h` released | 10 | **4** |

Two structural facts fall out, and both were invisible to a symbol count.

**The trace direction of each diagonal-shift triple is an exact CKM gauge.**

| triple | `\|J_F·1\|` | `\|J_M·1\|` | realised value | traceless? |
|--|--|--|--|--|
| `diag_shift_plus` | 2.2e-10 | 12.47 | `(+d, −d, 0)` | **yes** |
| `diag_shift_minus` | 8.0e-10 | 12.47 | `(0.1128, −0.0647, −0.0480)` | **yes** |

A uniform shift within a block moves the masses and moves **no** flavor
observable, so it cannot have been selected using flavor data. That is why each
triple measures rank 2 rather than 3. And both realised triples are traceless
to `~1e-10`, with `diag_shift_plus` collapsing further to a one-parameter
`(+d, −d, 0)` family.

**The v3 knobs alone already measure flavor rank 4** at the v4 lock — though
the mass-preserving version of that statement is marginal (`cond ≈ 3.2e5`, with
the fourth singular value at `8.7e-7`), so it is recorded as a curiosity rather
than as a claim. The robust statement is the full-chart one above.

---

## The "+3 parameters for +5 independent observables" claim

It cannot hold, for a reason that needs no computation.

> **A unitary 3×3 CKM has exactly four physical parameters.**

The nine quoted flavor-CP observables — `|V_us|`, `|V_cb|`, `|V_ub|`, `|V_td|`,
`|V_ts|`, `J`, `β`, `γ`, `sin δ` — are all functions of those four, so at most
four of them are independent. "+5 independent observables" exceeds the ceiling.

| | claimed | measured |
|--|--|--|
| new parameters | 3 | 9 symbols, **4** independent flavor directions |
| new independent observables | 5 | **≤ 4** (the ceiling) |
| net predictive surplus | +2 | **≤ 0** |

---

## What this does not say

The v4 numbers are not wrong. The lock reproduces nine observables at `≤ 1%`
and that is a real property of the realisation, unaffected by anything here.

What the rank shows is that the agreement is **not evidence for the
Hamiltonian**, because the same Hamiltonian could have reproduced any *locally
neighbouring* CKM equally well at the same masses. A fit that could have
accommodated any nearby outcome does not gain credit for accommodating the
observed one.

**Scope.** This is a **local, first-order** statement at the v4 lock. A rank
says nothing about how far the mass-preserving surface extends before
nonlinearity or positivity bites. And the excursions are not small — an
arbitrary unit `δy_F` needs `|δx| ≈ 50–80` in these coordinates — so
"reachable" is a statement about *directions*, not about comfortable distances.

---

## The ledger

| claim | verdict | evidence |
|--|--|--|
| the v4 CKM realization is a prediction of the Hamiltonian | **WITHDRAWN** | `rank K = 4` over all coordinates, **and `rank K_v4 = 4`** over the v4 calibration group alone with the frozen v3 lock and `φ_h` both fixed |
| the first draft's headline used the right restricted map | **NO — THIS ROUND'S OWN GAP** | it reported `rank J_F[:,G]` (masses free) beside an unrestricted `K`; the restricted map gives the same rank 4, so the headline strengthens |
| the tangent-space construction is independent evidence | **NO — RELABELLED** | algebra on the same matrices; the independent check is the `O(ε²)` re-run |
| there is a first-order flavor relation to compare with data | **NONE EXISTS** | zero left-null vectors of `K` |
| the derived `φ_h = π/k₅` produces a flavor prediction | **NOT BY ITSELF** | holding it fixed leaves `rank K = 4`; efficient (×4.8) but not identifying |
| +3 parameters bought +5 independent observables, net +2 | **REFUTED** | the ceiling is 4; measured calibration dimension is 4; net ≤ 0 |
| the six diagonal-shift numbers are six fit freedoms | **MISCOUNTED** | each triple measures flavor rank 2; the trace is an exact CKM gauge |
| the masses and the CKM constrain the parameters jointly | **NOT AT FIRST ORDER** | `ker J_M` is 14-dimensional and `J_F` on it still reaches rank 4 |

---

## The recommendation, on this round's own terms

`rank K = 4`, so: **downgrade the CKM realization to a successful but locally
flexible realisation, stop the quark parameter archaeology, and return to the
trunk.**

The three audits (#271, #272, this one) have converged on one structural
statement worth carrying forward: **the quark sector's observables are too few
to identify its parameters, in both the mass and the flavor sectors.** No
further geometric derivation of an individual knob can change that, because the
ceiling is set by the observable count, not by the quality of the derivation.
