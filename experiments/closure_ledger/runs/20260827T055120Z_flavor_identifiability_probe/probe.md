# Is the v4 CKM realization a prediction, or a locally flexible fit?

**7/7 checks pass.**

The decisive object is not another pseudoinverse. Build two response maps over one common parameter chart `x` — `J_M = ∂y_M/∂x` for the mass ratios and `J_F = ∂y_F/∂x` for `(θ₁₂, θ₂₃, θ₁₃, δ)` — take the mass-preserving tangent space `N_M = ker J_M`, and form **`K = J_F N_M`**. Its rank counts the physically independent CKM directions reachable *without disturbing the masses*.

## F0 — four coordinates, not ten redundant ones

| observable | rad | deg |
|--|--|--|
| `theta12` | `0.227526` | `13.036°` |
| `theta23` | `0.041932` | `2.403°` |
| `theta13` | `0.003681` | `0.211°` |
| `delta` | `1.150760` | `65.934°` |

Jarlskog `J = 3.093578e-05`, CKM unitary to `6.7e-16`, and the extraction is invariant under random rephasing of both eigenbases to `3.3e-16` — it reads only `|V_ij|` and `J`, never LAPACK's phase convention.

> A unitary 3x3 CKM modulo rephasing has exactly four real parameters, so no observable set can supply more than four independent flavor rows.

*J, alpha, beta, gamma, |V_td|, |V_ts| and the remaining moduli are functions of these four and are NOT independent observable rows.*

## F2 — the calibration-DOF census, measured

The measured calibration dimension of a knob group is `rank J_F` restricted to it — how many independent physical flavor directions it can actually move. **Symbol count and DOF count are different numbers.**

| group | symbols | measured flavor rank |
|--|--|--|
| v3 knobs | 8 | **4** |
| v4 targeted etas | 3 | **3** |
| eta_k3k5_minus retune | 1 | **1** |
| diag_shift_plus | 3 | **2** |
| diag_shift_minus | 3 | **2** |
| all diagonal shifts | 6 | **3** |
| phi_h | 1 | **1** |
| v4 additions, phi_h fixed | 9 | **4** |
| v4 additions, phi_h released | 10 | **4** |
| everything | 18 | **4** |

**The trace direction of each diagonal-shift triple is an exact CKM gauge:**

| triple | `\|J_F·1\|` | `\|J_M·1\|` | realised value | traceless? |
|--|--|--|--|--|
| `diag_shift_plus` | `2.2e-10` | `12.47` | `(0.0018961946, -0.0018961941, -4e-10)` | **yes** |
| `diag_shift_minus` | `8.0e-10` | `12.47` | `(0.1127516149, -0.0647265474, -0.0480250675)` | **yes** |

So a uniform shift within a block moves masses and **cannot have been selected on flavor data**. `diag_shift_plus` and diag_shift_minus carry three symbols each and measure rank 2 each, because the trace direction moves masses but no flavor observable and so cannot have been selected on flavor data. And the realised (+d, -d, 0) sits in a ONE-parameter family inside the two-dimensional traceless subspace.

## F3 — the decisive object: `rank K`

**`rank K = 4`**, against a physical flavor space of dimension 4. `ker J_M` is `14`-dimensional and verified to `5.3e-15`.

Singular values of `K`: `[2.611023, 0.468495, 0.008274, 0.006887]` — spread `379×`, no near-degeneracy.

Rank is invariant under invertible reparameterisation of `x`, which is exactly why it is the right object here: no metric, no pseudoinverse, no chart is being smuggled in. And it is stable across the whole grid:

| step | rcond | dim ker | rank K |
|--|--|--|--|
| `1e-05` | `1e-06` | 14 | **4** |
| `1e-05` | `1e-08` | 14 | **4** |
| `1e-05` | `1e-10` | 14 | **4** |
| `1e-06` | `1e-06` | 14 | **4** |
| `1e-06` | `1e-08` | 14 | **4** |
| `1e-06` | `1e-10` | 14 | **4** |
| `1e-07` | `1e-06` | 14 | **4** |
| `1e-07` | `1e-08` | 14 | **4** |
| `1e-07` | `1e-10` | 14 | **4** |

**Direct construction** — hit an arbitrary target `δy_F` using only mass-preserving directions:

| trial | flavor miss | mass disturbance | `\|δx\|` |
|--|--|--|--|
| 0 | `3.6e-14` | `4.1e-14` | `67.0` |
| 1 | `1.2e-14` | `5.7e-15` | `79.2` |
| 2 | `6.3e-15` | `2.0e-14` | `52.0` |

> **The mass-preserving parameter freedom spans the entire physical flavor space. Fitting the CKM is a successful realisation, not a prediction.**

Rank K = 4 leaves no w with w^T K = 0, so the model predicts no first-order flavor relation w^T dy_F = 0. The round was built to extract one if it existed.

## F4 — the `φ_h` A/B test

The library treats `φ_h = π/k₅` as *derived* and as the source of CP structure. If holding it fixed dropped the flavor rank from 4 to 3, and the missing direction were the observed CP relation, that would be evidence that topology removes one calibration freedom.

| case | coordinates | dim ker `J_M` | **rank K** | singular values |
|--|--|--|--|--|
| `phi_h_released` | 18 | 14 | **4** | `[2.611023, 0.468495, 0.008274, 0.006887]` |
| `phi_h_fixed` | 17 | 13 | **4** | `[0.542415, 0.365819, 0.007932, 0.00357]` |

> **Fixing the derived phase does NOT lower the flavor rank: with phi_h held at pi/k_5 the other fitted matrix elements still span all four physical flavor directions at fixed masses. The derived phase is the single most EFFICIENT CP handle — releasing it multiplies the leading singular value by 4.8 — but efficiency is not identifiability, and it produces no flavor prediction by itself.**

*PR #173 found that ADDING phi_h as an input left the observable rank unchanged. This is the stronger statement: with phi_h held FIXED, the mass-preserving freedom still covers the whole physical flavor space.*

## F5 — the counting claim

| | claimed | measured |
|--|--|--|
| new parameters | 3 | 9 symbols, **4** independent flavor directions |
| new independent observables | 5 | **≤ 4** (the ceiling) |
| net predictive surplus | +2 | **0** |

> **'+5 independent observables' exceeds the ceiling: a unitary 3x3 CKM has exactly four physical parameters, so at most four of the nine quoted flavor-CP observables are independent. Against that ceiling the v4 additions measure 4 independent calibration directions with phi_h fixed, so the net predictive surplus is at most zero, not +2.**

## F6 — the ledger

| claim | verdict | evidence |
|--|--|--|
| the v4 CKM realization is a prediction of the Hamiltonian | **WITHDRAWN** | rank K = 4 = the full physical flavor dimension; the mass-preserving freedom generates arbitrary CKM variation at fixed masses |
| there is a first-order flavor relation to compare with data | **NONE EXISTS** | rank K = 4 leaves 0 left-null vectors, so no w^T dy_F = 0 invariant is predicted |
| the derived phi_h = pi/k_5 produces a flavor prediction | **NOT BY ITSELF** | holding it fixed leaves rank K = 4; it is the most efficient CP handle (leading sv x4.8) but not an identifying one |
| +3 parameters bought +5 independent observables, net +2 | **REFUTED** | a unitary 3x3 CKM has exactly 4 physical parameters, so '+5 independent' exceeds the ceiling; the measured calibration dimension of the v4 additions is 4 with phi_h fixed, giving net surplus <= 0 |
| the six diagonal-shift numbers are six fit freedoms | **MISCOUNTED** | each triple measures flavor rank 2, not 3: the trace direction is an exact CKM gauge that moves masses only, and both realised triples are traceless to ~1e-10 |
| the masses and the CKM constrain the parameters jointly | **NOT AT FIRST ORDER** | ker J_M is 14-dimensional and J_F restricted to it still reaches rank 4; the two sectors do not intersect in a constraining way here |

**Headline.** Rank k = 4: the mass-preserving parameter freedom spans the whole physical flavor space, so the v4 ckm agreement is a successful but locally flexible realisation, not a prediction — and holding the derived phi_h fixed does not change that.

**What this does not say.** The v4 numbers are not wrong and the <= 1% agreement across nine observables is real. What the rank shows is that the agreement is not evidence FOR the Hamiltonian, because the same Hamiltonian could have reproduced any other CKM equally well at the same masses

**Scope.** A local, first-order statement at the v4 lock. Rank says nothing about how far the mass-preserving surface extends before nonlinearity or positivity bites, and the parameter excursions required are not small

**The recommendation.** Per the round's own terms: rank K = 4, so downgrade the CKM realization to a successful but locally flexible realisation, stop the quark parameter archaeology, and return to the trunk
