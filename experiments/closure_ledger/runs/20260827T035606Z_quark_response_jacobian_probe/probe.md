# The local geometry of the quark fit manifold

**7/7 checks pass.**

PR #272 measured one knob at a time and named the gap between individually-right and jointly-wrong residuals a *missing correlation*, guessing it lived between `transport` and `resistance`. **That guess was wrong, and the object is not a scalar relation.** It is the response map `J_ia = ∂ln(m_i/m_d)/∂ln p_a` and its SVD.

## R0 — three apparent nulls, three different objects

| direction | status | mechanism | check |
|--|--|--|--|
| `action_base` | **EXACT INVARIANCE (all orders)** | `H(a) = H(0) + a·I`, killed by the spectrum-zero subtraction | `1.8e-13` |
| `phase` | **ANTIUNITARY Z2 OF THE SPECTRUM, QUADRATIC** | `H(−φ, p) = H(φ, p)*` for *arbitrary* `p`; `H` Hermitian ⟹ isospectral | `0.0` |
| `partition_mixing` | **UNITARY-CONJUGATION Z2, QUADRATIC** | `H(−p) = D H(p) D†`, `D = diag(σ)` | `0.0` |

`action_base` is removed **upstream of the `d`-anchor** — the zero-point subtraction does it, and `spectrum_zero_mode="action_base"` kills it too. So it does *not* reappear in an absolute-scale observable of this model, and it is dropped from the first-order parameter space.

**Corrected from the first draft.** That draft called the `phase` `Z₂` *a property of the lock*, because the **matrix** stops satisfying `H(−φ) = H(φ)` once `partition_mixing ≠ 0` (`0.577` at `φ = 0.37, p = 0.3`). But the **spectrum** is even at every mixing — `0.0` on the eigenvalues *and* on the anchored masses. Matrix inequality is not spectral asymmetry.

> phase is an antiunitary (complex-conjugation) symmetry of the spectrum; partition_mixing is a unitary-conjugation symmetry. Both give a vanishing first derivative and a live second one, for different reasons.

Both are handled by the chart `q = x²`, under which `∂ln m/∂q` is finite and constant:

| knob | `∂ln m/∂q` | spread over two decades in `q` |
|--|--|--|
| `phase` | `5.0944` | `2.0e-06` |
| `partition_mixing` | `156.5314` | `8.9e-05` |

*`action_base` is flat to all orders and carries no information at any displacement; phase and partition_mixing are flat only at the symmetric point and carry ordinary quadratic response away from it.*

## R1 — the Jacobian converges

`extract_physical_spectrum` reassigns species by adiabatic continuation, so a relabelling inside the difference stencil would show up as step-dependent noise. None appears — every column norm is stable across three decades of step size.

| knob | column norm | relative spread, `1e-3`…`1e-6` |
|--|--|--|
| `beta` | `1.2486` | `6.5e-08` |
| `gamma_q` | `2.4874` | `4.2e-07` |
| `transport` | `1.0530` | `5.9e-06` |
| `pinhole` | `5.7272` | `1.6e-05` |
| `resistance` | `0.5388` | `2.2e-06` |
| `uplift_asymmetry` | `19.0165` | `1.2e-04` |
| `chi_q_k3` | `5.5046` | `1.5e-05` |
| `eta_k3k5_minus` | `0.1078` | `2.6e-07` |

The largest entry in the whole matrix is `∂ln m_b/∂ln uplift_asymmetry = −18.96`: **the shell-index axiom `ε = 1 − 1/k₅² = 24/25` is by far the most load-bearing number in the quark ladder**, and it is the one PR #272 classified as *exactly invariant* — correctly, since it reads no radial operator. Being safe from the operator correction is not the same as being unimportant, and its elasticity had never been measured.

## R2 — effective rank: four observables, eleven knobs

`rank J = 4` against `8` first-order knobs, so **4 directions in these coordinates move no observable at all**; `JᵀJ` carries `4` exact zeros. The identifiable tangent space is 4-dimensional — the *model* has more local non-identifiability still, since the two quadratic coordinates are held fixed here.

| `a` | `σ_a` | `σ_a/σ_1` |
|--|--|--|
| 1 | `19.2190` | `1.0000` |
| 2 | `7.5739` | `0.3941` |
| 3 | `2.5349` | `0.1319` |
| 4 | `0.8519` | `0.0443` |

The four **nonzero** singular values span only `22.6×` — no long hierarchy among the identifiable directions. That is *not* the same as "not a sloppy model": the four nonzero singular values span only 22.6x, but J^T J carries 4 exact zeros in these 8 coordinates and the full 11-coordinate model carries more; the problem is rank, not conditioning.

*Scope: Euclidean in the 8 positive knobs in unit log coordinates, quadratic coordinates held fixed. Rescaling a column changes every right-space quantity below.*

**Right singular vectors** (parameter combinations):

| knob | `v1` | `v2` | `v3` | `v4` |
|--|--|--|--|--|
| `beta` | `+0.0393` | `+0.0453` | `+0.2526` | `+0.7971` |
| `gamma_q` | `-0.0607` | `+0.0180` | `-0.8641` | `+0.1043` |
| `transport` | `+0.0254` | `-0.0100` | `+0.3666` | `-0.0459` |
| `pinhole` | `-0.0954` | `-0.7140` | `+0.1173` | `-0.3795` |
| `resistance` | `-0.0150` | `-0.0319` | `-0.1522` | `+0.0013` |
| `uplift_asymmetry` | `-0.9880` | `+0.1344` | `+0.0741` | `+0.0179` |
| `chi_q_k3` | `+0.0930` | `+0.6845` | `+0.1119` | `-0.4554` |
| `eta_k3k5_minus` | `+0.0025` | `+0.0123` | `+0.0096` | `-0.0022` |

**Left singular vectors** (species response):

| species | `u1` | `u2` | `u3` | `u4` |
|--|--|--|--|--|
| `s` | `-0.0138` | `-0.8530` | `+0.5216` | `-0.0144` |
| `c` | `-0.0130` | `+0.3584` | `+0.5652` | `-0.7429` |
| `b` | `+0.9992` | `+0.0064` | `+0.0373` | `+0.0140` |
| `t` | `-0.0357` | `+0.3794` | `+0.6380` | `+0.6691` |

The quadratic directions add no new observable direction — with rank 4 and four observables the column space is already all of `ℝ⁴`, so *any* 4-vector projects exactly. The content is **which** combination each mimics:

| knob | dominant equivalent first-order knob |
|--|--|
| `phase` | `gamma_q` |
| `partition_mixing` | `gamma_q` |

## R3 — which singular directions actually repair the masses

Reporting `σ_a` alone says nothing. `c_a = u_aᵀr` and `c_a/σ_a` say everything: a tiny `σ` is harmless if the residual has no projection on it, and a moderate one can dominate.

| `a` | `σ_a` | `c_a = u_aᵀr` | `c_a/σ_a` | share of `\|δx\|²` |
|--|--|--|--|--|
| 1 | `19.2190` | `-0.007223` | `-0.000376` | `0.0003` |
| 2 | `7.5739` | `+0.015230` | `+0.002011` | `0.0077` |
| 3 | `2.5349` | `-0.005377` | `-0.002121` | `0.0086` |
| 4 | `0.8519` | `+0.019377` | `+0.022745` | `0.9835` |

**98.3% of the repair rides `v4` — the *weakest* direction.** *(Scope: the share decomposition and delta_x_min are Euclidean in unit log coordinates; rescaling a column changes them. The exact re-run does not depend on that choice.)*

The correction has `|δ ln p| = 0.022935`, its largest single knob moving `1.78%`, and an **exact nonlinear re-run** — `exp(δ)` applied to the lock, spectrum re-solved — gives `0.0179%` from the lock's `1.58%`.

> The frozen Hamiltonian has enough local freedom to absorb essentially the whole residual: 1.58% -> 0.0179% under a displacement whose largest single knob moves 1.78%. That is a statement about parameter count, not about the physics.

## R4 — does the operator correction move toward the data?

Two different questions live here, and the first draft ran them together.

**(a) Which candidate displacement from the lock lands nearer the data?**

| operator | `\|δx_geom\|` | `cos(J·g, r)` | L2 log-residual | max rel err | null-space fraction |
|--|--|--|--|--|--|
| legacy | `0.014819` | `+0.4638` | `0.0548` | `3.19%` | `0.6630` |
| corrected | `0.007906` | `-0.6163` | `0.0433` | `3.80%` | `0.6922` |

**(b) Which way does the *correction itself* move?** That is `Δg = g_corrected − g_legacy`, compared against the residual the legacy triple leaves.

| `\|Δg\|` | `\|r − J·g_legacy\|` | `cos(J·Δg, r − J·g_legacy)` |
|--|--|--|
| `0.015537` | `0.054813` | **`+0.8734`** (`29.1°`) |

> **The first draft's headline -- 'the corrected geometry moves AWAY from what the masses want', from cos(J g_corrected, r) = -0.616 -- is a true statement about g_corrected as a displacement FROM THE LOCK, and an invalid basis for a claim about what correcting the operator does. Withdrawn.**

And the two objectives disagree, which is substance rather than bookkeeping: the L2 log-residual **improves** (`0.0548` → `0.0433`) while the repository's max-relative-error score **worsens** (`3.19%` → `3.80%`). Both are true; the direction of "improvement" is objective-dependent and no metric-free claim is available.

**What survives.** Per-knob proximity to a fitted value carries no information about the sign or size of the joint effect on the spectrum -- which was the point, and does not need the sign of any cosine.

## R5 — the holdout measures the regulariser, not the Hamiltonian

For each holdout `J_keep` is `3×8` with a five-dimensional kernel, and because the full `J` has rank 4 the held-out row is **not** in the span of the other three. So a `z ∈ ker(J_keep)` moving the withheld species always exists.

| held out | `dim ker` | `\|j_h·Z\|` | min-norm miss | after a kernel shift | norm cost |
|--|--|--|--|--|--|
| `s` | `5` | `4.2520` | `+1.912%` | **`+6.9e-16%`** | `1.02×` |
| `c` | `5` | `1.1095` | `+2.519%` | **`-3.5e-16%`** | `4.77×` |
| `b` | `5` | `17.6989` | `-9.705%` | **`-3.6e-15%`** | `1.03×` |
| `t` | `5` | `1.2102` | `-2.521%` | **`+0.0e+00%`** | `2.55×` |

> **This is a property of the minimum-log-norm regulariser, not of the Hamiltonian. The rank deficiency already establishes the local flexibility; reading these misses as failed predictions was an overclaim and is withdrawn.**

The first draft read the `−10.4%` as a failed prediction. It is not one: `δ + λz` fits the withheld mass to machine precision while holding the other three exact, at between `1.02×` and `4.77×` the minimum norm.

## R6 — the ledger

| claim | verdict | evidence |
|--|--|--|
| four scored masses can identify the current quark parameterisation | **WITHDRAWN** | rank J = 4 against 8 first-order knobs; J^T J carries 4 exact zeros in these coordinates and the full 11-knob model more |
| action_base is a free parameter of the mass spectrum | **WITHDRAWN — EXACT GAUGE** | H(a) = H(0) + a*I, removed by the spectrum-zero subtraction upstream of the d-anchor; flat to all orders |
| phase and partition_mixing are null directions | **MISCLASSIFIED** | two different Z2s -- phase antiunitary on the spectrum, partition_mixing unitary-conjugation -- both quadratic; the Jacobian sees them in q = x^2, not in x |
| the phase Z2 is a property of the lock, not the model | **WITHDRAWN — THIS ROUND'S OWN ERROR** | H(-phi, p) = conj(H(phi, p)) for arbitrary p, so the SPECTRUM is even at every mixing (0.0 on eigenvalues and anchored masses); the first draft mistook matrix inequality for spectral asymmetry |
| per-knob proximity to a fitted value determines the effect on the spectrum | **REFUTED** | the three radial residuals do not compose linearly into the ladder; individual agreement carries no information about the sign or size of the joint effect |
| the corrected geometry moves AWAY from what the masses want | **WITHDRAWN — THIS ROUND'S OWN ERROR** | that used g_corrected as a displacement from the lock. The move the correction makes is Delta g = g_corrected - g_legacy, and cos(J.Delta g, r - J.g_legacy) = +0.873 -- toward the residual, not away |
| one objective settles whether the correction helped | **NOT AVAILABLE** | L2 log-residual improves (0.0548 -> 0.0433) while max relative error worsens (3.19% -> 3.80%); the direction of improvement is objective-dependent |
| the v3 ladder's 1.61% floor is set by the functional form | **WITHDRAWN** | a displacement whose largest single knob moves 1.78% reaches 0.0179% on an exact nonlinear re-run |
| that repair is evidence for the model | **NOT ESTABLISHED** | 98.3% of it rides the weakest singular direction in the chosen log chart, and the rank deficiency makes some such displacement guaranteed to exist |
| leave-one-species-out shows the Hamiltonian fails to predict | **WITHDRAWN — THIS ROUND'S OWN OVERCLAIM** | ker(J_keep) is 5-dimensional and the held-out row is reachable from it in every case; a kernel shift fits the withheld mass to ~1e-15 at 1.02-4.77x the norm. It measures the regulariser, not the model |
| the quark model is 'sloppy' in the Sethna sense | **NOT THE RIGHT DIAGNOSIS** | the four nonzero singular values span only 22.6x -- no long hierarchy -- but the rank deficiency is exact; the pathology is structural non-identifiability, not ill-conditioning |
| N = 466 drifting is a defect of N, and any knob would drift | **REFRAMED, AND THE SECOND HALF IS FALSE** | identifiable share runs from 1.0000 (uplift_asymmetry) to 0.0003 (eta_k3k5_minus); which knobs drift is itself structure |
| the missing correlation is a scalar relation R = f(p, T) | **REFUTED** | the degeneracy is a linear subspace selected by the response map; the nearest pair is gamma_q/transport at 178.9 deg, not transport/resistance |

**Headline.** Four scored masses cannot identify the current quark parameterisation. The positive-knob Jacobian is numerically full row rank, and the radial residual triple's individual proximity to fitted knobs is not enough to infer its effect on the spectrum.

**What would settle it.** More observables, not more knobs. The v4 flavor-CP layer supplies CKM angles and J from the same Hamiltonian and inherits the v3 eigenvalues, so it adds observable rows -- but it is NOT automatically column-free: v4 QuarkParams carries phi_h, targeted eta couplings and per-shell diagonal shifts. A joint identifiability audit has to decide which of those are externally derived and which belong in the response columns.
