# The local geometry of the quark fit manifold

**7/7 checks pass.**

PR #272 measured one knob at a time and named the gap between individually-right and jointly-wrong residuals a *missing correlation*, guessing it lived between `transport` and `resistance`. **That guess was wrong, and the object is not a scalar relation.** It is the response map `J_ia = ∂ln(m_i/m_d)/∂ln p_a` and its SVD.

## R0 — three apparent nulls, three different objects

| direction | status | mechanism | check |
|--|--|--|--|
| `action_base` | **EXACT INVARIANCE (all orders)** | `H(a) = H(0) + a·I`, killed by the spectrum-zero subtraction | `1.8e-13` |
| `phase` | **Z2-EVEN, QUADRATIC** | enters only as `cos(φ·dk)` — *because* `partition_mixing = 0` | `0.0` |
| `partition_mixing` | **Z2-EVEN, QUADRATIC (unitary equivalence)** | `H(−p) = D H(p) D†`, `D = diag(σ)` — isospectral | `0.0` |

`action_base` is removed **upstream of the `d`-anchor** — the zero-point subtraction does it, and the alternative `spectrum_zero_mode="action_base"` kills it too. So it does *not* reappear in an absolute-scale observable of this model. It is dropped from the first-order parameter space.

The `phase` Z₂ is a property of the **lock**, not the model: switch `partition_mixing` on and the evenness breaks (`0.096`).

Both quadratic directions are handled by the chart `q = x²`, under which `∂ln m/∂q` is finite and constant:

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

`rank J = 4` against `8` first-order knobs, so **4 directions move no observable at all** — before any accidental degeneracy.

| `a` | `σ_a` | `σ_a/σ_1` |
|--|--|--|
| 1 | `19.2190` | `1.0000` |
| 2 | `7.5739` | `0.3941` |
| 3 | `2.5349` | `0.1319` |
| 4 | `0.8519` | `0.0443` |

Condition number `22.6`. **This is not a sloppy model** in the Sethna sense — a sloppy spectrum spans 10³–10⁶. The identifiable subspace is well conditioned; the degeneracy is *dimensional*: 8 first-order knobs (plus 2 quadratic and 1 exactly gauge) against 4 scored masses.

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

**98.3% of the repair rides `v4` — the *weakest* direction.** The min-norm repair is not a coherent physical combination; it is dominated by the softest singular direction — by the model's least-constrained freedom.

The correction itself has `|δ ln p| = 0.022935` and drives the linear residual to `1.6e-17` (exact re-run: `0.018%`, from the lock's `1.61%`). **That is not a success of the physics.** It shows the frozen Hamiltonian has enough local freedom to absorb essentially the whole residual — the old `1.61%` floor was set by *where the lock sits*, not by the functional form.

## R4 — does the geometry move where the data want?

The decisive test. Project the operator-induced displacement `δx_geom` onto the mass-optimal direction `δx_min = J⁺r`.

| operator | `\|δx_geom\|` | `cos Θ` (parameter) | `cos Θ` (observable) | invisible fraction |
|--|--|--|--|--|
| legacy | `0.014819` | `+0.2870` | **`+0.4638`** (`62.4°`) | `0.6630` |
| corrected | `0.007906` | `-0.2785` | **`-0.6163`** (`128.0°`) | `0.6922` |

> **Case three of the trichotomy: the corrected geometry moves the Hamiltonian AWAY from the mass-optimal direction (cos = -0.616) where the legacy geometry had partial overlap with it (cos = +0.464). PR #272's three per-knob improvements were not merely uninformative about the ladder — they were misleading about it.**

And About two thirds of each geometric displacement lies in the null space of J and moves no observable at all — so even the part that points somewhere mostly points nowhere observable.

## R5 — leave-one-species-out

Fit the minimum-norm correction on three masses; evaluate the fourth with **nothing readjusted**.

| held out | `\|δx\|` | fitted three (exact) | **held-out (exact)** |
|--|--|--|--|
| `s` | `0.022498` | `0.015%` | **`+1.905%`** |
| `c` | `0.004813` | `0.041%` | **`+2.562%`** |
| `b` | `0.022198` | `0.059%` | **`-10.410%`** |
| `t` | `0.008997` | `0.002%` | **`-2.519%`** |

Local flexibility, not structure: the four-parameter repair that drives the full ladder to 0.018% mispredicts a withheld species by up to 10.4%. The fitted species land at `0.002–0.06%` every time; the withheld one does not come along.

## R6 — the ledger

| claim | verdict | evidence |
|--|--|--|
| action_base is a free parameter of the mass spectrum | **WITHDRAWN — EXACT GAUGE** | H(a) = H(0) + a*I, removed by the spectrum-zero subtraction upstream of the d-anchor; flat to all orders |
| phase and partition_mixing are null directions | **MISCLASSIFIED** | Z2-even at the lock, quadratic away from it; the Jacobian cannot see them in x, only in q = x^2 |
| pinhole, transport, resistance are independently constrained by the quark masses | **WITHDRAWN** | rank J = 4 against 8 first-order knobs, so 4 directions move no observable; the masses fix at most four combinations of everything |
| the v3 ladder's 1.61% floor is set by the functional form | **WITHDRAWN** | a min-norm displacement of |dln p| = 0.0229 reaches 0.018%; the floor was where the lock sits, not what the model can do |
| that 0.018% fit is evidence for the model | **NOT ESTABLISHED** | 98.3% of the repair rides the weakest singular direction, and leave-one-out mispredicts by up to 10.4% |
| PR #272's three per-knob improvements moved the ladder toward the data | **REFUTED** | cos(Theta) in observable space is -0.616 (corrected) against +0.464 (legacy): the correction moves AWAY from what the masses want |
| the quark model is 'sloppy' in the Sethna sense | **NOT SUPPORTED** | condition number 22.6 over the identifiable subspace — well conditioned; the degeneracy is dimensional, not ill-conditioning |
| N = 466 drifting is a defect of N | **REFRAMED** | it is the visible symptom of a four-dimensional observable space; any knob would drift along the unconstrained directions |
| the missing correlation is a scalar relation R = f(p, T) | **REFUTED** | PR #272's conjecture; the degeneracy is a linear subspace selected by the response map, and the nearest pair is gamma_q/transport at 178.9 deg, not transport/resistance |

**Headline.** The quark masses constrain four combinations of eleven knobs; the residual is removable, the repair is a compensator, and the corrected geometry pushes against the data rather than with it.

**What would settle it.** More observables, not more knobs — the v4 flavor-CP layer already supplies CKM angles and J from the same Hamiltonian, which would raise the achievable rank above four and make these directions identifiable for the first time.
