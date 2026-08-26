# Scalar operator audit — correcting it, and pricing the correction

**Question.** the radial potential was short of the minimally coupled scalar master operator by 3A^2/(4r^2) -- which published claims are algebraically untouched, which keep their meaning with different digits, and which no longer say what they said?

**Answer.** the eigenvalues move at the 1e-3 level and less as l rises, and the cross-l operator is exactly invariant; but the barrier sums move enough that the two gamma statements in the tree swap places -- the canonical l = 1..5 claim improves threefold to -0.75% with nothing tuned, while the claim that the l = 0 channel closes the gap overshoots to +1.50% and is withdrawn

**8/8 checks pass.**

| id | check | result |
|----|-------|--------|
| T0 | the correction, with three independent confirmations | PASS |
| T1 | the eigenvalues barely move, and less as l rises | PASS |
| T2 | *** the barrier sums move, and the gamma story swaps *** | PASS |
| T3 | the cross-l operator survives exactly; its elements drift | PASS |
| T4 | ratios absorb most of the shift | PASS |
| T5 | actions move more than eigenvalues | PASS |
| T6 | *** the ledger: invariant / shifted / reinterpreted *** | PASS |
| T7 | *** the one narrow downstream re-derivation: gamma is the selector *** | PASS |

## T0 — the correction

| | operator |
|--|--|
| legacy | `A[l(l+2)/r^2 + 3 r_h^2/r^4]` |
| corrected | `A[l(l+n-1)/r^2 + n(n-2)A/(4r^2) + n A'/(2r)]` |
| gap | `3 A^2 / (4 r^2)` |

Three independent confirmations: the gap matches that closed form to `2.4e-15`; it carries no `ℓ` to `2.4e-15`; and the flat limit reproduces the Bessel form to `2.2e-16`, which is what settles *which* operator is the scalar one.

`V_scalar_tangherlini` agrees bitwise with `dynamics.master_potential`: `True`.

**The old name implied the canonical minimally coupled scalar operator, and the implementation was short of it by an l-independent term.**

## T1 — the eigenvalues barely move

| `ℓ` | legacy `ω_{ℓ,0}` | corrected | shift | min overlap |
|--|--|--|--|--|
| 0 | `1.00065891` | `1.00198000` | `+0.1320%` | `0.999998` |
| 1 | `1.05472694` | `1.05582653` | `+0.1043%` | `0.999998` |
| 2 | `1.13156946` | `1.13239953` | `+0.0734%` | `0.999998` |
| 3 | `1.21908274` | `1.21966785` | `+0.0480%` | `0.999999` |
| 4 | `1.30869618` | `1.30909388` | `+0.0304%` | `0.999999` |
| 5 | `1.39597349` | `1.39624205` | `+0.0192%` | `0.999999` |

The monotone fall with `ℓ` is not a coincidence: an eigenvalue averages the potential against a bound state, so an l-independent shift matters least where the centrifugal term already dominates.

## T2 — the barrier sums, and where the meaning moves

| channels | legacy sum | residual | corrected sum | residual | `R_OUTER` legacy | corrected |
|--|--|--|--|--|--|--|
| `l = 1..5` | `22.00824` | `-2.19%` | `22.33119` | `-0.75%` | `1.28737` | `1.26788` |
| `l = 0..5` | `22.45268` | `-0.21%` | `22.83642` | `+1.50%` | `1.26227` | `1.24614` |

The two `γ` statements move in **opposite directions**. The canonical README claim improves from `-2.19%` to `-0.75%` — and the canonical l = 1..5 statement moves the right way on its own; nothing was tuned. The `ℓ = 0…5` claim goes from `-0.21%` to `+1.50%`, and the sum closest to `22.5` swaps from **l = 0..5** to **l = 1..5**.

**What has to be reopened:** the claim that adding the l = 0 5D channel closes the gamma discrepancy -- under the corrected operator it overshoots, and the l = 1..5 sum is the one near 22.5.

## T3 — what survives exactly

`ΔV` carries no `ℓ`, so `V_{ℓ+2} − V_ℓ` is unchanged to `1.8e-15` — algebraically exact. Its matrix elements are not:

| `ℓ` | element legacy | corrected | drift |
|--|--|--|--|
| 0 | `3.223657e-01` | `3.212593e-01` | `-0.343%` |
| 1 | `4.442945e-01` | `4.428735e-01` | `-0.320%` |
| 2 | `5.186418e-01` | `5.171943e-01` | `-0.279%` |
| 3 | `5.488586e-01` | `5.476245e-01` | `-0.225%` |

**The operator v_{l+2} - v_l survives algebraically exactly; its matrix elements drift because the eigenfunctions do.**

## T4 — ratios absorb most of the shift

| `ℓ` | `α_q` legacy | corrected | drift |
|--|--|--|--|
| 0 | `+0.954557` | `+0.954802` | `+0.026%` |
| 1 | `+1.000000` | `+1.000000` | `+0.000%` |
| 2 | `+1.066078` | `+1.065786` | `-0.027%` |
| 3 | `+1.141558` | `+1.140895` | `-0.058%` |
| 4 | `+1.219166` | `+1.218249` | `-0.075%` |
| 5 | `+1.294669` | `+1.293542` | `-0.087%` |

*A ratio is not automatically safe -- the common part cancels, the differential does not, and this is measured rather than assumed.*

## T5 — actions move more than eigenvalues

| `ℓ` | action legacy | corrected | drift |
|--|--|--|--|
| 1 | `2.770577` | `2.736322` | `-1.236%` |
| 2 | `2.496428` | `2.494327` | `-0.084%` |
| 3 | `2.421247` | `2.421205` | `-0.002%` |

The action integrates sqrt(w^2 - v) against the potential directly, without the bound state's averaging.

## T6 — the ledger

Not whether the old tests stay green, but which published claims are algebraically untouched, which keep their meaning with different digits, and which no longer say what they said.

| claim | verdict | evidence |
|--|--|--|
| cross-l perturbation operator V_{l+2} - V_l | **EXACTLY INVARIANT** | unchanged to 1.8e-15 |
| Hopf fibration, Pin- structure, odd-k ladder, antipodal parity | **EXACTLY INVARIANT** | no dependence on the radial operator; not re-run, and proximity is not dependence |
| alpha_q(l,0) throat-flux ratios | **NUMERICALLY SHIFTED** | largest drift 0.087% |
| omega_{l,n} radial eigenfrequencies and eigenfunctions | **NUMERICALLY SHIFTED** | ground shifts 0.0192% to 0.1320%; overlaps > 0.99 |
| cross-l transport matrix elements | **NUMERICALLY SHIFTED** | operator exact, elements drift up to 0.343% |
| closed-orbit / WKB radial actions | **NUMERICALLY SHIFTED** | largest drift 1.236% |
| the 1.054 factor, omega(1,0) at the gamma-locked geometry | **NUMERICALLY SHIFTED** | 1.054727 -> 1.055827; the quoted 1.054 becomes 1.056, which exceeds the 0.04% Compton-bridge tolerance and so needs re-quoting, not re-deriving |
| pinhole gamma = Sum V_max[1..5] vs the locked 22.5 | **NUMERICALLY SHIFTED** | -2.19% -> -0.75% -- the canonical README claim IMPROVES, and nothing was tuned |
| R_OUTER as the gamma = 22.5 fixed point | **NUMERICALLY SHIFTED** | l=0..5 root 1.26227 -> 1.24614; l=1..5 root 1.28737 -> 1.26788 |
| 'adding the l = 0 5D channel closes the gamma discrepancy' | **INTERPRETATION CHANGED** | -0.21% -> +1.50%; the l=0..5 sum now OVERSHOOTS and the l=1..5 sum is the one near 22.5 -- the closest channel set swaps |
| any generation or mass result whose chain runs through gamma, R_OUTER, or a radial eigenvalue | **INTERPRETATION CHANGED** | inherits the gamma reopening above; each such chain must be re-derived from the corrected operator before its number is quoted again |

| verdict | count |
|--|--|
| EXACTLY INVARIANT | 2 |
| NUMERICALLY SHIFTED | 7 |
| INTERPRETATION CHANGED | 2 |

**Not re-run, and why.** Hopf, pin-, odd-k winding and antipodal parity have no dependence on the radial scalar operator; proximity is not dependence.

**Still open.** The gamma narrative. the l = 0 closure claim is withdrawn, not replaced -- the corrected l = 1..5 sum lands nearer 22.5 than the old one did, but that is an observation, not a derivation of why 22.5.

## T7 — the one narrow downstream re-derivation

A: R_OUTER = 1.26 fixed; B: enforce Sum[1..5] = 22.5; C: enforce Sum[0..5] = 22.5, each passed through the **locked** lepton Hamiltonian with **nothing retuned**.

| case | `R_OUTER` | `γ` | `m_μ` err | `m_τ` err |
|--|--|--|--|--|
| baseline, gamma = 22.5 | — | `22.50000` | `-0.04%` | `+0.12%` |
| legacy R=1.26, gamma[0..5] | `1.26000` | `22.45268` | `+3.78%` | `+4.04%` |
| legacy R=1.26, gamma[1..5] | `1.26000` | `22.00824` | `+63.78%` | `+65.58%` |
| A corrected R=1.26, gamma[1..5] | `1.26000` | `22.33119` | `+15.16%` | `+15.71%` |
| A corrected R=1.26, gamma[0..5] | `1.26000` | `22.83642` | `-20.52%` | `-20.89%` |
| B corrected root [1..5] = 22.5 | `1.26788` | `22.50000` | `-0.04%` | `+0.12%` |
| C corrected root [0..5] = 22.5 | `1.24614` | `22.50000` | `-0.04%` | `+0.12%` |

**B and C are bit-identical.** Compute_knotted_lepton_spectrum discards r_outer outright; the locked block consumes the geometry only through the scalar hard_pinhole_gamma. So the channel-set choice leaves no trace in any observable once `γ` is enforced — and the comparison cannot decide it.

> **`γ = 22.5` is the selector; `R_OUTER` is downstream of it.**

Fixing `R_OUTER` and letting `γ` float is what breaks the ladder: `+15.16%` and `-20.52%`, against the legacy geometry's `+3.78%`. So the correction **weakens** the geometry-supplies-`γ` story even while improving the `1..5` residual in isolation.

The reason is sensitivity: `d ln m_μ / d ln γ = -16.57`, so a sub-percent geometric residual is **not** a small residual in this chain.

**What this does not settle:** the channel-set question is not decidable by the lepton observables -- it has to be settled by what gamma MEANS geometrically, because the masses only ever saw the scalar.
