# Does negative-coupling EGB actually work?

**7/7 checks pass — the branch closes, on the classical field equations rather than the matter.**

Frozen in `docs/negative_egb_prereg.md` before this module existed.

> ## The step the previous round missed

> `α_GB` is a **coupling constant in the action**, so the same value acts in the exterior the throat is glued into. PR #277 analysed the throat in isolation, and it should not have. The NEC then pins one coupling — and at exactly that coupling the linearised classical operator degenerates.

> **Negative-coupling EGB closes on STRUCTURAL grounds: the NEC pins one coupling, and at exactly that coupling the principal symbol of the linearised classical equations loses its leading time derivative.**

---

## E1 — the exterior constrains `α_GB` in the opposite direction

For `ℝ × S⁴_R`: `R_kk = 3.000000 = 3/R²` — **positive**, because the exterior holds ordinary matter rather than flaring out — and `H_kk = 12.000000 = 12/R⁴`. The ratio `H_kk/R_kk = 4μ/f⁴ = 4.000000 = 4/R²` is **independent of `χ`** (`True`), so the exterior's constraint is a single number.

| `α_GB` | `8πG₅T_kk` outside | exterior NEC? |
|--|--|--|
| `+0.000000` | `+3.000000` | **yes** |
| `-0.125000` | `+1.500000` | **yes** |
| `-0.250000` | `+0.000000` | **yes** |
| `-0.500000` | `-3.000000` | **no** |

```
exterior needs  α_GB ≥ -0.250000
throat   needs  α_GB ≤ -0.250000
```

> In the throat R_kk < 0 because the neck must flare out; in the exterior R_kk = +3/R^2 because it holds ordinary matter. The NEC therefore demands opposite signs of the SAME bracket 1 + alpha H_kk/R_kk, and the two bounds land on one number.

## E2 — searching for a coupling that works everywhere

A 4001-point scan over `[-2.0, 2.0]`. A surviving *interval* would reopen the branch, so this looks for one rather than asserting none.

Survivors: **`[-0.25]`**, a set of width `0.0e+00`.

> **There is no open set of couplings for which this spacetime satisfies the matter NEC everywhere. A fundamental constant of the action would have to take one exact value, tuned to the radius of the universe.**

## E3 — the two bounds meet by continuity

| side of the seam | `H_kk/R_kk` |
|--|--|
| throat, at `s = S` | `4.000000000` |
| exterior, at `χ = a` | `4.000000000` |

Jump `5.3e-15`, and at the critical coupling the shared bracket is `1.3e-15` — exactly zero.

> mu/f^4 is 1/R^2 on both sides of the seam -- the same Misner-Sharp continuity finite_mouth's P2 established -- so the shared bracket 1 + 4 alpha mu/f^4 is continuous there. It is required to be <= 0 from the throat side and >= 0 from the exterior side, so it must be exactly 0. The two constraints meet because they are two one-sided limits of one continuous function, not because two independent numbers coincided.

## E4 — what the one surviving coupling costs

| `α_GB` | `8πG₅ρ` | `8πG₅p` | `ρ+p` | `w = p/ρ` |
|--|--|--|--|--|
| `+0.000000` | `+6.000000` | `-3.000000` | `+3.000000` | `-0.5000` |
| `-0.125000` | `+4.500000` | `-3.000000` | `+1.500000` | `-0.6667` |
| `-0.250000` | `+3.000000` | `-3.000000` | `+0.000000` | `-1.0000` |
| `-0.375000` | `+1.500000` | `-3.000000` | `-1.500000` | `-2.0000` |

> H^i_j = 0 for ANY ultrastatic product -- the spatial block is the 4D Euler tensor, topological in D = 4 -- so the entire Gauss-Bonnet correction lands in rho and none of it in p. The equation of state is therefore driven to -1 by shrinking rho onto a fixed p, not by adjusting both.

> **At the one surviving coupling the exterior -- which is supposed to BE the observed closed universe -- is pure vacuum energy with no ordinary matter anywhere, at half the Einstein-gravity density. Push alpha more negative and the EXTERIOR violates the NEC, which relocates the exotic matter from the throat into the universe around it rather than removing it.**

## E4b — but the throat matter there is not exotic

> **An earlier draft of this round said the exotic matter was 'merely relocated' at the critical coupling. It is not. rho + p_s = 3A(q-1) >= 0 with equality only at the mouth, rho + p_Omega = A(3q+1) > 0, and rho exceeds both |p_s| and |p_Omega|, so NEC, WEC and DEC all hold along the throat. Relocation happens only for alpha strictly below the critical value.**

`q` runs `1.0000` (mouth) to `131.1147` (neck). Minima along the throat: `ρ+p_s = +2.7e-15`, `ρ+p_Ω = +4.0000`, `ρ−|p_s| = +2.7e-15`. NEC **True**, WEC **True**, DEC **True**.

## E4c — the closure: the classical principal symbol

**This is a classical PDE test, not a quantum one.** Linearise the classical field equations `G_AB + α H_AB` about this background, take `h_AB` transverse-traceless, and read off the highest-derivative operator. The background is a product, not a maximally symmetric spacetime, so the textbook coefficient does not apply and the symbol is derived here.

| `α_GB` | `C_t` (`ω²`) | predicted | `C_s` (`κ²`) | `deg_ω P` | hyperbolic | `c²` |
|--|--|--|--|--|--|--|
| `-0.00000` | `-0.4999995` | `-0.5000000` | `+0.4999997` | `2` | yes | `1.00` |
| `-0.05000` | `-0.3999997` | `-0.4000000` | `+0.4999998` | `2` | yes | `1.25` |
| `-0.12500` | `-0.2500000` | `-0.2500000` | `+0.4999998` | `2` | yes | `2.00` |
| `-0.20000` | `-0.1000003` | `-0.1000000` | `+0.4999999` | `2` | yes | `5.00` |
| `-0.24000` | `-0.0200005` | `-0.0200000` | `+0.5000000` | `2` | yes | `25.00` |
| `-0.25000` | `-0.0000005` | `-0.0000000` | `+0.5000000` | `0` | **no** | ∞ |

```
P(omega, kappa) = C_t omega^2 + C_s kappa^2
C_t = -(1/2)(1 + 4 alpha/R^2),   C_s = +1/2
c^2 = 1/(1 + 4 alpha/R^2)
```

That `P` is this quadratic form is **measured, not assumed** — directions off the coordinate axes reproduce `C_t d_t² + C_s |d_space|²`:

| propagation direction | measured | predicted | error |
|--|--|--|--|
| spatial, along x1 | `+0.4999998` | `+0.4999998` | `0.0e+00` |
| spatial, along x4 | `+0.4999998` | `+0.4999998` | `7.2e-16` |
| spatial, 45 deg in the (x1, x4) plane | `+0.5000001` | `+0.4999998` | `2.7e-07` |
| mixed, 45 deg in the (t, x1) plane | `+0.1250000` | `+0.1249999` | `6.7e-08` |
| mixed, 30 deg from t toward x4 | `-0.0625000` | `-0.0625000` | `5.0e-08` |

> No quantisation anywhere. This is the principal symbol of the linearised classical field equations, and the question is well-posedness of the Cauchy problem for classical metric perturbations. An earlier draft of this round called it a 'graviton kinetic term', which imported a quantum-gravity framing this programme rejects; the wording is retracted, the calculation is unchanged.

> The familiar 1 + 2 alpha (D-3)(D-4) K coefficient is derived for a maximally symmetric SPACETIME. R x S^4_R is a product and is not maximally symmetric, so that formula does not apply here and the coefficient was computed by linearising the full field equations on this background.

> **At alpha = -R^2/4 the omega^2 coefficient vanishes while the kappa^2 coefficient is untouched, so the principal symbol drops from degree 2 to degree 0 in omega: the linearised equation loses its leading time derivative and degenerates from an evolution equation into a constraint. The lower-order term stays finite, so this is a degeneration of the principal part specifically, not the whole equation vanishing or an overall factor that could be divided out.**

> **c^2 = 1/(1 + 4 alpha/R^2) exceeds 1 for every -R^2/4 < alpha < 0 and diverges at the critical coupling, so the tensor characteristic cone lies outside the metric null cone well before the degeneration. That is a causality failure rather than ill-posedness -- the operator is still hyperbolic there -- so the branch is in trouble on an interval, for one reason, and at the endpoint for a second and worse one.**

## E5 — the ledger

| claim | verdict | evidence |
|--|--|--|
| the throat's constraint on alpha_GB can be read alone | **NO -- alpha_GB IS A COUPLING CONSTANT** | the same value acts in the exterior, where R_kk = +3.0000 > 0 and the NEC needs alpha >= -0.250000 -- the opposite direction |
| the two bounds coinciding is a coincidence | **NO -- ONE CONTINUOUS BRACKET** | T_kk = R_kk(1 + 4 alpha mu/f^4) on both sides, with R_kk < 0 inside and > 0 outside, and mu/f^4 continuous at the seam at 4.000000/4 |
| the NEC alone leaves an open set of couplings | **NO -- IT PINS ONE** | a 4001-point scan leaves width 0.0e+00 at -0.250000; equivalently the field equations select R^2 = -4 alpha_GB |
| at that coupling the throat matter is exotic | **NO -- NEC, WEC AND DEC ALL HOLD** | rho = 3Aq, p_s = -3A, p_Omega = A with q >= 1, so rho+p_s = 3A(q-1) >= 2.7e-15 (saturated at the mouth) and rho exceeds both pressures in magnitude. An earlier draft's 'relocated exotic matter' applies only below critical |
| at that coupling the observed universe must be empty | **OVERREACHED -- IT IS A 5D BULK STATEMENT** | w = -1.0 is the total stress supporting R x S^4_R, not the S^3 equator's; and a homogeneous -Lambda g_ab can be moved to the gravitational side. The defensible claim is a vacuum-energy-like 5D exterior |
| the linearised tensor sector is well posed at the critical coupling | **NO -- THE PRINCIPAL SYMBOL DEGENERATES** | C_t = -(1/2)(1 + 4 alpha/R^2) derived by linearising on THIS product background, matching to 5e-07; it is -5.0e-07 at criticality while C_s is untouched, so P drops from degree 2 to degree 0 in omega and the Cauchy problem stops being an evolution problem |
| the symbol being a quadratic form is an assumption | **NO -- IT IS MEASURED OFF-AXIS** | propagation directions off the coordinate axes, including mixed time-space ones, reproduce C_t d_t^2 + C_s |d_space|^2 to 3e-07 |
| the trouble starts only at the critical point | **NO -- THE CONE LEAVES THE NULL CONE FIRST** | c^2 = 1/(1 + 4 alpha/R^2) exceeds 1 for every -R^2/4 < alpha < 0, reaching 25 at 96% of the critical value before diverging. There the operator is still hyperbolic -- it is a causality failure, not ill-posedness; the two are distinct and an earlier draft ran them together |

**What this closes.** Constant-coupling EGB on this geometry, on STRUCTURAL rather than matter-content grounds. The NEC pins a unique coupling; the principal symbol of the linearised classical equations degenerates at exactly that coupling; and the tensor characteristic cone is already outside the metric null cone on the approach.

**What remains untested.** Whether an admissible SOURCE realises the throat's full anisotropic stress -- this round determines the stress the metric requires, not fields obeying their own equations that supply it. Also the scalar and vector perturbation sectors -- this round settles the tensor sector only -- dilatonic alpha(phi) L_GB where the scalar's own stress enters and where the heterotic term lives, f(R), and a DIFFERENT exterior: the constraint is derived for the round S^4_R completion this programme assumes.

**The remaining branches.**

- **1 accept the horizon** — the Tangherlini branch N(0) = 0 as the particle, abandoning MTY traversability
- **2 ghost** — a wrong-sign field, with its stability problem
- **3 quantum stress** — Casimir-type support, so the geometry is no longer classical
- **4 reinterpret** — particle exchange needs no traversable throat
- **5 dilatonic EGB or f(R)** — untested here, and the place where the heterotic term actually lives
