# Does negative-coupling EGB actually work?

**5/5 checks pass — the branch closes.**

Frozen in `docs/negative_egb_prereg.md` before this module existed.

> ## The step the previous round missed

> `α_GB` is a **coupling constant in the action**, so the same value acts in the exterior the throat is glued into. PR #277 analysed the throat in isolation, and it should not have.

> **Negative-coupling EGB closes on PHYSICAL grounds: it survives at one exact coupling, and there the observed universe must be empty of ordinary matter.**

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

> H^i_j = 0 on a maximally symmetric spatial slice, so the entire Gauss-Bonnet correction lands in rho and none of it in p. The equation of state is therefore driven to -1 by shrinking rho onto a fixed p, not by adjusting both.

> **At the one surviving coupling the exterior -- which is supposed to BE the observed closed universe -- is pure vacuum energy with no ordinary matter anywhere, at half the Einstein-gravity density. Push alpha more negative and the EXTERIOR violates the NEC, which relocates the exotic matter from the throat into the universe around it rather than removing it.**

## E5 — the ledger

| claim | verdict | evidence |
|--|--|--|
| the throat's constraint on alpha_GB can be read alone | **NO -- alpha_GB IS A COUPLING CONSTANT** | the same value acts in the exterior, where R_kk = +3.0000 > 0 and the NEC needs alpha >= -0.250000 -- the opposite direction |
| some open interval of alpha_GB works everywhere | **NO -- THE SURVIVING SET IS MEASURE ZERO** | a 4001-point scan over [-2.0, 2.0] leaves a set of width 0.0e+00, at the single value -0.250000 |
| the two bounds coinciding is a coincidence | **NO -- ONE CONTINUOUS BRACKET** | T_kk = R_kk(1 + 4 alpha mu/f^4) on both sides, with R_kk < 0 inside and > 0 outside, and mu/f^4 continuous at the seam at 4.000000/4 |
| Gauss-Bonnet can adjust the exterior pressure | **NO -- H^i_j = 0** | on a maximally symmetric slice the whole correction lands in rho; p = -3/R^2 whatever alpha is |
| the surviving coupling is physically usable | **NO -- IT EMPTIES THE UNIVERSE** | w = -1.0 exactly, rho = 3.000000 = half the Einstein value: pure vacuum energy, no ordinary matter anywhere |

**What this closes.** Global existence for constant-coupling EGB on this geometry. The previous round narrowed the branch by showing a negative coupling was needed; this one closes it by noting that the coupling is not the throat's to choose.

**What remains untested.** Stability and the graviton kinetic term at the critical coupling -- the second half of PR #277's open item, and not addressed here. Also dilatonic alpha(phi) L_GB, where the scalar's own stress enters the null contraction and known 5D solutions exist, and f(R). And a DIFFERENT exterior: the constraint is derived for the round S^4_R completion this programme assumes, so a different bulk completion would need its own version of this calculation.

**The remaining branches.**

- **1 accept the horizon** — the Tangherlini branch N(0) = 0 as the particle, abandoning MTY traversability
- **2 ghost** — a wrong-sign field, with its stability problem
- **3 quantum stress** — Casimir-type support, so the geometry is no longer classical
- **4 reinterpret** — particle exchange needs no traversable throat
- **5 dilatonic EGB or f(R)** — untested here, and the place where the heterotic term actually lives
