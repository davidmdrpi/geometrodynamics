# The whole-throat S-matrix of a supported traversable 5D throat

**7/7 checks pass.**

`transaction/network.py` (PR #216) replaced `handshake.py`'s phenomenological advanced confirmation with a Morris–Thorne–Yurtsever mechanism — everything propagates forward in local time, and the mouths' clock offset carries the return into the exterior past. That was a real relocation of the postulate. But the geometry was never supplied: `MouthPort.t`, `r_in`, `r_out` are inputs, and `closure_offset` *solves* for the `Δ` that makes the loop close.

This round supplies the geometry and recomputes the closure without tuning `Δ` to the answer.

> ## Three results

> **1.** Traversability costs exactly `-0.187500/G₅` = `−3/(16G₅a)` in the complete null-ray energy integral.

> **2.** No single clock offset closes the loop: phase closure and group closure disagree by up to `6.78` over the sampled band.

> **3.** `Λ = 1` — PR #216's completed transaction — cannot occur at any finite frequency, and the deficit closes only exponentially, `1 − |T|² ~ exp(-4.25 ω)`.

---

## Why not Tangherlini — an argument, not a computation

`supp G_ret(c,s) ⊂ J⁺(s)` and `supp G_adv(c,d) ⊂ J⁻(d)`, so a nonzero product `G_ret(c,s)G_adv(c,d)` needs `c ∈ J⁺(s) ∩ J⁻(d)`. But then `s → c → d` is a causal chain, so `d ∈ J⁺(s)`. Contrapositive: if `d ∉ J⁺(s)` the product vanishes for **every** `c`. **An advanced leg by itself does not evade ER non-traversability.** PR #275's `T_ℓ` is exterior→*horizon*, and is exactly that zero for the cross-mouth channel.

## G0 — the geometry, derived

`f(s) = sqrt(s^2 + a^2), the scalar-flat solution of f'^2 = 1 - a^2/f^2`, giving `V_l = l(l+2)/f^2 + (3/4)(s^2 + 2a^2)/f^4`.

| `ℓ` | `V_ℓ(0)` | exact `[ℓ(ℓ+2)+3/2]/a²` | tail coeff | exact `(ℓ+1)²−¼` |
|--|--|--|--|--|
| 0 | `1.500000` | `1.500000` | `0.750000` | `0.750000` |
| 1 | `4.500000` | `4.500000` | `3.750000` | `3.750000` |
| 2 | `9.500000` | `9.500000` | `8.750000` | `8.750000` |
| 3 | `16.500000` | `16.500000` | `15.750000` | `15.750000` |

**V is a single smooth symmetric maximum at s = 0, with no interior cavity and no separated interfaces. The Fabry-Perot form t_A t_B / (1 - r_inA r_inB e^{2 i w tau}) that PR #216 assumes is therefore a modelling choice the geometry does not supply, and this module does not impose it: the whole throat is solved as one scattering problem.**

V -> [(l+1)^2 - 1/4]/s^2 at BOTH ends -- numerically the same centrifugal tail as the Tangherlini exterior of PR #275 -- so that round's validated Riccati-Hankel boundary condition is imported rather than rewritten.

*V_l(0) = [l(l+2) + 3/2]/a^2 carries the same l(l+2) + 3/2 that PR #275 derived as the exact int V dr* on Tangherlini.*

## G1 — the GR price

| `s` | `8πG₅ρ` | `8πG₅p_s` | `8πG₅p_Ω` |
|--|--|--|--|
| `0` | `0.000000` | `-3.000000` | `+1.000000` |
| `0.5` | `0.000000` | `-1.920000` | `+0.640000` |
| `1` | `0.000000` | `-0.750000` | `+0.250000` |
| `2` | `0.000000` | `-0.120000` | `+0.040000` |
| `5` | `0.000000` | `-0.004438` | `+0.001479` |

`ρ = 0` identically, and the radial NEC `ρ + p_s < 0` **everywhere**. Along a complete radial null geodesic with `k^t̂ = 1`:

```
∫ T_ab k^a k^b dλ = -0.187500/G₅ = -3/(16 G5 a)   (exact)
```

The exoticity scales as `1/(G₅a)` — a wider throat is cheaper.

> **This is the support required by the STATIC traversable throat. It is NOT the energy cost of the clock offset Delta_BA. A frozen final metric can carry a large accumulated offset without any local stress tensor proportional to it -- the offset is produced by the mouths' differential-aging HISTORY. Costing Delta requires moving-mouth dynamics or an explicit gravitational-well history, and this module does not attempt it.**

## G2 — the whole-throat S-matrix

| `ω` | `\|T\|` | `\|R\|` | `\|R\|²+\|T\|²−1` |
|--|--|--|--|
| `0.05` | `0.00098705` | `0.99999951` | `-1.6e-14` |
| `0.1` | `0.00399871` | `0.99999201` | `+1.8e-14` |
| `0.2` | `0.01665893` | `0.99986123` | `+9.8e-15` |
| `0.5` | `0.12651406` | `0.99196481` | `+2.3e-14` |
| `1` | `0.64142490` | `0.76718583` | `-1.4e-14` |
| `2` | `0.99657150` | `0.08273606` | `-3.7e-14` |
| `5` | `0.99999999` | `0.00011608` | `-2.6e-14` |
| `10` | `1.00000000` | `0.00000000` | `-6.6e-14` |

Worst unitarity residual `6.6e-14`, imposed nowhere. V is even in s to machine precision, so R_L = R_R and T_LR = T_RL identically. Reciprocity here is a property of the geometry, not a numerical coincidence to be verified.

## G3 — the threshold law, and the false regression avoided

> A flux-normalised transmission amplitude between asymptotically flat four-dimensional spatial ends vanishes at threshold, from the asymptotic radial normalisation. g = pi^2 a^2 is a static boundary response; asserting T_0(0) = g would be a false regression. What the static solution controls is the COEFFICIENT of the low-frequency expansion.

Static: `I₃ = 2.0000`, `g = 9.869604 = π²a²`. The dynamical law is instead:

```
|T_0(w)| -> (pi/8)(a w)^2 ,  measured/predicted -> i
```

| `ω` | `\|T₀\|/(π/8)(aω)²` | `arg(ratio)/π` |
|--|--|--|
| `0.01` | `1.000294` | `+0.5000` |
| `0.015` | `1.000618` | `+0.5001` |
| `0.02` | `1.001043` | `+0.5001` |
| `0.03` | `1.002168` | `+0.5002` |
| `0.05` | `1.005398` | `+0.5006` |

> The oracle was stated up to the overall phase convention for the asymptotic waves, and the measured ratio is a CONSTANT +pi/2 to within 5e-3 across the sampled band -- i.e. exactly a factor of i, from this module's Riccati-Hankel normalisation to e^{+-ix}. The magnitude law, which is the physical content, is confirmed.

## G4 — no single clock offset closes the loop

PR #216 sets `Δ_BA = −(d_A + d_B + τ_th)`, which *solves* for closure. Once the throat is dispersive there is no unique `τ_th`: the geometry supplies `δ_ℓ(ω) = arg T_ℓ(ω)`, and phase closure `Φ = 2πn` and group closure `dΦ/dω = 0` demand **different** offsets.

With `d_A + d_B = 3.1416` (antipodal on the unit `S³`):

| `ω` | `δ_ℓ` | `dδ/dω` (Wigner) | `Δ_phase` | `Δ_group` | gap |
|--|--|--|--|--|--|
| `0.5` | `-2.9260` | `+0.9323` | `+2.7105` | `-4.0739` | **`6.7844`** |
| `1` | `-2.1523` | `+1.9921` | `-0.9893` | `-5.1337` | **`4.1445`** |
| `1.5` | `-1.3339` | `+1.1324` | `-2.2523` | `-4.2740` | **`2.0217`** |
| `2` | `-0.9400` | `+0.5424` | `-2.6716` | `-3.6840` | **`1.0124`** |
| `3` | `-0.6032` | `+0.2117` | `-2.9405` | `-3.3533` | **`0.4128`** |
| `5` | `-0.3563` | `+0.0724` | `-3.0703` | `-3.2140` | **`0.1437`** |

> One clock offset can phase-close a monochromatic carrier while failing to return a localised packet to the same event. Delta_n and Delta_g are the two demands; where they differ, no single offset does both.

d^2 delta / dw^2 sets how narrow the packet must be for the monochromatic answer to survive: the group delay varies across the packet by delta''(w0) * bandwidth.

**closure_offset(d_A, d_B, tau_th) = -(d_A + d_B + tau_th), which SOLVES for the offset that makes the loop close. Here the offset is a derived demand of the geometry, and it is frequency dependent.**

## G5 — can the transaction actually complete?

| `ω` | `1 − \|T\|²` |
|--|--|
| `1.50` | `7.7582e-02` |
| `2.87` | `1.2601e-04` |
| `4.24` | `3.3516e-07` |
| `5.61` | `1.0817e-09` |
| `6.97` | `3.8178e-12` |

> **|Lambda| = |t_net| |eta_topo| and a barrier with V > 0 everywhere has |T| < 1 at every finite frequency. So 1 - Lambda cannot vanish and G_eff = G_0/(1 - Lambda) has no true pole -- unless some element outside the throat supplies gain. PR #216's completed transaction is a high-Q near-resonance, not an exact one.**

> 1 - |T|^2 ~ exp(-4.25 w), so the quality factor of the near-transaction grows exponentially with frequency and the frequency needed for a given Q grows only logarithmically. Exact closure is approached in the UV, which is a chronology-horizon concern rather than the benign resonance the phenomenological ports suggest.

| `D_loop` | regime |
|--|--|
| `D_loop > 0` | ordinary delayed feedback |
| `D_loop = 0` | closed-null / marginal chronology condition |
| `D_loop < 0` | return before emission (CTC regime) |

## G6 — the ledger

| claim | verdict | evidence |
|--|--|--|
| an advanced leg alone evades ER non-traversability | **NO -- PROVED, NOT COMPUTED** | supp G_ret(c,s) in J+(s) and supp G_adv(c,d) in J-(d), so a nonzero product needs c in J+(s) cap J-(d), whence s -> c -> d and d in J+(s). If d is not in J+(s) the product vanishes for every c |
| PR #216's throat transfer is supplied by a geometry | **NOW IT IS** | T_l(w), R_l(w) computed from the whole traversable throat instead of MouthPort.t, r_in, r_out being inputs |
| the throat factorises as two mouths plus a cavity | **NOT BY THIS GEOMETRY** | V_l is a single smooth symmetric barrier; the Fabry-Perot form is a modelling choice, and is not imposed here |
| the static conductance g is T_0(0) | **NO -- IT IS THE THRESHOLD COEFFICIENT** | |T_0| -> (pi/8)(a w)^2, magnitude ratio 1.0003 at w = 0.01, with a constant phase convention factor of i |
| the traversable throat is free | **NO** | rho = 0 but p_s = -3a^2/(8 pi G5 f^4), radial NEC violated everywhere, and the complete null-ray integral is exactly -0.187500/G5 = -3/(16 G5 a) |
| that NEC integral is the cost of the clock offset | **NO -- DIFFERENT QUESTION** | it is the support of the static throat; Delta_BA comes from the mouths' aging history, and a frozen metric can carry a large offset with no local stress proportional to it |
| one clock offset closes the loop | **NOT FOR A PACKET** | phase closure and group closure demand different offsets; worst disagreement 6.7844 over the sampled band |
| Lambda = 1 (the completed transaction) can occur | **NOT AT ANY FINITE FREQUENCY** | |T| < 1 for a positive barrier, so 1 - Lambda cannot vanish; but 1 - |T|^2 ~ exp(-4.25 w), so exact closure is approached only in the UV |

**What this round establishes.** PR #216 relocated the advanced wave from a bare postulate into a Morris-Thorne-Yurtsever geometry, which was real progress. This round supplies the geometry it was missing and finds that the relocation has a price (NEC violation, exactly -3/(16 G5 a)), that the closure it assumed is frequency dependent rather than a single tuned constant, and that its completed transaction is a high-Q limit rather than an attainable pole.

**What remains postulated.** That the throat is traversable at all, and that the mouths carry the required clock offset. The first now has an exact price; the second needs moving-mouth dynamics this round does not do.

**Still open.** The history that produces Delta_BA, and whether any BAM exchange kernel is meant to approximate the whole-throat T_l derived here.
