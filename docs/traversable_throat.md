# The whole-throat S-matrix of a supported traversable 5D throat

*Module `geometrodynamics/tangherlini/traversable_throat.py`, probe
`experiments/closure_ledger/traversable_throat_probe.py` (7/7), tests
`tests/test_traversable_throat.py` (28/28).*

---

## What this round is for

`transaction/network.py` (PR #216) did something genuinely good: it replaced
`handshake.py`'s advanced confirmation — which its own docstring calls
*"phenomenological … a conjugated ('advanced') phase assigned by fiat"* — with a
Morris–Thorne–Yurtsever mechanism. Every wave propagates forward in local time;
the mouths' clock offset `Δ_BA` carries the return into the exterior past. The
postulate moved from a bare nonlocal rule to a geometric one.

But the geometry was never supplied. `MouthPort.t`, `r_in`, `r_out` are inputs,
`tau_th` is a constant, and

```python
def closure_offset(d_A, d_B, tau_th):
    """The clock offset Delta_BA for which the network path returns
    exactly to the emission event."""
    return -(d_A + d_B + tau_th)
```

*solves* for the offset that makes the loop close. This round supplies the
missing physical input and recomputes the closure without tuning `Δ` to the
answer.

---

## Three results

**1. Traversability has an exact price.**

```
∫ T_ab k^a k^b dλ  =  −3/(16 G₅ a)
```

along a complete radial null geodesic — exact, and scaling as `1/(G₅a)`.

**2. No single clock offset closes the loop.** Phase closure and group closure
demand offsets differing by up to **`6.78`** over the sampled band.

**3. `Λ = 1` cannot occur at any finite frequency**, and the deficit closes only
exponentially, `1 − |T|² ~ exp(−4.25 ω)`.

---

## Why not Tangherlini — an argument, not a computation

The cross-exterior channel on a maximally extended bridge vanishes for a reason
needing no numerics:

```
supp G_ret(c,s) ⊂ J⁺(s)      supp G_adv(c,d) ⊂ J⁻(d)
```

so `G_ret(c,s)·G_adv(c,d) ≠ 0` requires `c ∈ J⁺(s) ∩ J⁻(d)`. But then
`s → c → d` is a causal chain, so `d ∈ J⁺(s)`. Contrapositive:

> **If `d ∉ J⁺(s)`, the product vanishes for every `c`.**

**An advanced leg by itself does not evade ER non-traversability.** PR #275's
`T_ℓ` is exterior→*horizon*, not exterior→exterior — exactly the zero this
argument predicts for the cross-mouth channel. Computing `G_sym = ½(G_ret+G_adv)`
would produce another zero and arbitrate nothing.

So the network mechanism needs a genuinely traversable throat. That is a choice,
and this round prices it.

---

## The geometry, derived

`ds² = −dt² + ds² + (s²+a²)dΩ₃²`, with `f = √(s²+a²)` the scalar-flat solution
of `f′² = 1 − a²/f²`. With `Φ = e^{−iωt}Y_ℓ(Ω)u(s)` and `ψ = f^{3/2}u`, removing
the first-derivative term contributes `(3/4)(f′/f)² + (3/2)f″/f`, giving

```
ψ″ + [ω² − V_ℓ(s)]ψ = 0 ,   V_ℓ = ℓ(ℓ+2)/f² + (3/4)(s² + 2a²)/f⁴
```

| `ℓ` | `V_ℓ(0)` | exact `[ℓ(ℓ+2)+3/2]/a²` | tail coeff | exact `(ℓ+1)²−¼` |
|--|--|--|--|--|
| 0 | `1.500000` | `1.500000` | `0.750000` | `0.750000` |
| 1 | `4.500000` | `4.500000` | `3.750000` | `3.750000` |
| 2 | `9.500000` | `9.500000` | `8.750000` | `8.750000` |
| 3 | `16.500000` | `16.500000` | `15.750000` | `15.750000` |

Two structural facts, both load-bearing:

**It is one barrier, not two.** `V` is a single smooth symmetric maximum at
`s = 0`, positive everywhere, strictly decreasing outward. The Fabry–Pérot form
`t_A t_B/(1 − r_inA r_inB e^{2iωτ})` that PR #216 assumes is a *modelling
choice the geometry does not supply*, and this module does not impose it — the
whole throat is solved as one scattering problem over `s ∈ (−∞,∞)`. If a later
geometry really does contain two separated interfaces, the factorisation becomes
derived rather than assumed.

**The tail is PR #275's.** `V → [(ℓ+1)² − ¼]/s²` at *both* ends — numerically
the same centrifugal form as the Tangherlini exterior — so that round's
validated Riccati–Hankel Jost condition (order `ν = ℓ+1`, normalised to
`e^{±ix}`) is imported rather than rewritten. Plane waves would be the wrong
basis at low `ω` for exactly the same reason as there.

*A small echo: `V_ℓ(0) = [ℓ(ℓ+2)+3/2]/a²` carries the same `ℓ(ℓ+2)+3/2` that
#275 derived as the exact `∫V dr*` on Tangherlini.*

---

## The GR price

From the 5D Einstein tensor in an orthonormal frame — and `R₅ = 0`, though the
spacetime is not vacuum:

| `s` | `8πG₅ρ` | `8πG₅p_s` | `8πG₅p_Ω` |
|--|--|--|--|
| `0` | `0.000000` | `−3.000000` | `+1.000000` |
| `0.5` | `0.000000` | `−1.920000` | `+0.640000` |
| `1` | `0.000000` | `−0.750000` | `+0.250000` |
| `2` | `0.000000` | `−0.120000` | `+0.040000` |

```
ρ = 0    exactly
p_s = −3a²/(8πG₅ f⁴)      p_Ω = +a²/(8πG₅ f⁴)   (each of three)
```

The radial NEC `ρ + p_s < 0` fails **everywhere**, and the complete null-ray
integral is exactly `−3/(16G₅a)`. Exoticity scales as `1/(G₅a)`: a wider throat
is cheaper, which is the only knob available.

> **This is not the cost of the clock offset.** It is the support required by
> the *static* traversable throat. A frozen final metric can carry a large
> accumulated `Δ_BA` with no local stress tensor proportional to it — the offset
> is produced by the mouths' differential-aging *history*. Costing `Δ` needs
> moving-mouth dynamics, and this round does not attempt it.

---

## The scattering matrix

Piecewise-constant transfer matrix vectorised over `ω`, the method #275
validated, with Jost matching at both ends.

| `ω` | `\|T\|` | `\|R\|` | `\|R\|²+\|T\|²−1` |
|--|--|--|--|
| `0.05` | `0.00098705` | `0.99999951` | `−1.6e-14` |
| `0.20` | `0.01665893` | `0.99986123` | `+1.0e-14` |
| `1.00` | `0.64142490` | `0.76718583` | `−1.4e-14` |
| `2.00` | `0.99657150` | `0.08273606` | `−3.7e-14` |
| `10.0` | `1.00000000` | `0.00000000` | `−6.5e-14` |

Unitarity to `6.6e-14`, imposed nowhere. **Reciprocity is structural**: `V` is
even in `s` to machine precision, so `R_L = R_R` and `T_LR = T_RL` identically —
a property of the geometry, not a numerical coincidence.

---

## The threshold law, and the false regression avoided

The exact static monopole response has `I₃ = ∫ds/f³ = 2/a²` and
`g = 2π²/I₃ = π²a²`. It is tempting to regress `T₀(0)` against `g`. **That would
be wrong.** A flux-normalised transmission amplitude between asymptotically flat
four-dimensional spatial ends vanishes at threshold, from the asymptotic radial
normalisation. What the static solution controls is the *coefficient* of the
low-frequency expansion:

```
|T₀(ω)| → (π/8)(aω)²
```

| `ω` | `\|T₀\|/(π/8)(aω)²` | `arg(ratio)/π` |
|--|--|--|
| `0.010` | `1.000294` | `+0.5000` |
| `0.015` | `1.000618` | `+0.5001` |
| `0.020` | `1.001043` | `+0.5001` |
| `0.030` | `1.002168` | `+0.5002` |
| `0.050` | `1.005398` | `+0.5006` |

The magnitude law is confirmed and converges. The phase of the ratio is a
**constant** `+π/2` to within `5e-3` across the band — i.e. exactly a factor of
`i`, from this module's Riccati–Hankel normalisation to `e^{±ix}`. The oracle
was stated up to the asymptotic phase convention, and that is precisely what the
residual is.

Static conductance and dynamical threshold law are two faces of one calculation
— but they are not the same number, and the module says so.

---

## No single clock offset closes the loop

Once the throat is dispersive there is no unique `τ_th`. The geometry supplies
`δ_ℓ(ω) = arg T_ℓ(ω)`, and monochromatic closure reads

```
Φ_ℓ(ω) = ω(d_A + d_B + Δ_BA) + δ_ℓ(ω) + φ_topo = 2πn
```

while a wave packet additionally needs `dΦ/dω = 0`. So:

```
Δ_phase(ω) = −(d_A+d_B) − [δ_ℓ(ω) + φ_topo − 2πn]/ω
Δ_group(ω) = −(d_A+d_B) − dδ_ℓ/dω          ← the Wigner delay
```

With `d_A + d_B = π` (antipodal on the unit `S³`):

| `ω` | `δ_ℓ` | `dδ/dω` | `Δ_phase` | `Δ_group` | gap |
|--|--|--|--|--|--|
| `0.5` | `−2.9260` | `+0.9323` | `+2.7105` | `−4.0739` | **`6.7844`** |
| `1.0` | `−2.1523` | `+1.9921` | `−0.9893` | `−5.1337` | **`4.1445`** |
| `1.5` | `−1.3339` | `+1.1324` | `−2.2523` | `−4.2740` | **`2.0217`** |
| `2.0` | `−0.9400` | `+0.5424` | `−2.6716` | `−3.6840` | **`1.0124`** |
| `3.0` | `−0.6032` | `+0.2117` | `−2.9405` | `−3.3533` | **`0.4128`** |
| `5.0` | `−0.3563` | `+0.0724` | `−3.0703` | `−3.2140` | **`0.1437`** |

**One clock offset can phase-close a monochromatic carrier while failing to
return a localised packet to the same event.** The gap is largest where the
throat disperses most and dies monotonically as it becomes transparent — which
is what identifies it as a dispersion effect rather than a numerical artefact.
`d²δ/dω²` then sets how narrow the packet must be for the monochromatic answer
to survive.

This is the quantitative form of the undeclared monochromatic assumption: PR
#216's conventions line reads *"global monochromatic form `e^{−iωt}`"* while
`projected_kernel` multiplies frequency by frequency, and a transaction is
bilinear, `T̃(Ω) = ∫dω Ã_R(ω)Ã_A(Ω−ω)`.

---

## Can the transaction complete?

PR #216 carries `Λ_ℓ(ω) = t_net η_topo e^{iωD_loop}` and `G_eff = G₀/(1−Λ)`,
calling `Λ → 1` the completed transaction. With a derived `T_ℓ` that becomes
answerable.

| `ω` | `1 − \|T\|²` |
|--|--|
| `1.50` | `7.76e-02` |
| `2.87` | `1.26e-04` |
| `4.24` | `3.35e-07` |
| `5.61` | `1.08e-09` |
| `6.97` | `3.82e-12` |

**`|Λ| = 1` requires `|T_ℓ(ω)| = 1`, and a barrier with `V > 0` everywhere has
`|T| < 1` at every finite frequency.** So `1 − Λ` cannot vanish and `G_eff` has
no true pole, unless something outside the throat supplies gain. The completed
transaction is a high-`Q` near-resonance, not an attainable one.

But `1 − |T|² ~ exp(−4.25 ω)` over twelve decades, so the `Q` of the
near-transaction grows *exponentially* with frequency and the frequency needed
for a given `Q` grows only logarithmically. Exact closure is approached in the
UV — a chronology-horizon concern rather than the benign resonance the
phenomenological ports suggested.

| `D_loop` | regime |
|--|--|
| `> 0` | ordinary delayed feedback |
| `= 0` | closed-null / marginal chronology condition |
| `< 0` | return before emission (CTC regime) |

---

## Ledger

| claim | verdict | evidence |
|--|--|--|
| an advanced leg alone evades ER non-traversability | **NO — PROVED, NOT COMPUTED** | support composition forces `d ∈ J⁺(s)` |
| PR #216's throat transfer is supplied by a geometry | **NOW IT IS** | `T_ℓ`, `R_ℓ` from the whole throat instead of port inputs |
| the throat factorises as two mouths plus a cavity | **NOT BY THIS GEOMETRY** | one smooth symmetric barrier; not imposed here |
| the static conductance `g` is `T₀(0)` | **NO — IT IS THE THRESHOLD COEFFICIENT** | `\|T₀\| → (π/8)(aω)²`, ratio `1.0003` at `ω = 0.01` |
| the traversable throat is free | **NO** | NEC violated everywhere; null-ray integral exactly `−3/(16G₅a)` |
| that NEC integral is the cost of the clock offset | **NO — DIFFERENT QUESTION** | static support vs aging history |
| one clock offset closes the loop | **NOT FOR A PACKET** | phase and group closure differ by up to `6.78` |
| `Λ = 1` can occur | **NOT AT ANY FINITE FREQUENCY** | `\|T\| < 1`; deficit `~ exp(−4.25ω)` |

---

## What this round establishes, and what it does not

PR #216 relocated the advanced wave from a bare postulate into a
Morris–Thorne–Yurtsever geometry. That was real progress. This round supplies
the geometry it was missing and finds three things: the relocation has an exact
price, the closure it assumed is frequency-dependent rather than a single tuned
constant, and its completed transaction is a high-`Q` limit rather than an
attainable pole.

**What remains postulated:** that the throat is traversable at all, and that the
mouths carry the required clock offset. The first now has an exact price. The
second needs moving-mouth dynamics this round does not do.

**Still open:** the history that produces `Δ_BA`, and whether any BAM exchange
kernel is meant to approximate the whole-throat `T_ℓ` derived here.
